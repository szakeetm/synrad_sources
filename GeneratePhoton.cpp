//#include "Simulation.h"
#include "Random.h"
#include "GLApp/MathTools.h"
#include "GeneratePhoton.h"
#include "Distributions.h"
#include <vector>
#include <tuple>

extern Distribution2D /*polarization_distribution,*/integral_N_photons, integral_SR_power/*,g1h2_distribution*/;

GenPhoton GeneratePhoton(int pointId, Region_mathonly *current_region, int generation_mode,
	std::vector<std::vector<double>> &psi_distro, std::vector<std::vector<double>> &chi_distro,
	std::vector<std::vector<double>> &parallel_polarization, bool recalc) { //Generates a photon from point number 'pointId'

	/* interpolation between source points removed, wasn't useful and slowed things down
	//Interpolate source point
	Saturate(sourceId,0,(double)current_region->Points.size()-1.0000001); //make sure we stay within limits
	Trajectory_Point previousPoint=current_region->Points[(int)sourceId];
	Trajectory_Point nextPoint=current_region->Points[(int)sourceId+1];
	double overshoot=sourceId-(int)sourceId;
	Trajectory_Point source;
	source.position=Weigh(previousPoint.position,nextPoint.position,overshoot);
	source.direction=Weigh(previousPoint.direction,nextPoint.direction,overshoot);
	source.rho=Weigh(previousPoint.rho,nextPoint.rho,overshoot);
	*/

	Trajectory_Point *source = &(current_region->Points[pointId]);
	GenPhoton result;

	static double last_critical_energy, last_Bfactor, last_Bfactor_power; //to speed up calculation if critical energy didn't change
	static double last_average_ans;
	double average_photon_energy_for_region;

	if (recalc) last_critical_energy = 0.0; //force recalculation of B-factors (if, for example, region properties have changed)

	if (source->emittance_X == 0.0) { //Ideal beam
		result.offset_x = 0.0;
		result.offset_divx = 0.0;
	}
	else {
		//Generate particle emittance, around RMS emittance of the point
		double factorX, phaseX, x_unrotated, xprime_unrotated;

		//Generate point on ellipse whose axes are a,b precalculated previously
		factorX = sqrt(-2.0 * log(rnd())); //Rayleigh distribution, will be the size of the ellipse

		 phaseX = 2.0 * PI*rnd(); //We choose a point uniformly on the parametrized phase ellipse
		 x_unrotated = factorX * source->a_x * cos(phaseX); //Offset in [cm] if alpha=0
		 xprime_unrotated = factorX * source->b_x * sin(phaseX); //Divergence in [cm] if alpha=0

		//Now rotate by theta:
		result.offset_x = x_unrotated * cos(source->theta_X) - xprime_unrotated * sin(source->theta_X); //horizontal offset in [cm]
		result.offset_divx = x_unrotated * sin(source->theta_X) + xprime_unrotated * cos(source->theta_X); //horizontal divergence in [rad]
	}

	if (source->emittance_Y == 0.0) { //Ideal beam
		result.offset_y = 0.0;
		result.offset_divy = 0.0;
	}
	else {
		//Generate particle emittance, around RMS emittance of the point
		double factorY = sqrt(-2.0 * log(rnd()));

		//Generate point on ellipse whose axes are a,b precalculated previously
		double phaseY = 2 * PI*rnd(); //We choose a point uniformly on the parametrized phase ellipse
		double y_unrotated = factorY * source->a_y * cos(phaseY); //Offset in [cm] if alpha=0
		double yprime_unrotated = factorY * source->b_y * sin(phaseY); //Divergence in [cm] if alpha=0

		//Now rotate by theta:
		result.offset_y = y_unrotated * cos(source->theta_Y) - yprime_unrotated * sin(source->theta_Y); //vertical offset in [cm]
		result.offset_divy = y_unrotated * sin(source->theta_Y) + yprime_unrotated * cos(source->theta_Y); //vertical divergence in [rad]
	}

	Vector3d offset = source->X_local * result.offset_x; //apply dX offset
	offset = offset + source->Y_local * result.offset_y; //apply dY offset
	result.start_pos = source->position + offset;

	result.start_dir = source->Z_local; //choose orbit direction as original dir, then apply offset
	result.start_dir = Rotate(result.start_dir,Vector3d(0,0,0),source->Y_local, result.offset_divx);
	result.start_dir = Rotate(result.start_dir,Vector3d(0,0,0),source->X_local, - result.offset_divy);

	result.B = current_region->B(pointId, offset); //recalculate B at offset position
	result.B_par = Dot(result.start_dir, result.B);
	result.B_ort = sqrt(Sqr(result.B.Norme()) - Sqr(result.B_par));

	if (result.B_ort > VERY_SMALL)
		result.radius = current_region->params.E_GeV / 0.00299792458 / result.B_ort; //Energy in GeV divided by speed of light/1E9 converted to centimeters
	else { //Magnetic field very low, generate a 0-flux virtual photon tangent to the beam
		result.radius = 1.0E30;
	}

	result.critical_energy = 2.959E-5*pow(current_region->params.gamma, 3) / result.radius; //becomes ~1E-30 if radius is 1E30

	if (result.critical_energy == last_critical_energy) { //check if we've calculated the results already
		result.B_factor = last_Bfactor;
		result.B_factor_power = last_Bfactor_power;
		average_photon_energy_for_region = last_average_ans;
	}
	else {
		//what part of all photons we cover in our region [Emin,Emax]
		result.B_factor = (integral_N_photons.InterpolateY(log(current_region->params.energy_hi_eV / result.critical_energy),false)
			- integral_N_photons.InterpolateY(log(current_region->params.energy_low_eV / result.critical_energy),false))
			/ integral_N_photons.valuesY[integral_N_photons.size - 1];
		//what part of all power we cover in [Emin,Emax]
		result.B_factor_power = (integral_SR_power.InterpolateY(log(current_region->params.energy_hi_eV / result.critical_energy),false)
			- integral_SR_power.InterpolateY(log(current_region->params.energy_low_eV / result.critical_energy),false))
			/ integral_SR_power.valuesY[integral_SR_power.size - 1];
		average_photon_energy_for_region = Interval_Mean(current_region->params.energy_low_eV / result.critical_energy,
			current_region->params.energy_hi_eV / result.critical_energy);
	}
	
	if (result.B_factor < 1E-10) { //negligible photons inside region energy limits
		result.critical_energy = 1.0E-30;
		result.B_factor = result.B_factor_power = 0.0;
		result.natural_divx = result.natural_divy = 0.0;
		result.polarization = 0.5;
		result.SR_flux = result.SR_power = result.energy = 0.0;
		return result;
	}

	double average_photon_energy_whole_range = integral_SR_power.valuesY[integral_SR_power.size - 1] / integral_N_photons.valuesY[integral_N_photons.size - 1];

	double generated_energy = SYNGEN1(log(current_region->params.energy_low_eV / result.critical_energy), log(current_region->params.energy_hi_eV / result.critical_energy),
		generation_mode);

	int retries = 0;
	do {
		std::tie(result.natural_divy, result.polarization) = find_psi_and_polarization(generated_energy, psi_distro, parallel_polarization, current_region->params.polarizationCompIndex);
		result.natural_divy /= current_region->params.gamma;
	} while (result.natural_divy > current_region->params.psimaxY_rad || (++retries) > 1000);
	
	retries = 0;
	do {
		result.natural_divx = find_chi(result.natural_divy, current_region->params.gamma, chi_distro);
	} while (result.natural_divx > current_region->params.psimaxX_rad || (++retries) > 1000);

	//Symmetrize distribution
	if (rnd() < 0.5) result.natural_divx *= -1;
	if (rnd() < 0.5) result.natural_divy *= -1;

	//Flux and power
	double ratio_of_full_revolution = current_region->params.dL_cm / (result.radius * 2 * PI);
	result.SR_flux = ratio_of_full_revolution*current_region->params.gamma*4.13104E14*result.B_factor*current_region->params.current_mA*result.polarization;
	//Total flux per revolution for electrons: 8.084227E17*E[GeV]*I[mA] photons/sec
	//8.084227E17 * 0.000511GeV = 4.13104E14

	if (generation_mode == SYNGEN_MODE_POWERWISE)
		result.SR_flux *= average_photon_energy_for_region / generated_energy;

	if (result.B_factor > 0.0 && average_photon_energy_for_region > VERY_SMALL) {
		result.SR_power = result.SR_flux*(generated_energy*result.critical_energy)*average_photon_energy_whole_range / average_photon_energy_for_region; //flux already multiplied by current and polarization
		result.SR_power *= result.B_factor_power / result.B_factor; //correction for not generating on every possible energy, only between [E_min,E_max] of region
		result.SR_power *= 1.602189E-19; //eV per second -> J/second (Watts)
	}
	else result.SR_power = 0.0;

	//Store these calc results to speed up calculation for the same critical energy
	last_critical_energy = result.critical_energy;
	last_Bfactor = result.B_factor;
	last_Bfactor_power = result.B_factor_power;
	last_average_ans = average_photon_energy_for_region;

	//return values
	result.start_dir = Rotate(result.start_dir,Vector3d(0,0,0),source->Y_local, result.natural_divx);
	result.start_dir = Rotate(result.start_dir,Vector3d(0,0,0),source->X_local, - result.natural_divy);
	result.energy = generated_energy*result.critical_energy;

	return result;
}

double Interval_Mean(const double &min, const double &max) {
	//average of a cumulative distribution (differentiation included)

	/*int i_min = integral_N_photons.findXindex(log(min)) + 1;
	int i_max = integral_N_photons.findXindex(log(max)) + 1;
	if (i_min == 0) i_min++;
	if (i_min >= size - 1) {
	i_min = size - 2;
	i_max = size - 1;
	}
	if (i_max == 0) i_max++;

	double delta = log(max / min);
	int no_steps = i_max - i_min;
	if (i_max > i_min) delta = delta / (no_steps);

	double sum_power = 0.0;
	double sum_photons = 0.0;
	int i;
	for (i = 0; i < no_steps; i++){
	double x_lower = exp(log(min) + i*delta); //lower energy
	double x_middle = exp(log(min) + (i + 0.5)*delta); //middle energy
	double x_higher = exp(log(min) + (i + 1)*delta); //higher energy
	double exp_delta = x_higher - x_lower; //actual energy range for the next index

	double interval_dN = SYNRAD_FAST(x_middle)*exp_delta; //number of photons for the actual energy interval
	sum_photons += interval_dN;        //already calculated as valuesY[i]
	sum_power += interval_dN*x_middle; //already calculated as integral_SR_power.valuesY[i]
	}*/
	double sum_photons2 = integral_N_photons.InterpolateY(log(max),false) - integral_N_photons.InterpolateY(log(min), false);
	double sum_power2 = integral_SR_power.InterpolateY(log(max), false) - integral_SR_power.InterpolateY(log(min), false);
	if (sum_photons2 > VERY_SMALL) return sum_power2 / sum_photons2;
	else return 0.5*(min + max);
}