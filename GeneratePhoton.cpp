#include "Simulation.h"
#include "Random.h"
#include "GeneratePhoton.h"
extern Distribution2D /*polarization_distribution,*/integral_N_photons, integral_SR_power/*,g1h2_distribution*/;

GenPhoton GeneratePhoton(int pointId, Region_mathonly *current_region, int generation_mode,
	std::vector<std::vector<double>> &psi_distr, std::vector<std::vector<double>> &chi_distr, BOOL recalc) { //Generates a photon from point number 'pointId'

	/* interpolation between source points removed, wasn't useful and slowed things down
	//Interpolate source point
	SATURATE(sourceId,0,(double)current_region->Points.size()-1.0000001); //make sure we stay within limits
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

	result.offset_x = Gaussian(source->sigma_x);
	result.offset_y = Gaussian(source->sigma_y);
	result.offset_divx = Gaussian(source->sigma_x_prime);
	result.offset_divy = Gaussian(source->sigma_y_prime);

	Vector offset = Vector(0, 0, 0); //choose ideal beam as origin
	offset = Add(offset, ScalarMult(source->X_local, result.offset_x)); //apply dX offset
	offset = Add(offset, ScalarMult(source->Y_local, result.offset_y)); //apply dY offset
	result.start_pos = Add(source->position, offset);

	result.start_dir = source->Z_local; //choose orbit direction as original dir, then apply offset
	result.start_dir = result.start_dir.Rotate(source->Y_local, result.offset_divx);
	result.start_dir = result.start_dir.Rotate(source->X_local, result.offset_divy);

	result.B = current_region->B(pointId, offset); //recalculate B at offset position
	result.B_par = DotProduct(result.start_dir, result.B);
	result.B_ort = sqrt(Sqr(result.B.Norme()) - Sqr(result.B_par));

	if (result.B_ort > VERY_SMALL)
		result.radius = current_region->params.E / 0.00299792458 / result.B_ort; //Energy in GeV divided by speed of light/1E9 converted to centimeters
	else result.radius = 1.0E30;

	result.critical_energy = 2.959E-5*pow(current_region->params.gamma, 3) / result.radius; //becomes ~1E-30 if radius is 1E30


	if (result.critical_energy == last_critical_energy) { //check if we've calculated the results already
		result.B_factor = last_Bfactor;
		result.B_factor_power = last_Bfactor_power;
		average_photon_energy_for_region = last_average_ans;
	}
	else {
		//what part of all photons we cover in our region (Emin...Emax)
		result.B_factor = (/*exp(*/integral_N_photons.InterpolateY(log(current_region->params.energy_hi / result.critical_energy))/*)*/
			- /*exp(*/integral_N_photons.InterpolateY(log(current_region->params.energy_low / result.critical_energy)))/*)*/
			/ /*exp(*/integral_N_photons.valuesY[integral_N_photons.size - 1]/*)*/;
		//what part of all power we cover in Emin...Emax
		result.B_factor_power = (/*exp(*/integral_SR_power.InterpolateY(log(current_region->params.energy_hi / result.critical_energy))/*)*/
			- /*exp(*/integral_SR_power.InterpolateY(log(current_region->params.energy_low / result.critical_energy)))/*)*/
			/ /*exp(*/integral_SR_power.valuesY[integral_SR_power.size - 1]/*)*/;
		average_photon_energy_for_region = Interval_Mean(current_region->params.energy_low / result.critical_energy, current_region->params.energy_hi / result.critical_energy);
	}
	//double average_photon_energy_whole_range = integral_N_photons.sum_energy / integral_N_photons.sum_photons;
	double average_photon_energy_whole_range = integral_SR_power.valuesY[integral_SR_power.size - 1] / integral_N_photons.valuesY[integral_N_photons.size - 1];

	double generated_energy = SYNGEN1(log(current_region->params.energy_low / result.critical_energy), log(current_region->params.energy_hi / result.critical_energy),
		generation_mode);

	double ratio_of_full_revolution = current_region->params.dL / (result.radius * 2 * PI);
	result.SR_flux = ratio_of_full_revolution*current_region->params.gamma*4.13104E14*result.B_factor*current_region->params.current;
	//Total flux per revolution for electrons: 8.084227E17*E[GeV]*I[mA] photons/sec
	//8.084227E17 * 0.000511GeV = 4.13104E14

	if (generation_mode == SYNGEN_MODE_POWERWISE)
		result.SR_flux *= average_photon_energy_for_region / generated_energy;

	if (result.B_factor > 0.0 && average_photon_energy_for_region > VERY_SMALL) {
		result.SR_power = result.SR_flux*(generated_energy*result.critical_energy)*average_photon_energy_whole_range / average_photon_energy_for_region; //flux already multiplied by current
		result.SR_power *= result.B_factor_power / result.B_factor; //correction for not generating on every possible energy
		result.SR_power *= 1.602189E-19; //eV per second -> J/second (Watts)
	}
	else result.SR_power = 0.0;

	//}
	//double f=polarization_distribution.InterpolateY(generated_energy);
	/*result.g1h2=exp(g1h2_distribution.InterpolateY(generated_energy));
	double f_times_g1h2=f*result.g1h2;*/

	int retries = 0;
	/*
	//debug start
	size_t numSamples = 500000;
	FileWriter *fileOut = new FileWriter("histout.csv");
	fileOut->WriteDouble(generated_energy,"\n");
	for (size_t i = 0; i < numSamples; i++)
	fileOut->WriteDouble(find_psi(generated_energy, Sqr(current_region->gamma), f_times_g1h2,
	current_region->enable_par_polarization, current_region->enable_ort_polarization),"\n");
	delete(fileOut);
	throw Error("blabla");
	*/
	//debug end

	do {
		//result.natural_divy=find_psi(generated_energy,current_region->enable_par_polarization,current_region->enable_ort_polarization)/current_region->gamma;
		result.natural_divy = find_psi(generated_energy, psi_distr) / current_region->params.gamma;
	} while (result.natural_divy > current_region->params.psimaxY || (++retries) > 1000);

	retries = 0;
	do {
		/*result.natural_divx=find_chi(result.natural_divy,current_region->gamma,
			current_region->enable_par_polarization,current_region->enable_ort_polarization); //divided by sHandle->gamma inside the function*/
		result.natural_divx = find_chi(result.natural_divy, current_region->params.gamma,
			chi_distr);
	} while (result.natural_divx > current_region->params.psimaxX || (++retries) > 1000);

	//Symmetrize distribution
	if (rnd() < 0.5) result.natural_divx *= -1;
	if (rnd() < 0.5) result.natural_divy *= -1;

	//Store these calc results to speed up calculation for the same critical energy
	last_critical_energy = result.critical_energy;
	last_Bfactor = result.B_factor;
	last_Bfactor_power = result.B_factor_power;
	last_average_ans = average_photon_energy_for_region;

	//return values
	result.start_dir = result.start_dir.Rotate(source->Y_local, result.natural_divx);
	result.start_dir = result.start_dir.Rotate(source->X_local, result.natural_divy);
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
	double sum_photons2 = integral_N_photons.InterpolateY(log(max)) - integral_N_photons.InterpolateY(log(min));
	double sum_power2 = integral_SR_power.InterpolateY(log(max)) - integral_SR_power.InterpolateY(log(min));
	if (sum_photons2 > VERY_SMALL) return sum_power2 / sum_photons2;
	else return 0.5*(min + max);
}