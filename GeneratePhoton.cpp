#include "Simulation.h"
#include "Random.h"
extern Distribution2D polarization_distribution,integral_N_photons,integral_SR_power,g1h2_distribution;

GenPhoton GeneratePhoton(int pointId, Region_mathonly *current_region, int generation_mode,BOOL recalc) { //Generates a photon from point number 'pointId'

	/* interpolation removed, wasn't useful and slowed things down
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
	double average_;

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
		result.radius = current_region->E / 0.00299792458 / result.B_ort; //Energy in GeV divided by speed of light/1E9 converted to centimeters
	else result.radius = 1.0E30;

	result.critical_energy = 2.959E-5*pow(current_region->gamma, 3) / result.radius; //becomes ~1E-30 if radius is 1E30


	if (result.critical_energy == last_critical_energy) { //check if we've calculated the results already
		result.B_factor = last_Bfactor;
		result.B_factor_power = last_Bfactor_power;
		average_ = last_average_ans;
	} else {
		result.B_factor = (exp(integral_N_photons.InterpolateY(log(current_region->energy_hi / result.critical_energy)))
			- exp(integral_N_photons.InterpolateY(log(current_region->energy_low / result.critical_energy))))
			/ exp(integral_N_photons.valuesY[integral_N_photons.size - 1]);
		result.B_factor_power = (exp(integral_SR_power.InterpolateY(log(current_region->energy_hi / result.critical_energy)))
			- exp(integral_SR_power.InterpolateY(log(current_region->energy_low / result.critical_energy))))
			/ exp(integral_SR_power.valuesY[integral_SR_power.size - 1]);
		average_ = integral_N_photons.Interval_Mean(current_region->energy_low / result.critical_energy, current_region->energy_hi / result.critical_energy);
	}
	double average = integral_N_photons.average;

	double generated_energy = SYNGEN1(current_region->energy_low / result.critical_energy, current_region->energy_hi / result.critical_energy,
		generation_mode);

	result.SR_flux=current_region->dL/(result.radius*2*PI)*current_region->gamma*4.1289E14*result.B_factor*current_region->current;
	//Total flux per revolution for electrons: 8.08E17*E[GeV]*I[mA] photons/sec
	//8.08E17 * 0.000511GeV = 4.1289E14

	if (generation_mode==SYNGEN_MODE_POWERWISE)
		result.SR_flux=result.SR_flux/generated_energy*average_;

	if (result.B_factor>0.0 && average_>VERY_SMALL) {
		result.SR_power = result.SR_flux*generated_energy*result.critical_energy*1.602189E-19*average / average_*result.B_factor_power / result.B_factor; //flux already multiplied by current
	}
	else result.SR_power = 0.0;

	//}
	double f=polarization_distribution.InterpolateY(generated_energy);
	result.g1h2=exp(g1h2_distribution.InterpolateY(generated_energy));
	double f_times_g1h2=f*result.g1h2;
		
	int retries = 0;
	do {
		result.natural_divy=find_psi(generated_energy,Sqr(current_region->gamma),f_times_g1h2,
			current_region->enable_par_polarization,current_region->enable_ort_polarization)/current_region->gamma;
	} while (result.natural_divy>current_region->psimaxY || (++retries)>1000);
	retries = 0;
	do {
		result.natural_divx=find_chi(result.natural_divy,current_region->gamma,f_times_g1h2,
			current_region->enable_par_polarization,current_region->enable_ort_polarization); //divided by sHandle->gamma inside the function
	} while (result.natural_divx>current_region->psimaxX || (++retries)>1000);

	if (rnd()<0.5) result.natural_divx *= -1;
	if (rnd()<0.5) result.natural_divy *= -1;

	//Store these calc results to speed up calculation for the same critical energy
	last_critical_energy=result.critical_energy;
	last_Bfactor=result.B_factor;
	last_Bfactor_power=result.B_factor_power;
	last_average_ans=average_;

	//return values
	result.start_dir=result.start_dir.Rotate(source->Y_local,result.natural_divx);
	result.start_dir=result.start_dir.Rotate(source->X_local,result.natural_divy);
	result.energy=generated_energy*result.critical_energy;

	return result;
}