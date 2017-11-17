#pragma once
//Synrad stuff, distributions and interpolation

#include "Distributions.h"

#define NUMBER_OF_DISTRO_VALUES 100

#define UPPER_LIMIT log10(20)
#define LOWER_LIMIT log10(1E-10) //-10
#define NUMBER_OF_INTEGR_VALUES 5000
//sampling SR spectrum from [1E-10...20]*critical_energy in 5000 steps

#define INTEGRAL_MODE_N_PHOTONS 1
#define INTEGRAL_MODE_SR_POWER  2

#define SYNGEN_MODE_FLUXWISE  0
#define SYNGEN_MODE_POWERWISE 1

class Material { //2-variable interpolation
public:
	//ReflectivityTable(const ReflectivityTable &copy_src); //copy constructor
	//ReflectivityTable& operator= (const ReflectivityTable & other); //assignment op
	std::vector<double> BilinearInterpolate(const double &energy, const double &angle);
	std::vector<double> energyVals, angleVals; //energy and angle values (table headers)
	std::vector<std::vector<std::vector<double>>> reflVals; //actual table values, 3 components for forward/diffuse/back reflection
	void LoadMaterialCSV(FileReader *file);
	std::string name; //For display in the program
	int GetReflectionType(const double &energy, const double &angle,const double &rnd);
	int hasBackscattering; //has forward/diffuse/backscattering/transparent pass probabilities
	void InitAngles(std::vector<std::string> data);
};

/*
class Indexes {
public:
	double inferior_flux, inferior_power;
	double superior_flux, superior_power;
	int i1;
};
*/

//double g0ki(double x, double order, int kind);
double SYNRAD_FAST(const double &x);
//double Gi(double x,int order);
//double H(double x, int order);
//double calc_polarization_percentage(double energy,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
//double find_psi_and_polarization(double x,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
//double find_chi(double psi,double gamma,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
std::tuple<double,double> find_psi_and_polarization(const double& lambda_ratios,
	const std::vector<std::vector<double>> &psi_distr, const std::vector<std::vector<double>> &parallel_polarization, const size_t& polarizationComponent); //returns psi and parallel polarization ratio
double find_chi(const double& psi, const double& gamma, const std::vector<std::vector<double>> &chi_distr);
double SYNGEN1(const double& log10LoEnergyRatio, const double& log10HiEnergyRatio,
	double& interpFluxLo,double& interpFluxHi,double& interpPowerLo,double& interpPowerHi, const bool& calcInterpolates,
	const int& generation_mode);

//Distribution2D Generate_K_Distribution(double order);
//Distribution2D Generate_G1_H2_Distribution();
//Distribution2D Generate_Polarization_Distribution(bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
Distribution2D Generate_SR_spectrum(double log10_min, double log10_max, int mode);
/*
Distribution2D K_1_3_distribution=Generate_K_Distribution(1.0/3.0);
Distribution2D K_2_3_distribution=Generate_K_Distribution(2.0/3.0);
Distribution2D integral_N_photons=Generate_SR_spectrum(LOWER_LIMIT,UPPER_LIMIT,INTEGRAL_MODE_N_PHOTONS);
Distribution2D integral_SR_power=Generate_SR_spectrum(LOWER_LIMIT,UPPER_LIMIT,INTEGRAL_MODE_SR_POWER);
Distribution2D polarization_distribution=Generate_Polarization_Distribution(true,true);
Distribution2D g1h2_distribution=Generate_G1_H2_Distribution();
*/