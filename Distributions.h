#ifndef _DISTRIBUTIONS_
#define _DISTRIBUTIONS_
#define NUMBER_OF_DISTRO_VALUES 100
#define VERY_SMALL 1.0E-30

#define UPPER_LIMIT log(20)
#define LOWER_LIMIT log(1E-10)
#define NUMBER_OF_INTEGR_VALUES 5000
//sampling SR spectrum from [1E-10...20]*critical_energy in 5000 steps

#define INTEGRAL_MODE_N_PHOTONS 1
#define INTEGRAL_MODE_SR_POWER  2

#define SYNGEN_MODE_FLUXWISE  0
#define SYNGEN_MODE_POWERWISE 1

#include "File.h" //FileReader for LoadCSV
#include <vector>

class Averages {
public:
	double average;
	double average1;
	double average_;
};

class Distribution2D {
public:
	Distribution2D(int size);
	Distribution2D::Distribution2D(const Distribution2D &copy_src); //copy constructor
	Distribution2D& operator= (const Distribution2D & other); //assignment op
	~Distribution2D();
	void Resize(size_t N);
	double InterpolateY(const double &x); //interpolates the Y value corresponding to X (allows extrapolation)
	double InterpolateX(const double &y); //(no extrapolation, first/last X values are the output limits)
	double *valuesX, *valuesY;
	int findXindex(const double &x);
	int size;
	//double sum_energy,sum_photons;
};

class Material { //2-variable interpolation
public:
	//ReflectivityTable(const ReflectivityTable &copy_src); //copy constructor
	//ReflectivityTable& operator= (const ReflectivityTable & other); //assignment op
	double Interpolate(const double &energy, const double &angle);
	std::vector<double> energyVals, angleVals; //energy and angle values (table headers)
	std::vector<std::vector<double>> reflVals; //actual table values
	void LoadCSV(FileReader *file);
	std::string name; //For display in the program
};

class Indexes {
public:
	double inferior_flux, inferior_power;
	double superior_flux, superior_power;
	int i1;
};

double g0ki(double x, double order, int kind);
double SYNRAD_FAST(const double &x);
//double Gi(double x,int order);
//double H(double x, int order);
//double calc_polarization_percentage(double energy,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
//double find_psi(double x,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
//double find_chi(double psi,double gamma,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
double find_psi(double psi, std::vector<std::vector<double>> &psi_distro);
double find_chi(double psi, double gamma, std::vector<std::vector<double>> &chi_distro);
double SYNGEN1(double x_min, double x_max, int mode);

Distribution2D Generate_K_Distribution(double order);
//Distribution2D Generate_G1_H2_Distribution();
Distribution2D Generate_LN_Distribution(); //precalculated ln(x) values for the most used energies
//Distribution2D Generate_Polarization_Distribution(bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
Distribution2D Generate_Integral(double log_min, double log_max, int mode);
/*
Distribution2D K_1_3_distribution=Generate_K_Distribution(1.0/3.0);
Distribution2D K_2_3_distribution=Generate_K_Distribution(2.0/3.0);
Distribution2D integral_N_photons=Generate_Integral(LOWER_LIMIT,UPPER_LIMIT,INTEGRAL_MODE_N_PHOTONS);
Distribution2D integral_SR_power=Generate_Integral(LOWER_LIMIT,UPPER_LIMIT,INTEGRAL_MODE_SR_POWER);
Distribution2D polarization_distribution=Generate_Polarization_Distribution(true,true);
Distribution2D g1h2_distribution=Generate_G1_H2_Distribution();
*/

template <typename T> int binary_search(double key, T A, int size);

#endif