/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/
#include "SynradDistributions.h"
#include "SynradTypes.h"
#include "GLApp/MathTools.h"
#include <tuple>
#include "Random.h"
#include "File.h"
#include <sstream>
#include "Buffer_shared.h"

//Distribution2D K_1_3_distribution = Generate_K_Distribution(1.0 / 3.0);
//Distribution2D K_2_3_distribution = Generate_K_Distribution(2.0 / 3.0);
//DistributionND SR_spectrum_CDF = Generate_SR_spectrum(LOWER_LIMIT, UPPER_LIMIT);
Distribution2D integral_N_photons = Generate_SR_spectrum(LOWER_LIMIT, UPPER_LIMIT, INTEGRAL_MODE_N_PHOTONS);
Distribution2D integral_SR_power = Generate_SR_spectrum(LOWER_LIMIT, UPPER_LIMIT, INTEGRAL_MODE_SR_POWER);
//Distribution2D polarization_distribution=Generate_Polarization_Distribution(true,true);
//Distribution2D g1h2_distribution=Generate_G1_H2_Distribution();

/*int Distribution2D::findXindex(const double &x) {
	int superior_index;
	for (superior_index = 0; valuesX[superior_index] < x && superior_index < size; superior_index++); //replace by binary search
	return superior_index;
	}*/


//double g0ki(double x, double order, int kind) {
	/*ported from CALCF1.PAS:function g0ki(x,ord:realt1;kind:integer):realt1;

	{ Adapted from B. Diviacco, Sincrotrone Trieste, private comm. }
	{ Calculates the modified Bessel functions K1/3, and K2/3 for obtaining }
	{ the degree of polarization, and the integral of K5/3, used to obtain the }
	{ number of SR photons. In this case it's equivalent to SYNRAD_, but slower.}
	*/
/*
	double h1, g0, r1, q1, q2, s1, s2, t1, xs1; //absolutely no idea what these variables are or how this function works :(
	h1 = 0.5;
	g0 = 0.0;
	r1 = 0.0;
	do {
		r1 = r1 + 1.0;
		q1 = exp(r1*h1);
		q2 = exp(order*r1*h1);
		s1 = (q1 + 1.0 / q1) / 2.0;
		s2 = (q2 + 1.0 / q2) / 2.0;
		xs1 = x*s1;
		if (kind == 0) t1 = exp(-xs1)*s2; // kind=0 ---> calculates the functions 
		else t1 = exp(-xs1)*s2 / s1;     // kind=1 ---> calculates the integral  
		g0 = g0 + t1;
	} while (t1 > 1.0E-6);
	return h1*(exp(-x) / 2.0 + g0);
	*/
	/*double sum,r1,q1,q2,s1,s2,increment; //absolutely no idea what these variables are or how this function works :(
	sum=0.0;
	x=0.0;
	do {
	x+=1.0;
	q1=exp(0.5*x);
	q2=pow(q1,order);
	s1=(q1+1.0/q1)/2.0;
	s2=(q2+1.0/q2)/2.0;
	if (kind==0) increment=exp(-x*s1)*s2; // kind=0 ---> calculates the functions
	else increment=exp(-x*s1)*s2/s1;     // kind=1 ---> calculates the integral
	sum+=increment;
	} while (increment>1.0E-6);
	return 0.5*(exp(-x)/2.0+sum);*/
//}

/*
double Gi(double x,int order) {
//pascal code: Gi:=exp(ord*ln(x))*g0ki(x,5.0/3.0,1);
return pow(x,order)*g0ki(x,5.0/3.0,1);
}
*/
/*
double H(double x, int order) {
//pascal code: H:=exp(ord*ln(x))*Sqr(g0ki(x/2.0,2.0/3.0,0));
return pow(x,order)*pow(g0ki(x/2.0,2.0/3.0,0),2);
}
*/

/*
Distribution2D Generate_K_Distribution(double order){
	//Gives K_order[x] (previously gave the natural logarithm)
	Distribution2D result(NUMBER_OF_DISTRO_VALUES);
	//double stepDelta = log(UPPER_LIMIT / LOWER_LIMIT) / NUMBER_OF_DISTRO_VALUES;
	double stepDelta = (UPPER_LIMIT - LOWER_LIMIT) / NUMBER_OF_DISTRO_VALUES;
	for (int i = 0; i < NUMBER_OF_DISTRO_VALUES; i++) {
		//double x = LOWER_LIMIT*exp(i*stepDelta);
		double x = exp(LOWER_LIMIT + i*stepDelta);
		result.valuesX[i] = x;
		result.valuesY[i] = g0ki(x, order, 0); //previously took the logarithm
		//note: values go from 0 to 99 (in original code: 1 to 100)
	}
	return result;
}
*/

/*Distribution2D Generate_G1_H2_Distribution(){
Distribution2D result(NUMBER_OF_DISTRO_VALUES);
double stepDelta=log(UPPER_LIMIT/LOWER_LIMIT)/NUMBER_OF_DISTRO_VALUES;
for (int i=0;i<NUMBER_OF_DISTRO_VALUES;i++) {
double x=LOWER_LIMIT*exp(i*stepDelta);
result.valuesX[i]=x;
result.valuesY[i]=log(Max((Gi(x,1)/H(x,2)),VERY_SMALL));
//note: values go from 0 to 99 (in original code: 1 to 100)
}
return result;
}*/

Distribution2D Generate_SR_spectrum(double log10_min, double log10_max, int mode){
	//log10_min : Log10( Minimum_energy / Critical_energy )
	//log10_max : Log10( Maximum_energy / Critical_energy )

	/* ported from CALCF1.PAS:function integral(x1,x2:realt1):realt1;

	{ Calculates the two real vectors integ[i,1], for the numbers of SR photons,  }
	{ integ[i,2] for the SR power, 						      }
	{ and integ[i,3], for the corresponding energies, for generating SR photons   }
	{ according to the real distribution, see fig. 2 ref. G.K. Green for instance }
	*/

	//return value: Distribution2D a.k.a. pairs of double,double
	//Scale X: log(ratio of energy / critical energy), goes from -10 to +2
	//Scale Y: CDF (integral) of SR spectrum (flux or power) from log(E/E_crit)=-10 to X, goes from 0 to 1
	
	double delta, log10_delta, interval_dN, sum_photons, sum_power, mean_photons, x_lower, x_middle, x_higher;

	int i;
	Distribution2D result; result.Resize(NUMBER_OF_INTEGR_VALUES);
	//DistributionND result2;

	delta = (log10_max - log10_min) / (double)NUMBER_OF_INTEGR_VALUES;
	sum_photons = sum_power =  0.0;

	for (i = 0; i < NUMBER_OF_INTEGR_VALUES; i++) {
		x_lower = Pow10(log10_min + i*delta);          //lower energy of bin, in units of E_crit
		x_higher = Pow10(log10_min + (i + 1.0)*delta); //higher energy of bin, in units of E_crit
		x_middle = (x_lower + x_higher) / 2.0;     //average energy of bin
		log10_delta = x_higher - x_lower;            //bin energy range, in units of E_crit
		mean_photons = (SYNRAD_FAST(x_lower) + SYNRAD_FAST(x_higher)) / 2.0;

		interval_dN = mean_photons*log10_delta; //number of photons for the actual interval, different averaging
		
		if (mode == INTEGRAL_MODE_N_PHOTONS) {
			sum_photons += interval_dN; //total integrated flux from log10(E/E_crit)=-10 to this interval
			result.SetPair(i,log10(x_middle), sum_photons);
		}
		else {
			sum_power += interval_dN*x_middle; //total integrated power from log10(E/E_crit)=-10 to this interval (number of photons * average energy: energy of the interval)
			result.SetPair(i,log10(x_middle), sum_power);
		}
		//result.AddPair(log10(x_middle), { sum_photons,sum_power });
	}
	return result;
}

/*
Distribution2D Generate_Polarization_Distribution(bool calculate_parallel_polarization, bool calculate_orthogonal_polarization) {
//ported from CALCF1.PAS - procedure findpolarization;

//{ Calculates the polarization percentages vs. normalized photon energy }
//{ in the interval 1.0E-10 < (energy/critical_en.) < 100, as partially }
//{ reported in fig. 14, ref. G.K. Green.}
//{ An accuracy of better than 0.5% is usually obtained.}

double energy,delta;
Distribution2D result(NUMBER_OF_DISTRO_VALUES);

delta=log(UPPER_LIMIT/LOWER_LIMIT)/NUMBER_OF_DISTRO_VALUES;
for (int i=0;i<NUMBER_OF_DISTRO_VALUES;i++) {
energy=LOWER_LIMIT*exp(i*delta);
result.valuesX[i]=energy;
result.valuesY[i]=calc_polarization_percentage(energy,calculate_parallel_polarization,calculate_orthogonal_polarization); //F_pol percentage
}
result.valuesY[99]=0.997;
return result;
}
*/

/*
double calc_polarization_percentage(double energy,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization){
//ported from CALCF1.PAS: function find_polarization(x:realt1):realt1;
//Besself[6] removed as distribution is never used later

//{ Calculates the degree of polarization (linear and orthogonal) for a }
//{ given SR photon energy }
//{ Ref. G.K. Green, page 1-12, and fig. 15 }

//if (energy>20.0) __debugbreak();
double F_parallel,F_orthogonal;
double delta,csi,X1;
int index;

F_parallel=0.0;
F_orthogonal=VERY_SMALL;
delta=4.344*Gi(energy,1)/H(energy,2)/NUMBER_OF_DISTRO_VALUES;
index=0;

for (int i=0;i<NUMBER_OF_DISTRO_VALUES;i++) {
X1=i*delta;
csi=1.0+Sqr(X1);
csi=energy/2.0*csi*sqrt(csi);
//csi=log(csi);

if (calculate_parallel_polarization)
F_parallel+=Sqr((1.0+Sqr(X1))*exp(K_2_3_distribution.InterpolateY(csi)));
else F_parallel=0.0;

if (calculate_orthogonal_polarization)
F_orthogonal+=(1.0+Sqr(X1))*Sqr(X1*exp(K_1_3_distribution.InterpolateY(csi)));
else F_orthogonal=0.0;
}
return F_parallel/(F_parallel+F_orthogonal);
}
*/

/*double find_psi_and_polarization(double lambda_ratios,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization){
//vertical

//ported from FINDPRCS.PAS: function find_psi_and_polarization(x,maxangle:realt1):realt1;

//{ Derived from procedure find_polarization ... }

//RETURNS GAMMA*PSI

double F_parallel,F_orthogonal;
double delta,ksi,gamma_psi,factor_psi,seed;
Distribution2D local_polarization_integral(NUMBER_OF_DISTRO_VALUES+1);
F_parallel=0.0;
F_orthogonal=VERY_SMALL;
delta = -0.5*pow(lambda_ratios,0.35) / NUMBER_OF_DISTRO_VALUES; //At -3*lambdaratios^0.35 we encompass 10 orders of magnitude

//See eq. (2) and chapter 6 (polarization) of http://server2.phys.uniroma1.it/gr/lotus/Mariani_carlo/didattica/cap_01_Mobilio.pdf

//First we construct the integrated vertical distributions for a given energy (expressed as lambda_ratio),
//for orthogonal and parallel polarization, then we generate a random number according to this distribution
for (int i=1;i<=NUMBER_OF_DISTRO_VALUES;i++) {
gamma_psi = (i - 1)*delta; //actual gamma*psi
double one_plus_gammapsisquare = 1.0 + Sqr(gamma_psi); //1+(gamma*psi)^2
ksi = lambda_ratios / 2.0*pow(one_plus_gammapsisquare, 1.5); //ksi=(lambda_crit/lambda/2)(1+(gamma*psi)^2)^3/2

if (calculate_parallel_polarization)
F_parallel = Sqr(one_plus_gammapsisquare*K_2_3_distribution.InterpolateY(ksi)); //used to be exp(K_2_3)
else F_parallel=0.0;

if (calculate_orthogonal_polarization)
F_orthogonal = one_plus_gammapsisquare*Sqr(gamma_psi*K_1_3_distribution.InterpolateY(ksi)); //used to be exp(K_1_3)
//tested: one_plus_square_X1 shouldn't be squared (as opposed to typo in http://server2.phys.uniroma1.it/gr/lotus/Mariani_carlo/didattica/cap_01_Mobilio.pdf)
else F_orthogonal=0.0;

local_polarization_integral.valuesX[i - 1] = gamma_psi; //Will also set valuesX[0]=0
if (i>1) local_polarization_integral.valuesY[i]=local_polarization_integral.valuesY[i-1]+F_parallel+F_orthogonal; //integrate full polarization
else local_polarization_integral.valuesY[i]=F_parallel+F_orthogonal;
}
local_polarization_integral.valuesY[0]=0.0;
factor_psi=1.0;
seed=rnd()*factor_psi*local_polarization_integral.valuesY[NUMBER_OF_DISTRO_VALUES-1];
return local_polarization_integral.InterpolateX(seed);
}*/

std::tuple<double,double> find_psi_and_polarization(const double& lambda_ratios, const std::vector<std::vector<double>> &psi_distro,
	const std::vector<std::vector<double>> &parallel_polarization, const size_t& polarizationComponent) {
	
	//returns gamma*psi
	double lambda_relative = log10(lambda_ratios);
	double lambda_index = (lambda_relative + 10) / 0.1; //digitized for -10..+2 with delta=0.01
	int lambda_lower_index = (int)(lambda_index);
	double lambda_overshoot = lambda_index - (double)lambda_lower_index;

	double lookup = rnd();
	int foundAngle = 0;
	/*double interpolated_CDF;
	do { //to replace by binary search
	foundAngle++;
	interpolated_CDF = psi_distro[lambda_lower_index][foundAngle] + lambda_overshoot*(psi_distro[lambda_lower_index + 1][foundAngle] - psi_distro[lambda_lower_index][foundAngle]);

	} while (interpolated_CDF < lookup);*/

	//Binary search
	int imin = 0;
	int size = (int)psi_distro[0].size();
	int imax = size;
	double interpolated_CDF_lower, interpolated_CDF_higher;
	// continue searching while [imin,imax] is not empty
	while (imin <= imax)
	{
		// calculate the midpoint for roughly equal partition
		int imid = (imin + imax) / 2;
		interpolated_CDF_lower = Weigh(psi_distro[lambda_lower_index][imid], psi_distro[lambda_lower_index + 1][imid], lambda_overshoot);
		interpolated_CDF_higher = Weigh(psi_distro[lambda_lower_index][imid + 1], psi_distro[lambda_lower_index + 1][imid + 1], lambda_overshoot);
		if (imid == size - 1 || (interpolated_CDF_lower < lookup && lookup < interpolated_CDF_higher)) {
			// key found at index imid
			foundAngle = imid + 1; //will be lowered by 1
			break;
		}
		// determine which subarray to search
		else if (interpolated_CDF_lower < lookup)
			// change min index to search upper subarray
			imin = imid + 1;
		else
			// change max index to search lower subarray
			imax = imid - 1;
	}
	//TO DO: Treat not found errors
	int psi_lower_index = foundAngle - 1;
	//double previous_interpolated_CDF = psi_distro[lambda_lower_index][psi_lower_index] + lambda_overshoot*(psi_distro[lambda_lower_index + 1][psi_lower_index] - psi_distro[lambda_lower_index][psi_lower_index]);
	double psi_overshoot = (lookup - interpolated_CDF_lower) / (interpolated_CDF_higher - interpolated_CDF_lower);
	double psi = (((double)psi_lower_index + psi_overshoot)*0.005) * (4.0 / pow(lambda_ratios, 0.35)); //psi_relative=1 corresponds to psi=4/lambda_ratios^0.35
	
	double polarization;
	if (polarizationComponent == 0) {
		polarization = 1.0;
	}
	else {
		double polarization_pdf_lower = Weigh(parallel_polarization[lambda_lower_index][psi_lower_index], parallel_polarization[lambda_lower_index + 1][psi_lower_index], lambda_overshoot);
		double polarization_pdf_higher = Weigh(parallel_polarization[lambda_lower_index][psi_lower_index + 1], parallel_polarization[lambda_lower_index + 1][psi_lower_index + 1], lambda_overshoot);
		polarization = Weigh(polarization_pdf_lower, polarization_pdf_higher, psi_overshoot);
		if (polarizationComponent == 2) polarization = 1.0 - polarization; //orthogonal component
	}
	return std::tie(psi,polarization);
	//return 0.15;
}

/*double find_chi(double psi,double gamma_square,double f_times_g1h2,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization) {
//horizontal

//ported from FINDPRCS.PAS - function find_chi(psi:realt1):realt1;

//RETURNS ANGLE

double delta;
//double j1;
double F_parallel,F_orthogonal;
double chi,factor;
int i;
Distribution2D local_polarization_integral(NUMBER_OF_DISTRO_VALUES+1);

F_parallel=0.0;
F_orthogonal=VERY_SMALL;
//delta=4.344*f_times_g1h2/NUMBER_OF_DISTRO_VALUES/gamma_square;
delta = 10.0*f_times_g1h2 / NUMBER_OF_DISTRO_VALUES / gamma_square;

//{ delta:=pi/100.0/gamma_square;}
//{ delta gives the angular range (for traj.params.gamma*Psi) where the vertical }
//{ angle of emission, Psi, should be contained (actually, traj.params.gamma*Psi...) }

double cos_psi=cos(psi);
double sin_psi=sin(psi);
double sin_chi,cos_chi;
for (i=1;i<=NUMBER_OF_DISTRO_VALUES;i++){
chi=(i-0.5)*delta;
sin_chi=sin(chi);
cos_chi=cos(chi);

//{1  teta:=acos(cospsi*coschi);}
//{2  teta:=asin(sqrt(Sqr(sinpsi)+Sqr(cospsi*sinchi)));}
//{ Version 1 and 2 of teta=teta(psi,chi) generate round-off errors }

double theta=atan(sqrt(Sqr(sin_psi)+Sqr(cos_psi*sin_chi)))/(cos_psi*cos_chi);
double gamma_times_theta=gamma_square*theta;
double gamma_times_theta_square=Sqr(gamma_times_theta);
double sin_theta=sin(theta);
double sin_fi=sin_psi/sin_theta;
double cos_fi=cos_psi*sin_chi/sin_theta;
//factor = exp(-6.0*log(1.0 + gamma_times_theta_square));
factor = pow((1.0 + gamma_times_theta_square),-6.0);
//{      The 6.0 in the exp agrees with Coisson's results ... }
if (calculate_parallel_polarization)
F_parallel=Sqr(1.0-gamma_times_theta_square+2.0*Sqr(gamma_times_theta*sin_fi))*factor;
else F_parallel=0.0;

if (calculate_orthogonal_polarization)
F_orthogonal=Sqr(2.0*gamma_times_theta_square*sin_fi*cos_fi)*factor;
else F_orthogonal=0.0;

//here i goes from 1 to 100!
local_polarization_integral.valuesX[i]=chi;
if (i>1) local_polarization_integral.valuesY[i]=local_polarization_integral.valuesY[i-1]+F_parallel+F_orthogonal;
else local_polarization_integral.valuesY[i]=F_parallel+F_orthogonal;
}
local_polarization_integral.valuesX[0] = local_polarization_integral.valuesY[0] = 0.0;
double seed=rnd()*local_polarization_integral.valuesY[NUMBER_OF_DISTRO_VALUES-1];
return local_polarization_integral.InterpolateX(seed);
}*/

/*
double find_chi(double psi, double gamma, bool calculate_parallel_polarization, bool calculate_orthogonal_polarization) {
//horizontal

//ported from FINDPRCS.PAS - function find_chi(psi:realt1):realt1;

//RETURNS ANGLE
double delta;
//double j1;
double F_parallel, F_orthogonal;
double chi, factor;
int i;
Distribution2D local_polarization_integral(NUMBER_OF_DISTRO_VALUES + 1);

F_parallel = 0.0;
F_orthogonal = VERY_SMALL;
//delta=4.344*f_times_g1h2/NUMBER_OF_DISTRO_VALUES/gamma_square;
//delta = 4 * PI / gamma / NUMBER_OF_DISTRO_VALUES; //+- 4PI/gamma contains ten orders of magnitude
delta = PI / gamma / NUMBER_OF_DISTRO_VALUES; //+- 4PI/gamma contains ten orders of magnitude

double cos_psi = cos(psi);
double sin_psi = sin(psi);
double sin_chi, cos_chi;
for (i = 1; i <= NUMBER_OF_DISTRO_VALUES; i++){
if ((double)i < (0.9*(double)NUMBER_OF_DISTRO_VALUES)) {
chi = i*delta*0.5 / 0.9; //reserve 90% of points for the first half of the max. angle
}
else {
chi = (double)NUMBER_OF_DISTRO_VALUES*0.9*delta*0.5 / 0.9 + (i - (double)NUMBER_OF_DISTRO_VALUES*0.9)*0.5 / 0.1*delta; //reserve 10% for the second half
}
sin_chi = sin(chi);
cos_chi = cos(chi);

//{1  teta:=acos(cospsi*coschi);}
//{2  teta:=asin(sqrt(Sqr(sinpsi)+Sqr(cospsi*sinchi)));}
//{ Version 1 and 2 of teta=teta(psi,chi) generate round-off errors }

double theta = atan(sqrt(Sqr(sin_psi) + Sqr(cos_psi*sin_chi))) / (cos_psi*cos_chi);
double gamma_times_theta = gamma*theta;
double gamma_times_theta_square = Sqr(gamma_times_theta);
double sin_theta = sin(theta);
double sin_fi = sin_psi / sin_theta;
double cos_fi = cos_psi*sin_chi / sin_theta;
//factor = exp(-6.0*log(1.0 + gamma_times_theta_square));
factor = pow((1.0 + gamma_times_theta_square), -6.0);
//{      The 6.0 in the exp agrees with Coisson's results ... }
if (calculate_parallel_polarization)
F_parallel = Sqr(1.0 - gamma_times_theta_square + 2.0*Sqr(gamma_times_theta*sin_fi))*factor;
else F_parallel = 0.0;

if (calculate_orthogonal_polarization)
F_orthogonal = Sqr(2.0*gamma_times_theta_square*sin_fi*cos_fi)*factor;
else F_orthogonal = 0.0;

//here i goes from 1 to 100!
local_polarization_integral.valuesX[i] = chi;
if (i > 1) local_polarization_integral.valuesY[i] = local_polarization_integral.valuesY[i - 1] + F_parallel + F_orthogonal;
else local_polarization_integral.valuesY[i] = F_parallel + F_orthogonal;
}
local_polarization_integral.valuesX[0] = local_polarization_integral.valuesY[0] = 0.0;
double seed = rnd()*local_polarization_integral.valuesY[NUMBER_OF_DISTRO_VALUES - 1];
return local_polarization_integral.InterpolateX(seed);
}
*/
double find_chi(const double& psi, const double& gamma, const std::vector<std::vector<double>> &chi_distro) {

	double psi_index, chi_lower, chi_higher, chi;
	double psi_relative = log10(abs(psi)*(gamma / 10000.0)); //distributions are digitized for gamma=10000, and sampled logarithmically
	if (psi_relative < -7.0) {
		psi_index = 0; //use lowest angle
	}
	else {
		psi_index = (psi_relative + 7.0) / 0.02; //sampled from -7 to 0 with delta=0.02
	}

	int psi_lower_index = (int)(psi_index); //digitized for -2PI/10 .. +2PI/10 with delta=0.0025
	double psi_overshoot = psi_index - (double)psi_lower_index;

	double lookup = rnd();

	int foundAngle = 0;
	/*interpolated_CDF = 0.0;
	do { //to replace by binary search
	previous_interpolated_CDF = interpolated_CDF;
	foundAngle++;
	interpolated_CDF = Weigh(chi_distro[foundAngle][psi_lower_index], chi_distro[foundAngle][psi_lower_index + 1], psi_overshoot);
	} while (interpolated_CDF < lookup);*/

	//Binary search
	int imin = 0;
	int size = (int)chi_distro.size();
	int imax = size;
	double interpolated_CDF_lower, interpolated_CDF_higher;
	// continue searching while [imin,imax] is not empty
	while (imin <= imax)
	{
		// calculate the midpoint for roughly equal partition
		int imid = (imin + imax) / 2;
		interpolated_CDF_lower = Weigh(chi_distro[imid][psi_lower_index], chi_distro[imid][psi_lower_index + 1], psi_overshoot);
		interpolated_CDF_higher = Weigh(chi_distro[imid + 1][psi_lower_index], chi_distro[imid + 1][psi_lower_index + 1], psi_overshoot);
		if (imid == size - 1 || (interpolated_CDF_lower < lookup && lookup < interpolated_CDF_higher)) {
			// key found at index imid
			foundAngle = imid + 1; //will be lowered by 1
			break;
		}
		// determine which subarray to search
		else if (interpolated_CDF_lower < lookup)
			// change min index to search upper subarray
			imin = imid + 1;
		else
			// change max index to search lower subarray
			imax = imid - 1;
	}
	//TO DO: Treat not found errors

	int chi_lower_index = foundAngle - 1;

	double next_interpolated_CDF;
	if (foundAngle == chi_distro.size() - 1) next_interpolated_CDF = 1.0;
	else next_interpolated_CDF = Weigh(chi_distro[chi_lower_index + 2][psi_lower_index], chi_distro[chi_lower_index + 2][psi_lower_index + 1], psi_overshoot);
	//double chi_overshoot = (log10(lookup)-log10(previous_interpolated_CDF))/ (log10(interpolated_CDF) - log10(previous_interpolated_CDF));
	//double chi_overshoot = (lookup - previous_interpolated_CDF) / (interpolated_CDF - previous_interpolated_CDF); //linear scale

	if (chi_lower_index == 0) {
		chi_lower = 0;
		chi_higher = 1.0964782E-7;
		chi = Weigh(chi_lower, chi_higher, rnd()) / (gamma / 10000.0);
		//chi = 0;
	}
	else {

		double a = pow(10, -7.0 + ((double)chi_lower_index + 0.0)*0.02) / (gamma / 10000.0);
		double b = a*1.04712854805; /* pow(10, -7.0 + ((double)chi_lower_index + 1.0)*0.02) / (gamma / 10000.0);*/ //1.047=10^0.02
		double c = b*1.04712854805; /*pow(10, -7.0 + ((double)chi_lower_index + 2.0)*0.02) / (gamma / 10000.0);*/
		double FA = interpolated_CDF_lower;
		double FB = interpolated_CDF_higher;
		double FC = next_interpolated_CDF;
		
		chi = QuadraticInterpolateX(lookup, a, b, c, FA, FB, FC);
		//chi = Weigh(a, b,rnd()); //inverse linear interpolation

	}
	return chi;
}

double SYNGEN1(const double& log10LoEnergyRatio, const double& log10HiEnergyRatio,
	double& interpFluxLo, double& interpFluxHi, double& interpPowerLo, double& interpPowerHi, const bool& calcInterpolates,
	const int& generation_mode) {
	/*
	Originally called SYNGEN1.
	- Determines the CDF values belonging to log10_x_min and log10_x_max (they are expressed in E/E_crit)
	- Generates a random number between the two CDF values
	- Interpolates the energy belonging to the generated number

	{ Generates a random normalized SR photon energy in (xmin,xmax) from a        }
	{ cumulative distribution given by the Vector3d array integ[i,1] calculated by  }
	{ the procedure integral before calling SYNGEN1.                              }
	{ Makes a linear interpolation, on a log-log scale, of integ[i-1,1] and       }
	{ integ[i,1], in the energy interval (integ[i-1,3],integ[i,3])                }
	{ This function is generally slower than SYNGEN (ref. H. Burkhardt CERN/SL,   }
	{ LEP Note 632, 17 Dec. 1990, but allows an upper limit for SR photon energy. }
	{ For a fixed critical photon energy (i.e. constant radius), SYNGEN1 is faster}
	{ than SYNGEN. }
	*/

	/*
	double r1,m,value;
	int i;

	Indexes indexes=find_indexes(x_min,x_max,true); //subst. i1
	*/
	if (calcInterpolates) {
		interpFluxLo = integral_N_photons.InterpolateY(log10LoEnergyRatio, false);
		interpFluxHi = integral_N_photons.InterpolateY(log10HiEnergyRatio, false);
		interpPowerLo = integral_SR_power.InterpolateY(log10LoEnergyRatio, false);
		interpPowerHi = integral_SR_power.InterpolateY(log10HiEnergyRatio, false);
	}

	double generated_energy;
	if (generation_mode == SYNGEN_MODE_FLUXWISE) {
		double generated_flux = Weigh(interpFluxLo, interpFluxHi, rnd()); //uniform distribution between flux_min and flux_max
		generated_energy = Pow10(integral_N_photons.InterpolateX(generated_flux,false));
	}
	else { //Powerwise
		double generated_power = Weigh(interpPowerLo, interpPowerHi, rnd()); //uniform distribution between flux_min and flux_max
		generated_energy = Pow10(integral_SR_power.InterpolateX(generated_power,false));
	}
	return generated_energy;
}

double SYNRAD_FAST(const double &x) {
	/*
	{ Adapted from H.H. Umstaetter, CERN/PS/SM/81-13. }
	{ Works as g0ki(x,5/3,1), but about 2.5x faster }

	x = E/E_crit
	*/
	double Y, Z, A, B, P, Q;

	if (x < 6.0) {

		Z = Sqr(x) / 16.0 - 2.0;
		A = +0.0000000001;
		B = Z*A + 0.0000000023;
		A = Z*B - A + 0.0000000813;
		B = Z*A - B + 0.0000024575;
		A = Z*B - A + 0.0000618126;
		B = Z*A - B + 0.0012706638;
		A = Z*B - A + 0.0209121680;
		B = Z*A - B + 0.2688034606;
		A = Z*B - A + 2.6190218379;
		B = Z*A - B + 18.6525089687;
		A = Z*B - A + 92.9523266592;
		B = Z*A - B + 308.1591941313;
		A = Z*B - A + 644.8697965824;
		P = 0.5*Z*A - B + 414.5654364883;
		A = +0.0000000012;
B = Z*A + 0.0000000391;
A = Z*B - A + 0.0000011060;
B = Z*A - B + 0.0000258145;
A = Z*B - A + 0.0004876869;
B = Z*A - B + 0.0072845620;
A = Z*B - A + 0.0835793546;
B = Z*A - B + 0.7103136120;
A = Z*B - A + 4.2678026127;
B = Z*A - B + 17.0554078580;
A = Z*B - A + 41.8390348678;
Q = 0.5*Z*A - B + 28.4178737436;
Y = exp((2.0 / 3.0)*log(x));
//Y = pow(x, 2.0 / 3.0); //new
return (P / Y - Q*Y - 1.0)*1.8137993642;
	}
	else if (x < 80.0)
	{
		Z = 20.0 / x - 2.0;
		A = +0.0000000001;
		B = Z*A - 0.0000000004;
		A = Z*B - A + 0.0000000020;
		B = Z*A - B - 0.0000000110;
		A = Z*B - A + 0.0000000642;
		B = Z*A - B - 0.0000004076;
		A = Z*B - A + 0.0000028754;
		B = Z*A - B - 0.0000232125;
		A = Z*B - A + 0.0002250532;
		B = Z*A - B - 0.0028763680;
		A = Z*B - A + 0.0623959136;
		P = 0.5*Z*A - B + 1.0655239080;
		return P*sqrt(1.57079632679 / x) / exp(x);
	}
	else return 0.0;
}

void Material::LoadMaterialCSV(FileReader *file) {
	bool angleInit = true;
	bool energyInit = true;
	hasBackscattering = false;
	size_t lineCount = 0;
	while (!file->IsEof()) {
		std::string line = file->ReadLine(); lineCount++;
		std::stringstream   lineStream(line);
		std::string         cell;

		std::vector<std::string>    m_data;
		while (std::getline(lineStream, cell, ','))
		{
			m_data.push_back(cell);
		}
		if (angleInit) { //first line, load angles
			InitAngles(m_data);
			angleInit = false;
		}
		else if (m_data[0] == "diffuse" || m_data[0] == "back" || m_data[0] == "transparent") {
			hasBackscattering = true;
			energyInit = false;
			lineCount = 0;
			continue; //skip line
		}
		else { //regular row with values
			if (energyInit) {
				energyVals.push_back(atof(m_data[0].c_str()));
			}
			std::vector<std::vector<double>> newRow;
			for (size_t i = 1; i < m_data.size(); i++) {
				if (energyInit) { //Forward reflection, construct table
					std::vector<double> newCell;
					newCell.push_back(atof(m_data[i].c_str()));
					newRow.push_back(newCell);
				}
				else { //table already exist, add new components
					reflVals[lineCount - 1][i - 1].push_back(atof(m_data[i].c_str()));
				}
			}
			if (energyInit) reflVals.push_back(newRow);
		}
	}
}

void Material::InitAngles(std::vector<std::string> data) {
	angleVals = std::vector<double>(); //Reset
	for (size_t i = 1; i < data.size(); i++) {
		angleVals.push_back(std::stod(data[i]));
	}
}

std::vector<double> Material::BilinearInterpolate(const double &energy, const double &angle) {
	//Bilinear interpolation in a table where each cell is a list of probabilities for forward/diffuse/backscattering/transparent pass
	//Logarithmic values (both energy and angle), no extrapolation

	int angleLowerIndex = my_lower_bound(angle, angleVals);
	int energyLowerIndex = my_lower_bound(energy, energyVals);

	double angleOvershoot, angleDelta, energyOvershoot, energyDelta;
	bool interpolateAngle, interpolateEnergy;

	//Treat cases where out of table:
	if (angleLowerIndex == -1) {
		angleLowerIndex = 0;
		interpolateAngle = false;
	}
	else if (angleLowerIndex == (angleVals.size() - 1)) {
		interpolateAngle = false;
	}
	else {
		angleOvershoot = log10(angle) - log10(angleVals[angleLowerIndex]);
		angleDelta = log10(angleVals[angleLowerIndex + 1]) - log10(angleVals[angleLowerIndex]);
		interpolateAngle = true;
	}

	if (energyLowerIndex == -1) {
		energyLowerIndex = 0;
		interpolateEnergy = false;
	}
	else if (energyLowerIndex == (energyVals.size() - 1)) {
		interpolateEnergy = false;
	}
	else {
		energyOvershoot = log10(energy) - log10(energyVals[energyLowerIndex]);
		energyDelta = log10(energyVals[energyLowerIndex + 1]) - log10(energyVals[energyLowerIndex]);
		interpolateEnergy = true;
	}

	std::vector<double> interpRefl;
	for (size_t comp = 0; comp < reflVals[energyLowerIndex][angleLowerIndex].size(); comp++) {
		double interpolatedReflForLowerAngle,interpolatedReflForHigherAngle,componentReflValue;
		
		interpolatedReflForLowerAngle = 
			interpolateEnergy
			? Weigh(reflVals[energyLowerIndex][angleLowerIndex][comp], reflVals[energyLowerIndex + 1][angleLowerIndex][comp], energyOvershoot / energyDelta)
			: reflVals[energyLowerIndex][angleLowerIndex][comp];
		if (interpolateAngle) {
			interpolatedReflForHigherAngle = 
				interpolateEnergy
				? Weigh(reflVals[energyLowerIndex][angleLowerIndex + 1][comp], reflVals[energyLowerIndex + 1][angleLowerIndex + 1][comp], energyOvershoot / energyDelta)
				: reflVals[energyLowerIndex][angleLowerIndex + 1][comp];
			componentReflValue = Weigh(interpolatedReflForLowerAngle, interpolatedReflForHigherAngle, angleOvershoot / angleDelta);
		}
		else {
			componentReflValue = interpolatedReflForLowerAngle;
		}
		Saturate(componentReflValue, 0.0, 1.0);
		interpRefl.push_back(componentReflValue);
	}
	return interpRefl;
}

int Material::GetReflectionType(const double &energy, const double &angle, double const &rnd) {
	std::vector<double> components = BilinearInterpolate(energy, angle);
	if (rnd < components[0]) return REFL_FORWARD; //forward reflection
	else if (hasBackscattering) {
		if (rnd < (components[0] + components[1])) return REFL_DIFFUSE; //diffuse reflection
		else if (rnd < (components[0] + components[1] + components[2])) return REFL_BACK; //backscattering
		else if (rnd < (components[0] + components[1] + components[2] + components[3])) return REFL_TRANS; //transparent pass
	}
	return REFL_ABSORB; //absorption
}

double QuadraticInterpolateX(const double& y,
                             const double& a, const double& b, const double& c,
                             const double& FA, const double& FB, const double& FC) {
    double amb = a - b;
    double amc = a - c;
    double bmc = b - c;
    double amb_amc = amb*amc;
    double amc_bmc = amc*bmc;
    double divisor = (2 * (-(FA / (amb_amc)) + FB / (amb_amc)+FB / (amc_bmc)-FC / (amc_bmc)));

    if (fabs(divisor) < 1e-30) {
        //Divisor is 0 when the slope is 1st order (a straight line) i.e. (FC-FB)/(c-b) == (FB-FA)/(b-a)
        if ((FB - FA) < 1e-30) {
            //FA==FB, shouldn't happen
            return a + rnd()*(b - a);
        }
        else {
            //Inverse linear interpolation
            return a + (y - FA) * (b - a)/ (FB - FA);
        }
    }
    else {
        //(reverse interpolate y on a 2nd order polynomial fitted on {a,FA},{b,FB},{c,FC} where FA<y<FB):
        //Root of Lagrangian polynomials, solved by Mathematica
        return (FA / (amb)-(a*FA) / (amb_amc)-(b*FA) / (amb_amc)-FB / (amb)+(a*FB) / (amb_amc)+(b*FB) / (amb_amc)+(a*FB) / (amc_bmc)+(b*	FB) / (amc_bmc)-(a*FC) / (amc_bmc)
                -(b*FC) / (amc_bmc)-sqrt(Sqr(-(FA / (amb)) + (a*FA) / (amb_amc)+(b*FA) / (amb_amc)+FB / (amb)-(a*FB) / (amb_amc)-(b*FB) / (amb_amc)-(a*FB) / (amc_bmc)-(b*FB)
                                                                                                                                                                       / (amc_bmc)+(a*FC) / (amc_bmc)+(b*FC) / (amc_bmc)) - 4 * (-(FA / (amb_amc)) + FB / (amb_amc)+FB / (amc_bmc)-FC / (amc_bmc))*(-FA + (a*FA) / (amb)-(a*b*FA)
                                                                                                                                                                                                                                                                                                                         / (amb_amc)-(a*FB) / (amb)+(a*b*FB) / (amb_amc)+(a*b*FB) / (amc_bmc)-(a*b*FC) / (amc_bmc)+y))) / divisor;
    }
}