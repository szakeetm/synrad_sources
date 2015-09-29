#include <malloc.h>
#include <cstring>
#include <math.h>
#include "Distributions.h"
#include "Tools.h"
#include "Random.h"
#include "GLApp\GLTypes.h"
#include <math.h>

Distribution2D K_1_3_distribution = Generate_K_Distribution(1.0 / 3.0);
Distribution2D K_2_3_distribution = Generate_K_Distribution(2.0 / 3.0);
Distribution2D integral_N_photons = Generate_Integral(LOWER_LIMIT, UPPER_LIMIT, INTEGRAL_MODE_N_PHOTONS);
Distribution2D integral_SR_power = Generate_Integral(LOWER_LIMIT, UPPER_LIMIT, INTEGRAL_MODE_SR_POWER);
//Distribution2D polarization_distribution=Generate_Polarization_Distribution(true,true);
//Distribution2D g1h2_distribution=Generate_G1_H2_Distribution();

Distribution2D::Distribution2D(int N){
	if (!(N > 0) && (N < 10000000)) N = 1; //don't create 0 size distributions
	valuesX = (double*)malloc(N*sizeof(double));
	valuesY = (double*)malloc(N*sizeof(double));
	size = N;
	//sum_energy=sum_photons= 0.0;
}

Distribution2D::Distribution2D(const Distribution2D &copy_src){ //copy constructor to avoid shallow copy
	valuesX = (double*)malloc(copy_src.size*sizeof(double));
	valuesY = (double*)malloc(copy_src.size*sizeof(double));
	memcpy(valuesX, copy_src.valuesX, copy_src.size*sizeof(double));
	memcpy(valuesY, copy_src.valuesY, copy_src.size*sizeof(double));
	size = copy_src.size;
	//sum_energy = copy_src.sum_energy;
	//sum_photons = copy_src.sum_photons;
}

Distribution2D::~Distribution2D(){
	SAFE_FREE(valuesX);
	SAFE_FREE(valuesY);
}

Distribution2D& Distribution2D::operator= (const Distribution2D &copy_src) {
	if (this != &copy_src) // protect against invalid self-assignment
	{
		//*this=Distribution2D_tiny(copy_src.size);
		valuesX = (double*)malloc(copy_src.size*sizeof(double));
		valuesY = (double*)malloc(copy_src.size*sizeof(double));
		memcpy(valuesX, copy_src.valuesX, copy_src.size*sizeof(double));
		memcpy(valuesY, copy_src.valuesY, copy_src.size*sizeof(double));
		size = copy_src.size;
	}
	// by convention, always return *this
	return *this;
}

double Distribution2D::InterpolateY(const double &x) {
	int inferior_index, superior_index;
	double slope, overshoot;

	//for (superior_index = 0; valuesX[superior_index] < x && superior_index < size; superior_index++); //To replace by binary search
	inferior_index = binary_search(x, valuesX, size);
	if (inferior_index == size - 1) inferior_index--; //not found, x too large
	if (inferior_index == -1)    inferior_index++; //not found, x too small
	superior_index = inferior_index + 1;

	double diffX = valuesX[superior_index] - valuesX[inferior_index];
	double diffY = valuesY[superior_index] - valuesY[inferior_index];
	slope = diffY / diffX;
	overshoot = x - valuesX[inferior_index];

	return valuesY[inferior_index] + slope*overshoot;
}

double Distribution2D::InterpolateX(const double &y) {
	int inferior_index, superior_index;
	double slope, overshoot;

	//for (superior_index = 0; valuesY[superior_index] < y && superior_index < size; superior_index++); //To replace by binary search
	inferior_index = binary_search(y, valuesY, size);
	if (inferior_index == size-1) return valuesX[size - 1]; //not found, y too large
	if (inferior_index == -1)    return valuesX[0]; //not found, y too small
	superior_index = inferior_index + 1;

	double diffX = valuesX[superior_index] - valuesX[inferior_index];
	double diffY = valuesY[superior_index] - valuesY[inferior_index];
	slope = diffY / diffX;
	if (slope == 0.0) return valuesX[inferior_index];
	overshoot = y - valuesY[inferior_index];

	return valuesX[inferior_index] + overshoot / slope;
}

/*int Distribution2D::findXindex(const double &x) {
	int superior_index;
	for (superior_index = 0; valuesX[superior_index] < x && superior_index < size; superior_index++); //replace by binary search
	return superior_index;
	}*/

double g0ki(double x, double order, int kind) {
	/*ported from CALCF1.PAS:function g0ki(x,ord:realt1;kind:integer):realt1;

	{ Adapted from B. Diviacco, Sincrotrone Trieste, private comm. }
	{ Calculates the modified Bessel functions K1/3, and K2/3 for obtaining }
	{ the degree of polarization, and the integral of K5/3, used to obtain the }
	{ number of SR photons. In this case it's equivalent to SYNRAD_, but slower.}
	*/

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
}
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

Distribution2D Generate_Integral(double log_min, double log_max, int mode){
	/* ported from CALCF1.PAS:function integral(x1,x2:realt1):realt1;

	{ Calculates the two real vectors integ[i,1], for the numbers of SR photons,  }
	{ integ[i,2] for the SR power, 						      }
	{ and integ[i,3], for the corresponding energies, for generating SR photons   }
	{ according to the real distribution, see fig. 2 ref. G.K. Green for instance }
	*/

	//log_min: log(ratio of minimal energy / critical energy)
	//log_max: log(ratio of maximal energy / critical energy)

	double delta, exp_delta, interval_dN, interval_dN2, sum_photons, sum_power, temp2_power, x_middle, x_lower, x_higher;

	int i;
	Distribution2D result(NUMBER_OF_INTEGR_VALUES);

	delta = (log_max - log_min) / NUMBER_OF_INTEGR_VALUES;
	sum_photons = sum_power = /*result.sum_photons = result.sum_energy =*/ 0.0;

	for (i = 0; i < NUMBER_OF_INTEGR_VALUES; i++) {
		x_lower = exp(log_min + i*delta); //lower energy
		x_middle = exp(log_min + (i + 0.5)*delta); //middle energy
		x_higher = exp(log_min + (i + 1.0)*delta); //higher energy
		exp_delta = x_higher - x_lower; //actual energy range for the next index

		interval_dN = SYNRAD_FAST(x_middle)*exp_delta; //number of photons for the actual energy interval
		interval_dN2 = (SYNRAD_FAST(x_lower) + SYNRAD_FAST(x_higher)) / 2.0*exp_delta; //number of photons for the actual interval, different averaging
		sum_photons += interval_dN2; //flux increment
		sum_power += interval_dN2*x_middle; //2*number of photons * average energy: 2*energy of the interval
		result.valuesX[i] = log(x_middle);
		if (mode == INTEGRAL_MODE_N_PHOTONS) result.valuesY[i] = sum_photons; //used to be log(sum_flux)
		else if (mode == INTEGRAL_MODE_SR_POWER) result.valuesY[i] = sum_power; //used to be log(sum_power)

		//values are filled from 0 to 249!
		//result.sum_photons += interval_dN2;
		//result.sum_energy += interval_dN2*x_middle;
	}
	//i=NUMBER_OF_INTEGR_VALUES;
	//integral:=sum*7.7085761E16; //No return value in this implementation
	//{ integral(1e-10,100)=8.084227E17 ph/s/mA/GeV }
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

/*double find_psi(double lambda_ratios,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization){
//vertical

//ported from FINDPRCS.PAS: function find_psi(x,maxangle:realt1):realt1;

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

double find_psi(double lambda_ratios, std::vector<std::vector<double>> &psi_distro){
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
	int size = psi_distro[0].size();
	int imax = size;
	double interpolated_CDF_lower, interpolated_CDF_higher;
	// continue searching while [imin,imax] is not empty
	while (imin <= imax)
	{
		// calculate the midpoint for roughly equal partition
		int imid = (imin + imax) / 2;
		interpolated_CDF_lower = WEIGH(psi_distro[lambda_lower_index][imid],psi_distro[lambda_lower_index + 1][imid],lambda_overshoot);
		interpolated_CDF_higher = WEIGH(psi_distro[lambda_lower_index][imid+1], psi_distro[lambda_lower_index + 1][imid+1], lambda_overshoot);
		if (imid == size - 1 || (interpolated_CDF_lower < lookup && lookup < interpolated_CDF_higher)) {
			// key found at index imid
			foundAngle = imid+1; //will be lowered by 1
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
	return psi;
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
//{ delta gives the angular range (for traj.gamma*Psi) where the vertical }
//{ angle of emission, Psi, should be contained (actually, traj.gamma*Psi...) }


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
double find_chi(double psi, double gamma, std::vector<std::vector<double>> &chi_distro) {

	double psi_index, chi_lower, chi_higher, chi_higher2, chi;
	double psi_relative = log10(abs(psi)*(gamma / 10000.0)); //distributions are digitized for gamma=10000, and sampled logarithmically
	if (psi_relative < -7.0) {
		psi_index = 0; //use lowest angle
	}
	else {
		psi_index = (psi_relative + 7.0) / 0.04; //sampled from -7 to 0 with delta=0.04
	}

	int psi_lower_index = (int)(psi_index); //digitized for -2PI/10 .. +2PI/10 with delta=0.0025
	double psi_overshoot = psi_index - (double)psi_lower_index;

	double lookup = rnd();

	int foundAngle = 0;
	/*interpolated_CDF = 0.0;
	do { //to replace by binary search
		previous_interpolated_CDF = interpolated_CDF;
		foundAngle++;
		interpolated_CDF = WEIGH(chi_distro[foundAngle][psi_lower_index], chi_distro[foundAngle][psi_lower_index + 1], psi_overshoot);
	} while (interpolated_CDF < lookup);*/

	//Binary search
	int imin = 0;
	int size = chi_distro.size();
	int imax = size;
	double interpolated_CDF_lower, interpolated_CDF_higher;
	// continue searching while [imin,imax] is not empty
	while (imin <= imax)
	{
		// calculate the midpoint for roughly equal partition
		int imid = (imin + imax) / 2;
		interpolated_CDF_lower = WEIGH(chi_distro[imid][psi_lower_index], chi_distro[imid][psi_lower_index + 1], psi_overshoot);
		interpolated_CDF_higher = WEIGH(chi_distro[imid + 1][psi_lower_index], chi_distro[imid + 1][psi_lower_index + 1], psi_overshoot);
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


	//TO DO: Treat not found errors
	int chi_lower_index = foundAngle - 1;

	double next_interpolated_CDF;
	if (foundAngle == chi_distro.size() - 1) next_interpolated_CDF = 1.0;
	else next_interpolated_CDF = WEIGH(chi_distro[chi_lower_index + 2][psi_lower_index], chi_distro[chi_lower_index + 2][psi_lower_index + 1], psi_overshoot);
	//double chi_overshoot = (log10(lookup)-log10(previous_interpolated_CDF))/ (log10(interpolated_CDF) - log10(previous_interpolated_CDF));
	//double chi_overshoot = (lookup - previous_interpolated_CDF) / (interpolated_CDF - previous_interpolated_CDF); //linear scale

	if (chi_lower_index == 0) {
		chi_lower = 0;
		chi_higher = 1.0964782E-7;
		chi = WEIGH(chi_lower, chi_higher, rnd()) / (gamma / 10000.0);
		//chi = 0;
	}
	else {

		double a = pow(10, -7.0 + ((double)chi_lower_index + 0.0)*0.04) / (gamma / 10000.0);
		double b = a*1.0964782; /* pow(10, -7.0 + ((double)chi_lower_index + 1.0)*0.04) / (gamma / 10000.0);*/ //1.09=10^0.04
		double c = b*1.0964782; /*pow(10, -7.0 + ((double)chi_lower_index + 2.0)*0.04) / (gamma / 10000.0);*/
		double FA = interpolated_CDF_lower;
		double FB = interpolated_CDF_higher;
		double FC = next_interpolated_CDF;

		//inverse 2nd degree polynomial interpolation
		double V = lookup;
		double amb = a - b;
		double amc = a - c;
		double bmc = b - c;
		double amb_amc = amb*amc;
		double amc_bmc = amc*bmc;
		chi = (FA / (amb)-(a*FA) / (amb_amc)-(b*FA) / (amb_amc)-FB / (amb)+(a*FB) / (amb_amc)+(b*FB) / (amb_amc)+(a*FB) / (amc_bmc)+(b*	FB) / (amc_bmc)-(a*FC) / (amc_bmc)
			-(b*FC) / (amc_bmc)-sqrt(Sqr(-(FA / (amb)) + (a*FA) / (amb_amc)+(b*FA) / (amb_amc)+FB / (amb)-(a*FB) / (amb_amc)-(b*FB) / (amb_amc)-(a*FB) / (amc_bmc)-(b*FB)
			/ (amc_bmc)+(a*FC) / (amc_bmc)+(b*FC) / (amc_bmc)) - 4 * (-(FA / (amb_amc)) + FB / (amb_amc)+FB / (amc_bmc)-FC / (amc_bmc))*(-FA + (a*FA) / (amb)-(a*b*FA)
			/ (amb_amc)-(a*FB) / (amb)+(a*b*FB) / (amb_amc)+(a*b*FB) / (amc_bmc)-(a*b*FC) / (amc_bmc)+V))) / (2 * (-(FA / (amb_amc)) + FB / (amb_amc)+FB / (amc_bmc)-FC / (amc_bmc)));

		//inverse linear interpolation
		//chi = WEIGH(a, b,rnd());

	}
	return chi;
}

double SYNGEN1(double log_x_min, double log_x_max, int mode) {
	/*
	{ Generates a random normalized SR photon energy in (xmin,xmax) from a        }
	{ cumulative distribution given by the Vector array integ[i,1] calculated by  }
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

	double generated_energy;
	if (mode == SYNGEN_MODE_FLUXWISE) {
		double flux_min = /*exp(*/integral_N_photons.InterpolateY(/*log(*/log_x_min/*)*/)/*)*/;
		double flux_max = /*exp(*/integral_N_photons.InterpolateY(/*log(*/log_x_max/*)*/)/*)*/;
		double generated_flux = WEIGH(flux_min, flux_max, rnd());//uniform distribution between flux_min and flux_max
		generated_energy = exp(integral_N_photons.InterpolateX(/*log(*/generated_flux/*)*/));
	}
	else if (mode == SYNGEN_MODE_POWERWISE) {
		double power_min = /*exp(*/integral_SR_power.InterpolateY(/*log(*/log_x_min/*)*/)/*)*/;
		double power_max = /*exp(*/integral_SR_power.InterpolateY(/*log(*/log_x_max/*)*/)/*)*/;
		double generated_power = WEIGH(power_min, power_max, rnd()); //uniform distribution between flux_min and flux_max
		generated_energy = exp(integral_SR_power.InterpolateX(/*log(*/generated_power/*)*/));
	}
	return generated_energy;
}

double SYNRAD_FAST(const double &x) {
	/*
	{ Adapted from H.H. Umstaetter, CERN/PS/SM/81-13. }
	{ Works as g0ki(x,5/3,1), but about 2.5x faster }
	*/
	double Y, Z, A, B, P, Q;


	if (x < 6.0) {

		Z = pow(x, 2) / 16.0 - 2.0;
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

void Material::LoadCSV(FileReader *file){
	file->SeekForChar(","); //find first comma (skip A1 cell)
	file->ReadKeyword(",");
	double val;

	do { //store values and read subsequent angles
		val = file->ReadDouble();
		angleVals.push_back(val);
		if (!file->IsEol()) file->ReadKeyword(",");
	} while (!file->IsEol());

	do {
		std::vector<double> currentRow;
		val = file->ReadDouble();
		energyVals.push_back(val);
		for (int i = 0; i < (int)angleVals.size(); i++) {
			file->ReadKeyword(",");
			val = file->ReadDouble();
			currentRow.push_back(val);
		}
		reflVals.push_back(currentRow);
	} while (!file->IsEof());
}

double Material::Interpolate(const double &energy, const double &angle) {
	int angleLowerIndex, energyLowerIndex;
	//for (angleLowerIndex = 0; angleLowerIndex<((int)angleVals.size() - 1) && angle>angleVals[angleLowerIndex + 1]; angleLowerIndex++); //replace by binary search
	//for (energyLowerIndex = 0; energyLowerIndex<((int)energyVals.size() - 1) && energy>energyVals[energyLowerIndex + 1]; energyLowerIndex++); //replace by binary search

	angleLowerIndex = binary_search(angle, angleVals, angleVals.size());
	energyLowerIndex = binary_search(energy, energyVals, energyVals.size());

	if (angleLowerIndex == ((int)angleVals.size() - 1)) angleLowerIndex--; //if not in table
	if (energyLowerIndex == ((int)energyVals.size() - 1)) energyLowerIndex--; //if not in table

	double angleOvershoot = log(angle) - log(angleVals[angleLowerIndex]);
	double angleDelta = log(angleVals[angleLowerIndex + 1]) - log(angleVals[angleLowerIndex]);

	double energyOvershoot = log(energy) - log(energyVals[energyLowerIndex]);
	double energyDelta = log(energyVals[energyLowerIndex + 1]) - log(energyVals[energyLowerIndex]);

	double interpolatedReflForLowerAngle = WEIGH(reflVals[energyLowerIndex][angleLowerIndex], reflVals[energyLowerIndex + 1][angleLowerIndex], energyOvershoot / energyDelta);
	double interpolatedReflForHigherAngle = WEIGH(reflVals[energyLowerIndex][angleLowerIndex + 1], reflVals[energyLowerIndex + 1][angleLowerIndex + 1], energyOvershoot / energyDelta);

	double interpRefl = WEIGH(interpolatedReflForLowerAngle, interpolatedReflForHigherAngle, angleOvershoot / angleDelta);
	SATURATE(interpRefl, 0.0, 100.0);
	return interpRefl;
}

template <typename T> int binary_search(double key, T A, int size)
{
	int imin = 0;
	int imax = size - 1;
	// continue searching while [imin,imax] is not empty
	while (imin <= imax)
	{
		// calculate the midpoint for roughly equal partition
		int imid = (imin+imax)/2;
		if (imid==size-1 || (A[imid] < key && key < A[imid+1]))
			// key found at index imid
			return imid;
		// determine which subarray to search
		else if (A[imid] < key)
			// change min index to search upper subarray
			imin = imid + 1;
		else
			// change max index to search lower subarray
			imax = imid - 1;
	}
	// key was not found
	return -1;
}