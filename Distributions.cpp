#include <malloc.h>
#include <cstring>
#include <math.h>
#include "Distributions.h"
#include "Tools.h"
#include "Random.h"
#include "GLApp\GLTypes.h"
#include <math.h>

Distribution2D K_1_3_distribution=Generate_K_Distribution(1.0/3.0);
Distribution2D K_2_3_distribution=Generate_K_Distribution(2.0/3.0);
Distribution2D integral_N_photons=Generate_Integral(LOWER_LIMIT,UPPER_LIMIT,INTEGRAL_MODE_N_PHOTONS);
Distribution2D integral_SR_power=Generate_Integral(LOWER_LIMIT,UPPER_LIMIT,INTEGRAL_MODE_SR_POWER);
Distribution2D polarization_distribution=Generate_Polarization_Distribution(true,true);
Distribution2D g1h2_distribution=Generate_G1_H2_Distribution();

Distribution2D::Distribution2D(int N){
	if (!(N>0) && (N<10000000)) N=1; //don't create 0 size distributions
	valuesX=(double*)malloc(N*sizeof(double));
	valuesY=(double*)malloc(N*sizeof(double));
	size=N;
	average=0.0;
}

Distribution2D::Distribution2D(const Distribution2D &copy_src){ //copy constructor to avoid shallow copy
	valuesX=(double*)malloc(copy_src.size*sizeof(double));
	valuesY=(double*)malloc(copy_src.size*sizeof(double));
	memcpy(valuesX,copy_src.valuesX,copy_src.size*sizeof(double));
	memcpy(valuesY,copy_src.valuesY,copy_src.size*sizeof(double));
	size=copy_src.size;
	average=copy_src.average;
	average1=copy_src.average1;
}

Distribution2D::~Distribution2D(){
	SAFE_FREE(valuesX);
	SAFE_FREE(valuesY);
}

Distribution2D& Distribution2D::operator= (const Distribution2D &copy_src) {
	if (this != &copy_src) // protect against invalid self-assignment
	{
		//*this=Distribution2D_tiny(copy_src.size);
		valuesX=(double*)malloc(copy_src.size*sizeof(double));
		valuesY=(double*)malloc(copy_src.size*sizeof(double));
		memcpy(valuesX,copy_src.valuesX,copy_src.size*sizeof(double));
		memcpy(valuesY,copy_src.valuesY,copy_src.size*sizeof(double));
		size=copy_src.size;
	}
	// by convention, always return *this
	return *this;
}

double Distribution2D::InterpolateY(const double &x) {
	int inferior_index,superior_index;
	double slope, overshoot;

	for (superior_index=0;valuesX[superior_index]<x && superior_index<size;superior_index++);
	if (superior_index==size) superior_index--; //not found, x too large
	if (superior_index==0)    superior_index++; //not found, x too small
	inferior_index=superior_index-1;

	double diffX=valuesX[superior_index]-valuesX[inferior_index];
	double diffY=valuesY[superior_index]-valuesY[inferior_index];
	slope=diffY/diffX;
	overshoot=x-valuesX[inferior_index];

	return valuesY[inferior_index]+slope*overshoot;
}

double Distribution2D::InterpolateX(const double &y) {
	int inferior_index,superior_index;
	double slope, overshoot;

	for (superior_index=0;valuesY[superior_index]<y && superior_index<size;superior_index++);
	if (superior_index==size) return valuesX[size-1]; //not found, y too large
	if (superior_index==0)    return valuesX[0]; //not found, y too small
	inferior_index=superior_index-1;

	double diffX=valuesX[superior_index]-valuesX[inferior_index];
	double diffY=valuesY[superior_index]-valuesY[inferior_index];
	slope=diffY/diffX;
	if (slope==0.0) return valuesX[inferior_index];
	overshoot=y-valuesY[inferior_index];

	return valuesX[inferior_index]+overshoot/slope;
}

double Distribution2D::Interval_Mean(const double &x1,const double &x2) {
	//average of a cumulative distribution (differentiation included)
	
	int i1=integral_N_photons.findXindex(log(x1))+1;
	int i2=integral_N_photons.findXindex(log(x2))+1;
	if (i1==0) i1++;
	if (i1>=size-1) {
		i1=size-2;
		i2=size-1;
	}
	if (i2==0) i2=1;
	
	//double delta=log(b/a);
	double delta=log(x2/x1);
	int no_steps=i2-i1;
	if (i2>i1) delta=delta/(i2-i1);
	double temp=0.0;
	double tempp=0.0;
	double average=0.0;
	double average1=0.0;
	int i;
	for (i=0;i<=(i2-i1-1);i++){
		double xG1=log(x1)+i*delta;
		xG1=exp(xG1);
		double xG2=log(x1)+(i+1)*delta;
		xG2=exp(xG2);
		double delta1=xG2-xG1;
		double xG=log(x1)+(i+0.5)*delta;
		xG=exp(xG);
		double temp1=SYNRAD_FAST(xG)*delta1;
		average+=temp1*xG;
		average1+=temp1;

	}
	double temp3;
	if (average1>VERY_SMALL) temp3=average/average1;
	else temp3=(x1+x2)/2.0;

	return temp3;
}

int Distribution2D::findXindex(const double &x) {
	int superior_index;
	for (superior_index=0;valuesX[superior_index]<x && superior_index<size;superior_index++);
	return superior_index;
}

double g0ki(double x, double order, int kind) {
	/*ported from CALCF1.PAS:function g0ki(x,ord:realt1;kind:integer):realt1;

	{ Adapted from B. Diviacco, Sincrotrone Trieste, private comm. }
	{ Calculates the modified Bessel functions K1/3, and K2/3 for obtaining }
	{ the degree of polarization, and the integral of K5/3, used to obtain the }
	{ number of SR photons. In this case it's equivalent to SYNRAD_, but slower.}
	*/

	double h1,g0,r1,q1,q2,s1,s2,t1,xs1; //absolutely no idea what these variables are or how this function works :(
	h1=0.5;
	g0=0.0;
	r1=0.0;
	do {
		r1=r1+1.0;
		q1=exp(r1*h1);
		q2=exp(order*r1*h1);
		s1=(q1+1.0/q1)/2.0;
		s2=(q2+1.0/q2)/2.0;
		xs1=x*s1;
		if (kind==0) t1=exp(-xs1)*s2; // kind=0 ---> calculates the functions 
		else t1=exp(-xs1)*s2/s1;     // kind=1 ---> calculates the integral  
		g0=g0+t1;
	} while (t1>1.0E-6);
	return h1*(exp(-x)/2.0+g0);
}

double Gi(double x,int order) {
	//pascal code: Gi:=exp(ord*ln(x))*g0ki(x,5.0/3.0,1);
	return pow(x,order)*g0ki(x,5.0/3.0,1);
}

double H(double x, int order) {
	//pascal code: H:=exp(ord*ln(x))*Sqr(g0ki(x/2.0,2.0/3.0,0));
	return pow(x,order)*pow(g0ki(x/2.0,2.0/3.0,0),2);
}

Distribution2D Generate_K_Distribution(double order){
	Distribution2D result(NUMBER_OF_DISTRO_VALUES);
	double stepDelta=log(UPPER_LIMIT/LOWER_LIMIT)/NUMBER_OF_DISTRO_VALUES;
	for (int i=0;i<NUMBER_OF_DISTRO_VALUES;i++) {
		double x=LOWER_LIMIT*exp(i*stepDelta);
		result.valuesX[i]=x;
		result.valuesY[i]=log(Max(g0ki(x,order,0),VERY_SMALL));
		//note: values go from 0 to 99 (in original code: 1 to 100)
	}
	return result;
}

Distribution2D Generate_G1_H2_Distribution(){
	Distribution2D result(NUMBER_OF_DISTRO_VALUES);
	double stepDelta=log(UPPER_LIMIT/LOWER_LIMIT)/NUMBER_OF_DISTRO_VALUES;
	for (int i=0;i<NUMBER_OF_DISTRO_VALUES;i++) {
		double x=LOWER_LIMIT*exp(i*stepDelta);
		result.valuesX[i]=x;
		result.valuesY[i]=log(Max((Gi(x,1)/H(x,2)),VERY_SMALL));
		//note: values go from 0 to 99 (in original code: 1 to 100)
	}
	return result;
}

Distribution2D Generate_Integral(double x1,double x2,int mode){
	/* ported from CALCF1.PAS:function integral(x1,x2:realt1):realt1;

	{ Calculates the two real vectors integ[i,1], for the numbers of SR photons,  }
	{ integ[i,2] for the SR power, 						      }
	{ and integ[i,3], for the corresponding energies, for generating SR photons   }
	{ according to the real distribution, see fig. 2 ref. G.K. Green for instance }
	*/
	double delta,delta1,sum,temp1,temp2,sum_power,temp2_power,xG,xG1,xG2;
	
	int i;
	Distribution2D result(NUMBER_OF_INTEGR_VALUES);

	delta=log(x2/x1)/NUMBER_OF_INTEGR_VALUES;
	sum=0.0;
	sum_power=0.0;
	result.average=0.0;
	result.average1=0.0;
	for (i=0;i<NUMBER_OF_INTEGR_VALUES;i++) {
		xG1=log(x1)+i*delta;
		xG1=exp(xG1); //left discrete point
		xG2=log(x1)+(i+1)*delta; 
		xG2=exp(xG2); //right discrete point
		delta1=xG2-xG1; //distance between them
		xG=log(x1)+(i+0.5)*delta;
		xG=exp(xG); //middle point
		
		temp1=SYNRAD_FAST(xG)*delta1;
		temp2=(SYNRAD_FAST(xG1)+SYNRAD_FAST(xG2))*delta1;
		
		temp2_power=temp2*xG;
		sum+=temp2;
		sum_power+=temp2_power;
		result.valuesX[i]=log(xG);
		if (mode==INTEGRAL_MODE_N_PHOTONS) result.valuesY[i]=log(sum);
		else if (mode==INTEGRAL_MODE_SR_POWER) result.valuesY[i]=log(sum_power);
		//values are filled from 0 to 249!
		result.average+=temp1*xG;
		result.average1+=temp1;
	}
	//i=NUMBER_OF_INTEGR_VALUES;
	//integral:=sum*7.7085761E16; //No return value in this implementation
	//{ integral(1e-10,100)=8.084227E17 ph/s/mA/GeV }
	if (result.average1>0.0) result.average=result.average/result.average1;
	else result.average=x1;
	return result;
}

Distribution2D Generate_Polarization_Distribution(bool calculate_parallel_polarization, bool calculate_orthogonal_polarization) {
	/*ported from CALCF1.PAS - procedure findpolarization;

	{ Calculates the polarization percentages vs. normalized photon energy }
	{ in the interval 1.0E-10 < (energy/critical_en.) < 100, as partially }
	{ reported in fig. 14, ref. G.K. Green.}
	{ An accuracy of better than 0.5% is usually obtained.}
	*/

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

double calc_polarization_percentage(double energy,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization){
	/*ported from CALCF1.PAS: function find_polarization(x:realt1):realt1;
	Besself[6] removed as distribution is never used later

	{ Calculates the degree of polarization (linear and orthogonal) for a }
	{ given SR photon energy }
	{ Ref. G.K. Green, page 1-12, and fig. 15 }
	*/
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

double find_psi(double lambda_ratios,double gamma_square,double f_times_g1h2,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization){
	/*ported from FINDPRCS.PAS: function find_psi(x,maxangle:realt1):realt1;

	{ Derived from procedure find_polarization ... }
	*/

	//RETURNS ANGLE*GAMMA

	double F_parallel,F_orthogonal;
	double delta,csi,X1,factor_psi,seed;
	Distribution2D local_polarization_integral(NUMBER_OF_DISTRO_VALUES+1);
	F_parallel=0.0;
	F_orthogonal=VERY_SMALL;
	delta=4.344*f_times_g1h2/NUMBER_OF_DISTRO_VALUES;
	//index=1;
	/*
	{ delta gives the angular range (for traj.gamma*Psi) where the vertical }
	{ angle of emission, Psi, should be contained (actually, traj.gamma*Psi...) }
	*/

	for (int i=1;i<=NUMBER_OF_DISTRO_VALUES;i++) { 
		X1=(i-1)*delta;
		csi=1.0+Sqr(X1);
		csi=lambda_ratios/2.0*pow(csi,3.0/2.0);
		// Nem marad le egy LN???

		double one_plus_square_X1=1.0+Sqr(X1);

		if (calculate_parallel_polarization) 
			F_parallel=Sqr(one_plus_square_X1*exp(K_2_3_distribution.InterpolateY(csi)));
		else F_parallel=0.0;

		if (calculate_orthogonal_polarization) 
			F_orthogonal=one_plus_square_X1*Sqr(X1*exp(K_1_3_distribution.InterpolateY(csi)));
		else F_orthogonal=0.0;

		local_polarization_integral.valuesX[i-1]=X1;
		if (i>1) local_polarization_integral.valuesY[i]=local_polarization_integral.valuesY[i-1]+F_parallel+F_orthogonal; //integrate full polarization
		else local_polarization_integral.valuesY[i]=F_parallel+F_orthogonal;
	}
	local_polarization_integral.valuesY[0]=0.0;
	factor_psi=1.0;
	seed=rnd()*factor_psi*local_polarization_integral.valuesY[NUMBER_OF_DISTRO_VALUES-1];
	return local_polarization_integral.InterpolateX(seed);
}

double find_chi(double psi,double gamma_square,double f_times_g1h2,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization) {
	/*
	ported from FINDPRCS.PAS - function find_chi(psi:realt1):realt1;
	*/

	//RETURNS ANGLE

	double delta;
	//double j1;
	double F_parallel,F_orthogonal;
	double chi,factor;
	int i;
	Distribution2D local_polarization_integral(NUMBER_OF_DISTRO_VALUES);

	F_parallel=0.0;
	F_orthogonal=VERY_SMALL;
	delta=4.344*f_times_g1h2/NUMBER_OF_DISTRO_VALUES/gamma_square;
	/*
	{ delta:=pi/100.0/gamma_square;}
	{ delta gives the angular range (for traj.gamma*Psi) where the vertical }
	{ angle of emission, Psi, should be contained (actually, traj.gamma*Psi...) }
	*/

	double cos_psi=cos(psi);
	double sin_psi=sin(psi);
	double sin_chi,cos_chi;
	for (i=1;i<=NUMBER_OF_DISTRO_VALUES;i++){
		chi=(i-0.5)*delta;
		sin_chi=sin(chi);
		cos_chi=cos(chi);
		/*
		{1  teta:=acos(cospsi*coschi);}
		{2  teta:=asin(sqrt(Sqr(sinpsi)+Sqr(cospsi*sinchi)));}
		{ Version 1 and 2 of teta=teta(psi,chi) generate round-off errors }
		*/
		double theta=atan(sqrt(Sqr(sin_psi)+Sqr(cos_psi*sin_chi)))/(cos_psi*cos_chi);
		double gamma_times_theta=gamma_square*theta;
		double gamma_times_theta_square=Sqr(gamma_times_theta);
		double sin_theta=sin(theta);
		double sin_fi=sin_psi/sin_theta;
		double cos_fi=cos_psi*sin_chi/sin_theta;
		factor=exp(-6.0*log(1.0+gamma_times_theta_square));
		//{      The 6.0 in the exp agrees with Coisson's results ... }
		if (calculate_parallel_polarization)
			F_parallel=Sqr(1.0-gamma_times_theta_square+2.0*Sqr(gamma_times_theta*sin_fi))*factor;
		else F_parallel=0.0;
		
		if (calculate_orthogonal_polarization)
			F_orthogonal=Sqr(2.0*gamma_times_theta_square*sin_fi*cos_fi)*factor;
		else F_orthogonal=0.0;
		
		//here i goes from 1 to 100!
		local_polarization_integral.valuesX[i-1]=chi;
		if (i>1) local_polarization_integral.valuesY[i-1]=local_polarization_integral.valuesY[i-2]+F_parallel+F_orthogonal;
		else local_polarization_integral.valuesY[i-1]=F_parallel+F_orthogonal;
	}

	double seed=rnd()*local_polarization_integral.valuesY[NUMBER_OF_DISTRO_VALUES-1];
	return local_polarization_integral.InterpolateX(seed);
}

double SYNGEN1(double x_min,double x_max,int mode) {
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
	if (mode==SYNGEN_MODE_FLUXWISE) {
		double flux_min=exp(integral_N_photons.InterpolateY(log(x_min)));
		double flux_max=exp(integral_N_photons.InterpolateY(log(x_max)));
		double delta_flux=flux_max-flux_min;
		double generated_flux=flux_min+rnd()*delta_flux; //uniform distribution between flux_min and flux_max
		generated_energy=integral_N_photons.InterpolateX(log(generated_flux));
	}
	else if (mode==SYNGEN_MODE_POWERWISE) {
		double power_min=exp(integral_SR_power.InterpolateY(log(x_min)));
		double power_max=exp(integral_SR_power.InterpolateY(log(x_max)));
		double delta_power=power_max-power_min;
		double generated_power=power_min+rnd()*delta_power; //uniform distribution between flux_min and flux_max
		generated_energy=integral_SR_power.InterpolateX(log(generated_power));
	}
	return exp(generated_energy);
}

double SYNRAD_FAST(const double &x) {
/*
{ Adapted from H.H. Umstaetter, CERN/PS/SM/81-13. }
{ Works as g0ki(x,5/3,1), but about 2.5x faster }
*/
double Y,Z,A,B,P,Q;


  if (x<6.0) { 
   
      Z=pow(x,2)/16.0-2.0;
      A=         +0.0000000001;
      B=   Z*A   +0.0000000023;
      A=   Z*B-A +0.0000000813;
      B=   Z*A-B +0.0000024575;
      A=   Z*B-A +0.0000618126;
      B=   Z*A-B +0.0012706638;
      A=   Z*B-A +0.0209121680;
      B=   Z*A-B +0.2688034606;
      A=   Z*B-A +2.6190218379;
      B=   Z*A-B +18.6525089687;
      A=   Z*B-A +92.9523266592;
      B=   Z*A-B +308.1591941313;
      A=   Z*B-A +644.8697965824;
      P=0.5*Z*A-B+414.5654364883;
      A=         +0.0000000012;
      B=   Z*A   +0.0000000391;
      A=   Z*B-A +0.0000011060;
      B=   Z*A-B +0.0000258145;
      A=   Z*B-A +0.0004876869;
      B=   Z*A-B +0.0072845620;
      A=   Z*B-A +0.0835793546;
      B=   Z*A-B +0.7103136120;
      A=   Z*B-A +4.2678026127;
      B=   Z*A-B +17.0554078580;
      A=   Z*B-A +41.8390348678;
      Q=0.5*Z*A-B+28.4178737436;
      Y=exp((2.0/3.0)*log(x));
      return (P/Y-Q*Y-1.0)*1.8137993642;
  }
  else if (x<80.0)  
  {
      Z=20.0/x-2.0;
      A=         +0.0000000001;
      B=   Z*A   -0.0000000004;
      A=   Z*B-A +0.0000000020;
      B=   Z*A-B -0.0000000110;
      A=   Z*B-A +0.0000000642;
      B=   Z*A-B -0.0000004076;
      A=   Z*B-A +0.0000028754;
      B=   Z*A-B -0.0000232125;
      A=   Z*B-A +0.0002250532;
      B=   Z*A-B -0.0028763680;
      A=   Z*B-A +0.0623959136;
      P=0.5*Z*A-B+1.0655239080;
      return P*sqrt(1.57079632679/x)/exp(x);
  }
  else return 0.0;
}

void Material::LoadCSV(FileReader *file){
	file->SeekForChar(","); //find first comma (skip A1 cell)
	file->ReadKeyword(",");
	double val;
	
	do { //store values and read subsequent angles
		val=file->ReadDouble();
		angleVals.push_back(val);
		if (!file->IsEol()) file->ReadKeyword(",");
	} while (!file->IsEol());
	
	do {
		std::vector<double> currentRow;
		val=file->ReadDouble();
		energyVals.push_back(val);
		for (int i=0;i<(int)angleVals.size();i++) {
			file->ReadKeyword(",");
			val=file->ReadDouble();
			currentRow.push_back(val);
		}
		reflVals.push_back(currentRow);
	} while (!file->IsEof());
}

double Material::Interpolate(const double &energy,const double &angle) {
	int angleLowerIndex,energyLowerIndex;
	for (angleLowerIndex=0;angleLowerIndex<((int)angleVals.size()-1)&&angle>angleVals[angleLowerIndex+1];angleLowerIndex++);
	for (energyLowerIndex=0;energyLowerIndex<((int)energyVals.size()-1)&&energy>energyVals[energyLowerIndex+1];energyLowerIndex++);
	
	if (angleLowerIndex==((int)angleVals.size()-1)) angleLowerIndex--; //if not in table
	if (energyLowerIndex==((int)energyVals.size()-1)) energyLowerIndex--; //if not in table

	double angleOvershoot=log(angle)-log(angleVals[angleLowerIndex]);
	double angleDelta=log(angleVals[angleLowerIndex+1])-log(angleVals[angleLowerIndex]);

	double energyOvershoot=log(energy)-log(energyVals[energyLowerIndex]);
	double energyDelta=log(energyVals[energyLowerIndex+1])-log(energyVals[energyLowerIndex]);

	double interpolatedReflForLowerAngle=reflVals[energyLowerIndex][angleLowerIndex]+
		energyOvershoot/energyDelta*(reflVals[energyLowerIndex+1][angleLowerIndex]-reflVals[energyLowerIndex][angleLowerIndex]);

	double interpolatedReflForHigherAngle=reflVals[energyLowerIndex][angleLowerIndex+1]+
		energyOvershoot/energyDelta*(reflVals[energyLowerIndex+1][angleLowerIndex+1]-reflVals[energyLowerIndex][angleLowerIndex+1]);

	double interpRefl=interpolatedReflForLowerAngle+angleOvershoot/angleDelta*(interpolatedReflForHigherAngle-interpolatedReflForLowerAngle);
	SATURATE(interpRefl, 0.0, 100.0);
	return interpRefl;
}