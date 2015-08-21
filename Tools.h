#ifndef _TOOLS_
#define _TOOLS_

#include "Distributions.h"
typedef  unsigned __int64 LONGINT;

#define PI 3.14159265358979323846
#define proton_charge 1.602E-19 //in Coulombs
#define proton_mass 1.67E-27 //in kg

#define WEIGH(a,b,weigh) a+(b-a)*weigh

class Vector {
public:
	double x;
	double y;
	double z;
	double Norme();
	Vector();
	Vector(const Vector &src);
	Vector& operator=(const Vector &src);
	Vector(const double &x,const double &y,const double &z);
	//~Vector();
	Vector Rotate(Vector Axis_dir,double alpha);
	Vector Normalize();
};

class Trajectory_Point {
public:
	Vector position,direction,X_local,Y_local,Z_local,B;
	Vector rho;
	double critical_energy, emittance_X, emittance_Y, beta_X, beta_Y, eta, eta_prime, energy_spread;
	double sigma_x, sigma_y, sigma_x_prime, sigma_y_prime;
	
	double Critical_Energy(const double &gamma);
	double dAlpha(const double &dL);
};

class GenPhoton {
public:
	double natural_divx, natural_divy, offset_x, offset_y, offset_divx, offset_divy;
	double radius, critical_energy, B_ort, B_par, energy;
	double g1h2, B_factor, B_factor_power, SR_flux, SR_power;
	Vector start_pos, start_dir,B;
};

class Histogram {
private:
	int number_of_bins;
	double delta,min,max;
	double *counts;
public:
	Histogram(double min,double max,int number_of_bins,bool logscale);
	~Histogram();
	double GetCount(int index);
	double GetFrequency(int index);
	//double GetNormalized(int index);
	double GetX(int index);
	void Add(const double &x,const double &dY,const double &bandwidth=-1);
	bool logarithmic;
	double /*max_count,*/total_count;
	void ResetCounts();
};

double DotProduct(const Vector &a,const Vector &b);
Vector Crossproduct(const Vector &v1,const Vector &v2);
Vector ScalarMult(const Vector &v,const double &r);
Vector RandomPerpendicularVector(const Vector &v,const double &length);
Vector Add(const Vector &v1,const Vector &v2);
Vector Subtract(const Vector &v1,const Vector &v2);
Vector Weigh(const Vector &v1,const Vector &v2,const double &weigh);
double Min(const double &a,const double &b);
double Max(const double &a,const double &b);
int Min(const int &a,const int &b);
int Max(const int &a,const int &b);
double Sqr(const double &x);
double Gaussian(const double &sigma);

#endif
