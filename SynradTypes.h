#pragma once

#include "Vector.h"

class Trajectory_Point {
public:
	Vector3d position,direction,X_local,Y_local,Z_local,B;
	Vector3d rho;
	double critical_energy, emittance_X, emittance_Y, beta_X, beta_Y, eta, eta_prime, alpha_X, alpha_Y;
	double sigma_x, sigma_y, sigma_x_prime, sigma_y_prime, theta_X, theta_Y, gamma_X, gamma_Y;
	double a_x, b_x, a_y, b_y;
	
	double Critical_Energy(const double &gamma);
	double dAlpha(const double &dL);
};

class GenPhoton {
public:
	double natural_divx, natural_divy, offset_x, offset_y, offset_divx, offset_divy;
	double radius, critical_energy, B_ort, B_par, energy;
	double/* g1h2,*/ B_factor, B_factor_power, SR_flux, SR_power;
	Vector3d start_pos, start_dir,B;
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