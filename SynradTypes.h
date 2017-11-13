#pragma once

#include "Vector.h"
#include "Buffer_shared.h"
#include "GLApp/GLTypes.h"

// Hit type
#define HIT_DES   1
#define HIT_ABS   2
#define HIT_REF   3
#define HIT_TRANS 4
#define HIT_TELEPORT 5
#define HIT_LAST 6

//Reflection type
#define REFL_ABSORB 0
#define REFL_FORWARD 1
#define REFL_DIFFUSE 2
#define REFL_BACK 3
#define REFL_TRANS 4

// Reflection type
#define REF_DIFFUSE 0   // Diffuse (cosine law)
#define REF_MIRROR  1   // Mirror
#define REF_MATERIAL 10 //Real rough surface reflection

// Profile type

#define REC_NONE       0  // No recording
#define REC_PRESSUREU  1  // Pressure profile (U direction)
#define REC_PRESSUREV  2  // Pressure profile (V direction)
#define REC_ANGULAR    3  // Angular profile

// Density/Hit field stuff
#define HITMAX_INT64 18446744073709551615
#define HITMAX_DOUBLE 1E308

typedef struct {
	int      generation_mode; //fluxwise or powerwise
	bool	 lowFluxMode;
	double	 lowFluxCutoff;
} SHMODE;

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
	double radius, critical_energy, B_ort, B_par, energy, polarization;
	double/* g1h2,*/ B_factor, B_factor_power, SR_flux, SR_power;
	Vector3d start_pos, start_dir,B;
};

#define SPECTRUM_SIZE (size_t)100 //number of histogram bins
class Histogram {
private:
	int number_of_bins;
	double delta,min,max;
	double *counts;
public:
	Histogram(double min,double max,int number_of_bins,bool logscale);
	~Histogram();
	double GetCount(size_t index);
	double GetFrequency(size_t index);
	//double GetNormalized(int index);
	double GetX(size_t index);
	void Add(const double &x,const double &dY);
	bool logarithmic;
	double /*max_count,*/total_count;
	void ResetCounts();
};