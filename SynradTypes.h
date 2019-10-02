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
#pragma once

#include "Vector.h"
//#include "Buffer_shared.h"
#include "GLApp/GLTypes.h"

// Hit type
#define HIT_DES   1
#define HIT_ABS   2
#define HIT_REF   3
#define HIT_TRANS 4
#define HIT_TELEPORTSOURCE 5
#define HIT_TELEPORTDEST 6
#define HIT_LAST 10

//Reflection type
#define REFL_ABSORB 0
#define REFL_FORWARD 1
#define REFL_DIFFUSE 2
#define REFL_BACK 3
#define REFL_TRANS 4

// Reflection type
#define REFLECTION_DIFFUSE 0   // Diffuse (cosine law)
#define REFLECTION_SPECULAR  1   // Mirror
#define REFLECTION_MATERIAL 10 //Real rough surface reflection

// Profile type

#define PROFILE_NONE  0  // No recording
#define PROFILE_U  1  // Flux, power (U direction)
#define PROFILE_V  2  // Flux, power (V direction)
#define PROFILE_ANGULAR    3  // Angular profile

// Density/Hit field stuff
#define HITMAX_INT64 18446744073709551615
#define HITMAX_DOUBLE 1E308

//Texture modes
#define TEXTURE_MODE_MCHITS 0
#define TEXTURE_MODE_FLUX 1
#define TEXTURE_MODE_POWER 2

class ProfileSlice {
public:
	size_t count_incident;
	size_t count_absorbed;
	double flux_incident;
	double flux_absorbed;
	double power_incident;
	double power_absorbed;
	ProfileSlice& operator+=(const ProfileSlice& rhs);
    ProfileSlice& operator=(const ProfileSlice& rhs);

};

class TextureCell {
public:
	size_t count;
	double flux;
	double power;
	TextureCell& operator+=(const TextureCell& rhs);
    TextureCell& operator=(const TextureCell& rhs);
};

class Trajectory_Point {
public:
	Vector3d position,direction,X_local,Y_local,Z_local,B;
	Vector3d rho;
	double critical_energy, emittance_X, emittance_Y, beta_X, beta_Y, eta, eta_prime, alpha_X, alpha_Y;
	double sigma_x, sigma_y, sigma_x_prime, sigma_y_prime, theta_X, theta_Y, gamma_X, gamma_Y;
	double a_x, b_x, a_y, b_y;
	
	double Critical_Energy(const double &gamma);
	double dAlpha(const double &dL);

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(
                CEREAL_NVP(position),
                CEREAL_NVP(direction),
                CEREAL_NVP(X_local),CEREAL_NVP(Y_local),CEREAL_NVP(Z_local),
                CEREAL_NVP(rho),

                CEREAL_NVP(critical_energy),
                CEREAL_NVP(emittance_X),CEREAL_NVP(emittance_Y),
                CEREAL_NVP(beta_X),CEREAL_NVP(beta_Y),
                CEREAL_NVP(eta),CEREAL_NVP(eta_prime),
                CEREAL_NVP(alpha_X),CEREAL_NVP(alpha_Y),
                CEREAL_NVP(sigma_x),CEREAL_NVP(sigma_y),
                CEREAL_NVP(sigma_x_prime),CEREAL_NVP(sigma_y_prime),
                CEREAL_NVP(theta_X),CEREAL_NVP(theta_Y),
                CEREAL_NVP(gamma_X),CEREAL_NVP(gamma_Y),
                CEREAL_NVP(a_x),CEREAL_NVP(a_y),
                CEREAL_NVP(b_x),CEREAL_NVP(b_y)
        );
    }
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
	//ProfileSlice *counts;
    std::vector<ProfileSlice> counts;
public:
    Histogram();
    Histogram(double min,double max,int number_of_bins,bool logscale);
	~Histogram();
	ProfileSlice GetCounts(size_t index);
	//double GetNormalized(int index);
	double GetX(size_t index);
	void Add(const double &x,const ProfileSlice &increment);
	bool logarithmic;
	void ResetCounts();
};