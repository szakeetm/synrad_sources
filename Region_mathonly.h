//Generic region that can do the calculations but not additional methods like display, etc.

#ifndef _REGION_MATHONLY_
#define _REGION_MATHONLY_

#include <math.h>
#include <stdio.h>
#include <fstream>
#include "SynradTypes.h"
#include <vector>
//#include "File.h"
#include "Quadrupole.h"
#include "Distributions.h"
#include "GLApp\GLTypes.h"
//#include "GLApp\GLTypes.h"

#define B_MODE_CONSTANT        1
#define B_MODE_DIRECTION       2
#define B_MODE_ALONGBEAM       3
#define B_MODE_SINCOS          4
#define B_MODE_QUADRUPOLE      5
#define B_MODE_ANALYTIC        6
#define B_MODE_HELICOIDAL      7
#define B_MODE_ROTATING_DIPOLE 8

struct RegionParams {

	//Parameters
	double energy_low_eV;
	double energy_hi_eV;
	double dL_cm;
	double particleMass_GeV;
	double E_GeV;
	double current_mA;
	double emittance_cm;
	double betax_const_cm, betay_const_cm;
	double eta_x_const_cm, eta_x_prime_const;
	double psimaxX_rad, psimaxY_rad;
	double coupling_percent;
	double energy_spread_percent;
	
	double gamma;  // Calculated:   gamma=E/abs(particleMass)

	int Bx_mode, By_mode, Bz_mode, beta_kind; //switches between constant or file-based magnetic fields

	BOOL enable_ort_polarization, enable_par_polarization;
	Vector3d nbDistr_MAG; //MAG file points
	size_t nbDistr_BXY; //BXY file points
	Vector3d Bx_dir,By_dir,Bz_dir; //direction in which MAG files are oriented (in their second line)
	double Bx_period,By_period,Bz_period; //first line of length files
	double Bx_phase,By_phase,Bz_phase; //phase in case of helicoidal B file

	size_t nbPointsToCopy; //passes info about the number of points that needs to be read from the buffer on loading
	Vector3d startPoint, startDir, B_const, limits;//AABBmin,AABBmax
	Quadrupole quad_params;
	BOOL showPhotons; //Whether to include photons from this region in the hit cache
};

class Region_mathonly { //Beam trajectory
public:
	RegionParams params; //everything except distributions
	//References
	//std::string fileName,MAGXfileName,MAGYfileName,MAGZfileName,BXYfileName;

	//Loaded from files
	Distribution2D Bx_distr,By_distr,Bz_distr; //B field distribution (if file-based)
	DistributionND betaFunctions; //BetaX, BetaY, EtaX, EtaX', AlphaX, AlphaY

	//Calculated data
	std::vector<Trajectory_Point> Points;

	//Methods
	Region_mathonly();
	Region_mathonly(const Region_mathonly &src);
	Region_mathonly& operator=(const Region_mathonly &src);

	Vector3d B(int pointId,const Vector3d &offset); //returns the B field at a given point, allows interpolation between points and offset due to non-ideal beam
};

#endif