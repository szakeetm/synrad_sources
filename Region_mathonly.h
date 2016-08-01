//Generic region that can do the calculations but not additional methods like display, etc.

#ifndef _REGION_MATHONLY_
#define _REGION_MATHONLY_

#include <math.h>
#include <stdio.h>
#include <fstream>
#include "Tools.h"
#include <vector>
#include "File.h"
#include "Distributions.h"
#include "Quadrupole.h"
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
	double energy_low;
	double energy_hi;
	double dL;
	double particleMass;
	double E;
	double current;
	double emittance;
	double betax, betay;
	double eta, etaprime;
	double psimaxX, psimaxY;
	double coupling;
	double energy_spread;

	int Bx_mode, By_mode, Bz_mode, beta_kind; //switches between constant or file-based magnetic fields

	bool enable_ort_polarization, enable_par_polarization;
	Vector nbDistr_MAG; //MAG file points
	size_t nbDistr_BXY; //BXY file points
	Vector Bx_dir,By_dir,Bz_dir; //direction in which MAG files are oriented (in their second line)
	double Bx_period,By_period,Bz_period; //first line of length files
	double Bx_phase,By_phase,Bz_phase; //phase in case of helicoidal B file

	size_t nbPointsToCopy; //passes info about the number of points that needs to be read from the buffer on loading
	Vector startPoint, startDir, B_const, limits;//AABBmin,AABBmax
	Quadrupole quad;
	//int selectedPoint;
	double gamma;  //    =E/abs(particleMass)
};

class Region_mathonly { //Beam trajectory
public:
	RegionParams params; //everything except distributions

	//References
	//std::string fileName,MAGXfileName,MAGYfileName,MAGZfileName,BXYfileName;

	//Loaded from files
	Distribution2D Bx_distr,By_distr,Bz_distr; //B field distribution (if file-based)
	Distribution2D beta_x_distr,beta_y_distr,eta_distr,etaprime_distr,e_spread_distr; //beta field distribution (if file-based)

	//Calculated data
	std::vector<Trajectory_Point> Points;

	//Methods
	Region_mathonly();
	~Region_mathonly();
	Region_mathonly(const Region_mathonly &src);
	Region_mathonly& operator=(const Region_mathonly &src);

	Vector B(int pointId,const Vector &offset); //returns the B field at a given point, allows interpolation between points and offset due to non-ideal beam
};

#endif