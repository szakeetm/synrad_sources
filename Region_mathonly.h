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
#include <cereal/types/vector.hpp>

#define B_MODE_CONSTANT        1
#define B_MODE_DIRECTION       2
#define B_MODE_ALONGBEAM       3
#define B_MODE_SINCOS          4
#define B_MODE_QUADRUPOLE      5
#define B_MODE_ANALYTIC        6
#define B_MODE_HELICOIDAL      7
#define B_MODE_ROTATING_DIPOLE 8
#define B_MODE_COMBINED_FUNCTION 9

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

	bool enable_ort_polarization, enable_par_polarization;
	size_t polarizationCompIndex;

	Vector3d nbDistr_MAG; //MAG file points
	size_t nbDistr_BXY; //BXY file points
	Vector3d Bx_dir,By_dir,Bz_dir; //direction in which MAG files are oriented (in their second line)
	double Bx_period,By_period,Bz_period; //first line of length files
	double Bx_phase,By_phase,Bz_phase; //phase in case of helicoidal B file

	size_t nbPointsToCopy; //passes info about the number of points that needs to be read from the buffer on loading
	Vector3d startPoint, startDir, B_const, limits;//AABBmin,AABBmax
	Quadrupole quad_params;
	bool showPhotons; //Whether to include photons from this region in the hit cache
	size_t structureId; //Which structure generated photons belong to

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(
                CEREAL_NVP(energy_low_eV),
        CEREAL_NVP(energy_hi_eV),
        CEREAL_NVP(dL_cm),
        CEREAL_NVP(particleMass_GeV),
        CEREAL_NVP(E_GeV),
        CEREAL_NVP(current_mA),
        CEREAL_NVP(emittance_cm),
        CEREAL_NVP(betax_const_cm), CEREAL_NVP(betay_const_cm),
        CEREAL_NVP(eta_x_const_cm), CEREAL_NVP(eta_x_prime_const),
        CEREAL_NVP(psimaxX_rad), CEREAL_NVP(psimaxY_rad),
        CEREAL_NVP(coupling_percent),
        CEREAL_NVP(energy_spread_percent),

        CEREAL_NVP(gamma), // Calculated:   gamma=E/abs(particleMass)

        CEREAL_NVP(Bx_mode), CEREAL_NVP(By_mode), CEREAL_NVP(Bz_mode), CEREAL_NVP(beta_kind),//switches between constant or file-based magnetic fields

        CEREAL_NVP(enable_ort_polarization), CEREAL_NVP(enable_par_polarization),
                CEREAL_NVP(polarizationCompIndex),

        CEREAL_NVP(nbDistr_MAG),//MAG file points
        CEREAL_NVP(nbDistr_BXY),//BXY file points
        CEREAL_NVP(Bx_dir), CEREAL_NVP(By_dir), CEREAL_NVP(Bz_dir),//direction in which MAG files are oriented (in their second line)
        CEREAL_NVP(Bx_period), CEREAL_NVP(By_period), CEREAL_NVP(Bz_period),//first line of length files
        CEREAL_NVP(Bx_phase), CEREAL_NVP(By_phase), CEREAL_NVP(Bz_phase),//phase in case of helicoidal B file

                CEREAL_NVP(nbPointsToCopy),//passes info about the number of points that needs to be read from the buffer on loading
        CEREAL_NVP(startPoint), CEREAL_NVP(startDir), CEREAL_NVP(B_const), CEREAL_NVP(limits),//AABBmin,AABBmax
        CEREAL_NVP(quad_params),
        CEREAL_NVP(showPhotons),//Whether to include photons from this region in the hit cache
        CEREAL_NVP(structureId)//Which structure generated photons belong to
        );
    }
};

class Region_mathonly { //Beam trajectory
public:
	RegionParams params; //everything except distributions (Plain Old Data)
	//References
	//std::string fileName,MAGXfileName,MAGYfileName,MAGZfileName,BXYfileName;

	//Loaded from files
	Distribution2D Bx_distr,By_distr,Bz_distr; //B field distribution (if file-based)
	DistributionND latticeFunctions; //BetaX, BetaY, EtaX, EtaX', AlphaX, AlphaY

	//Calculated data
	std::vector<Trajectory_Point> Points;

	//Methods
	Region_mathonly();
	//Region_mathonly(const Region_mathonly &src);
	//Region_mathonly& operator=(const Region_mathonly &src);

	Vector3d B(size_t pointId,const Vector3d &offset); //returns the B field at a given point, allows interpolation between points and offset due to non-ideal beam

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(
                CEREAL_NVP(params),
                CEREAL_NVP(Points),
                CEREAL_NVP(Bx_distr),CEREAL_NVP(By_distr),CEREAL_NVP(Bx_distr),
                CEREAL_NVP(latticeFunctions)
        );
    }
};

#endif