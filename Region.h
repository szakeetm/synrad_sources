#ifndef _B_TRAJ_
#define _B_TRAJ_

#include <math.h>
#include <stdio.h>
#include <fstream>
#include "Tools.h"
#include <Process.h> //for getpid
#include <vector>
#include "File.h"
#include "Distributions.h"
#include "Quadrupole.h"
#include "GLApp\GLToolkit.h"
#include "GLApp\GLProgress.h"
//#include "GLApp\GLTypes.h"

#define B_MODE_CONSTANT        1
#define B_MODE_COORDVALUES     2
#define B_MODE_SINCOS          3
#define B_MODE_QUADRUPOLE      4
#define B_MODE_ANALYTIC        5
#define B_MODE_HELICOIDAL      6
#define B_MODE_ROTATING_DIPOLE 7

class Region { //Beam trajectory
public:
	//Parameters
	double energy_low;
	double energy_hi;
	double dL;
	double particleMass;
	double E;
	/*
	double x0,y0,z0;
	double teta0,alfa0;
	*/
	double emittance;
	double betax,betay;
	double eta,etaprime;
	double psimaxX,psimaxY;
	double coupling;
	double energy_spread;
	
	int generation_mode;     //default value, can change in program
	/*
	int no_scans;            //from old Synrad code, unused
	int nregions;            //from old Synrad code, unused
	int i_struct1;           //from old Synrad code, unused
	*/
	int Bx_mode,By_mode,Bz_mode,beta_kind; //switches between constant or file-based magnetic fields

	bool enable_ort_polarization,enable_par_polarization;
	
	//References
	std::string fileName,MAGXfileName,MAGYfileName,MAGZfileName,BXYfileName;

	//Loaded from files
	Vector nbDistr_MAG; //MAG file points
	int nbDistr_BXY; //BXY file points
	Distribution2D *Bx_distr,*By_distr,*Bz_distr; //pointer to the B field distribution (if file-based)
	Distribution2D *beta_x_distr,*beta_y_distr,*eta_distr,*etaprime_distr; //pointer to the beta field distribution (if file-based)
	Vector Bx_dir,By_dir,Bz_dir; //direction in which MAG files are oriented (in their second line)
	double Bx_period,By_period,Bz_period; //first line of length files
	double Bx_phase,By_phase,Bz_phase; //phase in case of helicoidal B file

	//Calculated data
	std::vector<Trajectory_Point> Points;
	int nbPoints;
	bool isLoaded;
	Vector startPoint,startDir,B_const,limits,AABBmin,AABBmax;
	Quadrupole quad;
	int selectedPoint;
	double gamma;  //    =E/abs(particleMass)
	double emittance_x,emittance_y;

	//Methods
	Region();
	~Region();
	Region(const Region &src); //used at push_back
	Region& operator=(const Region &src); //used in CopyGeometryBuffer
	
	std::vector<Trajectory_Point> CalculateTrajectory(double dL,double E,Vector limits,Vector pos0,Vector dir0,int max_steps);
	Trajectory_Point OneStep(Trajectory_Point p0,double dL,double E); //moves the beam by dL length
	bool isOutsideBoundaries(Vector a,BOOL recalcDirs);
	Vector B(Vector position); //returns the B field at position
	void LoadPAR(FileReader *file,GLProgress *prg);
	void LoadParam(FileReader *f,GLProgress *prg);
	void LoadRegion(FileReader *file,GLProgress *prg);
	Distribution2D LoadMAGFile(FileReader *file,Vector *dir,double *period,double *phase,int mode);
	int LoadBXY(FileReader *file,Distribution2D *beta_x_distr,Distribution2D *beta_y_distr,
		Distribution2D *eta_distr,Distribution2D *etaprime_distr);
	void Render(int dispNumTraj,GLMATERIAL *B_material,double vectorLength);
	void SelectTrajPoint(int x,int y);
	void ExportPoints(FileWriter *f,GLProgress *prg,int frequency=1,BOOL doFullScan=FALSE);
	void SaveParam(FileWriter *f);
	
};

#endif