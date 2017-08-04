#ifndef _REGION_FULL_
#define _REGION_FULL_

#include <stdio.h>
#include <fstream>
#include "Region_mathonly.h"
#include <Process.h> //for getpid

#include "File.h"

#include "GLApp\GLToolkit.h"
#include "GLApp\GLProgress.h"
//#include "GLApp\GLTypes.h"

class Region_full : public Region_mathonly { //Beam trajectory
public:
	
	//References
	std::string fileName,MAGXfileName,MAGYfileName,MAGZfileName,BXYfileName;
	bool isLoaded;
	Vector3d AABBmin,AABBmax;
	int selectedPoint;

	//Methods
	Region_full();
	~Region_full();
	//Region_full(const Region_full &src); //used at push_back
	//Region_full& operator=(const Region_full &src); //used in CopyGeometryBuffer
	//need converter between full<->mathonly?
	
	void CalculateTrajectory(int max_steps);
	Trajectory_Point OneStep(int pointId); //moves the beam by dL length from the i-th calculated point
	void CalcPointProperties(int pointId); //calculates magnetic field, sigma values, etc. for a trajectory point
	bool isOutsideBoundaries(Vector3d a,bool recalcDirs);
	void LoadPAR(FileReader *file);
	void LoadParam(FileReader *f);
	Distribution2D LoadMAGFile(FileReader *file,Vector3d *dir,double *period,double *phase,int mode);
	int LoadBXY(const std::string& fileName); //Throws error
	void Render(const int& regionId, const int& dispNumTraj, GLMATERIAL *B_material, const double& vectorLength);
	void SelectTrajPoint(int x,int y, size_t regionId);
	void SaveParam(FileWriter *f);
	
};

#endif