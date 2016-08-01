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



class Region_full : public Region_mathonly { //Beam trajectory0
public:
	
	//References
	std::string fileName,MAGXfileName,MAGYfileName,MAGZfileName,BXYfileName;
	bool isLoaded;
	Vector AABBmin,AABBmax;
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
	bool isOutsideBoundaries(Vector a,BOOL recalcDirs);
	void LoadPAR(FileReader *file);
	void LoadParam(FileReader *f);
	void LoadRegion(FileReader *file,GLProgress *prg);
	Distribution2D LoadMAGFile(FileReader *file,Vector *dir,double *period,double *phase,int mode);
	int LoadBXY(FileReader *file);
	void Render(int dispNumTraj,GLMATERIAL *B_material,double vectorLength);
	void SelectTrajPoint(int x,int y,int regionId);
	void ExportPoints(FileWriter *f,GLProgress *prg,int frequency=1,BOOL doFullScan=FALSE);
	void SaveParam(FileWriter *f);
	
};

#endif