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
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>

class Region_full : public Region_mathonly { //Beam trajectory
public:
	
	//References
	std::string fileName,MAGXfileName,MAGYfileName,MAGZfileName,BXYfileName;
	bool isLoaded;
	Vector3d AABBmin,AABBmax;
	int selectedPointId; //can be -1 if nothing's selected

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
	void Render(const size_t& regionId, const size_t& dispNumTraj, GLMATERIAL *B_material, const double& vectorLength);
	void SelectTrajPoint(int x,int y, size_t regionId);
	void SaveParam(FileWriter *f);

    /*template<class Archive>
    void serialize(Archive & archive)
    {

        archive(
                CEREAL_NVP(fileName),
                CEREAL_NVP(MAGXfileName),CEREAL_NVP(MAGYfileName),CEREAL_NVP(MAGZfileName),
                CEREAL_NVP(BXYfileName),
                CEREAL_NVP(isLoaded),
                CEREAL_NVP(AABBmin),
                CEREAL_NVP(AABBmax),
                CEREAL_NVP(selectedPointId),

                CEREAL_NVP(params),
                CEREAL_NVP(Points),
                CEREAL_NVP(Bx_distr),CEREAL_NVP(By_distr),CEREAL_NVP(Bx_distr),
                CEREAL_NVP(latticeFunctions)
        );
    }*/
};

#endif