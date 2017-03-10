/*
  File:        SynRad.cpp
  Description: Main application class (GUI management)
  Program:     SynRad+
  Author:      R. KERSEVAN / M ADY
  Copyright:   E.S.R.F / CERN

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "Interface.h"
#include "ExportDesorption.h"
#include "FacetMesh.h"
#include "FacetDetails.h"
#include "TrajectoryDetails.h"
#include "Viewer3DSettings.h"
#include "TextureSettings.h"
#include "GlobalSettings.h"
#include "ProfilePlotter.h"
#include "SpectrumPlotter.h"
#include "TexturePlotter.h"
#include "RegionInfo.h"
#include "RegionEditor.h"

class Worker;

class SynRad : public Interface
{
public:
    SynRad();

	void LoadFile(char *fName = NULL);
	void InsertGeometry(BOOL newStr, char *fName = NULL);
	void SaveFile();

    // Util functions
	//void SendHeartBeat(BOOL forced=FALSE);
    //void SetParam(GLTextField *txt,double value);
    void PlaceComponents();
	void LoadParam(char *fName=NULL,int position=-1);
	//void ExportDes(bool selectedOnly);
    void ClearFacetParams();
    void UpdateFacetParams(BOOL updateSelection=FALSE);
    void ApplyFacetParams();
    void StartStopSimulation();
    void SaveConfig();
    void LoadConfig();
	
    void UpdateFacetHits(BOOL all=FALSE);

	void QuickPipe();

	void UpdatePlotters();
	
	void ClearRegions();
	void RemoveRegion(int index);
	void NewRegion();

    // Recent files   	
	char *recentPARs[MAX_RECENT];
    int  nbRecentPAR;
    void AddRecentPAR(char *fileName);
    void RemoveRecentPAR(char *fileName);
	void UpdateRecentPARMenu();
	void PlaceScatteringControls(BOOL newReflectionMode);

    // Components

	GLTextField   *doseNumber;
    
	GLTextField   *facetTeleport;
    GLTextField   *facetSticking;
	GLTextField   *facetRMSroughness;
	GLTextField   *facetAutoCorrLength;
	GLCombo       *facetReflType;
	GLToggle      *facetDoScattering;
    GLTextField   *facetSuperDest;
    GLCombo       *facetRecType;
	GLToggle      *facetSpectrumToggle;
    GLButton      *facetTexBtn;
    GLLabel       *modeLabel;
	GLLabel       *doseLabel;
	GLLabel       *facetTPLabel;
	GLLabel       *facetRMSroughnessLabel;
	GLLabel       *facetAutoCorrLengthLabel;
    GLLabel       *facetLinkLabel;
    GLLabel       *facetStrLabel;
    GLTextField   *facetSILabel;
   // GLLabel       *facetTLabel;
    GLLabel       *facetRLabel;
    GLLabel       *facetReLabel;

	GLMenu        *PARloadToMenu;
	GLMenu        *PARremoveMenu;
	GLMenu		  *ShowHitsMenu;

	//Materials (photon reflection)
	vector<string> materialPaths;

	void RebuildPARMenus();

    //Dialog
	RegionInfo       *regionInfo;
	ExportDesorption *exportDesorption;
    FacetMesh        *facetMesh;
    FacetDetails     *facetDetails;
	TrajectoryDetails *trajectoryDetails;
    Viewer3DSettings  *viewer3DSettings;
    TextureSettings  *textureSettings;
	GlobalSettings	 *globalSettings;
    ProfilePlotter   *profilePlotter;
	SpectrumPlotter  *spectrumPlotter;
    TexturePlotter   *texturePlotter;
	RegionEditor     *regionEditor;

	//char *nbF;

    // Testing
    //int     nbSt;
    //void LogProfile();
    void BuildPipe(double ratio,int steps=0);
    //void BuildPipeStick(double s);
	
	void CrashHandler(Error *e);

protected:

    int  OneTimeSceneInit();
    int  RestoreDeviceObjects();
    int  InvalidateDeviceObjects();
	int  FrameMove();
    void ProcessMessage(GLComponent *src,int message);
	BOOL EvaluateVariable(VLIST * v, Worker * w, Geometry * geom);
};

