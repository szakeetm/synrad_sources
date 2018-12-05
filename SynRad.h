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
#include "Interface.h"
//#include "ExportDesorption.h"
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
	void InsertGeometry(bool newStr, char *fName = NULL);
	void SaveFile();

    // Util functions
	//void SendHeartBeat(bool forced=false);
    //void SetParam(GLTextField *txt,double value);
    void PlaceComponents();
	void LoadParam(char *fName=NULL,int position=-1);
	//void ExportDes(bool selectedOnly);
    void ClearFacetParams();
    void UpdateFacetParams(bool updateSelection=false);
    void ApplyFacetParams();
    void StartStopSimulation();
    void SaveConfig();
    void LoadConfig();
	
    void UpdateFacetHits(bool all=false);

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
	void PlaceScatteringControls(bool newReflectionMode);

	bool EvaluateVariable(VLIST * v);

    // Components

	GLTextField   *doseNumber;
    
	GLTextField   *facetTeleport;
    GLTextField   *facetSticking;
	GLTextField   *facetRMSroughness;
	GLTextField   *facetAutoCorrLength;
	GLCombo       *facetReflType;
	GLToggle      *facetDoScattering;
    GLTextField   *facetSuperDest;
    GLCombo       *facetProfileCombo;
	GLToggle	  *facetSpectrumToggle;
    GLButton      *facetTexBtn;
    GLLabel       *modeLabel;
	GLLabel       *doseLabel;
	GLLabel       *facetTPLabel;
	GLLabel       *facetRMSroughnessLabel;
	GLLabel       *facetAutoCorrLengthLabel;
    GLLabel       *facetLinkLabel;
    GLLabel       *facetStrLabel;
    GLTextField   *facetStructure;
   // GLLabel       *facetTLabel;
    GLLabel       *facetRLabel;
    GLLabel       *facetReLabel;
	GLLabel		*facetSpecLabel;

	GLMenu        *PARloadToMenu;
	GLMenu        *PARremoveMenu;
	GLMenu		  *ShowHitsMenu;

	//Materials (photon reflection)
	vector<string> materialPaths;

	void RebuildPARMenus();

    //Dialog
	RegionInfo       *regionInfo;
	//ExportDesorption *exportDesorption;
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
	void EmptyGeometry();
	
	void CrashHandler(Error *e);

protected:

    int  OneTimeSceneInit();
    int  RestoreDeviceObjects();
    int  InvalidateDeviceObjects();
	int  FrameMove();
    void ProcessMessage(GLComponent *src,int message);
};

