/*
  File:        SynRad.cpp
  Description: Main application class (GUI management)
  Program:     SynRad
  Author:      R. KERSEVAN / M SZAKACS
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

#include "GLApp/GLApp.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLList.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLMenuBar.h"

#include "Worker.h"
#include "GeometryViewer.h"
#include "FormulaSettings.h"
#include "CollapseSettings.h"
#include "MoveVertex.h"
#include "ScaleVertex.h"
#include "ScaleFacet.h"
#include "MoveFacet.h"
#include "ExportDesorption.h"
#include "MirrorFacet.h"
#include "RotateFacet.h"
#include "AlignFacet.h"
#include "AddVertex.h"
#include "FacetMesh.h"
#include "FacetCoordinates.h"
#include "VertexCoordinates.h"
#include "FacetDetails.h"
#include "Viewer3DSettings.h"
#include "TextureSettings.h"
#include "GlobalSettings.h"
#include "ProfilePlotter.h"
#include "SpectrumPlotter.h"
#include "ViewEditor.h"
#include "TexturePlotter.h"
//#include "OutgassingMap.h"
#include "SelectDialog.h"
#include "RegionInfo.h"
#include "RegionEditor.h"

#define MAX_FORMULA 10
#define MAX_VIEW    19
#define MAX_SELECTION 19
#define MAX_RECENT  10

#ifdef _DEBUG
#define APP_NAME "SynRad+ development version (Compiled "__DATE__" "__TIME__") DEBUG MODE"
#else
//#define APP_NAME "SynRad+ development version ("__DATE__")"
#define APP_NAME "Synrad+ 1.1.0.3 ("__DATE__")"
#endif

extern int changedSinceSave;

typedef struct {
  GLLabel     *name;
  GLTextField *value;
  GLButton    *setBtn;
  GLParser    *parser;
} FORMULA;

class SynRad : public GLApplication
{
public:
    SynRad();

    // Worker handle
    Worker worker;

    // Current directory
    void UpdateCurrentDir(char *fileName);
    char currentDir[512];
    void UpdateCurrentSelDir(char *fileName);
    char currentSelDir[512];

    // Simulation state
    float    lastUpdate;   // Last 'hit update' time
    double   hps;          // Hit per second
    double   dps;          // Desorption (or new particle) per second
    double   lastHps;      // hps measurement
    double   lastDps;      // dps measurement
    llong    lastNbHit;    // measurement
    llong    lastNbDes;    // measurement
    llong    nbDesStart;   // measurement
    llong    nbHitStart;   // measurement
    int      nbProc;       // Temporary var (use Worker::GetProcNumber)
	//float lastHeartBeat;   //last time a heartbeat was sent to the subprocesses
	float lastAppTime;

    float    lastMeasTime; // Last measurement time (for hps and dps)
	double tolerance; //Select coplanar tolerance

    // Util functions
	//void SendHeartBeat(BOOL forced=FALSE);
    void SetParam(GLTextField *txt,double value);
    char *FormatInt(llong v,char *unit);
    char *FormatPS(double v,char *unit);
    char *FormatSize(DWORD size);
    char *FormatTime(float t);
    void PlaceComponents();
    void LoadSelection(char *fName=NULL);
    void SaveSelection();
    void LoadFile(char *fName=NULL);
	void LoadParam(char *fName=NULL,int position=-1);
	void LoadFileMemory();
	void InsertGeometry(BOOL newStr,char *fName=NULL);
	void SaveFile();
    void SaveFileAs();
    void ExportSelection();
	void ExportTexture(int mode);
	//void ExportDes(bool selectedOnly);
    void ClearFacetParams();
    void UpdateFacetParams(BOOL updateSelection=FALSE);
    void ApplyFacetParams();
    void UpdateModelParams();
	void UpdateViewerFlags();
    void StartStopSimulation();
    void ResetSimulation(BOOL askConfirm=TRUE);
    void EditFacet();
    void SaveConfig();
    void LoadConfig();
    void UpdateStructMenu();
    void UpdateTitle();
    void UpdateFacetHits(BOOL all=FALSE);
    void AnimateViewerChange(int next,BOOL init=FALSE);
    void UpdateViewerParams();
    void SelectViewer(int s);
    void Place3DViewer();
	void QuickPipe();
	void UpdateFacetlistSelected();
	BOOL AskToSave();
	BOOL AskToReset(Worker *work=NULL);
	float GetAppTime();
	void ResetAutoSaveTimer();
	void AutoSave(BOOL crashSave=FALSE);
	void CheckForRecovery();
	void UpdateViewers();
	void SetFacetSearchPrg(BOOL visible,char *text);
	void DisplayCollapseDialog();
	void ClearTraj();
	void RemoveRegion(int index);
	void RenumberSelections(int startFacetId);
	void NewRegion();

    // Formula management
    void ProcessFormulaButtons(GLComponent *src);
    void UpdateFormula();
	void OffsetFormula(char* expression,int offset,int filter=0);
	void RenumberFormulas(int startId);
    void AddFormula(GLParser *f,BOOL doUpdate=TRUE);
	
    // Recent files
	float lastSaveTime;
	float lastSaveTimeSimu;
   
	char *recents[MAX_RECENT];
    int  nbRecent;
    void AddRecent(char *fileName);
    void RemoveRecent(char *fileName);
	
	char *recentPARs[MAX_RECENT];
    int  nbRecentPAR;
    void AddRecentPAR(char *fileName);
    void RemoveRecentPAR(char *fileName);

    // Components
    GLMenuBar     *menu;
    GeometryViewer *viewer[MAX_VIEWER];
    GLTextField   *geomNumber;
    GLToggle      *showNormal;
    GLToggle      *showRule;
    GLToggle      *showUV;
    GLToggle      *showLeak;
    GLToggle      *showHit;
    GLToggle      *showLine;
    GLToggle      *showVolume;
    GLToggle      *showTexture;
    GLToggle      *showFilter;
    GLToggle      *showIndex;
    GLToggle      *showVertex;
    GLButton      *showMoreBtn;
    GLButton      *startSimu;
    GLButton      *resetSimu;
    GLCombo       *modeCombo;
    GLTextField   *hitNumber;
    GLTextField   *desNumber;
    GLTextField   *leakNumber;
	GLTextField   *doseNumber;
    GLTextField   *sTime;
    GLMenu        *facetMenu;
	GLTextField   *facetTeleport;
    GLTextField   *facetSticking;
	GLTextField   *facetRoughness;
    GLTextField   *facetSuperDest;
    GLTextField   *facetOpacity;
	GLTextField   *facetArea;
	GLCombo       *facetSideType;
    GLCombo       *facetReflType;
    GLCombo       *facetRecType;
    GLButton      *facetApplyBtn;
    GLButton      *facetMoreBtn;
    GLButton      *facetCoordBtn;
    GLButton      *facetTexBtn;
    GLTitledPanel *facetPanel;
    GLList        *facetList;
    GLTitledPanel *togglePanel;
    GLLabel       *modeLabel;
	GLLabel       *facetAreaLabel;
    GLButton      *singleACBtn;
    GLLabel       *hitLabel;
    GLLabel       *desLabel;
    GLLabel       *leakLabel;
	GLLabel       *doseLabel;
    GLLabel       *sTimeLabel;
    GLTitledPanel *simuPanel;
	GLLabel       *facetTPLabel;
    GLLabel       *facetStickingLabel;
	GLLabel       *facetRoughnessLabel;
	GLLabel       *facetSideLabel;
    GLLabel       *facetLinkLabel;
    GLLabel       *facetStrLabel;
    GLTextField   *facetSILabel;
    GLLabel       *facetTLabel;
    GLLabel       *facetRLabel;
    GLLabel       *facetReLabel;
	GLLabel       *facetSpectrumLabel;
	GLCombo      *facetSpectrumCombo;
    GLMenu        *structMenu;
    GLMenu        *viewsMenu;
	GLMenu        *selectionsMenu;
	GLMenu        *memorizeSelectionsMenu;
	GLMenu        *memorizeViewsMenu;
	GLMenu        *clearSelectionsMenu;
	GLMenu        *clearViewsMenu;
	GLMenu        *PARloadToMenu;
	GLMenu        *PARremoveMenu;

    // Formulas
    FORMULA formulas[MAX_FORMULA];
    int nbFormula;
    void ClearFormula();
    void AddFormula(char *fName,char *formula);

	//Materials (photon reflection)
	vector<string> materialPaths;

    // Views
    void SelectView(int v);
	void AddView(char *selectionName,AVIEW v);
    void AddView();
    void ClearViewMenus();
	void ClearAllViews();
	void OverWriteView(int idOvr);
	void ClearView(int idClr);
	void RebuildViewMenus();
	
    // Selections
    void SelectSelection(int v);
    void AddSelection(char *selectionName,ASELECTION s);
    void AddSelection();
    void ClearSelectionMenus();
	void ClearAllSelections();
	void OverWriteSelection(int idOvr);
	void ClearSelection(int idClr);
	void RebuildSelectionMenus();

	void RebuildPARMenus();

	void AddStruct();
	void DeleteStruct();

    AVIEW   views[MAX_VIEW];
    int     nbView;
    int     idView;
    int     curViewer;
    int     modeSolo;

	ASELECTION selections[MAX_SELECTION];
	int nbSelection;
	int idSelection;

    //Dialog
    FormulaSettings  *formulaSettings;
    CollapseSettings *collapseSettings;
	RegionInfo       *regionInfo;
	MoveVertex		 *moveVertex;
	ScaleVertex      *scaleVertex;
	ScaleFacet       *scaleFacet;
	SelectDialog     *selectDialog;
	MoveFacet		 *moveFacet;
	ExportDesorption *exportDesorption;
	MirrorFacet	     *mirrorFacet;
	RotateFacet      *rotateFacet;
	AlignFacet       *alignFacet;
	AddVertex		 *addVertex;
    FacetMesh        *facetMesh;
    FacetDetails     *facetDetails;
    Viewer3DSettings  *viewer3DSettings;
    TextureSettings  *textureSettings;
	GlobalSettings	 *globalSettings;
    FacetCoordinates *facetCoordinates;
	VertexCoordinates *vertexCoordinates;
    ProfilePlotter   *profilePlotter;
	SpectrumPlotter  *spectrumPlotter;
    ViewEditor       *viewEditor;
    TexturePlotter   *texturePlotter;
	//OutgassingMap    *outgassingMap;
	RegionEditor     *regionEditor2;
	char *nbF;

    // Testing
    //int     nbSt;
    //void LogProfile();
    void BuildPipe(double ratio);
    //void BuildPipeStick(double s);
	
	void CrashHandler(Error *e);

protected:

    int  OneTimeSceneInit();
    int  RestoreDeviceObjects();
    int  InvalidateDeviceObjects();
    int  OnExit();
    int  FrameMove();
    void ProcessMessage(GLComponent *src,int message);
    int  Resize(DWORD width, DWORD height, BOOL forceWindowed);

};

