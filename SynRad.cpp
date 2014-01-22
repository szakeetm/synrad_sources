/*
File:        SynRad.cpp
Description: Main application class (GUI management)
Program:     SynRad
Author:      Roberto KERSEVAN / Marton ADY
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

#include <math.h>
#include <malloc.h>
#include "SynRad.h"
#include "File.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLSaveDialog.h"
#include "GLApp/GLInputBox.h"
#include "GLApp/GLFileBox.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "Utils.h" //for Remainder
#include "direct.h"
//#include <winsparkle.h>

float m_fTime;

static const char *fileLFilters = "All SynRad supported files\0*.txt;*.syn;*.syn7z;*.geo;*.geo7z;*.str;*.stl;*.ase\0SYN files\0*.syn;*.syn7z;\0GEO files\0*.geo;*.geo7z;\0TXT files\0*.txt\0STR files\0*.str\0STL files\0*.stl\0ASE files\0*.ase\0";
static const int   nbLFilter = 8;
static const char *fileInsFilters = "SYN files\0*.syn;*.syn7z;\0GEO files\0*.geo;*.geo7z;\0Text files\0*.txt\0STL files\0*.stl\0";
static const int   nbInsFilter = 4;
static const char *fileSFilters = "SYN files\0*.syn;*.syn7z;\0GEO files\0*.geo;*.geo7z;\0Text files\0*.txt\0All files\0*.*\0";
static const int   nbSFilter = 4;
static const char *fileSelFilters = "Selection files\0*.sel\0All files\0*.*\0";
static const int   nbSelFilter = 2;
static const char *fileTexFilters = "Text files\0*.txt\0Texture files\0*.tex\0All files\0*.*\0";
static const int   nbTexFilter = 3;
static const char *fileParFilters = "param files\0*.param\0PAR files\0*.par\0All files\0*.*\0";
static const int   nbParFilter = 3;
static const char *fileDesFilters = "DES files\0*.des\0All files\0*.*\0";
static const int   nbDesFilter = 2;

static const int   cWidth[] = {30,56,50,50,50};
static const char *cName[] = {"#","Hits","Flux","Power","Abs"};

BOOL EndsWithParam(const char* s);
BOOL changedSinceSave;
SynRad *theApp;
extern double autoSaveFrequency;
extern int checkForUpdates;
extern int compressSavedFiles;
extern int autoSaveSimuOnly;
extern int numCPU;

#define MENU_FILE_LOAD       11
#define MENU_FILE_SAVE       12
#define MENU_FILE_SAVEAS     13
#define MENU_FILE_INSERTGEO  140
#define MENU_FILE_INSERTGEO_NEWSTR  141
#define MENU_FILE_EXPORTMESH      15
#define MENU_FILE_EXPORT_DESORP 18

#define MENU_FILE_EXPORTTEXTURE_AREA 151
#define MENU_FILE_EXPORTTEXTURE_MCHITS 152
#define MENU_FILE_EXPORTTEXTURE_FLUX 153
#define MENU_FILE_EXPORTTEXTURE_POWER 154
#define MENU_FILE_EXPORTTEXTURE_FLUXPERAREA 155
#define MENU_FILE_EXPORTTEXTURE_POWERPERAREA 156

/*#define MENU_FILE_EXPORTTEXTURE_DENSITY 152
#define MENU_FILE_EXPORTTEXTURE_AVG_VELOCITY 153
#define MENU_FILE_EXPORTTEXTURE_VELOCITY 154
#define MENU_FILE_EXPORTTEXTURE_COUNT 155
#define MENU_FILE_EXPORTTEXTURE_HITS 156*/

#define MENU_FILE_LOADRECENT 110
#define MENU_FILE_EXIT       20

#define MENU_REGIONS_LOADPAR      60
#define MENU_REGIONS_LOADRECENT   600
#define MENU_REGIONS_CLEARALL    61
#define MENU_REGIONS_LOADTO       630
#define MENU_REGIONS_REMOVE       660
#define MENU_REGIONS_REGIONINFO 62
#define MENU_REGIONS_NEW        63

#define MENU_EDIT_3DSETTINGS   21
#define MENU_EDIT_TSCALING     22
#define MENU_EDIT_ADDFORMULA   23
#define MENU_EDIT_GLOBALSETTINGS 24

#define MENU_FACET_COLLAPSE    300
#define MENU_FACET_SWAPNORMAL  301
#define MENU_FACET_SHIFTVERTEX 302
#define MENU_FACET_COORDINATES 303
#define MENU_FACET_PROFPLOTTER 304
#define MENU_FACET_SPECTRUMPLOTTER 305
#define MENU_FACET_DETAILS     306
#define MENU_FACET_MESH        307
#define MENU_FACET_TEXPLOTTER  308
#define MENU_FACET_REMOVESEL   309
#define MENU_FACET_EXPLODE     310
#define MENU_FACET_SELECTALL   311
#define MENU_FACET_SELECTSTICK 312
#define MENU_FACET_SELECTDES   313
#define MENU_FACET_SELECTABS   314
#define MENU_FACET_SELECTTRANS 315
#define MENU_FACET_SELECTREFL  316
#define MENU_FACET_SELECT2SIDE 317
#define MENU_FACET_SELECTTEXT  318
#define MENU_FACET_SELECTPROF  319
#define MENU_FACET_SELECTSPECTRUM 320
#define MENU_FACET_SELECTDEST  321
#define MENU_FACET_SELECTVOL   322
#define MENU_FACET_SELECTERR   323
#define MENU_FACET_SELECTHITS  324
#define MENU_FACET_SELECTTP    325
#define MENU_FACET_SAVESEL     330
#define MENU_FACET_LOADSEL     331
#define MENU_FACET_INVERTSEL   332
#define MENU_FACET_MOVE		   333
#define MENU_FACET_SCALE       334
#define MENU_FACET_MIRROR	   335
#define MENU_FACET_ROTATE	   336
#define MENU_FACET_ALIGN       337
//#define MENU_FACET_OUTGASSINGMAP 338
#define MENU_FACET_CREATE_DIFFERENCE 3381

#define MENU_SELECTION_ADDNEW             339
#define MENU_SELECTION_CLEARALL           340

#define MENU_SELECTION_MEMORIZESELECTIONS   3300
#define MENU_SELECTION_SELECTIONS           3400
#define MENU_SELECTION_CLEARSELECTIONS      3500

#define MENU_SELECTION_SELECTFACETNUMBER 360

#define MENU_VERTEX_SELECTALL   701
#define MENU_VERTEX_UNSELECTALL 702
#define MENU_VERTEX_SELECT_ISOLATED 703
#define MENU_VERTEX_CREATE_POLY_CONVEX   7040
#define MENU_VERTEX_CREATE_POLY_ORDER    7041
#define MENU_VERTEX_SELECT_COPLANAR   705
#define MENU_VERTEX_MOVE   706
#define MENU_VERTEX_SCALE  707
#define MENU_VERTEX_ADD	   708
#define MENU_VERTEX_REMOVE 709
#define MENU_VERTEX_COORDINATES 710


#define MENU_VIEW_STRUCTURE       4000
#define MENU_VIEW_STRUCTURE_P     40
#define MENU_VIEW_NEWSTRUCT       401
#define MENU_VIEW_DELSTRUCT       402
#define MENU_VIEW_PREVSTRUCT	  403		
#define MENU_VIEW_NEXTSTRUCT	  404
#define MENU_VIEW_FULLSCREEN      41

#define MENU_VIEW_ADDNEW             431
#define MENU_VIEW_CLEARALL           432

#define MENU_VIEW_MEMORIZEVIEWS   4300
#define MENU_VIEW_VIEWS           4400
#define MENU_VIEW_CLEARVIEWS      4500

#define MENU_TEST_PIPE0001        50
#define MENU_TEST_PIPE1           51
#define MENU_TEST_PIPE10          52
#define MENU_TEST_PIPE100         53
#define MENU_TEST_PIPE1000        54
#define MENU_TEST_PIPE10000       55

#define MENU_QUICKPIPE            56

//-----------------------------------------------------------------------------
// Name: WinMain()
// Desc: Entry point to the program. Initializes everything, and goes into a
//       message-processing loop. Idle time is used to render the scene.
//-----------------------------------------------------------------------------

INT WINAPI WinMain( HINSTANCE hInst, HINSTANCE, LPSTR, INT )
{
	theApp = new SynRad();
	if( !theApp->Create( 1024 , 768 , FALSE ) ) {
		char *logs = GLToolkit::GetLogs();
#ifdef WIN32
		if(logs) MessageBox(NULL,logs,"Synrad [Fatal error]",MB_OK);
#else
		if( logs ) {
			printf("Synrad [Fatal error]\n");
			printf(logs);
		}
#endif
		SAFE_FREE(logs);
		delete theApp;
		return -1;
	}
	try {
		theApp->Run();
	} catch(Error &e) {
		theApp->CrashHandler(&e);
	}
	delete theApp;
	return 0;
}



//-----------------------------------------------------------------------------
// Name: SynRad()
// Desc: Application constructor. Sets default attributes for the app.
//-----------------------------------------------------------------------------
SynRad::SynRad()
{
	//Get number of cores
	SYSTEM_INFO sysinfo;
	GetSystemInfo( &sysinfo );
	numCPU = (int) sysinfo.dwNumberOfProcessors;

	lastSaveTime=0.0f;
	lastSaveTimeSimu=0.0f;
	changedSinceSave=FALSE;
	nbDesStart=0;
	nbHitStart=0;
	lastUpdate=0.0;
	nbFormula = 0;
	nbRecent = 0;
	nbRecentPAR = 0;
	nbView = 0;
	nbSelection = 0;
	idView = 0;
	idSelection = 0;
	nbProc = numCPU;
	curViewer = 0;
	strcpy(currentDir,".");
	strcpy(currentSelDir,".");
	memset(formulas,0,sizeof formulas);
	formulaSettings = NULL;
	collapseSettings = NULL;
	moveVertex = NULL;
	scaleVertex = NULL;
	scaleFacet = NULL;
	regionInfo = NULL;
	selectDialog = NULL;
	moveFacet = NULL;
	exportDesorption = NULL;
	mirrorFacet = NULL;
	rotateFacet = NULL;
	alignFacet = NULL;
	addVertex = NULL;
	facetMesh = NULL;
	facetDetails = NULL;
	viewer3DSettings = NULL;
	textureSettings = NULL;
	globalSettings = NULL;
	facetCoordinates = NULL;
	vertexCoordinates = NULL;
	profilePlotter = NULL;
	spectrumPlotter = NULL;
	viewEditor = NULL;
	texturePlotter = NULL;
	//outgassingMap = NULL;
	m_strWindowTitle = APP_NAME;
	wnd->SetBackgroundColor(212,208,200);
	m_bResizable = TRUE;
	m_minScreenWidth  = 800;
	m_minScreenHeight = 600;
	tolerance=1e-8;
}

//-----------------------------------------------------------------------------
// Name: OneTimeSceneInit()
// Desc: Called during initial app startup, this function performs all the
//       permanent initialization.
//-----------------------------------------------------------------------------
int SynRad::OneTimeSceneInit()
{

	GLToolkit::SetIcon32x32("images/app_icon.png");

	for(int i=0;i<MAX_VIEWER;i++) {
		viewer[i] = new GeometryViewer(i);
		Add(viewer[i]);
	}
	modeSolo = TRUE;
	//nbSt = 0;

	LoadConfig();

	//Materials
	materialPaths.push_back("param\\copper.csv");
	materialPaths.push_back("param\\aluminium.csv");
	materialPaths.push_back("param\\al2o3.csv");

	menu = new GLMenuBar(0);
	wnd->SetMenuBar(menu);
	menu->Add("File");
	menu->GetSubMenu("File")->Add("&Load",MENU_FILE_LOAD,SDLK_o,CTRL_MODIFIER);
	menu->GetSubMenu("File")->Add("&Save",MENU_FILE_SAVE,SDLK_s,CTRL_MODIFIER);
	menu->GetSubMenu("File")->Add("Save as",MENU_FILE_SAVEAS);
	menu->GetSubMenu("File")->Add("Export DES file (deprecated)",MENU_FILE_EXPORT_DESORP);
	menu->GetSubMenu("File")->Add(NULL); //separator
	menu->GetSubMenu("File")->Add("&Insert geometry");
	menu->GetSubMenu("File")->GetSubMenu("Insert geometry")->Add("&To current structure",MENU_FILE_INSERTGEO);
	menu->GetSubMenu("File")->GetSubMenu("Insert geometry")->Add("&To new structure",MENU_FILE_INSERTGEO_NEWSTR);
	menu->GetSubMenu("File")->Add("&Export selected facets",MENU_FILE_EXPORTMESH);

	menu->GetSubMenu("File")->Add("Export selected textures");
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Element Area",MENU_FILE_EXPORTTEXTURE_AREA);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("MC hits",MENU_FILE_EXPORTTEXTURE_MCHITS);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Flux(ph/sec)/Beam current(mA)",MENU_FILE_EXPORTTEXTURE_FLUX);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Power(W)/Beam Current(mA)",MENU_FILE_EXPORTTEXTURE_POWER);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Flux/Area/Beam current (ph/sec/cm2/mA)",MENU_FILE_EXPORTTEXTURE_FLUXPERAREA);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Power/Area/Beam current (W/sec/cm2/mA)",MENU_FILE_EXPORTTEXTURE_POWERPERAREA);
	/*
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Average Velocity",MENU_FILE_EXPORTTEXTURE_AVG_VELOCITY);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Velocity",MENU_FILE_EXPORTTEXTURE_VELOCITY);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Count",MENU_FILE_EXPORTTEXTURE_COUNT);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Hits",MENU_FILE_EXPORTTEXTURE_HITS);
	
	
	TO DO FOR SYNRAD!
	*/

	menu->GetSubMenu("File")->Add(NULL); // Separator
	menu->GetSubMenu("File")->Add("Load recent");
	for(int i=nbRecent-1;i>=0;i--)
		menu->GetSubMenu("File")->GetSubMenu("Load recent")->Add(recents[i],MENU_FILE_LOADRECENT+i);
	menu->GetSubMenu("File")->Add(NULL); // Separator
	menu->GetSubMenu("File")->Add("E&xit",MENU_FILE_EXIT);

	menu->GetSubMenu("File")->SetIcon(MENU_FILE_LOAD,65,24);
	menu->GetSubMenu("File")->SetIcon(MENU_FILE_SAVE,83,24);
	menu->GetSubMenu("File")->SetIcon(MENU_FILE_SAVEAS,101,24);

	menu->Add("Regions");
	menu->GetSubMenu("Regions")->Add("New...",MENU_REGIONS_NEW);
	menu->GetSubMenu("Regions")->Add("Load...",MENU_REGIONS_LOADPAR);
	menu->GetSubMenu("Regions")->Add("Load recent",MENU_REGIONS_LOADRECENT);
	for(int i=nbRecentPAR-1;i>=0;i--)
		menu->GetSubMenu("Regions")->GetSubMenu("Load recent")->Add(recentPARs[i],MENU_REGIONS_LOADRECENT+i);
	menu->GetSubMenu("Regions")->Add("Remove all",MENU_REGIONS_CLEARALL);
	menu->GetSubMenu("Regions")->Add(NULL); // Separator
	menu->GetSubMenu("Regions")->Add("Load to");
	PARloadToMenu=menu->GetSubMenu("Regions")->GetSubMenu("Load to");
	menu->GetSubMenu("Regions")->Add("Remove");
	PARremoveMenu=menu->GetSubMenu("Regions")->GetSubMenu("Remove");
	menu->GetSubMenu("Regions")->Add(NULL); // Separator
	menu->GetSubMenu("Regions")->Add("Region Info ...",MENU_REGIONS_REGIONINFO);

	menu->Add("Selection");
	menu->GetSubMenu("Selection")->Add("Select All Facets",MENU_FACET_SELECTALL,SDLK_a,CTRL_MODIFIER);
	menu->GetSubMenu("Selection")->Add("Select by Facet Number...",MENU_SELECTION_SELECTFACETNUMBER,SDLK_n,ALT_MODIFIER);
	menu->GetSubMenu("Selection")->Add("Select Sticking",MENU_FACET_SELECTSTICK);
	menu->GetSubMenu("Selection")->Add("Select Transparent",MENU_FACET_SELECTTRANS);
	menu->GetSubMenu("Selection")->Add("Select 2 sided",MENU_FACET_SELECT2SIDE);
	menu->GetSubMenu("Selection")->Add("Select Texture",MENU_FACET_SELECTTEXT);
	menu->GetSubMenu("Selection")->Add("Select Profile",MENU_FACET_SELECTPROF);
	menu->GetSubMenu("Selection")->Add("Select Spectrum",MENU_FACET_SELECTSPECTRUM);
	menu->GetSubMenu("Selection")->Add(NULL); // Separator
	menu->GetSubMenu("Selection")->Add("Select Abs > 0",MENU_FACET_SELECTABS);
	menu->GetSubMenu("Selection")->Add("Select Hit > 0",MENU_FACET_SELECTHITS);
	menu->GetSubMenu("Selection")->Add(NULL); // Separator
	menu->GetSubMenu("Selection")->Add("Select teleport",MENU_FACET_SELECTTP);
	menu->GetSubMenu("Selection")->Add("Select link facet",MENU_FACET_SELECTDEST);
	menu->GetSubMenu("Selection")->Add("Select volatile facet",MENU_FACET_SELECTVOL);
	menu->GetSubMenu("Selection")->Add("Select non simple",MENU_FACET_SELECTERR);
	menu->GetSubMenu("Selection")->Add("Invert selection",MENU_FACET_INVERTSEL,SDLK_i, CTRL_MODIFIER);
	menu->GetSubMenu("Selection")->Add(NULL); // Separator 

	menu->GetSubMenu("Selection")->Add("Memorize selection to");
	memorizeSelectionsMenu = menu->GetSubMenu("Selection")->GetSubMenu("Memorize selection to");
	memorizeSelectionsMenu->Add("Add new...",MENU_SELECTION_ADDNEW,SDLK_w,CTRL_MODIFIER);
	memorizeSelectionsMenu->Add(NULL); // Separator

	menu->GetSubMenu("Selection")->Add("Select memorized");
	selectionsMenu = menu->GetSubMenu("Selection")->GetSubMenu("Select memorized");

	menu->GetSubMenu("Selection")->Add("Clear memorized",MENU_SELECTION_CLEARSELECTIONS);
	clearSelectionsMenu = menu->GetSubMenu("Selection")->GetSubMenu("Clear memorized");
	clearSelectionsMenu->Add("Clear All",MENU_SELECTION_CLEARALL);
	clearSelectionsMenu->Add(NULL); // Separator

	menu->GetSubMenu("Selection")->Add(NULL); // Separator
	menu->GetSubMenu("Selection")->Add(NULL); // Separator
	menu->GetSubMenu("Selection")->Add("Select all vertex",MENU_VERTEX_SELECTALL);
	menu->GetSubMenu("Selection")->Add("Unselect all vertex",MENU_VERTEX_UNSELECTALL);
	menu->GetSubMenu("Selection")->Add("Select Isolated vertex",MENU_VERTEX_SELECT_ISOLATED);
	menu->GetSubMenu("Selection")->Add("Select Coplanar vertex (visible on screen)",MENU_VERTEX_SELECT_COPLANAR);


	menu->Add("Tools");
	menu->GetSubMenu("Tools")->Add("3D Settings ..."   ,MENU_EDIT_3DSETTINGS,SDLK_b,CTRL_MODIFIER);
	menu->GetSubMenu("Tools")->Add("Texture scaling...",MENU_EDIT_TSCALING,SDLK_d,CTRL_MODIFIER);
	menu->GetSubMenu("Tools")->Add("Add formula ..."   ,MENU_EDIT_ADDFORMULA);
	menu->GetSubMenu("Tools")->Add("Global Settings ..."   ,MENU_EDIT_GLOBALSETTINGS);
	menu->GetSubMenu("Tools")->Add(NULL); // Separator
	menu->GetSubMenu("Tools")->Add("Texture Plotter ...",MENU_FACET_TEXPLOTTER,SDLK_t,ALT_MODIFIER);
	menu->GetSubMenu("Tools")->Add("Profile Plotter ...",MENU_FACET_PROFPLOTTER,SDLK_p,ALT_MODIFIER);
	menu->GetSubMenu("Tools")->Add("Spectrum Plotter ...",MENU_FACET_SPECTRUMPLOTTER);


	menu->GetSubMenu("Tools")->SetIcon(MENU_EDIT_3DSETTINGS,119,24);
	menu->GetSubMenu("Tools")->SetIcon(MENU_EDIT_TSCALING,137,24);
	menu->GetSubMenu("Tools")->SetIcon(MENU_EDIT_ADDFORMULA,155,24);

	menu->Add("Facet");
	menu->GetSubMenu("Facet")->Add("Collapse ...",MENU_FACET_COLLAPSE);
	menu->GetSubMenu("Facet")->Add("Swap normal",MENU_FACET_SWAPNORMAL,SDLK_n,CTRL_MODIFIER);
	menu->GetSubMenu("Facet")->Add("Shift vertex",MENU_FACET_SHIFTVERTEX,SDLK_h,CTRL_MODIFIER);
	menu->GetSubMenu("Facet")->Add("Edit coordinates ...",MENU_FACET_COORDINATES);
	menu->GetSubMenu("Facet")->Add("Move ...",MENU_FACET_MOVE);
	menu->GetSubMenu("Facet")->Add("Scale ...",MENU_FACET_SCALE);
	menu->GetSubMenu("Facet")->Add("Mirror ...",MENU_FACET_MIRROR);
	menu->GetSubMenu("Facet")->Add("Rotate ...",MENU_FACET_ROTATE);
	menu->GetSubMenu("Facet")->Add("Align ...",MENU_FACET_ALIGN);
	menu->GetSubMenu("Facet")->Add("Remove selected",MENU_FACET_REMOVESEL,SDLK_DELETE,CTRL_MODIFIER);
	menu->GetSubMenu("Facet")->Add("Explode selected",MENU_FACET_EXPLODE);
	menu->GetSubMenu("Facet")->Add("Create difference of 2",MENU_FACET_CREATE_DIFFERENCE);
	menu->GetSubMenu("Facet")->Add(NULL); // Separator

	menu->GetSubMenu("Facet")->Add("Facet Details ...",MENU_FACET_DETAILS);
	menu->GetSubMenu("Facet")->Add("Facet Mesh ...",MENU_FACET_MESH);

	facetMenu = menu->GetSubMenu("Facet");
	facetMenu->SetEnabled(MENU_FACET_MESH,FALSE);

	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_COLLAPSE,173,24);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_SWAPNORMAL,191,24);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_SHIFTVERTEX,90,77);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_COORDINATES,209,24);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_PROFPLOTTER,227,24);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_DETAILS,54,77);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_MESH,72,77);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_TEXPLOTTER,108,77);

	menu->Add("Vertex");
	menu->GetSubMenu("Vertex")->Add("Create Facet from Selected");
	menu->GetSubMenu("Vertex")->GetSubMenu("Create Facet from Selected")->Add("Convex Hull",MENU_VERTEX_CREATE_POLY_CONVEX,SDLK_v,ALT_MODIFIER);
	menu->GetSubMenu("Vertex")->GetSubMenu("Create Facet from Selected")->Add("Keep selection order",MENU_VERTEX_CREATE_POLY_ORDER);
	menu->GetSubMenu("Vertex")->Add("Remove selected",MENU_VERTEX_REMOVE);
	menu->GetSubMenu("Vertex")->Add("Vertex coordinates...",MENU_VERTEX_COORDINATES);
	menu->GetSubMenu("Vertex")->Add("Move selected...",MENU_VERTEX_MOVE);
	menu->GetSubMenu("Vertex")->Add("Scale selected...",MENU_VERTEX_SCALE);
	menu->GetSubMenu("Vertex")->Add("Add new...",MENU_VERTEX_ADD);

	menu->Add("View");

	menu->GetSubMenu("View")->Add("Structure",MENU_VIEW_STRUCTURE_P);
	structMenu = menu->GetSubMenu("View")->GetSubMenu("Structure");
	UpdateStructMenu();

	menu->GetSubMenu("View")->Add("Full Screen",MENU_VIEW_FULLSCREEN);

	menu->GetSubMenu("View")->Add(NULL); // Separator 

	menu->GetSubMenu("View")->Add("Memorize view to");
	memorizeViewsMenu = menu->GetSubMenu("View")->GetSubMenu("Memorize view to");
	memorizeViewsMenu->Add("Add new...",MENU_VIEW_ADDNEW,SDLK_q,CTRL_MODIFIER);
	memorizeViewsMenu->Add(NULL); // Separator

	menu->GetSubMenu("View")->Add("Select memorized");
	viewsMenu = menu->GetSubMenu("View")->GetSubMenu("Select memorized");

	menu->GetSubMenu("View")->Add("Clear memorized",MENU_VIEW_CLEARVIEWS);
	clearViewsMenu = menu->GetSubMenu("View")->GetSubMenu("Clear memorized");
	clearViewsMenu->Add("Clear All",MENU_VIEW_CLEARALL);

	menu->GetSubMenu("View")->SetIcon(MENU_VIEW_STRUCTURE_P,0,77);
	menu->GetSubMenu("View")->SetIcon(MENU_VIEW_FULLSCREEN,18,77);

	menu->Add("T&est");
	menu->GetSubMenu("Test")->Add("Pipe (L/R=0.0001)",MENU_TEST_PIPE0001);
	menu->GetSubMenu("Test")->Add("Pipe (L/R=1)",MENU_TEST_PIPE1);
	menu->GetSubMenu("Test")->Add("Pipe (L/R=10)",MENU_TEST_PIPE10);
	menu->GetSubMenu("Test")->Add("Pipe (L/R=100)",MENU_TEST_PIPE100);
	menu->GetSubMenu("Test")->Add("Pipe (L/R=1000)",MENU_TEST_PIPE1000);
	menu->GetSubMenu("Test")->Add("Pipe (L/R=10000)",MENU_TEST_PIPE10000);
	//Quick test pipe
	menu->GetSubMenu("Test")->Add(NULL);
	menu->GetSubMenu("Test")->Add("Quick Pipe",MENU_QUICKPIPE,SDLK_q,ALT_MODIFIER);

	geomNumber = new GLTextField(0,NULL);
	geomNumber->SetEditable(FALSE);
	Add(geomNumber);

	togglePanel = new GLTitledPanel("3D Viewer settings");
	togglePanel->SetClosable(TRUE);
	Add(togglePanel);

	showNormal = new GLToggle(0,"Normals");
	togglePanel->Add(showNormal);

	showRule = new GLToggle(0,"Rules");
	togglePanel->Add(showRule);

	showUV = new GLToggle(0,"\201,\202");
	togglePanel->Add(showUV);

	showLeak = new GLToggle(0,"Leaks");
	togglePanel->Add(showLeak);

	showHit = new GLToggle(0,"Hits");
	togglePanel->Add(showHit);

	showLine = new GLToggle(0,"Lines");
	togglePanel->Add(showLine);

	showVolume = new GLToggle(0,"Volume");
	togglePanel->Add(showVolume);

	showTexture = new GLToggle(0,"Texture");
	togglePanel->Add(showTexture);

	showFilter = new GLToggle(0,"Filtering");
	togglePanel->Add(showFilter);

	showIndex = new GLToggle(0,"Indices");
	togglePanel->Add(showIndex);

	showVertex = new GLToggle(0,"Vertices");
	togglePanel->Add(showVertex);

	showMoreBtn = new GLButton(0,"More ...");
	togglePanel->Add(showMoreBtn);

	simuPanel = new GLTitledPanel("Simulation");
	simuPanel->SetClosable(TRUE);
	Add(simuPanel);

	startSimu = new GLButton(0,"Start/Stop");
	simuPanel->Add(startSimu);

	resetSimu = new GLButton(0,"Reset");
	simuPanel->Add(resetSimu);

	modeLabel = new GLLabel("Mode");
	simuPanel->Add(modeLabel);

	modeCombo = new GLCombo(0);
	modeCombo->SetEditable(TRUE);
	modeCombo->SetSize(2);
	modeCombo->SetValueAt(0,"Fluxwise");
	modeCombo->SetValueAt(1,"Powerwise");
	modeCombo->SetSelectedIndex(1);
	simuPanel->Add(modeCombo);

	singleACBtn = new GLButton(0,"1");
	singleACBtn->SetEnabled(FALSE);
	//simuPanel->Add(singleACBtn);

	hitLabel = new GLLabel("Hits");
	simuPanel->Add(hitLabel);

	hitNumber = new GLTextField(0,NULL);
	hitNumber->SetEditable(FALSE);
	simuPanel->Add(hitNumber);

	desLabel = new GLLabel("Gen.");
	simuPanel->Add(desLabel);

	desNumber = new GLTextField(0,NULL);
	desNumber->SetEditable(FALSE);
	simuPanel->Add(desNumber);

	leakLabel = new GLLabel("Leaks");
	simuPanel->Add(leakLabel);

	leakNumber = new GLTextField(0,NULL);
	leakNumber->SetEditable(FALSE);
	simuPanel->Add(leakNumber);

	doseLabel = new GLLabel("Dose");
	simuPanel->Add(doseLabel);

	doseNumber = new GLTextField(0,NULL);
	doseNumber->SetEditable(FALSE);
	simuPanel->Add(doseNumber);

	sTimeLabel = new GLLabel("Time");
	simuPanel->Add(sTimeLabel);

	sTime = new GLTextField(0,NULL);
	sTime->SetEditable(FALSE);
	simuPanel->Add(sTime);

	facetPanel = new GLTitledPanel("Selected Facet");
	facetPanel->SetClosable(TRUE);
	Add(facetPanel);


	facetStickingLabel = new GLLabel("Sticking factor:");
	facetPanel->Add(facetStickingLabel);
	facetSticking = new GLTextField(0,NULL);
	facetPanel->Add(facetSticking);

	facetRoughnessLabel = new GLLabel("Material roughness:");
	facetRoughnessLabel->SetVisible(FALSE);
	facetPanel->Add(facetRoughnessLabel);
	facetRoughness = new GLTextField(0,NULL);
	facetRoughness->SetVisible(FALSE);
	facetPanel->Add(facetRoughness);

	facetSideLabel = new GLLabel("Sides:");
	facetPanel->Add(facetSideLabel);

	facetSideType = new GLCombo(0);
	facetSideType->SetSize(2);
	facetSideType->SetValueAt(0,"1 Sided");
	facetSideType->SetValueAt(1,"2 Sided");
	facetPanel->Add(facetSideType);

	facetTLabel = new GLLabel("Opacity:");
	facetPanel->Add(facetTLabel);
	facetOpacity = new GLTextField(0,NULL);
	facetPanel->Add(facetOpacity);

	facetRLabel = new GLLabel("Reflection:");
	facetPanel->Add(facetRLabel);
	facetReflType = new GLCombo(0);
	facetReflType->SetSize(5); //TODO
	facetReflType->SetValueAt(0,"Diffuse",REF_DIFFUSE);
	facetReflType->SetValueAt(1,"Mirror",REF_MIRROR);
	facetReflType->SetValueAt(2,"Copper",REF_MATERIAL+0);
	facetReflType->SetValueAt(3,"Aluminium",REF_MATERIAL+1);
	facetReflType->SetValueAt(4,"Al2O3",REF_MATERIAL+2);
	facetPanel->Add(facetReflType);

	facetAreaLabel = new GLLabel("Area (cm\262):");
	facetPanel->Add(facetAreaLabel);
	facetArea = new GLTextField(0,NULL);
	facetPanel->Add(facetArea);

	facetTPLabel = new GLLabel("Teleport to facet  #");
	facetPanel->Add(facetTPLabel);
	facetTeleport = new GLTextField(0,NULL);
	facetPanel->Add(facetTeleport);

	facetLinkLabel = new GLLabel("Link:");
	facetPanel->Add(facetLinkLabel);
	facetSILabel = new GLTextField(0,"");
	facetSILabel->SetEditable(TRUE);

	facetPanel->Add(facetSILabel);

	facetStrLabel = new GLLabel("Structure:");
	facetPanel->Add(facetStrLabel);
	facetSuperDest = new GLTextField(0,NULL);
	facetPanel->Add(facetSuperDest);

	facetReLabel = new GLLabel("Profile:");
	facetPanel->Add(facetReLabel);
	facetRecType = new GLCombo(0);
	facetRecType->SetSize(4);
	facetRecType->SetValueAt(0,"None");
	facetRecType->SetValueAt(1,"Along \201");
	facetRecType->SetValueAt(2,"Along \202");
	facetRecType->SetValueAt(3,"Angular");
	facetPanel->Add(facetRecType);

	facetSpectrumLabel = new GLLabel("Record spectrum:");
	facetPanel->Add(facetSpectrumLabel);
	facetSpectrumCombo = new GLCombo(0);
	facetSpectrumCombo->SetSize(2);
	facetSpectrumCombo->SetValueAt(0,"Off");
	facetSpectrumCombo->SetValueAt(1,"On");
	facetPanel->Add(facetSpectrumCombo);

	facetTexBtn = new GLButton(0,"Mesh...");
	facetTexBtn->SetEnabled(FALSE);
	facetPanel->Add(facetTexBtn);


	facetMoreBtn = new GLButton(0,"Details...");
	facetPanel->Add(facetMoreBtn);

	facetCoordBtn = new GLButton(0,"Coord...");
	facetPanel->Add(facetCoordBtn);

	facetApplyBtn = new GLButton(0,"Apply");
	facetApplyBtn->SetEnabled(FALSE);
	facetPanel->Add(facetApplyBtn);

	facetList = new GLList(0);
	facetList->SetWorker(&worker);
	facetList->SetGrid(TRUE);
	facetList->SetSelectionMode(MULTIPLE_ROW);
	facetList->SetSize(5,1);
	facetList->SetColumnWidths((int*)cWidth);
	facetList->SetColumnLabels((char **)cName);
	facetList->SetColumnLabelVisible(TRUE);
	facetList->Sortable = TRUE;
	Add(facetList);

	ClearFacetParams();
	UpdateViewerParams();
	PlaceComponents();

	//LoadFile();
	try {
		for (int i=0;i<(int)materialPaths.size();i++) //load materials
			worker.AddMaterial(&(materialPaths[i]));
	} catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"Failed to load material library\n%s",e.GetMsg());
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
	}

	try {
		worker.SetProcNumber(nbProc);

	} catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"Failed to start working sub-process(es), simulation not available\n%s",e.GetMsg());
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
	}

	//AnimateViewerChange(0,TRUE);

	PlaceComponents();
	SelectViewer(0);

	viewer[0]->Paint();

	if (checkForUpdates) {
		//Launch updater tool
		char command[2048];
		char CWD [MAX_PATH];
		_getcwd( CWD, MAX_PATH );
		//sprintf(tmp5,"%s\\synrad_updater_tmp.exe",CWD);
		if(FileUtils::Exist("synrad_updater_tmp.exe")) { //rename after new installation
			sprintf(command,"move \"%s\\synrad_updater_tmp.exe\" \"%s\\synrad_updater.exe\"",CWD,CWD);
			system(command);
		}
		//win_sparkle_set_appcast_url("http://test.xml");
		//win_sparkle_init();
		//win_sparkle_check_update_with_ui();
		if(FileUtils::Exist("synrad_updater.exe"))
			StartProc_background("synrad_updater.exe");
		else GLMessageBox::Display("synrad_updater.exe not found. You will not receive updates to Synrad."
			"\n(You can disable checking for updates in Tools/Global Settings)","Updater module missing.",GLDLG_OK,GLDLG_ICONINFO);
	}
	return GL_OK;
}

//-----------------------------------------------------------------------------
// Name: Resize()
// Desc: Called when the window is resized
//-----------------------------------------------------------------------------
int SynRad::Resize(DWORD width, DWORD height, BOOL forceWindowed) {
	int r = GLApplication::Resize(width,height,forceWindowed);
	PlaceComponents();
	return r;
}

//-----------------------------------------------------------------------------
// Name: PlaceComponents()
// Desc: Place components on screen
//-----------------------------------------------------------------------------
void SynRad::PlaceComponents() {

	int sx = m_screenWidth-205;  
	int sy = 30;

	Place3DViewer();

	// ---------------------------------------------------------
	geomNumber->SetBounds(sx,3,202,18);

	// Viewer settings ----------------------------------------
	togglePanel->SetBounds(sx,sy,202,112);

	togglePanel->SetCompBounds(showRule,5,20,60,18);
	togglePanel->SetCompBounds(showNormal,70,20,60,18);
	togglePanel->SetCompBounds(showUV,135,20,60,18);

	togglePanel->SetCompBounds(showLine,5,42,60,18);
	togglePanel->SetCompBounds(showLeak,70,42,60,18);
	togglePanel->SetCompBounds(showHit,135,42,60,18);

	togglePanel->SetCompBounds(showVolume,5,64,60,18);
	togglePanel->SetCompBounds(showTexture,70,64,60,18);
	togglePanel->SetCompBounds(showFilter,135,64,60,18);

	togglePanel->SetCompBounds(showVertex,5,86,60,18);
	togglePanel->SetCompBounds(showIndex,70,86,60,18);
	togglePanel->SetCompBounds(showMoreBtn,137,86,55,19);

	sy += (togglePanel->GetHeight() + 5);

	// Selected facet -----------------------------------------
	facetPanel->SetBounds(sx,sy,202,265);

	facetPanel->SetCompBounds(facetStickingLabel     ,7  ,15 ,100 ,18);
	facetPanel->SetCompBounds(facetSticking   ,110 ,15 ,82,18);

	facetPanel->SetCompBounds(facetRoughnessLabel     ,7  ,15 ,100 ,18);
	facetPanel->SetCompBounds(facetRoughness   ,110 ,15 ,82,18);

	facetPanel->SetCompBounds(facetRLabel     ,7  ,40,60 ,18);
	facetPanel->SetCompBounds(facetReflType   ,65 ,40,130,18);		

	facetPanel->SetCompBounds(facetSideLabel,7,65,50,18);
	facetPanel->SetCompBounds(facetSideType   ,65,65 ,130 ,18);

	facetPanel->SetCompBounds(facetTLabel     ,7  ,90 ,100 ,18);
	facetPanel->SetCompBounds(facetOpacity    ,110 ,90 ,82 ,18);

	facetPanel->SetCompBounds(facetAreaLabel     ,7,115,100 ,18);
	facetPanel->SetCompBounds(facetArea       ,110,115,82 ,18);

	facetPanel->SetCompBounds(facetTPLabel     ,7,140,100,18);
	facetPanel->SetCompBounds(facetTeleport   ,110,140,82,18);

	facetPanel->SetCompBounds(facetStrLabel   ,7  ,165,55 ,18); //Structure:
	facetPanel->SetCompBounds(facetSILabel    ,65 ,165,42 ,18); //Editable Textfield
	facetPanel->SetCompBounds(facetLinkLabel  ,115,165,18 ,18); //Link
	facetPanel->SetCompBounds(facetSuperDest  ,148,165,42 ,18); //Textfield

	facetPanel->SetCompBounds(facetReLabel    ,7  ,190,60 ,18);
	facetPanel->SetCompBounds(facetRecType    ,65 ,190,130,18);

	facetPanel->SetCompBounds(facetSpectrumLabel    ,7  ,215,70 ,18);
	facetPanel->SetCompBounds(facetSpectrumCombo    ,100 ,215,50,18);

	facetPanel->SetCompBounds(facetMoreBtn    ,5  ,240,45 ,18);
	facetPanel->SetCompBounds(facetCoordBtn   ,53 ,240,44 ,18);
	facetPanel->SetCompBounds(facetTexBtn     ,101,240,50 ,18);
	facetPanel->SetCompBounds(facetApplyBtn   ,155,240,40 ,18);


	sy += (facetPanel->GetHeight() + 5);

	// Simulation ---------------------------------------------
	simuPanel->SetBounds(sx,sy,202,194);

	simuPanel->SetCompBounds(startSimu,5,20,92,19);
	simuPanel->SetCompBounds(resetSimu,102,20,93,19);
	simuPanel->SetCompBounds(modeLabel,5,45,30,18);
	simuPanel->SetCompBounds(modeCombo,40,45,85,18);
	simuPanel->SetCompBounds(hitLabel,5,70,30,18);
	simuPanel->SetCompBounds(hitNumber,40,70,155,18);
	simuPanel->SetCompBounds(desLabel,5,95,30,18);
	simuPanel->SetCompBounds(desNumber,40,95,155,18);
	simuPanel->SetCompBounds(leakLabel,5,120,30,18);
	simuPanel->SetCompBounds(leakNumber,40,120,155,18);
	simuPanel->SetCompBounds(doseLabel,5,145,30,18);
	simuPanel->SetCompBounds(doseNumber,40,145,155,18);
	simuPanel->SetCompBounds(sTimeLabel,5,170,30,18);
	simuPanel->SetCompBounds(sTime,40,170,155,18);

	sy += (simuPanel->GetHeight() + 5);

	// ---------------------------------------------------------
	int lg = m_screenHeight - (nbFormula*25 + 23);

	facetList->SetBounds(sx,sy,202,lg-sy);

	// ---------------------------------------------------------

	for(int i=0;i<nbFormula;i++) {
		formulas[i].name->SetBounds(sx,lg+5,95,18);
		formulas[i].value->SetBounds(sx+90,lg+5,87,18);
		formulas[i].setBtn->SetBounds(sx+182,lg+5,20,18);
		lg+=25;
	}

}

void SynRad::Place3DViewer() {

	int sx = m_screenWidth-205;

	// 3D Viewer ----------------------------------------------
	int fWidth = m_screenWidth-215;
	int fHeight = m_screenHeight-27;
	int Width2 = fWidth/2-1;
	int Height2 = fHeight/2-1;

	if( modeSolo ) {
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->SetVisible(FALSE);
		viewer[curViewer]->SetBounds(3,3,fWidth,fHeight);
		viewer[curViewer]->SetVisible(TRUE);
	} else {
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->SetVisible(TRUE);
		viewer[0]->SetBounds(3       ,3        ,Width2,Height2);
		viewer[1]->SetBounds(6+Width2,3        ,Width2,Height2);
		viewer[2]->SetBounds(3       ,6+Height2,Width2,Height2);
		viewer[3]->SetBounds(6+Width2,6+Height2,Width2,Height2);
	}

}

void SynRad::UpdateViewers() {
	for(int i=0;i<MAX_VIEWER;i++)
		viewer[i]->UpdateMatrix();
}

void SynRad::SetFacetSearchPrg(BOOL visible,char *text) {
	for(int i=0;i<MAX_VIEWER;i++) {
		viewer[i]->facetSearchState->SetVisible(visible);
		viewer[i]->facetSearchState->SetText(text);
	}
	GLWindowManager::Repaint();
}

void SynRad::UpdateViewerParams() {

	showNormal->SetCheck(viewer[curViewer]->showNormal);
	showRule->SetCheck(viewer[curViewer]->showRule);
	showUV->SetCheck(viewer[curViewer]->showUV);
	showLeak->SetCheck(viewer[curViewer]->showLeak);
	showHit->SetCheck(viewer[curViewer]->showHit);
	showVolume->SetCheck(viewer[curViewer]->showVolume);
	showLine->SetCheck(viewer[curViewer]->showLine);
	showTexture->SetCheck(viewer[curViewer]->showTexture);
	showFilter->SetCheck(viewer[curViewer]->showFilter);
	showVertex->SetCheck(viewer[curViewer]->showVertex);
	showIndex->SetCheck(viewer[curViewer]->showIndex);

	// Force all views to have the same showColormap
	viewer[1]->showColormap = viewer[0]->showColormap;
	viewer[2]->showColormap = viewer[0]->showColormap;
	viewer[3]->showColormap = viewer[0]->showColormap;
	worker.GetGeometry()->texColormap = viewer[0]->showColormap;

}

//-----------------------------------------------------------------------------
// Name: ClearFacetParams()
// Desc: Reset selected facet parameters.
//-----------------------------------------------------------------------------

void SynRad::ClearFacetParams() 
{
	facetPanel->SetTitle("Selected Facet (none)");
	facetSticking->Clear();
	facetSticking->SetEditable(FALSE);
	facetRoughness->Clear();
	facetRoughnessLabel->SetVisible(FALSE);
	facetRoughness->SetVisible(FALSE);
	facetTeleport->Clear();
	facetTeleport->SetEditable(FALSE);
	facetArea->SetEditable(FALSE);
	facetArea->Clear();
	facetSuperDest->Clear();
	facetSuperDest->SetEditable(FALSE);
	facetSILabel->Clear();
	facetSILabel->SetEditable(FALSE);
	facetOpacity->Clear();
	facetOpacity->SetEditable(FALSE);
	facetSideType->SetSelectedValue("");
	facetSideType->SetEditable(FALSE);
	facetReflType->SetSelectedValue("");
	facetReflType->SetEditable(FALSE);
	facetRecType->SetSelectedValue("");
	facetRecType->SetEditable(FALSE); 
	facetSpectrumCombo->SetSelectedValue("");
	facetSpectrumCombo->SetEditable(FALSE);
}

//-----------------------------------------------------------------------------
// Name: ApplyFacetParams()
// Desc: Apply facet parameters.
//-----------------------------------------------------------------------------

void SynRad::ApplyFacetParams() {
	if (!AskToReset()) return;
	changedSinceSave = TRUE;
	Geometry *geom = worker.GetGeometry();
	int nbFacet = geom->GetNbFacet();

	// Sticking
	double sticking;
	BOOL doSticking = FALSE;
	if( facetSticking->GetNumber(&sticking) ) {
		if( sticking<0.0 || sticking>1.0 ) {
			GLMessageBox::Display("Sticking must be in the range [0,1]","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
		doSticking = TRUE;
	} else {
		if( strcmp(facetSticking->GetText(),"..." )==0 ) doSticking = FALSE;
		else {
			GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
	}

	// Roughness
	double roughness;
	BOOL doRoughness = FALSE;
	if( facetRoughness->GetNumber(&roughness) ) {
		if( sticking<0.0) {
			GLMessageBox::Display("Roughness must be non-negative","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
		doRoughness = TRUE;
	} else {
		if( strcmp(facetRoughness->GetText(),"..." )==0 ) doRoughness = FALSE;
		else {
			GLMessageBox::Display("Invalid roughness number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
	}

	// teleport
	int teleport;
	BOOL doTeleport = FALSE;

	if( facetTeleport->GetNumberInt(&teleport) ) {
		if( teleport<0 || teleport>nbFacet ) {
			GLMessageBox::Display("Invalid teleport destination\n(If no teleport: set number to 0)","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		} else if (teleport>0 && geom->GetFacet(teleport-1)->selected) {
			char tmp[256];
			sprintf(tmp,"The teleport destination of facet #%d can't be itself!",teleport);
			GLMessageBox::Display(tmp,"Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
		doTeleport = TRUE;
	} else {
		if( strcmp(facetTeleport->GetText(),"..." )==0 ) doTeleport = FALSE;
		else {
			GLMessageBox::Display("Invalid teleport destination\n(If no teleport: set number to 0)","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
	}


	// opacity
	double opacity;
	BOOL doOpacity = FALSE;
	if( facetOpacity->GetNumber(&opacity) ) {
		if( opacity<0.0 || opacity>1.0 ) {
			GLMessageBox::Display("Opacity must be in the range [0,1]","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
		doOpacity = TRUE;
	} else {
		if( strcmp(facetOpacity->GetText(),"..." )==0 ) doOpacity = FALSE;
		else {
			GLMessageBox::Display("Invalid opacity number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
	}


	// Superstructure
	int superStruct;
	BOOL doSuperStruct = FALSE;
	if( sscanf(facetSILabel->GetText(),"%d",&superStruct)>0 && superStruct>0 && superStruct<=geom->GetNbStructure()) doSuperStruct = TRUE;
	else {
		if( strcmp(facetSILabel->GetText(),"..." )==0 ) doSuperStruct = FALSE;
		else{
			GLMessageBox::Display("Invalid superstructre number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
	}


	// Super structure destination
	int superDest;
	BOOL doSuper = FALSE;
	if( strcmp(facetSuperDest->GetText(),"none")==0 || strcmp(facetSuperDest->GetText(),"no")==0 || strcmp(facetSuperDest->GetText(),"0")==0 ) {
		doSuper = TRUE;
		superDest=0;
	} else if( sscanf(facetSuperDest->GetText(),"%d",&superDest)>0) {
		if (superDest==superStruct) {
			GLMessageBox::Display("Link and superstructure can't be the same","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		} else if (superDest>0 && superDest<=geom->GetNbStructure()) doSuper = TRUE;
	} else {
		if( strcmp(facetSuperDest->GetText(),"..." )==0 ) doSuper = FALSE;
		GLMessageBox::Display("Invalid superstructre destination","Error",GLDLG_OK,GLDLG_ICONERROR);
		UpdateFacetParams();
		return;
	}


	// Reflection type
	int reflType = facetReflType->GetSelectedIndex();

	// Record type
	int rType = facetRecType->GetSelectedIndex();

	// Spectrum recording
	int specType = facetSpectrumCombo->GetSelectedIndex();

	// 2sided
	int is2Sided = facetSideType->GetSelectedIndex();

	BOOL structChanged=FALSE; //if a facet gets into a new structure, we have to re-render the geometry
	// Update facets (local)
	for(int i=0;i<nbFacet;i++) {
		Facet *f = geom->GetFacet(i);
		if( f->selected ) {
			if(doSticking) f->sh.sticking = sticking;
			if(doRoughness) f->sh.roughness = roughness;
			if(doTeleport) f->sh.teleportDest = teleport;
			if(doOpacity) f->sh.opacity = opacity;
			if(reflType>=0) f->sh.reflectType = reflType;
			if(rType>=0) {
				f->sh.profileType = rType;
				f->sh.isProfile = (rType!=REC_NONE);
				
			} if (profilePlotter) profilePlotter->Refresh();
			if(specType>=0) {
				f->sh.hasSpectrum = specType;
				
			} if (spectrumPlotter) spectrumPlotter->Refresh();
			if(is2Sided>=0) f->sh.is2sided = is2Sided;
			if(doSuperStruct) {
				if (f->sh.superIdx != (superStruct-1)) {
					f->sh.superIdx = superStruct-1;
					structChanged=TRUE;
				}
			}
			if(doSuper) {
				f->sh.superDest = superDest;
				if(superDest) f->sh.opacity = 1.0; // Force opacity for link facet
			}
			f->UpdateFlags();
		}
	}
	if (structChanged) geom->RebuildLists();
	// Send to sub process
	try { worker.Reload(); } catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	UpdateFacetParams();
}

void SynRad::UpdateFacetlistSelected() {
	int nbSelected=0;
	Geometry *geom = worker.GetGeometry();
	int nbFacet=geom->GetNbFacet();
	int* selection=(int*)malloc(nbFacet*sizeof(int));
	for (int i=0;i<nbFacet;i++) {
		if (geom->GetFacet(i)->selected) {
			selection[nbSelected]=i;
			nbSelected++;
		}
	}
	if (nbSelected>1000) { //Index lookup would be very slow for 1000+ selected lines
		facetList->ReOrder();
		facetList->SetSelectedRows(selection,nbSelected,FALSE);
	} else {
		facetList->SetSelectedRows(selection,nbSelected,TRUE);
	}
	SAFE_FREE(selection);
}

//-----------------------------------------------------------------------------
// Name: UpdateFacetParams()
// Desc: Update selected facet parameters.
//-----------------------------------------------------------------------------

void SynRad::UpdateFacetParams(BOOL updateSelection) {

	char tmp[256];
	int sel0 = -1;

	// Update params
	Geometry *geom = worker.GetGeometry();
	int nbSel = geom->GetNbSelected();
	if(nbSel>0) {

		Facet *f0;
		Facet *f;

		// Get list of selected facet
		int *selection = (int *)malloc(nbSel*sizeof(int));
		int count=0;
		for(int i=0;i<geom->GetNbFacet();i++)
			if( geom->GetFacet(i)->selected )
				selection[count++] = i;

		sel0 = selection[0];
		f0 = geom->GetFacet(selection[0]);

		double sticking = f0->sh.sticking;
		double roughness = f0->sh.roughness;
		int    teleport = f0->sh.teleportDest;
		double opacity = f0->sh.opacity;
		double area = f0->sh.area;
		int    superDest = f0->sh.superDest;
		int    superIdx = f0->sh.superIdx;
		int    reflectType = f0->sh.reflectType;
		int    recType = f0->sh.profileType;
		int    is2sided = f0->sh.is2sided;
		int    hasSpectrum = f0->sh.hasSpectrum;

		BOOL stickingE = TRUE;
		BOOL roughnessE = TRUE;
		BOOL teleportE = TRUE;
		BOOL opacityE = TRUE;
		BOOL superDestE = TRUE;
		BOOL superIdxE = TRUE;
		BOOL reflectTypeE = TRUE;
		BOOL recordE = TRUE;
		BOOL is2sidedE = TRUE;
		BOOL hasSpectrumE = TRUE;
		BOOL isMaterialE = TRUE;

		for(int i=1;i<count;i++) {
			f = geom->GetFacet(selection[i]);
			stickingE = stickingE && (abs(f0->sh.sticking - f->sh.sticking)<1e-7);
			roughnessE = roughnessE && (abs(f0->sh.roughness - f->sh.roughness)<1e-7);
			isMaterialE = isMaterialE && ((f0->sh.reflectType<2) == (f->sh.reflectType<2));
			teleportE = teleportE && (f0->sh.teleportDest == f->sh.teleportDest);
			opacityE = opacityE && (abs(f0->sh.opacity - f->sh.opacity)<1e-7);
			superDestE = superDestE && (f0->sh.superDest == f->sh.superDest);
			superIdxE = superIdxE && (f0->sh.superIdx == f->sh.superIdx);
			is2sidedE = is2sidedE && (f0->sh.is2sided == f->sh.is2sided);
			reflectTypeE = reflectTypeE && (f0->sh.reflectType == f->sh.reflectType);
			recordE = recordE && (f0->sh.profileType == f->sh.profileType);
			hasSpectrumE = hasSpectrumE && (f0->sh.hasSpectrum == f->sh.hasSpectrum);
			if (f->sh.area>0) area+=f->sh.area;
		}

		if( nbSel==1 )
			sprintf(tmp,"Selected Facet (#%d)",selection[0]+1);
		else
			sprintf(tmp,"Selected Facet (%d selected)",count);

		// Old STR compatibility
		if( stickingE && f0->sh.superDest ) stickingE = FALSE;

		facetPanel->SetTitle(tmp);
		if (count>1) facetAreaLabel->SetText("Sum Area (cm\262):");
		else facetAreaLabel->SetText("Area (cm\262):");
		sprintf(tmp,"%g",area);
		facetArea->SetText(tmp);
		if(stickingE && isMaterialE) SetParam(facetSticking,sticking); else facetSticking->SetText("...");
		if(roughnessE && isMaterialE) SetParam(facetRoughness,roughness); else facetRoughness->SetText("...");
		if(teleportE) SetParam(facetTeleport,teleport); else facetTeleport->SetText("...");
		if(opacityE) SetParam(facetOpacity,opacity); else facetOpacity->SetText("...");
		if(is2sidedE) facetSideType->SetSelectedIndex(f0->sh.is2sided); else facetSideType->SetSelectedValue("...");
		if(reflectTypeE) facetReflType->SetSelectedIndex(f0->sh.reflectType); else facetReflType->SetSelectedValue("...");
		if(recordE) facetRecType->SetSelectedIndex(f0->sh.profileType); else facetRecType->SetSelectedValue("...");
		if(hasSpectrumE) facetSpectrumCombo->SetSelectedIndex(f0->sh.hasSpectrum); else facetSpectrumCombo->SetSelectedValue("...");
		if(superDestE) {
			if( f0->sh.superDest==0 ) {
				facetSuperDest->SetText("no");
			} else {
				sprintf(tmp,"%d",f0->sh.superDest);
				facetSuperDest->SetText(tmp);
			}
		} else {
			facetSuperDest->SetText("...");
		}
		if(superIdxE) {
			sprintf(tmp,"%d",f0->sh.superIdx+1);
			facetSILabel->SetText(tmp);
		} else {
			facetSILabel->SetText("...");
		}
		if( updateSelection ) {
			if (nbSel>1000 || geom->GetNbFacet()>50000) { //If it would take too much time to look up every selected facet in the list
				facetList->ReOrder();
				facetList->SetSelectedRows(selection,nbSel,FALSE);
			} else {
				facetList->SetSelectedRows(selection,nbSel,TRUE);
			}
			facetList->lastRowSel=-1;
		}

		free(selection);

		facetSticking->SetEditable(isMaterialE);
		facetRoughness->SetEditable(isMaterialE);
		BOOL isMaterial=facetReflType->GetSelectedIndex()>=2;
		facetSticking->SetVisible(!isMaterial);
		facetStickingLabel->SetVisible(!isMaterial);
		facetRoughness->SetVisible(isMaterial);
		facetRoughnessLabel->SetVisible(isMaterial);
		facetTeleport->SetEditable(TRUE);
		facetOpacity->SetEditable(TRUE);
		facetSuperDest->SetEditable(TRUE);
		facetSILabel->SetEditable(TRUE);
		facetSideType->SetEditable(TRUE);
		facetReflType->SetEditable(TRUE);
		facetRecType->SetEditable(TRUE);
		facetSpectrumCombo->SetEditable(TRUE);

		facetApplyBtn->SetEnabled(FALSE);
		facetMenu->SetEnabled(MENU_FACET_MESH,TRUE);
		facetTexBtn->SetEnabled(TRUE);

	} else {
		ClearFacetParams();
		facetMenu->SetEnabled(MENU_FACET_MESH,FALSE);
		facetTexBtn->SetEnabled(FALSE);
		if( updateSelection ) facetList->ClearSelection();
	}

	if( facetDetails ) facetDetails->Update();
	if( facetCoordinates ) facetCoordinates->UpdateFromSelection();
	if( texturePlotter ) texturePlotter->Update(m_fTime,TRUE);
	//if( outgassingMap ) outgassingMap->Update(m_fTime,TRUE);
}
//-----------------------------------------------------------------------------
// Name: Evaluate formula
// Desc: Update formula
//-----------------------------------------------------------------------------

int getVariable(char *name,char *preffix) {

	char tmp[256];
	int  idx;
	int lgthP = (int)strlen(preffix);
	int lgthN = (int)strlen(name);

	if( lgthP>=lgthN ) {
		return -1;
	} else {
		strcpy(tmp,name);
		tmp[lgthP]=0;
		if( _stricmp(tmp,preffix)==0 ) {
			strcpy(tmp,name+lgthP);
			int conv = sscanf(tmp,"%d",&idx);
			if(conv) {
				return idx;
			} else {
				return -1;
			}
		}
	}
	return -1;

}

void SynRad::UpdateFormula() {

	char tmp[256];

	int idx;
	Geometry *geom = worker.GetGeometry();
	int nbFacet = geom->GetNbFacet();

	for(int i=0;i<nbFormula;i++) {

		GLParser *f = formulas[i].parser;

		// Variables
		int nbVar = f->GetNbVariable();
		BOOL ok = TRUE;
		BOOL iName = FALSE;
		for(int j=0;j<nbVar && ok;j++) {

			VLIST *v = f->GetVariableAt(j);
			if( (idx=getVariable(v->name,"A"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = (double)geom->GetFacet(idx-1)->sh.counter.nbAbsorbed;
			} else if( (idx=getVariable(v->name,"D"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = (double)geom->GetFacet(idx-1)->sh.counter.nbDesorbed;
			} else if( (idx=getVariable(v->name,"H"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = (double)geom->GetFacet(idx-1)->sh.counter.nbHit;
			} else if( (idx=getVariable(v->name,"_A"))>0 ) {
				ok = (idx<=nbFacet);
			} else if( (idx=getVariable(v->name,"AR"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = geom->GetFacet(idx-1)->sh.area;
			} else if( _stricmp(v->name,"SUMDES")==0 ) {
				v->value = (double)worker.nbDesorption;
			} else if( _stricmp(v->name,"SUMABS")==0 ) {
				v->value = (double)worker.nbAbsorption;
			} else if( _stricmp(v->name,"SUMHIT")==0 ) {
				v->value = (double)worker.nbHit;
			} else if( _stricmp(v->name,"MPP")==0 ) {
				v->value = worker.distTraveledTotal/(double)worker.nbDesorption;
			} else if( _stricmp(v->name,"MFP")==0 ) {
				v->value = worker.distTraveledTotal/(double)worker.nbHit;
			}  else if( _stricmp(v->name,"ABSAR")==0 ) {
				double sumArea = 0.0;
				for(int i=0;i<geom->GetNbFacet();i++) {
					Facet *f = geom->GetFacet(i);
					if(f->sh.sticking>0.0) sumArea += f->sh.area*f->sh.opacity*(f->sh.is2sided?2.0:1.0);
				}
				v->value = sumArea;
			} else if( _stricmp(v->name,"KB")==0 ) {
				v->value = 1.3806504e-23;
			} else if( _stricmp(v->name,"R")==0 ) {
				v->value = 8.314472;
			} else if( _stricmp(v->name,"Na")==0 ) {
				v->value = 6.02214179e23;
			} else {
				formulas[i].value->SetText("Invalid var name");
				ok = FALSE;
				iName = TRUE;
			}

			if(!ok && !iName) formulas[i].value->SetText("Invalid var index");

		}

		// Evaluation
		if( ok ) {
			double r;
			if( f->Evaluate(&r) ) {
				sprintf(tmp,"%g",r);
				formulas[i].value->SetText(tmp);
			} else {
				formulas[i].value->SetText(f->GetErrorMsg());
			}
		}
	}
}

void SynRad::OffsetFormula(char *expression,int offset,int filter) {
	//will increase or decrease facet numbers in a formula
	//only applies to facet numbers larger than "filter" parameter

	vector<string> prefixes;
	prefixes.push_back("A");
	prefixes.push_back("D");
	prefixes.push_back("H");
	prefixes.push_back("AR");
	prefixes.push_back(","); //for sum formulas

	string expr=expression; //convert char array to string

	size_t pos=0; //analyzed until this position
	while (pos<expr.size()) { //while not end of expression
		
		vector<size_t> location; //for each prefix, we store where it was found

		for (int j=0;j<(int)prefixes.size();j++) { //try all expressions
			location.push_back(expr.find(prefixes[j],pos));
		}
		size_t minPos=string::npos;
		size_t maxLength=0;
		for (int j=0;j<(int)prefixes.size();j++)  //try all expressions, find first prefix location
			if (location[j]<minPos) minPos=location[j];
		for (int j=0;j<(int)prefixes.size();j++)  //try all expressions, find longest prefix at location
			if (location[j]==minPos && prefixes[j].size()>maxLength) maxLength=prefixes[j].size();
		int digitsLength=0;
		if (minPos!=string::npos) { //found expression, let's find tailing facet number digits
			while ((minPos+maxLength+digitsLength)<expr.length() && expr[minPos+maxLength+digitsLength]>='0' && expr[minPos+maxLength+digitsLength]<='9')
				digitsLength++;
			if (digitsLength>0) { //there was a digit after the prefix
				int facetNumber;
				if (sscanf(expr.substr(minPos+maxLength,digitsLength).c_str(),"%d",&facetNumber)){
					if (facetNumber>filter) {
						char tmp[10];
						sprintf(tmp,"%d",facetNumber+=offset);
						expr.replace(minPos+maxLength,digitsLength,tmp);
					}
				}
			}
		}
		if (minPos!=string::npos) pos=minPos+maxLength+digitsLength;
		else pos=minPos;
	}
	strcpy(expression,expr.c_str());
}

void SynRad::RenumberFormulas(int startId) {
	for (int i=0;i<nbFormula;i++) {
		char expression[1024];
		strcpy(expression,this->formulas[i].parser->GetExpression());
		OffsetFormula(expression,-1,startId);
		this->formulas[i].parser->SetExpression(expression);
		this->formulas[i].parser->Parse();
	}
}

static float lastWrite=0.0f;
static llong lastNbD=0;

void InitLog(char *tmp) {
	FILE *f = fopen("D:\\c++\\synrad\\tests\\data\\scan.txt","w");
	fprintf(f,tmp);
	fprintf(f,"\n");
	fclose(f);
}

void Log(char *tmp) {
	FILE *f = fopen("D:\\c++\\synrad\\tests\\data\\scan.txt","a");
	fprintf(f,tmp);
	fprintf(f,"\n");
	fclose(f);
}

// ----------------------------------------------------------------------------
void SynRad::ResetAutoSaveTimer() {
	if (autoSaveSimuOnly) lastSaveTimeSimu=worker.simuTime+(m_fTime-worker.startTime);
	else lastSaveTime=m_fTime;
}

//-----------------------------------------------------------------------------
void SynRad::AutoSave(BOOL crashSave) {
	if (!changedSinceSave) return;
	GLProgress *progressDlg2 = new GLProgress("Peforming autosave...","Please wait");
	progressDlg2->SetProgress(0.0);
	progressDlg2->SetVisible(TRUE);
	//GLWindowManager::Repaint();
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	char filename[2048];
	sprintf(filename,"%s\\Synrad_AutoSave.syn7z",CWD);
	try {
		ResetAutoSaveTimer();
		worker.SaveGeometry(filename,progressDlg2,FALSE,FALSE,TRUE,crashSave);
	} catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),worker.GetFileName());
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
	}
	//lastSaveTime=(worker.simuTime+(m_fTime-worker.startTime));
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
}

//-------------------------------------------------------------------------------
void SynRad::CheckForRecovery() {
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	char filename[2048];
	sprintf(filename,"%s\\Synrad_AutoSave.syn7z",CWD);
	if (FileUtils::Exist(filename)) {
		int rep = GLMessageBox::Display("Autosave file found. Load it now?\nIf you click CANCEL the file will be discarded.","Autosave recovery",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);
		if( rep == GLDLG_OK ) {
			LoadFile(filename);
			RemoveRecent(filename);
		}
		return;
	}
	sprintf(filename,"%s\\Synrad_AutoSave.syn",CWD);
	if (FileUtils::Exist(filename)) {
		int rep = GLMessageBox::Display("Autosave file found. Load it now?\nIf you click CANCEL the file will be discarded.","Autosave recovery",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);
		if( rep == GLDLG_OK ) {
			LoadFile(filename);
		}
	}
}


//-----------------------------------------------------------------------------
// Name: FrameMove()
// Desc: Called once per frame, the call is the entry point for animating
//       the scene.
//-----------------------------------------------------------------------------
int SynRad::FrameMove()
{
	char tmp[256];
	Geometry *geom = worker.GetGeometry();

	//Autosave routines	
	if (geom->IsLoaded()) {
		if (autoSaveSimuOnly) {
			if (worker.running) {
				if(((worker.simuTime+(m_fTime-worker.startTime))-lastSaveTimeSimu)>=(float)autoSaveFrequency*60.0f) {
					AutoSave(); 
				}
			}
		} else {
			if((m_fTime-lastSaveTime)>=(float)autoSaveFrequency*60.0f) {
				AutoSave();
			}
		}
	}

	// Simulation monitoring
	if(globalSettings) globalSettings->SMPUpdate(m_fTime);
	if(profilePlotter) profilePlotter->Update(m_fTime);
	if(spectrumPlotter) spectrumPlotter->Update(m_fTime);
	if(texturePlotter) texturePlotter->Update(m_fTime);
	if(worker.running) {
		if( m_fTime-lastUpdate>=1.0f ) {

			// Update hits
			worker.Update(m_fTime);
			if( textureSettings ) textureSettings->Update();
			lastUpdate = m_fTime;

			// Update timing measurements
			if( worker.nbHit!=lastNbHit || worker.nbDesorption!=lastNbDes ) {
				double dTime = (double)(m_fTime-lastMeasTime);
				hps = (double)(worker.nbHit-lastNbHit)/dTime;
				dps = (double)(worker.nbDesorption-lastNbDes)/dTime;
				if( lastHps!=0.0 ) {
					hps = 0.2*(hps) + 0.8*lastHps;
					dps = 0.2*(dps) + 0.8*lastDps;
				}
				lastHps = hps;
				lastDps = dps;
				lastNbHit = worker.nbHit;
				lastNbDes = worker.nbDesorption;
				lastMeasTime = m_fTime;
			}

		}
		sprintf(tmp,"Running: %s",FormatTime(worker.simuTime+(m_fTime-worker.startTime)));
		sTime->SetText( tmp );
	} else {
		if( worker.simuTime>0.0 ) {
			hps = (double)(worker.nbHit-nbHitStart)/worker.simuTime;
			dps = (double)(worker.nbDesorption-nbDesStart)/worker.simuTime;
		} else {
			hps = 0.0;
			dps = 0.0;
		}
		sprintf(tmp,"Stopped: %s",FormatTime(worker.simuTime));
		sTime->SetText( tmp );
	}

	if( (m_fTime-worker.startTime <= 2.0f) && worker.running ) {
		hitNumber->SetText("Starting...");
		desNumber->SetText("Starting...");
		doseNumber->SetText("Starting...");
	} else {
		sprintf(tmp,"%s (%s)",FormatInt(worker.nbHit,"hit"),FormatPS(hps,"hit"));
		hitNumber->SetText(tmp);
		sprintf(tmp,"%s (%s)",FormatInt(worker.nbDesorption,"gen"),FormatPS(dps,"gen"));
		desNumber->SetText(tmp);
	}


	if( worker.nbLeakTotal ) {
		sprintf(tmp,"%g (%.4f%%)",(double)worker.nbLeakTotal,(double)(worker.nbLeakTotal*100)/(double)worker.nbDesorption);
		leakNumber->SetText(tmp);
	} else {
		leakNumber->SetText("None");
	}

	if( worker.nbDesorption ) {
		double no_scans=(double)worker.nbDesorption/
			((worker.nbTrajPoints>0)?(double)worker.nbTrajPoints:1.0);
		sprintf(tmp,"Scn:%.1f F=%.3g P=%.3g",(no_scans==1.0)?0.0:no_scans,worker.totalFlux/no_scans,worker.totalPower/no_scans);
		if (worker.nbTrajPoints>0) doseNumber->SetText(tmp);
	} else {
		doseNumber->SetText("");
	}

	resetSimu->SetEnabled(!worker.running&&worker.nbDesorption>0);

	if (worker.running) {
		startSimu->SetText("Pause");
	} else if (worker.nbHit>0) {
		startSimu->SetText("Resume");
	} else {
		startSimu->SetText("Begin");
	}

	// Facet parameters and hits
	if( viewer[0]->SelectionChanged() ||
		viewer[1]->SelectionChanged() ||
		viewer[2]->SelectionChanged() ||
		viewer[3]->SelectionChanged() ) {
			UpdateFacetParams(TRUE);
	}

	UpdateFacetHits();

	// Formulas
	UpdateFormula();

	// Sleep a bit to avoid unwanted CPU load
	if( viewer[0]->IsDragging() || 
		viewer[1]->IsDragging() || 
		viewer[2]->IsDragging() || 
		viewer[3]->IsDragging() || !worker.running )
		SDL_Delay(32);
	else
		SDL_Delay(60);

	return GL_OK;
}

// ----------------------------------------------------------------

void SynRad::UpdateFacetHits(BOOL all) {
	char tmp[256];
	Geometry *geom = worker.GetGeometry();
	Worker *worker=&(theApp->worker);
	double no_scans;
	if (worker->nbTrajPoints==0 || worker->nbDesorption==0) no_scans=1.0;
	else no_scans=(double)worker->nbDesorption/(double)worker->nbTrajPoints;   
	try{
		// Facet list
		if(geom->IsLoaded()) {

			int sR,eR;
			facetList->GetVisibleRows(&sR,&eR);
			if (all) {sR=0;eR=geom->GetNbFacet()-1;}
			for(int i=sR;i<=eR;i++) {
				int facetId=facetList->GetValueInt(i,0)-1;
				if (facetId==-2) facetId=i;
				if (i>=geom->GetNbFacet()) {
					char errMsg[512];
					sprintf(errMsg,"Synrad::UpdateFacetHits()\nError while updating facet hits. Was looking for facet #%d in list.\nSynrad will now autosave and crash.",i+1);
					GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
					AutoSave(TRUE);
				}
				Facet *f = geom->GetFacet(facetId);
				sprintf(tmp,"%d",facetId+1);
				facetList->SetValueAt(0,i,tmp);
				facetList->SetColumnLabel(1,"Hits");
				sprintf(tmp,"%I64d",f->sh.counter.nbHit);
				facetList->SetValueAt(1,i,tmp);
				sprintf(tmp,"%.3g",f->sh.counter.fluxAbs/no_scans);
				facetList->SetValueAt(2,i,tmp);
				sprintf(tmp,"%.3g",f->sh.counter.powerAbs/no_scans);
				facetList->SetValueAt(3,i,tmp);
				sprintf(tmp,"%I64d",f->sh.counter.nbAbsorbed);
				facetList->SetValueAt(4,i,tmp);

			}

		}
	} 	 catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"%s\nError while updating facet hits",e.GetMsg());
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
	}

}

//-----------------------------------------------------------------------------
// Name: SetParam()
// Desc: print the specified param
//-----------------------------------------------------------------------------
void SynRad::SetParam(GLTextField *txt,double value)
{
	char tmp[256];
	sprintf(tmp,"%g",value);
	txt->SetText(tmp);
}

//-----------------------------------------------------------------------------
// Name: FormatInt()
// Desc: Format an integer in K,M,G,..
//-----------------------------------------------------------------------------
char *SynRad::FormatInt(llong v,char *unit)
{

	double x = (double)v;

	static char ret[64];
	if(x<1E3) {
		sprintf(ret,"%g %s",(double)x,unit);
	} else if(x<1E6) {
		sprintf(ret,"%.1f K%s",x/1E3,unit);
	} else if(x<1E9) {
		sprintf(ret,"%.2f M%s",x/1E6,unit);
	} else if(x<1E12) {
		sprintf(ret,"%.2f G%s",x/1E9,unit);
	} else {
		sprintf(ret,"%.2f T%s",x/1E12,unit);
	}

	return ret;

}

//-----------------------------------------------------------------------------
// Name: FormatTime()
// Desc: Format time in HH:MM:SS
//-----------------------------------------------------------------------------
char *SynRad::FormatTime(float t) {
	static char ret[64];
	int nbSec = (int)(t+0.5f);
	sprintf(ret,"%02d:%02d:%02d",nbSec/3600,(nbSec%3600)/60,nbSec%60);
	return ret;
}

//-----------------------------------------------------------------------------
// Name: FormatPS()
// Desc: Format a double in K,M,G,.. per sec
//-----------------------------------------------------------------------------
char *SynRad::FormatPS(double v,char *unit)
{

	static char ret[64];
	if(v<1000.0) {
		sprintf(ret,"%.1f %s/s",v,unit);
	} else if(v<1000000.0) {
		sprintf(ret,"%.1f K%s/s",v/1000.0,unit);
	} else if(v<1000000000.0) {
		sprintf(ret,"%.1f M%s/s",v/1000000.0,unit);
	} else {
		sprintf(ret,"%.1f G%s/s",v/1000000000.0,unit);
	}

	return ret;

}

//-----------------------------------------------------------------------------
// Name: FormatSize()
// Desc: Format a double in K,M,G,.. per sec
//-----------------------------------------------------------------------------
char *SynRad::FormatSize(DWORD size)
{

	static char ret[64];
	if(size<1024UL) {
		sprintf(ret,"%d Bytes",size);
	} else if(size<1048576UL) {
		sprintf(ret,"%.1f KB",(double)size/1024.0);
	} else if(size<1073741824UL) {
		sprintf(ret,"%.1f MB",(double)size/1048576.0);
	} else {
		sprintf(ret,"%.1f GB",(double)size/1073741824.0);
	}

	return ret;

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Name: RestoreDeviceObjects()
// Desc: Initialize scene objects.
//-----------------------------------------------------------------------------
int SynRad::RestoreDeviceObjects()
{
	Geometry *geom = worker.GetGeometry();
	geom->RestoreDeviceObjects();
	//worker.Update(0.0f);

	// Restore dialog which are not displayed
	// Those which are displayed are invalidated by the window manager
	RVALIDATE_DLG(formulaSettings);
	RVALIDATE_DLG(collapseSettings);
	RVALIDATE_DLG(regionInfo);
	RVALIDATE_DLG(moveVertex);
	RVALIDATE_DLG(scaleVertex);
	RVALIDATE_DLG(scaleFacet);
	RVALIDATE_DLG(selectDialog);
	RVALIDATE_DLG(moveFacet);
	RVALIDATE_DLG(exportDesorption);
	RVALIDATE_DLG(mirrorFacet);
	RVALIDATE_DLG(rotateFacet);
	RVALIDATE_DLG(alignFacet);
	RVALIDATE_DLG(addVertex);
	RVALIDATE_DLG(facetMesh);
	RVALIDATE_DLG(facetDetails);
	RVALIDATE_DLG(viewer3DSettings);
	RVALIDATE_DLG(textureSettings);
	RVALIDATE_DLG(globalSettings);
	RVALIDATE_DLG(facetCoordinates);
	RVALIDATE_DLG(vertexCoordinates);
	RVALIDATE_DLG(profilePlotter);
	RVALIDATE_DLG(spectrumPlotter);
	RVALIDATE_DLG(viewEditor);
	RVALIDATE_DLG(texturePlotter);
	//RVALIDATE_DLG(outgassingMap);

	UpdateTitle();
	return GL_OK;
}

//-----------------------------------------------------------------------------
// Name: InvalidateDeviceObjects()
// Desc: Free all alocated resource
//-----------------------------------------------------------------------------

int SynRad::InvalidateDeviceObjects()
{
	Geometry *geom = worker.GetGeometry();
	geom->InvalidateDeviceObjects();

	// Invalidate dialog which are not displayed
	// Those which are displayed are invalidated by the window manager
	IVALIDATE_DLG(formulaSettings);
	IVALIDATE_DLG(collapseSettings);
	IVALIDATE_DLG(regionInfo);
	IVALIDATE_DLG(moveVertex);
	IVALIDATE_DLG(scaleVertex);
	IVALIDATE_DLG(scaleFacet);
	IVALIDATE_DLG(selectDialog);
	IVALIDATE_DLG(moveFacet);
	IVALIDATE_DLG(exportDesorption);
	IVALIDATE_DLG(mirrorFacet);
	IVALIDATE_DLG(rotateFacet);
	IVALIDATE_DLG(alignFacet);
	IVALIDATE_DLG(addVertex);
	IVALIDATE_DLG(facetDetails);
	IVALIDATE_DLG(viewer3DSettings);
	IVALIDATE_DLG(textureSettings);
	IVALIDATE_DLG(globalSettings);
	IVALIDATE_DLG(facetCoordinates);
	IVALIDATE_DLG(vertexCoordinates);
	IVALIDATE_DLG(profilePlotter);
	IVALIDATE_DLG(spectrumPlotter);
	IVALIDATE_DLG(viewEditor);
	IVALIDATE_DLG(texturePlotter);
	//IVALIDATE_DLG(outgassingMap);

	return GL_OK;
}

//-----------------------------------------------------------------------------
// Name: InvalidateDeviceObjects()
// Desc: Called before exiting
//-----------------------------------------------------------------------------
int SynRad::OnExit() {
	SaveConfig();
	worker.Exit();
	remove("Synrad_AutoSave.syn");
	remove("Synrad_AutoSave.syn7z");
	//empty TMP directory
	char tmp[2048];
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	sprintf(tmp,"del /Q \"%s\\tmp\\*.*\"",CWD);
	system(tmp);
	//win_sparkle_cleanup();
	return GL_OK;
}

//-----------------------------------------------------------------------------

void SynRad::UpdateCurrentDir(char *fileName) {

	strcpy(currentDir,fileName);
	char *dp = strrchr(currentDir,'\\');
	if(!dp) dp = strrchr(currentDir,'/');
	if(dp) *dp=0;

}

//-----------------------------------------------------------------------------

void SynRad::UpdateCurrentSelDir(char *fileName) {

	strcpy(currentSelDir,fileName);
	char *dp = strrchr(currentSelDir,'\\');
	if(!dp) dp = strrchr(currentSelDir,'/');
	if(dp) *dp=0;

}

//-----------------------------------------------------------------------------

void SynRad::UpdateTitle() {

	static char title[128];

	Geometry *geom = worker.GetGeometry();

	if( !geom->IsLoaded() ) {
		sprintf(title,"%s",APP_NAME);
	} else {
		if( geom->viewStruct<0 ) {
			sprintf(title,"%s [%s]",APP_NAME,worker.GetShortFileName());
		} else {
			sprintf(title,"%s [%s: Struct #%d %s]",APP_NAME,worker.GetShortFileName(),geom->viewStruct+1,geom->GetStructureName(geom->viewStruct));
		}
	}

	SetTitle(title);

}

//-----------------------------------------------------------------------------

void SynRad::SaveFileAs() {

	FILENAME *fn = GLFileBox::SaveFile(currentDir,worker.GetShortFileName(),"Save File",fileSFilters,nbSFilter);

	GLProgress *progressDlg2 = new GLProgress("Saving file...","Please wait");
	progressDlg2->SetProgress(0.0);
	progressDlg2->SetVisible(TRUE);
	//GLWindowManager::Repaint();  
	if( fn ) {
		try {
			worker.SaveGeometry(fn->fullName,progressDlg2);
			ResetAutoSaveTimer();
			changedSinceSave=FALSE;
			UpdateCurrentDir(worker.fullFileName);
			UpdateTitle();
			AddRecent(worker.fullFileName);
		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
			RemoveRecent(fn->fullName);
		}
	}

	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
}

//-----------------------------------------------------------------------------

void SynRad::ExportSelection() {

	Geometry *geom = worker.GetGeometry();
	if(geom->GetNbSelected()==0) {
		GLMessageBox::Display("Empty selection","Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	FILENAME *fn = GLFileBox::SaveFile(currentDir,worker.GetShortFileName(),"Export selection",fileSFilters,nbSFilter);
	GLProgress *progressDlg2 = new GLProgress("Saving file...","Please wait");
	progressDlg2->SetProgress(0.0);
	progressDlg2->SetVisible(TRUE);
	//GLWindowManager::Repaint();
	if( fn ) {

		try {
			worker.SaveGeometry(fn->fullName,progressDlg2,TRUE,TRUE);
			UpdateCurrentDir(fn->fullName);
			UpdateTitle();
		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		}

	}

	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
}

//-----------------------------------------------------------------------------

void SynRad::ExportTexture(int mode) {

	Geometry *geom = worker.GetGeometry();
	if(geom->GetNbSelected()==0) {
		GLMessageBox::Display("Empty selection","Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	FILENAME *fn = GLFileBox::SaveFile(currentDir,NULL,"Save Texture File",fileTexFilters,nbTexFilter);

	if( fn ) {

		try {
			worker.ExportTextures(fn->fullName,mode,TRUE,TRUE);
			//UpdateCurrentDir(fn->fullName);
			//UpdateTitle();
		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		}

	}

}

void SynRad::SaveFile() {
	if(strlen(worker.fullFileName)>0){

		GLProgress *progressDlg2 = new GLProgress("Saving...","Please wait");
		progressDlg2->SetProgress(0.5);
		progressDlg2->SetVisible(TRUE);
		//GLWindowManager::Repaint();

		try {
			worker.SaveGeometry(worker.fullFileName,progressDlg2,FALSE);
			ResetAutoSaveTimer();
		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),worker.GetFileName());
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		}
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		changedSinceSave=FALSE;

	} else SaveFileAs();
}

//-----------------------------------------------------------------------------

void SynRad::SaveSelection() {

	FileWriter *f = NULL;
	Geometry *geom = worker.GetGeometry();
	if(geom->GetNbSelected()==0) return;
	GLProgress *progressDlg2 = new GLProgress("Saving file","Please wait");
	progressDlg2->SetProgress(0.5);
	progressDlg2->SetVisible(TRUE);
	//GLWindowManager::Repaint();

	FILENAME *fn = GLFileBox::SaveFile(currentSelDir,worker.GetShortFileName(),"Save selection",fileSelFilters,nbSelFilter);

	if( fn ) {

		try {


			char *ext = fn->fullName+strlen(fn->fullName)-4;

			if(!(*ext=='.') ) {
				sprintf(fn->fullName,"%s.sel",fn->fullName); //set to default SEL format
				ext = strrchr(fn->fullName,'.');
			}
			ext++;

			f = new FileWriter(fn->fullName);
			int nbSelected = geom->GetNbSelected();
			int nbFacet = geom->GetNbFacet();
			for(int i=0;i<nbFacet;i++) {
				if(geom->GetFacet(i)->selected) f->WriteInt(i,"\n");
			}

		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		}

		SAFE_DELETE(f);

	}
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
	changedSinceSave=FALSE;
}

//-----------------------------------------------------------------------------

void SynRad::LoadSelection(char *fName) {

	char fullName[512];
	strcpy(fullName,"");
	FileReader *f=NULL;

	if(fName==NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentSelDir,NULL,"Load Selection",fileSelFilters,nbSelFilter);
		if( fn )
			strcpy(fullName,fn->fullName);
	} else {
		strcpy(fullName,fName);
	}

	if(strlen(fullName)==0) return;

	try {

		Geometry *geom = worker.GetGeometry();
		geom->Unselect();
		int nbFacet = geom->GetNbFacet();

		f = new FileReader(fullName);
		while(!f->IsEof()) {
			int s = f->ReadInt();
			if( s>=0 && s<nbFacet ) geom->Select(s);
		}
		geom->UpdateSelection();

		UpdateFacetParams(TRUE);
		UpdateCurrentSelDir(fullName);

	} catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fullName);
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);

	}

	SAFE_DELETE(f);
	changedSinceSave=FALSE;

}

//-----------------------------------------------------------------------------

void SynRad::LoadFile(char *fName) {


	char fullName[512];
	char shortName[512];
	strcpy(fullName,"");

	if(fName==NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir,NULL,"Open File",fileLFilters,nbLFilter);
		if( fn )
			strcpy(fullName,fn->fullName);
	} else {
		strcpy(fullName,fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Preparing to load file...","Please wait");
	progressDlg2->SetVisible(TRUE);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if(strlen(fullName)==0) {
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		return;
	}


	char *lPart = strrchr(fullName,'\\');
	if(lPart) strcpy(shortName,lPart+1);
	else strcpy(shortName,fullName);

	try {
		worker.ClearRegions();
		worker.LoadGeometry(fullName);

		Geometry *geom = worker.GetGeometry();

		// Default initialisation
		viewer[0]->SetWorker(&worker);
		viewer[1]->SetWorker(&worker);
		viewer[2]->SetWorker(&worker);
		viewer[3]->SetWorker(&worker);
		startSimu->SetEnabled(TRUE);
		ClearFacetParams();
		nbDesStart = worker.nbDesorption;
		nbHitStart = worker.nbHit;
		AddRecent(fullName);
		geom->viewStruct = -1;

		RebuildPARMenus();
		UpdateStructMenu();
		if(profilePlotter) profilePlotter->Reset();
		if(spectrumPlotter) spectrumPlotter->Reset();
		UpdateCurrentDir(fullName);

		// Check non simple polygon
		progressDlg2->SetMessage("Checking for non simple polygons...");

		geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();
		// Set up view
		// Default
		viewer[0]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[0]->ToFrontView();
		viewer[1]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[1]->ToTopView();
		viewer[2]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[2]->ToSideView();
		viewer[3]->SetProjection(PERSPECTIVE_PROJ);
		viewer[3]->ToFrontView();
		SelectViewer(0);

		ResetAutoSaveTimer();
		if(profilePlotter) profilePlotter->Refresh();
		if(spectrumPlotter) spectrumPlotter->Refresh();
		if(texturePlotter) texturePlotter->Update(m_fTime,TRUE);
		//if(outgassingMap) outgassingMap->Update(m_fTime,TRUE);
		if(facetDetails) facetDetails->Update();
		if(facetCoordinates) facetCoordinates->UpdateFromSelection();
		if(vertexCoordinates) vertexCoordinates->Update();

	} catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),shortName);
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		RemoveRecent(fName);

	}
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
	changedSinceSave=FALSE;
}

//-----------------------------------------------------------------------------

void SynRad::InsertGeometry(BOOL newStr,char *fName) {
	if (!AskToReset()) return;
	ResetSimulation(FALSE);

	char fullName[512];
	char shortName[512];
	strcpy(fullName,"");

	if(fName==NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir,NULL,"Open File",fileInsFilters,nbInsFilter);
		if( fn )
			strcpy(fullName,fn->fullName);
	} else {
		strcpy(fullName,fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Loading file...","Please wait");
	progressDlg2->SetVisible(TRUE);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if(strlen(fullName)==0) {
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		return;
	}


	char *lPart = strrchr(fullName,'\\');
	if(lPart) strcpy(shortName,lPart+1);
	else strcpy(shortName,fullName);

	try {

		worker.InsertGeometry(newStr,fullName);
		Geometry *geom = worker.GetGeometry();

		startSimu->SetEnabled(TRUE);

		AddRecent(fullName);
		geom->viewStruct = -1;

		RebuildPARMenus();
		UpdateStructMenu();
		if(profilePlotter) profilePlotter->Reset();
		if(spectrumPlotter) spectrumPlotter->Reset();

		geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();

		if(profilePlotter) profilePlotter->Refresh();
		if(spectrumPlotter) spectrumPlotter->Refresh();
		if(texturePlotter) texturePlotter->Update(m_fTime,TRUE);
		//if(outgassingMap) outgassingMap->Update(m_fTime,TRUE);
		if(facetDetails) facetDetails->Update();
		if(facetCoordinates) facetCoordinates->UpdateFromSelection();
		if(vertexCoordinates) vertexCoordinates->Update();
		if (regionInfo) regionInfo->Update();

	} catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),shortName);
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		RemoveRecent(fName);
	}
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
	changedSinceSave=TRUE;
}

//-----------------------------------------------------------------------------
// Name: UpdateModelParams()
// Desc: Update displayed model parameter on geometry ghange
//-----------------------------------------------------------------------------
void SynRad::UpdateModelParams() {

	Geometry *geom = worker.GetGeometry();
	char tmp[256];
	double sumArea=0;
	facetList->SetSize(5,geom->GetNbFacet(),TRUE);
	facetList->SetColumnWidths((int*)cWidth);
	facetList->SetColumnLabels((char **)cName);
	UpdateFacetHits(TRUE);
	AABB bb = geom->GetBB();

	for (int i=0;i<geom->GetNbFacet();i++) {
		Facet *f=geom->GetFacet(i);
		if (f->sh.area>0) sumArea+=f->sh.area;
	}

	sprintf(tmp,"V:%d F:%d Dim:(%g,%g,%g) Area:%g",geom->GetNbVertex(),geom->GetNbFacet(),
		(bb.max.x-bb.min.x),(bb.max.y-bb.min.y),(bb.max.z-bb.min.z),sumArea);
	geomNumber->SetText(tmp);

}

//-----------------------------------------------------------------------------
// Name: AddFormula()
// Desc: Add a formula
//-----------------------------------------------------------------------------
void SynRad::AddFormula(GLParser *f,BOOL doUpdate) {

	if( f ) {
		if( nbFormula<MAX_FORMULA ) {
			formulas[nbFormula].parser = f;
			formulas[nbFormula].name = new GLLabel(f->GetName());
			Add(formulas[nbFormula].name);
			formulas[nbFormula].value = new GLTextField(0,"");
			formulas[nbFormula].value->SetEditable(FALSE);
			Add(formulas[nbFormula].value);
			formulas[nbFormula].setBtn = new GLButton(0,"...");
			Add(formulas[nbFormula].setBtn);
			nbFormula++;
			PlaceComponents();
			if(doUpdate) UpdateFormula();
		} else {
			SAFE_DELETE(f);
		}
	}

}

void SynRad::ClearFormula() {

	for(int i=0;i<nbFormula;i++) {
		wnd->PostDelete(formulas[i].name);
		wnd->PostDelete(formulas[i].value);
		wnd->PostDelete(formulas[i].setBtn);
		formulas[i].name=NULL;
		formulas[i].value=NULL;
		formulas[i].setBtn=NULL;
		SAFE_DELETE(formulas[i].parser);
	}
	nbFormula=0;
	PlaceComponents();

}

void SynRad::AddFormula(char *fName,char *formula) {

	GLParser *f = new GLParser();
	f->SetExpression(formula);
	f->SetName(fName);
	f->Parse();
	AddFormula(f,FALSE);

}

//-----------------------------------------------------------------------------
// Name: ProcessFormulaButtons()
// Desc: Handle forumla button event
//-----------------------------------------------------------------------------
void SynRad::ProcessFormulaButtons(GLComponent *src) {

	// Search formula buttons
	BOOL found = FALSE;
	int i=0;
	while(!found && i<nbFormula) {
		found = (src == formulas[i].setBtn);
		if(!found) i++;
	}
	if( found ) {
		if( !formulaSettings ) formulaSettings = new FormulaSettings();
		if( formulaSettings->EditFormula(formulas[i].parser) ) {
			// Apply change
			formulas[i].name->SetText(formulas[i].parser->GetName());
			UpdateFormula();
		} else {
			// Delete
			wnd->PostDelete(formulas[i].name);
			wnd->PostDelete(formulas[i].value);
			wnd->PostDelete(formulas[i].setBtn);
			formulas[i].name=NULL;
			formulas[i].value=NULL;
			formulas[i].setBtn=NULL;
			SAFE_DELETE(formulas[i].parser);
			for(int j=i;j<nbFormula-1;j++)
				formulas[j] = formulas[j+1];
			nbFormula--;
			PlaceComponents();
			UpdateFormula();
		}
	}

}

//-----------------------------------------------------------------------------

void SynRad::StartStopSimulation() {
	if ((int)worker.regions.size()==0) {
		GLMessageBox::Display("No regions loaded","Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	worker.StartStop(m_fTime,0);
	if(profilePlotter) profilePlotter->Update(m_fTime,TRUE);
	if(spectrumPlotter) spectrumPlotter->Update(m_fTime,TRUE);
	if(texturePlotter) texturePlotter->Update(m_fTime,TRUE);

	// Frame rate measurement
	lastMeasTime = m_fTime;
	dps = 0.0;
	hps = 0.0;
	lastHps = hps;
	lastDps = dps;
	lastNbHit = worker.nbHit;
	lastNbDes = worker.nbDesorption;
	lastUpdate = 0.0;

}

//-----------------------------------------------------------------------------

void SynRad::ResetSimulation(BOOL askConfirm) {

	BOOL ok = TRUE;
	if( askConfirm ) 
		ok = GLMessageBox::Display("Reset simulation ?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK;

	if(ok) {
		worker.Reset(m_fTime);
		nbDesStart = 0;
		nbHitStart = 0;
	}
	//resetSimu->SetEnabled(FALSE);
	if(profilePlotter) profilePlotter->Update(m_fTime,TRUE);
	if(spectrumPlotter) spectrumPlotter->Update(m_fTime,TRUE);
	if(texturePlotter) texturePlotter->Update(m_fTime,TRUE);

}

//-----------------------------------------------------------------------------

void SynRad::SelectViewer(int s) {

	curViewer = s;
	for(int i=0;i<MAX_VIEWER;i++) viewer[i]->SetSelected(i==curViewer);
	UpdateViewerParams();

}

//-----------------------------------------------------------------------------
// Name: EventProc()
// Desc: Message proc function to handle key and mouse input
//-----------------------------------------------------------------------------
void SynRad::ProcessMessage(GLComponent *src,int message)
{
	Geometry *geom = worker.GetGeometry();

	switch(message) {

		//MENU --------------------------------------------------------------------
	case MSG_MENU:
		switch(src->GetId()) {
		case MENU_FILE_LOAD:
			if (AskToSave()) {
				if (worker.running) worker.Stop_Public();
				LoadFile();
			}
			break;
		case MENU_REGIONS_NEW:
			NewRegion();					
			break;
		case MENU_REGIONS_LOADPAR:
			//if (AskToSave()) {
			if (worker.running) worker.Stop_Public();
			LoadParam();
			//}
			break;

		case MENU_REGIONS_CLEARALL:
			if( GLMessageBox::Display("Remove all magnetic regions?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				if (worker.running) worker.Stop_Public();
				ClearTraj();
			}
			break;
		case MENU_FILE_INSERTGEO:
			if( geom->IsLoaded() ) {
				if (worker.running) worker.Stop_Public();
				InsertGeometry(FALSE);
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FILE_INSERTGEO_NEWSTR:
			if( geom->IsLoaded() ) {
				if (worker.running) worker.Stop_Public();
				InsertGeometry(TRUE);
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FILE_SAVEAS:
			if( geom->IsLoaded() ) {
				SaveFileAs();
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FILE_EXPORTMESH:
			ExportSelection();
			break;
		case MENU_FILE_EXPORTTEXTURE_AREA:
			ExportTexture(0);
			break;
		case MENU_FILE_EXPORTTEXTURE_MCHITS:
			ExportTexture(1);
			break;
		case MENU_FILE_EXPORTTEXTURE_FLUX:
			ExportTexture(2);
			break;
		case MENU_FILE_EXPORTTEXTURE_POWER:
			ExportTexture(3);
			break;
		case MENU_FILE_EXPORTTEXTURE_FLUXPERAREA:
			ExportTexture(4);
			break;
		case MENU_FILE_EXPORTTEXTURE_POWERPERAREA:
			ExportTexture(5);
			break;
		case MENU_FILE_EXPORT_DESORP:
			if (!geom->IsLoaded()) {
				GLMessageBox::Display("No geometry loaded.","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
			if( !exportDesorption ) exportDesorption = new ExportDesorption(geom,&worker);
			exportDesorption->SetVisible(TRUE);
			break;
		case MENU_FILE_SAVE:
			if( geom->IsLoaded() ) SaveFile();
			else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FILE_EXIT:
			if (AskToSave()) Exit();
			break;
		case MENU_EDIT_3DSETTINGS:
			if( !viewer3DSettings ) viewer3DSettings = new Viewer3DSettings();
			viewer3DSettings->Display(geom,viewer[curViewer]);
			UpdateViewerParams();
			break;
		case MENU_EDIT_TSCALING:
			if( !textureSettings ) textureSettings = new TextureSettings();
			textureSettings->Display(&worker,viewer);
			break;
		case MENU_EDIT_ADDFORMULA:
			if( !formulaSettings ) formulaSettings = new FormulaSettings();
			AddFormula(formulaSettings->NewFormula());
			break;
		case MENU_EDIT_GLOBALSETTINGS:
			if( !globalSettings ) globalSettings = new GlobalSettings();
			globalSettings->Display(&worker);
			break;
		case MENU_FACET_COLLAPSE:
			if( geom->IsLoaded() ) {
				DisplayCollapseDialog();
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FACET_SWAPNORMAL:
			if (AskToReset()) {
				geom->SwapNormal();
				// Send to sub process
				try { worker.Reload(); } catch(Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
				}
			}
			break;
		case MENU_FACET_SHIFTVERTEX:
			if (AskToReset()) {
				geom->ShiftVertex();
				// Send to sub process
				try { worker.Reload(); } catch(Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
				}
			}
			break;
		case MENU_FACET_COORDINATES:

			if(!facetCoordinates) facetCoordinates = new FacetCoordinates();
			facetCoordinates->Display(&worker);
			break;
		case MENU_FACET_MOVE:
			if(!moveFacet) moveFacet = new MoveFacet(geom,&worker);
			moveFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_SCALE:
			if( geom->IsLoaded() ) {
				if( !scaleFacet ) scaleFacet = new ScaleFacet(geom,&worker);

				scaleFacet->SetVisible(TRUE);

			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FACET_MIRROR:
			if(!mirrorFacet) mirrorFacet = new MirrorFacet(geom,&worker);
			mirrorFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_ROTATE:
			if(!rotateFacet) rotateFacet = new RotateFacet(geom,&worker);
			rotateFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_ALIGN:
			if(!alignFacet) alignFacet = new AlignFacet(geom,&worker);
			alignFacet->MemorizeSelection();
			alignFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_PROFPLOTTER:
			if(!profilePlotter) profilePlotter = new ProfilePlotter();
			profilePlotter->Display(&worker);
			break;
		case MENU_FACET_SPECTRUMPLOTTER:
			if(!spectrumPlotter) spectrumPlotter = new SpectrumPlotter();
			spectrumPlotter->Display(&worker);
			break;
		case MENU_FACET_MESH:
			if( !facetMesh ) facetMesh = new FacetMesh();
			facetMesh->EditFacet(&worker);
			UpdateFacetParams();
			break;          
		case MENU_FACET_TEXPLOTTER:
			if( !texturePlotter ) texturePlotter = new TexturePlotter();
			texturePlotter->Display(&worker);
			break;
		/*case MENU_FACET_OUTGASSINGMAP:
			if( !outgassingMap ) outgassingMap = new OutgassingMap();
			outgassingMap->Display(&worker);
			break;*/
		case MENU_FACET_REMOVESEL:
			if( GLMessageBox::Display("Remove selected facets?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				if (AskToReset()) {
					if (worker.running) worker.Stop_Public();
					geom->RemoveSelected();
					//geom->CheckIsolatedVertex();
					UpdateModelParams();
					if (vertexCoordinates) vertexCoordinates->Update();
					if (facetCoordinates) facetCoordinates->UpdateFromSelection();
					if (profilePlotter) profilePlotter->Refresh();
					if (spectrumPlotter) spectrumPlotter->Refresh();
					//if (pressureEvolution) pressureEvolution->Refresh();
					//if (timewisePlotter) timewisePlotter->Refresh();
					// Send to sub process
					try { worker.Reload(); } catch(Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
					}
				}
			}
			break;
		case MENU_FACET_EXPLODE:
			if( GLMessageBox::Display("Explode selected facet?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				if (AskToReset()) {
					int err = geom->ExplodeSelected();
					if( err==-1 ) {
						GLMessageBox::Display("Empty selection","Error",GLDLG_OK,GLDLG_ICONERROR);
					} else if ( err==-2 ) {
						GLMessageBox::Display("All selected facets must have a mesh with boudary correction enabled","Error",GLDLG_OK,GLDLG_ICONERROR);
					} else if( err==0 ) {

						UpdateModelParams();
						UpdateFacetParams(TRUE);
						// Send to sub process
						try { worker.Reload(); } catch(Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
						}
					}
				}
			}
			break;
		case MENU_FACET_DETAILS:
			if( facetDetails==NULL ) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
			break;

		case MENU_REGIONS_REGIONINFO:
			if ((int)worker.regions.size()>0) {
				if( regionInfo==NULL ) regionInfo = new RegionInfo(&worker);
				regionInfo->SetVisible(TRUE);
			} else {
				GLMessageBox::Display("No regions loaded","Error",GLDLG_OK,GLDLG_ICONERROR);
			}
			break;
		case MENU_FACET_SELECTALL:
			geom->SelectAll();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTSTICK:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.sticking != 0.0 && !geom->GetFacet(i)->IsLinkFacet() )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTTRANS:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.opacity != 1.0 && geom->GetFacet(i)->sh.opacity != 2.0 )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTREFL:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++) {
				Facet *f = geom->GetFacet(i);
				if( f->sh.sticking==0.0 && f->sh.opacity>0.0 )
					geom->Select(i);
			}
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECT2SIDE:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.is2sided )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTVOL:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.isVolatile )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTTEXT:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.isTextured!=NULL )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTPROF:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.isProfile!=NULL )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTSPECTRUM:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.hasSpectrum )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTERR:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.sign==0.0 )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTTP:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.teleportDest!=0 )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTDEST:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.superDest!=0 )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTABS:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.counter.nbAbsorbed > 0 )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTHITS:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.counter.nbHit > 0 )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_INVERTSEL:
			for(int i=0;i<geom->GetNbFacet();i++)
				geom->GetFacet(i)->selected = !geom->GetFacet(i)->selected;
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_SELECTION_SELECTFACETNUMBER:
			if( !selectDialog ) selectDialog = new SelectDialog(&worker);
			selectDialog->SetVisible(TRUE);
			break;
		case MENU_FACET_SAVESEL:
			SaveSelection();
			break;
		case MENU_FACET_LOADSEL:
			LoadSelection();
			break;
		case MENU_SELECTION_ADDNEW:
			AddSelection();
			break;
			break;
		case  MENU_SELECTION_CLEARALL:
			if( GLMessageBox::Display("Clear all selections ?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				ClearAllSelections();
			}
			break;
		case MENU_VERTEX_UNSELECTALL:
			geom->UnselectAllVertex();
			break;
		case MENU_VERTEX_SELECTALL:
			geom->SelectAllVertex();
			break;
		case MENU_VERTEX_SELECT_ISOLATED:
			geom->SelectIsolatedVertices();
			break;
		case MENU_VERTEX_CREATE_POLY_CONVEX:
				if (AskToReset()) {
					try {
						geom->CreatePolyFromVertices_Convex();
					} catch (Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(),"Error creating polygon",GLDLG_OK,GLDLG_ICONERROR);
					}
					//UpdateModelParams();
					try {
						worker.Reload();
					} catch(Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
					}
				}
				break;
		case MENU_VERTEX_CREATE_POLY_ORDER:
				if (AskToReset()) {
					try {
						geom->CreatePolyFromVertices_Order();
					} catch (Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(),"Error creating polygon",GLDLG_OK,GLDLG_ICONERROR);
					}
					//UpdateModelParams();
					try {
						worker.Reload();
					} catch(Error &e) {

						GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
					}
				}
				break;
		case MENU_FACET_CREATE_DIFFERENCE:
			if (geom->IsLoaded()){
				try { 
					if (AskToReset()) {
						geom->CreateDifference(); }} catch (Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(),"Error creating polygon",GLDLG_OK,GLDLG_ICONERROR);
						}
						//UpdateModelParams();
						try { worker.Reload(); } catch(Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
						}
			}else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;

		case MENU_VERTEX_SELECT_COPLANAR:
			char *input;
			char tmp[128];
			if( geom->IsLoaded() ) {
				if (geom->GetNbSelectedVertex()!=3) {
					GLMessageBox::Display("Select exactly 3 vertices","Can't define plane",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				sprintf(tmp,"%g",tolerance);
				//sprintf(title,"Pipe L/R = %g",L/R);
				input = GLInputBox::GetInput(tmp,"Tolerance (cm)","Select coplanar vertices");
				if( !input ) return;	
				if(( sscanf(input,"%lf",&tolerance)<=0 )||(tolerance<=0.0)) {
					GLMessageBox::Display("Invalid number","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				try { viewer[curViewer]->SelectCoplanar(tolerance);  } catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error selecting coplanar vertices",GLDLG_OK,GLDLG_ICONERROR);
				}
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_VERTEX_MOVE:
			if( geom->IsLoaded() ) {
				if( !moveVertex ) moveVertex = new MoveVertex(geom,&worker);

				moveVertex->SetVisible(TRUE);

			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_VERTEX_SCALE:
			if( geom->IsLoaded() ) {
				if( !scaleVertex ) scaleVertex = new ScaleVertex(geom,&worker);

				scaleVertex->SetVisible(TRUE);

			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_VERTEX_COORDINATES:

			if(!vertexCoordinates) vertexCoordinates = new VertexCoordinates();
			vertexCoordinates->Display(&worker);
			break;

		case MENU_VERTEX_ADD:
			if( geom->IsLoaded() ) {
				if( !addVertex ) addVertex = new AddVertex(geom,&worker);
				addVertex->SetVisible(TRUE);
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;

		case MENU_VERTEX_REMOVE:
			if( geom->IsLoaded() ) {
				if( GLMessageBox::Display("Remove Selected vertices?\nNote: It will also affect facets that contain them!","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK )  {        
									if (AskToReset()) {
					if (worker.running) worker.Stop_Public();
					geom->RemoveSelectedVertex();
					geom->Rebuild(); //Will recalculate facet parameters
					UpdateModelParams();
					if (vertexCoordinates) vertexCoordinates->Update();
					if (facetCoordinates) facetCoordinates->UpdateFromSelection();
					if (profilePlotter) profilePlotter->Refresh();	
					if (spectrumPlotter) spectrumPlotter->Refresh();
					//if (pressureEvolution) pressureEvolution->Refresh();
					//if (timewisePlotter) timewisePlotter->Refresh();
					// Send to sub process
					try { worker.Reload(); } catch(Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
					}
				}
				}

			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;

		case MENU_VIEW_FULLSCREEN:
			if( m_bWindowed ) {
				ToggleFullscreen();
				PlaceComponents();
			} else {
				Resize(1024,768,TRUE);
			}
			menu->GetSubMenu("View")->SetCheck(MENU_VIEW_FULLSCREEN,!m_bWindowed);
			break;

		case MENU_VIEW_ADDNEW:
			AddView();
			break;
			break;
		case  MENU_VIEW_CLEARALL:
			if( GLMessageBox::Display("Clear all views ?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				ClearAllViews();
			}
			break;

		case MENU_TEST_PIPE0001:
			if (AskToSave()) BuildPipe(0.0001);
			break;
		case MENU_TEST_PIPE1:
			if (AskToSave()) BuildPipe(1.0);
			break;
		case MENU_TEST_PIPE10:
			if (AskToSave()) BuildPipe(10.0);
			break;
		case MENU_TEST_PIPE100:
			if (AskToSave()) BuildPipe(100.0);
			break;
		case MENU_TEST_PIPE1000:
			if (AskToSave()) BuildPipe(1000.0);
			break;
		case MENU_TEST_PIPE10000:
			if (AskToSave()) BuildPipe(10000.0);
			break;
		case MENU_QUICKPIPE:
			if (AskToSave()) QuickPipe();
			break;
		}
		// Load recent menu
		if( src->GetId()>=MENU_FILE_LOADRECENT && src->GetId()<MENU_FILE_LOADRECENT+nbRecent ) {
			if (AskToSave()) {
				if (worker.running) worker.Stop_Public();
				LoadFile(recents[src->GetId()-MENU_FILE_LOADRECENT]);
			}
		}
		// Load recent PAR menu
		if( src->GetId()>=MENU_REGIONS_LOADRECENT && src->GetId()<MENU_REGIONS_LOADRECENT+nbRecentPAR ) {
			if (worker.running) worker.Stop_Public();
			LoadParam(recentPARs[src->GetId()-MENU_REGIONS_LOADRECENT]);
		}
		// Show structure menu
		if( src->GetId()>=MENU_VIEW_STRUCTURE && src->GetId()<MENU_VIEW_STRUCTURE+geom->GetNbStructure() ) {
			geom->viewStruct = src->GetId()-MENU_VIEW_STRUCTURE-1;
			if(src->GetId()>MENU_VIEW_STRUCTURE) geom->Unselect();
			UpdateStructMenu();
		}
		// Load PAR to... menu
		if( src->GetId()>=MENU_REGIONS_LOADTO && src->GetId()<MENU_REGIONS_LOADTO+nbRecentPAR ) {
			if (worker.running) worker.Stop_Public();
			LoadParam(NULL,src->GetId()-MENU_REGIONS_LOADTO);
		}
		// Remove region menu
		if( src->GetId()>=MENU_REGIONS_REMOVE && src->GetId()<MENU_REGIONS_REMOVE+nbRecentPAR ) {
			if (worker.running) worker.Stop_Public();
			RemoveRegion(src->GetId()-MENU_REGIONS_REMOVE);
		}
		if( src->GetId()==MENU_VIEW_NEWSTRUCT) {
			AddStruct();
			UpdateStructMenu();
		}
		if( src->GetId()==MENU_VIEW_DELSTRUCT) {
			DeleteStruct();
			UpdateStructMenu();
		}
		if( src->GetId()==MENU_VIEW_PREVSTRUCT) {
			geom->viewStruct=Remainder(geom->viewStruct-1,geom->GetNbStructure());
			geom->Unselect();
			UpdateStructMenu();
		}
		if( src->GetId()==MENU_VIEW_NEXTSTRUCT) {
			geom->viewStruct=Remainder(geom->viewStruct+1,geom->GetNbStructure());
			geom->Unselect();
			UpdateStructMenu();
		}


		// Select selection
		if( src->GetId()>=MENU_SELECTION_SELECTIONS && src->GetId()<MENU_SELECTION_SELECTIONS+nbSelection ) {
			SelectSelection( src->GetId() - MENU_SELECTION_SELECTIONS );
		}
		// Clear selection
		if( src->GetId()>=MENU_SELECTION_CLEARSELECTIONS && src->GetId()<MENU_SELECTION_CLEARSELECTIONS+nbSelection ) {
			char tmpname[256];
			sprintf(tmpname,"Clear %s?",selections[src->GetId() - MENU_SELECTION_CLEARSELECTIONS].name);
			if( GLMessageBox::Display(tmpname,"Confirmation",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				ClearSelection( src->GetId() - MENU_SELECTION_CLEARSELECTIONS );
			}
		}
		// Memorize selection
		if( src->GetId()>=MENU_SELECTION_MEMORIZESELECTIONS && src->GetId()<MENU_SELECTION_MEMORIZESELECTIONS+nbSelection ) {
			OverWriteSelection( src->GetId() - MENU_SELECTION_MEMORIZESELECTIONS );
		}

		// Select view
		if( src->GetId()>=MENU_VIEW_VIEWS && src->GetId()<MENU_VIEW_VIEWS+nbView ) {
			SelectView( src->GetId() - MENU_VIEW_VIEWS );
		}
		// Clear view
		if( src->GetId()>=MENU_VIEW_CLEARVIEWS && src->GetId()<MENU_VIEW_CLEARVIEWS+nbView ) {
			char tmpname[256];
			sprintf(tmpname,"Clear %s?",views[src->GetId() - MENU_VIEW_CLEARVIEWS].name);
			if( GLMessageBox::Display(tmpname,"Confirmation",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				ClearView( src->GetId() - MENU_VIEW_CLEARVIEWS );
			}
		}
		// Memorize view
		if( src->GetId()>=MENU_VIEW_MEMORIZEVIEWS && src->GetId()<MENU_VIEW_MEMORIZEVIEWS+nbView ) {
			OverWriteView( src->GetId() - MENU_VIEW_MEMORIZEVIEWS );
		}
		break;

		//TEXT --------------------------------------------------------------------
	case MSG_TEXT_UPD:
		if( src == facetSticking ) {
			facetApplyBtn->SetEnabled(TRUE);
		} else if( src == facetRoughness ) {
			facetApplyBtn->SetEnabled(TRUE);
		}else if ( src == facetTeleport ) {
			facetApplyBtn->SetEnabled(TRUE);
		} else if ( src == facetOpacity ) {
			facetApplyBtn->SetEnabled(TRUE);
		} else if ( src == facetSuperDest || src == facetSILabel) {
			facetApplyBtn->SetEnabled(TRUE);
		}
		break;

	case MSG_TEXT:
		if( src == facetSticking || src == facetRoughness || src == facetTeleport 
			|| src == facetOpacity || src == facetSuperDest || src==facetSILabel) {
				ApplyFacetParams();
		}
		break;

		//COMBO -------------------------------------------------------------------
	case MSG_COMBO:
		if ( src == facetReflType ) {
			facetApplyBtn->SetEnabled(TRUE);
			BOOL isMaterial=facetReflType->GetSelectedIndex()>=2;
			facetSticking->SetEditable(TRUE);
			facetRoughness->SetEditable(TRUE);
			facetSticking->SetVisible(!isMaterial);
			facetStickingLabel->SetVisible(!isMaterial);
			facetRoughness->SetVisible(isMaterial);
			facetRoughnessLabel->SetVisible(isMaterial);
		} else if ( src == facetRecType || src == facetSideType || src == facetSpectrumCombo) {
			facetApplyBtn->SetEnabled(TRUE);
		} else if ( src == modeCombo ) {
			changedSinceSave = TRUE;
			worker.generation_mode=modeCombo->GetSelectedIndex(); //fluxwise or powerwise
			worker.Reload();
			UpdateFacetHits();
		}
		break;

		//TOGGLE ------------------------------------------------------------------
	case MSG_TOGGLE:
		// Update viewer flags
		UpdateViewerFlags();
		break;

		//LIST --------------------------------------------------------------------
	case MSG_LIST:
		if( src == facetList && geom->IsLoaded()) {
			int *sels = (int *)malloc((geom->GetNbFacet())*sizeof(int));
			int nbSel;
			facetList->GetSelectedRows(&sels,&nbSel,TRUE);
			geom->Unselect();
			for(int i=0;i<nbSel;i++)
				geom->Select(sels[i]);
			geom->UpdateSelection();
			UpdateFacetParams();
			SAFE_FREE(sels);
		}
		break;

		//GEOMVIEWER ------------------------------------------------------------------
	case MSG_GEOMVIEWER_MAXIMISE:
		if( src==viewer[0] ) {
			AnimateViewerChange(0);
		} else if( src==viewer[1] ) {
			AnimateViewerChange(1);
		} else if( src==viewer[2] ) {
			AnimateViewerChange(2);
		} else if( src==viewer[3] ) {
			AnimateViewerChange(3);
		}
		Place3DViewer();
		break;

	case MSG_GEOMVIEWER_SELECT: {
		SelectViewer(src->GetId());
								}break;

		//BUTTON ------------------------------------------------------------------
	case MSG_BUTTON:
		if( src == startSimu ) {
			changedSinceSave=TRUE;
			StartStopSimulation();
			resetSimu->SetEnabled(!worker.running);
		} else if ( src == resetSimu ) {
			changedSinceSave=TRUE;
			ResetSimulation();
		} else if ( src == facetApplyBtn ) {
			changedSinceSave=TRUE;
			ApplyFacetParams();
		} else if ( src == facetMoreBtn ) {
			if( facetDetails==NULL ) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
		} else if ( src == facetCoordBtn ) {
			if(!facetCoordinates) facetCoordinates = new FacetCoordinates();
			facetCoordinates->Display(&worker);
		} else if ( src == facetTexBtn ) {
			if( !facetMesh ) facetMesh = new FacetMesh();
			facetMesh->EditFacet(&worker);
			changedSinceSave=TRUE;
			UpdateFacetParams();
		} else if ( src == showMoreBtn ) {
			if( !viewer3DSettings ) viewer3DSettings = new Viewer3DSettings();
			viewer3DSettings->Display(geom,viewer[curViewer]);
			UpdateViewerParams();
		} else {
			ProcessFormulaButtons(src);
		}
		break;

		//Panel open/close ---------------------------------------------------------
	case MSG_PANELR:
		PlaceComponents();
		break;

	}

}

void SynRad::AnimateViewerChange(int next,BOOL init) {

	double xs1,ys1,xs2,ys2;
	double xe1,ye1,xe2,ye2;
	int sx = m_screenWidth-205;  
	int fWidth  = m_screenWidth-215;
	int fHeight = m_screenHeight-27;
	int Width2  = fWidth/2-1;
	int Height2 = fHeight/2-1;

	// Reset to layout and make all visible

	if (!init) {
	for(int i=0;i<MAX_VIEWER;i++)  viewer[i]->SetVisible(TRUE);
	viewer[0]->SetBounds(3       ,3        ,Width2,Height2);
	viewer[1]->SetBounds(6+Width2,3        ,Width2,Height2);
	viewer[2]->SetBounds(3       ,6+Height2,Width2,Height2);
	viewer[3]->SetBounds(6+Width2,6+Height2,Width2,Height2);

	if( modeSolo ) {

		// Go from single to layout
		xs1 = (double)3;
		ys1 = (double)3;
		xs2 = (double)fWidth+xs1;
		ys2 = (double)fHeight+ys1;

		switch( next ) {
		case 0:
			xe1 = (double)(3);
			ye1 = (double)(3);
			break;
		case 1:
			xe1 = (double)(5+Width2);
			ye1 = (double)(3);
			break;
		case 2:
			xe1 = (double)(3);
			ye1 = (double)(5+Height2);
			break;
		case 3:
			xe1 = (double)(5+Width2);
			ye1 = (double)(5+Height2);
			break;
		}

		xe2 = (double)(Width2)+xe1;
		ye2 = (double)(Height2)+ye1;

	} else {

		// Go from layout to single
		xe1 = (double)3;
		ye1 = (double)3;
		xe2 = (double)fWidth+xe1;
		ye2 = (double)fHeight+ye1;

		switch( next ) {
		case 0:
			xs1 = (double)(3);
			ys1 = (double)(3);
			break;
		case 1:
			xs1 = (double)(5+Width2);
			ys1 = (double)(3);
			break;
		case 2:
			xs1 = (double)(3);
			ys1 = (double)(5+Height2);
			break;
		case 3:
			xs1 = (double)(5+Width2);
			ys1 = (double)(5+Height2);
			break;
		}

		xs2 = (double)(Width2)+xs1;
		ys2 = (double)(Height2)+ys1;

	}

	double t0=(double)SDL_GetTicks()/1000.0;
	double t1 = t0;
	double T = 0.15;

	while((t1-t0)<T) {
		double t = (t1-t0) / T;
		int x1 = (int)(xs1 + t*(xe1-xs1) + 0.5);
		int y1 = (int)(ys1 + t*(ye1-ys1) + 0.5);
		int x2 = (int)(xs2 + t*(xe2-xs2) + 0.5);
		int y2 = (int)(ys2 + t*(ye2-ys2) + 0.5);
		viewer[next]->SetBounds(x1,y1,x2-x1,y2-y1);
		if (!init) wnd->Paint();
		// Overides moving component
		if (!init) viewer[next]->Paint();
		// Paint modeless
		int n;
		if (!init) n = GLWindowManager::GetNbWindow();
		if (!init) GLWindowManager::RepaintRange(1,n);
		t1 = (double)SDL_GetTicks()/1000.0;
	}

	} else {
		wnd->Paint();
		viewer[0]->Paint();
		int n = GLWindowManager::GetNbWindow();
		GLWindowManager::RepaintRange(1,n);
	}

	modeSolo = !modeSolo;
	SelectViewer(next);

}


BOOL SynRad::AskToSave() {
	if (!changedSinceSave) return TRUE;
	int ret = GLSaveDialog::Display("Save current geometry first?","File not saved",GLDLG_SAVE | GLDLG_DISCARD | GLDLG_CANCEL_S,GLDLG_ICONINFO);
	if (ret==GLDLG_SAVE) {
		FILENAME *fn = GLFileBox::SaveFile(currentDir,worker.GetShortFileName(),"Save File",fileSFilters,nbSFilter);
		if( fn ) {
			GLProgress *progressDlg2 = new GLProgress("Saving file...","Please wait");
			progressDlg2->SetVisible(TRUE);
			progressDlg2->SetProgress(0.0);
			//GLWindowManager::Repaint();
			try {
				worker.SaveGeometry(fn->fullName,progressDlg2);
				changedSinceSave=FALSE;
				UpdateCurrentDir(fn->fullName);
				UpdateTitle();
				AddRecent(fn->fullName);
			} catch (Error &e) {
				char errMsg[512];
				sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
				GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
				RemoveRecent(fn->fullName);
			}
			progressDlg2->SetVisible(FALSE);
			SAFE_DELETE(progressDlg2);
			return TRUE;
		} else return FALSE;
	} else if (ret==GLDLG_DISCARD) return TRUE;
	return FALSE;
}

BOOL SynRad::AskToReset(Worker *work) {
	SynRad *mApp = (SynRad *)theApp;
	if (work==NULL) work=&worker;
	if (work->nbHit>0) {
		int rep = GLMessageBox::Display("This will reset simulation data.","Geometry change",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);
		if( rep == GLDLG_OK ) {
			work->Reset(m_fTime);
			mApp->nbDesStart = 0;
			mApp->nbHitStart = 0;

			//resetSimu->SetEnabled(FALSE);
			if(mApp->profilePlotter) mApp->profilePlotter->Update(m_fTime,TRUE);
			if(mApp->spectrumPlotter) mApp->spectrumPlotter->Update(m_fTime,TRUE);
			if(mApp->texturePlotter) mApp->texturePlotter->Update(m_fTime,TRUE);
			return TRUE;
		} else return FALSE;
	} else return TRUE;
}

void SynRad::QuickPipe() {

	Geometry *geom = worker.GetGeometry();
	char tmp[256];
	double R = 1.0;
	double L = 5.0;
	int    step=5;
	ResetSimulation(FALSE);
	worker.ClearRegions();
	geom->BuildPipe(L,R,0,step);
	worker.nbDesorption = 0;
	sprintf(tmp,"L|R %g",L/R);
	nbDesStart = 0;
	nbHitStart = 0;
	for(int i=0;i<MAX_VIEWER;i++)
		viewer[i]->SetWorker(&worker);
	startSimu->SetEnabled(TRUE);
	ClearFacetParams();
	if( nbFormula==0 ) {
		GLParser *f = new GLParser();
		f->SetExpression("A2/SUMDES");
		f->SetName("Trans. Prob.");
		f->Parse();
		AddFormula(f);
	}
	UpdateStructMenu();
	// Send to sub process
	try { worker.Reload(); } catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	UpdateTitle();
	changedSinceSave=FALSE;
	ResetAutoSaveTimer();
}

void SynRad::BuildPipe(double ratio) {

	char tmp[128];
	Geometry *geom = worker.GetGeometry();

	double R = 1.0;
	double L = ratio * R;
	int    step;

	sprintf(tmp,"100");
	char *nbF = GLInputBox::GetInput(tmp,"Number of facet","Build Pipe");
	if( !nbF ) return;
	if(( sscanf(nbF,"%d",&step)<=0 )||(step<3)) {
		GLMessageBox::Display("Invalid number","Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	ResetSimulation(FALSE);
	worker.ClearRegions();
	geom->BuildPipe(L,R,0,step);
	worker.nbDesorption = 0;
	sprintf(tmp,"L|R %g",L/R);
	nbDesStart = 0;
	nbHitStart = 0;
	for(int i=0;i<MAX_VIEWER;i++)
		viewer[i]->SetWorker(&worker);
	startSimu->SetEnabled(TRUE);
	ClearFacetParams();
	if(profilePlotter) profilePlotter->Refresh();
	if(spectrumPlotter) spectrumPlotter->Refresh();
	if(texturePlotter) texturePlotter->Update(m_fTime,TRUE);
		//if(outgassingMap) outgassingMap->Update(m_fTime,TRUE);
	if(facetDetails) facetDetails->Update();
	if(facetCoordinates) facetCoordinates->UpdateFromSelection();
	if(vertexCoordinates) vertexCoordinates->Update();
	if( nbFormula==0 ) {
		GLParser *f = new GLParser();
		f->SetExpression("A2/SUMDES");
		f->SetName("Trans. Prob.");
		f->Parse();
		AddFormula(f);
	}
	UpdateStructMenu();
	// Send to sub process
	try { worker.Reload(); } catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	UpdateTitle();
	changedSinceSave=FALSE;
	ResetAutoSaveTimer();
}

//-----------------------------------------------------------------------------

void SynRad::UpdateStructMenu() {

	char tmp[128];
	Geometry *geom = worker.GetGeometry();

	structMenu->Clear();
	structMenu->Add("New structure...",MENU_VIEW_NEWSTRUCT);
	structMenu->Add("Delete structure...",MENU_VIEW_DELSTRUCT);
	structMenu->Add(NULL); //Separator
	structMenu->Add("Show all",MENU_VIEW_STRUCTURE,SDLK_F1,CTRL_MODIFIER);
	structMenu->Add("Show previous",MENU_VIEW_PREVSTRUCT,SDLK_F11,CTRL_MODIFIER);
	structMenu->Add("Show next",MENU_VIEW_NEXTSTRUCT,SDLK_F12,CTRL_MODIFIER);
	structMenu->Add(NULL); //Separator

	for(int i=0;i<geom->GetNbStructure();i++) {
		sprintf(tmp,"Show #%d (%s)",i+1,geom->GetStructureName(i));
		if( i<10 )
			structMenu->Add(tmp,MENU_VIEW_STRUCTURE+(i+1),SDLK_F1+i+1,CTRL_MODIFIER);
		else
			structMenu->Add(tmp,MENU_VIEW_STRUCTURE+(i+1));
	}

	structMenu->SetCheck(MENU_VIEW_STRUCTURE+geom->viewStruct+1,TRUE);

	UpdateTitle();
}

//SELECTIONS
//-----------------------------------------------------------------------------

void SynRad::SelectView(int v) {
	viewer[curViewer]->SetCurrentView(views[v]);
}

//-----------------------------------------------------------------------------

void SynRad::SelectSelection(int v) {
	Geometry *geom = worker.GetGeometry();
	geom->SetSelection((&selections[v].selection),&(selections[v].nbSel));
}

//-----------------------------------------------------------------------------
void SynRad::ClearSelectionMenus() {
	memorizeSelectionsMenu->Clear();
	memorizeSelectionsMenu->Add("Add new...",MENU_SELECTION_ADDNEW,SDLK_w,CTRL_MODIFIER);
	memorizeSelectionsMenu->Add(NULL); // Separator
	clearSelectionsMenu->Clear();
	clearSelectionsMenu->Add("Clear All",MENU_SELECTION_CLEARALL);
	clearSelectionsMenu->Add(NULL); // Separator
	selectionsMenu->Clear();
}

void SynRad::RebuildSelectionMenus() {
	ClearSelectionMenus();
	for(int i=0;i<nbSelection;i++){
		if (i<=8) {
			selectionsMenu->Add(selections[i].name,MENU_SELECTION_SELECTIONS+i,SDLK_1+i,ALT_MODIFIER);
		} else {
			selectionsMenu->Add(selections[i].name,MENU_SELECTION_SELECTIONS+i); //no place for ALT+shortcut
		}
		clearSelectionsMenu->Add(selections[i].name,MENU_SELECTION_CLEARSELECTIONS+i);
		memorizeSelectionsMenu->Add(selections[i].name,MENU_SELECTION_MEMORIZESELECTIONS+i);
	}
}

void SynRad::AddSelection(char *selectionName,ASELECTION s) {

	if(nbSelection<MAX_SELECTION) {
		selections[nbSelection] = s;
		selections[nbSelection].name = _strdup(selectionName);
		nbSelection++;
	} else {
		SAFE_FREE(selections[0].name);
		for(int i=0;i<MAX_SELECTION-1;i++) selections[i] = selections[i+1];
		selections[MAX_SELECTION-1] = s;
		selections[MAX_SELECTION-1].name = _strdup(selectionName);
	}
	RebuildSelectionMenus();
}

void SynRad::ClearSelection(int idClr) {
	SAFE_FREE(selections[0].name);
	for(int i=0;i<nbSelection-1;i++) selections[i] = selections[i+1];
	nbSelection--;
	RebuildSelectionMenus();
}

void SynRad::ClearAllSelections() {
	for(int i=0;i<nbSelection;i++) SAFE_FREE(selections[i].name);
	nbSelection=0;
	ClearSelectionMenus();
}

void SynRad::OverWriteSelection(int idOvr) {
	Geometry *geom = worker.GetGeometry();
	char *selectionName = GLInputBox::GetInput(selections[idOvr].name,"Selection name","Enter selection name");
	if( !selectionName ) return;

	geom->GetSelection(&(selections[idOvr].selection),&(selections[idOvr].nbSel));
	selections[idOvr].name = _strdup(selectionName);
	RebuildSelectionMenus();
}

void SynRad::AddSelection() {
	Geometry *geom = worker.GetGeometry();
	char tmp[32];
	sprintf(tmp,"Selection #%d",nbSelection+1);
	char *selectionName = GLInputBox::GetInput(tmp,"Selection name","Enter selection name");
	if( !selectionName ) return;

	if(nbSelection<MAX_SELECTION) {
		geom->GetSelection(&(selections[nbSelection].selection),&(selections[nbSelection].nbSel));
		selections[nbSelection].name = _strdup(selectionName);
		nbSelection++;
	} else {
		SAFE_FREE(selections[0].name);
		for(int i=0;i<MAX_SELECTION-1;i++) selections[i] = selections[i+1];
		geom->GetSelection(&(selections[MAX_SELECTION-1].selection),&(selections[MAX_SELECTION-1].nbSel));
		selections[MAX_SELECTION-1].name = _strdup(selectionName);
	}
	RebuildSelectionMenus();
}

//VIEWS
//-----------------------------------------------------------------------------
void SynRad::ClearViewMenus() {
	memorizeViewsMenu->Clear();
	memorizeViewsMenu->Add("Add new...",MENU_VIEW_ADDNEW,SDLK_q,CTRL_MODIFIER);
	memorizeViewsMenu->Add(NULL); // Separator
	clearViewsMenu->Clear();
	clearViewsMenu->Add("Clear All",MENU_VIEW_CLEARALL);
	clearViewsMenu->Add(NULL); // Separator
	viewsMenu->Clear();
}

void SynRad::RebuildViewMenus() {
	ClearViewMenus();
	for(int i=0;i<nbView;i++){
		int id=i;
		if (nbView>=10) id=i-nbView+8;
		if (id>=0 && id<=8) {
			viewsMenu->Add(views[i].name,MENU_VIEW_VIEWS+i,SDLK_F1+id,ALT_MODIFIER);
		} else {
			viewsMenu->Add(views[i].name,MENU_VIEW_VIEWS+i);
		}
		clearViewsMenu->Add(views[i].name,MENU_VIEW_CLEARVIEWS+i);
		memorizeViewsMenu->Add(views[i].name,MENU_VIEW_MEMORIZEVIEWS+i);
	}
}

void SynRad::RebuildPARMenus() {
	PARloadToMenu->Clear();
	PARremoveMenu->Clear();
	for(int i=0;i<(int)worker.regions.size();i++){
		char tmp[256];
		sprintf(tmp,"Region %d",i+1);
		if (worker.regions[i].fileName.length()>0)
			sprintf(tmp,"%s (%s)",tmp,worker.regions[i].fileName.c_str());
		PARloadToMenu->Add(tmp,MENU_REGIONS_LOADTO+i);
		PARremoveMenu->Add(tmp,MENU_REGIONS_REMOVE+i);
	}
}

void SynRad::AddView(char *viewName,AVIEW v) {

	if(nbView<MAX_VIEW) {
		views[nbView] = v;
		views[nbView].name = _strdup(viewName);
		nbView++;
	} else {
		SAFE_FREE(views[0].name);
		for(int i=0;i<MAX_VIEW-1;i++) views[i] = views[i+1];
		views[MAX_VIEW-1] = v;
		views[MAX_VIEW-1].name = _strdup(viewName);
	}
	RebuildViewMenus();
}

void SynRad::ClearView(int idClr) {
	SAFE_FREE(views[0].name);
	for(int i=0;i<nbView-1;i++) views[i] = views[i+1];
	nbView--;
	RebuildViewMenus();
}

void SynRad::ClearAllViews() {
	for(int i=0;i<nbView;i++) SAFE_FREE(views[i].name);
	nbView=0;
	ClearViewMenus();
}

void SynRad::OverWriteView(int idOvr) {
	Geometry *geom = worker.GetGeometry();
	char *viewName = GLInputBox::GetInput(views[idOvr].name,"View name","Enter view name");
	if( !viewName ) return;

	views[idOvr] = viewer[curViewer]->GetCurrentView();
	views[idOvr].name = _strdup(viewName);
	RebuildViewMenus();
}

void SynRad::AddView() {
	Geometry *geom = worker.GetGeometry();
	char tmp[32];
	sprintf(tmp,"View #%d",nbView+1);
	char *viewName = GLInputBox::GetInput(tmp,"View name","Enter view name");
	if( !viewName ) return;

	if(nbView<MAX_VIEW) {
		views[nbView] = viewer[curViewer]->GetCurrentView();
		views[nbView].name = _strdup(viewName);
		nbView++;
	} else {
		SAFE_FREE(views[0].name);
		for(int i=0;i<MAX_VIEW-1;i++) views[i] = views[i+1];
		views[MAX_VIEW-1] = viewer[curViewer]->GetCurrentView();
		views[MAX_VIEW-1].name = _strdup(viewName);
	}
	RebuildViewMenus();
}

//-----------------------------------------------------------------------------

void SynRad::RemoveRecent(char *fileName) {

	if( !fileName ) return;

	BOOL found = FALSE;
	int i = 0;
	while(!found && i<nbRecent) {
		found = strcmp(fileName,recents[i])==0;
		if(!found) i++;
	}
	if(!found) return;

	SAFE_FREE(recents[i]);
	for(int j=i;j<nbRecent-1;j++)
		recents[j] = recents[j+1];
	nbRecent--;

	// Update menu
	GLMenu *m = menu->GetSubMenu("File")->GetSubMenu("Load recent");
	m->Clear();
	for(i=nbRecent-1;i>=0;i--)
		m->Add(recents[i],MENU_FILE_LOADRECENT+i);
	SaveConfig();
}

//-----------------------------------------------------------------------------

void SynRad::AddRecent(char *fileName) {

	// Check if already exists
	BOOL found = FALSE;
	int i = 0;
	while(!found && i<nbRecent) {
		found = strcmp(fileName,recents[i])==0;
		if(!found) i++;
	}
	if(found) {
		for (int j=i;j<nbRecent-1;j++) {
			recents[j]=recents[j+1];
		}
		recents[nbRecent-1]=_strdup(fileName);
		// Update menu
		GLMenu *m = menu->GetSubMenu("File")->GetSubMenu("Load recent");
		m->Clear();
		for(int i=nbRecent-1;i>=0;i--)
			m->Add(recents[i],MENU_FILE_LOADRECENT+i);
		SaveConfig();
		return;
	}

	// Add the new recent file
	if(nbRecent<MAX_RECENT) {
		recents[nbRecent] = _strdup(fileName);
		nbRecent++;
	} else {
		// Shift
		SAFE_FREE(recents[0]);
		for(int i=0;i<MAX_RECENT-1;i++)
			recents[i] = recents[i+1];
		recents[MAX_RECENT-1]=_strdup(fileName);
	}

	// Update menu
	GLMenu *m = menu->GetSubMenu("File")->GetSubMenu("Load recent");
	m->Clear();
	for(int i=nbRecent-1;i>=0;i--)
		m->Add(recents[i],MENU_FILE_LOADRECENT+i);
	SaveConfig();
}

void SynRad::RemoveRecentPAR(char *fileName) {

	if( !fileName ) return;

	BOOL found = FALSE;
	int i = 0;
	while(!found && i<nbRecentPAR) {
		found = strcmp(fileName,recentPARs[i])==0;
		if(!found) i++;
	}
	if(!found) return;

	SAFE_FREE(recentPARs[i]);
	for(int j=i;j<nbRecentPAR-1;j++)
		recentPARs[j] = recentPARs[j+1];
	nbRecentPAR--;

	// Update menu
	GLMenu *m = menu->GetSubMenu("Regions")->GetSubMenu("Load recent");
	m->Clear();
	for(i=nbRecentPAR-1;i>=0;i--)
		m->Add(recentPARs[i],MENU_REGIONS_LOADRECENT+i);
	SaveConfig();
}

//-----------------------------------------------------------------------------

void SynRad::AddRecentPAR(char *fileName) {

	// Check if already exists
	BOOL found = FALSE;
	int i = 0;
	while(!found && i<nbRecentPAR) {
		found = strcmp(fileName,recentPARs[i])==0;
		if(!found) i++;
	}
	if(found) {
		for (int j=i;j<nbRecentPAR-1;j++) {
			recentPARs[j]=recentPARs[j+1];
		}
		recentPARs[nbRecentPAR-1]=_strdup(fileName);
		// Update menu
		GLMenu *m = menu->GetSubMenu("Regions")->GetSubMenu("Load recent");
		m->Clear();
		for(int i=nbRecentPAR-1;i>=0;i--)
			m->Add(recentPARs[i],MENU_REGIONS_LOADRECENT+i);
		SaveConfig();
		return;
	}

	// Add the new recent file
	if(nbRecentPAR<MAX_RECENT) {
		recentPARs[nbRecentPAR] = _strdup(fileName);
		nbRecentPAR++;
	} else {
		// Shift
		SAFE_FREE(recentPARs[0]);
		for(int i=0;i<MAX_RECENT-1;i++)
			recentPARs[i] = recentPARs[i+1];
		recentPARs[MAX_RECENT-1]=_strdup(fileName);
	}

	// Update menu
	GLMenu *m = menu->GetSubMenu("Regions")->GetSubMenu("Load recent");
	m->Clear();
	for(int i=nbRecentPAR-1;i>=0;i--)
		m->Add(recentPARs[i],MENU_REGIONS_LOADRECENT+i);
	SaveConfig();
}


//-----------------------------------------------------------------------------

void SynRad::LoadConfig() {

	FileReader *f = NULL;
	char *w;
	nbRecent = 0;

	try {

		f = new FileReader("synrad.cfg");
		Geometry *geom = worker.GetGeometry();

		f->ReadKeyword("showRules");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showRule = f->ReadInt();
		f->ReadKeyword("showNormals");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showNormal = f->ReadInt();
		f->ReadKeyword("showUV");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showUV = f->ReadInt();
		f->ReadKeyword("showLines");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showLine = f->ReadInt();
		f->ReadKeyword("showLeaks");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showLeak = f->ReadInt();
		f->ReadKeyword("showHits");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showHit = f->ReadInt();
		f->ReadKeyword("showVolume");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showVolume = f->ReadInt();
		f->ReadKeyword("showTexture");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showTexture = f->ReadInt();
		f->ReadKeyword("showFilter");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showFilter = f->ReadInt();
		f->ReadKeyword("showIndices");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showIndex = f->ReadInt();
		f->ReadKeyword("showVertices");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showVertex = f->ReadInt();
		f->ReadKeyword("showMode");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showBack = f->ReadInt();
		f->ReadKeyword("showMesh");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showMesh = f->ReadInt();
		f->ReadKeyword("showHidden");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showHidden = f->ReadInt();
		f->ReadKeyword("showHiddenVertex");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showHiddenVertex = f->ReadInt();
		f->ReadKeyword("texColormap");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showColormap = f->ReadInt();
		f->ReadKeyword("translation");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->transStep = f->ReadDouble();
		f->ReadKeyword("dispNumLines");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->dispNumHits = f->ReadInt();
		f->ReadKeyword("dispNumLeaks");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->dispNumLeaks = f->ReadInt();
		f->ReadKeyword("dispNumTraj");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->dispNumTraj = f->ReadInt();
		f->ReadKeyword("angle");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->angleStep = f->ReadDouble();
		f->ReadKeyword("autoScale");f->ReadKeyword(":");
		geom->texAutoScale = f->ReadInt();
		f->ReadKeyword("texMin_MC");f->ReadKeyword(":");
		geom->texMin_MC = (llong)f->ReadDouble();
		f->ReadKeyword("texMax_MC");f->ReadKeyword(":");
		geom->texMax_MC = (llong)f->ReadDouble();
		f->ReadKeyword("texMin_flux");f->ReadKeyword(":");
		geom->texMin_flux = f->ReadDouble();
		f->ReadKeyword("texMax_flux");f->ReadKeyword(":");
		geom->texMax_flux = f->ReadDouble();
		f->ReadKeyword("texMin_power");f->ReadKeyword(":");
		geom->texMin_power = f->ReadDouble();
		f->ReadKeyword("texMax_power");f->ReadKeyword(":");
		geom->texMax_power = f->ReadDouble();
		f->ReadKeyword("processNum");f->ReadKeyword(":");
		nbProc = f->ReadInt();
		if(nbProc<=0) nbProc=1;
		f->ReadKeyword("recents");f->ReadKeyword(":");f->ReadKeyword("{");
		w = f->ReadString();
		while( strcmp(w,"}")!=0 && nbRecent<MAX_RECENT)  {
			recents[nbRecent] = _strdup(w);
			nbRecent++;
			w = f->ReadString();
		}
		f->ReadKeyword("recentPARs");f->ReadKeyword(":");f->ReadKeyword("{");
		w = f->ReadString();
		while( strcmp(w,"}")!=0 && nbRecentPAR<MAX_RECENT)  {
			recentPARs[nbRecentPAR] = _strdup(w);
			nbRecentPAR++;
			w = f->ReadString();
		}
		f->ReadKeyword("cdir");f->ReadKeyword(":");
		strcpy(currentDir,f->ReadString());
		f->ReadKeyword("cseldir");f->ReadKeyword(":");
		strcpy(currentSelDir,f->ReadString());
		f->ReadKeyword("autonorme");f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("centernorme");f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		f->ReadKeyword("normeratio");f->ReadKeyword(":");
		geom->SetNormeRatio((float)(f->ReadDouble()));
		f->ReadKeyword("showDirection");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showDir = f->ReadInt();
		f->ReadKeyword("shadeLines");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->shadeLines = f->ReadInt();
		f->ReadKeyword("showTP");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showTP = f->ReadInt();
		f->ReadKeyword("autoSaveFrequency");f->ReadKeyword(":");
		autoSaveFrequency = f->ReadDouble();
		f->ReadKeyword("autoSaveSimuOnly");f->ReadKeyword(":");
		autoSaveSimuOnly = f->ReadInt();
		f->ReadKeyword("checkForUpdates");f->ReadKeyword(":");
		checkForUpdates = f->ReadInt();
		f->ReadKeyword("compressSavedFiles");f->ReadKeyword(":");
		compressSavedFiles = f->ReadInt();

	} catch (Error &err) {
		printf("Warning, load config file (one or more feature not supported) %s\n",err.GetMsg());
	}

	SAFE_DELETE(f);


}

//-----------------------------------------------------------------------------

#define WRITEI(name,var) {             \
	f->Write(name);                      \
	f->Write(":");                       \
	for(int i=0;i<MAX_VIEWER;i++)        \
	f->WriteInt(viewer[i]->var," ");   \
	f->Write("\n");                      \
}                                      \

#define WRITED(name,var) {             \
	f->Write(name);                      \
	f->Write(":");                       \
	for(int i=0;i<MAX_VIEWER;i++)        \
	f->WriteDouble(viewer[i]->var," ");\
	f->Write("\n");                      \
}

//----------------------------------------------------------------------------
void SynRad::UpdateViewerFlags() {
	viewer[curViewer]->showNormal=showNormal->IsChecked();
	viewer[curViewer]->showRule=showRule->IsChecked();
	viewer[curViewer]->showUV=showUV->IsChecked();
	viewer[curViewer]->showLeak=showLeak->IsChecked();
	viewer[curViewer]->showHit=showHit->IsChecked();
	viewer[curViewer]->showLine=showLine->IsChecked();
	viewer[curViewer]->showVolume=showVolume->IsChecked();
	viewer[curViewer]->showTexture=showTexture->IsChecked();
	viewer[curViewer]->showFilter=showFilter->IsChecked();
	viewer[curViewer]->showVertex=showVertex->IsChecked();
	viewer[curViewer]->showIndex=showIndex->IsChecked();
	//worker.Update(0.0);
}

//-----------------------------------------------------------------------------

void SynRad::SaveConfig() {

	FileWriter *f = NULL;

	try {

		f = new FileWriter("synrad.cfg");
		Geometry *geom = worker.GetGeometry();

		// Save flags
		WRITEI("showRules",showRule);
		WRITEI("showNormals",showNormal);
		WRITEI("showUV",showUV);
		WRITEI("showLines",showLine);
		WRITEI("showLeaks",showLeak);
		WRITEI("showHits",showHit);
		WRITEI("showVolume",showVolume);
		WRITEI("showTexture",showTexture);
		WRITEI("showFilter",showFilter);
		WRITEI("showIndices",showIndex);
		WRITEI("showVertices",showVertex);
		WRITEI("showMode",showBack);
		WRITEI("showMesh",showMesh);
		WRITEI("showHidden",showHidden);
		WRITEI("showHiddenVertex",showHiddenVertex);
		WRITEI("texColormap",showColormap);
		WRITED("translation",transStep);
		WRITEI("dispNumLines",dispNumHits);
		WRITEI("dispNumLeaks",dispNumLeaks);
		WRITEI("dispNumTraj",dispNumTraj);
		WRITED("angle",angleStep);
		f->Write("autoScale:");f->WriteInt(geom->texAutoScale,"\n");
		f->Write("texMin_MC:");f->WriteLLong(geom->texMin_MC,"\n");
		f->Write("texMax_MC:");f->WriteLLong(geom->texMax_MC,"\n");
		f->Write("texMin_flux:");f->WriteDouble(geom->texMin_flux,"\n");
		f->Write("texMax_flux:");f->WriteDouble(geom->texMax_flux,"\n");
		f->Write("texMin_power:");f->WriteDouble(geom->texMin_power,"\n");
		f->Write("texMax_power:");f->WriteDouble(geom->texMax_power,"\n");
		f->Write("processNum:");f->WriteInt(worker.GetProcNumber(),"\n");
		f->Write("recents:{\n");
		for(int i=0;i<nbRecent;i++) {
			f->Write("\"");
			f->Write(recents[i]);
			f->Write("\"\n");
		}
		f->Write("}\n");
		f->Write("recentPARs:{\n");
		for(int i=0;i<nbRecentPAR;i++) {
			f->Write("\"");
			f->Write(recentPARs[i]);
			f->Write("\"\n");
		}
		f->Write("}\n");
		f->Write("cdir:\"");f->Write(currentDir);f->Write("\"\n");
		f->Write("cseldir:\"");f->Write(currentSelDir);f->Write("\"\n");
		f->Write("autonorme:");f->WriteInt(geom->GetAutoNorme(),"\n");
		f->Write("centernorme:");f->WriteInt(geom->GetCenterNorme(),"\n");
		f->Write("normeratio:");f->WriteDouble((double)(geom->GetNormeRatio()),"\n");
		WRITEI("showDirection",showDir);f->Write("\n");
		WRITEI("shadeLines",showDir);f->Write("\n");
		WRITEI("showTP",showDir);f->Write("\n");
		f->Write("autoSaveFrequency:");f->WriteDouble(autoSaveFrequency,"\n");
		f->Write("autoSaveSimuOnly:");f->WriteInt(autoSaveSimuOnly,"\n");
		f->Write("checkForUpdates:");f->WriteInt(checkForUpdates,"\n");
		f->Write("compressSavedFiles:");f->WriteInt(compressSavedFiles,"\n");

	} catch (Error &err) {
		printf("Warning, failed to save config file %s\n",err.GetMsg());
	}

	SAFE_DELETE(f);

}

void SynRad::CrashHandler(Error *e) {
	char tmp[1024];
	sprintf(tmp,"Well, that's emberassing. Synrad crashed and will exit now.\nBefore that, an autosave will be attempted.\nHere is the error info:\n\n%s",(char *)e->GetMsg());
	GLMessageBox::Display(tmp,"Main crash handler",GLDLG_OK,GLDGL_ICONDEAD);
	try {
		theApp->AutoSave(TRUE); //crashSave
		GLMessageBox::Display("Good news, autosave worked!","Main crash handler",GLDLG_OK,GLDGL_ICONDEAD);
	} catch(Error &e) {
		e.GetMsg();
		GLMessageBox::Display("Sorry, I couldn't even autosave.","Main crash handler",GLDLG_OK,GLDGL_ICONDEAD);
	}
}

void SynRad::AddStruct() {
	Geometry *geom = worker.GetGeometry();
	char tmp[32];
	sprintf(tmp,"Structure #%d",geom->GetNbStructure()+1);
	char *structName = GLInputBox::GetInput(tmp,"Structure name","Enter name of new structure");
	if( !structName ) return;
	geom->AddStruct(structName);
	// Send to sub process
	try { worker.Reload(); } catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
}

void SynRad::DeleteStruct() {
	Geometry *geom = worker.GetGeometry();
	char *structNum = GLInputBox::GetInput("","Structure number","Number of structure to delete:");
	if( !structNum ) return;
	int structNumInt;
	if ( !sscanf(structNum,"%d",&structNumInt)) {
		GLMessageBox::Display("Invalid structure number");
		return;
	}
	if (structNumInt<1 || structNumInt>geom->GetNbStructure()) {
		GLMessageBox::Display("Invalid structure number");
		return;
	}
	BOOL hasFacets=FALSE;
	for (int i=0;i<geom->GetNbFacet() && !hasFacets;i++) {
		if (geom->GetFacet(i)->sh.superIdx==(structNumInt-1)) hasFacets=TRUE;
	}
	if (hasFacets) {
		int rep = GLMessageBox::Display("This structure has facets. They will be deleted with the structure.","Structure delete",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);
		if( rep != GLDLG_OK ) return;
	}
	if (!AskToReset()) return;
	geom->DelStruct(structNumInt-1);
	// Send to sub process
	try { worker.Reload(); } catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
}

void SynRad::DisplayCollapseDialog() {
	Geometry *geom = worker.GetGeometry();
	if( !collapseSettings ) collapseSettings = new CollapseSettings();
	collapseSettings->SetGeometry(geom,&worker);
	collapseSettings->SetVisible(TRUE);
}

void SynRad::LoadParam(char *fName,int position) {

	if (!worker.GetGeometry()->IsLoaded()) {
		GLMessageBox::Display("You have to load a geometry before adding regions.","Add magnetic region",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	char fullName[512];
	char shortName[512];
	strcpy(fullName,"");

	if(fName==NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir,NULL,"Open File",fileParFilters,nbParFilter);
		if( fn )
			strcpy(fullName,fn->fullName);
	} else {
		strcpy(fullName,fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Preparing to load file...","Please wait");
	progressDlg2->SetVisible(TRUE);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if(strlen(fullName)==0) {
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		return;
	}


	char *lPart = strrchr(fullName,'\\');
	if(lPart) strcpy(shortName,lPart+1);
	else strcpy(shortName,fullName);

	try {
		AddRecentPAR(fullName);
		worker.AddRegion(fullName,position);
		if (regionInfo) regionInfo->Update();
	} catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),shortName);
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		RemoveRecentPAR(fullName);
	}
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
	changedSinceSave=FALSE;
	RebuildPARMenus();
}

void SynRad::ClearTraj(){
	worker.ClearRegions();
	if (regionInfo) {
		regionInfo->SetVisible(FALSE);
		SAFE_DELETE(regionInfo);
	}
	RebuildPARMenus();
}

void SynRad::RemoveRegion(int index){
	worker.nbTrajPoints-=(int)worker.regions[index].Points.size();
	worker.regions.erase(worker.regions.begin() + index);
	if (regionInfo) regionInfo->Update();
	RebuildPARMenus();
	worker.GetGeometry()->InitializeGeometry(); //recalculate bounding box
	worker.Reload();
}


void SynRad::RenumberSelections(int startFacetId){
	for (int i=0;i<nbSelection;i++) {
		BOOL found=FALSE;
		for (int j=0;j<selections[i].nbSel && !found;j++) { //remove from selection
			if (selections[i].selection[j]==startFacetId) {

				for (int k=j;k<(selections[i].nbSel-1);k++)
					selections[i].selection[k]=selections[i].selection[k+1];
				selections[i].nbSel--;
				if (selections[i].nbSel==0) ClearSelection(i);
				found=TRUE;
			} 
		}
		for (int j=0;j<selections[i].nbSel;j++) { //renumber selection
			if (selections[i].selection[j]>startFacetId) {
				//decrease number by 1
				selections[i].selection[j]--; 
			}
		}
	}
}

void SynRad::NewRegion(){
	if (!worker.GetGeometry()->IsLoaded()) {
		GLMessageBox::Display("You have to load a geometry before adding regions.","Add magnetic region",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	if (worker.running) worker.Stop_Public();
	FILENAME *fn = GLFileBox::SaveFile(NULL,NULL,"Save Region","param files\0*.param\0All files\0*.*\0",2);
	if (!fn || !fn->fullName) return;
	if (!EndsWithParam(fn->fullName))
		sprintf(fn->fullName,"%s.param",fn->fullName); //append .param extension
	try {	
		Region_full newreg;
		newreg.fileName.assign(fn->fullName);
		//worker.regions.push_back(newreg);
		FileWriter *file=new FileWriter(fn->fullName);
		newreg.SaveParam(file);
		SAFE_DELETE(file);
		AddRecentPAR(fn->fullName);
		worker.AddRegion(fn->fullName);
		RebuildPARMenus();
	} catch(Error &e) {
		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		RemoveRecentPAR(fn->fullName);
		return;
	}
	if( regionInfo ) regionInfo->Update();
	regionEditor2 = new RegionEditor();
	regionEditor2->Display(&worker,(int)worker.regions.size()-1);
	regionEditor2->DoModal();
	SAFE_DELETE(regionEditor2);
}

/*BOOL EndsWithParam(const char* s)
{
int ret = 0;

if (s != NULL)
{
size_t size = strlen(s);

if (size >= 6 &&
s[size-6] == '.' &&
s[size-5] == 'p' &&
s[size-4] == 'a' &&
s[size-3] == 'r' &&
s[size-2] == 'a' &&
s[size-1] == 'm')
{
ret = 1;
}
}

return ret;
}*/