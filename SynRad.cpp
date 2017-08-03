/*
File:        SynRad.cpp
Description: Main application class (GUI management)
Program:     SynRad+
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
#include "GLApp/GLMenuBar.h"
#include "RecoveryDialog.h"
#include "Facet.h"
#include "SynradGeometry.h"
#include "GLApp\MathTools.h"
//for Remainder
#include "direct.h"
#include <vector>
#include <string>
#include <io.h>
#include <numeric> //std::iota
//#include <winsparkle.h>

#include "FacetCoordinates.h"
#include "SmartSelection.h"
#include "VertexCoordinates.h"
#include "FormulaEditor.h"

static const char *fileLFilters = "All SynRad supported files\0*.xml;*.zip;*.txt;*.syn;*.syn7z;*.geo;*.geo7z;*.str;*.stl;*.ase\0All files\0*.*\0";
const char *fileSFilters = "SYN files\0*.syn;*.syn7z;\0GEO files\0*.geo;*.geo7z;\0Text files\0*.txt\0All files\0*.*\0";
//static const char *fileSelFilters = "Selection files\0*.sel\0All files\0*.*\0";
//static const char *fileTexFilters = "Text files\0*.txt\0Texture files\0*.tex\0All files\0*.*\0";
static const char *fileParFilters = "param files\0*.param\0PAR files\0*.par\0All files\0*.*\0";
static const char *fileDesFilters = "DES files\0*.des\0All files\0*.*\0";

int   cSize = 5;
int   cWidth[] = { 30, 56, 50, 50, 50 };
char *cName[] = { "#", "Hits", "Flux", "Power", "Abs" };

int appVersion = 1413;
std::string appId = "Synrad";
#ifdef _DEBUG
std::string appName = "SynRad+ development version (Compiled " __DATE__ " " __TIME__ ") DEBUG MODE";
#else
std::string appName = "Synrad+ 1.4.13 (" __DATE__ ")";
#endif

std::vector<string> formulaPrefixes = { "H","A","F","P","AR","h","a","f","p","ar","," };
std::string formulaSyntax =
R"(MC Variables: An (Absorption on facet n), Hn (Hit on facet n)
Fn (Flux absorbed on facet n), Pn (Power absorbed on facet n)
SUMABS (total absorbed), SUMDES (total desorbed), SUMHIT (total hit)
SUMFLUX (total gen. flux), SUMPOWER (total gen. power), SCANS (no. of scans)

Area variables: ARn (Area of facet n), ABSAR (total absorption area)

Math functions: sin(), cos(), tan(), sinh(), cosh(), tanh(),
                   asin(), acos(), atan(), exp(), ln(), pow(x,y)
                   log2(), log10(), inv(), sqrt(), abs()

Utils functions: ci95(p,N) 95% confidence interval (p=prob,N=count)
                  SUM(prefix,i,j) sum variables ex: sum(F,1,10)=F1+F2+...+F10
                  SUM(prefix,Si) sum of variables on selection group i  ex: sum(F,S1)=all flux on group 1
                  SUM(prefix,SEL) calculates the sum of hits on the current selection. (Works with H,A,F,P,AR)

Constants: Kb (Boltzmann's constant), R (Gas constant)
             Na (Avogadro's number [mol]), PI

Expression example:
  (A1+A45)/(D23+D12)
  sqrt(A1^2+A45^2)*DESAR/SUMDES
)";
int formulaSyntaxHeight = 320;

float m_fTime;
SynRad *mApp;

//Menu elements, Synrad specific
#define MENU_FILE_EXPORT_DESORP 140

#define MENU_FILE_EXPORTTEXTURE_AREA 151
#define MENU_FILE_EXPORTTEXTURE_MCHITS 152
#define MENU_FILE_EXPORTTEXTURE_FLUX 153
#define MENU_FILE_EXPORTTEXTURE_POWER 154
#define MENU_FILE_EXPORTTEXTURE_FLUXPERAREA 155
#define MENU_FILE_EXPORTTEXTURE_POWERPERAREA 156
#define MENU_FILE_EXPORTTEXTURE_ANSYS_POWER 157

#define MENU_FILE_EXPORTTEXTURE_AREA_COORD 171
#define MENU_FILE_EXPORTTEXTURE_MCHITS_COORD 172
#define MENU_FILE_EXPORTTEXTURE_FLUX_COORD 173
#define MENU_FILE_EXPORTTEXTURE_POWER_COORD 174
#define MENU_FILE_EXPORTTEXTURE_FLUXPERAREA_COORD 175
#define MENU_FILE_EXPORTTEXTURE_POWERPERAREA_COORD 176
#define MENU_FILE_EXPORTTEXTURE_ANSYS_POWER_COORD 177

#define MENU_REGIONS_NEW        901
#define MENU_REGIONS_LOADPAR      902
#define MENU_REGIONS_CLEARALL    903
#define MENU_REGIONS_REGIONINFO 904

#define MENU_REGIONS_LOADRECENT		910
#define MENU_REGIONS_LOADTO			930
#define MENU_REGIONS_REMOVE			950
#define MENU_REGIONS_SHOWHITS		970
#define MENU_REGIONS_SHOWNONE		990
#define MENU_REGIONS_SHOWALL		991

#define MENU_TOOLS_SPECTRUMPLOTTER 403

#define MENU_FACET_MESH        360

#define MENU_FACET_SELECTSPECTRUM 361

//-----------------------------------------------------------------------------
// Name: WinMain()
// Desc: Entry point to the program. Initializes everything, and goes into a
//       message-processing loop. Idle time is used to render the scene.
//-----------------------------------------------------------------------------

INT WINAPI WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, INT)
{

	SynRad *mApp = new SynRad();

	if (!mApp->Create(1024, 800, false)) {
		char *logs = GLToolkit::GetLogs();
#ifdef WIN
		if (logs) MessageBox(NULL, logs, "Synrad [Fatal error]", MB_OK);
#else
		if (logs) {
			printf("Synrad [Fatal error]\n");
			printf(logs);
		}
#endif
		SAFE_FREE(logs);
		delete mApp;
		return -1;
	}
	try {
		mApp->Run();
	}
	catch (Error &e) {
		mApp->CrashHandler(&e);
	}
	delete mApp;

	return 0;
}



//-----------------------------------------------------------------------------
// Name: SynRad()
// Desc: Application constructor. Sets default attributes for the app.
//-----------------------------------------------------------------------------
SynRad::SynRad()
{
	mApp = this; //to refer to the app as extern variable

	nbRecentPAR = 0;

	//Different Synrad implementation:
	facetMesh = NULL;
	facetDetails = NULL;
	viewer3DSettings = NULL;
	textureSettings = NULL;
	formulaEditor = NULL;
	globalSettings = NULL;
	profilePlotter = NULL;
	texturePlotter = NULL;

	//Synrad only:
	regionInfo = NULL;
	exportDesorption = NULL;
	trajectoryDetails = NULL;
	spectrumPlotter = NULL;
	regionEditor = NULL;

	materialPaths = std::vector<std::string>();
}

//-----------------------------------------------------------------------------
// Name: OneTimeSceneInit()
// Desc: Called during initial app startup, this function performs all the
//       permanent initialization.
//-----------------------------------------------------------------------------
int SynRad::OneTimeSceneInit()
{

	OneTimeSceneInit_shared();

	menu->GetSubMenu("File")->Add("Export DES file (deprecated)", MENU_FILE_EXPORT_DESORP);

	menu->GetSubMenu("File")->Add("Export selected textures");
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Facet by facet");
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Element Area", MENU_FILE_EXPORTTEXTURE_AREA);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("MC hits", MENU_FILE_EXPORTTEXTURE_MCHITS);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Flux(ph/sec)", MENU_FILE_EXPORTTEXTURE_FLUX);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Power(W)", MENU_FILE_EXPORTTEXTURE_POWER);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Flux density (ph/sec/cm\262)", MENU_FILE_EXPORTTEXTURE_FLUXPERAREA);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Power density (W/mm\262)", MENU_FILE_EXPORTTEXTURE_POWERPERAREA);

	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("By X,Y,Z coordinates");
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Element Area", MENU_FILE_EXPORTTEXTURE_AREA_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("MC hits", MENU_FILE_EXPORTTEXTURE_MCHITS_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Flux(ph/sec)", MENU_FILE_EXPORTTEXTURE_FLUX_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Power(W)", MENU_FILE_EXPORTTEXTURE_POWER_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Flux density (ph/sec/cm\262)", MENU_FILE_EXPORTTEXTURE_FLUXPERAREA_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Power density (W/mm\262)", MENU_FILE_EXPORTTEXTURE_POWERPERAREA_COORD);

	menu->GetSubMenu("File")->Add(NULL); // Separator
	menu->GetSubMenu("File")->Add("E&xit", MENU_FILE_EXIT);  //Moved here from OnetimeSceneinit_shared to assert it's the last menu item

	menu->Add("Regions");
	menu->GetSubMenu("Regions")->Add("New...", MENU_REGIONS_NEW);
	menu->GetSubMenu("Regions")->Add("Load...", MENU_REGIONS_LOADPAR);
	menu->GetSubMenu("Regions")->Add("Load recent", MENU_REGIONS_LOADRECENT); //Will add recent files after loading config
	menu->GetSubMenu("Regions")->Add("Remove all", MENU_REGIONS_CLEARALL);
	menu->GetSubMenu("Regions")->Add(NULL); // Separator

	PARloadToMenu = menu->GetSubMenu("Regions")->Add("Load to");
	PARremoveMenu = menu->GetSubMenu("Regions")->Add("Remove");
	menu->GetSubMenu("Regions")->Add(NULL); // Separator
	ShowHitsMenu = menu->GetSubMenu("Regions")->Add("Show lines, hits");
	menu->GetSubMenu("Regions")->Add(NULL); // Separator
	menu->GetSubMenu("Regions")->Add("Region Info ...", MENU_REGIONS_REGIONINFO);

	menu->GetSubMenu("Selection")->Add(NULL); // Separator
	menu->GetSubMenu("Selection")->Add("Select volatile facet", MENU_FACET_SELECTVOL);
	menu->GetSubMenu("Selection")->Add("Select Spectrum", MENU_FACET_SELECTSPECTRUM);

	menu->GetSubMenu("Tools")->Add("Spectrum Plotter ...", MENU_TOOLS_SPECTRUMPLOTTER);

	showFilter = new GLToggle(0, "Filtering");
	togglePanel->Add(showFilter);

	showMoreBtn = new GLButton(0, "<< View");
	togglePanel->Add(showMoreBtn);

	/*
	shortcutPanel = new GLTitledPanel("Shortcuts");
	shortcutPanel->SetClosable(true);
	shortcutPanel->Close();
	Add(shortcutPanel);

	profilePlotterBtn = new GLButton(0, "Profile pl.");
	shortcutPanel->Add(profilePlotterBtn);

	texturePlotterBtn = new GLButton(0, "Texture pl.");
	shortcutPanel->Add(texturePlotterBtn);

	textureScalingBtn = new GLButton(0, "Tex.scaling");
	shortcutPanel->Add(textureScalingBtn);
	*/

	modeLabel = new GLLabel("Mode");
	simuPanel->Add(modeLabel);

	modeCombo = new GLCombo(0);
	modeCombo->SetEditable(true);
	modeCombo->SetSize(2);
	modeCombo->SetValueAt(0, "Fluxwise");
	modeCombo->SetValueAt(1, "Powerwise");
	modeCombo->SetSelectedIndex(worker.generation_mode);
	simuPanel->Add(modeCombo);

	/*globalSettingsBtn = new GLButton(0, "<< Sim");
	simuPanel->Add(globalSettingsBtn);*/

	doseLabel = new GLLabel("Dose");
	simuPanel->Add(doseLabel);

	doseNumber = new GLTextField(0, NULL);
	doseNumber->SetEditable(false);
	simuPanel->Add(doseNumber);

	//Reflection materials
	//Find material files in param directory
	intptr_t file;
	_finddata_t filedata;
	file = _findfirst("param\\Materials\\*.csv", &filedata);
	if (file != -1)
	{
		do
		{
			materialPaths.push_back(filedata.name);
		} while (_findnext(file, &filedata) == 0);
	}
	_findclose(file);

	facetRLabel = new GLLabel("Refl:");
	facetPanel->Add(facetRLabel);
	facetReflType = new GLCombo(0);
	facetReflType->SetSize((int)materialPaths.size() + 2);
	facetReflType->SetValueAt(0, "Diffuse, Sticking->", REF_DIFFUSE);
	facetReflType->SetValueAt(1, "Mirror, Sticking->", REF_MIRROR);
	for (int i = 0; i < (int)materialPaths.size(); i++) {
		size_t lastindex = materialPaths[i].find_last_of("."); //cut extension
		facetReflType->SetValueAt(i + 2, materialPaths[i].substr(0, lastindex).c_str(), REF_MATERIAL + i);
	}
	facetPanel->Add(facetReflType);

	facetSticking = new GLTextField(0, "");
	facetPanel->Add(facetSticking);

	facetDoScattering = new GLToggle(0, "Rough surface scattering");
	facetDoScattering->SetState(0);
	facetPanel->Add(facetDoScattering);

	facetRMSroughnessLabel = new GLLabel("sigma (nm):");
	facetPanel->Add(facetRMSroughnessLabel);
	facetRMSroughness = new GLTextField(0, NULL);
	facetRMSroughness->SetEditable(false);
	facetPanel->Add(facetRMSroughness);

	facetAutoCorrLengthLabel = new GLLabel("T (nm):");
	facetPanel->Add(facetAutoCorrLengthLabel);
	facetAutoCorrLength = new GLTextField(0, NULL);
	facetAutoCorrLength->SetEditable(false);
	facetPanel->Add(facetAutoCorrLength);

	facetTPLabel = new GLLabel("Teleport to facet  #");
	facetPanel->Add(facetTPLabel);
	facetTeleport = new GLTextField(0, NULL);
	facetPanel->Add(facetTeleport);

	facetLinkLabel = new GLLabel("Link:");
	facetPanel->Add(facetLinkLabel);
	facetSILabel = new GLTextField(0, "");
	facetSILabel->SetEditable(true);

	facetPanel->Add(facetSILabel);

	facetStrLabel = new GLLabel("Structure:");
	facetPanel->Add(facetStrLabel);
	facetSuperDest = new GLTextField(0, NULL);
	facetPanel->Add(facetSuperDest);

	facetReLabel = new GLLabel("Profile:");
	facetPanel->Add(facetReLabel);
	facetRecType = new GLCombo(0);
	facetRecType->SetSize(4);
	facetRecType->SetValueAt(0, "None");
	facetRecType->SetValueAt(1, "Along \201");
	facetRecType->SetValueAt(2, "Along \202");
	facetRecType->SetValueAt(3, "Angular");
	facetPanel->Add(facetRecType);

	facetSpectrumToggle = new GLToggle(0, "Record Spectrum");
	facetPanel->Add(facetSpectrumToggle);

	facetTexBtn = new GLButton(0, "Mesh...");
	facetTexBtn->SetEnabled(false);
	facetPanel->Add(facetTexBtn);

	facetList = new GLList(0);
	facetList->SetWorker(&worker);
	facetList->SetGrid(true);
	facetList->SetSelectionMode(MULTIPLE_ROW);
	facetList->SetSize(5, 1);
	facetList->SetColumnWidths((int*)cWidth);
	facetList->SetColumnLabels((char **)cName);
	facetList->SetColumnLabelVisible(true);
	facetList->Sortable = true;
	Add(facetList);



	ClearFacetParams();
	LoadConfig();
	UpdateRecentMenu();
	UpdateRecentPARMenu();
	UpdateViewerPanel();
	PlaceComponents();
	CheckNeedsTexture();

	try {
		worker.SetProcNumber(nbProc);
	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "Failed to start working sub-process(es), simulation not available\n%s", e.GetMsg());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

	int index;
	try {
		for (index = 0; index < (int)materialPaths.size(); index++) //load materials
			worker.AddMaterial(&(materialPaths[index]));
	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "Failed to load material reflection file:\n%s\n%s", materialPaths[index].c_str(), e.GetMsg());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

	FileReader *f;
	try {
		worker.chi_distros.resize(3);

		f = new FileReader("param\\Distributions\\sum_psi_distr_0to4perlambdar_0.35_delta5E-3_full_pol.csv");
		//vertical (psi) distribution for different e_crit/e values
		//each row is for a logarithm of lambda_ratio, starting from -10 to +2
		//each column is for a psi angle, starting from 0 going to 1, with a delta of 0.005
		//where 1 corresponds to 4/lambda_ratio^0.35
		worker.psi_distro = worker.ImportCSV(f);
		SAFE_DELETE(f);
		f = new FileReader("param\\Distributions\\psi_chi_gamma10000_logsampled_-7to0_delta0.02_full_pol.csv");
		//each column corresponds to a Log10[PSI*(gamma/10000)] value. First column: -99, second column: -7, delta: 0.02, max: 0
		//each row corresponds to a    Log10[CHI*(gamma/10000)] value. First column: -7,                     delta: 0.02, max: 0
		worker.chi_distros[0] = worker.ImportCSV(f);
		SAFE_DELETE(f);

		//Parallel polarization
		f = new FileReader("param\\Distributions\\psi_chi_gamma10000_logsampled_-7to0_delta0.02_par_pol.csv");
		worker.chi_distros[1] = worker.ImportCSV(f);
		SAFE_DELETE(f);

		//Orthogonal polarization
		f = new FileReader("param\\Distributions\\psi_chi_gamma10000_logsampled_-7to0_delta0.02_ort_pol.csv");
		worker.chi_distros[2] = worker.ImportCSV(f);
		SAFE_DELETE(f);

		f = new FileReader("param\\Distributions\\degree_of_parallel_polarization_0to4perlambdar_0.35_delta5E-3.csv");
		worker.parallel_polarization = worker.ImportCSV(f);
		SAFE_DELETE(f);
	}
	catch (Error &e) {
		SAFE_DELETE(f);
		char errMsg[512];
		sprintf(errMsg, "Failed to load angular distribution file.\nIt should be in the param\\Distributions directory.\n%s", e.GetMsg());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

	if (checkForUpdates) {
		//Launch updater tool
		char command[2048];
		char CWD[MAX_PATH];
		_getcwd(CWD, MAX_PATH);
		//sprintf(tmp5,"%s\\synrad_updater_tmp.exe",CWD);
		if (FileUtils::Exist("synrad_updater_tmp.exe")) { //rename after new installation
			sprintf(command, "move \"%s\\synrad_updater_tmp.exe\" \"%s\\synrad_updater.exe\"", CWD, CWD);
			system(command);
		}
		//win_sparkle_set_appcast_url("http://test.xml");
		//win_sparkle_init();
		//win_sparkle_check_update_with_ui();
		if (FileUtils::Exist("synrad_updater.exe"))
			StartProc("synrad_updater.exe", STARTPROC_BACKGROUND);
		else GLMessageBox::Display("synrad_updater.exe not found. You will not receive updates to Synrad."
			"\n(You can disable checking for updates in Tools/Global Settings)", "Updater module missing.", GLDLG_OK, GLDLG_ICONINFO);
	}
	return GL_OK;
}

//-----------------------------------------------------------------------------
// Name: PlaceComponents()
// Desc: Place components on screen
//-----------------------------------------------------------------------------
void SynRad::PlaceComponents() {

	int sx = m_screenWidth - 205;
	int sy = 30;

	Place3DViewer();

	// ---------------------------------------------------------
	geomNumber->SetBounds(sx, 3, 202, 18);

	// Viewer settings ----------------------------------------
	togglePanel->SetBounds(sx, sy, 202, 112);

	togglePanel->SetCompBounds(showRule, 5, 20, 60, 18);
	togglePanel->SetCompBounds(showNormal, 70, 20, 60, 18);
	togglePanel->SetCompBounds(showUV, 135, 20, 60, 18);

	togglePanel->SetCompBounds(showLine, 5, 42, 60, 18);
	togglePanel->SetCompBounds(showLeak, 70, 42, 60, 18);
	togglePanel->SetCompBounds(showHit, 135, 42, 60, 18);

	togglePanel->SetCompBounds(showVolume, 5, 64, 60, 18);
	togglePanel->SetCompBounds(showTexture, 70, 64, 60, 18);
	togglePanel->SetCompBounds(showFilter, 135, 64, 60, 18);

	togglePanel->SetCompBounds(showVertex, 70, 86, 60, 18);
	togglePanel->SetCompBounds(showIndex, 137, 86, 60, 18);
	togglePanel->SetCompBounds(showMoreBtn, 5, 86, 55, 19);

	sy += (togglePanel->GetHeight() + 5);

	// Selected facet -----------------------------------------
	facetPanel->SetBounds(sx, sy, 202, 280);

	facetPanel->SetCompBounds(facetRLabel, 7, 15, 35, 18);
	facetPanel->SetCompBounds(facetReflType, 48, 15, 115, 18);

	facetPanel->SetCompBounds(facetSticking, 165, 15, 32, 18);

	facetPanel->SetCompBounds(facetDoScattering, 5, 36, 150, 18);

	facetPanel->SetCompBounds(facetRMSroughnessLabel, 7, 55, 35, 18);
	facetPanel->SetCompBounds(facetRMSroughness, 65, 55, 45, 18);

	facetPanel->SetCompBounds(facetAutoCorrLengthLabel, 113, 55, 35, 18);
	facetPanel->SetCompBounds(facetAutoCorrLength, 150, 55, 45, 18);

	PlaceScatteringControls(worker.newReflectionModel);

	facetPanel->SetCompBounds(facetSideLabel, 7, 80, 50, 18);
	facetPanel->SetCompBounds(facetSideType, 65, 80, 130, 18);

	facetPanel->SetCompBounds(facetTLabel, 7, 105, 100, 18);
	facetPanel->SetCompBounds(facetOpacity, 110, 105, 82, 18);

	facetPanel->SetCompBounds(facetAreaLabel, 7, 130, 100, 18);
	facetPanel->SetCompBounds(facetArea, 110, 130, 82, 18);

	facetPanel->SetCompBounds(facetTPLabel, 7, 155, 100, 18);
	facetPanel->SetCompBounds(facetTeleport, 110, 155, 82, 18);

	facetPanel->SetCompBounds(facetStrLabel, 7, 180, 55, 18); //Structure:
	facetPanel->SetCompBounds(facetSILabel, 65, 180, 42, 18); //Editable Textfield
	facetPanel->SetCompBounds(facetLinkLabel, 115, 180, 18, 18); //Link
	facetPanel->SetCompBounds(facetSuperDest, 148, 180, 42, 18); //Textfield

	facetPanel->SetCompBounds(facetReLabel, 7, 205, 60, 18);
	facetPanel->SetCompBounds(facetRecType, 65, 205, 130, 18);

	facetPanel->SetCompBounds(facetSpectrumToggle, 5, 230, 150, 18);

	facetPanel->SetCompBounds(facetDetailsBtn, 5, 255, 45, 18);
	facetPanel->SetCompBounds(facetCoordBtn, 53, 255, 44, 18);
	facetPanel->SetCompBounds(facetTexBtn, 101, 255, 50, 18);
	facetPanel->SetCompBounds(facetApplyBtn, 155, 255, 40, 18);


	sy += (facetPanel->GetHeight() + 5);

	// Simulation ---------------------------------------------
	simuPanel->SetBounds(sx, sy, 202, 219);
	int height = 20;
	simuPanel->SetCompBounds(startSimu, 5, height, 92, 19);
	simuPanel->SetCompBounds(resetSimu, 102, height, 93, 19);
	height += 25;
	simuPanel->SetCompBounds(autoFrameMoveToggle, 5, height, 65, 19);
	simuPanel->SetCompBounds(forceFrameMoveButton, 128, height, 66, 19);
	height += 25;
	simuPanel->SetCompBounds(modeLabel, 5, height, 30, 18);
	simuPanel->SetCompBounds(modeCombo, 40, height, 85, 18);
	height += 25;
	simuPanel->SetCompBounds(hitLabel, 5, height, 30, 18);
	simuPanel->SetCompBounds(hitNumber, 40, height, 155, 18);
	height += 25;
	simuPanel->SetCompBounds(desLabel, 5, height, 30, 18);
	simuPanel->SetCompBounds(desNumber, 40, height, 155, 18);
	height += 25;
	simuPanel->SetCompBounds(leakLabel, 5, height, 30, 18);
	simuPanel->SetCompBounds(leakNumber, 40, height, 155, 18);
	height += 25;
	simuPanel->SetCompBounds(doseLabel, 5, height, 30, 18);
	simuPanel->SetCompBounds(doseNumber, 40, height, 155, 18);
	height += 25;
	simuPanel->SetCompBounds(sTimeLabel, 5, height, 30, 18);
	simuPanel->SetCompBounds(sTime, 40, height, 155, 18);

	sy += (simuPanel->GetHeight() + 5);

	// ---------------------------------------------------------
	int lg = m_screenHeight - 23 /*- (nbFormula * 25)*/;

	facetList->SetBounds(sx, sy, 202, lg - sy);

	// ---------------------------------------------------------

	/*
	for (int i = 0; i < nbFormula; i++) {
	formulas[i].name->SetBounds(sx, lg + 5, 95, 18);
	formulas[i].value->SetBounds(sx + 90, lg + 5, 87, 18);
	formulas[i].setBtn->SetBounds(sx + 182, lg + 5, 20, 18);
	lg += 25;
	}
	*/
}

//-----------------------------------------------------------------------------
// Name: ClearFacetParams()
// Desc: Reset selected facet parameters.
//-----------------------------------------------------------------------------

void SynRad::ClearFacetParams()
{
	facetPanel->SetTitle("Selected Facet (none)");
	facetReflType->SetSelectedValue("");
	facetReflType->SetEditable(false);
	facetSticking->Clear();
	facetSticking->SetVisible(false);
	facetDoScattering->SetState(0);
	facetDoScattering->SetEnabled(false);
	facetRMSroughness->Clear();
	facetRMSroughness->SetEditable(false);
	facetAutoCorrLength->Clear();
	facetAutoCorrLength->SetEditable(false);
	facetTeleport->Clear();
	facetTeleport->SetEditable(false);
	facetArea->SetEditable(false);
	facetArea->Clear();
	facetSuperDest->Clear();
	facetSuperDest->SetEditable(false);
	facetSILabel->Clear();
	facetSILabel->SetEditable(false);
	facetOpacity->Clear();
	facetOpacity->SetEditable(false);
	facetSideType->SetSelectedValue("");
	facetSideType->SetEditable(false);

	facetRecType->SetSelectedValue("");
	facetRecType->SetEditable(false);
	facetSpectrumToggle->SetState(0);
	facetSpectrumToggle->SetEnabled(false);
}

//-----------------------------------------------------------------------------
// Name: ApplyFacetParams()
// Desc: Apply facet parameters.
//-----------------------------------------------------------------------------

void SynRad::ApplyFacetParams() {
	if (!AskToReset()) return;
	changedSinceSave = true;
	Geometry *geom = worker.GetGeometry();
	size_t nbFacet = geom->GetNbFacet();

	//Reflection type and sticking
	double sticking;
	bool doSticking = false;
	int reflType = facetReflType->GetSelectedIndex();

	if (reflType == 0 || reflType == 1) { //Diffuse or mirror
		if (facetSticking->GetNumber(&sticking)) {
			if (sticking<0.0 || sticking>1.0) {
				GLMessageBox::Display("Sticking must be in the range [0,1]", "Error", GLDLG_OK, GLDLG_ICONERROR);
				UpdateFacetParams();
				return;
			}
			doSticking = true;
		}
		else { //Not a double number
			if (strcmp(facetSticking->GetText(), "...") == 0) doSticking = false;
			else {
				GLMessageBox::Display("Invalid sticking number", "Error", GLDLG_OK, GLDLG_ICONERROR);
				UpdateFacetParams();
				return;
			}
		}
	}

	//Scattering, roughness, autocorr.length
	bool doScattering = false;
	bool doRoughness = false;
	bool doCorrLength = false;
	double roughness, corrLength;

	if (facetDoScattering->GetState() < 2) {
		doScattering = true;
		if (facetDoScattering->GetState() == 1) { //Do scattering, read roughness and autocorr.length

			// Roughness
			if (facetRMSroughness->GetNumber(&roughness)) {
				if (roughness < 0.0) {
					GLMessageBox::Display("Roughness must be non-negative", "Error", GLDLG_OK, GLDLG_ICONERROR);
					UpdateFacetParams();
					return;
				}
				doRoughness = true;
			}
			else {
				if (strcmp(facetRMSroughness->GetText(), "...") == 0) doRoughness = false;
				else {
					GLMessageBox::Display("Invalid roughness number", "Error", GLDLG_OK, GLDLG_ICONERROR);
					UpdateFacetParams();
					return;
				}
			}

			// Autocorrelation length
			if (facetAutoCorrLength->GetNumber(&corrLength)) {
				if (corrLength <= 0.0) {
					GLMessageBox::Display("Autocorr.length must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
					UpdateFacetParams();
					return;
				}
				doCorrLength = true;
			}
			else {
				if (strcmp(facetAutoCorrLength->GetText(), "...") == 0) doCorrLength = false;
				else {
					GLMessageBox::Display("Invalid autocorrelation length", "Error", GLDLG_OK, GLDLG_ICONERROR);
					UpdateFacetParams();
					return;
				}
			}
		}
	}


	// teleport
	int teleport;
	bool doTeleport = false;

	if (facetTeleport->GetNumberInt(&teleport)) {
		if (teleport<-1 || teleport>nbFacet) {
			GLMessageBox::Display("Invalid teleport destination\n(If no teleport: set number to 0)", "Error", GLDLG_OK, GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
		else if (teleport > 0 && geom->GetFacet(teleport - 1)->selected) {
			char tmp[256];
			sprintf(tmp, "The teleport destination of facet #%d can't be itself!", teleport);
			GLMessageBox::Display(tmp, "Error", GLDLG_OK, GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
		doTeleport = true;
	}
	else {
		if (strcmp(facetTeleport->GetText(), "...") == 0) doTeleport = false;
		else {
			GLMessageBox::Display("Invalid teleport destination\n(If no teleport: set number to 0)", "Error", GLDLG_OK, GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
	}


	// opacity
	double opacity;
	bool doOpacity = false;
	if (facetOpacity->GetNumber(&opacity)) {
		if (opacity<0.0 || opacity>1.0) {
			GLMessageBox::Display("Opacity must be in the range [0,1]", "Error", GLDLG_OK, GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
		doOpacity = true;
	}
	else {
		if (strcmp(facetOpacity->GetText(), "...") == 0) doOpacity = false;
		else {
			GLMessageBox::Display("Invalid opacity number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
	}


	// Superstructure
	int superStruct;
	bool doSuperStruct = false;
	if (sscanf(facetSILabel->GetText(), "%d", &superStruct) > 0 && superStruct > 0 && superStruct <= geom->GetNbStructure()) doSuperStruct = true;
	else {
		if (strcmp(facetSILabel->GetText(), "...") == 0) doSuperStruct = false;
		else {
			GLMessageBox::Display("Invalid superstructre number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
	}


	// Super structure destination (link)
	int superDest;
	bool doLink = false;
	if (strcmp(facetSuperDest->GetText(), "none") == 0 || strcmp(facetSuperDest->GetText(), "no") == 0 || strcmp(facetSuperDest->GetText(), "0") == 0) {
		doLink = true;
		superDest = 0;
	}
	else if (sscanf(facetSuperDest->GetText(), "%d", &superDest) > 0) {
		if (superDest == superStruct) {
			GLMessageBox::Display("Link and superstructure can't be the same", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		else if (superDest < 0 || superDest > geom->GetNbStructure()) {
			GLMessageBox::Display("Link destination points to a structure that doesn't exist", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		else
			doLink = true;
	}
	else if (strcmp(facetSuperDest->GetText(), "...") == 0) doLink = false;
	else {
		GLMessageBox::Display("Invalid superstructure destination", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	// Record type
	int rType = facetRecType->GetSelectedIndex();

	// Spectrum recording
	int specType = facetSpectrumToggle->GetState();

	// 2sided
	int is2Sided = facetSideType->GetSelectedIndex();

	bool structChanged = false; //if a facet gets into a new structure, we have to re-render the geometry
	// Update facets (local)
	for (int i = 0; i < nbFacet; i++) {
		Facet *f = geom->GetFacet(i);
		if (f->selected) {
			if (reflType >= 0) {
				if (reflType >= 2)
					f->sh.reflectType = reflType + 8; //Material reflections: 10, 11, 12...
				else //Diffuse or Mirror
					f->sh.reflectType = reflType;
			}
			if (doSticking) f->sh.sticking = sticking;
			if (doScattering) f->sh.doScattering = facetDoScattering->GetState();
			if (worker.newReflectionModel) { //Apply values directly
				if (doRoughness) f->sh.rmsRoughness = roughness*1E-9; //nm->m
				if (doCorrLength) f->sh.autoCorrLength = corrLength*1E-9; //nm->m
			}
			else { //Convert roughness ratio
				if (doRoughness) {
					f->sh.autoCorrLength = 1E-5; //10000 nm
					f->sh.rmsRoughness = roughness*f->sh.autoCorrLength; //roughness = roughness_ratio * autocorr_length
				}
			}
			if (doTeleport) f->sh.teleportDest = teleport;
			if (doOpacity) f->sh.opacity = opacity;

			if (rType >= 0) {
				f->sh.profileType = rType;
				f->sh.isProfile = (rType != REC_NONE);

			} if (profilePlotter) profilePlotter->Refresh();
			if (specType < 2) { //Not mixed state
				f->sh.hasSpectrum = specType;

			} if (spectrumPlotter) spectrumPlotter->Refresh();
			if (is2Sided >= 0) f->sh.is2sided = is2Sided;
			if (doSuperStruct) {
				if (f->sh.superIdx != (superStruct - 1)) {
					f->sh.superIdx = superStruct - 1;
					structChanged = true;
				}
			}
			if (doLink) {
				f->sh.superDest = superDest;
				if (superDest) f->sh.opacity = 1.0; // Force opacity for link facet
			}
			f->UpdateFlags();
		}
	}
	if (structChanged) geom->BuildGLList();
	// Send to sub process
	try { worker.Reload(); }
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	UpdateFacetParams();
}

//-----------------------------------------------------------------------------
// Name: UpdateFacetParams()
// Desc: Update selected facet parameters.
//-----------------------------------------------------------------------------

void SynRad::UpdateFacetParams(bool updateSelection) {

	char tmp[256];
	int sel0 = -1;

	// Update params
	Geometry *geom = worker.GetGeometry();
	size_t nbSel = geom->GetNbSelectedFacets();
	if (nbSel > 0) {

		Facet *f0;
		Facet *f;

		// Get list of selected facet
		std::vector<size_t> selectedFacets = geom->GetSelectedFacets();

		f0 = geom->GetFacet(selectedFacets[0]);

		double f0Area = f0->GetArea();
		double sumArea = f0Area; //sum facet area

		bool reflectTypeE = true;
		bool stickingE = true;
		bool rmsRoughnessE = true;
		bool autoCorrLengthE = true;
		bool teleportE = true;
		bool opacityE = true;
		bool superDestE = true;
		bool superIdxE = true;
		bool recordE = true;
		bool is2sidedE = true;
		bool hasSpectrumE = true;
		bool doScatteringE = true;
		bool isAllRegular = (f0->sh.reflectType < 2); //All facets are (mirror OR diffuse)
		bool isAllDiffuse = (f0->sh.reflectType == 0); //All facets are diffuse (otherwise enable choice for rough surface scattering)

		for (size_t sel = 1; sel < selectedFacets.size(); sel++) {
			f = geom->GetFacet(selectedFacets[sel]);
			double fArea = f->GetArea();
			reflectTypeE = reflectTypeE && (f0->sh.reflectType == f->sh.reflectType);
			stickingE = stickingE && (abs(f0->sh.sticking - f->sh.sticking) < 1e-7);
			if (worker.newReflectionModel) { //RMS roughness compare
				rmsRoughnessE = rmsRoughnessE && (abs(f0->sh.rmsRoughness - f->sh.rmsRoughness) < 1e-15);
			}
			else { //Roughness ratio compare
				rmsRoughnessE = rmsRoughnessE && (abs(f0->sh.rmsRoughness / f0->sh.autoCorrLength - f->sh.rmsRoughness / f->sh.autoCorrLength) < 1e-10);
			}
			autoCorrLengthE = autoCorrLengthE && (abs(f0->sh.autoCorrLength - f->sh.autoCorrLength) < 1e-15);
			doScatteringE = doScatteringE && (f0->sh.doScattering == f->sh.doScattering);
			isAllRegular = isAllRegular && (f->sh.reflectType < 2);
			isAllDiffuse = isAllDiffuse && (f->sh.reflectType == 0);
			teleportE = teleportE && (f0->sh.teleportDest == f->sh.teleportDest);
			opacityE = opacityE && (abs(f0->sh.opacity - f->sh.opacity) < 1e-7);
			superDestE = superDestE && (f0->sh.superDest == f->sh.superDest);
			superIdxE = superIdxE && (f0->sh.superIdx == f->sh.superIdx);
			is2sidedE = is2sidedE && (f0->sh.is2sided == f->sh.is2sided);

			recordE = recordE && (f0->sh.profileType == f->sh.profileType);
			hasSpectrumE = hasSpectrumE && (f0->sh.hasSpectrum == f->sh.hasSpectrum);
			sumArea += fArea;
		}

		if (nbSel == 1)
			sprintf(tmp, "Selected Facet (#%zd)", selectedFacets[0] + 1);
		else
			sprintf(tmp, "Selected Facet (%zd selected)", selectedFacets.size());

		// Old STR compatibility
		if (stickingE && f0->sh.superDest) stickingE = false;

		facetPanel->SetTitle(tmp);
		if (selectedFacets.size() > 1) facetAreaLabel->SetText("Sum Area (cm\262):");
		else facetAreaLabel->SetText("Area (cm\262):");
		facetArea->SetText(sumArea);

		if (reflectTypeE) {
			if (isAllRegular) facetReflType->SetSelectedIndex(f0->sh.reflectType); //Diffuse or Mirror
			else if (f0->sh.reflectType >= 10 && f0->sh.reflectType < 10 + worker.materials.size())
				facetReflType->SetSelectedIndex(f0->sh.reflectType - 8); //Map 10,11,12 to 2,3,4...
			else
				facetReflType->SetSelectedValue("Invalid material");
		}
		else facetReflType->SetSelectedValue("...");

		if (isAllRegular) {
			if (stickingE) facetSticking->SetText(f0->sh.sticking);
			else facetSticking->SetText("...");
		}
		else {
			facetSticking->SetText("");
		}
		if (doScatteringE)
			facetDoScattering->SetState(f0->sh.doScattering);
		else
			facetDoScattering->SetState(2); //Mixed state
		facetDoScattering->AllowMixedState(!doScatteringE);

		if (isAllDiffuse) {
			facetDoScattering->SetEnabled(false);
			facetRMSroughness->SetText("");
			facetAutoCorrLength->SetText("");
		}
		else {
			facetDoScattering->SetEnabled(true);
			if (rmsRoughnessE) {
				if (worker.newReflectionModel) { //RMS roughness (nm)
					facetRMSroughness->SetText(f0->sh.rmsRoughness*1E9);
				}
				else { //Roughness ratio
					facetRMSroughness->SetText(f0->sh.rmsRoughness / f0->sh.autoCorrLength);
				}
			}
			else facetRMSroughness->SetText("..."); //m->nm
			if (autoCorrLengthE) facetAutoCorrLength->SetText(f0->sh.autoCorrLength*1E9); else facetAutoCorrLength->SetText("..."); //m->nm
		}

		if (teleportE) facetTeleport->SetText(f0->sh.teleportDest); else facetTeleport->SetText("...");
		if (opacityE)facetOpacity->SetText(f0->sh.opacity); else facetOpacity->SetText("...");
		if (is2sidedE) facetSideType->SetSelectedIndex(f0->sh.is2sided); else facetSideType->SetSelectedValue("...");

		if (recordE) facetRecType->SetSelectedIndex(f0->sh.profileType); else facetRecType->SetSelectedValue("...");
		if (hasSpectrumE) facetSpectrumToggle->SetState(f0->sh.hasSpectrum); else facetSpectrumToggle->SetState(2);
		facetSpectrumToggle->AllowMixedState(!hasSpectrumE);
		if (superDestE) {
			if (f0->sh.superDest == 0) {
				facetSuperDest->SetText("no");
			}
			else {
				sprintf(tmp, "%zd", f0->sh.superDest);
				facetSuperDest->SetText(tmp);
			}
		}
		else {
			facetSuperDest->SetText("...");
		}
		if (superIdxE) {
			sprintf(tmp, "%zd", f0->sh.superIdx + 1);
			facetSILabel->SetText(tmp);
		}
		else {
			facetSILabel->SetText("...");
		}
		if (updateSelection) {
			if (nbSel > 1000 || geom->GetNbFacet() > 50000) { //If it would take too much time to look up every selected facet in the list
				facetList->ReOrder();
				facetList->SetSelectedRows(selectedFacets, false);
			}
			else {
				facetList->SetSelectedRows(selectedFacets, true);
			}
			facetList->lastRowSel = -1;
		}

		facetReflType->SetEditable(true);
		facetSticking->SetVisible(isAllRegular);
		facetDoScattering->SetEnabled(!isAllDiffuse);
		facetRMSroughness->SetEditable(!doScatteringE || f0->sh.doScattering);
		facetAutoCorrLength->SetEditable(!doScatteringE || f0->sh.doScattering);

		facetTeleport->SetEditable(true);
		facetOpacity->SetEditable(true);
		facetSuperDest->SetEditable(true);
		facetSILabel->SetEditable(true);
		facetSideType->SetEditable(true);

		facetRecType->SetEditable(true);
		facetSpectrumToggle->SetEnabled(true);

		facetApplyBtn->SetEnabled(false);
		menu->GetSubMenu("Facet")->SetEnabled(MENU_FACET_MESH, true);
		facetTexBtn->SetEnabled(true);

	}
	else {
		ClearFacetParams();
		menu->GetSubMenu("Facet")->SetEnabled(MENU_FACET_MESH, false);
		facetTexBtn->SetEnabled(false);
		if (updateSelection) facetList->ClearSelection();
	}

	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (texturePlotter) texturePlotter->Update(m_fTime, true); //Selected facet change
	//if( outgassingMap ) outgassingMap->Update(m_fTime,true);
}


bool SynRad::EvaluateVariable(VLIST *v) {
	bool ok = true;
	Geometry* geom = worker.GetGeometry();
	size_t nbFacet = geom->GetNbFacet();
	size_t idx;

	if ((idx = GetVariable(v->name, "A")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->counterCache.nbAbsorbed;
	}
	else if ((idx = GetVariable(v->name, "H")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->counterCache.nbHit;
	}
	else if ((idx = GetVariable(v->name, "F")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->counterCache.fluxAbs / worker.no_scans;
	}
	else if ((idx = GetVariable(v->name, "P")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->counterCache.powerAbs / worker.no_scans;
	}
	else if ((idx = GetVariable(v->name, "AR")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) v->value = geom->GetFacet(idx - 1)->sh.area;
	}
	else if (_stricmp(v->name, "SCANS") == 0) {
		v->value = worker.no_scans;
	}
	else if (_stricmp(v->name, "SUMDES") == 0) {
		v->value = (double)worker.nbDesorption;
	}
	else if (_stricmp(v->name, "SUMABS") == 0) {
		v->value = (double)worker.nbAbsorption;
	}
	else if (_stricmp(v->name, "SUMHIT") == 0) {
		v->value = (double)worker.nbHit;
	}
	else if (_stricmp(v->name, "SUMFLUX") == 0) {
		v->value = worker.totalFlux / worker.no_scans;
	}
	else if (_stricmp(v->name, "SUMPOWER") == 0) {
		v->value = worker.totalPower / worker.no_scans;
	}
	else if (_stricmp(v->name, "MPP") == 0) {
		v->value = worker.distTraveledTotal / (double)worker.nbDesorption;
	}
	else if (_stricmp(v->name, "MFP") == 0) {
		v->value = worker.distTraveledTotal / (double)worker.nbHit;
	}
	else if (_stricmp(v->name, "ABSAR") == 0) {
		double sumArea = 0.0;
		for (int i = 0; i < geom->GetNbFacet(); i++) {
			Facet *f = geom->GetFacet(i);
			if (f->sh.sticking > 0.0) sumArea += f->sh.area*f->sh.opacity*(f->sh.is2sided ? 2.0 : 1.0);
		}
		v->value = sumArea;
	}
	else if (_stricmp(v->name, "KB") == 0) {
		v->value = 1.3806504e-23;
	}
	else if (_stricmp(v->name, "R") == 0) {
		v->value = 8.314472;
	}
	else if (_stricmp(v->name, "Na") == 0) {
		v->value = 6.02214179e23;
	}
	else if ((beginsWith(v->name, "SUM(") || beginsWith(v->name, "sum(")) /*|| (beginsWith(v->name, "AVG("))*/ && endsWith(v->name, ")")) {
		//bool avgMode = beginsWith(v->name, "AVG("); //else SUM mode
		std::string inside = v->name; inside.erase(0, 4); inside.erase(inside.size() - 1, 1);
		std::vector<std::string> tokens = SplitString(inside, ',');
		if (!Contains({ 2,3 }, tokens.size()))
			return false;
		/*if (avgMode) {
			if (!Contains({ "P","DEN","Z" }, tokens[0]))
				return false;
		}*/
		else {
			if (!Contains({ "H","A","F","P","AR","h","a","f","p","ar" }, tokens[0]))
				return false;
		}
		std::vector<size_t> facetsToSum;
		if (tokens.size() == 3) { // Like SUM(H,3,6) = H3 + H4 + H5 + H6
			size_t startId, endId, pos;
			try {
				startId = std::stol(tokens[1], &pos); if (pos != tokens[1].size() || startId > geom->GetNbFacet() || startId == 0) return false;
				endId = std::stol(tokens[2], &pos); if (pos != tokens[2].size() || endId > geom->GetNbFacet() || endId == 0) return false;
			}
			catch (...) {
				return false;
			}
			if (!(startId < endId)) return false;
			facetsToSum = std::vector<size_t>(endId - startId + 1);
			std::iota(facetsToSum.begin(), facetsToSum.end(), startId - 1);
		}
		else { //Selection group
			if (!(beginsWith(tokens[1], "S") || beginsWith(tokens[1], "s"))) return false;
			std::string selIdString = tokens[1]; selIdString.erase(0, 1);
			if (Contains({ "EL","el" }, selIdString)) { //Current selections
				facetsToSum = geom->GetSelectedFacets();
			}
			else {
				size_t selGroupId, pos;
				try {
					selGroupId = std::stol(selIdString, &pos); if (pos != selIdString.size() || selGroupId > selections.size() || selGroupId == 0) return false;
				}
				catch (...) {
					return false;
				}
				facetsToSum = selections[selGroupId - 1].selection;
			}
		}
		llong sumLL = 0;
		double sumD = 0.0;
		double sumArea = 0.0; //We average by area
		for (auto sel : facetsToSum) {
			if (Contains({ "H","h" }, tokens[0])) {
				sumLL += geom->GetFacet(sel)->counterCache.nbHit;
			}
			else if (Contains({ "A","a" }, tokens[0])) {
				sumLL += geom->GetFacet(sel)->counterCache.nbAbsorbed;
			}
			else if (Contains({ "AR","ar" }, tokens[0])) {
				sumArea += geom->GetFacet(sel)->GetArea();
			}
			else if (Contains({ "F","f" }, tokens[0])) {
				sumD += geom->GetFacet(sel)->counterCache.fluxAbs /*/ worker.no_scans*/;
			}
			else if (Contains({ "P","p" }, tokens[0])) {
				sumD += geom->GetFacet(sel)->counterCache.powerAbs/* / worker.no_scans*/;
			}
			else return false;
		}
		//if (avgMode) v->value = sumD * worker.GetMoleculesPerTP(worker.displayedMoment)*1E4 / sumArea;
		/*else*/ if (Contains({ "AR" , "ar"},tokens[0])) v->value = sumArea;
		else if (Contains({ "F","P","f","p" }, tokens[0])) v->value = sumD / worker.no_scans;
		else v->value = (double)sumLL;
	}
	else ok = false;
	return ok;
}

void SynRad::UpdatePlotters()
{
	if (mApp->profilePlotter) mApp->profilePlotter->Update(m_fTime, true);
	if (mApp->spectrumPlotter) mApp->spectrumPlotter->Update(m_fTime, true);
	if (mApp->texturePlotter) mApp->texturePlotter->Update(m_fTime, true);
}


//-----------------------------------------------------------------------------
// Name: FrameMove()
// Desc: Called once per frame, the call is the entry point for animating
//       the scene.
//-----------------------------------------------------------------------------
int SynRad::FrameMove()
{
	if (worker.running && ((m_fTime - lastUpdate) >= 1.0f)) {
		if (textureSettings) textureSettings->Update();
	}
	Interface::FrameMove(); //might reset lastupdate
	char tmp[256];
	if (globalSettings) globalSettings->SMPUpdate();


	if ((m_fTime - worker.startTime <= 2.0f) && worker.running) {
		hitNumber->SetText("Starting...");
		desNumber->SetText("Starting...");
		doseNumber->SetText("Starting...");
	}
	else {
		sprintf(tmp, "%s (%s)", FormatInt(worker.nbHit, "hit"), FormatPS(hps, "hit"));
		hitNumber->SetText(tmp);
		sprintf(tmp, "%s (%s)", FormatInt(worker.nbDesorption, "des"), FormatPS(dps, "des"));
		desNumber->SetText(tmp);
	}

	if (worker.no_scans) {
		sprintf(tmp, "Scn:%.1f F=%.3g P=%.3g", (worker.no_scans == 1.0) ? 0.0 : worker.no_scans,
			worker.totalFlux / worker.no_scans, worker.totalPower / worker.no_scans);
		if (worker.nbTrajPoints > 0) doseNumber->SetText(tmp);
	}
	else {
		doseNumber->SetText("");
	}
	return GL_OK;
}

// ----------------------------------------------------------------

void SynRad::UpdateFacetHits(bool all) {
	char tmp[256];
	Geometry *geom = worker.GetGeometry();

	try {
		// Facet list
		if (geom->IsLoaded()) {


			int sR, eR;
			if (all)
			{
				sR = 0;
				eR = (int)facetList->GetNbRow() - 1;
			}
			else
			{
				facetList->GetVisibleRows(&sR, &eR);
			}

			for (int i = sR; i <= eR; i++) {
				int facetId = facetList->GetValueInt(i, 0) - 1;
				if (facetId == -2) facetId = i;
				if (i >= geom->GetNbFacet()) {
					char errMsg[512];
					sprintf(errMsg, "Synrad::UpdateFacetHits()\nError while updating facet hits. Was looking for facet #%d in list.\nSynrad will now autosave and crash.", i + 1);
					GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
					AutoSave(true);
				}
				Facet *f = geom->GetFacet(facetId);
				sprintf(tmp, "%d", facetId + 1);
				facetList->SetValueAt(0, i, tmp);
				facetList->SetColumnLabel(1, "Hits");
				sprintf(tmp, "%I64d", f->counterCache.nbHit);
				facetList->SetValueAt(1, i, tmp);
				sprintf(tmp, "%.3g", f->counterCache.fluxAbs / worker.no_scans);
				facetList->SetValueAt(2, i, tmp);
				sprintf(tmp, "%.3g", f->counterCache.powerAbs / worker.no_scans);
				facetList->SetValueAt(3, i, tmp);
				sprintf(tmp, "%I64d", f->counterCache.nbAbsorbed);
				facetList->SetValueAt(4, i, tmp);

			}

		}
	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "%s\nError while updating facet hits", e.GetMsg());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Name: RestoreDeviceObjects()
// Desc: Initialize scene objects.
//-----------------------------------------------------------------------------
int SynRad::RestoreDeviceObjects()
{
	RestoreDeviceObjects_shared();

	//Different SynRad implementations
	RVALIDATE_DLG(facetMesh);
	RVALIDATE_DLG(facetDetails);
	RVALIDATE_DLG(smartSelection);
	RVALIDATE_DLG(viewer3DSettings);
	RVALIDATE_DLG(textureSettings);
	RVALIDATE_DLG(globalSettings);
	RVALIDATE_DLG(profilePlotter);
	RVALIDATE_DLG(texturePlotter);

	//Synrad only
	RVALIDATE_DLG(regionInfo);
	RVALIDATE_DLG(exportDesorption);
	RVALIDATE_DLG(spectrumPlotter);
	RVALIDATE_DLG(trajectoryDetails);

	return GL_OK;
}

//-----------------------------------------------------------------------------
// Name: InvalidateDeviceObjects()
// Desc: Free all alocated resource
//-----------------------------------------------------------------------------

int SynRad::InvalidateDeviceObjects()
{
	InvalidateDeviceObjects_shared();

	//Different SynRad implementations
	IVALIDATE_DLG(facetMesh);
	IVALIDATE_DLG(facetDetails);
	IVALIDATE_DLG(smartSelection);
	IVALIDATE_DLG(viewer3DSettings);
	IVALIDATE_DLG(textureSettings);
	IVALIDATE_DLG(globalSettings);
	IVALIDATE_DLG(profilePlotter);
	IVALIDATE_DLG(texturePlotter);

	//Synrad only
	IVALIDATE_DLG(regionInfo);
	IVALIDATE_DLG(exportDesorption);
	IVALIDATE_DLG(spectrumPlotter);
	IVALIDATE_DLG(trajectoryDetails);

	return GL_OK;
}

/*
void SynRad::SaveFileAs() {

	FILENAME *fn = GLFileBox::SaveFile(currentDir, worker.GetShortFileName(), "Save File", fileSFilters, 0);

	GLProgress *progressDlg2 = new GLProgress("Saving file...", "Please wait");
	progressDlg2->SetProgress(0.0);
	progressDlg2->SetVisible(true);
	//GLWindowManager::Repaint();
	if (fn) {
		try {
			worker.SaveGeometry(fn->fullName, progressDlg2);
			ResetAutoSaveTimer();
			changedSinceSave = false;
			UpdateCurrentDir(worker.fullFileName);
			UpdateTitle();
			AddRecent(worker.fullFileName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
			RemoveRecent(fn->fullName);
		}
	}

	progressDlg2->SetVisible(false);
	SAFE_DELETE(progressDlg2);
}
*/

/*
void SynRad::ExportTextures(int grouping, int mode) {

	Geometry *geom = worker.GetGeometry();
	if (geom->GetNbSelectedFacets() == 0) {
		GLMessageBox::Display("Empty selection", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	if (!worker.IsDpInitialized()) {
		GLMessageBox::Display("Worker Dataport not initialized yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	FILENAME *fn = GLFileBox::SaveFile(currentDir, NULL, "Save Texture File", fileTexFilters, nbTexFilter);

	if (fn) {

		try {
			worker.ExportTextures(fn->fullName, grouping, mode, true, true);
			//UpdateCurrentDir(fn->fullName);
			//UpdateTitle();
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}

	}

}
*/

void SynRad::SaveFile() {
	if (strlen(worker.fullFileName) > 0) {

		GLProgress *progressDlg2 = new GLProgress("Saving...", "Please wait");
		progressDlg2->SetProgress(0.5);
		progressDlg2->SetVisible(true);
		//GLWindowManager::Repaint();

		try {
			worker.SaveGeometry(worker.fullFileName, progressDlg2, false);
			ResetAutoSaveTimer();
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), worker.GetFileName());
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
		progressDlg2->SetVisible(false);
		SAFE_DELETE(progressDlg2);
		changedSinceSave = false;

	}
	else SaveFileAs();
}

void SynRad::LoadFile(char *fName) {


	char fullName[512];
	char shortName[512];
	strcpy(fullName, "");

	if (fName == NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Open File", fileLFilters, 0);
		if (fn)
			strcpy(fullName, fn->fullName);
	}
	else {
		strcpy(fullName, fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Preparing to load file...", "Please wait");
	progressDlg2->SetVisible(true);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if (strlen(fullName) == 0) {
		progressDlg2->SetVisible(false);
		SAFE_DELETE(progressDlg2);
		return;
	}


	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);

	try {
		ClearFormulas();
		ClearAllSelections();
		ClearAllViews();
		ClearRegions();
		worker.LoadGeometry(fullName);

		Geometry *geom = worker.GetGeometry();

		// Default initialisation
		viewer[0]->SetWorker(&worker);
		viewer[1]->SetWorker(&worker);
		viewer[2]->SetWorker(&worker);
		viewer[3]->SetWorker(&worker);
		startSimu->SetEnabled(true);
		ClearFacetParams();
		nbDesStart = worker.nbDesorption;
		nbHitStart = worker.nbHit;
		AddRecent(fullName);
		geom->viewStruct = -1;

		RebuildPARMenus();
		UpdateStructMenu();
		if (profilePlotter) profilePlotter->Reset();
		if (spectrumPlotter) spectrumPlotter->Reset();
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
		UpdatePlotters();
		if (textureSettings) textureSettings->Update();
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();

	}
	catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(fName);

	}
	progressDlg2->SetVisible(false);
	SAFE_DELETE(progressDlg2);
	changedSinceSave = false;
}

//-----------------------------------------------------------------------------

void SynRad::InsertGeometry(bool newStr, char *fName) {
	if (!AskToReset()) return;
	ResetSimulation(false);

	char fullName[512];
	char shortName[512];
	strcpy(fullName, "");

	if (fName == NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Open File", fileLFilters, 0);
		if (fn)
			strcpy(fullName, fn->fullName);
	}
	else {
		strcpy(fullName, fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Loading file...", "Please wait");
	progressDlg2->SetVisible(true);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if (strlen(fullName) == 0) {
		progressDlg2->SetVisible(false);
		SAFE_DELETE(progressDlg2);
		return;
	}


	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);

	try {

		worker.LoadGeometry(fullName, true, newStr);
		Geometry *geom = worker.GetGeometry();

		startSimu->SetEnabled(true);

		AddRecent(fullName);
		geom->viewStruct = -1;

		RebuildPARMenus();
		UpdateStructMenu();
		if (profilePlotter) profilePlotter->Reset();
		if (spectrumPlotter) spectrumPlotter->Reset();

		geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();

		UpdatePlotters();
		//if(outgassingMap) outgassingMap->Update(m_fTime,true);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();
		if (regionInfo) regionInfo->Update();

	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(fName);
	}
	progressDlg2->SetVisible(false);
	SAFE_DELETE(progressDlg2);
	changedSinceSave = true;
}

void SynRad::StartStopSimulation() {
	if ((int)worker.regions.size() == 0) {
		GLMessageBox::Display("No regions loaded", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	worker.StartStop(m_fTime, 0);
	UpdatePlotters();
	if (autoUpdateFormulas && formulaEditor && formulaEditor->IsVisible()) formulaEditor->ReEvaluate();

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
// Name: EventProc()
// Desc: Message proc function to handle key and mouse input
//-----------------------------------------------------------------------------
void SynRad::ProcessMessage(GLComponent *src, int message)
{
	if (ProcessMessage_shared(src, message)) return; //Already processed by common interface

	Geometry *geom = worker.GetGeometry();

	switch (message) {

		//MENU --------------------------------------------------------------------
	case MSG_MENU:
		switch (src->GetId()) {

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
			if (GLMessageBox::Display("Remove all magnetic regions?", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
				if (worker.running) worker.Stop_Public();
				ClearRegions();
			}
			break;

		case MENU_FILE_EXPORTTEXTURE_AREA:
			ExportTextures(0, 0);
			break;
		case MENU_FILE_EXPORTTEXTURE_MCHITS:
			ExportTextures(0, 1);
			break;
		case MENU_FILE_EXPORTTEXTURE_FLUX:
			ExportTextures(0, 2);
			break;
		case MENU_FILE_EXPORTTEXTURE_POWER:
			ExportTextures(0, 3);
			break;
		case MENU_FILE_EXPORTTEXTURE_FLUXPERAREA:
			ExportTextures(0, 4);
			break;
		case MENU_FILE_EXPORTTEXTURE_POWERPERAREA:
			ExportTextures(0, 5);
			break;

		case MENU_FILE_EXPORTTEXTURE_AREA_COORD:
			ExportTextures(1, 0);
			break;
		case MENU_FILE_EXPORTTEXTURE_MCHITS_COORD:
			ExportTextures(1, 1);
			break;
		case MENU_FILE_EXPORTTEXTURE_FLUX_COORD:
			ExportTextures(1, 2);
			break;
		case MENU_FILE_EXPORTTEXTURE_POWER_COORD:
			ExportTextures(1, 3);
			break;
		case MENU_FILE_EXPORTTEXTURE_FLUXPERAREA_COORD:
			ExportTextures(1, 4);
			break;
		case MENU_FILE_EXPORTTEXTURE_POWERPERAREA_COORD:
			ExportTextures(1, 5);
			break;

		case MENU_FILE_EXPORT_DESORP:
			if (!geom->IsLoaded()) {
				GLMessageBox::Display("No geometry loaded.", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (!exportDesorption) exportDesorption = new ExportDesorption(geom, &worker);
			exportDesorption->SetVisible(true);
			break;

		case MENU_EDIT_TSCALING:
			if (!textureSettings || !textureSettings->IsVisible()) {
				SAFE_DELETE(textureSettings);
				textureSettings = new TextureSettings();
				textureSettings->Display(&worker, viewer);
			}
			break;
		case MENU_EDIT_GLOBALSETTINGS:
			if (!globalSettings) globalSettings = new GlobalSettings();
			globalSettings->Display(&worker);
			break;
		case MENU_TOOLS_PROFPLOTTER:
			if (!profilePlotter) profilePlotter = new ProfilePlotter();
			profilePlotter->Display(&worker);
			break;
		case MENU_TOOLS_TEXPLOTTER:
			if (!texturePlotter) texturePlotter = new TexturePlotter();
			texturePlotter->Display(&worker);
			break;
		case MENU_TOOLS_SPECTRUMPLOTTER:
			if (!spectrumPlotter) spectrumPlotter = new SpectrumPlotter();
			spectrumPlotter->Display(&worker);
			break;

		case MENU_FACET_MESH:
			if (!facetMesh) facetMesh = new FacetMesh();
			facetMesh->EditFacet(&worker);
			UpdateFacetParams();
			break;
			/*case MENU_FACET_OUTGASSINGMAP:
				if( !outgassingMap ) outgassingMap = new OutgassingMap();
				outgassingMap->Display(&worker);
				break;*/
		case MENU_FACET_REMOVESEL:
			if (GLMessageBox::Display("Remove selected facets?", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
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
					try { worker.Reload(); }
					catch (Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(), "Error reloading worker", GLDLG_OK, GLDLG_ICONERROR);
					}
				}
			}
			break;

		case MENU_FACET_DETAILS:
			if (facetDetails == NULL) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
			break;

		case MENU_REGIONS_REGIONINFO:
			if ((int)worker.regions.size() > 0) {
				if (regionInfo == NULL) regionInfo = new RegionInfo(&worker);
				regionInfo->SetVisible(true);
			}
			else {
				GLMessageBox::Display("No regions loaded", "Error", GLDLG_OK, GLDLG_ICONERROR);
			}
			break;

		case MENU_FACET_SELECTSTICK:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.sticking != 0.0 && !geom->GetFacet(i)->IsTXTLinkFacet())
					geom->SelectFacet(i);
			geom->UpdateSelection();
			UpdateFacetParams(true);
			break;

		case MENU_FACET_SELECTREFL:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++) {
				Facet *f = geom->GetFacet(i);
				if (f->sh.sticking == 0.0 && f->sh.opacity > 0.0)
					geom->SelectFacet(i);
			}
			geom->UpdateSelection();
			UpdateFacetParams(true);
			break;

		case MENU_FACET_SELECTVOL:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.isVolatile)
					geom->SelectFacet(i);
			geom->UpdateSelection();
			UpdateFacetParams(true);
			break;


		case MENU_FACET_SELECTSPECTRUM:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.hasSpectrum)
					geom->SelectFacet(i);
			geom->UpdateSelection();
			UpdateFacetParams(true);
			break;

		case MENU_VERTEX_REMOVE:
			if (geom->IsLoaded()) {
				if (GLMessageBox::Display("Remove Selected vertices?\nNote: It will also affect facets that contain them!", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
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
						try { worker.Reload(); }
						catch (Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(), "Error reloading worker", GLDLG_OK, GLDLG_ICONERROR);
						}
					}
				}

			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;

		case MENU_REGIONS_SHOWALL:
			for (size_t i = 0; i < worker.regions.size(); i++)
				worker.regions[i].params.showPhotons = true;
			worker.ChangeSimuParams();
			for (size_t i = 0; i < worker.regions.size(); i++)
				ShowHitsMenu->SetCheck(MENU_REGIONS_SHOWHITS + i, worker.regions[i].params.showPhotons);
			break;

		case MENU_REGIONS_SHOWNONE:
			for (size_t i = 0; i < worker.regions.size(); i++)
				worker.regions[i].params.showPhotons = false;
			worker.ChangeSimuParams();
			for (size_t i = 0; i < worker.regions.size(); i++)
				ShowHitsMenu->SetCheck(MENU_REGIONS_SHOWHITS + i, worker.regions[i].params.showPhotons);
			break;
		}


		// Load recent PAR menu
		if (src->GetId() >= MENU_REGIONS_LOADRECENT && src->GetId() < MENU_REGIONS_LOADRECENT + nbRecentPAR) {
			if (AskToReset()) {
				if (worker.running) worker.Stop_Public();
				LoadParam(recentPARs[src->GetId() - MENU_REGIONS_LOADRECENT]);
			}
		}

		// Load PAR to... menu
		else if (src->GetId() >= MENU_REGIONS_LOADTO && src->GetId() < MENU_REGIONS_LOADTO + Min((int)worker.regions.size(), 19)) {
			if (AskToReset()) {
				if (worker.running) worker.Stop_Public();
				LoadParam(NULL, src->GetId() - MENU_REGIONS_LOADTO);
			}
		}
		// Remove region menu
		else if (src->GetId() >= MENU_REGIONS_REMOVE && src->GetId() < MENU_REGIONS_REMOVE + Min((int)worker.regions.size(), 19)) {
			if (AskToReset()) {
				if (worker.running) worker.Stop_Public();
				RemoveRegion(src->GetId() - MENU_REGIONS_REMOVE);
			}
		}

		//View hits menu
		else if (src->GetId() >= MENU_REGIONS_SHOWHITS && src->GetId() < MENU_REGIONS_SHOWHITS + Min((int)worker.regions.size(), 19)) {
			size_t index = src->GetId() - MENU_REGIONS_SHOWHITS;
			worker.regions[index].params.showPhotons = !ShowHitsMenu->GetCheck(MENU_REGIONS_SHOWHITS + index);
			worker.ChangeSimuParams();
			ShowHitsMenu->SetCheck(MENU_REGIONS_SHOWHITS + index, worker.regions[index].params.showPhotons);
		}

		//TEXT --------------------------------------------------------------------
	case MSG_TEXT_UPD:
		if (src == facetSticking) {
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetRMSroughness || src == facetAutoCorrLength) {
			facetDoScattering->SetState(1);
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetTeleport) {
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetOpacity) {
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetSuperDest || src == facetSILabel) {
			facetApplyBtn->SetEnabled(true);
		}
		break;

	case MSG_TEXT:
		if (src == facetSticking || src == facetRMSroughness || src == facetAutoCorrLength
			|| src == facetTeleport || src == facetOpacity || src == facetSuperDest || src == facetSILabel) {
			ApplyFacetParams();
		}
		break;

		//COMBO -------------------------------------------------------------------
	case MSG_COMBO:
		if (src == facetReflType) {
			facetApplyBtn->SetEnabled(true);
			bool isMaterial = facetReflType->GetSelectedIndex() >= 2;
			if (isMaterial) facetSticking->SetText("");
			facetSticking->SetVisible(!isMaterial);
			bool isDiffuse = facetReflType->GetSelectedIndex() == 0;
			if (isDiffuse) facetDoScattering->SetState(0);  //Diffuse surface overrides rough scattering model
			facetDoScattering->SetEnabled(!isDiffuse); //Diffuse surface overrides rough scattering model
			bool scatter = (facetDoScattering->GetState() != 0);
			facetRMSroughness->SetEditable(!isDiffuse && scatter);
			facetAutoCorrLength->SetEditable(!isDiffuse && scatter);
		}
		else if (src == facetRecType || src == facetSideType) {
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == modeCombo) {
			/*if (!AskToReset()) {
				modeCombo->SetSelectedIndex(worker.generation_mode);
				return;
			}
			changedSinceSave = true;
			worker.generation_mode = modeCombo->GetSelectedIndex(); //fluxwise or powerwise
			// Send to sub process
			worker.Reload();
			UpdateFacetHits();
			*/
			worker.generation_mode = modeCombo->GetSelectedIndex(); //fluxwise or powerwise
			worker.ChangeSimuParams();
		}
		break;

		//TOGGLE ------------------------------------------------------------------
	case MSG_TOGGLE:
		// Update viewer flags
		if (src == facetDoScattering) {
			bool scatter = (facetDoScattering->GetState() != 0);
			facetRMSroughness->SetEditable(scatter);
			facetAutoCorrLength->SetEditable(scatter);
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetSpectrumToggle)
			facetApplyBtn->SetEnabled(true);
		else if (src == autoFrameMoveToggle) {
			autoFrameMove = autoFrameMoveToggle->GetState();
			forceFrameMoveButton->SetEnabled(!autoFrameMove);
		}
		else UpdateViewerFlags(); //Viewer flags clicked
		break;

		//BUTTON ------------------------------------------------------------------
	case MSG_BUTTON:
		if (src == startSimu) {
			changedSinceSave = true;
			StartStopSimulation();
			resetSimu->SetEnabled(!worker.running);
		}

		else if (src == facetApplyBtn) {
			changedSinceSave = true;
			ApplyFacetParams();
		}
		else if (src == facetDetailsBtn) {
			if (facetDetails == NULL) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
		}
		else if (src == facetCoordBtn) {
			if (!facetCoordinates) facetCoordinates = new FacetCoordinates();
			facetCoordinates->Display(&worker);
		}
		else if (src == facetTexBtn) {
			if (!facetMesh) facetMesh = new FacetMesh();
			facetMesh->EditFacet(&worker);
			changedSinceSave = true;
			UpdateFacetParams();
		}
		else if (src == showMoreBtn) {
			if (!viewer3DSettings) viewer3DSettings = new Viewer3DSettings();
			viewer3DSettings->SetVisible(!viewer3DSettings->IsVisible());
			viewer3DSettings->Reposition();
			viewer3DSettings->Refresh(geom, viewer[curViewer]);
		}

		/*else {
			ProcessFormulaButtons(src);
		}*/
		break;
	}
}

void SynRad::QuickPipe() {

	BuildPipe(5.0, 5);
}

void SynRad::BuildPipe(double ratio, int steps) {

	char tmp[128];
	SynradGeometry *geom = worker.GetSynradGeometry();

	double R = 1.0;
	double L = ratio * R;
	int    step;

	if (steps) step = steps; //Quick Pipe
	else {
		sprintf(tmp, "100");
		char *nbF = GLInputBox::GetInput(tmp, "Number of facet", "Build Pipe");
		if (!nbF) return;
		if ((sscanf(nbF, "%d", &step) <= 0) || (step < 3)) {
			GLMessageBox::Display("Invalid number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}

	ResetSimulation(false);
	ClearFormulas();
	ClearAllSelections();
	ClearAllViews();
	ClearRegions();
	geom->BuildPipe(L, R, 0, step);
	worker.nbDesorption = 0;
	sprintf(tmp, "L|R %g", L / R);
	nbDesStart = 0;
	nbHitStart = 0;
	for (int i = 0; i < MAX_VIEWER; i++)
		viewer[i]->SetWorker(&worker);
	startSimu->SetEnabled(true);
	ClearFacetParams();
	//UpdatePlotters();
	if (profilePlotter) profilePlotter->Reset();
	if (spectrumPlotter) spectrumPlotter->Reset();

	if (textureSettings) textureSettings->Update();
	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (vertexCoordinates) vertexCoordinates->Update();
	UpdateStructMenu();
	// Send to sub process
	try { worker.Reload(); }
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	UpdateTitle();
	changedSinceSave = false;
	ResetAutoSaveTimer();
}

void SynRad::RebuildPARMenus() {
	PARloadToMenu->Clear();
	PARremoveMenu->Clear();
	ShowHitsMenu->Clear(); ShowHitsMenu->Add("Show All", MENU_REGIONS_SHOWALL); ShowHitsMenu->Add("Show None", MENU_REGIONS_SHOWNONE); ShowHitsMenu->Add(NULL);
	for (int i = 0; i < Min((int)worker.regions.size(), 19); i++) {
		std::ostringstream tmp;
		tmp << "Region " << (i + 1);
		if (worker.regions[i].fileName.length() > 0) {
			tmp << "(" << worker.regions[i].fileName << ")";
		}
		PARloadToMenu->Add(tmp.str().c_str(), MENU_REGIONS_LOADTO + i);
		PARremoveMenu->Add(tmp.str().c_str(), MENU_REGIONS_REMOVE + i);
		ShowHitsMenu->Add(tmp.str().c_str(), MENU_REGIONS_SHOWHITS + i);
		ShowHitsMenu->SetCheck(MENU_REGIONS_SHOWHITS + i, worker.regions[i].params.showPhotons);
	}
}

void SynRad::RemoveRecentPAR(char *fileName) {

	if (!fileName) return;

	bool found = false;
	int i = 0;
	while (!found && i < nbRecentPAR) {
		found = strcmp(fileName, recentPARs[i]) == 0;
		if (!found) i++;
	}
	if (!found) return;

	SAFE_FREE(recentPARs[i]);
	for (int j = i; j < nbRecentPAR - 1; j++)
		recentPARs[j] = recentPARs[j + 1];
	nbRecentPAR--;

	// Update menu
	GLMenu *m = menu->GetSubMenu("Regions")->GetSubMenu("Load recent");
	m->Clear();
	for (i = nbRecentPAR - 1; i >= 0; i--)
		m->Add(recentPARs[i], MENU_REGIONS_LOADRECENT + i);
	SaveConfig();
}

void SynRad::AddRecentPAR(char *fileName) {

	// Check if already exists
	bool found = false;
	int i = 0;
	while (!found && i < nbRecentPAR) {
		found = strcmp(fileName, recentPARs[i]) == 0;
		if (!found) i++;
	}
	if (found) {
		for (int j = i; j < nbRecentPAR - 1; j++) {
			recentPARs[j] = recentPARs[j + 1];
		}
		recentPARs[nbRecentPAR - 1] = _strdup(fileName);
		UpdateRecentPARMenu();
		SaveConfig();
		return;
	}

	// Add the new recent file
	if (nbRecentPAR < MAX_RECENT) {
		recentPARs[nbRecentPAR] = _strdup(fileName);
		nbRecentPAR++;
	}
	else {
		// Shift
		SAFE_FREE(recentPARs[0]);
		for (int i = 0; i < MAX_RECENT - 1; i++)
			recentPARs[i] = recentPARs[i + 1];
		recentPARs[MAX_RECENT - 1] = _strdup(fileName);
	}

	UpdateRecentPARMenu();
	SaveConfig();
}


//-----------------------------------------------------------------------------

void SynRad::LoadConfig() {

	FileReader *f = NULL;
	char *w;
	nbRecent = 0;

	try {

		f = new FileReader("synrad.cfg");
		SynradGeometry *geom = worker.GetSynradGeometry();

		f->ReadKeyword("showRules"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showRule = f->ReadInt();
		f->ReadKeyword("showNormals"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showNormal = f->ReadInt();
		f->ReadKeyword("showUV"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showUV = f->ReadInt();
		f->ReadKeyword("showLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showLine = f->ReadInt();
		f->ReadKeyword("showLeaks"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showLeak = f->ReadInt();
		f->ReadKeyword("showHits"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHit = f->ReadInt();
		f->ReadKeyword("showVolume"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showVolume = f->ReadInt();
		f->ReadKeyword("showTexture"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showTexture = f->ReadInt();
		f->ReadKeyword("showFilter"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showFilter = f->ReadInt();
		f->ReadKeyword("showIndices"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showIndex = f->ReadInt();
		f->ReadKeyword("showVertices"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showVertex = f->ReadInt();
		f->ReadKeyword("showMode"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showBack = f->ReadInt();
		f->ReadKeyword("showMesh"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showMesh = f->ReadInt();
		f->ReadKeyword("showHidden"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHidden = f->ReadInt();
		f->ReadKeyword("showHiddenVertex"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHiddenVertex = f->ReadInt();
		f->ReadKeyword("texColormap"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			//viewer[i]->showColormap = 
			f->ReadInt();
		f->ReadKeyword("translation"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->transStep = f->ReadDouble();
		f->ReadKeyword("dispNumLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumHits = f->ReadLLong();
		f->ReadKeyword("dispNumLeaks"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumLeaks = f->ReadLLong();
		f->ReadKeyword("dispNumTraj"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumTraj = f->ReadInt();
		f->ReadKeyword("angle"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->angleStep = f->ReadDouble();
		f->ReadKeyword("autoScale"); f->ReadKeyword(":");
		geom->texAutoScale = f->ReadInt();
		f->ReadKeyword("texMin_MC"); f->ReadKeyword(":");
		geom->texMin_MC = (llong)f->ReadDouble();
		f->ReadKeyword("texMax_MC"); f->ReadKeyword(":");
		geom->texMax_MC = (llong)f->ReadDouble();
		f->ReadKeyword("texMin_flux"); f->ReadKeyword(":");
		geom->texMin_flux = f->ReadDouble();
		f->ReadKeyword("texMax_flux"); f->ReadKeyword(":");
		geom->texMax_flux = f->ReadDouble();
		f->ReadKeyword("texMin_power"); f->ReadKeyword(":");
		geom->texMin_power = f->ReadDouble();
		f->ReadKeyword("texMax_power"); f->ReadKeyword(":");
		geom->texMax_power = f->ReadDouble();
		f->ReadKeyword("processNum"); f->ReadKeyword(":");
		nbProc = f->ReadInt();
#ifdef _DEBUG
		nbProc = 1;
#endif
		if (nbProc <= 0) nbProc = 1;
		f->ReadKeyword("recents"); f->ReadKeyword(":"); f->ReadKeyword("{");
		w = f->ReadString();
		while (strcmp(w, "}") != 0 && nbRecent < MAX_RECENT) {
			recents[nbRecent] = _strdup(w);
			nbRecent++;
			w = f->ReadString();
		}

		f->ReadKeyword("recentPARs"); f->ReadKeyword(":"); f->ReadKeyword("{");
		w = f->ReadString();
		while (strcmp(w, "}") != 0 && nbRecentPAR < MAX_RECENT) {
			recentPARs[nbRecentPAR] = _strdup(w);
			nbRecentPAR++;
			w = f->ReadString();
		}
		f->ReadKeyword("cdir"); f->ReadKeyword(":");
		strcpy(currentDir, f->ReadString());
		f->ReadKeyword("cseldir"); f->ReadKeyword(":");
		strcpy(currentSelDir, f->ReadString());
		f->ReadKeyword("autonorme"); f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("centernorme"); f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		f->ReadKeyword("normeratio"); f->ReadKeyword(":");
		geom->SetNormeRatio((float)(f->ReadDouble()));
		f->ReadKeyword("showDirection"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showDir = f->ReadInt();
		f->ReadKeyword("shadeLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->shadeLines = f->ReadInt();
		f->ReadKeyword("showTP"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showTP = f->ReadInt();
		f->ReadKeyword("autoSaveFrequency"); f->ReadKeyword(":");
		autoSaveFrequency = f->ReadDouble();
		f->ReadKeyword("autoSaveSimuOnly"); f->ReadKeyword(":");
		autoSaveSimuOnly = f->ReadInt();
		f->ReadKeyword("checkForUpdates"); f->ReadKeyword(":");
		checkForUpdates = f->ReadInt();
		f->ReadKeyword("autoUpdateFormulas"); f->ReadKeyword(":");
		autoUpdateFormulas = f->ReadInt();
		f->ReadKeyword("compressSavedFiles"); f->ReadKeyword(":");
		compressSavedFiles = f->ReadInt();
		f->ReadKeyword("lowFluxMode"); f->ReadKeyword(":");
		worker.lowFluxMode = f->ReadInt();
		f->ReadKeyword("lowFluxCutoff"); f->ReadKeyword(":");
		worker.lowFluxCutoff = f->ReadDouble();
		f->ReadKeyword("textureLogScale"); f->ReadKeyword(":");
		geom->texLogScale = f->ReadInt();
		f->ReadKeyword("newReflectionModel"); f->ReadKeyword(":");
		worker.newReflectionModel = f->ReadInt();
		PlaceScatteringControls(worker.newReflectionModel);
		f->ReadKeyword("hideLot"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->hideLot = f->ReadInt();
	}
	catch (Error &err) {
		printf("Warning, load config file (one or more feature not supported) %s\n", err.GetMsg());
	}

	SAFE_DELETE(f);


}

//-----------------------------------------------------------------------------

#define WRITEI(name,var) {             \
	f->Write(name);                      \
	f->Write(":");                       \
	for(int i=0;i<MAX_VIEWER;i++)        \
	f->Write(viewer[i]->var," ");   \
	f->Write("\n");                      \
}                                      \

#define WRITED(name,var) {             \
	f->Write(name);                      \
	f->Write(":");                       \
	for(int i=0;i<MAX_VIEWER;i++)        \
	f->Write(viewer[i]->var," ");\
	f->Write("\n");                      \
}


void SynRad::SaveConfig(bool increaseSessionCount) {

	FileWriter *f = NULL;

	try {

		f = new FileWriter("synrad.cfg");
		SynradGeometry *geom = worker.GetSynradGeometry();

		// Save flags
		WRITEI("showRules", showRule);
		WRITEI("showNormals", showNormal);
		WRITEI("showUV", showUV);
		WRITEI("showLines", showLine);
		WRITEI("showLeaks", showLeak);
		WRITEI("showHits", showHit);
		WRITEI("showVolume", showVolume);
		WRITEI("showTexture", showTexture);
		WRITEI("showFilter", showFilter);
		WRITEI("showIndices", showIndex);
		WRITEI("showVertices", showVertex);
		WRITEI("showMode", showBack);
		WRITEI("showMesh", showMesh);
		WRITEI("showHidden", showHidden);
		WRITEI("showHiddenVertex", showHiddenVertex);
		//WRITEI("texColormap", showColormap);
		f->Write("texColormap:1 1 1 1\n");
		WRITED("translation", transStep);
		WRITEI("dispNumLines", dispNumHits);
		WRITEI("dispNumLeaks", dispNumLeaks);
		WRITEI("dispNumTraj", dispNumTraj);
		WRITED("angle", angleStep);
		f->Write("autoScale:"); f->Write(geom->texAutoScale, "\n");
		f->Write("texMin_MC:"); f->Write(geom->texMin_MC, "\n");
		f->Write("texMax_MC:"); f->Write(geom->texMax_MC, "\n");
		f->Write("texMin_flux:"); f->Write(geom->texMin_flux, "\n");
		f->Write("texMax_flux:"); f->Write(geom->texMax_flux, "\n");
		f->Write("texMin_power:"); f->Write(geom->texMin_power, "\n");
		f->Write("texMax_power:"); f->Write(geom->texMax_power, "\n");
#ifdef _DEBUG
		f->Write("processNum:"); f->Write(numCPU, "\n");
#else
		f->Write("processNum:"); f->Write(worker.GetProcNumber(), "\n");
#endif
		f->Write("recents:{\n");
		for (int i = 0; i < nbRecent; i++) {
			f->Write("\"");
			f->Write(recents[i]);
			f->Write("\"\n");
		}
		f->Write("}\n");
		f->Write("recentPARs:{\n");
		for (int i = 0; i < nbRecentPAR; i++) {
			f->Write("\"");
			f->Write(recentPARs[i]);
			f->Write("\"\n");
		}
		f->Write("}\n");
		f->Write("cdir:\""); f->Write(currentDir); f->Write("\"\n");
		f->Write("cseldir:\""); f->Write(currentSelDir); f->Write("\"\n");
		f->Write("autonorme:"); f->Write(geom->GetAutoNorme(), "\n");
		f->Write("centernorme:"); f->Write(geom->GetCenterNorme(), "\n");
		f->Write("normeratio:"); f->Write((double)(geom->GetNormeRatio()), "\n");
		WRITEI("showDirection", showDir); f->Write("\n");
		WRITEI("shadeLines", showDir); f->Write("\n");
		WRITEI("showTP", showDir); f->Write("\n");
		f->Write("autoSaveFrequency:"); f->Write(autoSaveFrequency, "\n");
		f->Write("autoSaveSimuOnly:"); f->Write(autoSaveSimuOnly, "\n");
		f->Write("checkForUpdates:"); f->Write(checkForUpdates, "\n");
		f->Write("autoUpdateFormulas:"); f->Write(autoUpdateFormulas, "\n");
		f->Write("compressSavedFiles:"); f->Write(compressSavedFiles, "\n");
		f->Write("lowFluxMode:"); f->Write(worker.lowFluxMode, "\n");
		f->Write("lowFluxCutoff:"); f->Write(worker.lowFluxCutoff, "\n");
		f->Write("textureLogScale:"); f->Write(geom->texLogScale, "\n");
		f->Write("newReflectionModel:"); f->Write(worker.newReflectionModel, "\n");

		WRITEI("hideLot", hideLot);
	}
	catch (Error &err) {
		printf("Warning, failed to save config file %s\n", err.GetMsg());
	}

	SAFE_DELETE(f);

}

void SynRad::CrashHandler(Error *e) {
	char tmp[1024];
	sprintf(tmp, "Well, that's emberassing. Synrad crashed and will exit now.\nBefore that, an autosave will be attempted.\nHere is the error info:\n\n%s", (char *)e->GetMsg());
	GLMessageBox::Display(tmp, "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	try {
		AutoSave(true); //crashSave
		GLMessageBox::Display("Good news, autosave worked!", "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	}
	catch (Error &e) {
		e.GetMsg();
		GLMessageBox::Display("Sorry, I couldn't even autosave.", "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	}
}

void SynRad::LoadParam(char *fName, int position) {

	if (!worker.GetGeometry()->IsLoaded()) {
		GLMessageBox::Display("You have to load a geometry before adding regions.", "Add magnetic region", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	std::vector<FILENAME> files;
	if (!AskToReset()) return;
	if (fName == NULL) {
		if (position == -1)
			files = GLFileBox::OpenMultipleFiles(fileParFilters, "Add magnetic region(s)");
		else { //To given position, allow only one file
			FILENAME *file = GLFileBox::OpenFile(currentDir, NULL, "Add magnetic region", fileParFilters, 0);
			files.push_back(*file);
		}
	}
	else { //Filename already defined
		char shortName[256];
		char *lPart = strrchr(fName, '\\');
		if (lPart) strcpy(shortName, lPart + 1);
		else strcpy(shortName, fName);
		FILENAME file;
		strcpy(file.file, shortName);
		strcpy(file.fullName, fName);
		files.push_back(file);
	}

	GLProgress *progressDlg2 = new GLProgress("Preparing to load file...", "Please wait");
	progressDlg2->SetVisible(true);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if (files.size() == 0) { //Nothing selected
		progressDlg2->SetVisible(false);
		SAFE_DELETE(progressDlg2);
		return;
	}
	for (size_t i = 0; i < files.size(); i++) {
		char tmp[256];
		sprintf(tmp, "Adding %s...", files[i].file);
		progressDlg2->SetMessage(tmp);
		progressDlg2->SetProgress((double)i / (double)files.size());
		try {
			worker.AddRegion(files[i].fullName, position);
			AddRecentPAR(files[i].fullName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), files[i].file);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
		if (regionInfo) regionInfo->Update();
	}
	progressDlg2->SetVisible(false);
	SAFE_DELETE(progressDlg2);
	changedSinceSave = false;
	worker.GetGeometry()->RecalcBoundingBox(); //recalculate bounding box
	RebuildPARMenus();
}

void SynRad::ClearRegions() {
	if (!AskToReset()) return;
	worker.ClearRegions();
	changedSinceSave = true;
	worker.Reload();
	if (regionInfo) {
		regionInfo->SetVisible(false);
		SAFE_DELETE(regionInfo);
	}
	worker.GetGeometry()->RecalcBoundingBox(); //recalculate bounding box
	RebuildPARMenus();
}

void SynRad::RemoveRegion(int index) {
	if (!AskToReset()) return;
	if (regionEditor != NULL && (index <= regionEditor->GetRegionId())) { //Editing a region that changes
		regionEditor->SetVisible(false);
		SAFE_DELETE(regionEditor);
	}
	worker.RemoveRegion(index);
	changedSinceSave = true;
	worker.Reload();
	if (regionInfo) regionInfo->Update();
	RebuildPARMenus();

}

void SynRad::NewRegion() {
	if (!worker.GetGeometry()->IsLoaded()) {
		GLMessageBox::Display("You have to load a geometry before adding regions.", "Add magnetic region", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (worker.running) worker.Stop_Public();
	FILENAME *fn = GLFileBox::SaveFile(NULL, NULL, "Save Region", "param files\0*.param\0All files\0*.*\0", 2);
	if (!fn || !fn->fullName) return;
	if (!(FileUtils::GetExtension(fn->fullName) == "param"))
		sprintf(fn->fullName, "%s.param", fn->fullName); //append .param extension
	try {
		Region_full newreg;
		newreg.fileName.assign(fn->fullName);
		//worker.regions.push_back(newreg);
		FileWriter *file = new FileWriter(fn->fullName);
		newreg.SaveParam(file);
		SAFE_DELETE(file);
		AddRecentPAR(fn->fullName);
		worker.AddRegion(fn->fullName);
		RebuildPARMenus();
	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecentPAR(fn->fullName);
		return;
	}
	if (regionInfo) regionInfo->Update();
	if (mApp->regionEditor == NULL) regionEditor = new RegionEditor();
	regionEditor->Display(&worker, (int)worker.regions.size() - 1);
}

/*bool EndsWithParam(const char* s)
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

void SynRad::UpdateRecentPARMenu() {
	// Update menu
	GLMenu *m = menu->GetSubMenu("Regions")->GetSubMenu("Load recent");
	m->Clear();
	for (int i = nbRecentPAR - 1; i >= 0; i--)
		m->Add(recentPARs[i], MENU_REGIONS_LOADRECENT + i);
}

void SynRad::PlaceScatteringControls(bool newReflectionMode) {
	facetRMSroughnessLabel->SetText(newReflectionMode ? "sigma (nm):" : "Roughness ratio:");
	facetPanel->SetCompBounds(facetRMSroughness, newReflectionMode ? 65 : 100, 55, newReflectionMode ? 45 : 70, 18);
	facetAutoCorrLengthLabel->SetVisible(newReflectionMode);
	facetAutoCorrLength->SetVisible(newReflectionMode);
	UpdateFacetParams();
}