/*
File:        Worker.cpp
Description: Sub processes handling
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

#include "Worker.h"
#include "GLApp/GLApp.h"
#include "GLApp/GLMessageBox.h"
#include <math.h>
#include <stdlib.h>
#include <Process.h>
#include "GLApp/GLUnitDialog.h"
#include "Synrad.h"
#include <direct.h>

extern GLApplication *theApp;
extern HANDLE compressProcessHandle;
extern double totalOutgassing;
extern int autoSaveSimuOnly;
extern int needsReload;
extern int compressSavedFiles;
BOOL EndsWithPar(const char* s);

// -------------------------------------------------------------

Worker::Worker() {
	pid = _getpid();
	sprintf(ctrlDpName,"SRDCTRL%d",pid);
	sprintf(loadDpName,"SRDLOAD%d",pid);
	sprintf(hitsDpName,"SRDHITS%d",pid);
	sprintf(materialsDpName,"SRDMATS%d",pid);
	nbProcess = 0;
	maxDesorption = 0;
	ResetWorkerStats();
	geom = new Geometry();
	regions = std::vector<Region>();
	dpControl = NULL;
	dpHit = NULL;
	nbHHit = 0;
	nbHit = 0;
	startTime = 0.0f;
	stopTime = 0.0f;
	simuTime = 0.0f;
	nbTrajPoints=0;
	running = FALSE;
	calcAC = FALSE;
	strcpy(fullFileName,"");
}

// -------------------------------------------------------------

Worker::~Worker() {
	CLOSEDP(dpHit);
	Exit();
	delete geom;
}

// -------------------------------------------------------------

Geometry *Worker::GetGeometry() {
	return geom;
}

// -------------------------------------------------------------

char *Worker::GetFileName() {
	return fullFileName;
}

char *Worker::GetShortFileName() {

	static char ret[256];
	char *r = strrchr(fullFileName,'/');
	if(!r) r = strrchr(fullFileName,'\\');
	if(!r) strcpy(ret,fullFileName);
	else   {
		r++;
		strcpy(ret,r);
	}

	return ret;

}

// -------------------------------------------------------------

void Worker::SetFileName(char *fileName) {
	strcpy(fullFileName,fileName);
}

/*
// Execute command-line command and return the screen output
std::string execCMD(char* cmd);
*/
// -------------------------------------------------------------

void Worker::SaveGeometry(char *fileName,GLProgress *prg,BOOL askConfirm,BOOL saveSelected,BOOL autoSave,BOOL crashSave) {

	if (needsReload&&(!crashSave)) RealReload();
	char tmp[10000]; //compress.exe command line
	char fileNameWithGeo[2048]; //file name with .geo extension (instead of .geo7z)
	char fileNameWithGeo7z[2048];
	char fileNameWithSyn[2048]; //file name with .syn extension (instead of .syn7z)
	char fileNameWithSyn7z[2048];
	char fileNameWithoutExtension[2048]; //file name without extension
	//char *ext = fileName+strlen(fileName)-4;
	char *ext,*dir;
	SynRad *mApp = (SynRad *)theApp;

	dir = strrchr(fileName,'\\');
	ext = strrchr(fileName,'.');

	if(!(ext) || !(*ext=='.') || ((dir)&&(dir>ext)) ) { 
		sprintf(fileName,"%s.syn7z",fileName); //set to default SYN7Z format
		ext = strrchr(fileName,'.');
	}

	ext++;

	// Read a file
	BOOL ok = TRUE;
	FileWriter *f = NULL;
	BOOL isTXT = _stricmp(ext,"txt")==0;
	BOOL isSTR = _stricmp(ext,"str")==0;
	BOOL isGEO = _stricmp(ext,"geo")==0;
	BOOL isGEO7Z = _stricmp(ext,"geo7z")==0;
	BOOL isSYN = _stricmp(ext,"syn")==0;
	BOOL isSYN7Z = _stricmp(ext,"syn7z")==0;

	if(isTXT || isGEO || isGEO7Z || isSYN || isSYN7Z || isSTR) {

		if (WAIT_TIMEOUT==WaitForSingleObject(compressProcessHandle,0)) {
			GLMessageBox::Display("Compressing a previous save file is in progress. Wait until that finishes"
				"or close process \"compress.exe\"\nIf this was an autosave attempt,"
				"you have to lower the autosave frequency.","Can't save right now.",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
		if (isGEO || isSYN) {
			memcpy(fileNameWithoutExtension,fileName,sizeof(char)*(strlen(fileName)-4));
			fileNameWithoutExtension[strlen(fileName)-4]='\0';
		} else if (isGEO7Z || isSYN7Z) { //geo7z
			memcpy(fileNameWithoutExtension,fileName,sizeof(char)*(strlen(fileName)-6));
			fileNameWithoutExtension[strlen(fileName)-6]='\0';
		}

		if (isGEO) {
			sprintf(fileNameWithGeo7z,"%s7z",fileName);
			memcpy(fileNameWithGeo,fileName,(strlen(fileName)+1)*sizeof(char));

			if(!autoSave && FileUtils::Exist(fileNameWithGeo7z) && compressSavedFiles) {
				sprintf(tmp,"A .geo7z file of the same name exists. Overwrite that file ?\n%s",fileNameWithGeo7z);
				ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
			}
		} 
		if (isGEO7Z) {
			memcpy(fileNameWithGeo,fileName,sizeof(char)*(strlen(fileName)-2));
			fileNameWithGeo[strlen(fileName)-2]='\0';
			memcpy(fileNameWithGeo7z,fileName,(1+strlen(fileName))*sizeof(char));
			sprintf(tmp,"A .geo file of the same name exists. Overwrite that file ?\n%s",fileNameWithGeo);
			if(!autoSave && FileUtils::Exist(fileNameWithGeo) ) {
				ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
			}
		} 		
		if (isSYN) {
			sprintf(fileNameWithSyn7z,"%s7z",fileName);
			memcpy(fileNameWithSyn,fileName,(strlen(fileName)+1)*sizeof(char));

			if(!autoSave && FileUtils::Exist(fileNameWithSyn7z) && compressSavedFiles) {
				sprintf(tmp,"A .syn7z file of the same name exists. Overwrite that file ?\n%s",fileNameWithSyn7z);
				ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
			}
		} 
		if (isSYN7Z) {
			memcpy(fileNameWithSyn,fileName,sizeof(char)*(strlen(fileName)-2));
			fileNameWithSyn[strlen(fileName)-2]='\0';
			memcpy(fileNameWithSyn7z,fileName,(1+strlen(fileName))*sizeof(char));
			sprintf(tmp,"A .syn file of the same name exists. Overwrite that file ?\n%s",fileNameWithSyn);
			if(!autoSave && FileUtils::Exist(fileNameWithSyn) ) {
				ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
			}
		} 
		if(!autoSave && FileUtils::Exist(fileName) ) {
			sprintf(tmp,"Overwrite existing file ?\n%s",fileName);
			if (askConfirm) ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
		}

		if( ok ) {
			if( isSTR ) {
				geom->SaveSTR(dpHit,saveSelected);
			} else {
				try {
					if (isGEO7Z) {
						memcpy(fileNameWithGeo,fileName,sizeof(char)*(strlen(fileName)-2));
						fileNameWithGeo[strlen(fileName)-2]='\0';
						f = new FileWriter(fileNameWithGeo);
					} else if (isSYN7Z) {
						memcpy(fileNameWithSyn,fileName,sizeof(char)*(strlen(fileName)-2));
						fileNameWithSyn[strlen(fileName)-2]='\0';
						f = new FileWriter(fileNameWithSyn);
					} else
						f = new FileWriter(fileName);
				} catch(Error &e) {
					SAFE_DELETE(f);
					GLMessageBox::Display((char*)e.GetMsg(),"Error in Worker::SaveGeometry()",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				geom->tNbDesorptionMax = maxDesorption;
				if( isTXT ) geom->SaveTXT(f,dpHit,saveSelected);
				else if( isGEO || isGEO7Z ) {
					// Retrieve leak cache
					int nbLeakSave,nbHHitSave;
					LEAK pLeak[NBHLEAK];
					if (!crashSave) GetLeak(pLeak,&nbLeakSave);
					// Retrieve hit cache (lines and dots)
					HIT pHits[NBHHIT];
					if (!crashSave) GetHHit(pHits,&nbHHitSave);
					geom->SaveGEO(f,prg,dpHit,saveSelected,pLeak,&nbLeakSave,pHits,&nbHHitSave,crashSave);
				} else if( isSYN || isSYN7Z ) {
					// Retrieve leak cache
					int nbLeakSave,nbHHitSave;
					LEAK pLeak[NBHLEAK];
					if (!crashSave) GetLeak(pLeak,&nbLeakSave);
					else nbLeakSave=0;
					// Retrieve hit cache (lines and dots)
					HIT pHits[NBHHIT];
					if (!crashSave) GetHHit(pHits,&nbHHitSave);
					else nbHHitSave=0;

					//Save regions
					for (int i=0;i<(int)regions.size();i++) {
						if (EndsWithPar((char*)regions[i].fileName.c_str())) { //fileName doesn't end with .param, save a param version
							sprintf(tmp,"%sam",regions[i].fileName.c_str());
							regions[i].fileName=tmp;
						}	
						SaveRegion((char*)regions[i].fileName.c_str(),i,TRUE); //save with forced overwrite
					}
					geom->SaveSYN(f,prg,dpHit,saveSelected,pLeak,&nbLeakSave,pHits,&nbHHitSave,crashSave);
				}
			}
			if (!autoSave) strcpy(fullFileName,fileName);
			if (!autoSave) {
				remove("Synrad_AutoSave.syn");
				remove("Synrad_AutoSave.syn7z");
			}
		}
	} else {
		SAFE_DELETE(f);
		throw Error("SaveGeometry(): Invalid file extension [only geo,txt,str]");
	}

	SAFE_DELETE(f);
	if (ok && isGEO || isGEO7Z) {

		if (compressSavedFiles) {
			if (FileUtils::Exist("compress.exe")) { //compress GEO file to GEO7Z using 7-zip launcher "compress.exe"
				sprintf(tmp,"compress.exe \"%s\" Geometry.geo",fileNameWithGeo);
				int procId = StartProc_background(tmp);
				compressProcessHandle=OpenProcess(PROCESS_ALL_ACCESS, TRUE, procId);
				fileName=fileNameWithGeo7z;
			} else {
				GLMessageBox::Display("compress.exe (part of Molfow) not found.\n Will save as uncompressed GEO file.","Compressor not found",GLDLG_OK,GLDLG_ICONERROR);
				fileName=fileNameWithGeo;
			}
		} else fileName=fileNameWithGeo;
		if (!autoSave) {
			SetFileName(fileName);
			mApp->UpdateTitle();
		}
	} else if (ok && isSYN || isSYN7Z) {

		if (compressSavedFiles) {
			if (FileUtils::Exist("compress.exe")) { //compress SYN file to SYN7Z using 7-zip launcher "compress.exe"
				sprintf(tmp,"compress.exe \"%s\" Geometry.syn",fileNameWithSyn);
				for (int i=0;i<(int)regions.size();i++) {
					sprintf(tmp,"%s \"%s\"",tmp,regions[i].fileName.c_str());
					if (!regions[i].MAGXfileName.empty())
						sprintf(tmp,"%s \"%s\"",tmp,regions[i].MAGXfileName.c_str());
					if (!regions[i].MAGYfileName.empty())
						sprintf(tmp,"%s \"%s\"",tmp,regions[i].MAGYfileName.c_str());
					if (!regions[i].MAGZfileName.empty())
						sprintf(tmp,"%s \"%s\"",tmp,regions[i].MAGZfileName.c_str());
					if (!regions[i].BXYfileName.empty())
						sprintf(tmp,"%s \"%s\"",tmp,regions[i].BXYfileName.c_str());
				}
				int procId = StartProc_background(tmp);
				compressProcessHandle=OpenProcess(PROCESS_ALL_ACCESS, TRUE, procId);
				fileName=fileNameWithSyn7z;
			} else {
				GLMessageBox::Display("compress.exe (part of Molfow) not found.\n Will save as uncompressed SYN file.","Compressor not found",GLDLG_OK,GLDLG_ICONERROR);
				fileName=fileNameWithSyn;
			}
		} else fileName=fileNameWithSyn;
		if (!autoSave) {
			SetFileName(fileName);
			mApp->UpdateTitle();
		}
	}
}

void Worker::ExportTextures(char *fileName,int mode,BOOL askConfirm,BOOL saveSelected) {

	char tmp[512];

	// Read a file
	FILE *f = NULL;
	BOOL ok = TRUE;
	if( askConfirm ) {
		if( FileUtils::Exist(fileName) ) {
			sprintf(tmp,"Overwrite existing file ?\n%s",fileName);
			ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
		}
	}
	if( ok ) {
		f=fopen(fileName,"w");
		geom->SaveTexture(f,mode,dpHit,saveSelected);
		fclose(f);
	}
}

void Worker::ExportRegionPoints(char *fileName,GLProgress *prg,int regionId,int exportFrequency,BOOL doFullScan){
	char tmp[512];
	// Read a file
	FileWriter *f = NULL;
	BOOL ok = TRUE;
	char *ext,*dir;

	dir = strrchr(fileName,'\\');
	ext = strrchr(fileName,'.');

	if(!(ext) || !(*ext=='.') || ((dir)&&(dir>ext)) ) { 
		sprintf(fileName,"%s.csv",fileName); //set to default CSV extension
		ext = strrchr(fileName,'.');
	}
	if( FileUtils::Exist(fileName) ) {
		sprintf(tmp,"Overwrite existing file ?\n%s",fileName);
		ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
	}
	if( ok ) {
		f= new FileWriter(fileName);
		regions[regionId].ExportPoints(f,prg,exportFrequency,doFullScan);
		SAFE_DELETE(f);
	}
}

void Worker::ExportDesorption(char *fileName,bool selectedOnly,int mode,double eta0,double alpha,Distribution2D *distr) {

	// Read a file
	FILE *f = NULL;
	BOOL ok = TRUE;

	char *ext,*dir;

	dir = strrchr(fileName,'\\');
	ext = strrchr(fileName,'.');

	if(!(ext) || !(*ext=='.') || ((dir)&&(dir>ext)) ) { 
		sprintf(fileName,"%s.des",fileName); //set to default DES extension
		ext = strrchr(fileName,'.');
	}
	char tmp[512];
	sprintf(tmp,"A .geo file of the same name exists. Overwrite that file ?\n%s",fileName);
	if(FileUtils::Exist(fileName) ) {
		ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
	}
	if (!ok) return;
	f=fopen(fileName,"w");
	if (needsReload) RealReload();
	geom->SaveDesorption(f,dpHit,selectedOnly,mode,eta0,alpha,distr);
	fclose(f);

}

// -------------------------------------------------------------

void Worker::LoadGeometry(char *fileName) {
	needsReload=TRUE;
	SynRad *mApp = (SynRad *)theApp;
	char *ext,*filebegin;
	BOOL isGEO7Z,isSYN7Z;
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );

	ext = strrchr(fileName,'.');

	if( ext==NULL )
		throw Error("LoadGeometry(): Invalid file extension [Only syn,syn7z,txt,stl,str,geo,geo7z or ase]");
	ext++;

	// Read a file
	FileReader *f = NULL;
	GLProgress *progressDlg = new GLProgress("Reading file...","Please wait");
	progressDlg->SetVisible(TRUE);
	progressDlg->SetProgress(0.0);
	isGEO7Z=(_stricmp(ext,"geo7z")==0);
	isSYN7Z=(_stricmp(ext,"syn7z")==0);
	if(_stricmp(ext,"txt")==0) {

		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			geom->LoadTXT(f,progressDlg);
			nbHit = geom->tNbHit;
			nbDesorption = geom->tNbDesorption;
			//totalFlux = geom->tFlux;
			//totalPower = geom->tPower;
			maxDesorption = geom->tNbDesorptionMax;
			nbLeak = geom->tNbLeak;
			//RealReload();
			strcpy(fullFileName,fileName);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"stl")==0) {
		try {
			int ret = GLUnitDialog::Display("","Choose STL file units:",GLDLG_MM | GLDLG_CM| GLDLG_M| GLDLG_INCH| GLDLG_FOOT | GLDLG_CANCEL_U,GLDLG_ICONNONE);
			double scaleFactor=1.0;
			switch(ret) {
			case GLDLG_MM:
				scaleFactor=0.1;
				break;
			case GLDLG_CM:
				scaleFactor=1.0;
				break;
			case GLDLG_M:
				scaleFactor=100;
				break;
			case GLDLG_INCH:
				scaleFactor=2.54;
				break;
			case GLDLG_FOOT:
				scaleFactor=30.48;
				break;
			}
			if (ret!=GLDLG_CANCEL_U) {
				progressDlg->SetMessage("Resetting worker...");
				ResetWorkerStats();
				f = new FileReader(fileName);
				progressDlg->SetVisible(TRUE);
				progressDlg->SetMessage("Reading geometry...");
				geom->LoadSTL(f,progressDlg,scaleFactor);
				//progressDlg->SetMessage("Reloading worker...");
				//RealReload();
				strcpy(fullFileName,fileName);
				mApp->DisplayCollapseDialog();
			}
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;

		}

	} else if(_stricmp(ext,"str")==0) {

		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadSTR(f,progressDlg);
			//RealReload();
			strcpy(fullFileName,fileName);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"syn")==0 || isSYN7Z) {
		char tmp2[2048];
		int version;
		progressDlg->SetVisible(TRUE);
		try {
			if (isSYN7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[2048];
				char fileOnly[512];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				//execCMD(tmp);
				system(tmp);

				filebegin= strrchr(fileName,'\\');
				if (filebegin) filebegin++;
				else filebegin=fileName;
				memcpy(fileOnly,filebegin,sizeof(char)*(strlen(filebegin)-2));
				fileOnly[strlen(filebegin)-2]='\0';
				sprintf(tmp2,"%s\\tmp\\Geometry.syn",CWD);
				f = new FileReader(tmp2);
			}
			progressDlg->SetMessage("Resetting worker...");
			ResetWorkerStats();
			if (!isSYN7Z) f = new FileReader(fileName);

			//leaks
			LEAK pLeak[NBHLEAK];
			//hits
			HIT pHits[NBHHIT];

			PARfileList regionsToLoad=geom->LoadSYN(f,progressDlg,pLeak,&nbLeak,pHits,&nbHHit,&version);
			nbLeak = geom->tNbLeak;
			nbHit = geom->tNbHit;
			nbDesorption = geom->tNbDesorption;
			maxDesorption = geom->tNbDesorptionMax;
			totalFlux = geom->tFlux;
			totalPower = geom->tPower;

			//Load regions
			if (regionsToLoad.nbFiles>0) {
				char tmp[256];
				sprintf(tmp,"This geometry refers to %d regions. Load them now?",regionsToLoad.nbFiles);
				BOOL loadThem = ( GLMessageBox::Display(tmp,"File load",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK );
				//BOOL loadThem = TRUE;
				if (loadThem) {
					for (int i=0;i<regionsToLoad.nbFiles;i++) {
						std::string toLoad;
						if (isSYN7Z) { //PAR file to load in tmp dir (just extracted)
							toLoad=CWD;
							toLoad.append("\\tmp\\");
						} else { //PAR file in same dir as SYN file
							char tmp[512];
							filebegin= strrchr(fileName,'\\'); //throw out absolute path, keep only filename
							if (filebegin) filebegin++;
							else filebegin=fileName;
							memcpy(tmp,fileName,filebegin-fileName);
							tmp[(int)(filebegin-fileName)]=NULL;
							toLoad=tmp;
						}
						toLoad.append(regionsToLoad.fileNames[i]);
						char *toLoadChar=_strdup(toLoad.c_str());
						AddRegion(toLoadChar,-1);
					}
				}
			}

			progressDlg->SetMessage("Reloading worker with new geometry...");
			RealReload(); //for the loading of textures
			geom->LoadProfileSYN(f,dpHit);
			geom->LoadSpectrumSYN(f,dpHit);
			SHGHITS *gHits = (SHGHITS *)dpHit->buff;
			SetLeak(pLeak,&nbLeak,gHits);
			SetHHit(pHits,&nbHHit,gHits);
			SAFE_DELETE(f);
			progressDlg->SetMessage("Loading textures...");
			loadTextures((isSYN7Z)?tmp2:fileName,version);
			strcpy(fullFileName,fileName);

			//if (isSYN7Z) remove(tmp2);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			//if (isSYN7Z) remove(tmp2);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"geo")==0 || isGEO7Z) {
		char tmp2[2048];
		int version;
		progressDlg->SetVisible(TRUE);
		try {
			if (isGEO7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[2048];
				char fileOnly[512];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				//execCMD(tmp);
				system(tmp);

				char *filebegin= strrchr(fileName,'\\');
				if (filebegin) filebegin++;
				else filebegin=fileName;
				memcpy(fileOnly,filebegin,sizeof(char)*(strlen(filebegin)-2));
				fileOnly[strlen(filebegin)-2]='\0';
				sprintf(tmp2,"%s\\tmp\\Geometry.geo",CWD);
				if (!FileUtils::Exist(tmp2)) //fall back to old geo7z format
					sprintf(tmp2,"%s\\tmp\\%s",CWD,fileOnly);
				f = new FileReader(tmp2);
			}
			progressDlg->SetMessage("Resetting worker...");
			ResetWorkerStats();
			if (!isGEO7Z) f = new FileReader(fileName);

			//leaks
			LEAK pLeak[NBHLEAK];
			//hits
			HIT pHits[NBHHIT];

			geom->LoadGEO(f,progressDlg,pLeak,&nbLeak,pHits,&nbHHit,&version);
			nbLeak = 0;
			nbHit = 0;
			nbDesorption = 0;
			maxDesorption = geom->tNbDesorptionMax;

			progressDlg->SetMessage("Reloading worker with new geometry...");
			RealReload(); //for the loading of textures
			//if (version>=8) geom->LoadProfileGEO(f);
			//SHGHITS *gHits = (SHGHITS *)dpHit->buff;
			//SetLeak(pLeak,&nbLeak,gHits);
			//SetHHit(pHits,&nbHHit,gHits);
			SAFE_DELETE(f);
			//progressDlg->SetMessage("Loading textures...");
			//loadTextures((isGEO7Z)?tmp2:fileName,version);
			strcpy(fullFileName,fileName);
			//if (isGEO7Z) remove(tmp2);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			//if (isGEO7Z) remove(tmp2);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"ase")==0) {

		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadASE(f,progressDlg);
			//RealReload();
			strcpy(fullFileName,fileName);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else {

		SAFE_DELETE(f);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("LoadGeometry(): Invalid file extension [Only txt,geo,geo7z,ase,stl or str]");

	}
	//geom->CalcTotalOutGassing();
	
	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
	SAFE_DELETE(f);

}

// -------------------------------------------------------------
// -------------------------------------------------------------

void Worker::InsertGeometry(BOOL newStr,char *fileName) {
	SynRad *mApp = (SynRad *)theApp;
	if (needsReload) RealReload();
	char *ext,*filebegin;
	BOOL isGEO7Z,isSYN7Z;
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );

	ext = strrchr(fileName,'.');

	if( ext==NULL )
		throw Error("InsertGeometry(): Invalid file extension [Only txt,stl,geo,geo7z,syn,syn7z]");
	ext++;

	// Read a file
	FileReader *f = NULL;

	GLProgress *progressDlg = new GLProgress("Loading file...","Please wait");
	progressDlg->SetProgress(0.0);
	progressDlg->SetVisible(TRUE);
	isGEO7Z=(_stricmp(ext,"geo7z")==0);
	isSYN7Z=(_stricmp(ext,"syn7z")==0);
	if(_stricmp(ext,"txt")==0) {

		try {
			f=new FileReader(fileName);
			changedSinceSave = TRUE;
			geom->InsertTXT(f,progressDlg,newStr);
			//geom->sh.nbSuper++;
			//nbHit = geom->tNbHit;
			//nbDesorption = geom->tNbDesorption;
			//maxDesorption = geom->tNbDesorptionMax;
			//nbLeak = geom->tNbLeak;
			//nbHit = 0;
			//nbDesorption = 0;
			//maxDesorption = 0;
			//nbLeak = 0;
			Reload();
			//strcpy(fullFileName,fileName);
		} catch(Error &e) {
			//geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"stl")==0) {
		try {
			int ret = GLUnitDialog::Display("","Choose STL file units:",GLDLG_MM | GLDLG_CM| GLDLG_M| GLDLG_INCH| GLDLG_FOOT | GLDLG_CANCEL_U,GLDLG_ICONNONE);
			double scaleFactor=1.0;
			switch(ret) {
			case GLDLG_MM:
				scaleFactor=0.1;
				break;
			case GLDLG_CM:
				scaleFactor=1.0;
				break;
			case GLDLG_M:
				scaleFactor=100;
				break;
			case GLDLG_INCH:
				scaleFactor=2.54;
				break;
			case GLDLG_FOOT:
				scaleFactor=30.48;
				break;
			}
			if (ret!=GLDLG_CANCEL_U) {
				changedSinceSave = TRUE;
				ResetWorkerStats();
				f = new FileReader(fileName);
				geom->InsertSTL(f,progressDlg,newStr);
				//geom->sh.nbSuper++;
				//nbHit = geom->tNbHit;
				//nbDesorption = geom->tNbDesorption;
				//maxDesorption = geom->tNbDesorptionMax;
				//nbLeak = geom->tNbLeak;
				//nbHit = 0;
				//nbDesorption = 0;
				//maxDesorption = 0;
				//nbLeak = 0;
				Reload();}
			//strcpy(fullFileName,fileName);
		} catch(Error &e) {
			//geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	} else if(_stricmp(ext,"str")==0) {
		/*
		try {
		ResetWorkerStats();
		f = new FileReader(fileName);
		progressDlg->SetVisible(TRUE);
		geom->LoadSTR(f,progressDlg);
		Reload();
		strcpy(fullFileName,fileName);
		} catch(Error &e) {
		geom->Clear();
		SAFE_DELETE(f);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw e;
		}
		*/

	} else if(_stricmp(ext,"syn")==0 || isSYN7Z) {
		char tmp2[2048];
		
		progressDlg->SetVisible(TRUE);
		try {
			if (isSYN7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[2048];
				char fileOnly[512];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				//execCMD(tmp);
				system(tmp);

				filebegin= strrchr(fileName,'\\');
				if (filebegin) filebegin++;
				else filebegin=fileName;
				memcpy(fileOnly,filebegin,sizeof(char)*(strlen(filebegin)-2));
				fileOnly[strlen(filebegin)-2]='\0';
				sprintf(tmp2,"%s\\tmp\\Geometry.syn",CWD);
				f = new FileReader(tmp2);
			}
			//progressDlg->SetMessage("Resetting worker...");
			//ResetWorkerStats();
			if (!isSYN7Z) f = new FileReader(fileName);


			PARfileList regionsToLoad=geom->InsertSYN(f,progressDlg,newStr);
			//nbLeak = 0;
			//nbHit = 0;
			//nbDesorption = 0;
			//maxDesorption = geom->tNbDesorptionMax;

			progressDlg->SetMessage("Reloading worker with new geometry...");
			RealReload(); //for the loading of textures
			//geom->LoadProfileSYN(f,dpHit);
			//geom->LoadSpectrumSYN(f,dpHit);
			//SHGHITS *gHits = (SHGHITS *)dpHit->buff;
			//SetLeak(pLeak,&nbLeak,gHits);
			//SetHHit(pHits,&nbHHit,gHits);
			SAFE_DELETE(f);
			//progressDlg->SetMessage("Loading textures...");
			//loadTextures((isSYN7Z)?tmp2:fileName,version);
			strcpy(fullFileName,fileName);
			//Load regions
			if (regionsToLoad.nbFiles>0) {
				char tmp[256];
				sprintf(tmp,"This geometry refers to %d regions. Load them now?",regionsToLoad.nbFiles);
				BOOL loadThem = ( GLMessageBox::Display(tmp,"File load",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK );
				if (loadThem) {
					for (int i=0;i<regionsToLoad.nbFiles;i++) {
						std::string toLoad;
						if (isSYN7Z) { //PAR file to load in tmp dir (just extracted)
							toLoad=CWD;
							toLoad.append("\\tmp\\");
						} else { //PAR file in same dir as SYN file
							char tmp[512];
							filebegin= strrchr(fileName,'\\');
							if (filebegin) filebegin++;
							else filebegin=fileName;
							memcpy(tmp,fileName,filebegin-fileName);
							tmp[(int)(filebegin-fileName)]=NULL;
							toLoad=tmp;
						}
						toLoad.append(regionsToLoad.fileNames[i]);
						char *toLoadChar=_strdup(toLoad.c_str());
						AddRegion(toLoadChar,-1);
					}
				}
			}
			//if (isSYN7Z) remove(tmp2);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			//if (isSYN7Z) remove(tmp2);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	} else if(isGEO7Z || _stricmp(ext,"geo")==0) {
		try {
			char tmp2[2048];
			if (isGEO7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[2048];
				
				char fileOnly[512];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				//execCMD(tmp);
				system(tmp);

				char *filebegin= strrchr(fileName,'\\');
				if (filebegin) filebegin++;
				else filebegin=fileName;
				memcpy(fileOnly,filebegin,sizeof(char)*(strlen(filebegin)-2));
				fileOnly[strlen(filebegin)-2]='\0';
				sprintf(tmp2,"%s\\tmp\\Geometry.geo",CWD);
				if (!FileUtils::Exist(tmp2)) //fall back to old geo7z format
					sprintf(tmp2,"%s\\tmp\\%s",CWD,fileOnly);
				f = new FileReader(tmp2);
			}
			progressDlg->SetMessage("Resetting worker...");
			ResetWorkerStats();
			if (!isGEO7Z) f = new FileReader(fileName);
			changedSinceSave = TRUE;
			geom->InsertGEO(f,progressDlg,newStr);
			//if (isGEO7Z) remove(tmp2);
			//geom->sh.nbSuper++;
			//nbHit = geom->tNbHit;
			//nbDesorption = geom->tNbDesorption;
			//maxDesorption = geom->tNbDesorptionMax;
			//nbLeak = geom->tNbLeak;
			//nbHit = 0;
			//nbDesorption = 0;
			//maxDesorption = 0;
			//nbLeak = 0;
			Reload();
			//strcpy(fullFileName,fileName);
		} catch(Error &e) {
			//geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	} else if(_stricmp(ext,"ase")==0) {
	} else {

		SAFE_DELETE(f);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("InsertGeometry(): Invalid file extension [Only txt,stl,geo,geo7z,ase,syn,syn7z or str]");

	}
	
	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
	SAFE_DELETE(f);
	mApp->UpdateFacetlistSelected();	
	mApp->UpdateViewers();

}

// -------------------------------------------------------------

void Worker::GetLeak(LEAK *buffer,int *nb) {

	*nb=0;
	if( dpHit ) {
		memcpy(buffer,leakCache,sizeof(LEAK)*NBHLEAK);
		*nb = MIN(nbLeak,NBHLEAK);
	}

}

// -------------------------------------------------------------

void Worker::SetLeak(LEAK *buffer,int *nb,SHGHITS *gHits) { //When loading from file

	if( nb>0 ) {
		memcpy(leakCache,buffer,sizeof(LEAK)*MIN(NBHLEAK,*nb));
		memcpy(gHits->pLeak,buffer,sizeof(LEAK)*MIN(NBHLEAK,*nb));
		gHits->nbLeak = *nb;
	}

}

// -------------------------------------------------------------
void Worker::loadTextures(char *fileName,int version) {
	// Read a file
	char *ext;

	ext = strrchr(fileName,'.');
	FileReader *f = NULL;

	GLProgress *progressDlg = new GLProgress("Loading textures","Please wait");
	progressDlg->SetProgress(0.0);


	if(_stricmp(ext,".syn")==0) {

		try {

			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->loadTextures(f,progressDlg,dpHit,version);
			if( AccessDataport(dpHit) ) {
				BYTE *buffer = (BYTE *)dpHit->buff;
				geom->BuildTexture(buffer);
				ReleaseDataport(dpHit);
			}
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);

		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	}

	else {

		SAFE_DELETE(f);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
	}

	SAFE_DELETE(f);

}



// -------------------------------------------------------------

void Worker::GetHHit(HIT *buffer,int *nb) {

	*nb=0;
	if( dpHit ) {
		memcpy(buffer,hhitCache,sizeof(HIT)*NBHHIT);
		*nb = nbHHit;
	}

}

// -------------------------------------------------------------

void Worker::SetHHit(HIT *buffer,int *nb,SHGHITS *gHits) {

	if( nb>0 ) {
		memcpy(hhitCache,buffer,sizeof(HIT)*MIN(*nb,NBHHIT));
		memcpy(gHits->pHits,buffer,sizeof(LEAK)*MIN(NBHLEAK,*nb));
		gHits->nbHHit = *nb;
	}

}

// -------------------------------------------------------------

void Worker::InnerStop(float appTime) {

	stopTime =appTime;
	simuTime+=appTime-startTime;
	running  = FALSE;
	calcAC = FALSE;

}





void Worker::Stop_Public() {
	// Stop
	InnerStop(m_fTime);
	try {
		Stop();
		Update(m_fTime);
	} catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error (Stop)",GLDLG_OK,GLDLG_ICONERROR);
	}
}

// -------------------------------------------------------------

void Worker::StartStop(float appTime,int mode) {

	if( running )  {

		// Stop
		InnerStop(appTime);
		try {
			Stop();
			Update(appTime);
		} catch(Error &e) {
			GLMessageBox::Display((char *)e.GetMsg(),"Error (Stop)",GLDLG_OK,GLDLG_ICONERROR);
		}

	} else {

		// Start
		if (needsReload) RealReload();
		startTime = appTime;
		running = TRUE;
		calcAC = FALSE;
		this->mode = mode;
		try {
			Start();
		} catch(Error &e) {
			running = FALSE;
			GLMessageBox::Display((char *)e.GetMsg(),"Error (Start)",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}

		// Particular case when simulation ends before getting RUN state
		if(allDone) {
			Update(appTime);
			GLMessageBox::Display("Max desorption reached","Information (Start)",GLDLG_OK,GLDLG_ICONINFO);
		}

	}

}

// -------------------------------------------------------------

void Worker::Update(float appTime) {
	if (needsReload) RealReload();
	//if (!needsReload) {
		// Check calculation ending
		BOOL done = TRUE;
		BOOL error = TRUE;
		if( dpControl ) {
			if( AccessDataport(dpControl) ) {
				int i = 0;
				SHMASTER *master =(SHMASTER *)dpControl->buff;
				for(int i=0;i<nbProcess && done;i++) {
					done = done && (master->states[i]==PROCESS_DONE);
					error = error && (master->states[i]==PROCESS_ERROR);
				}
				ReleaseDataport(dpControl);
			}
		}

		// End of simulation reached (Stop GUI)
		if( (error || done) && running && appTime!=0.0f )
			InnerStop(appTime);

		// Retrieve hit count recording from the shared memory
		if( dpHit ) {

			if( AccessDataport(dpHit) ) {
				BYTE *buffer = (BYTE *)dpHit->buff;

				changedSinceSave = TRUE;
				// Globals
				SHGHITS *gHits = (SHGHITS *)buffer;

				// Global hits and leaks
				nbHit = gHits->total.nbHit;
				nbAbsorption = gHits->total.nbAbsorbed;
				nbDesorption = gHits->total.nbDesorbed;
				totalFlux = gHits->total.fluxAbs;
				totalPower = gHits->total.powerAbs;
				nbLeak = gHits->nbLeak;
				nbHHit = gHits->nbHHit;
				memcpy(hhitCache,gHits->pHits,sizeof(HIT)*NBHHIT);
				memcpy(leakCache,gHits->pLeak,sizeof(LEAK)*NBHHIT);

				// Facets
				int nbFacet = geom->GetNbFacet();
				for(int i=0;i<nbFacet;i++) {    
					Facet *f = geom->GetFacet(i);
					memcpy(&(f->sh.counter),buffer+f->sh.hitOffset,sizeof(SHHITS));
				}
				geom->BuildTexture(buffer);
				ReleaseDataport(dpHit);
			}
		}
	//}
}

// -------------------------------------------------------------

void Worker::SendHits() {
	//if (!needsReload) {
	if(dpHit) {
		if( AccessDataport(dpHit) ) {

			// Store initial hit counts in the shared memory
			BYTE *pBuff = (BYTE *)dpHit->buff;
			memset(pBuff,0,geom->GetHitsSize());

			SHGHITS *gHits = (SHGHITS *)pBuff;
			gHits->total.nbHit = nbHit;
			gHits->total.nbDesorbed = nbDesorption;
			gHits->total.nbAbsorbed = nbAbsorption;
			gHits->total.fluxAbs = totalFlux;
			gHits->total.powerAbs = totalPower;

			int nbFacet = geom->GetNbFacet();
			for(int i=0;i<nbFacet;i++) {
				Facet *f = geom->GetFacet(i);
				memcpy(pBuff+f->sh.hitOffset,&(f->sh.counter),sizeof(SHHITS));
			}
			ReleaseDataport(dpHit);

		} else {
			throw Error("Failed to initialise 'hits' dataport");
		}
	}
	//}
}

// -------------------------------------------------------------

void  Worker::ReleaseHits() {

	ReleaseDataport(dpHit);

}

// -------------------------------------------------------------

BYTE *Worker::GetHits() {

	if( dpHit )
		if( AccessDataport(dpHit) ) 
			return (BYTE *)dpHit->buff;

	return NULL;

}

// -------------------------------------------------------------

void Worker::ThrowSubProcError(char *message) {

	char errMsg[1024];
	if( !message )
		sprintf(errMsg,"Bad response from sub process(es):\n%s",GetErrorDetails());
	else
		sprintf(errMsg,"%s\n%s",message,GetErrorDetails());
	throw Error(errMsg);

}

void Worker::Reload (){
	needsReload=true;
}

// -------------------------------------------------------------

void Worker::RealReload() { //Sharing geometry with workers

	if(nbProcess==0) return;

	GLProgress *progressDlg = new GLProgress("Asking subprocesses to clear geometry...","Passing Geometry to workers");
	progressDlg->SetVisible(TRUE);
	progressDlg->SetProgress(0.0);

	// Clear geometry
	CLOSEDP(dpHit);
	if( !ExecuteAndWait(COMMAND_CLOSE,PROCESS_READY) )
		ThrowSubProcError();

	if(!geom->IsLoaded()) return;

	progressDlg->SetMessage("Creating dataport...");
	// Create the temporary geometry shared structure
	int loadSize = geom->GetGeometrySize(&regions,&materials);
	Dataport *loader = CreateDataport(loadDpName,loadSize);
	if( !loader )
		throw Error("Failed to create 'loader' dataport");
	progressDlg->SetMessage("Accessing dataport...");
	AccessDataportTimed(loader,3000+nbProcess*(int)((double)loadSize/10000.0));
	progressDlg->SetMessage("Assembling geometry and regions to pass...");
	geom->CopyGeometryBuffer((BYTE *)loader->buff,&regions,&materials);
	progressDlg->SetMessage("Releasing dataport...");
	ReleaseDataport(loader);

	int hitSize = geom->GetHitsSize();
	dpHit = CreateDataport(hitsDpName,hitSize);
	if( !dpHit ) {
		CLOSEDP(loader);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("Failed to create 'hits' dataport");
	}

	// Compute number of max desorption per process
	if(AccessDataportTimed(dpControl,3000+nbProcess*(int)((double)loadSize/10000.0))) {
		SHMASTER *m = (SHMASTER *)dpControl->buff;
		llong common = maxDesorption / (llong)nbProcess;
		int remain = (int)(maxDesorption % (llong)nbProcess);
		for(int i=0;i<nbProcess;i++) {
			m->cmdParam2[i] = common;
			if(i<remain) m->cmdParam2[i]++;
		}
		ReleaseDataport(dpControl);
	}

	// Load geometry
	progressDlg->SetMessage("Waiting for subprocesses to load geometry...");
	if( !ExecuteAndWait(COMMAND_LOAD,PROCESS_READY,loadSize,progressDlg) ) {
		CLOSEDP(loader);
		char errMsg[1024];
		sprintf(errMsg,"Failed to send geometry to sub process:\n%s",GetErrorDetails());
		GLMessageBox::Display(errMsg,"Warning (Load)",GLDLG_OK,GLDLG_ICONWARNING);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		return;
	}
	progressDlg->SetProgress(1.0);
	progressDlg->SetMessage("Sending hits...");
	//Send hit counts
	try {
		SendHits();
	} catch(Error &e) {
		// Geometry not loaded !
		CLOSEDP(dpHit);
		CLOSEDP(loader);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw e;
	}
	progressDlg->SetMessage("Closing dataport...");
	CLOSEDP(loader);
	needsReload=false;
	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
}

// -------------------------------------------------------------

void Worker::SetMaxDesorption(llong max) {

	try {
		Reset(0.0);
		maxDesorption = max;
		Reload();
	} catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
	}

}

// -------------------------------------------------------------

char *Worker::GetErrorDetails() {

	static char err[1024];
	strcpy(err,"");

	AccessDataport(dpControl);
	SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
	for(int i=0;i<nbProcess;i++) {
		char tmp[256];
		if( pID[i]!=0 ) {
			int st = shMaster->states[i];
			if(st==PROCESS_ERROR) {
				sprintf(tmp,"[#%d] Process [PID %d] %s: %s\n",i,pID[i],prStates[st],shMaster->statusStr[i]);
			} else {
				sprintf(tmp,"[#%d] Process [PID %d] %s\n",i,pID[i],prStates[st]);
			}
		} else {
			sprintf(tmp,"[#%d] Process [PID ???] Not started\n",i);
		}
		strcat(err,tmp);
	}
	ReleaseDataport(dpControl);

	return err;
}

// -------------------------------------------------------------

void Worker::ClearHits() {
	if (needsReload) RealReload();
	if(dpHit) {
		AccessDataport(dpHit);
		memset(dpHit->buff,0,geom->GetHitsSize());
		ReleaseDataport(dpHit);
	}

}

// -------------------------------------------------------------

BOOL Worker::Wait(int waitState,int timeout,GLProgress *prg) {

	BOOL ok = FALSE;
	BOOL error = FALSE;
	int t = 0;
	int nbReady = 0;
	double initialProgress = 0.0;
	if (prg) initialProgress = prg->GetProgress();
	int nbError = 0;
	allDone = TRUE;

	// Wait for completion
	while(!ok && t<timeout) {

		ok = TRUE;
		AccessDataport(dpControl);
		SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
		nbReady=nbError=0;
		for(int i=0;i<nbProcess;i++) {
			if (shMaster->states[i]==waitState) nbReady++;
			ok = ok & (shMaster->states[i]==waitState || shMaster->states[i]==PROCESS_ERROR || shMaster->states[i]==PROCESS_DONE);
			if( shMaster->states[i]==PROCESS_ERROR ) {
				error = TRUE;
				nbError++;
			}
			allDone = allDone & (shMaster->states[i]==PROCESS_DONE);
		}
		ReleaseDataport(dpControl);

		if(!ok) {
			if (prg) prg->SetProgress(double(nbReady)/(double)nbProcess);
			Sleep(500);
			t+=500;
		}

	}

	if (t>=timeout) {
		if ((prg) && ((double)nbReady/(double)nbProcess)>initialProgress) //progress advanced, wait more
			return Wait(waitState,timeout,prg);
		char tmp[512];
		sprintf(tmp,"Total workers : %d\n"
			"%d are ready, %d reported errors\n"
			"Do you want to wait a bit more?\n"
			"(Loading continues while this dialog is visible)\n",nbProcess,nbReady,nbError,nbProcess);
		int waitmore = GLMessageBox::Display(tmp,"Info",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK;
		if (waitmore) {
			t=0;
			return Wait(waitState,timeout,prg);
		}
	}
	return ok && !error;

}

// -------------------------------------------------------------

BOOL Worker::ExecuteAndWait(int command,int waitState,int param,GLProgress *prg) {

	if(!dpControl) return FALSE;

	// Send command
	AccessDataport(dpControl);
	SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
	for(int i=0;i<nbProcess;i++) {
		shMaster->states[i]=command;
		shMaster->cmdParam[i]=param;
	}
	ReleaseDataport(dpControl);

	Sleep(100);
	return Wait(waitState,8000,prg);
}

void Worker::ResetWorkerStats() {

	nbAbsorption = 0;
	nbDesorption = 0;
	nbHit = 0;
	nbLeak = 0;
	totalFlux = 0.0;
	totalPower = 0.0;

}

// -------------------------------------------------------------

void Worker::Reset(float appTime) {

	if( calcAC ) {
		GLMessageBox::Display("Reset not allowed while calculating AC","Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	stopTime = 0.0f;
	startTime = 0.0f;
	simuTime = 0.0f;
	running = FALSE;
	if( nbProcess==0 ) 
		return;

	try {
		ResetWorkerStats();
		if( !ExecuteAndWait(COMMAND_RESET,PROCESS_READY) )
			ThrowSubProcError();
		ClearHits();
		Update(appTime);
	} catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
	}

}

// -------------------------------------------------------------

void Worker::Start() {
	/*
	// Check that at least one desortion facet exists
	BOOL found = FALSE;
	int nbF = geom->GetNbFacet();
	int i=0;
	while(i<nbF && !found) {
	found = (geom->GetFacet(i)->sh.desorbType!=DES_NONE);
	if(!found) i++;
	}

	if( !found )
	throw Error("No desorption facet found");

	if(!(totalOutgassing>0.0))
	throw Error("Total outgassing is zero.");
	*/
	if( nbProcess==0 )
		throw Error("No sub process found. (Simulation not available)");

	if( !ExecuteAndWait(COMMAND_START,PROCESS_RUN,mode) )
		ThrowSubProcError();

}

// -------------------------------------------------------------

void Worker::Stop() {

	if( nbProcess==0 )
		throw Error("No sub process found. (Simulation not available)");

	if( !ExecuteAndWait(COMMAND_PAUSE,PROCESS_READY) )
		ThrowSubProcError();

}

// -------------------------------------------------------------

void Worker::KillAll() {

	if( dpControl && nbProcess>0 ) {
		if( !ExecuteAndWait(COMMAND_EXIT,PROCESS_KILLED) ) {
			AccessDataport(dpControl);
			SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
			for(int i=0;i<nbProcess;i++) 
				if(shMaster->states[i]==PROCESS_KILLED) pID[i]=0;
			ReleaseDataport(dpControl);
			// Force kill
			for(int i=0;i<nbProcess;i++)
				if(pID[i]) KillProc(pID[i]);
		}
		CLOSEDP(dpHit);
	}
	nbProcess = 0;

}

// -------------------------------------------------------------

void Worker::Exit() {

	if( dpControl && nbProcess>0 ) {
		KillAll();
		CLOSEDP(dpControl);
	}

}

// -------------------------------------------------------------

void Worker::GetProcStatus(int *states,char **status) {

	if(nbProcess==0) return;

	AccessDataport(dpControl);
	SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
	memcpy(states,shMaster->states,MAX_PROCESS*sizeof(int));
	memcpy(status,shMaster->statusStr,MAX_PROCESS*64);
	ReleaseDataport(dpControl);

}

void Worker::SetProcNumber(int n) {

	char cmdLine[512];

	// Kill all sub process
	KillAll();

	// Create new control dataport
	if( !dpControl ) 
		dpControl = CreateDataport(ctrlDpName,sizeof(SHMASTER));
	if( !dpControl )
		throw Error("Failed to create 'control' dataport");
	AccessDataport(dpControl);
	memset(dpControl->buff,0,sizeof(SHMASTER));
	ReleaseDataport(dpControl);

	// Launch n subprocess
	for(int i=0;i<n;i++) {
		sprintf(cmdLine,"synradSub.exe %d %d",pid,i);
		pID[i] = StartProc(cmdLine);
		Sleep(25); // Wait a bit
		if( pID[i]==0 ) {
			nbProcess = 0;
			throw Error(cmdLine);
		}
	}

	nbProcess = n;

	if( !Wait(PROCESS_READY,3000) )
		ThrowSubProcError("Sub process(es) starting failure");

}

int Worker::GetProcNumber() {
	return nbProcess;
}

DWORD Worker::GetPID(int prIdx) {
	return pID[prIdx];
}

void Worker::AddRegion(char *fileName,int position) {
	//if (!geom->IsLoaded()) throw Error("Load geometry first!");
	needsReload=TRUE;
	SynRad *mApp = (SynRad *)theApp;
	char *ext;

	ext = strrchr(fileName,'.');

	if( ext==NULL )
		throw Error("LoadParam(): Invalid file extension [Only PAR or param]");
	ext++;

	// Read a file
	FileReader *f = NULL;

	GLProgress *progressDlg = new GLProgress("Reading file...","Please wait");
	progressDlg->SetVisible(TRUE);
	progressDlg->SetProgress(0.0);

	BOOL isPAR=(_stricmp(ext,"par")==0);
	BOOL isParam=(_stricmp(ext,"param")==0);

	if(isPAR || isParam) {

		try {
			f = new FileReader(fileName);
			Region newtraj;
			if (isPAR) newtraj.LoadPAR(f,progressDlg);
			else newtraj.LoadParam(f,progressDlg);
			newtraj.fileName=fileName;
			if (position==-1) regions.push_back(newtraj);
			else {
				nbTrajPoints-=(int)regions[position].Points.size();
				regions[position]=newtraj;
				regions[position].Points=newtraj.Points; //need to force because of operator=
			}
			geom->InitializeGeometry(-1,TRUE); //recalculate bounding box
			nbTrajPoints+=(int)newtraj.Points.size();
		} catch(Error &e) {
			//geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	}  else {
		SAFE_DELETE(f);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("LoadParam(): Invalid file extension [Only par or param]");
	}
	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
	SAFE_DELETE(f);
}

void Worker::RecalcRegion(int regionId) {
	
	needsReload=TRUE;
	SynRad *mApp = (SynRad *)theApp;
	try {
			
			Region newtraj;
			newtraj=regions[regionId]; //copy all except the points
			newtraj.Points=newtraj.CalculateTrajectory(newtraj.dL,newtraj.E,newtraj.limits,
				newtraj.startPoint,newtraj.startDir,1000000); //points calculated
				nbTrajPoints-=(int)regions[regionId].Points.size();
				regions[regionId]=newtraj;
				regions[regionId].Points=newtraj.Points; //need to force because of operator=
			geom->InitializeGeometry(-1,TRUE); //recalculate bounding box
			nbTrajPoints+=(int)newtraj.Points.size();
		} catch(Error &e) {
			throw e;
		}

}

void Worker::ClearRegions() {
	//if ((int)regions.size()>0) regions.clear();
	regions=std::vector<Region>();
	geom->InitializeGeometry(); //recalculate bounding box
	nbTrajPoints=0;
	Reload();
}

void Worker::AddMaterial(string *fileName){
	Material result;
	char tmp[512];
	strcpy(tmp,fileName->c_str());
	FileReader *f=new FileReader(tmp);
	result.LoadCSV(f);
	materials.push_back(result);
}


void Worker::SaveRegion(char *fileName,int position,BOOL overwrite) {
	char tmp[512];

	// Read a file
	FILE *f = NULL;
	BOOL ok = TRUE;
	
	if( FileUtils::Exist(fileName) && !overwrite) {
		sprintf(tmp,"Overwrite existing file ?\n%s",fileName);
		ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
	}
	
	if( ok || overwrite ) {
		FileWriter *file=new FileWriter(fileName);
		regions[position].SaveParam(file);
		SAFE_DELETE(file);
	}
}

BOOL EndsWithPar(const char* s)
{
  int ret = 0;

  if (s != NULL)
  {
    size_t size = strlen(s);

    if (size >= 4 &&
        s[size-4] == '.' &&
        s[size-3] == 'p' &&
        s[size-2] == 'a' &&
        s[size-1] == 'r')
    {
      ret = 1;
    }
  }

  return ret;
}