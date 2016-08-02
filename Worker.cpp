/*
File:        Worker.cpp
Description: Sub processes handling
Program:     SynRad
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


#include "Worker.h"
#include "GLApp/GLApp.h"
#include "GLApp/GLMessageBox.h"
#include <math.h>
#include <stdlib.h>
#include <Process.h>
#include "GLApp/GLUnitDialog.h"
#include "Synrad.h"
#include <direct.h>

#include "ZipUtils/zip.h"
#include "ZipUtils/unzip.h"

using namespace pugi;
extern SynRad *mApp;

BOOL EndsWithPar(const char* s);
Worker::Worker() {

	pid = _getpid();
	sprintf(ctrlDpName,"SRDCTRL%d",pid);
	sprintf(loadDpName,"SRDLOAD%d",pid);
	sprintf(hitsDpName,"SRDHITS%d",pid);
	sprintf(materialsDpName,"SRDMATS%d",pid);
	nbProcess = 0;
	maxDesorption = 0;
	distTraveledTotal=0.0;
	lowFluxCutoff = 1E-7;
	lowFluxMode = FALSE;
	ResetWorkerStats();
	geom = new Geometry();
	regions = std::vector<Region_full>();
	generation_mode=SYNGEN_MODE_POWERWISE;
	dpControl = NULL;
	dpHit = NULL;
	nbHHit = 0;
	nbHit = 0;
	nbLeakTotal = 0;
	nbLastLeaks = 0;
	startTime = 0.0f;
	stopTime = 0.0f;
	simuTime = 0.0f;
	nbTrajPoints=0;
	running = FALSE;

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

BOOL Worker::IsDpInitialized(){

	return (dpHit != NULL);
}
// -------------------------------------------------------------

char *Worker::GetFileName() {
	return fullFileName;
}

char *Worker::GetShortFileName() {

	static char ret[512];
	char *r = strrchr(fullFileName,'/');
	if(!r) r = strrchr(fullFileName,'\\');
	if(!r) strcpy(ret,fullFileName);
	else   {
		r++;
		strcpy(ret,r);
	}

	return ret;

}

char *Worker::GetShortFileName(char* longFileName) {

	static char ret[512];
	char *r = strrchr(longFileName, '/');
	if (!r) r = strrchr(longFileName, '\\');
	if (!r) strcpy(ret, longFileName);
	else   {
		r++;
		strcpy(ret, r);
	}

	return ret;

}

// -------------------------------------------------------------

void Worker::SetFileName(char *fileName) {

	strcpy(fullFileName,fileName);
}








void Worker::SaveGeometry(char *fileName,GLProgress *prg,BOOL askConfirm,BOOL saveSelected,BOOL autoSave,BOOL crashSave) {

	try {
		if (needsReload&&(!crashSave && !saveSelected)) RealReload();


	} catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"Error reloading worker. Trying crash save:\n%s",e.GetMsg());
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		crashSave=TRUE;
	} 
	char tmp[65536]; //compress.exe command line
	/*char fileNameWithGeo[2048]; //file name with .geo extension (instead of .geo7z)
	char fileNameWithGeo7z[2048];*/
	char fileNameWithSyn[2048]; //file name with .syn extension (instead of .syn7z)
	char fileNameWithSyn7z[2048];
	char fileNameWithoutExtension[2048]; //file name without extension
	//char *ext = fileName+strlen(fileName)-4;
	char *ext,*dir;

	dir = strrchr(fileName,'\\');
	ext = strrchr(fileName,'.');

	if(!(ext) || !(*ext=='.') || ((dir)&&(dir>ext)) ) { 
		sprintf(fileName, mApp->compressSavedFiles ? "%s.syn7z" : "%s.syn", fileName); //set to default SYN/SYN7Z format
		ext = strrchr(fileName,'.');
	}

	ext++;

	// Read a file
	BOOL ok = TRUE;
	FileWriter *f = NULL;
	BOOL isTXT = _stricmp(ext,"txt")==0;
	BOOL isSTR = _stricmp(ext,"str")==0;
	/*BOOL isGEO = _stricmp(ext,"geo")==0;
	BOOL isGEO7Z = _stricmp(ext,"geo7z")==0;*/
	BOOL isSYN = _stricmp(ext,"syn")==0;
	BOOL isSYN7Z = _stricmp(ext,"syn7z")==0;


	if(isTXT || isSYN || isSYN7Z || isSTR) {

		if (WAIT_TIMEOUT==WaitForSingleObject(mApp->compressProcessHandle,0)) {
			GLMessageBox::Display("Compressing a previous save file is in progress. Wait until that finishes"
				"or close process \"compress.exe\"\nIf this was an autosave attempt,"
				"you have to lower the autosave frequency.","Can't save right now.",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
		if (isSYN) {
			memcpy(fileNameWithoutExtension,fileName,sizeof(char)*(strlen(fileName)-4));
			fileNameWithoutExtension[strlen(fileName)-4]='\0';
			sprintf(fileNameWithSyn7z,"%s7z",fileName);
			memcpy(fileNameWithSyn,fileName,(strlen(fileName)+1)*sizeof(char));

		} else if (isSYN7Z) {
			memcpy(fileNameWithoutExtension,fileName,sizeof(char)*(strlen(fileName)-6));
			fileNameWithoutExtension[strlen(fileName)-6]='\0';
			memcpy(fileNameWithSyn,fileName,sizeof(char)*(strlen(fileName)-2));
			fileNameWithSyn[strlen(fileName)-2]='\0';
			memcpy(fileNameWithSyn7z,fileName,(1+strlen(fileName))*sizeof(char));
			sprintf(tmp,"A .syn file of the same name exists. Overwrite that file ?\n%s",fileNameWithSyn);
			if(!autoSave && FileUtils::Exist(fileNameWithSyn) ) {

				ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
			}
		}

		if(!autoSave && ok && FileUtils::Exist(fileName) ) {
			sprintf(tmp,"Overwrite existing file ?\n%s",fileName);
			if (askConfirm) ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
		}

		if( ok ) {
			if( isSTR ) {
				geom->SaveSTR(dpHit,saveSelected);


			} else {
				try {

					if (isSYN7Z) {
						/*memcpy(fileNameWithSyn,fileName,sizeof(char)*(strlen(fileName)-2));
						fileNameWithSyn[strlen(fileName)-2]='\0';*/
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
				else if( isSYN || isSYN7Z ) {
					// Retrieve leak cache
					int nbLeakSave,nbHHitSave;
					LEAK pLeak[NBHLEAK];
					if (!crashSave && !saveSelected) GetLeak(pLeak,&nbLeakSave);
					else nbLeakSave=0;
					// Retrieve hit cache (lines and dots)
					HIT pHits[NBHHIT];
					if (!crashSave && !saveSelected) GetHHit(pHits,&nbHHitSave);

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
			if (!autoSave && !saveSelected) {
				strcpy(fullFileName,fileName);
				remove("Synrad_AutoSave.syn");

				remove("Synrad_AutoSave.syn7z");
			}
		}
	} else {
		SAFE_DELETE(f);
		throw Error("SaveGeometry(): Invalid file extension [only syn,txt,str]");
	}

	SAFE_DELETE(f);
	/*if (ok && isGEO || isGEO7Z) {



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
	} else*/ if (ok && isSYN7Z) {
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
				mApp->compressProcessHandle=OpenProcess(PROCESS_ALL_ACCESS, TRUE, procId);
				fileName=fileNameWithSyn7z;
			} else {
				GLMessageBox::Display("compress.exe (part of Synrad) not found.\n Will save as uncompressed SYN file.","Compressor not found",GLDLG_OK,GLDLG_ICONERROR);
				fileName=fileNameWithSyn;
			}
		} else if (ok && isSYN) fileName=fileNameWithSyn;
		if (!autoSave && !saveSelected) {
			SetFileName(fileName);
			mApp->UpdateTitle();
		}
	



}



void Worker::ExportTextures(char *fileName,int grouping,int mode,BOOL askConfirm,BOOL saveSelected) {

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
		if (!f) {



			throw Error("Couldn't open file for writing");
		}
		geom->ExportTextures(f,grouping,mode,no_scans,dpHit,saveSelected);
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
	try {
	if (needsReload) RealReload();
	geom->SaveDesorption(f,dpHit,selectedOnly,mode,eta0,alpha,distr);




	} catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"Save error:\n%s",e.GetMsg());
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
	}
	fclose(f);

}

void Worker::LoadGeometry(char *fileName, BOOL insert, BOOL newStr) {
	if (!insert) {
		needsReload=TRUE;
	}
	else { //insert
		if (needsReload) RealReload();
	}
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );

	std::string ext = FileUtils::GetExtension(fileName);

	if( ext=="" )
		throw Error("LoadGeometry(): No file extension, can't determine type");

	// Read a file
	FileReader *f = NULL;
	GLProgress *progressDlg = new GLProgress("Reading file...","Please wait");
	progressDlg->SetVisible(TRUE);
	progressDlg->SetProgress(0.0);

	if (!insert) {
		//Clear hits and leaks cache
		memset(hhitCache, 0, sizeof(HIT)*NBHHIT);
		memset(leakCache, 0, sizeof(LEAK)*NBHLEAK);
		
	}

	if(ext=="txt" || ext=="TXT") {

		try {
			if (!insert) ResetWorkerStats();
			f = new FileReader(fileName);
			if (!insert) {
				geom->LoadTXT(f, progressDlg);
				strcpy(fullFileName, fileName);
			}
			else { //insert
				mApp->changedSinceSave = TRUE;
				geom->InsertTXT(f, progressDlg, newStr);
				Reload();
			}

		} catch(Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}


	} else if(ext=="stl") {
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
				progressDlg->SetVisible(TRUE);
				ResetWorkerStats();				
				progressDlg->SetMessage("Reading geometry...");
				f = new FileReader(fileName);
				if (!insert) {					
					geom->LoadSTL(f, progressDlg, scaleFactor);
					strcpy(fullFileName, fileName);
					mApp->DisplayCollapseDialog();
				}
				else { //insert
					mApp->changedSinceSave = TRUE;
					geom->InsertSTL(f, progressDlg, scaleFactor, newStr);
					Reload();
				}
			}
		} catch(Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;

		}

	} else if(ext=="str") {
		if (insert) throw Error("STR file inserting is not supported.");
		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadSTR(f,progressDlg);

			strcpy(fullFileName,fileName);
		} 
		catch(Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else if (ext == "syn" || ext == "syn7z") { //Synrad file
		int version;
		progressDlg->SetVisible(TRUE);
		try {
			if (ext == "syn7z") {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				system(tmp);
				f = new FileReader((std::string)CWD + "\\tmp\\Geometry.syn"); //Open extracted file
			} else f = new FileReader(fileName); //syn file, open it directly
			std::vector<std::string> regionsToLoad;
			LEAK pLeak[NBHLEAK];
			HIT pHits[NBHHIT];
			if (!insert) {
				progressDlg->SetMessage("Resetting worker...");
				ResetWorkerStats();

				regionsToLoad = geom->LoadSYN(f, progressDlg, pLeak, &nbLastLeaks, pHits, &nbHHit, &version);
				//copy temp values from geom to worker:
				nbLeakTotal = geom->tNbLeak;
				nbHit = geom->tNbHit;
				nbDesorption = geom->tNbDesorption;
				nbAbsorption = geom->tNbAbsorption;
				distTraveledTotal = geom->distTraveledTotal;
				maxDesorption = geom->tNbDesorptionMax;
				totalFlux = geom->tFlux;
				totalPower = geom->tPower;
				no_scans = geom->loaded_no_scans;
			}
			else { //insert
				regionsToLoad = geom->InsertSYN(f, progressDlg, newStr);
			}
			
			//Load regions
			if (regionsToLoad.size()>0) {
				char tmp[256];
				sprintf(tmp,"This geometry refers to %d regions. Load them now?",regionsToLoad.size());
				BOOL loadThem = ( GLMessageBox::Display(tmp,"File load",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK );
				if (loadThem) {
					progressDlg->SetMessage("Loading regions");
					int i = 0;
					for (auto&& regionFileName:regionsToLoad) {
						std::string toLoad;
						if (ext=="syn7z") { //PAR file to load in tmp dir (just extracted)
							toLoad=(std::string)CWD+"\\tmp\\";
							toLoad += FileUtils::GetFilename(regionFileName);
						} else { //PAR file in same dir as SYN file
							/*char tmp[512];
							filebegin= strrchr(fileName,'\\'); //throw out absolute path, keep only filename
							if (filebegin) filebegin++;
							else filebegin=fileName;
							memcpy(tmp,fileName,filebegin-fileName);
							tmp[(int)(filebegin-fileName)]=NULL;
							toLoad=tmp;*/
							toLoad = FileUtils::GetPath(fileName);
							toLoad+=regionFileName;
						}
						progressDlg->SetMessage("Adding "+regionFileName+"...");
						progressDlg->SetProgress((double)i++ / (double)regionsToLoad.size());
						AddRegion(toLoad.c_str(), -1);
					}
				}
			}
			if (!insert) {
				progressDlg->SetMessage("Reloading worker with new geometry...");
				RealReload(); //for the loading of textures
				geom->LoadProfileSYN(f, dpHit);
				geom->LoadSpectrumSYN(f, dpHit);
				SHGHITS *gHits = (SHGHITS *)dpHit->buff;
				SetLeak(pLeak, &nbLastLeaks, gHits);
				SetHHit(pHits, &nbHHit, gHits);
				progressDlg->SetMessage("Loading texture values...");
				LoadTexturesSYN((ext=="syn7z") ? ((std::string)CWD+"\\tmp\\Geometry.syn").c_str() : fileName, version);
				strcpy(fullFileName, fileName);
			}
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
		SAFE_DELETE(f);
	} else if(ext=="geo" || ext=="geo7z") {
		std::string toOpen;
		int version;
		progressDlg->SetVisible(TRUE);
		try {
			if (ext=="geo7z") {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				//execCMD(tmp);
				system(tmp);

				/*filebegin= strrchr(fileName,'\\');
				if (filebegin) filebegin++;
				else filebegin=fileName;
				memcpy(fileOnly,filebegin,sizeof(char)*(strlen(filebegin)-2));
				fileOnly[strlen(filebegin)-2]='\0';
				sprintf(tmp2,"%s\\tmp\\Geometry.geo",CWD);*/
				toOpen = (std::string)CWD + "\\tmp\\Geometry.geo"; //newer geo7z format: contain Geometry.geo
				if (!FileUtils::Exist(toOpen)) toOpen = ((std::string)fileName).substr(0, strlen(fileName) - 2); //Inside the zip, try original filename with extension changed from geo7z to geo
				f = new FileReader(toOpen);
			}
			else { //not geo7z
				toOpen = fileName;
				f = new FileReader(fileName); //geo file, open it directly
			}
			progressDlg->SetMessage("Resetting worker...");
			ResetWorkerStats();
			if (!insert) {
				//leaks
				LEAK pLeak[NBHLEAK];
				//hits
				HIT pHits[NBHHIT];

				geom->LoadGEO(f, progressDlg, pLeak, &nbLastLeaks, pHits, &nbHHit, &version);
				SAFE_DELETE(f);
				maxDesorption = 0;
				progressDlg->SetMessage("Reloading worker with new geometry...");
				RealReload();
				strcpy(fullFileName,fileName);
			}
			else { //insert
				mApp->changedSinceSave = TRUE;
				geom->InsertGEO(f,progressDlg,newStr);
				Reload();
			}
			
		} catch(Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			//if (isGEO7Z) remove(tmp2);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else if (ext == "xml" || ext == "zip") { //XML file, optionally in ZIP container
		xml_document loadXML;
		xml_parse_result parseResult;
		progressDlg->SetVisible(TRUE);
		try {
			if (ext == "zip") { //compressed in ZIP container
				//decompress file
				progressDlg->SetMessage("Decompressing file...");

				HZIP hz = OpenZip(fileName, 0);
				if (!hz) {
					throw Error("Can't open ZIP file");
				}
				ZIPENTRY ze; GetZipItem(hz, -1, &ze); int numitems = ze.index;
				BOOL notFoundYet = TRUE;
				for (int i = 0; i < numitems && notFoundYet; i++) { //extract first xml file found in ZIP archive
					GetZipItem(hz, i, &ze);
					std::string zipFileName = ze.name;

					if (FileUtils::GetExtension(zipFileName) == "xml") { //if it's an .xml file
						notFoundYet = FALSE;
						std::string tmpFileName = "tmp/" + zipFileName;
						UnzipItem(hz, i, tmpFileName.c_str()); //unzip it to tmp directory
						CloseZip(hz);
						progressDlg->SetMessage("Reading and parsing XML file...");
						parseResult = loadXML.load_file(tmpFileName.c_str()); //load and parse it
					}
				}
				if (notFoundYet) {
					CloseZip(hz);
					throw Error("Didn't find any XML file in the ZIP file.");
				}
			}
			else parseResult = loadXML.load_file(fileName); //parse xml file directly

			ResetWorkerStats();
			if (!parseResult) {
				//Parse error
				std::stringstream err;
				err << "XML parsed with errors.\n";
				err << "Error description: " << parseResult.description() << "\n";
				err << "Error offset: " << parseResult.offset << "\n";
				throw Error(err.str().c_str());
			}

			progressDlg->SetMessage("Building geometry...");
			if (!insert) {
				geom->LoadXML_geom(loadXML, this, progressDlg);
				maxDesorption = 0;
				progressDlg->SetMessage("Reloading worker with new geometry...");
				RealReload();
				geom->UpdateName(fileName);
			}
			else { //insert
				mApp->changedSinceSave = TRUE;
				geom->InsertXML(loadXML, this, progressDlg, newStr);
				Reload();
			}
		}
		catch (Error &e) {
			if (!insert) geom->Clear();
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else if (ext == "ase") {
		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadASE(f,progressDlg);
			//RealReload();
			strcpy(fullFileName,fileName);

		} catch(Error &e) {
			if (!insert) geom->Clear();
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
	
	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
	SAFE_DELETE(f);
}

/*void Worker::InsertGeometry(BOOL newStr,char *fileName) {
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
	BOOL isXML = _stricmp(ext, "xml") == 0;
	BOOL isXMLzip = _stricmp(ext, "zip") == 0;

	if(_stricmp(ext,"txt")==0) {

		try {
			f=new FileReader(fileName);
			mApp->changedSinceSave = TRUE;
			geom->InsertTXT(f,progressDlg,newStr);
			Reload();
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
				mApp->changedSinceSave = TRUE;
				ResetWorkerStats();
				f = new FileReader(fileName);
				geom->InsertSTL(f,progressDlg,scaleFactor,newStr);
				Reload();}
		} catch(Error &e) {
			//geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"str")==0) {
		throw Error("STR file inserting is not supported.");

	} else if(_stricmp(ext,"syn")==0 || isSYN7Z) {

		char tmp2[1024];
		
		progressDlg->SetVisible(TRUE);
		try {
			if (isSYN7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
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


			std::vector<std::string> regionsToLoad = geom->InsertSYN(f, progressDlg, newStr);

			//nbLeak = 0;
			//nbHit = 0;
			//nbDesorption = 0;
			//maxDesorption = geom->tNbDesorptionMax;

			//progressDlg->SetMessage("Reloading worker with new geometry...");

			//RealReload(); //for the loading of textures
			//geom->LoadProfileSYN(f,dpHit);
			//geom->LoadSpectrumSYN(f,dpHit);
			//SHGHITS *gHits = (SHGHITS *)dpHit->buff;
			//SetLeak(pLeak,&nbLeak,gHits);
			//SetHHit(pHits,&nbHHit,gHits);
			//SAFE_DELETE(f);
			//progressDlg->SetMessage("Loading textures...");
			//loadTextures((isSYN7Z)?tmp2:fileName,version);
			//strcpy(fullFileName,fileName);
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
			char tmp2[1024];
			if (isGEO7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				
				char fileOnly[512];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				//execCMD(tmp);
				system(tmp);

				/*char *filebegin= strrchr(fileName,'\\');
				if (filebegin) filebegin++;
				else filebegin=fileName;
				memcpy(fileOnly,filebegin,sizeof(char)*(strlen(filebegin)-2));
				fileOnly[strlen(filebegin)-2]='\0';
				sprintf(tmp2,"%s\\tmp\\Geometry.geo",CWD);*//*
				toOpen = (std::string)CWD + "\\tmp\\Geometry.geo"; //newer geo7z format: contain Geometry.geo
				if (!FileUtils::Exist(toOpen)) toOpen = ((std::string)fileName).substr(0, strlen(fileName) - 2); //Inside the zip, try original filename with extension changed from geo7z to geo
				f = new FileReader(toOpen);
			}
			progressDlg->SetMessage("Resetting worker...");
			ResetWorkerStats();
			if (!isGEO7Z) f = new FileReader(fileName);
			mApp->changedSinceSave = TRUE;
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
	} 
	else if (isXML || isXMLzip) {
		xml_document loadXML;
		xml_parse_result parseResult;
		progressDlg->SetVisible(TRUE);
		try {
			if (isXMLzip) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");

				HZIP hz = OpenZip(fileName, 0);
				if (!hz) {
					throw Error("Can't open ZIP file");
				}
				ZIPENTRY ze; GetZipItem(hz, -1, &ze); int numitems = ze.index;
				BOOL notFoundYet = TRUE;
				for (int i = 0; i < numitems && notFoundYet; i++) { //extract first XML file found in ZIP archive
					GetZipItem(hz, i, &ze);
					std::string fileName = ze.name;
					if (fileName.length() >= 4 && !fileName.substr(fileName.length() - 4, 4).compare(".xml")) { //if it's an .xml file
						notFoundYet = FALSE;
						std::string tmpFileName = "tmp/" + fileName;
						UnzipItem(hz, i, tmpFileName.c_str()); //unzip it to tmp directory
						CloseZip(hz);
						progressDlg->SetMessage("Reading and parsing XML file...");
						parseResult = loadXML.load_file(tmpFileName.c_str()); //parse it
					}
				}
				if (notFoundYet) {
					throw Error("No XML file in the ZIP file.");
				}
			}
			ResetWorkerStats();
			progressDlg->SetMessage("Reading and parsing XML file...");
			if (!isXMLzip) parseResult = loadXML.load_file(fileName); //parse it
			if (!parseResult) {
				//Parse error
				std::stringstream err;
				err << "XML parsed with errors.\n";
				err << "Error description: " << parseResult.description() << "\n";
				err << "Error offset: " << parseResult.offset << "\n";
				throw Error(err.str().c_str());
			}

			progressDlg->SetMessage("Building geometry...");
			geom->InsertXML(loadXML, this, progressDlg, newStr);
			mApp->changedSinceSave = TRUE;
			nbHit = 0;
			nbDesorption = 0;
			maxDesorption = 0;
			nbLeakTotal = 0;
			Reload();
		}
		catch (Error &e) {
			//geom->Clear();
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else if(_stricmp(ext,"ase")==0) {
		throw Error("ASE file inserting is not supported.");


	} else {

		SAFE_DELETE(f);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("InsertGeometry(): Invalid file extension [Only txt,stl,geo,geo7z,ase,syn,syn7z,str,xml or zip]");

	}
	
	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
	SAFE_DELETE(f);
	mApp->UpdateFacetlistSelected();	
	mApp->UpdateViewers();




}*/

// -------------------------------------------------------------

void Worker::GetLeak(LEAK *buffer,int *nb) {

	*nb=0;
	if( dpHit ) {
		memcpy(buffer,leakCache,sizeof(LEAK)*NBHLEAK);
		*nb = (int)MIN(nbLastLeaks,NBHLEAK);
	}

}

// -------------------------------------------------------------

void Worker::SetLeak(LEAK *buffer,int *nb,SHGHITS *gHits) { //When loading from file

	if( nb>0 ) {
		memcpy(leakCache,buffer,sizeof(LEAK)*MIN(NBHLEAK,*nb));
		memcpy(gHits->pLeak,buffer,sizeof(LEAK)*MIN(NBHLEAK,*nb));
		gHits->nbLastLeaks = *nb;
	}

}

// -------------------------------------------------------------
void Worker::LoadTexturesSYN(const char *fileName,int version) {

	if (FileUtils::GetExtension(fileName) == "syn") {
		GLProgress *progressDlg = new GLProgress("Loading texture values", "Please wait");
		progressDlg->SetProgress(0.0);
		FileReader *f = NULL;
		try {
			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadTextures(f, progressDlg, dpHit, version);
			RebuildTextures();
		}
		catch (Error &e) {
			char tmp[256];
			sprintf(tmp, "Couldn't load some textures. To avoid continuing a partially loaded state, it is recommended to reset the simulation.\n%s", e.GetMsg());
			GLMessageBox::Display(tmp, "Error while loading textures.", GLDLG_OK, GLDLG_ICONWARNING);
		}
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		SAFE_DELETE(f);
	}
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
		try {
			if (needsReload) RealReload();
			startTime = appTime;
			running = TRUE;

			this->mode = mode;
		
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
		if( (error || done) && running && appTime!=0.0f ) {
			InnerStop(appTime);
			if (error) ThrowSubProcError();
		}

		// Retrieve hit count recording from the shared memory
		if( dpHit ) {

			if( AccessDataport(dpHit) ) {
				BYTE *buffer = (BYTE *)dpHit->buff;

				mApp->changedSinceSave = TRUE;
				// Globals
				SHGHITS *gHits = (SHGHITS *)buffer;

				// Global hits and leaks
				nbHit = gHits->total.nbHit;
				nbAbsorption = gHits->total.nbAbsorbed;
				distTraveledTotal = gHits->distTraveledTotal;


				nbDesorption = gHits->total.nbDesorbed;
				if (nbDesorption && nbTrajPoints) no_scans = (double) nbDesorption / (double) nbTrajPoints;
				else no_scans=1.0;
				totalFlux = gHits->total.fluxAbs;
				totalPower = gHits->total.powerAbs;
				nbLeakTotal = gHits->nbLeakTotal;
				nbHHit = gHits->nbHHit;
				nbLastLeaks = gHits->nbLastLeaks;
				memcpy(hhitCache,gHits->pHits,sizeof(HIT)*NBHHIT);
				memcpy(leakCache,gHits->pLeak,sizeof(LEAK)*NBHHIT);


				// Facets
				int nbFacet = geom->GetNbFacet();
				for(int i=0;i<nbFacet;i++) {    
					Facet *f = geom->GetFacet(i);
					memcpy(&(f->sh.counter),buffer+f->sh.hitOffset,sizeof(SHHITS));
				}
				try {
					geom->BuildTexture(buffer);
				}
				catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(), "Error building texture", GLDLG_OK, GLDLG_ICONERROR);
					ReleaseDataport(dpHit);
					return;
				}
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
			gHits->nbLeakTotal = nbLeakTotal;
			gHits->total.nbDesorbed = nbDesorption;
			gHits->total.nbAbsorbed = nbAbsorption;

			gHits->total.fluxAbs = totalFlux;
			gHits->total.powerAbs = totalPower;
			gHits->distTraveledTotal = distTraveledTotal;

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

	if (!geom->IsLoaded()) {
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		return;
	}

	progressDlg->SetMessage("Creating dataport...");
	// Create the temporary geometry shared structure
	size_t loadSize = geom->GetGeometrySize(&regions, &materials, psi_distr, chi_distr);
	Dataport *loader = CreateDataport(loadDpName,loadSize);
	if (!loader) {
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("Failed to create 'loader' dataport.\nMost probably out of memory.\nReduce number of subprocesses or texture size.");
	}
	progressDlg->SetMessage("Accessing dataport...");
	AccessDataportTimed(loader,3000+nbProcess*(int)((double)loadSize/10000.0));
	progressDlg->SetMessage("Assembling geometry and regions to pass...");
	geom->CopyGeometryBuffer((BYTE *)loader->buff,&regions,&materials,psi_distr,chi_distr,generation_mode,lowFluxMode,lowFluxCutoff);
	progressDlg->SetMessage("Releasing dataport...");
	ReleaseDataport(loader);

	int hitSize = geom->GetHitsSize();
	dpHit = CreateDataport(hitsDpName,hitSize);
	if( !dpHit ) {
		CLOSEDP(loader);


		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("Failed to create 'hits' dataport: out of memory");
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
	try {
		if (needsReload) RealReload();
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
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
			"(Loading continues while this dialog is visible)\n",nbProcess,nbReady,nbError);
		int waitmore = GLMessageBox::Display(tmp,"Info",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK;
		if (waitmore) {
			t=0;
			return Wait(waitState,timeout,prg);
		}
	}

	return ok && !error;

}

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
	return Wait(waitState,3000 + nbProcess * 500,prg);
}

void Worker::ResetWorkerStats() {

	nbAbsorption = 0;
	nbDesorption = 0;
	nbHit = 0;
	nbLeakTotal = 0;

	totalFlux = 0.0;
	totalPower = 0.0;
	distTraveledTotal = 0.0;
	no_scans=1.0;

}

// -------------------------------------------------------------

void Worker::Reset(float appTime) {

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

void Worker::ImportCSV(FileReader *file,std::vector<std::vector<double>>& table){
	table = std::vector<std::vector<double>>(); //reset table
	do {
		std::vector<double> currentRow;
		do {
			currentRow.push_back(file->ReadDouble());
			if (!file->IsEol()) file->ReadKeyword(",");
		} while (!file->IsEol());
		table.push_back(currentRow);
	} while (!file->IsEof());
}

int Worker::GetProcNumber() {
	return nbProcess;
}



DWORD Worker::GetPID(int prIdx) {
	return pID[prIdx];
}


void Worker::AddRegion(const char *fileName,int position) {
	//if (!geom->IsLoaded()) throw Error("Load geometry first!");
	needsReload=TRUE;

	std::string ext = FileUtils::GetExtension(fileName);
	FileReader *f = NULL;
	if(ext=="par" || ext=="PAR" || ext=="param") {

			f = new FileReader(fileName);
			Region_full newtraj;
			if (ext=="par" || ext=="PAR") newtraj.LoadPAR(f);
			else newtraj.LoadParam(f);
			SAFE_DELETE(f);
			newtraj.fileName=fileName;
			if (position==-1) regions.push_back(newtraj);
			else {
				nbTrajPoints-=(int)regions[position].Points.size();
				regions[position]=newtraj;
				//regions[position].Points=newtraj.Points; //need to force because of operator=
			}
			geom->InitializeGeometry(-1,TRUE); //recalculate bounding box
			nbTrajPoints+=(int)newtraj.Points.size();
	}  else {
		SAFE_DELETE(f);
		throw Error("LoadParam(): Invalid file extension [Only par or param]");
	}
}


void Worker::RecalcRegion(int regionId) {
	needsReload=TRUE;
	try {
			//Region_full newtraj;
			//newtraj=regions[regionId]; //copy all except the points
			nbTrajPoints-=(int)regions[regionId].Points.size();
			regions[regionId].Points = std::vector<Trajectory_Point>(); //clear points
			regions[regionId].CalculateTrajectory(1000000); //points calculated
			nbTrajPoints+=(int)regions[regionId].Points.size();
			//regions[regionId]=newtraj;
			//regions[regionId].Points=newtraj.Points; //need to force because of operator=
			geom->InitializeGeometry(-1,TRUE); //recalculate bounding box
		} catch(Error &e) {
			throw e;
		}
}

void Worker::ClearRegions() {
	//if ((int)regions.size()>0) regions.clear();
	regions=std::vector<Region_full>();
	geom->InitializeGeometry(-1,TRUE); //recalculate bounding box
	nbTrajPoints=0;
	Reload();
}

void Worker::RemoveRegion(int index) {
	nbTrajPoints -= (int)regions[index].Points.size();
	/*
	//Explicit removal as Region_full doesn't copy Points for some reason
	for (size_t i = index; i < (regions.size() - 1); i++) { //Copy next
		regions[i].Points = regions[i + 1].Points; //Points
		regions[i] = regions[i + 1]; //Everything else
	}
	regions.erase(regions.end() - 1); //delete last
	*/
	regions.erase(regions.begin() + index);
	geom->InitializeGeometry(-1, TRUE); //recalculate bounding box
}

void Worker::AddMaterial(std::string *fileName){
	Material result;
	char tmp[512];
	sprintf(tmp,"param\\Materials\\%s",fileName->c_str());
	FileReader *f=new FileReader(tmp);
	result.LoadMaterialCSV(f);
	int lastindex = fileName->find_last_of(".");
	result.name = fileName->substr(0, lastindex);
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

void Worker::RebuildTextures() {
	if (AccessDataport(dpHit)) {
		BYTE *buffer = (BYTE *)dpHit->buff;
		try{ geom->BuildTexture(buffer); }
		catch (Error &e) {
			ReleaseDataport(dpHit);
			throw e;
		}
		ReleaseDataport(dpHit);
	}
}