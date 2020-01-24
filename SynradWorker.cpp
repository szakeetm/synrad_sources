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
#include "Worker.h"
#include "Facet_shared.h"
#include "SynradGeometry.h"
#include "SynradDistributions.h"
#include "SynradFacet.h"
#include "GLApp/GLApp.h"
#include "GLApp/GLMessageBox.h"
#include <math.h>
#include <stdlib.h>
#include <Process.h>
#include "GLApp/GLUnitDialog.h"
#include "File.h" //GetFileName
#include "GLApp/MathTools.h"
#ifdef MOLFLOW
#include "MolFlow.h"
#endif

#ifdef SYNRAD
#include "SynRad.h"
#endif

#include <direct.h>

#include "ziplib/ZipArchive.h"
#include "ziplib/ZipArchiveEntry.h"
#include "ziplib/ZipFile.h"
#include "File.h" //File utils (Get extension, etc)
#include <cereal/types/vector.hpp>
using namespace pugi;

#ifdef MOLFLOW
extern MolFlow *mApp;
#endif

#ifdef SYNRAD
extern SynRad*mApp;
#endif

Worker::Worker() {

	pid = _getpid();
	sprintf(ctrlDpName,"SNRDCTRL%d",pid);
	sprintf(loadDpName,"SNRDLOAD%d",pid);
	sprintf(hitsDpName,"SNRDHITS%d",pid);
	sprintf(logDpName, "SNRDLOG%d", pid);
	
	globalHitCache.distTraveledTotal=0.0;
	ontheflyParams.nbProcess = 0;
	ontheflyParams.enableLogging = false;
	ontheflyParams.desorptionLimit = 0;
	ontheflyParams.lowFluxCutoff = 1E-7;
	ontheflyParams.lowFluxMode = false;
	ontheflyParams.generation_mode=SYNGEN_MODE_FLUXWISE;
	wp.newReflectionModel = false;
	ResetWorkerStats();
	geom = new SynradGeometry();
	regions = std::vector<Region_full>();
	
	dpControl = NULL;
	dpHit = NULL;
	dpLog = NULL;
	globalHitCache.hitCacheSize = 0;
    globalHitCache.leakCacheSize = 0;

    globalHitCache.nbLeakTotal = 0;
	startTime = 0.0f;
	stopTime = 0.0f;
	simuTime = 0.0f;
	wp.nbTrajPoints=0;
	isRunning = false;
	calcAC = false; //not used, reserved for shared function

	strcpy(fullFileName,"");

}

SynradGeometry *Worker::GetSynradGeometry() {
	return geom;
}

void Worker::SaveGeometry(std::string fileName,GLProgress *prg,bool askConfirm,bool saveSelected,bool autoSave,bool crashSave) {

	try {
		if (needsReload&&(!crashSave && !saveSelected)) RealReload();

	} catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"Error reloading worker. Trying crash save:\n%s",e.GetMsg());
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		crashSave=true;
	}
    std::string compressCommandLine; //compress.exe command line
	/*char fileNameWithGeo[2048]; //file name with .geo extension (instead of .geo7z)
	char fileNameWithGeo7z[2048];*/
    std::string fileNameWithSyn; //file name with .syn extension (instead of .syn7z)
    std::string fileNameWithSyn7z;
    std::string fileNameWithoutExtension; //file name without extension
	//char *ext = fileName+strlen(fileName)-4;

    std::string ext = FileUtils::GetExtension(fileName);
    std::string dir = FileUtils::GetPath(fileName);

    bool ok = true;
    if (ext.empty()) {
        fileName = fileName + (mApp->compressSavedFiles ? ".syn7z" : ".syn");
        ext = FileUtils::GetExtension(fileName);
        if (!autoSave && FileUtils::Exist(fileName)) {
            char tmp[1024];
            sprintf(tmp, "Overwrite existing file ?\n%s", fileName.c_str());
            if (askConfirm) ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
        }
    }

	// Read a file
	FileWriter *f = NULL;

	/*bool isGEO = _stricmp(ext,"geo")==0;
	bool isGEO7Z = _stricmp(ext,"geo7z")==0;*/

    bool isTXT = Contains({ "txt","TXT" },ext);
    bool isSTR = Contains({ "str","STR" }, ext);
    bool isSYN = ext == "syn";
    bool isSYN7Z = ext == "syn7z";

	if(isTXT || isSYN || isSYN7Z || isSTR) {
#ifdef _WIN32
		if (WAIT_TIMEOUT==WaitForSingleObject(mApp->compressProcessHandle,0)) {
			GLMessageBox::Display("Compressing a previous save file is in progress. Wait until that finishes "
				"or close process \"compress.exe\"\nIf this was an autosave attempt,"
				"you have to lower the autosave frequency.","Can't save right now.",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
#endif

        if (isSYN) {
            fileNameWithoutExtension = fileName.substr(0,fileName.length()-4);
            fileNameWithSyn7z = fileName + "7z";

        } else if (isSYN7Z) {
            fileNameWithoutExtension = fileName.substr(0, fileName.length() - 6);
            fileNameWithSyn = fileName.substr(0, fileName.length() - 2);
            fileNameWithSyn7z = fileName;
            std::ostringstream tmp;
            tmp << "A .syn file of the same name exists. Overwrite that file ?\n" << fileNameWithSyn;
            if (!autoSave && FileUtils::Exist(fileNameWithSyn)) {
                ok = (GLMessageBox::Display(tmp.str().c_str(), "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
            }
		}

		if(!autoSave && ok && FileUtils::Exist(fileName) ) {
            char tmp[1024];
            sprintf(tmp,"Overwrite existing file ?\n%s",fileName.c_str());
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
				geom->loaded_desorptionLimit = ontheflyParams.desorptionLimit;
				if( isTXT ) geom->SaveTXT(f,dpHit,saveSelected);
				else if( isSYN || isSYN7Z ) {
					//Save regions
					for (int i=0;i<(int)regions.size();i++) {
						if (FileUtils::GetExtension(regions[i].fileName)=="par") { //fileName doesn't end with .param, save a param version
                            char tmp[1024];
                            sprintf(tmp,"%sam",regions[i].fileName.c_str());
							regions[i].fileName=tmp;
						}	

						SaveRegion((char*)regions[i].fileName.c_str(),i,true); //save with forced overwrite
					}

					geom->SaveSYN(f,prg,dpHit,saveSelected,globalHitCache.leakCache,(&globalHitCache.leakCacheSize),globalHitCache.hitCache,&(globalHitCache.hitCacheSize),crashSave);
				}
			}
			if (!autoSave && !saveSelected) {
                strcpy(fullFileName, fileName.c_str());
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

				compressProcessHandle=OpenProcess(PROCESS_ALL_ACCESS, true, procId);
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
				
				std::vector<std::string> paths;
				std::vector<std::string> fileNames;
				std::vector<size_t> regionIds;
				bool foundConflict = false;

				for (size_t i = 0; !foundConflict && i < regions.size(); i++) {
					Region_full& r = regions[i];

					foundConflict = CheckFilenameConflict(r.fileName, i, paths, fileNames, regionIds);

					if (!foundConflict && r.MAGXfileName != "") {
						foundConflict = CheckFilenameConflict(r.MAGXfileName, i, paths, fileNames, regionIds);
					}
					if (!foundConflict && r.MAGYfileName != "") {
						foundConflict = CheckFilenameConflict(r.MAGYfileName, i, paths, fileNames, regionIds);
					}
					if (!foundConflict && r.MAGZfileName != "") {
						foundConflict = CheckFilenameConflict(r.MAGZfileName, i, paths, fileNames, regionIds);
					}
					if (!foundConflict && r.BXYfileName != "") {
						foundConflict = CheckFilenameConflict(r.BXYfileName, i, paths, fileNames, regionIds);
					}
				}

				std::ostringstream commandLine;
				if (!foundConflict) {
					
					commandLine << "compress.exe \"" << fileNameWithSyn << "\" Geometry.syn";
					if (paths.size() > 0) {
						std::ofstream includeFile;
						includeFile.open("compressionInclude.txt");
						if (!includeFile.is_open()) {
							std::ostringstream msg;
							msg << "compressionInclude.txt cannot be opened for writing. Maybe already in use?\n\n";
							msg << "File will be saved as uncompressed syn file\n\n";
							GLMessageBox::Display(msg.str().c_str(), "Name conflict", GLDLG_OK, GLDLG_ICONWARNING);
							fileName = fileNameWithSyn;
							foundConflict = true;
						}
						else {
							for (auto path : paths) {
								includeFile << path << "\n";
							}
							includeFile.close();
							commandLine << " @compressionInclude.txt";
						}
					}
				}

				if (!foundConflict) {
					int procId = StartProc(commandLine.str().c_str(), STARTPROC_BACKGROUND);
					mApp->compressProcessHandle = OpenProcess(PROCESS_ALL_ACCESS, true, procId);
					fileName = fileNameWithSyn7z;
				}
				else {
					fileName=fileNameWithSyn;
				}
			} else {
				GLMessageBox::Display("compress.exe (part of Synrad) not found.\n Will save as uncompressed SYN file.","Compressor not found",GLDLG_OK,GLDLG_ICONERROR);
				fileName=fileNameWithSyn;
			}
		} else if (ok && isSYN) fileName=fileNameWithSyn;
		if (!autoSave && !saveSelected) {
			SetCurrentFileName(fileName.c_str());
			mApp->UpdateTitle();
		}
}

/*
void Worker::ExportDesorption(char *fileName,bool selectedOnly,int mode,double eta0,double alpha,const Distribution2D &distr) {

	// Read a file

	FILE *f = NULL;
	bool ok = true;

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
*/

void Worker::LoadGeometry(const std::string& fileName, bool insert, bool newStr) {
	if (!insert) {
		needsReload=true;
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
	progressDlg->SetVisible(true);
	progressDlg->SetProgress(0.0);

    ResetWorkerStats();

	if (!insert) {
		//Clear hits and leaks cache
		//Reset leak and hit cache
		globalHitCache.leakCacheSize = 0;
        //SetLeakCache(globalHitCache.leakCache, &(globalHitCache.leakCacheSize), dpHit); //will only write leakCacheSize
        globalHitCache.hitCacheSize = 0;
		//SetHitCache(globalHitCache.hitCache, (&globalHitCache.hitCacheSize), dpHit); //will only write hitCacheSize
	}

	if(ext=="txt" || ext=="TXT") {

		try {
			//if (!insert) ResetWorkerStats();
			f = new FileReader(fileName);
			if (!insert) {
				geom->LoadTXT(f, progressDlg,this);
				strcpy(fullFileName, fileName.c_str());
			}
			else { //insert
				mApp->changedSinceSave = true;
				geom->InsertTXT(f, progressDlg, newStr);
				Reload();
			}
			SAFE_DELETE(f);
		} catch(Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(false);
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
				progressDlg->SetVisible(true);
				ResetWorkerStats();				
				progressDlg->SetMessage("Reading geometry...");
				f = new FileReader(fileName);
				if (!insert) {					
					geom->LoadSTL(f, progressDlg, scaleFactor);
					strcpy(fullFileName, fileName.c_str());
					mApp->DisplayCollapseDialog();
				}
				else { //insert
					mApp->changedSinceSave = true;
					geom->InsertSTL(f, progressDlg, scaleFactor, newStr);
					Reload();
				}
				SAFE_DELETE(f);
			}
		} catch(Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);
			throw e;

		}

	} else if(ext=="str") {
		if (insert) throw Error("STR file inserting is not supported.");
		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			progressDlg->SetVisible(true);
			geom->LoadSTR(f,progressDlg);
			SAFE_DELETE(f);
			strcpy(fullFileName,fileName.c_str());
		}
		catch(Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else if (ext == "syn" || ext == "syn7z") { //Synrad file
		int version;
		progressDlg->SetVisible(true);
		try {
            if (ext == "syn7z") {
                //decompress file
                progressDlg->SetMessage("Decompressing file...");
                f = ExtractFrom7zAndOpen(fileName, "Geometry.syn");
            }
            else {
                f = new FileReader(fileName);  //original file opened
            }

			std::vector<std::string> regionsToLoad;
			
			LEAK loaded_leakCache[LEAKCACHESIZE];
			size_t loaded_nbLeak;
			HIT hitCache[HITCACHESIZE];
			//size_t loaded_nbHit;
			if (!insert) {
				progressDlg->SetMessage("Resetting worker...");
				ResetWorkerStats();
				
				regionsToLoad = geom->LoadSYN(f, progressDlg, &version, this);
				//copy temp values from geom to worker. They will be sent to shared memory in LoadTextures() which connects to dpHit
				globalHitCache.nbLeakTotal = geom->loaded_nbLeak;
                globalHitCache.globalHits.hit.nbMCHit = geom->loaded_nbMCHit;
                globalHitCache.globalHits.hit.nbHitEquiv = geom->loaded_nbHitEquiv;
                globalHitCache.globalHits.hit.nbDesorbed = geom->loaded_nbDesorption;
                globalHitCache.globalHits.hit.nbAbsEquiv = geom->loaded_nbAbsEquiv;
                globalHitCache.distTraveledTotal = geom->loaded_distTraveledTotal;
				ontheflyParams.desorptionLimit = geom->loaded_desorptionLimit;
                globalHitCache.globalHits.hit.fluxAbs = geom->loaded_totalFlux;
                globalHitCache.globalHits.hit.powerAbs = geom->loaded_totalPower;
				no_scans = geom->loaded_no_scans;
			}
			else { //insert
				regionsToLoad = geom->InsertSYN(f, progressDlg, newStr);
			}
			
			//Load regions
			if (regionsToLoad.size()>0) {
				char tmp[256];
				sprintf(tmp,"This geometry refers to %zd regions. Load them now?",regionsToLoad.size());
				bool loadThem = ( GLMessageBox::Display(tmp,"File load",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK );
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
				geom->LoadProfileSYN(f, dpHit, version);
				geom->LoadSpectrumSYN(f, dpHit, version);
				//SetLeakCache(loaded_leakCache, &loaded_nbLeak, dpHit);
				//SetHitCache(hitCache, &(globalHitCache.hitCacheSize), dpHit);
				progressDlg->SetMessage("Loading texture values...");
				LoadTexturesSYN(f, version);
				strcpy(fullFileName, fileName.c_str());
			}
			SAFE_DELETE(f);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			
			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	} else if(ext=="geo" || ext=="geo7z") {
		std::string toOpen;
		int version;
		progressDlg->SetVisible(true);
		try {
			if (ext=="geo7z") {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName.c_str());
				system(tmp);
				toOpen = (std::string)CWD + "\\tmp\\Geometry.geo"; //newer geo7z format: contains Geometry.geo
				if (!FileUtils::Exist(toOpen)) toOpen = (fileName).substr(0, strlen(fileName.c_str()) - 2); //Inside the zip, try original filename with extension changed from geo7z to geo
				f = new FileReader(toOpen);
			}
			else { //not geo7z
				toOpen = fileName;
				f = new FileReader(fileName); //geo file, open it directly
			}
			progressDlg->SetMessage("Resetting worker...");
			ResetWorkerStats();
			if (!insert) {
				geom->LoadGEO(f, progressDlg, &version, this);
				SAFE_DELETE(f);
				ontheflyParams.desorptionLimit = 0;
				progressDlg->SetMessage("Reloading worker with new geometry...");
				RealReload();
                SendToHitBuffer(); //Global hit counters and hit/leak cache
                SendFacetHitCounts(); // From facetHitCache to dpHit's const.flow counter

                strcpy(fullFileName,fileName.c_str());
			}
			else { //insert
				mApp->changedSinceSave = true;
				geom->InsertGEO(f,progressDlg,newStr);
				SAFE_DELETE(f);
				Reload();
			}
		} catch(Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			//if (isGEO7Z) remove(tmp2);
			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	}
	else if (ext == "xml" || ext == "zip") { //XML file, optionally in ZIP container
		xml_document loadXML;
		xml_parse_result parseResult;
		progressDlg->SetVisible(true);
		try {
			if (ext == "zip") { //compressed in ZIP container
				//decompress file
				progressDlg->SetMessage("Decompressing file...");

                ZipArchive::Ptr zip = ZipFile::Open(fileName);
                if (zip == nullptr) {
                    throw Error("Can't open ZIP file");
                }
                size_t numitems = zip->GetEntriesCount();
                bool notFoundYet = true;
                for (int i = 0; i < numitems && notFoundYet; i++) { //extract first xml file found in ZIP archive
                    auto zipItem = zip->GetEntry(i);
                    std::string zipFileName = zipItem->GetName();

					if (FileUtils::GetExtension(zipFileName) == "xml") { //if it's an .xml file
						notFoundYet = false;

                        FileUtils::CreateDir("tmp");// If doesn't exist yet

                        std::string tmpFileName = "tmp/" + zipFileName;
                        ZipFile::ExtractFile(fileName, zipFileName, tmpFileName);
						progressDlg->SetMessage("Reading and parsing XML file...");
						parseResult = loadXML.load_file(tmpFileName.c_str()); //load and parse it
					}
				}
				if (notFoundYet) {
					throw Error("Didn't find any XML file in the ZIP file.");
				}
			}
			else parseResult = loadXML.load_file(fileName.c_str()); //parse xml file directly

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
				ontheflyParams.desorptionLimit = 0;
				progressDlg->SetMessage("Reloading worker with new geometry...");
				RealReload();
                strcpy(fullFileName, fileName.c_str());

				geom->UpdateName(fileName.c_str());
			}
			else { //insert
				mApp->changedSinceSave = true;
				geom->InsertXML(loadXML, this, progressDlg, newStr);
				Reload();
			}
		}
		catch (Error &e) {
			if (!insert) geom->Clear();
			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	}
	else if (ext == "ase") {
		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			progressDlg->SetVisible(true);
			geom->LoadASE(f,progressDlg);
			SAFE_DELETE(f);
			strcpy(fullFileName,fileName.c_str());

		} catch(Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);
			throw e;
		}
		
	} else {
		progressDlg->SetVisible(false);
		SAFE_DELETE(progressDlg);
		throw Error("LoadGeometry(): Invalid file extension [Only txt,geo,geo7z,ase,stl or str]");
	}
	progressDlg->SetVisible(false);
	SAFE_DELETE(progressDlg);
}

void Worker::LoadTexturesSYN(FileReader* f,int version) {
		GLProgress *progressDlg = new GLProgress("Loading texture values", "Please wait");
		progressDlg->SetProgress(0.0);
		try {
			progressDlg->SetVisible(true);
			geom->LoadTextures(f, progressDlg, dpHit, version);
			RebuildTextures();
		}
		catch (Error &e) {
			char tmp[256];
			sprintf(tmp, "Couldn't load some textures. To avoid continuing a partially loaded state, it is recommended to reset the simulation.\n%s", e.GetMsg());
			GLMessageBox::Display(tmp, "Error while loading textures.", GLDLG_OK, GLDLG_ICONWARNING);
		}
		progressDlg->SetVisible(false);
		SAFE_DELETE(progressDlg);
}

void Worker::InnerStop(float appTime) {

	stopTime =appTime;
	simuTime+=appTime-startTime;
	isRunning  = false;

}

void Worker::StartStop(float appTime) {

	if( isRunning )  {

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
			isRunning = true;

			//this->mode = mode;
		
			Start();
		} catch(Error &e) {
			isRunning = false;
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

void Worker::RealReload(bool sendOnly) { //Sharing geometry with workers

	if(ontheflyParams.nbProcess==0) return;

	GLProgress *progressDlg = new GLProgress("Asking subprocesses to clear geometry...","Passing Geometry to workers");
	progressDlg->SetVisible(true);
	progressDlg->SetProgress(0.0);

    //this->regions = std::vector<Region_full>();
    //this->materials = std::vector<Material>();
    //this->psi_distro = std::vector<std::vector<double>>();
    //this->chi_distros = std::vector<std::vector<std::vector<double>>>();
    //this->parallel_polarization = std::vector<std::vector<double>>();

	// Clear geometry
	CLOSEDP(dpHit);
	CLOSEDP(dpLog);
	std::string lastError;
	if (!ExecuteAndWait(COMMAND_CLOSE, PROCESS_READY)) {
		progressDlg->SetVisible(false);
		SAFE_DELETE(progressDlg);
		ThrowSubProcError(lastError);
	}

	if (!geom->IsLoaded()) {
		progressDlg->SetVisible(false);
		SAFE_DELETE(progressDlg);
		return;
	}

	progressDlg->SetMessage("Creating dataport...");
	// Create the temporary geometry shared structure
	//size_t loadSize = geom->GetGeometrySize(regions, materials, psi_distro, chi_distros, parallel_polarization);
	if (ontheflyParams.enableLogging) {
		dpLog = CreateDataport(logDpName, sizeof(size_t) + sizeof(ParticleLoggerItem)*ontheflyParams.logLimit);
		if (!dpLog)
			throw Error("Failed to create 'dpLog' dataport.\nMost probably out of memory.\nReduce number of logged particles in Particle Logger.");
		//*((size_t*)dpLog->buff) = 0; //Automatic 0-filling
	}

    std::string loaderString = SerializeForLoader().str();

    //size_t loadSize = geom->GetGeometrySize();
    //Dataport *loader = CreateDataport(loadDpName, loadSize);

    size_t loadSize = loaderString.size();

	Dataport *loader = CreateDataport(loadDpName,loadSize);
	if (!loader) {
		progressDlg->SetVisible(false);
		SAFE_DELETE(progressDlg);
		throw Error("Failed to create 'loader' dataport.\nMost probably out of memory.\nReduce number of subprocesses or texture size.");
	}
	progressDlg->SetMessage("Accessing dataport...");
	AccessDataportTimed(loader,(DWORD)(3000.0+(double)(ontheflyParams.nbProcess*loadSize)/10000.0));
	progressDlg->SetMessage("Assembling geometry and regions to pass...");
	//this->ontheflyParams;
	/*geom->CopyGeometryBuffer((BYTE *)loader->buff,regions,materials,psi_distro,chi_distros,
		parallel_polarization,wp.newReflectionModel,ontheflyParams);
	*/

    BYTE* buffer = (BYTE*)loader->buff;
    //memcpy(loader->buff, loaderString.c_str(), loadSize);
    std::copy(loaderString.begin(), loaderString.end(), buffer);
	 progressDlg->SetMessage("Releasing dataport...");
	ReleaseDataport(loader);

    if (!sendOnly) {
        size_t hitSize = geom->GetHitsSize();
        dpHit = CreateDataport(hitsDpName, hitSize);
        if (!dpHit) {
            CLOSEDP(loader);

            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            throw Error("Failed to create 'hits' dataport: out of memory");
        }
    }
	// Compute number of max desorption per process
	if(AccessDataportTimed(dpControl, (DWORD)(3000.0 + (double)(ontheflyParams.nbProcess*loadSize) / 10000.0))) {
		SHCONTROL *m = (SHCONTROL *)dpControl->buff;
		size_t common = ontheflyParams.desorptionLimit / ontheflyParams.nbProcess;
		size_t remain = (ontheflyParams.desorptionLimit % ontheflyParams.nbProcess);
		for(size_t i=0;i<ontheflyParams.nbProcess;i++) {
			m->cmdParam2[i] = common;
			if(i<remain) m->cmdParam2[i]++;
		}
		ReleaseDataport(dpControl);
	}

	// Load geometry
	progressDlg->SetMessage("Waiting for subprocesses to load geometry...");
	std::string errorMsg;
	if( !ExecuteAndWait(COMMAND_LOAD,PROCESS_READY,loadSize) ) {
		CLOSEDP(loader);
		char errMsg[1024];
		sprintf(errMsg,"Failed to send geometry to sub process:\n%s",GetErrorDetails());
		GLMessageBox::Display(errMsg,"Warning (Load)",GLDLG_OK,GLDLG_ICONWARNING);

		progressDlg->SetVisible(false);
		SAFE_DELETE(progressDlg);
		return;
	}

	progressDlg->SetProgress(1.0);
	progressDlg->SetMessage("Sending hits...");
	//Send hit counts
	try {
		WriteHitBuffer();

	} catch(Error &e) {
		// Geometry not loaded !
		CLOSEDP(dpHit);
		CLOSEDP(loader);
		progressDlg->SetVisible(false);
		SAFE_DELETE(progressDlg);
		throw e;
	}

	progressDlg->SetMessage("Closing dataport...");
	CLOSEDP(loader);
	needsReload=false;
	progressDlg->SetVisible(false);
	SAFE_DELETE(progressDlg);
}

void Worker::ClearHits(bool noReload) {
	try {
		if (!noReload && needsReload) RealReload();
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (dpHit) {
		AccessDataport(dpHit);
		memset(dpHit->buff, 0, geom->GetHitsSize());
		ReleaseDataport(dpHit);
	}

}

/**
* \brief Resets workers global hit cache
* Reset function mainly used for initialisation / reload procedures
*/
void Worker::ResetWorkerStats() {

    memset(&globalHitCache, 0, sizeof(GlobalHitBuffer));

	no_scans=1.0;

}

void Worker::Start() {

	if(ontheflyParams.nbProcess==0 )
		throw Error("No sub process found. (Simulation not available)");

	if( !ExecuteAndWait(COMMAND_START,PROCESS_RUN,ontheflyParams.generation_mode) ) //here the param is ignored
		ThrowSubProcError();
}

void Worker::AddRegion(const char *fileName,int position) {
	//if (!geom->IsLoaded()) throw Error("Load geometry first!");
	needsReload=true;

	std::string ext = FileUtils::GetExtension(fileName);
	if(ext=="par" || ext=="PAR" || ext=="param") {
        Region_full newtraj;

        FileReader* fr = new FileReader(fileName);

        if (ext == "par" || ext == "PAR")
            newtraj.LoadPAR(fr);
        else {
            newtraj.LoadParam(fr);
        }
        SAFE_DELETE(fr);

        newtraj.fileName=fileName;
			if (position==-1) regions.push_back(newtraj);
			else {
				wp.nbTrajPoints-=(int)regions[position].Points.size();
				regions[position]=newtraj;
				//regions[position].Points=newtraj.Points; //need to force because of operator=
			}
			geom->RecalcBoundingBox(); //recalculate bounding box
			wp.nbTrajPoints+=(int)newtraj.Points.size();
	}  else {
		throw Error("LoadParam(): Invalid file extension [Only par or param]");
	}
}

void Worker::RecalcRegion(int regionId) {
	needsReload=true;
	try {
			//Region_full newtraj;
			//newtraj=regions[regionId]; //copy all except the points
			wp.nbTrajPoints-=(int)regions[regionId].Points.size();
			regions[regionId].Points = std::vector<Trajectory_Point>(); //clear points
			regions[regionId].CalculateTrajectory(1000000); //points calculated
			wp.nbTrajPoints+=(int)regions[regionId].Points.size();
			//regions[regionId]=newtraj;
			//regions[regionId].Points=newtraj.Points; //need to force because of operator=
			geom->RecalcBoundingBox(); //recalculate bounding box
		} catch(Error &e) {
			throw e;
		}
}

void Worker::ClearRegions() {
	//if ((int)regions.size()>0) regions.clear();
	regions=std::vector<Region_full>();
	geom->RecalcBoundingBox(); //recalculate bounding box
	wp.nbTrajPoints=0;
	Reload();
}

void Worker::RemoveRegion(int index) {
	wp.nbTrajPoints -= (int)regions[index].Points.size();
	/*
	//Explicit removal as Region_full doesn't copy Points for some reason
	for (size_t i = index; i < (regions.size() - 1); i++) { //Copy next
		regions[i].Points = regions[i + 1].Points; //Points
		regions[i] = regions[i + 1]; //Everything else
	}
	regions.erase(regions.end() - 1); //delete last
	*/
	regions.erase(regions.begin() + index);
	geom->RecalcBoundingBox(); //recalculate bounding box
}

void Worker::AddMaterial(std::string *fileName){
	Material result;
	char tmp[512];
	sprintf(tmp,"param\\Materials\\%s",fileName->c_str());
	FileReader *f=new FileReader(tmp);
	result.LoadMaterialCSV(f);
	delete f;
	size_t lastindex = fileName->find_last_of(".");
	result.name = fileName->substr(0, lastindex);
	materials.push_back(result);
}

void Worker::SaveRegion(const char *fileName,int position,bool overwrite) {
	char tmp[512];

	// Read a file
	FILE *f = NULL;
	bool ok = true;
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

bool Worker::CheckFilenameConflict(const std::string & newPath, const size_t & regionId, std::vector<std::string>& paths, std::vector<std::string>& fileNames, std::vector<size_t>& regionIds)
{
	std::string newFileName = FileUtils::GetFilename(newPath);
	size_t firstFound = FirstIndex(fileNames, newFileName);

	if (firstFound<fileNames.size()) { //There's already a file like this
		if (paths[firstFound] != newPath) {
			std::ostringstream msg;
			msg << "You have multiple files with same names but different locations.\nThey can't be packed in the same syn7z (compressed) file.\n\n";
			msg << "Region " << regionIds[firstFound] + 1 << ": " << paths[firstFound] << "\n";
			msg << "Region " << regionId + 1 << ": " << newPath << "\n\n";
			msg << "Until the name conflict is resolved, the file will be saved as an uncompressed syn file.";
			GLMessageBox::Display(msg.str().c_str(), "Name conflict", GLDLG_OK, GLDLG_ICONWARNING);
			return true;
		}
		else {
			//Same file, simply don't add again
		}
	}
	else {
		paths.push_back(newPath);
		fileNames.push_back(newFileName);
		regionIds.push_back(regionId);
	}
	return false;
}

bool EndsWithPar(const char* s)
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

void Worker::WriteHitBuffer() {
	//if (!needsReload) {
	if (dpHit) {
		if (AccessDataport(dpHit)) {

			// Store initial hit counts in the shared memory
			BYTE *pBuff = (BYTE *)dpHit->buff;
			//memset(pBuff, 0, geom->GetHitsSize());

			GlobalHitBuffer *gHits = (GlobalHitBuffer *)pBuff;

			gHits->globalHits.hit.nbMCHit = globalHitCache.globalHits.hit.nbMCHit;
			gHits->globalHits.hit.nbHitEquiv = globalHitCache.globalHits.hit.nbHitEquiv;
			gHits->nbLeakTotal = globalHitCache.nbLeakTotal;
			gHits->globalHits.hit.nbDesorbed = globalHitCache.globalHits.hit.nbDesorbed;
			gHits->globalHits.hit.nbAbsEquiv = globalHitCache.globalHits.hit.nbAbsEquiv;
			gHits->distTraveledTotal = globalHitCache.distTraveledTotal;
			gHits->globalHits.hit.fluxAbs = globalHitCache.globalHits.hit.fluxAbs;
			gHits->globalHits.hit.powerAbs = globalHitCache.globalHits.hit.powerAbs;

			size_t nbFacet = geom->GetNbFacet();
			for (int i = 0; i<nbFacet; i++) {
				Facet *f = geom->GetFacet(i);
				memcpy(pBuff + f->sh.hitOffset, &(f->facetHitCache), sizeof(FacetHitBuffer));
			}
			ReleaseDataport(dpHit);
		}
		else {
			throw Error("Failed to initialise 'hits' dataport");
		}

	}
}


/**
* \brief Function that updates the global hit counter with the cached value + releases the mutex
* \param skipFacetHits TODO: check if not necessary anymore
*
* Send total and facet hit counts to subprocesses
*/
void Worker::SendToHitBuffer() {
    if (dpHit) {
        if (AccessDataport(dpHit)) {

            GlobalHitBuffer *gHits = (GlobalHitBuffer *)dpHit->buff;

            *gHits = globalHitCache;

            ReleaseDataport(dpHit);

        }
        else {
            throw Error("Failed to initialize 'hits' dataport");
        }
    }
}

/**
* \brief Saves current facet hit counter from cache to results
*/
void Worker::SendFacetHitCounts() {
    // Store initial hit counts in the shared memory
    BYTE *pBuff = (BYTE *)dpHit->buff;

    size_t nbFacet = geom->GetNbFacet();
    for (size_t i = 0; i < nbFacet; i++) {
        Facet *f = geom->GetFacet(i);
        *((FacetHitBuffer*)(pBuff + f->sh.hitOffset )) = f->facetHitCache; //Only const.flow
        memcpy(pBuff + f->sh.hitOffset, &(f->facetHitCache), sizeof(FacetHitBuffer));
    }
}

/*void Worker::SendLeakCache(Dataport* dpHit) { //From worker.globalhitCache to dpHit
    if (dpHit) {
        AccessDataport(dpHit);
        GlobalHitBuffer *gHits = (GlobalHitBuffer *)dpHit->buff;
        size_t nbCopy = Min(LEAKCACHESIZE, globalHitCache.leakCacheSize);
        gHits->leakCache = globalHitCache.leakCache;
        gHits->lastLeakIndex = nbCopy-1;
        gHits->leakCacheSize = globalHitCache.leakCacheSize;
        ReleaseDataport(dpHit);
    }
}

void Worker::SendHitCache(Dataport* dpHit) { //From worker.globalhitCache to dpHit
    if (dpHit) {
        AccessDataport(dpHit);
        GlobalHitBuffer *gHits = (GlobalHitBuffer *)dpHit->buff;
        size_t nbCopy = Min(HITCACHESIZE, globalHitCache.hitCacheSize);
        gHits->hitCache = globalHitCache.hitCache;
        gHits->lastHitIndex = nbCopy - 1;
        gHits->hitCacheSize = globalHitCache.hitCacheSize;
        ReleaseDataport(dpHit);
    }
}*/

/**
* \brief Serialization function for a binary cereal archive for the worker attributes
* \return output string stream containing the result of the archiving
*/
std::ostringstream Worker::SerializeForLoader()
{
    //std::ofstream is("data.json");
    //cereal::JSONOutputArchive outputarchive( is );
    std::ostringstream result;
    cereal::BinaryOutputArchive outputarchive(result);

    //sh.nbRegion = regions.size();
    //sh.newReflectionModel = newReflectionModel;

    wp.nbRegion = regions.size();
    //wp.newReflectionModel = newReflectionModel;

    // Build shared buffer for trajectory (see Shared.h)
    for (auto& reg : regions) {
        reg.params.nbDistr_MAG = Vector3d((int)reg.Bx_distr.GetSize(), (int)reg.By_distr.GetSize(), (int)reg.Bz_distr.GetSize());
        reg.params.nbPointsToCopy = reg.Points.size();
    }

    outputarchive(
            CEREAL_NVP(wp),
            CEREAL_NVP(ontheflyParams),
            CEREAL_NVP(regions), // full
            CEREAL_NVP(materials),
            CEREAL_NVP(psi_distro),
            CEREAL_NVP(chi_distros),
            CEREAL_NVP(parallel_polarization)
    ); //Worker

    geom->SerializeForLoader(outputarchive);

    return result;
}