/*
File:        ExportDesorption.cpp
Description: Do the photon dose -> outgassing conversion before exporting
Program:     SynRad

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*/

#include "ExportDesorption.h"
#include "Facet.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "SynRad.h"
#include "GLApp/GLFileBox.h"

extern SynRad *mApp;

ExportDesorption::ExportDesorption(Geometry *g,Worker *w):GLWindow() {

	int wD = 500;
	int hD = 220;
	fileName = NULL;
	fileLoaded = false;
	mode=2; //default mode: use equation
	conversionDistr = new Distribution2D(1); //placeholder

	SetTitle("Export desorption map");

	toggle1=new GLToggle(0,"No conversion");
	toggle1->SetBounds(10,5,170,18);
	Add(toggle1);

	toggle2=new GLToggle(0,"Eta =           * Dose ^");
	toggle2->SetBounds(10,40,170,18);
	toggle2->SetState(true);
	Add(toggle2);

	eta0Field = new GLTextField(0,"1");
	eta0Field->SetBounds(60,40,30,18);
	Add(eta0Field);

	alphaField = new GLTextField(0,"-1");
	alphaField->SetBounds(145,40,30,18);
	Add(alphaField);

	toggle3 = new GLToggle(0,"File");
	toggle3->SetBounds(10,75,40,18);
	Add(toggle3);

	browseButton = new GLButton(0,"Browse...");
	browseButton->SetBounds(60,73,60,21);
	Add(browseButton);

	fileNameLabel = new GLLabel("No file loaded.");
	fileNameLabel->SetBounds(60,100,170,18);
	Add(fileNameLabel);

	fileInfoLabel = new GLLabel("");
	fileInfoLabel->SetBounds(60,125,170,18);
	Add(fileInfoLabel);

	selectedToggle = new GLToggle(0,"Export selected facets only");
	selectedToggle->SetBounds(10,150,170,18);
	Add(selectedToggle);

	exportButton = new GLButton(0,"Export");
	exportButton->SetBounds(10,hD-44,85,21);
	Add(exportButton);

	cancelButton = new GLButton(0,"Cancel");
	cancelButton->SetBounds(100,hD-44,85,21);
	Add(cancelButton);

	// Center dialog
	int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;
	SetBounds(xD,yD,wD,hD);

	RestoreDeviceObjects();

	geom = g;
	work = w;

}

void ExportDesorption::ProcessMessage(GLComponent *src,int message) {

	switch(message) {
	case MSG_BUTTON:

		if(src==cancelButton) {

			GLWindow::ProcessMessage(NULL,MSG_CLOSE);

		} else if (src==exportButton) {
			//do checks
			if (selectedToggle->GetState() && geom->GetNbSelectedFacets()==0) {
				GLMessageBox::Display("No facets selected","Nothing to export",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
			if( mode==2 && !eta0Field->GetNumber(&eta0) ) {
				GLMessageBox::Display("Invalid eta0 value","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
			if( mode==2 && !alphaField->GetNumber(&alpha) ) {
				GLMessageBox::Display("Invalid alpha value","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}

			if( mode==3 && !fileLoaded ) {
				GLMessageBox::Display("No conversion file loaded","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}

			bool meshExists = false;
			for (int i=0;i<geom->GetNbFacet() && !meshExists;i++) {
				if (geom->GetFacet(i)->hasMesh) meshExists = true;
			}
			if (!meshExists) {
				GLMessageBox::Display("There are no textures to convert to desorption","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}

			FILENAME *fn = GLFileBox::SaveFile(mApp->currentDir,NULL,"Save Desorption File","DES files\0*.des\0All files\0*.*\0",2);

			if( fn ) {
				
				try {
					work->ExportDesorption(fn->fullName,selectedToggle->GetState(),mode,eta0,alpha,conversionDistr);
					GLWindow::ProcessMessage(NULL,MSG_CLOSE);
				} catch (Error &e) {
					char errMsg[512];
					sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
					GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
				}

			}
		} else if (src==browseButton) {
			toggle1->SetState(false);
			toggle2->SetState(false);
			toggle3->SetState(true);
			//load file dialog
			FILENAME *convFile=GLFileBox::OpenFile(mApp->currentDir,NULL,"Open conversion file","CONV files\0*.conv\0All files\0*.*\0",2);
			if (!convFile) return;
			fileName=convFile->fullName;
			if (!fileName) return;
			//load file
			int nbPt=LoadConvFile(fileName);
			if (nbPt>0) {
			fileNameLabel->SetText(fileName);
			char tmp[256];
			sprintf(tmp,"%d points loaded.",nbPt);
			fileInfoLabel->SetText(tmp);
			} else {
				GLMessageBox::Display("Couldn't load conversion file","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
		}
		break;
	case MSG_TOGGLE:
		if (src!=selectedToggle) {
			toggle1->SetState(src==toggle1);
			toggle2->SetState(src==toggle2);
			toggle3->SetState(src==toggle3);
			if (src==toggle1) mode=1;
			if (src==toggle2) mode=2;
			if (src==toggle3) mode=3;
		}
		break;
	}

	GLWindow::ProcessMessage(src,message);
}

int ExportDesorption::LoadConvFile(char *fileName) {
	int nbPoints=0;
	FileReader *f = NULL;
	try {
		f=new FileReader(fileName);
		//count number of entries
		while (!f->IsEof()) {
			f->ReadDouble();
			f->ReadDouble();
			nbPoints++;
		}
		if (!(nbPoints>0)) throw Error("Invalid number of entries in file");
		//conversionDistr=new Distribution2D(nbPoints);
		
		f->SeekStart(); //restart from the beginning
		
		/*for (int i=0;i<nbPoints;i++) {
			conversionDistr->valuesX[i]=f->ReadDouble();
			conversionDistr->valuesY[i]=f->ReadDouble();
		}*/
		SAFE_DELETE(f);
	}  catch (Error &e) {
		SAFE_DELETE(f);
		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fileName);
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
	}
	fileLoaded=true;
	return nbPoints;
}