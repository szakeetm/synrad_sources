/*
File:        FacetMesh.cpp
Description: Facet mesh configuration dialog
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

#include <math.h>

#include "FacetMesh.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Utils.h" 
#include "SynRad.h"

//-----------------------------------------------------------------------------

FacetMesh::FacetMesh():GLWindow() {

	int wD = 275;
	int hD = 365;

	SetTitle("Facet Mesh");

	iPanel = new GLTitledPanel("Facet Dimension");
	iPanel->SetBounds(5,10,wD-10,45);
	Add(iPanel);

	GLLabel *l1 = new GLLabel("\201 length");
	l1->SetBounds(10,30,50,18);
	Add(l1);

	uLength = new GLTextField(0,"");
	uLength->SetBounds(60,30,70,18);
	uLength->SetEditable(FALSE);
	Add(uLength);

	GLLabel *l2 = new GLLabel("\202 length");
	l2->SetBounds(140,30,50,18);
	Add(l2);

	vLength = new GLTextField(0,"");
	vLength->SetBounds(190,30,70,18);
	vLength->SetEditable(FALSE);
	Add(vLength);

	GLTitledPanel *aPanel = new GLTitledPanel("Mesh properties");
	aPanel->SetBounds(5,60,wD-10,130);
	Add(aPanel);

	enableBtn = new GLToggle(0,"Enable");
	enableBtn->SetBounds(10,80,55,18);
	enableBtn->SetState(FALSE);
	Add(enableBtn);

	boundaryBtn = new GLToggle(0,"Boundary correction");
	boundaryBtn->SetBounds(10,100,100,18);
	boundaryBtn->SetEnabled(FALSE);
	boundaryBtn->SetTextColor(110,110,110);
	boundaryBtn->SetState(TRUE);
	Add(boundaryBtn);

	recordAbsBtn = new GLToggle(0,"Count absorption");
	recordAbsBtn->SetBounds(10,150,100,18);
	recordAbsBtn->SetState(FALSE);
	Add(recordAbsBtn);

	recordReflBtn = new GLToggle(0,"Count reflection");
	recordReflBtn->SetBounds(120,130,110,18);
	recordReflBtn->SetState(FALSE);
	Add(recordReflBtn);

	recordTransBtn = new GLToggle(0,"Count transparent pass");
	recordTransBtn->SetBounds(120,150,110,18);
	recordTransBtn->SetState(FALSE);
	Add(recordTransBtn);

	recordDirBtn = new GLToggle(0,"Record direction");
	recordDirBtn->SetBounds(120,170,110,18);
	recordDirBtn->SetState(FALSE);
	Add(recordDirBtn);

	GLLabel *l5 = new GLLabel("Resolution (Sample/cm)");
	l5->SetBounds(140,80,110,18);
	Add(l5);

	resolutionText = new GLTextField(0,"");
	resolutionText->SetBounds(140,100,50,18);
	Add(resolutionText);

	vPanel = new GLTitledPanel("View Settings");
	vPanel->SetBounds(5,195,wD-10,45);
	Add(vPanel);

	showTexture = new GLToggle(0,"Show texture");
	showTexture->SetBounds(10,215,55,18);
	showTexture->SetState(TRUE);
	Add(showTexture);

	showVolume = new GLToggle(0,"Show volume");
	showVolume->SetBounds(100,215,55,18);
	showVolume->SetState(TRUE);
	showVolume->SetVisible(TRUE); //not working yet
	Add(showVolume);

	quickApply = new GLButton(0,"Apply View");  //Apply View Settings without stopping the simulation
	quickApply->SetBounds(190,215,72,19);
	Add(quickApply);

	GLTitledPanel *mPanel = new GLTitledPanel("Memory/Cell");
	mPanel->SetBounds(5,245,wD-10,72);
	Add(mPanel);

	GLLabel *l7 = new GLLabel("Memory");
	l7->SetBounds(10,265,70,18);
	Add(l7);

	ramText = new GLTextField(0,"");
	ramText->SetBounds(80,265,100,18);
	Add(ramText);

	GLLabel *l8 = new GLLabel("Cells");
	l8->SetBounds(10,290,70,18);
	Add(l8);

	cellText = new GLTextField(0,"");
	cellText->SetBounds(80,290,100,18);
	Add(cellText);

	applyButton = new GLButton(0,"Apply mesh");
	applyButton->SetBounds(wD-200,hD-43,95,21);
	Add(applyButton);

	cancelButton = new GLButton(0,"Cancel");
	cancelButton->SetBounds(wD-100,hD-43,95,21);
	Add(cancelButton);

	// Center dialog
	int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;
	SetBounds(xD,yD,wD,hD);

	RestoreDeviceObjects();

}



void FacetMesh::UpdateSize() {

	char tmp[64];

	if( enableBtn->GetState() ) {

		llong ram = 0;
		llong cell = 0;
		int nbFacet = geom->GetNbFacet();

			for(int i=0;i<nbFacet;i++) {
				Facet *f = geom->GetFacet(i);
				cell += (llong)f->GetNbCell();
				ram += (llong)f->GetTexRamSize();
			}

		
		ramText->SetText(FormatMemoryLL(ram));
		sprintf(tmp,"%d",(int)cell);
		cellText->SetText(tmp);

	} else {

		ramText->SetText("0 bytes");
		cellText->SetText("0");

	}

}

void FacetMesh::UpdateSizeForRatio() {

	double ratio;
	char tmp[64];
	BOOL boundMap = boundaryBtn->GetState();
	BOOL recordDir = recordDirBtn->GetState();

	if( !enableBtn->GetState() ) {
		ramText->SetText(FormatMemory(0));
		cellText->SetText("0");
		return;
	}

	if( sscanf(resolutionText->GetText(),"%lf",&ratio)==0 ) {
		ramText->SetText("");
		cellText->SetText("");
		return;
	}

	llong ram = 0;
	llong cell = 0;
	int nbFacet = geom->GetNbFacet();

		for(int i=0;i<nbFacet;i++) {
			Facet *f = geom->GetFacet(i);
			if(f->selected) {
				cell += (llong)f->GetNbCellForRatio(ratio);
				ram += (llong)f->GetTexRamSizeForRatio(ratio,boundMap,recordDir);
			} else {
				cell += (llong)f->GetNbCell();
				ram += (llong)f->GetTexRamSize();
			}
		}

	

	ramText->SetText(FormatMemoryLL(ram));
	sprintf(tmp,"%d",(int)cell);
	cellText->SetText(tmp);

}

//-----------------------------------------------------------------------------

void FacetMesh::EditFacet(Worker *w) {

	char tmp[128];
	double maxU=0.0;
	double maxV=0.0;
	double minU=1.0e100;
	double minV=1.0e100;

	worker = w;
	geom   = w->GetGeometry();

	int nbS=0;
	int nbF=geom->GetNbFacet();
	int sel=0;
	BOOL allEnabled = TRUE;
	BOOL allBound = TRUE;
	BOOL ratioE = TRUE;
	BOOL allCountAbs = TRUE;
	BOOL allCountRefl = TRUE;
	BOOL allCountTrans = TRUE;
	BOOL allCountDir = TRUE;
	BOOL allTexVisible = TRUE;
	BOOL allVolVisible = TRUE;
	double tRatio = 1e100;

	for(int i=0;i<nbF;i++) {

		Facet *f = geom->GetFacet(i);
		if( f->selected ) {
			double nU = Norme(&(f->sh.U));
			double nV = Norme(&(f->sh.V));
			maxU = MAX(maxU,nU);
			maxV = MAX(maxV,nV);
			minU = MIN(minU,nU);
			minV = MIN(minV,nV);
			sel = i;
			allEnabled = allEnabled && f->sh.isTextured;
			allBound = allBound && (f->mesh!=NULL);
			allCountAbs = allCountAbs && f->sh.countAbs;
			allCountRefl = allCountRefl && f->sh.countRefl;
			allCountTrans = allCountTrans && f->sh.countTrans;
			allCountDir = allCountDir && f->sh.countDirection;
			allTexVisible = allTexVisible && f->textureVisible;
			allVolVisible = allVolVisible && f->volumeVisible;
			if( tRatio == 1e100 ) tRatio = f->tRatio;
			ratioE = ratioE & (tRatio == f->tRatio);
			nbS++;
		}

	}

	if( nbS==1 ) {
		Facet *f = geom->GetFacet(sel);
		sprintf(tmp,"Facet Info (#%d)",sel+1);
		iPanel->SetTitle(tmp);
		sprintf(tmp,"%g",maxU);
		uLength->SetText(tmp);
		sprintf(tmp,"%g",maxV);
		vLength->SetText(tmp);
	} else {
		sprintf(tmp,"Facet Info (%d selected)",nbS);
		iPanel->SetTitle(tmp);
		sprintf(tmp,"%g (MAX)",maxU);
		uLength->SetText(tmp);
		sprintf(tmp,"%g (MAX)",maxV);
		vLength->SetText(tmp);
	}

	enableBtn->SetState(allEnabled);
	//boundaryBtn->SetState(allBound);
	boundaryBtn->SetState(TRUE);
	recordAbsBtn->SetState(allCountAbs);
	recordReflBtn->SetState(allCountRefl);
	recordTransBtn->SetState(allCountTrans);
	recordDirBtn->SetState(allCountDir);
	showTexture->SetState(allTexVisible);
	showVolume->SetState(allVolVisible);

	if( allEnabled && ratioE ) {
		sprintf(tmp,"%g",tRatio);
		resolutionText->SetText(tmp);
	} else {
		resolutionText->SetText("...");
	}

	UpdateSize();
	DoModal();

}


//-----------------------------------------------------------------------------

BOOL FacetMesh::Apply() {
	extern GLApplication *theApp;
	SynRad *mApp = (SynRad *)theApp;
	if (!mApp->AskToReset(worker)) return FALSE;
	BOOL boundMap = boundaryBtn->GetState();
	double nbSelected = (double)geom->GetNbSelected();
	double nbPerformed = 0.0;

	if( enableBtn->GetState() ) {

		// Check counting mode
		if( !recordAbsBtn->GetState() && 
			!recordReflBtn->GetState() && !recordTransBtn->GetState() && 
			!recordDirBtn->GetState() ) {
				GLMessageBox::Display("Please select counting mode","Error",GLDLG_OK,GLDLG_ICONERROR);
				return FALSE;
		}

		// Auto resolution
		double ratio;
		if( sscanf(resolutionText->GetText(),"%lf",&ratio)==0 ) {
			GLMessageBox::Display("Invalid number format for sample/cm","Error",GLDLG_OK,GLDLG_ICONERROR);
			return FALSE;
		}

		progressDlg = new GLProgress("Applying mesh settings","Please wait");
		progressDlg->SetVisible(TRUE);
		progressDlg->SetProgress(0.0);
		int count=0;
		for(int i=0;i<geom->GetNbFacet();i++) {
			Facet *f = geom->GetFacet(i);
			if( f->selected ) {
				f->sh.countAbs = recordAbsBtn->GetState();
				f->sh.countRefl = recordReflBtn->GetState();
				f->sh.countTrans = recordTransBtn->GetState();
				f->sh.countDirection = recordDirBtn->GetState();
				f->textureVisible = showTexture->GetState();
				f->volumeVisible = showVolume->GetState();
				try {
					geom->SetFacetTexture(i,ratio,boundMap);
				} catch (Error &e) {
					char errMsg[512];
					sprintf(errMsg,"Error setting textures:\n%s",e.GetMsg());
					GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
					progressDlg->SetVisible(FALSE);
					SAFE_DELETE(progressDlg);
					return FALSE;
				}
				nbPerformed+=1.0;
				progressDlg->SetProgress(nbPerformed/nbSelected);
			}

		}

	} else {
		// Disable texture
		progressDlg = new GLProgress("Applying mesh settings","Please wait");
		progressDlg->SetVisible(TRUE);
		progressDlg->SetProgress(0.0);


		for(int i=0;i<geom->GetNbFacet();i++) {

			Facet *f = geom->GetFacet(i);
			if( f->selected ) {
				geom->SetFacetTexture(i,0.0,FALSE);
				f->textureVisible = showTexture->GetState();
				f->volumeVisible = showVolume->GetState();
				nbPerformed+=1.0;
				progressDlg->SetProgress(nbPerformed/nbSelected);
			}
		}

	}

	// Send to sub process
	try {
		worker->Reload();

	} catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
	}
	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
	return TRUE;

}

//-----------------------------------------------------------------------------
void FacetMesh::QuickApply() {
	//Apply view settings without stopping the simulation


	double nbSelected = (double)geom->GetNbSelected();
	double nbPerformed = 0.0;

	for(int i=0;i<geom->GetNbFacet();i++) {

		Facet *f = geom->GetFacet(i);
		if( f->selected ) {

			f->textureVisible = showTexture->GetState();
			f->volumeVisible = showVolume->GetState();

			nbPerformed+=1.0;
			progressDlg->SetProgress(nbPerformed/nbSelected);
		}

	}
	geom->RebuildLists();
}

//-----------------------------------------------------------------------------

void FacetMesh::UpdateToggle(GLComponent *src) {

	if(src==enableBtn) {
		//boundaryBtn->SetState(enableBtn->GetState());
	} else if(src==recordAbsBtn ) {
		enableBtn->SetState(TRUE);
		boundaryBtn->SetState(TRUE);
	} else if(src==recordReflBtn ) {
		enableBtn->SetState(TRUE);
		boundaryBtn->SetState(TRUE);
	} else if(src==recordTransBtn ) {
		enableBtn->SetState(TRUE);
		boundaryBtn->SetState(TRUE);
	} else if(src==recordDirBtn ) {
		enableBtn->SetState(TRUE);
		boundaryBtn->SetState(TRUE);
	}
	UpdateSizeForRatio();
}

//-----------------------------------------------------------------------------

void FacetMesh::ProcessMessage(GLComponent *src,int message) {

	switch(message) {

		// -------------------------------------------------------------
	case MSG_BUTTON:
		if(src==cancelButton) {

			GLWindow::ProcessMessage(NULL,MSG_CLOSE);

		} else if (src==applyButton) {


			//if (worker->running) worker->Stop_Public();
			if( Apply() )
				GLWindow::ProcessMessage(NULL,MSG_CLOSE);

		} else if (src==quickApply) {

			progressDlg = new GLProgress("Applying view settings","Please wait");
			progressDlg->SetVisible(TRUE);
			progressDlg->SetProgress(0.5);

			QuickApply();
			GLWindow::ProcessMessage(NULL,MSG_CLOSE);

			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);

		} 
		break;

		// -------------------------------------------------------------
	case MSG_TEXT_UPD:
		enableBtn->SetState(TRUE);
		UpdateSizeForRatio();
		break;

		// -------------------------------------------------------------
	case MSG_TOGGLE:
		UpdateToggle(src);
		break;
	}

	GLWindow::ProcessMessage(src,message);
}

