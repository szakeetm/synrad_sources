/*
File:        ScaleFacet.cpp
Description: Mirror facet to plane dialog
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

#define XYZMODE 0
#define VERTEXMODE 1
#define FACETMODE 2

#define UNIFORMMODE 0
#define DISTORTMODE 1

#include "ScaleFacet.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "SynRad.h"

extern SynRad *theApp;

ScaleFacet::ScaleFacet(Geometry *g,Worker *w):GLWindow() {

	int wD = 290;
	int hD = 220;
	invariantMode=XYZMODE;
	scaleMode=UNIFORMMODE;

	SetTitle("Scale selected vertices");

	iPanel = new GLTitledPanel("Invariant point definiton mode");
	iPanel->SetBounds(5,3,wD-10,97);
	Add(iPanel);

	l1 = new GLToggle(0,"");
	l1->SetBounds(10,20,20,18);
	l1->SetCheck(TRUE);
	iPanel->Add(l1);

	GLLabel *xLabel = new GLLabel("X=");
	xLabel->SetBounds(30,20,20,18);
	iPanel->Add(xLabel);

	xText = new GLTextField(0,"0");
	xText->SetBounds(45,20,50,18);
	//xText->SetEditable(FALSE);
	iPanel->Add(xText);

	GLLabel *yLabel = new GLLabel("Y=");
	yLabel->SetBounds(100,20,20,18);
	iPanel->Add(yLabel);

	yText = new GLTextField(0,"0");
	yText->SetBounds(115,20,50,18);
	//yText->SetEditable(FALSE);
	iPanel->Add(yText);

	GLLabel* zLabel = new GLLabel("Z=");
	zLabel->SetBounds(170,20,20,18);
	iPanel->Add(zLabel);

	zText = new GLTextField(0,"0");
	zText->SetBounds(185,20,50,18);
	//zText->SetEditable(FALSE);
	iPanel->Add(zText);

	l2 = new GLToggle(0,"Selected Vertex");
	l2->SetBounds(10,45,100,18);
	iPanel->Add(l2);

	l3 = new GLToggle(0,"Center of Facet #");
	l3->SetBounds(10,70,100,18);
	iPanel->Add(l3);

	facetNumber = new GLTextField(0,"0");
	facetNumber->SetBounds(115,70,50,18);
	facetNumber->SetEditable(FALSE);
	iPanel->Add(facetNumber);

	sPanel = new GLTitledPanel("Scale factor");
	sPanel->SetBounds(5,101,wD-10,65);
	Add(sPanel);

	uniform = new GLToggle(0,"Uniform");
	uniform->SetBounds(10,115,100,18);
	uniform->SetCheck(TRUE);
	Add(uniform);

	factorNumber = new GLTextField(0,"1");
	factorNumber->SetBounds(115,115,40,18);
	//factorNumber->SetEditable(FALSE);
	Add(factorNumber);

	distort = new GLToggle(0,"Distorted");
	distort->SetBounds(10,140,100,18);
	Add(distort);

	GLLabel* x2Label = new GLLabel("X:");
	x2Label->SetBounds(103,141,20,18);
	Add(x2Label);

	factorNumberX = new GLTextField(0,"1");
	factorNumberX->SetBounds(115,140,40,18);
	factorNumberX->SetEditable(FALSE);
	Add(factorNumberX);

	GLLabel* y2Label = new GLLabel("Y:");
	y2Label->SetBounds(163,141,20,18);
	Add(y2Label);

	factorNumberY = new GLTextField(0,"1");
	factorNumberY->SetBounds(175,140,40,18);
	factorNumberY->SetEditable(FALSE);
	Add(factorNumberY);

	GLLabel* z2Label = new GLLabel("Z:");
	z2Label->SetBounds(223,141,20,18);
	Add(z2Label);

	factorNumberZ = new GLTextField(0,"1");
	factorNumberZ->SetBounds(235,140,40,18);
	factorNumberZ->SetEditable(FALSE);
	Add(factorNumberZ);

	scaleButton = new GLButton(0,"Scale");
	scaleButton->SetBounds(5,hD-44,85,21);
	Add(scaleButton);

	copyButton = new GLButton(0,"Copy");
	copyButton->SetBounds(95,hD-44,85,21);
	Add(copyButton);

	cancelButton = new GLButton(0,"Dismiss");
	cancelButton->SetBounds(185,hD-44,85,21);
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

void ScaleFacet::ProcessMessage(GLComponent *src,int message) {
	SynRad *mApp = (SynRad *)theApp;
	double x,y,z,factor,factorX,factorY,factorZ;
	int facetNum;

	switch(message) {
		// -------------------------------------------------------------
	case MSG_TOGGLE:
		UpdateToggle(src);
		break;

	case MSG_BUTTON:

		if(src==cancelButton) {

			GLWindow::ProcessMessage(NULL,MSG_CLOSE);

		} else if (src==scaleButton || src==copyButton) {
			if (geom->GetNbSelected()==0) {
				GLMessageBox::Display("No facets selected","Nothing to scale",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}

			VERTEX3D invariant;
			BOOL found;

			switch (invariantMode) {
			case XYZMODE:
				if( !(xText->GetNumber(&x))) {
					GLMessageBox::Display("Invalid X coordinate","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				if( !(yText->GetNumber(&y))) {
					GLMessageBox::Display("Invalid Y coordinate","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				if( !(zText->GetNumber(&z))) {
					GLMessageBox::Display("Invalid Z coordinate","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				invariant.x=x;
				invariant.y=y;
				invariant.z=z;
				break;
			case FACETMODE:
				if( !(facetNumber->GetNumberInt(&facetNum))||facetNum<1||facetNum>geom->GetNbFacet() ) {
					GLMessageBox::Display("Invalid facet number","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				invariant=geom->GetFacet(facetNum-1)->sh.center;
				break;
			case VERTEXMODE:
				if( !(geom->GetNbSelectedVertex()==1) ) {
					GLMessageBox::Display("Select exactly one vertex","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				found=FALSE;
				for (int i=0;!found && i<geom->GetNbVertex();i++) {
					if (geom->GetVertex(i)->selected)
						invariant=*(geom->GetVertex(i));
				}
				break;
			default:
				GLMessageBox::Display("Select a plane definition mode.","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}

			if (scaleMode==UNIFORMMODE) {
				if( !(factorNumber->GetNumber(&factor))) {
					GLMessageBox::Display("Invalid scale factor number","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
			} else {
				if( !(factorNumberX->GetNumber(&factorX))) {
					GLMessageBox::Display("Invalid X scale factor number","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}

				if( !(factorNumberY->GetNumber(&factorY))) {
					GLMessageBox::Display("Invalid Y scale factor number","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}

				if( !(factorNumberZ->GetNumber(&factorZ))) {
					GLMessageBox::Display("Invalid Z scale factor number","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
			}
			if (mApp->AskToReset()) {
				if (scaleMode==UNIFORMMODE) factorX=factorY=factorZ=factor;
				geom->ScaleSelectedFacets(invariant,factorX,factorY,factorZ,src==copyButton,work);
				theApp->UpdateModelParams();
				try { work->Reload(); } catch(Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
				}  
				theApp->UpdateFacetlistSelected();
				mApp->UpdateViewers();
				//GLWindowManager::FullRepaint();
			}
		}
		break;
	}

	GLWindow::ProcessMessage(src,message);
}

void ScaleFacet::UpdateToggle(GLComponent *src) {
	
	if (src==l1 || src==l2 || src==l3) {

		l1->SetCheck(src==l1);
		l2->SetCheck(src==l2);
		l3->SetCheck(src==l3);

		facetNumber->SetEditable(src==l3);
		xText->SetEditable(src==l1);
		yText->SetEditable(src==l1);
		zText->SetEditable(src==l1);

		if (src==l1) invariantMode=XYZMODE;
		if (src==l2) invariantMode=VERTEXMODE;
		if (src==l3) invariantMode=FACETMODE;
	} else {
		uniform->SetCheck(src==uniform);
		distort->SetCheck(src==distort);

		factorNumber->SetEditable(src==uniform);
		factorNumberX->SetEditable(src==distort);
		factorNumberY->SetEditable(src==distort);
		factorNumberZ->SetEditable(src==distort);

		if (src==uniform) scaleMode=UNIFORMMODE;
		if (src==distort) scaleMode=DISTORTMODE;
	}
}