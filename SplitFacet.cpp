/*
  File:        SplitFacet.cpp
  Description: Split facet by plane dialog
  Program:     MolFlow


  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#define PLANEEQ_MODE 0
#define FACET_MODE 1
#define VERTEX_MODE 2

#include "SplitFacet.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "SynRad.h"

extern SynRad *mApp;

SplitFacet::SplitFacet(Geometry *g,Worker *w):GLWindow() {
	
	int wD = 356;
	int hD = 276;
	planeDefPanel = new GLTitledPanel("Plane definition mode");
	planeDefPanel->SetBounds(12, 12, 326, 179);
	Add(planeDefPanel);
	label1 = new GLLabel("*X +");
	planeDefPanel->SetCompBounds(label1, 54, 45, 27, 13);
	planeDefPanel->Add(label1);

	eqmodeCheckbox = new GLToggle(0, "By equation");
	planeDefPanel->SetCompBounds(eqmodeCheckbox, 6, 19, 82, 17);
	planeDefPanel->Add(eqmodeCheckbox);

	XZplaneButton = new GLButton(0, "XZ plane");
	planeDefPanel->SetCompBounds(XZplaneButton, 168, 68, 75, 23);
	planeDefPanel->Add(XZplaneButton);

	YZplaneButton = new GLButton(0, "YZ plane");
	planeDefPanel->SetCompBounds(YZplaneButton, 87, 68, 75, 23);
	planeDefPanel->Add(YZplaneButton);

	XYplaneButton = new GLButton(0, "XY plane");
	planeDefPanel->SetCompBounds(XYplaneButton, 6, 68, 75, 23);
	planeDefPanel->Add(XYplaneButton);

	dTextbox = new GLTextField(0, "");
	planeDefPanel->SetCompBounds(dTextbox, 249, 42, 42, 20);
	planeDefPanel->Add(dTextbox);

	label4 = new GLLabel("= 0");
	planeDefPanel->SetCompBounds(label4, 297, 45, 22, 13);
	planeDefPanel->Add(label4);

	cTextbox = new GLTextField(0, "");
	planeDefPanel->SetCompBounds(cTextbox, 168, 42, 42, 20);
	planeDefPanel->Add(cTextbox);

	label3 = new GLLabel("*Z +");
	planeDefPanel->SetCompBounds(label3, 216, 45, 27, 13);
	planeDefPanel->Add(label3);

	bTextbox = new GLTextField(0, "");
	planeDefPanel->SetCompBounds(bTextbox, 87, 42, 42, 20);
	planeDefPanel->Add(bTextbox);

	label2 = new GLLabel("*Y +");
	planeDefPanel->SetCompBounds(label2, 135, 45, 27, 13);
	planeDefPanel->Add(label2);

	aTextbox = new GLTextField(0, "");
	planeDefPanel->SetCompBounds(aTextbox, 6, 42, 42, 20);
	planeDefPanel->Add(aTextbox);

	vertexModeCheckbox = new GLToggle(0, "By 3 selected vertex");
	planeDefPanel->SetCompBounds(vertexModeCheckbox, 6, 155, 122, 17);
	planeDefPanel->Add(vertexModeCheckbox);

	facetIdTextbox = new GLTextField(0, "");
	planeDefPanel->SetCompBounds(facetIdTextbox, 107, 111, 55, 20);
	planeDefPanel->Add(facetIdTextbox);

	facetmodeCheckbox = new GLToggle(0, "Plane of facet:");
	planeDefPanel->SetCompBounds(facetmodeCheckbox, 6, 114, 95, 17);
	planeDefPanel->Add(facetmodeCheckbox);

	resultLabel = new GLLabel("resultLabel");
	resultLabel->SetBounds(15, 194, 58, 13);
	Add(resultLabel);

	splitButton = new GLButton(0, "Split");
	splitButton->SetBounds(18, 221, 129, 23);
	Add(splitButton);

	undoButton = new GLButton(0, "Undo");
	undoButton->SetBounds(200, 221, 129, 23);
	Add(undoButton);

	SetTitle("SplitFacet");
	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);
  RestoreDeviceObjects();
	
  resultLabel->SetText("");
  undoButton->SetEnabled(FALSE);
  planeMode = PLANEEQ_MODE;
  geom = g;
  work = w;

}

SplitFacet::~SplitFacet() {
	ClearUndoFacets();
}

void SplitFacet::ClearUndoFacets() {
	//Destroy old undo facets
	for (auto delFacet : deletedFacetList)
		delete delFacet.f;
	deletedFacetList.clear();
	resultLabel->SetText("");
}

void SplitFacet::ProcessMessage(GLComponent *src,int message) {
  double a,b,c,d;
  int facetNum;

  switch(message) {
	// -------------------------------------------------------------
    case MSG_TOGGLE:
		if (src == eqmodeCheckbox) planeMode = PLANEEQ_MODE;
		else if (src == facetmodeCheckbox) planeMode = FACET_MODE;
		else if (src == vertexModeCheckbox) planeMode = VERTEX_MODE;
		EnableDisableControls(planeMode);	  
      break;

	case MSG_BUTTON:

	if (src == XYplaneButton) {
		aTextbox->SetText(0);
		bTextbox->SetText(0);
		cTextbox->SetText(1);
		dTextbox->SetText(0);
		planeMode = PLANEEQ_MODE;
		EnableDisableControls(planeMode);
		}
	else if (src == YZplaneButton) {
		aTextbox->SetText(1);
		bTextbox->SetText(0);
		cTextbox->SetText(0);
		dTextbox->SetText(0);
		planeMode = PLANEEQ_MODE;
		EnableDisableControls(planeMode);
	}
	else if (src == XZplaneButton) {
		aTextbox->SetText(0);
		bTextbox->SetText(1);
		cTextbox->SetText(0);
		dTextbox->SetText(0);
		planeMode = PLANEEQ_MODE;
		EnableDisableControls(planeMode);
	}
    else if(src==undoButton) {
		if (nbFacet == geom->GetNbFacet()) { //Assume no change since the split operation
			std::vector<size_t> newlyCreatedList;
			for (size_t index = (geom->GetNbFacet() - nbCreated);index < geom->GetNbFacet();index++) {
				newlyCreatedList.push_back(index);
			}
			geom->RemoveFacets(newlyCreatedList);
			geom->RestoreFacets(deletedFacetList,FALSE); //Restore to original position
		} else {
			int answer = GLMessageBox::Display("Geometry changed since split, restore to end without deleting the newly created facets?", "Split undo", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO);
			geom->RestoreFacets(deletedFacetList, TRUE); //Restore to end
		}
		deletedFacetList.clear();
		undoButton->SetEnabled(FALSE);
		resultLabel->SetText("");
		//Renumberformula
		work->Reload();
		mApp->UpdateModelParams();
		mApp->UpdateFacetlistSelected();
		mApp->UpdateViewers();
    } else if (src==splitButton) {
		if (geom->GetNbSelected()==0) {
			GLMessageBox::Display("No facets selected","Nothing to split",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
			//Calculate the plane
			VERTEX3D P0,N;
			int nbSelectedVertex;
			double nN2;
			int *vIdx = (int *)malloc(geom->GetNbVertex()*sizeof(int));
			memset(vIdx,0xFF,geom->GetNbVertex()*sizeof(int));
			switch (planeMode) {
			case FACET_MODE:
				if( !(facetIdTextbox->GetNumberInt(&facetNum))||facetNum<1||facetNum>geom->GetNbFacet() ) {
					GLMessageBox::Display("Invalid facet number","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				P0=*geom->GetVertex(geom->GetFacet(facetNum-1)->indices[0]);
				N=geom->GetFacet(facetNum-1)->sh.N;
				break;
			case VERTEX_MODE:
				if (geom->GetNbSelectedVertex()!=3) {
						        GLMessageBox::Display("Select exactly 3 vertices","Can't define plane",GLDLG_OK,GLDLG_ICONERROR);
								return;
				}
				nbSelectedVertex = 0;

				for(int i=0;i<geom->GetNbVertex()&&nbSelectedVertex<geom->GetNbSelectedVertex();i++ ) {
				//VERTEX3D *v = GetVertex(i);
					if( geom->GetVertex(i)->selected ) {
						vIdx[nbSelectedVertex] = i;
						nbSelectedVertex++;
				 }
				}

				VERTEX3D U2 = *geom->GetVertex(vIdx[0]) - *geom->GetVertex(vIdx[1]);
				Normalize(&U2); // Normalize U2

				VERTEX3D V2 = *geom->GetVertex(vIdx[0]) - *geom->GetVertex(vIdx[2]);
				Normalize(&V2); // Normalize V2

				VERTEX3D N2 = CrossProduct(V2,U2); //We have a normal vector
				nN2 = Norme(N2);
				if (nN2<1e-8) {
					GLMessageBox::Display("The 3 selected vertices are on a line.","Can't define plane",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				; // Normalize N2
				N=N2*(1.0/nN2);
				P0=*(geom->GetVertex(vIdx[0]));
				break;
			case PLANEEQ_MODE:
				if( !(aTextbox->GetNumber(&a)) ) {
					GLMessageBox::Display("Invalid A coefficient","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				if( !(bTextbox->GetNumber(&b)) ) {
					GLMessageBox::Display("Invalid B coefficient","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				if( !(cTextbox->GetNumber(&c)) ) {
					GLMessageBox::Display("Invalid C coefficient","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				if( !(dTextbox->GetNumber(&d)) ) {
					GLMessageBox::Display("Invalid D coefficient","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				if ((a==0.0)&&(b==0.0)&&(c==0.0)) {
					GLMessageBox::Display("A, B, C are all zero. That's not a plane.","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				N.x=a;N.y=b;N.z=c;
				P0.x=0.0;P0.y=0;P0.z=0;
				if (!a==0) P0.x=-d/a;
				else if (!b==0) P0.y=-d/b;
				else if (!c==0) P0.z=-d/c;
			break;
			default:
					GLMessageBox::Display("Select a plane definition mode.","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
			}
			SAFE_FREE(vIdx);
			if (mApp->AskToReset()) {
				GLProgress *prg = new GLProgress("Splitting facets", "Facet split");
				
				ClearUndoFacets();
				nbCreated = 0;
				deletedFacetList=geom->SplitSelectedFacets(P0, N, &nbCreated,prg);
				nbFacet = geom->GetNbFacet();
				SAFE_DELETE(prg);
				std::stringstream tmp;
				tmp << deletedFacetList.size() << " facets split, creating " << nbCreated <<" new.";
				resultLabel->SetText(tmp.str().c_str());
				if (deletedFacetList.size() > 0) undoButton->SetEnabled(TRUE);
				work->Reload();
				mApp->UpdateModelParams();
				mApp->UpdateFacetlistSelected();
				mApp->UpdateViewers();
	       		//GLWindowManager::FullRepaint();
			}
    }
    break;
	case MSG_TEXT_UPD:
		if (src == aTextbox || src == bTextbox || src == cTextbox || src == dTextbox) {
			planeMode = PLANEEQ_MODE;
			EnableDisableControls(planeMode);
		}
		else if (src == facetIdTextbox) {
			planeMode = FACET_MODE;
			EnableDisableControls(planeMode);
		}
		break;
  }

  GLWindow::ProcessMessage(src,message);
}

void SplitFacet::EnableDisableControls(int mode) {
	eqmodeCheckbox->SetState(mode == PLANEEQ_MODE);
	/*aTextbox->SetEditable(mode == PLANEEQ_MODE);
	bTextbox->SetEditable(mode == PLANEEQ_MODE);
	cTextbox->SetEditable(mode == PLANEEQ_MODE);
	dTextbox->SetEditable(mode == PLANEEQ_MODE);*/

	facetmodeCheckbox->SetState(mode == FACET_MODE);
	//facetIdTextbox->SetEditable(mode == FACET_MODE);

	vertexModeCheckbox->SetState(mode == VERTEX_MODE);
}