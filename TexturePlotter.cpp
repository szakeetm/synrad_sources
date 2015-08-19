/*
File:        TexturePlotter.cpp
Description: Texture plotter dialog
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

#include "TexturePlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLFileBox.h"
#include "SynRad.h"

extern SynRad *mApp;

static const char *fileFilters = "Text files\0*.txt";
static const int   nbFilter = sizeof(fileFilters) / (2*sizeof(char *));

// --------------------------------------------------------------------

TexturePlotter::TexturePlotter():GLWindow() {

	int wD = 500;
	int hD = 300;
	lastUpdate = 0.0f;
	strcpy(currentDir,".");

	SetTitle("Texture plotter");
	SetResizable(TRUE);
	SetIconfiable(TRUE);
	SetMinimumSize(wD,hD);

	mapList = new GLList(0);
	mapList->SetColumnLabelVisible(TRUE);
	mapList->SetRowLabelVisible(TRUE);
	mapList->SetAutoColumnLabel(TRUE);
	mapList->SetAutoRowLabel(TRUE);
	mapList->SetRowLabelMargin(20);
	mapList->SetGrid(TRUE);
	mapList->SetSelectionMode(BOX_CELL);
	mapList->SetCornerLabel("\202\\\201");
	Add(mapList);

	viewLabel = new GLLabel("View");
	Add(viewLabel);
	viewCombo = new GLCombo(0);
	viewCombo->SetSize(4);
	viewCombo->SetValueAt(0,"Elem.area (cm\262)");
	viewCombo->SetValueAt(1,"SR Flux/cm\262");
	viewCombo->SetValueAt(2,"SR Power/cm\262");
	viewCombo->SetValueAt(3,"MC Hits");



	viewCombo->SetSelectedIndex(1);
	Add(viewCombo);

	saveButton = new GLButton(0,"Save");
	Add(saveButton);

	sizeButton = new GLButton(0,"Auto size");
	Add(sizeButton);

	maxButton = new GLButton(0,"Find Max.");
	Add(maxButton);

	cancelButton = new GLButton(0,"Dismiss");
	Add(cancelButton);

	autoSizeOnUpdate = new GLToggle(0,"Autosize on every update");
	autoSizeOnUpdate->SetState(TRUE);
	Add(autoSizeOnUpdate);

	// Center dialog
	int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;
	SetBounds(xD,yD,wD,hD);

	RestoreDeviceObjects();

	worker = NULL;

}

// --------------------------------------------------------------------

void TexturePlotter::PlaceComponents() {

	mapList->SetBounds(5,5,width-15,height-80);
	saveButton->SetBounds(10,height-70,70,19);
	sizeButton->SetBounds(10,height-45,70,19);
	autoSizeOnUpdate->SetBounds(90,height-45,120,19);
	maxButton->SetBounds(90,height-70,70,19);
	viewLabel->SetBounds(330,height-70,30,19);
	viewCombo->SetBounds(360,height-70,120,19);
	cancelButton->SetBounds(width-90,height-45,80,19);

}

// -----------------------------------------------------------------

void TexturePlotter::SetBounds(int x,int y,int w,int h) {

	GLWindow::SetBounds(x,y,w,h);
	PlaceComponents();

}

// --------------------------------------------------------------------

void TexturePlotter::GetSelected() {

	if(!worker) return;

	Geometry *geom = worker->GetGeometry();
	selFacet = NULL;
	int i = 0;
	int nb = geom->GetNbFacet();
	while(!selFacet && i<nb) {
		if( geom->GetFacet(i)->selected ) selFacet = geom->GetFacet(i);
		if(!selFacet) i++;
	}

	char tmp[32];
	sprintf(tmp,"Texture plotter #%d",i+1);
	SetTitle(tmp);

}

// --------------------------------------------------------------------

void TexturePlotter::Update(float appTime,BOOL force) {

	if(!IsVisible()) return;

	if(force) {
		UpdateTable();
		lastUpdate = appTime;
		return;
	}

	if( (appTime-lastUpdate>1.0f) ) {
		if(worker->running) UpdateTable();
		lastUpdate = appTime;
	}

}

// --------------------------------------------------------------------

void TexturePlotter::UpdateTable() {

	maxValue=0.0f;
	GetSelected();
	if( !selFacet || !selFacet->mesh) {
		mapList->Clear();
		return;
	}

	SHELEM *mesh = selFacet->mesh;
	if( mesh ) {

		char tmp[256];
		int w = selFacet->sh.texWidth;
		int h = selFacet->sh.texHeight;
		int textureSize_double=w*h*sizeof(double);
		int textureSize_llong=w*h*sizeof(llong);
		int profile_memory=PROFILE_SIZE*(2*sizeof(double)+sizeof(llong));
		mapList->SetSize(w,h);
		mapList->SetColumnAlign(ALIGN_CENTER);


		int mode = viewCombo->GetSelectedIndex();

		switch(mode) {

		case 0: {// Cell area
			for(int i=0;i<w;i++) {
				for(int j=0;j<h;j++) {
					float val=selFacet->mesh[i+j*w].area;
					sprintf(tmp,"%g",val);
					if (val>maxValue) {
						maxValue=val;
						maxX=i;maxY=j;
					}
					mapList->SetValueAt(i,j,tmp);
				}
			}
			break; }

		case 1:  {// Flux

			// Lock during update
			BYTE *buffer = worker->GetHits();

			try {
				if(buffer) {
					SHGHITS *shGHit = (SHGHITS *)buffer;
					int profSize = (selFacet->sh.isProfile)?profile_memory:0;
					double *hits_flux = (double *)((BYTE *)buffer + (selFacet->sh.hitOffset + sizeof(SHHITS) + profSize + textureSize_llong));




					for(int i=0;i<w;i++) {
						for(int j=0;j<h;j++) {
							double val=hits_flux[i+j*w]/worker->no_scans; //already divided by area
							if (val>maxValue) {
								maxValue=(float)val;
								maxX=i;maxY=j;
							}
							sprintf(tmp,"%g",val);
							mapList->SetValueAt(i,j,tmp);
						}
					}
					worker->ReleaseHits();
				}
			} catch (...) { //incorrect hits reference
				worker->ReleaseHits();
			}
			break;}




		case 2:  {// Power




			// Lock during update
			BYTE *buffer = worker->GetHits();
			try{
				if(buffer) {


					SHGHITS *shGHit = (SHGHITS *)buffer;
					int profSize = (selFacet->sh.isProfile)?profile_memory:0;
					double *hits_power = (double *)((BYTE *)buffer + (selFacet->sh.hitOffset + sizeof(SHHITS) + profSize + textureSize_llong + textureSize_double));

					for(int i=0;i<w;i++) {
						for(int j=0;j<h;j++) {
							double val=hits_power[i+j*w]/worker->no_scans;

							if (val>maxValue) {
								maxValue=(float)val;
								maxX=i;maxY=j;
							}
							sprintf(tmp,"%g",val);
							mapList->SetValueAt(i,j,tmp);
						}
					}
					worker->ReleaseHits();

				}
			} catch (...) {
				worker->ReleaseHits();

			}

			break;}

		case 3: {// MC Hits
			// Lock during update
			BYTE *buffer = worker->GetHits();
			try{
				if(buffer) {
					SHGHITS *shGHit = (SHGHITS *)buffer;
					int profSize = (selFacet->sh.isProfile)?profile_memory:0;
					llong *hits_MC = (llong *)((BYTE *)buffer + (selFacet->sh.hitOffset + sizeof(SHHITS) + profSize));
					float dCoef = 1.0f;
					//if( shGHit->mode == MC_MODE ) dCoef=1.0;//dCoef = (float)totalOutgassing / (float)shGHit->total.nbDesorbed;

					for(int i=0;i<w;i++) {
						for(int j=0;j<h;j++) {


							llong val=hits_MC[i+j*w];
							if (val>maxValue) {
								maxValue=(float)val;
								maxX=i;maxY=j;
							}
							sprintf(tmp,"%g",(double)val);



							mapList->SetValueAt(i,j,tmp);
						}
					}
					worker->ReleaseHits();
				}
			} catch (...) {
				worker->ReleaseHits();

			}

			break;}
		}

	}

	if (autoSizeOnUpdate->GetState()) mapList->AutoSizeColumn();
}

// --------------------------------------------------------------------

void TexturePlotter::Display(Worker *w) {

	worker = w;
	UpdateTable();
	SetVisible(TRUE);

}

// --------------------------------------------------------------------

void TexturePlotter::Close() {
	worker = NULL;
	if(selFacet) selFacet->UnselectElem();
	mapList->Clear();
}

// --------------------------------------------------------------------

void TexturePlotter::SaveFile() {

	if(!selFacet) return;

	FILENAME *fn = GLFileBox::SaveFile(currentDir,NULL,"Save File",fileFilters,nbFilter);

	if( fn ) {

		int u,v,wu,wv;
		if( !mapList->GetSelectionBox(&u,&v,&wu,&wv) ) {
			u=0;
			v=0;
			wu = mapList->GetNbRow();
			wv = mapList->GetNbColumn();
		}

		// Save tab separated text
		FILE *f = fopen(fn->fullName,"w");

		if( f==NULL ) {
			char errMsg[512];
			sprintf(errMsg,"Cannot open file\nFile:%s",fn->fullName);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}

		for(int i=u;i<u+wu;i++) {
			for(int j=v;j<v+wv;j++) {
				char *str = mapList->GetValueAt(j,i);
				if( str ) fprintf(f,"%s",str);
				if( j<v+wv-1 ) 
					fprintf(f,"\t");
			}
			fprintf(f,"\r\n");
		}
		fclose(f);

	}

}

// --------------------------------------------------------------------

void TexturePlotter::ProcessMessage(GLComponent *src,int message) {

	switch(message) {

	case MSG_CLOSE:
		Close();
		break;

	case MSG_BUTTON:
		if(src==cancelButton) {
			Close();
			GLWindow::ProcessMessage(NULL,MSG_CLOSE);
		} else if (src==sizeButton) {
			mapList->AutoSizeColumn();
		} else if (src==saveButton) {
			SaveFile();
		} else if (src==maxButton) {
			int u,v,wu,wv;
			mapList->SetSelectedCell(maxX,maxY);
			if( mapList->GetSelectionBox(&v,&u,&wv,&wu) )
				selFacet->SelectElem(u,v,wu,wv);
		}
		break;

	case MSG_LIST:
		if(src==mapList) {
			int u,v,wu,wv;
			if( mapList->GetSelectionBox(&v,&u,&wv,&wu) )
				selFacet->SelectElem(u,v,wu,wv);
		}
		break;

	case MSG_COMBO:
		if(src==viewCombo) {
			UpdateTable();
			maxButton->SetEnabled(TRUE);
			//maxButton->SetEnabled(viewCombo->GetSelectedIndex()!=2);
		}
		break;

	}

	GLWindow::ProcessMessage(src,message);
}
