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
#include "TexturePlotter.h"
#include "Facet_shared.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLFileBox.h"
#include "SynRad.h"

extern SynRad *mApp;

static const char *fileFilters = "Text files\0*.txt";
static const int   nbFilter = sizeof(fileFilters) / (2*sizeof(char *));

TexturePlotter::TexturePlotter():GLWindow() {

	int wD = 500;
	int hD = 300;
	lastUpdate = 0.0f;
	strcpy(currentDir,".");

	SetTitle("Texture plotter");
	SetResizable(true);
	SetIconfiable(true);
	SetMinimumSize(wD,hD);

	mapList = new GLList(0);
	mapList->SetColumnLabelVisible(true);
	mapList->SetRowLabelVisible(true);
	mapList->SetAutoColumnLabel(true);
	mapList->SetAutoRowLabel(true);
	mapList->SetRowLabelMargin(20);
	mapList->SetGrid(true);
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
	autoSizeOnUpdate->SetState(true);
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

void TexturePlotter::SetBounds(int x,int y,int w,int h) {

	GLWindow::SetBounds(x,y,w,h);
	PlaceComponents();

}

void TexturePlotter::GetSelected() {

	if(!worker) return;

	Geometry *geom = worker->GetGeometry();
	selFacet = NULL;
	size_t i = 0;
	size_t nb = geom->GetNbFacet();
	while(!selFacet && i<nb) {
		if( geom->GetFacet(i)->selected ) selFacet = geom->GetFacet(i);
		if(!selFacet) i++;
	}

	char tmp[32];
	sprintf(tmp,"Texture plotter #%zd",i+1);
	SetTitle(tmp);

}

void TexturePlotter::Update(float appTime,bool force) {

	if(!IsVisible()) return;

	if(force) {
		UpdateTable();
		lastUpdate = appTime;
		return;
	}

	if( (appTime-lastUpdate>1.0f) ) {
		if(worker->isRunning) UpdateTable();
		lastUpdate = appTime;
	}

}

void TexturePlotter::UpdateTable() {

	maxValue=0.0f;
	GetSelected();
	if( !selFacet || !selFacet->cellPropertiesIds) {
		mapList->Clear();
		return;
	}

	//SHELEM *mesh = selFacet->mesh;
	if(selFacet->cellPropertiesIds) {

		char tmp[256];
		size_t w = selFacet->sh.texWidth;
		size_t h = selFacet->sh.texHeight;
		size_t textureSize=w*h*sizeof(TextureCell);
		size_t profile_memory=PROFILE_SIZE*sizeof(ProfileSlice);
		mapList->SetSize(w,h);
		mapList->SetAllColumnAlign(ALIGN_CENTER);

		int mode = viewCombo->GetSelectedIndex();

		switch(mode) {

		case 0: {// Cell area
			for(size_t i=0;i<w;i++) {
				for(size_t j=0;j<h;j++) {
					float val=selFacet->GetMeshArea(i+j*w);
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
					GlobalHitBuffer *shGHit = (GlobalHitBuffer *)buffer;
					size_t profSize = (selFacet->sh.isProfile)?profile_memory:0;
					TextureCell *texture = (TextureCell *)((BYTE *)buffer + (selFacet->sh.hitOffset + sizeof(FacetHitBuffer) + profSize));

					for(size_t i=0;i<w;i++) {
						for(size_t j=0;j<h;j++) {
							double val=texture[i+j*w].flux/worker->no_scans; //already divided by area
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

					GlobalHitBuffer *shGHit = (GlobalHitBuffer *)buffer;
					size_t profSize = (selFacet->sh.isProfile)?profile_memory:0;
					TextureCell *texture = (TextureCell *)((BYTE *)buffer + (selFacet->sh.hitOffset + sizeof(FacetHitBuffer) + profSize));

					for(size_t i=0;i<w;i++) {
						for(size_t j=0;j<h;j++) {
							double val=texture[i+j*w].power/worker->no_scans;

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
					GlobalHitBuffer *shGHit = (GlobalHitBuffer *)buffer;
					size_t profSize = (selFacet->sh.isProfile)?profile_memory:0;
					TextureCell *texture = (TextureCell *)((BYTE *)buffer + (selFacet->sh.hitOffset + sizeof(FacetHitBuffer) + profSize));
					float dCoef = 1.0f;
					//if( shGHit->mode == MC_MODE ) dCoef=1.0;//dCoef = (float)totalOutgassing / (float)shGHit->total.nbDesorbed;

					for(size_t i=0;i<w;i++) {
						for(size_t j=0;j<h;j++) {

							llong val=texture[i+j*w].count;
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

void TexturePlotter::Display(Worker *w) {

	worker = w;
	UpdateTable();
	SetVisible(true);

}

void TexturePlotter::Close() {
	worker = NULL;
	if(selFacet) selFacet->UnselectElem();
	mapList->Clear();
}

void TexturePlotter::SaveFile() {

	if(!selFacet) return;

	FILENAME *fn = GLFileBox::SaveFile(currentDir,NULL,"Save File",fileFilters,nbFilter);

	if( fn ) {

		size_t u,v,wu,wv;
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

		for(size_t i=u;i<u+wu;i++) {
			for(size_t j=v;j<v+wv;j++) {
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
			size_t u,v,wu,wv;
			mapList->SetSelectedCell(maxX,maxY);
			if( mapList->GetSelectionBox(&v,&u,&wv,&wu) )
				selFacet->SelectElem(u,v,wu,wv);
		}
		break;

	case MSG_LIST:
		if(src==mapList) {
			size_t u,v,wu,wv;
			if( mapList->GetSelectionBox(&v,&u,&wv,&wu) )
				selFacet->SelectElem(u,v,wu,wv);
		}
		break;

	case MSG_COMBO:
		if(src==viewCombo) {
			UpdateTable();
			maxButton->SetEnabled(true);
			//maxButton->SetEnabled(viewCombo->GetSelectedIndex()!=2);
		}
		break;

	}

	GLWindow::ProcessMessage(src,message);
}
