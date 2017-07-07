/*
File:        FacetDetails.cpp
Description: Facet details window
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

#include "FacetDetails.h"
#include "Facet.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/MathTools.h" //FormatMemory
#include "GLApp/GLMessageBox.h"

#include "SynRad.h"

extern SynRad *mApp;

typedef struct {

	char *name;
	int   width;
	int   align;

} COLUMN;

static COLUMN allColumn[] = {
	{"#"             , 40 , ALIGN_CENTER} ,
	{"Sticking"      , 80 , ALIGN_CENTER} ,
	{"Opacity"       , 80 , ALIGN_CENTER} ,
	{"Structure"     , 80 , ALIGN_CENTER} ,
	{"Link"          , 40 , ALIGN_CENTER} ,
	{"Reflection"    , 80 , ALIGN_CENTER} ,
	{"Roughness (nm)", 80 , ALIGN_CENTER } ,
	{"Autocorr.l (nm)", 80 , ALIGN_CENTER } ,
	{"2 Sided"       , 60 , ALIGN_CENTER} ,
	{"Vertex"        , 80 , ALIGN_CENTER} ,

	{"Area"          , 80 , ALIGN_CENTER} ,
	{"Facet 2D Box"  , 100, ALIGN_CENTER} ,
	{"Texture (u,v)" , 120, ALIGN_CENTER} ,
	{"Sample/cm"   , 80 , ALIGN_CENTER} ,
	{"Count"         , 80 , ALIGN_CENTER} ,
	{"Memory"        , 80 , ALIGN_CENTER} ,
	{"Planarity"      , 80 , ALIGN_CENTER} ,
	{"Profile"       , 80 , ALIGN_CENTER} ,
	{"Spectrum"      , 80 , ALIGN_CENTER} ,
	{"Hits"        , 80 , ALIGN_CENTER} ,

	{"Abs."        , 80 , ALIGN_CENTER} ,
	{"Flux"        , 80 , ALIGN_CENTER} ,
	{"Power"        , 80 , ALIGN_CENTER} ,
};

static const char *profStr[] = {
	"None",
	"Pressure (u)",
	"Pressure (v)",
	"Angular"
};

static const char *ynStr[] = {
	"No",
	"Yes"
};


// -----------------------------------------------------------------

FacetDetails::FacetDetails():GLWindow() {

	int wD = 502;
	int hD = 400;
	worker = NULL;

	SetTitle("Facets details");
	SetIconfiable(true);
	SetResizable(true);
	SetMinimumSize(502,200);

	checkAllButton = new GLButton(0,"Check All");
	Add(checkAllButton);
	uncheckAllButton = new GLButton(0,"Uncheck All");
	Add(uncheckAllButton);
	updateButton = new GLButton(0,"Update");
	Add(updateButton);
	dismissButton = new GLButton(0,"Dismiss");
	Add(dismissButton);

	facetListD = new GLList(0);
	facetListD->SetColumnLabelVisible(true);
	facetListD->SetGrid(true);
	Add(facetListD);

	sPanel = new GLTitledPanel("Show column");
	sPanel->SetClosable(true);
	Add(sPanel);

	show[0] = new GLToggle(0,"#");
	show[0]->SetState(true);  // Always visible (not displayed)

	show[1] = new GLToggle(1,"Sticking");
	show[1]->SetState(true);
	sPanel->Add(show[1]);
	show[2] = new GLToggle(2,"Opacity");
	show[2]->SetState(true);
	sPanel->Add(show[2]);
	show[3] = new GLToggle(3,"Structure");
	show[3]->SetState(true);
	sPanel->Add(show[3]);
	show[4] = new GLToggle(4,"Link");
	show[4]->SetState(true);
	sPanel->Add(show[4]);
	show[5] = new GLToggle(5,"Reflection");
	show[5]->SetState(true);
	sPanel->Add(show[5]);
	show[6] = new GLToggle(6,"Roughness");
	show[6]->SetState(true);
	sPanel->Add(show[6]);
	show[7] = new GLToggle(7, "AutoCorr");
	show[7]->SetState(true);
	sPanel->Add(show[7]);
	show[8] = new GLToggle(8,"2 Sided");
	show[8]->SetState(true);
	sPanel->Add(show[8]);
	show[9] = new GLToggle(9,"Vertex nb");
	show[9]->SetState(true);
	sPanel->Add(show[9]);
	show[10] = new GLToggle(10,"Area");
	show[10]->SetState(true);
	sPanel->Add(show[10]);

	show[11] = new GLToggle(11,"2D Box");
	show[11]->SetState(true);
	sPanel->Add(show[11]);
	show[12] = new GLToggle(12,"Texture UV");
	show[12]->SetState(true);
	sPanel->Add(show[12]);

	show[13] = new GLToggle(13,"Sample/cm");
	show[13]->SetState(true);
	sPanel->Add(show[13]);
	show[14] = new GLToggle(14,"Count mode");
	show[14]->SetState(true);
	sPanel->Add(show[14]);
	show[15] = new GLToggle(14,"Memory");
	show[15]->SetState(true);
	sPanel->Add(show[15]);
	show[16] = new GLToggle(16,"Planarity");
	show[16]->SetState(true);
	sPanel->Add(show[16]);
	show[17] = new GLToggle(17,"Profile");
	show[17]->SetState(true);
	sPanel->Add(show[17]);
	show[18] = new GLToggle(18,"Spectrum");
	show[18]->SetState(true);
	sPanel->Add(show[18]);
	show[19] = new GLToggle(19,"Hits");
	show[19]->SetState(true);
	sPanel->Add(show[19]);
	show[20] = new GLToggle(20,"Abs.");
	show[20]->SetState(true);
	sPanel->Add(show[20]);

	show[21] = new GLToggle(21,"Flux");
	show[21]->SetState(true);
	sPanel->Add(show[21]);
	show[22] = new GLToggle(22,"Power");
	show[22]->SetState(true);
	sPanel->Add(show[22]);

	// Center dialog
	int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;
	SetBounds(xD,yD,wD,hD);

	RestoreDeviceObjects();

}

// -----------------------------------------------------------------

void FacetDetails::PlaceComponents() {

	// Show toggle panel
	int nbW = (width-20)/80;
	int nbL = (NB_FDCOLUMN-1)/nbW + (((NB_FDCOLUMN-1)%nbW)?1:0);
	int hp;
	if(!sPanel->IsClosed())
		hp = 20*(nbL+1);
	else
		hp = 20;
	sPanel->SetBounds(5,height-(hp+52),width-10,hp);
	for(int i=0;i<NB_FDCOLUMN;i++)
		sPanel->SetCompBounds(show[i],5+80*((i-1)%nbW),18+20*((i-1)/nbW),85,19);

	facetListD->SetBounds(5,5,width-10,height-(62+hp));

	checkAllButton->SetBounds(5,height-45,90,19);
	uncheckAllButton->SetBounds(100,height-45,90,19);
	updateButton->SetBounds(195,height-45,90,19);
	dismissButton->SetBounds(width-100,height-45,90,19);

}

// -----------------------------------------------------------------

void FacetDetails::SetBounds(int x,int y,int w,int h) {

	GLWindow::SetBounds(x,y,w,h);
	PlaceComponents();

}

// -----------------------------------------------------------------

char *FacetDetails::GetCountStr(Facet *f) {
	static char ret[128];
	strcpy(ret,"");
	if(f->sh.countAbs) if(strlen(ret)==0) strcat(ret,"ABS"); else strcat(ret,"+ABS");
	if(f->sh.countRefl) if(strlen(ret)==0) strcat(ret,"REFL"); else strcat(ret,"+REFL");
	if(f->sh.countTrans) if(strlen(ret)==0) strcat(ret,"TRANS"); else strcat(ret,"+TRANS");
	return ret;
}

// -----------------------------------------------------------------

char *FacetDetails::FormatCell(int idx,Facet *f,int mode) {

	Worker *worker=&(mApp->worker);     

	static char ret[256];
	strcpy(ret,"");

	switch(mode) {
	case 0:
		sprintf(ret,"%d",idx+1);
		break;
	case 1:
		sprintf(ret,"%g",f->sh.sticking);
		break;
	case 2:
		sprintf(ret,"%g",f->sh.opacity);
		break;
	case 3:
		sprintf(ret,"%d",f->sh.superIdx);
		break;
	case 4:
		sprintf(ret,"%d",f->sh.superDest);
		break;
	case 5:
		if (f->sh.reflectType == 0) sprintf(ret, "Diffuse");
		else if (f->sh.reflectType == 1) sprintf(ret, "Mirror");
		else sprintf(ret,"%s",worker->materials[f->sh.reflectType-10].name.c_str());
		break;
	case 6:
		sprintf(ret,"%g",f->sh.rmsRoughness*1E9);
		break;
	case 7:
		sprintf(ret,"%g",f->sh.autoCorrLength*1E9);
		break;
	case 8:
		sprintf(ret,"%s",ynStr[f->sh.is2sided]);      
		break;
	case 9:
		sprintf(ret,"%d",f->sh.nbIndex);
		break;
	case 10:
		sprintf(ret,"%g",f->sh.area);
		break;

	case 11:
		sprintf(ret,"%g x %g",f->sh.U.Norme(),f->sh.V.Norme());
		break;
	case 12:
		if( f->sh.isTextured ) {
			sprintf(ret,"%dx%d (%g x %g)",f->sh.texWidth,f->sh.texHeight,f->sh.texWidthD,f->sh.texHeightD);
		} else {
			sprintf(ret,"None");
		}
		break;
	case 13:
		sprintf(ret,"%g",f->tRatio);
		break;
	case 14:
		sprintf(ret,"%s",GetCountStr(f));
		break;
	case 15:
		sprintf(ret,"%s",FormatMemory(f->GetTexRamSize()));
		break;
	case 16:
		sprintf(ret,"%f",f->err);
		break;
	case 17:
		sprintf(ret,"%s",profStr[f->sh.profileType]);
		break;
	case 18:
		sprintf(ret,f->sh.hasSpectrum==1?"Yes":"No");
		break;
	case 19:
		sprintf(ret,"%d",f->counterCache.nbHit);
		break;
	case 20:
		sprintf(ret,"%d",f->counterCache.nbAbsorbed);
		break;

	case 21:
		sprintf(ret,"%g",f->counterCache.fluxAbs/worker->no_scans);
		break;
	case 22:
		sprintf(ret,"%g",f->counterCache.powerAbs/worker->no_scans);
		break;
	}

	return ret;

}

// -----------------------------------------------------------------

void FacetDetails::UpdateTable() {

	Geometry *s = worker->GetGeometry();
	int nbFacet = s->GetNbFacet();
	int nbS = s->GetNbSelectedFacets();
	static char ret[256];
	strcpy(ret,"");

	/*
	//SUM Counters
	double sumArea=0;
	DWORD sumMemory=0;
	*/

	char *tmpName[NB_FDCOLUMN];
	int  tmpWidth[NB_FDCOLUMN];
	int  tmpAlign[NB_FDCOLUMN];


	int  nbCol = 1;
	tmpName[0]  = allColumn[0].name;
	tmpWidth[0] = allColumn[0].width;
	tmpAlign[0] = allColumn[0].align;
	shown[0] = 0;

	for(int i=1;i<NB_FDCOLUMN;i++) {
		if(show[i]->GetState()) {
			tmpName[nbCol]  = allColumn[i].name;
			tmpWidth[nbCol] = allColumn[i].width;
			tmpAlign[nbCol] = allColumn[i].align;
			shown[nbCol] = i;
			nbCol++;
		}
	}

	facetListD->SetSize(nbCol,nbS);
	facetListD->SetColumnWidths(tmpWidth);
	facetListD->SetColumnLabels(tmpName);
	facetListD->SetColumnAligns(tmpAlign);

	nbS = 0;
	for(int i=0;i<nbFacet;i++) {
		Facet *f = s->GetFacet(i);
		if(f->selected) {
			for(int j=0;j<nbCol;j++)
				facetListD->SetValueAt(j,nbS,FormatCell(i,f,shown[j]));

			/*
			sumArea+=f->sh.area;
			sumMemory+=f->GetTexRamSize();
			*/
			nbS++;
		}

	}
	/*
	//SUM values
	for(int j=0;j<nbCol;j++) {
	switch (shown[j]) {
	case 0:
	facetListD->SetValueAt(j,nbS,"SUM");
	break;
	case 9:
	sprintf(ret,"%g",sumArea);
	facetListD->SetValueAt(j,nbS,ret);
	break;
	case 14:
	sprintf(ret,"%s",FormatMemory(sumMemory));
	facetListD->SetValueAt(j,nbS,ret);
	break;
	default:
	facetListD->SetValueAt(j,nbS,"");
	}

	}
	*/
}

// -----------------------------------------------------------------

void FacetDetails::Update() {

	if(!worker) return;
	if(!IsVisible()) return;

	Geometry *s = worker->GetGeometry();
	int nbS = s->GetNbSelectedFacets();

	if(nbS==0) {
		facetListD->Clear();
		return;
	}

	UpdateTable();

}

// -----------------------------------------------------------------

void FacetDetails::Display(Worker *w) {

	worker = w;
	SetVisible(true);
	Update();

}

// -----------------------------------------------------------------

void FacetDetails::ProcessMessage(GLComponent *src,int message) {

	switch(message) {

	case MSG_BUTTON:
		if(src==dismissButton) {
			SetVisible(false);
		} else if (src==checkAllButton) {
			for(int i=0;i<NB_FDCOLUMN;i++) show[i]->SetState(true);
			UpdateTable();
		} else if (src==uncheckAllButton) {
			for(int i=0;i<NB_FDCOLUMN;i++) show[i]->SetState(false);
			UpdateTable();
		}	else if (src==updateButton) {
			UpdateTable();
		}
		break;

	case MSG_TOGGLE:
		UpdateTable();
		break;

	case MSG_LIST_COL:
		if( src==facetListD ) {
			// Save column width
			int c = facetListD->GetDraggedCol();
			allColumn[shown[c]].width = facetListD->GetColWidth(c);
		}
		break;

	case MSG_PANELR:
		PlaceComponents();
		break;

	}

	GLWindow::ProcessMessage(src,message);

}

