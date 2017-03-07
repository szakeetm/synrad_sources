/*
File:        ProfilePlotter.cpp
Description: Profile plotter window
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

#include "ProfilePlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Facet.h"

#include <math.h>
#include "Synrad.h"

extern SynRad *mApp;

static const char*profType[] = {"None","\201","\202","Angle"};
static const char*profMode[] = {"MC Hits","SR Flux/cm\262","SR Power/cm\262"};

ProfilePlotter::ProfilePlotter():GLWindow() {

	int wD = 750;
	int hD = 400;

	SetTitle("Profile plotter");
	SetIconfiable(TRUE);
	nbView = 0;
	worker = NULL;
	lastUpdate = 0.0f;

	nbColors = 8;
	colors[0] = new GLCColor(); colors[0]->r = 255; colors[0]->g = 000; colors[0]->b = 055; //red
	colors[1] = new GLCColor(); colors[1]->r = 000; colors[1]->g = 000; colors[1]->b = 255; //blue
	colors[2] = new GLCColor(); colors[2]->r = 000; colors[2]->g = 204; colors[2]->b = 051; //green
	colors[3] = new GLCColor(); colors[3]->r = 000; colors[3]->g = 000; colors[3]->b = 000; //black
	colors[4] = new GLCColor(); colors[4]->r = 255; colors[4]->g = 153; colors[4]->b = 051; //orange
	colors[5] = new GLCColor(); colors[5]->r = 153; colors[5]->g = 204; colors[5]->b = 255; //light blue
	colors[6] = new GLCColor(); colors[6]->r = 153; colors[6]->g = 000; colors[6]->b = 102; //violet
	colors[7] = new GLCColor(); colors[7]->r = 255; colors[7]->g = 230; colors[7]->b = 005; //yellow

	chart = new GLChart(0);
	chart->SetBorder(BORDER_BEVEL_IN);
	chart->GetY1Axis()->SetGridVisible(TRUE);
	chart->GetXAxis()->SetGridVisible(TRUE);
	chart->GetY1Axis()->SetAutoScale(TRUE);
	chart->GetY2Axis()->SetAutoScale(TRUE);
	chart->GetY1Axis()->SetAnnotation(VALUE_ANNO);
	chart->GetXAxis()->SetAnnotation(VALUE_ANNO);
	Add(chart);

	dismissButton = new GLButton(0,"Dismiss");
	Add(dismissButton);

	selButton = new GLButton(0,"Show Facet");
	Add(selButton);

	addButton = new GLButton(0,"Add curve");
	Add(addButton);

	removeButton = new GLButton(0,"Remove curve");
	Add(removeButton);

	resetButton = new GLButton(0,"Remove all");
	Add(resetButton);

	profCombo = new GLCombo(0);
	profCombo->SetEditable(TRUE);
	Add(profCombo);

	logToggle = new GLToggle(0,"Log Y scale");
	logToggle->SetState(FALSE);
	Add(logToggle);

	normToggle = new GLToggle(0,"Normalize");
	normToggle->SetState(FALSE);
	Add(normToggle);

	formulaText = new GLTextField(0,"");
	formulaText->SetEditable(TRUE);
	Add(formulaText);

	formulaBtn = new GLButton(0,"-> Plot");
	Add(formulaBtn);

	// Center dialog
	int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;
	SetBounds(xD,yD,wD,hD);
	SetResizable(TRUE);
	SetMinimumSize(wD,220);

	RestoreDeviceObjects();

}

void ProfilePlotter::SetBounds(int x,int y,int w,int h) {

	chart->SetBounds(7,5,w-15,h-85);
	profCombo->SetBounds(7,h-70,117,19);
	selButton->SetBounds(130,h-70,80,19);
	addButton->SetBounds(215,h-70,80,19);
	removeButton->SetBounds(300,h-70,80,19);
	resetButton->SetBounds(385,h-70,80,19);
	logToggle->SetBounds(487,h-70,50,19);
	normToggle->SetBounds(537,h-70,105,19);
	formulaText->SetBounds(7,h-45,525,19);
	formulaBtn->SetBounds(537,h-45,80,19);
	dismissButton->SetBounds(w-100,h-45,90,19);

	GLWindow::SetBounds(x,y,w,h);

}

void ProfilePlotter::Refresh() {

	if(!worker) return;

	//Rebuild selection combo box
	Geometry *geom = worker->GetGeometry();
	int nb = geom->GetNbFacet();
	int nbProf = 0;
	for(int i=0;i<nb;i++)
		if(geom->GetFacet(i)->sh.isProfile) nbProf+=3;
	profCombo->Clear(); profCombo->SetSelectedIndex(0);
	if(nbProf) profCombo->SetSize(nbProf);
	nbProf=0;
	for (int i = 0; i < nb; i++) {
		Facet *f = geom->GetFacet(i);
		if (f->sh.isProfile) {
			char tmp[128];
			for (int mode = 0; mode < 3; mode++) { //MC hits, flux, power
				sprintf(tmp, "F#%d %s %s", i + 1, profType[f->sh.profileType], profMode[mode]);
				profCombo->SetValueAt(nbProf, tmp, i * 3 + mode);
				profCombo->SetSelectedIndex(0);
				nbProf++;
			}
		}
	}
	//Remove profiles that aren't present anymore
	for (int v = 0; v < nbView; v++)
		if (views[v]->userData1 >= geom->GetNbFacet() || !geom->GetFacet(views[v]->userData1)->sh.isProfile) {
			chart->GetY1Axis()->RemoveDataView(views[v]);
			SAFE_DELETE(views[v]);
			for (int j = v; j < nbView - 1; j++) views[j] = views[j + 1];
			nbView--;
		}

	//Update values
	refreshViews();

}

void ProfilePlotter::Display(Worker *w) {

	worker = w;
	Refresh();
	SetVisible(TRUE);

}

void ProfilePlotter::Update(float appTime,BOOL force) {

	if(!IsVisible() || IsIconic()) return;  

	if(force) {
		refreshViews();
		lastUpdate = appTime;
		return;
	}

	if( (appTime-lastUpdate>1.0f || force) && nbView ) {
		if(worker->running) refreshViews();
		lastUpdate = appTime;
	}

}

void ProfilePlotter::plot() {

	GLParser *parser = new GLParser();
	parser->SetExpression( formulaText->GetText() );
	if( !parser->Parse() ) {
		GLMessageBox::Display(parser->GetErrorMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	int nbVar = parser->GetNbVariable();
	if( nbVar==0 ) {
		GLMessageBox::Display("Variable 'x' not found","Error",GLDLG_OK,GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	if( nbVar>1 ) {
		GLMessageBox::Display("Too much variables or unknown constant","Error",GLDLG_OK,GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	VLIST *var = parser->GetVariableAt(0);
	if(_stricmp(var->name,"x")!=0) {
		GLMessageBox::Display("Variable 'x' not found","Error",GLDLG_OK,GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	Geometry *geom = worker->GetGeometry();
	GLDataView *v;

	// Check that view is not already added
	BOOL found = FALSE;
	int i = 0; 
	while(i<nbView && !found) {
		found = (views[i]->userData1 == -1);
		if(!found) i++;
	}

	if( found ) {
		v = views[i];
		v->SetName(formulaText->GetText());
		v->Reset();
	} else {
		if(nbView<32) {
			v = new GLDataView();
			v->SetName(formulaText->GetText());
			v->userData1 = -1;
			chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}
	}

	// Plot
	for(int i=0;i<1000;i++) {
		double x=(double)i;
		double y;
		var->value = x;
		parser->Evaluate(&y);
		v->Add(x,y,FALSE);
	}
	v->CommitChange();

	delete parser;

}

void ProfilePlotter::refreshViews() {

	// Lock during update
	BYTE *buffer = worker->GetHits();
	int normalize = normToggle->GetState();

	if(!buffer) return;

	Geometry *geom = worker->GetGeometry();
	SHGHITS *gHits = (SHGHITS *)buffer;
	double nbAbs = (double)gHits->total.nbAbsorbed;
	double nbDes = (double)gHits->total.nbDesorbed;
	double nbHit = (double)gHits->total.nbHit;
	
	for(int i=0;i<nbView;i++) {

		GLDataView *v = views[i];
		if( v->userData1>=0 && v->userData1<geom->GetNbFacet()) {
			Facet *f = geom->GetFacet(v->userData1);
			int mode=v->userData2;
			v->Reset();
			llong   *profilePtr_MC = (llong *)(buffer + f->sh.hitOffset + sizeof(SHHITS));
			double  *profilePtr_flux = (double *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + PROFILE_SIZE*sizeof(llong));
			double  *profilePtr_power= (double *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + PROFILE_SIZE*(sizeof(llong)+sizeof(double)));
			llong max_MC;
			double max_flux,max_power;
			double elemArea=f->sh.area/PROFILE_SIZE;

			switch(normalize) {
			case 0: //linear no normalization
				if (mode==0) { //mc hits
					for(int j=0;j<PROFILE_SIZE;j++)
						v->Add((double)j,(double)profilePtr_MC[j],FALSE);
				} else if (mode==1) { //flux
					for(int j=0;j<PROFILE_SIZE;j++)
						v->Add((double)j,profilePtr_flux[j]/worker->no_scans/elemArea,FALSE);
				} else if (mode==2) { //power
					for(int j=0;j<PROFILE_SIZE;j++)
						v->Add((double)j,profilePtr_power[j]/worker->no_scans/elemArea,FALSE);
				}
				break;

			case 1: //normalize max. value to 1
				if (mode==0) { //mc hits
					max_MC=0;
					for(int j=0;j<PROFILE_SIZE;j++)
						if (profilePtr_MC[j]>max_MC) max_MC=profilePtr_MC[j];
					if (max_MC>0)
						for(int j=0;j<PROFILE_SIZE;j++)
							v->Add((double)j,((double)profilePtr_MC[j]/(double)max_MC),FALSE);
				} else if (mode==1) { //flux
					max_flux=0.0;
					for(int j=0;j<PROFILE_SIZE;j++)
						if (profilePtr_flux[j]>max_flux) max_flux=profilePtr_flux[j];
					if (max_flux>0.0)
						for(int j=0;j<PROFILE_SIZE;j++)
							v->Add((double)j,profilePtr_flux[j]/max_flux,FALSE);
				} else if (mode==2) { //power
					max_power=0.0;
					for(int j=0;j<PROFILE_SIZE;j++)
						if (profilePtr_power[j]>max_power) max_power=profilePtr_power[j];
					if (max_power>0.0)
						for(int j=0;j<PROFILE_SIZE;j++)
							v->Add((double)j,profilePtr_power[j]/max_power,FALSE);
				}
				break;
			}
			v->CommitChange();
		} else {
			if( v->userData1==-2 && nbDes!=0.0 ) {

				// Volatile profile
				v->Reset();
				int nb = geom->GetNbFacet();
				for(int j=0;j<nb;j++) {
					Facet *f = geom->GetFacet(j);
					if( f->sh.isVolatile ) {
						SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
						double z = geom->GetVertex(f->indices[0])->z;
						v->Add(z,(double)(fCount->nbAbsorbed)/nbDes,FALSE);
					}
				}
				// Last
				Facet *f = geom->GetFacet(28);
				SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
				double fnbAbs = (double)fCount->nbAbsorbed;
				v->Add(1000.0,fnbAbs/nbDes,FALSE);
				v->CommitChange();
			}
		}

		}

		worker->ReleaseHits();

	}


	void ProfilePlotter::addView(int facet,int mode) {

		char tmp[128];
		Geometry *geom = worker->GetGeometry();

		// Check that view is not already added
		BOOL found = FALSE;
		int i = 0; 
		while(i<nbView && !found) {
			found = (views[i]->userData1 == facet && views[i]->userData2 == mode);
			if(!found) i++;
		}
		if( found ) {
			GLMessageBox::Display("Profile already plotted","Error",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
		if(nbView<32) {
			Facet *f = geom->GetFacet(facet);
			GLDataView *v = new GLDataView();
			sprintf(tmp,"F#%d %s %s",facet+1,profType[f->sh.profileType],profMode[mode]);
			v->SetName(tmp);
			v->SetColor(*colors[nbView%nbColors]);
			v->SetMarkerColor(*colors[nbView%nbColors]);
			v->SetLineWidth(2);
			v->userData1 = facet;
			v->userData2 = mode;
			chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}

	}

	void ProfilePlotter::remView(int facet,int mode) {

		Geometry *geom = worker->GetGeometry();

		BOOL found = FALSE;
		int i = 0; 
		while(i<nbView && !found) {
			found = (views[i]->userData1 == facet && views[i]->userData2 == mode);
			if(!found) i++;
		}
		if( !found ) {
			GLMessageBox::Display("Profile not plotted","Error",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
		chart->GetY1Axis()->RemoveDataView(views[i]);
		SAFE_DELETE(views[i]);
		for(int j=i;j<nbView-1;j++) views[j] = views[j+1];
		nbView--;

	}

	void ProfilePlotter::Reset() {

		chart->GetY1Axis()->ClearDataView();
		for(int i=0;i<nbView;i++) SAFE_DELETE(views[i]);
		nbView=0;

	}

	void ProfilePlotter::ProcessMessage(GLComponent *src,int message) {
		Geometry *geom = worker->GetGeometry();
		
		switch(message) {
		case MSG_BUTTON:
			if(src==dismissButton) {
				SetVisible(FALSE);
			} else if(src==selButton) {
				int idx = profCombo->GetSelectedIndex();
				if (idx>=0) {
					geom->UnselectAll();
					int facetId=(int)((double)profCombo->GetUserValueAt(idx)/3.0);
					geom->GetFacet(facetId)->selected = TRUE;
					mApp->UpdateFacetParams(TRUE);
					geom->UpdateSelection();
					mApp->facetList->SetSelectedRow(facetId);
					mApp->facetList->ScrollToVisible(facetId,1,TRUE);
				}
			} else if(src==addButton) {
				int idx = profCombo->GetSelectedIndex();
				int facetId=(int)((double)profCombo->GetUserValueAt(idx)/3.0);
				if(idx>=0) addView(facetId,profCombo->GetUserValueAt(idx)%3);
				refreshViews();
			} else if(src==removeButton) {
				int idx = profCombo->GetSelectedIndex();
				int facetId=(int)((double)profCombo->GetUserValueAt(idx)/3.0);
				if(idx>=0) remView(facetId,profCombo->GetUserValueAt(idx)%3);
				refreshViews();
			} else if(src==resetButton) {
				Reset();
			} else if(src==formulaBtn) {
				plot();
			}
		case MSG_TOGGLE:
			if( src==normToggle ) {
				refreshViews();
			} else if (src==logToggle) {
				chart->GetY1Axis()->SetScale(logToggle->GetState()?LOG_SCALE:LINEAR_SCALE);
			}
		}

		GLWindow::ProcessMessage(src,message);

	}

void ProfilePlotter::SetViews(std::vector<int> views,int mode) {
	Reset();
	for (int view : views)
		addView(view,mode);
	Refresh();
}

std::vector<int> ProfilePlotter::GetViews() {
	std::vector<int>v;
	v.reserve(nbView);
	for (int i = 0; i < nbView; i++)
		v.push_back(views[i]->userData1);
	return v;
}

BOOL ProfilePlotter::IsLogScaled() {
	return chart->GetY1Axis()->GetScale();
}
void ProfilePlotter::SetLogScaled(BOOL logScale){
	chart->GetY1Axis()->SetScale(logScale);
	logToggle->SetState(logScale);
}