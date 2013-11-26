/*
File:        SpectrumPlotter.cpp
Description: Spectrum plotter window
Program:     SynRad
Author:      R. KERSEVAN / M SZAKACS
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

#include "SpectrumPlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Utils.h"
#include <math.h>
#include "Synrad.h"

extern GLApplication *theApp;

static const char*specMode[] = {"Flux (ph/sec/.1%)","Power (W/.1%)"};

SpectrumPlotter::SpectrumPlotter():GLWindow() {

	int wD = 750;
	int hD = 350;

	SetTitle("Spectrum plotter");
	SetIconfiable(TRUE);
	nbView = 0;
	worker = NULL;
	lastUpdate = 0.0f;

	chart = new GLChart(0);
	chart->SetBorder(BORDER_BEVEL_IN);
	chart->GetY1Axis()->SetGridVisible(TRUE);
	chart->GetXAxis()->SetGridVisible(TRUE);
	chart->GetY1Axis()->SetAutoScale(TRUE);
	chart->GetY2Axis()->SetAutoScale(TRUE);
	chart->GetY1Axis()->SetAnnotation(VALUE_ANNO);
	chart->GetXAxis()->SetAnnotation(VALUE_ANNO);
	chart->GetXAxis()->SetScale(LOG_SCALE); //logarithmic X scale
	chart->GetY1Axis()->SetScale(LOG_SCALE); //logarithmic Y scale
	Add(chart);
	delta=1.0;

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

	specCombo = new GLCombo(0);
	specCombo->SetEditable(TRUE);
	Add(specCombo);

	logToggle = new GLToggle(0,"Log Y scale");
	logToggle->SetCheck(TRUE);
	Add(logToggle);

	normToggle = new GLToggle(0,"Normalize");
	normToggle->SetCheck(TRUE);
	Add(normToggle);

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

void SpectrumPlotter::SetBounds(int x,int y,int w,int h) {
	chart->SetBounds(7,5,w-15,h-60);
	specCombo->SetBounds(7,h-45,117,19);
	selButton->SetBounds(130,h-45,80,19);
	addButton->SetBounds(215,h-45,80,19);
	removeButton->SetBounds(300,h-45,80,19);
	resetButton->SetBounds(385,h-45,80,19);
	logToggle->SetBounds(487,h-45,50,19);
	normToggle->SetBounds(537,h-45,105,19);
	dismissButton->SetBounds(w-100,h-45,90,19);
	GLWindow::SetBounds(x,y,w,h);
}

void SpectrumPlotter::Refresh() {

	if(!worker) return;

	Geometry *geom = worker->GetGeometry();
	int nb = geom->GetNbFacet();
	int nbSpec = 0;
	for(int i=0;i<nb;i++)
		if(geom->GetFacet(i)->sh.hasSpectrum) nbSpec+=2;
	specCombo->Clear();
	if(nbSpec) specCombo->SetSize(nbSpec);
	nbSpec=0;
	for(int i=0;i<nb;i++) {
		Facet *f = geom->GetFacet(i);
		if(f->sh.hasSpectrum) {
			char tmp[128];
			for (int mode=0;mode<2;mode++) { //Flux, power
				sprintf(tmp,"F#%d %s",i+1,specMode[mode]);
				specCombo->SetValueAt(nbSpec,tmp,i*2+mode);
				specCombo->SetSelectedIndex(0);
				nbSpec++;
			}
		}
	}

	refreshViews();

}

void SpectrumPlotter::SetScale() {
	if (worker) {
	if ((int)worker->regions.size()>0) {
		Region *traj=&(worker->regions[0]); //scale axis X
		if (traj->isLoaded) {
			chart->GetXAxis()->SetMinimum(traj->energy_low);
			chart->GetXAxis()->SetMaximum(traj->energy_hi);
			delta=(log(traj->energy_hi)-log(traj->energy_low))/SPECTRUM_SIZE;
		}
	}
	}
}

void SpectrumPlotter::Display(Worker *w) {

	worker = w;
	Refresh();
	SetVisible(TRUE);
	SetScale();
}

void SpectrumPlotter::Update(float appTime,BOOL force) {

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

void SpectrumPlotter::refreshViews() {

	// Lock during update
	BYTE *buffer = worker->GetHits();
	int normalize = normToggle->IsChecked();

	if(!buffer) return;

	Geometry *geom = worker->GetGeometry();	
	double no_scans;
	if (worker->nbTrajPoints==0 || worker->nbDesorption==0) no_scans=1.0;
	else no_scans=(double)worker->nbDesorption/(double)worker->nbTrajPoints;

	for(int i=0;i<nbView;i++) {

		GLDataView *v = views[i];
		if( v->userData1>=0 && v->userData1<geom->GetNbFacet()) {
			Facet *f = geom->GetFacet(v->userData1);
			int mode=v->userData2;
			v->Reset();

			int profileSize=(f->sh.isProfile)?(PROFILE_SIZE*(sizeof(llong)+2*sizeof(double))):0;
			int textureSize=(f->sh.isTextured)?(f->sh.texWidth*f->sh.texHeight*(2*sizeof(double)+sizeof(llong))):0;
			int directionSize=(f->sh.countDirection)?(f->sh.texWidth*f->sh.texHeight*sizeof(VHIT)):0;

			double *shSpectrum_fluxwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize 
				+ textureSize + directionSize));
			double *shSpectrum_powerwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize 
				+ textureSize + directionSize + SPECTRUM_SIZE*sizeof(double)));

			double max_flux,max_power;

			switch(normalize) {
			case 0: //no normalization
				if (mode==0) { //flux
					for(int j=0;j<PROFILE_SIZE;j++)
						v->Add(exp(log(chart->GetXAxis()->GetMinimum())+j*delta),shSpectrum_fluxwise[j]/no_scans,FALSE);
				} else if (mode==1) { //power
					for(int j=0;j<PROFILE_SIZE;j++)
						v->Add(exp(log(chart->GetXAxis()->GetMinimum())+j*delta),shSpectrum_powerwise[j]/no_scans,FALSE);
				}
				break;
			case 1: //normalize max. value to 1
				if (mode==0) { //flux
					max_flux=0.0;
					for(int j=0;j<PROFILE_SIZE;j++)
						if (shSpectrum_fluxwise[j]>max_flux) max_flux=shSpectrum_fluxwise[j];
					if (max_flux>0.0)
						for(int j=0;j<PROFILE_SIZE;j++)
							v->Add(exp(log(chart->GetXAxis()->GetMinimum())+j*delta),shSpectrum_fluxwise[j]/max_flux,FALSE);
				} else if (mode==1) { //power
					max_power=0.0;
					for(int j=0;j<PROFILE_SIZE;j++)
						if (shSpectrum_powerwise[j]>max_power) max_power=shSpectrum_powerwise[j];
					if (max_power>0.0)
						for(int j=0;j<PROFILE_SIZE;j++)
							v->Add(exp(log(chart->GetXAxis()->GetMinimum())+j*delta),shSpectrum_powerwise[j]/max_power,FALSE);
				}
				break;
			}
			v->CommitChange();
		}

	}

	worker->ReleaseHits();

}


void SpectrumPlotter::addView(int facet,int mode) {

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
		GLMessageBox::Display("Spectrum already plotted","Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	if(nbView<32) {
		Facet *f = geom->GetFacet(facet);
		GLDataView *v = new GLDataView();
		sprintf(tmp,"F#%d %s",facet+1,specMode[mode]);
		v->SetName(tmp);
		v->userData1 = facet;
		v->userData2 = mode;
		chart->GetY1Axis()->AddDataView(v);
		views[nbView] = v;
		nbView++;
	}

}

void SpectrumPlotter::remView(int facet,int mode) {

	Geometry *geom = worker->GetGeometry();

	BOOL found = FALSE;
	int i = 0; 
	while(i<nbView && !found) {
		found = (views[i]->userData1 == facet && views[i]->userData2 == mode);
		if(!found) i++;
	}
	if( !found ) {
		GLMessageBox::Display("Spectrum not plotted","Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	chart->GetY1Axis()->RemoveDataView(views[i]);
	SAFE_DELETE(views[i]);
	for(int j=i;j<nbView-1;j++) views[j] = views[j+1];
	nbView--;

}

void SpectrumPlotter::Reset() {

	chart->GetY1Axis()->ClearDataView();
	for(int i=0;i<nbView;i++) SAFE_DELETE(views[i]);
	nbView=0;

}

void SpectrumPlotter::ProcessMessage(GLComponent *src,int message) {
	Geometry *geom = worker->GetGeometry();
	SynRad *mApp = (SynRad *)theApp;
	switch(message) {
	case MSG_BUTTON:
		if(src==dismissButton) {
			SetVisible(FALSE);
		} else if(src==selButton) {
			int idx = specCombo->GetSelectedIndex();
			geom->UnSelectAll();
			int facetId=(int)((double)specCombo->GetUserValueAt(idx)/2.0);
			geom->GetFacet(facetId)->selected = TRUE;
			mApp->UpdateFacetParams(TRUE);
			geom->UpdateSelection();
			mApp->facetList->SetSelectedRow(facetId);
			mApp->facetList->ScrollToVisible(facetId,1,TRUE);
		} else if(src==addButton) {
			int idx = specCombo->GetSelectedIndex();
			int facetId=(int)((double)specCombo->GetUserValueAt(idx)/2.0);
			if(idx>=0) addView(facetId,specCombo->GetUserValueAt(idx)%2);
			refreshViews();
		} else if(src==removeButton) {
			int idx = specCombo->GetSelectedIndex();
			int facetId=(int)((double)specCombo->GetUserValueAt(idx)/2.0);
			if(idx>=0) remView(facetId,specCombo->GetUserValueAt(idx)%2);
			refreshViews();
		} else if(src==resetButton) {
			Reset();
		}
	case MSG_TOGGLE:
		if( src==normToggle ) {
			refreshViews();
		} else if (src==logToggle) {
			chart->GetY1Axis()->SetScale(logToggle->IsChecked()?LOG_SCALE:LINEAR_SCALE);
		}
	}

	GLWindow::ProcessMessage(src,message);

}

