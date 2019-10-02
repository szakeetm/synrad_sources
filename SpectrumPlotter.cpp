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
#include "SpectrumPlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/MathTools.h"
#include "Facet_shared.h"

#include <math.h>
#include "Synrad.h"

extern GLApplication *theApp;

static const char*specMode[] = { "MC hits (inc.)","MC hits (abs.)","Flux inc. (ph/sec/.1%BW)","Flux abs. (ph/sec/.1%BW)","Power inc. (W/.1%BW)","Power abs. (W/.1%BW)" };

SpectrumPlotter::SpectrumPlotter():GLWindow() {

	int wD = 750;
	int hD = 350;

	SetTitle("Spectrum plotter");
	SetIconfiable(true);
	nbView = 0;
	worker = NULL;
	lastUpdate = 0.0f;

	nbColors = 8;
	colors[0] = new GLColor(); colors[0]->r = 255; colors[0]->g = 000; colors[0]->b = 055; //red
	colors[1] = new GLColor(); colors[1]->r = 000; colors[1]->g = 000; colors[1]->b = 255; //blue
	colors[2] = new GLColor(); colors[2]->r = 000; colors[2]->g = 204; colors[2]->b = 051; //green
	colors[3] = new GLColor(); colors[3]->r = 000; colors[3]->g = 000; colors[3]->b = 000; //black
	colors[4] = new GLColor(); colors[4]->r = 255; colors[4]->g = 153; colors[4]->b = 051; //orange
	colors[5] = new GLColor(); colors[5]->r = 153; colors[5]->g = 204; colors[5]->b = 255; //light blue
	colors[6] = new GLColor(); colors[6]->r = 153; colors[6]->g = 000; colors[6]->b = 102; //violet
	colors[7] = new GLColor(); colors[7]->r = 255; colors[7]->g = 230; colors[7]->b = 005; //yellow

	chart = new GLChart(0);
	chart->SetBorder(BORDER_BEVEL_IN);
	chart->GetY1Axis()->SetGridVisible(true);
	chart->GetXAxis()->SetGridVisible(true);
	chart->GetY1Axis()->SetAutoScale(true);
	chart->GetY2Axis()->SetAutoScale(true);
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
	specCombo->SetEditable(true);
	Add(specCombo);

	logToggle = new GLToggle(0,"Log Y scale");
	logToggle->SetState(true);
	Add(logToggle);

	normToggle = new GLToggle(0,"Normalize");
	normToggle->SetState(true);
	Add(normToggle);

	// Center dialog
	int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;
	SetBounds(xD,yD,wD,hD);
	SetResizable(true);
	SetMinimumSize(wD,220);

	RestoreDeviceObjects();

}

void SpectrumPlotter::SetBounds(int x,int y,int w,int h) {
	chart->SetBounds(7,5,w-15,h-60);
	specCombo->SetBounds(7,h-45,167,19);
	selButton->SetBounds(180,h-45,80,19);
	addButton->SetBounds(265,h-45,80,19);
	removeButton->SetBounds(350,h-45,80,19);
	resetButton->SetBounds(435,h-45,80,19);
	logToggle->SetBounds(537,h-45,50,19);
	normToggle->SetBounds(587,h-45,105,19);
	dismissButton->SetBounds(w-100,h-45,90,19);
	GLWindow::SetBounds(x,y,w,h);
}

void SpectrumPlotter::Refresh() {

	if(!worker) return;

	Geometry *geom = worker->GetGeometry();
	size_t nb = geom->GetNbFacet();
	size_t nbSpec = 0;
	for(size_t i=0;i<nb;i++)
		if(geom->GetFacet(i)->sh.recordSpectrum) nbSpec+=6;
	specCombo->Clear();
	specCombo->SetSelectedIndex(-1);
	if(nbSpec) specCombo->SetSize(nbSpec);
	nbSpec=0;
	for(size_t i=0;i<nb;i++) {
		Facet *f = geom->GetFacet(i);
		if(f->sh.recordSpectrum) {
			char tmp[128];
			for (size_t mode=0;mode<6;mode++) { //Flux, power
				sprintf(tmp,"F#%zd %s",i+1,specMode[mode]);
				specCombo->SetValueAt(nbSpec,tmp,(int)(i*6+mode));
				specCombo->SetSelectedIndex(0);
				nbSpec++;
			}
		}
	}

	refreshViews();

}

void SpectrumPlotter::SetScale() {
	if (worker) {
	if (worker->regions.size()>0) {
		Region_full *traj=&(worker->regions[0]); //scale axis X
		if (traj->isLoaded) {
			chart->GetXAxis()->SetMinimum(traj->params.energy_low_eV);
			chart->GetXAxis()->SetMaximum(traj->params.energy_hi_eV);
			delta=(log10(traj->params.energy_hi_eV)-log10(traj->params.energy_low_eV))/SPECTRUM_SIZE;
		}
	}
	}
}

void SpectrumPlotter::Display(Worker *w) {

	worker = w;
	Refresh();
	SetVisible(true);
	SetScale();
}

void SpectrumPlotter::Update(float appTime,bool force) {

	if(!IsVisible() || IsIconic()) return;  

	if(force) {
		refreshViews();
		lastUpdate = appTime;
		return;
	}

	if( (appTime-lastUpdate>1.0f || force) && nbView ) {
		if(worker->isRunning) refreshViews();
		lastUpdate = appTime;
	}

}

void SpectrumPlotter::refreshViews() {

	// Lock during update
	BYTE *buffer = worker->GetHits();
	int normalize = normToggle->GetState();

	if(!buffer) return;

	Geometry *geom = worker->GetGeometry();

	for(int i=0;i<nbView;i++) {

		GLDataView *v = views[i];
		if( v->userData1>=0 && v->userData1<geom->GetNbFacet()) {
			Facet *f = geom->GetFacet(v->userData1);
			int mode=v->userData2;
			v->Reset();

			size_t profileSize=f->sh.isProfile?PROFILE_SIZE*sizeof(ProfileSlice):0;
			size_t textureSize=f->sh.isTextured?f->sh.texWidth*f->sh.texHeight*sizeof(TextureCell):0;
			size_t directionSize=f->sh.countDirection?f->sh.texWidth*f->sh.texHeight*sizeof(DirectionCell):0;

			ProfileSlice* spectrum = (ProfileSlice*)(buffer + (f->sh.hitOffset + sizeof(FacetHitBuffer) + profileSize + textureSize + directionSize));
			
			double max;			

			switch(normalize) {
			case 0: //no normalization
			{
				// Bandwidth correction
				// Preferred bin width: 0.1% * Bin_lower_limit
				// Recorded bin width: Bin_higher_limit - Bin_lower_limit = Bin_lower_limit*10^delta - Bin_lower_limit = Bin_lower_limit*(10^delta - 1)
				// Preferred / recorded: 0.1% / (10^delta - 1)
				double log10_min = log10(chart->GetXAxis()->GetMinimum());
				double log10_max = log10(chart->GetXAxis()->GetMaximum());
				double bandwidthCorrection = 0.001 / (Pow10(delta) - 1.0);

				if (mode == 0) {
					for (size_t j = 0; j < SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10_min + (0.5 + (double)j)*delta), (double)spectrum[j].count_incident);
				}
				else if (mode == 1) {
					for (size_t j = 0; j < SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10_min + (0.5 + (double)j)*delta), (double)spectrum[j].count_absorbed);
				}
				else if (mode == 2) {
					for (size_t j = 0; j < SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10_min + (0.5 + (double)j)*delta), spectrum[j].flux_incident / worker->no_scans * bandwidthCorrection, false); //0.5: point should be center of bin
				}
				else if (mode == 3) {
					for (size_t j = 0; j < SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10_min + (0.5 + (double)j)*delta), spectrum[j].flux_absorbed / worker->no_scans * bandwidthCorrection, false); //0.5: point should be center of bin
				}
				else if (mode == 4) {
					for (size_t j = 0; j < SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10_min + (0.5 + (double)j)*delta), spectrum[j].power_incident / worker->no_scans * bandwidthCorrection, false); //0.5: point should be center of bin
				}
				else if (mode == 5) {
					for (size_t j = 0; j < SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10_min + (0.5 + (double)j)*delta), spectrum[j].power_absorbed / worker->no_scans * bandwidthCorrection, false); //0.5: point should be center of bin
				}
				break;
			}
			case 1: //normalize max. value to 1
				if (mode == 0) {
					max = 0.0;
					for (size_t j = 0; j < SPECTRUM_SIZE; j++) {
						if ((double)spectrum[j].count_incident > max) max = (double)spectrum[j].count_incident;
					}
					for(size_t j=0;j<SPECTRUM_SIZE;j++)
						v->Add(Pow10(log10(chart->GetXAxis()->GetMinimum()) + (0.5 + (double)j)*delta), (double)spectrum[j].count_incident / max, false);
				}
				else if (mode == 1) {
					max = 0.0;
					for (size_t j = 0; j < SPECTRUM_SIZE; j++) {
						if ((double)spectrum[j].count_absorbed > max) max = (double)spectrum[j].count_absorbed;
					}
					for (size_t j = 0; j<SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10(chart->GetXAxis()->GetMinimum()) + (0.5 + (double)j)*delta), (double)spectrum[j].count_absorbed / max, false);
				}
				else if (mode == 2) {
					max = 0.0;
					for (size_t j = 0; j < SPECTRUM_SIZE; j++) {
						if (spectrum[j].flux_incident > max) max = spectrum[j].flux_incident;
					}
					for (size_t j = 0; j<SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10(chart->GetXAxis()->GetMinimum()) + (0.5 + (double)j)*delta), spectrum[j].flux_incident / max, false);
				}
				else if (mode == 3) {
					max = 0.0;
					for (size_t j = 0; j < SPECTRUM_SIZE; j++) {
						if (spectrum[j].flux_absorbed > max) max = spectrum[j].flux_absorbed;
					}
					for (size_t j = 0; j<SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10(chart->GetXAxis()->GetMinimum()) + (0.5 + (double)j)*delta), spectrum[j].flux_absorbed / max, false);
				}
				else if (mode == 4) {
					max = 0.0;
					for (size_t j = 0; j < SPECTRUM_SIZE; j++) {
						if (spectrum[j].power_incident > max) max = spectrum[j].power_incident;
					}
					for (size_t j = 0; j<SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10(chart->GetXAxis()->GetMinimum()) + (0.5 + (double)j)*delta), spectrum[j].power_incident / max, false);
				}
				else if (mode == 5) {
					max = 0.0;
					for (size_t j = 0; j < SPECTRUM_SIZE; j++) {
						if (spectrum[j].power_absorbed > max) max = spectrum[j].power_absorbed;
					}
					for (size_t j = 0; j<SPECTRUM_SIZE; j++)
						v->Add(Pow10(log10(chart->GetXAxis()->GetMinimum()) + (0.5 + (double)j)*delta), spectrum[j].power_absorbed / max, false);
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
	bool found = false;
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

void SpectrumPlotter::remView(int facet,int mode) {

	Geometry *geom = worker->GetGeometry();

	bool found = false;
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
			SetVisible(false);
		} else if(src==selButton) {
			int idx = specCombo->GetSelectedIndex();
			if (idx>=0) {
				geom->UnselectAll();
				int facetId=(int)((double)specCombo->GetUserValueAt(idx)/6.0);
				geom->GetFacet(facetId)->selected = true;
				mApp->UpdateFacetParams(true);
				geom->UpdateSelection();
				mApp->facetList->SetSelectedRow(facetId);
				mApp->facetList->ScrollToVisible(facetId,1,true);
			}
		} else if(src==addButton) {
			int idx = specCombo->GetSelectedIndex();
			if (idx >= 0) {
				int facetId = (int)((double)specCombo->GetUserValueAt(idx) / 6.0);
				addView(facetId, specCombo->GetUserValueAt(idx) % 6);
				refreshViews();
			}
		} else if(src==removeButton) {
			int idx = specCombo->GetSelectedIndex();
			if (idx >= 0) {
				int facetId = (int)((double)specCombo->GetUserValueAt(idx) / 6.0);
				remView(facetId, specCombo->GetUserValueAt(idx) % 6);
				refreshViews();
			}
		} else if(src==resetButton) {
			Reset();
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

