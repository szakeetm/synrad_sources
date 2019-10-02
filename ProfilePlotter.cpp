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
#include "ProfilePlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Facet_shared.h"
#include "GLApp/MathTools.h"

#include <math.h>
#include "Synrad.h"

extern SynRad *mApp;

static const char*profType[] = { "None","\201","\202","Angle" };
static const char*profMode[] = { "MC Hits (incident)","MC Hits (absorbed)","SR Flux/cm\262 (incident)","SR Flux/cm\262 (absorbed)","SR Power/cm\262 (incident)","SR Power/cm\262 (absorbed)" };

ProfilePlotter::ProfilePlotter() :GLWindow() {

	int wD = 750;
	int hD = 400;

	SetTitle("Profile plotter");
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
	Add(chart);

	dismissButton = new GLButton(0, "Dismiss");
	Add(dismissButton);

	selButton = new GLButton(0, "Show Facet");
	Add(selButton);

	addButton = new GLButton(0, "Add curve");
	Add(addButton);

	removeButton = new GLButton(0, "Remove curve");
	Add(removeButton);

	resetButton = new GLButton(0, "Remove all");
	Add(resetButton);

	profCombo = new GLCombo(0);
	profCombo->SetEditable(true);
	Add(profCombo);

	logToggle = new GLToggle(0, "Log Y scale");
	logToggle->SetState(false);
	Add(logToggle);

	normToggle = new GLToggle(0, "Normalize");
	normToggle->SetState(false);
	Add(normToggle);

	formulaText = new GLTextField(0, "");
	formulaText->SetEditable(true);
	Add(formulaText);

	formulaBtn = new GLButton(0, "-> Plot");
	Add(formulaBtn);

	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);
	SetResizable(true);
	SetMinimumSize(wD, 220);

	RestoreDeviceObjects();

}

void ProfilePlotter::SetBounds(int x, int y, int w, int h) {

	chart->SetBounds(7, 5, w - 15, h - 85);
	profCombo->SetBounds(7, h - 70, 177, 19);
	selButton->SetBounds(190, h - 70, 80, 19);
	addButton->SetBounds(275, h - 70, 80, 19);
	removeButton->SetBounds(360, h - 70, 80, 19);
	resetButton->SetBounds(445, h - 70, 80, 19);
	logToggle->SetBounds(547, h - 70, 50, 19);
	normToggle->SetBounds(597, h - 70, 105, 19);
	formulaText->SetBounds(7, h - 45, 525, 19);
	formulaBtn->SetBounds(537, h - 45, 80, 19);
	dismissButton->SetBounds(w - 100, h - 45, 90, 19);

	GLWindow::SetBounds(x, y, w, h);

}

void ProfilePlotter::Refresh() {

	if (!worker) return;

	//Rebuild selection combo box
	Geometry *geom = worker->GetGeometry();
	size_t nb = geom->GetNbFacet();
	size_t nbProf = 0;
	size_t nbProfileModes = 6;
	for (size_t i = 0; i < nb; i++)
		if (geom->GetFacet(i)->sh.isProfile) nbProf += nbProfileModes;
	profCombo->Clear(); profCombo->SetSelectedIndex(-1);
	if (nbProf) profCombo->SetSize(nbProf);
	nbProf = 0;
	for (size_t i = 0; i < nb; i++) {
		Facet *f = geom->GetFacet(i);
		if (f->sh.isProfile) {
			char tmp[128];
			for (size_t mode = 0; mode < nbProfileModes; mode++) { //MC hits, flux, power
				sprintf(tmp, "F#%zd %s %s", i + 1, profType[f->sh.profileType], profMode[mode]);
				profCombo->SetValueAt(nbProf, tmp, (int)(i * nbProfileModes + mode));
				profCombo->SetSelectedIndex(0);
				nbProf++;
			}
		}
	}
	//Remove profiles that aren't present anymore
	for (size_t v = 0; v < nbView; v++) {
		if (views[v]->userData1 >= geom->GetNbFacet() || !geom->GetFacet(views[v]->userData1)->sh.isProfile) {
			chart->GetY1Axis()->RemoveDataView(views[v]);
			SAFE_DELETE(views[v]);
			for (size_t j = v; j < (int)nbView - 1; j++) views[j] = views[j + 1];
			nbView--;
		}
	}

	//Update values
	refreshViews();

}

void ProfilePlotter::Display(Worker *w) {

	worker = w;
	Refresh();
	SetVisible(true);

}

void ProfilePlotter::Update(float appTime, bool force) {

	if (!IsVisible() || IsIconic()) return;

	if (force) {
		refreshViews();
		lastUpdate = appTime;
		return;
	}

	if ((appTime - lastUpdate > 1.0f || force) && nbView) {
		if (worker->isRunning) refreshViews();
		lastUpdate = appTime;
	}

}

void ProfilePlotter::plot() {

	GLParser *parser = new GLParser();
	parser->SetExpression(formulaText->GetText().c_str());
	if (!parser->Parse()) {
		GLMessageBox::Display(parser->GetErrorMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	int nbVar = parser->GetNbVariable();
	if (nbVar == 0) {
		GLMessageBox::Display("Variable 'x' not found", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	if (nbVar > 1) {
		GLMessageBox::Display("Too much variables or unknown constant", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	VLIST *var = parser->GetVariableAt(0);
	if (_stricmp(var->name, "x") != 0) {
		GLMessageBox::Display("Variable 'x' not found", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	Geometry *geom = worker->GetGeometry();
	GLDataView *v;

	// Check that view is not already added
	bool found = false;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData1 == -1);
		if (!found) i++;
	}

	if (found) {
		v = views[i];
		v->SetName(formulaText->GetText().c_str());
		v->Reset();
	}
	else {
		if (nbView < 50) {
			v = new GLDataView();
			v->SetName(formulaText->GetText().c_str());
			v->userData1 = -1;
			chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}
	}

	// Plot
	for (int i = 0; i < 1000; i++) {
		double x = (double)i;
		double y;
		var->value = x;
		parser->Evaluate(&y);
		v->Add(x, y, false);
	}
	v->CommitChange();

	delete parser;

}

void ProfilePlotter::refreshViews() {

	// Lock during update
	BYTE *buffer = worker->GetHits();
	int normalize = normToggle->GetState();

	if (!buffer) return;

	Geometry *geom = worker->GetGeometry();
	GlobalHitBuffer *gHits = (GlobalHitBuffer *)buffer;
	//double nbAbs = gHits->globalHits.hit.nbAbsEquiv;
	double nbDes = (double)gHits->globalHits.hit.nbDesorbed;
	//double nbHit = (double)gHits->globalHits.hit.nbMCHit;

	for (int i = 0; i < nbView; i++) {

		GLDataView *v = views[i];
		if (v->userData1 >= 0 && v->userData1 < geom->GetNbFacet()) {
			Facet *f = geom->GetFacet(v->userData1);
			int mode = v->userData2;
			v->Reset();
			ProfileSlice   *profile = (ProfileSlice *)(buffer + f->sh.hitOffset + sizeof(FacetHitBuffer));
			
			size_t max_MC;
			double max_flux, max_power;
			
			double elemArea = f->sh.area / PROFILE_SIZE;

			switch (normalize) {
			case 0: //linear no normalization
				 //mc hits, put raw value
				if (mode == 0) {
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, (double)profile[j].count_incident, false);
				}
				else if (mode == 1) { //mc hits, put raw value
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, (double)profile[j].count_absorbed, false);
				}
				//flux or power, normalize by no_scans and elemarea
				else if (mode == 2) {
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].flux_incident / worker->no_scans / elemArea, false);
				}
				else if (mode == 3) {
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].flux_absorbed / worker->no_scans / elemArea, false);
				}
				else if (mode == 4) {
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].power_incident / worker->no_scans / elemArea, false);
				}
				else if (mode == 5) {
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].power_absorbed / worker->no_scans / elemArea, false);
				}
				break;

			case 1: //normalize max. value to 1
				if (mode == 0) { //mc hits
					max_MC = 0;
					for (int j = 0; j < PROFILE_SIZE; j++)
						if (profile[j].count_incident > max_MC) max_MC = profile[j].count_incident;
					if (max_MC > 0)
						for (int j = 0; j < PROFILE_SIZE; j++)
							v->Add((double)j, ((double)profile[j].count_incident / (double)max_MC), false);
				}
				else if (mode == 1) { //mc hits
					max_MC = 0;
					for (int j = 0; j < PROFILE_SIZE; j++)
						if (profile[j].count_absorbed > max_MC) max_MC = profile[j].count_absorbed;
					if (max_MC > 0)
						for (int j = 0; j < PROFILE_SIZE; j++)
							v->Add((double)j, ((double)profile[j].count_absorbed / (double)max_MC), false);
				}
				else if (mode == 2) { //flux
					max_flux = 0.0;
					for (int j = 0; j < PROFILE_SIZE; j++)
						if (profile[j].flux_incident > max_flux) max_flux = profile[j].flux_incident;
					if (max_flux > 0.0)
						for (int j = 0; j < PROFILE_SIZE; j++)
							v->Add((double)j, profile[j].flux_incident / max_flux, false);
				}
				else if (mode == 3) { //flux
					max_flux = 0.0;
					for (int j = 0; j < PROFILE_SIZE; j++)
						if (profile[j].flux_absorbed > max_flux) max_flux = profile[j].flux_absorbed;
					if (max_flux > 0.0)
						for (int j = 0; j < PROFILE_SIZE; j++)
							v->Add((double)j, profile[j].flux_absorbed / max_flux, false);
				}
				else if (mode == 4) { //power
					max_power = 0.0;
					for (int j = 0; j < PROFILE_SIZE; j++)
						if (profile[j].power_incident > max_power) max_power = profile[j].power_incident;
					if (max_power > 0.0)
						for (int j = 0; j < PROFILE_SIZE; j++)
							v->Add((double)j, profile[j].power_incident / max_power, false);
				}
				else if (mode == 5) { //power
					max_power = 0.0;
					for (int j = 0; j < PROFILE_SIZE; j++)
						if (profile[j].power_absorbed > max_power) max_power = profile[j].power_absorbed;
					if (max_power > 0.0)
						for (int j = 0; j < PROFILE_SIZE; j++)
							v->Add((double)j, profile[j].power_absorbed / max_power, false);
				}
				break;
			}
			v->CommitChange();
		}
		else {
			if (v->userData1 == -2 && nbDes != 0.0) {

				// Volatile profile
				v->Reset();
				size_t nb = geom->GetNbFacet();
				for (size_t j = 0; j < nb; j++) {
					Facet *f = geom->GetFacet(j);
					if (f->sh.isVolatile) {
						FacetHitBuffer *fCount = (FacetHitBuffer *)(buffer + f->sh.hitOffset);
						double z = geom->GetVertex(f->indices[0])->z;
						v->Add(z, fCount->hit.nbAbsEquiv / nbDes, false);
					}
				}
				// Last
				Facet *f = geom->GetFacet(28);
				FacetHitBuffer *fCount = (FacetHitBuffer *)(buffer + f->sh.hitOffset);
				double fnbAbs = fCount->hit.nbAbsEquiv;
				v->Add(1000.0, fnbAbs / nbDes, false);
				v->CommitChange();
			}
		}
	}

	worker->ReleaseHits();

}

void ProfilePlotter::addView(int facet, int mode) {

	char tmp[128];
	Geometry *geom = worker->GetGeometry();

	// Check that view is not already added
	bool found = false;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData1 == facet && views[i]->userData2 == mode);
		if (!found) i++;
	}
	if (found) {
		GLMessageBox::Display("Profile already plotted", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (nbView < 50) {
		Facet *f = geom->GetFacet(facet);
		GLDataView *v = new GLDataView();
		sprintf(tmp, "F#%d %s %s", facet + 1, profType[f->sh.profileType], profMode[mode]);
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

void ProfilePlotter::remView(int facet, int mode) {

	Geometry *geom = worker->GetGeometry();

	bool found = false;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData1 == facet && views[i]->userData2 == mode);
		if (!found) i++;
	}
	if (!found) {
		GLMessageBox::Display("Profile not plotted", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	chart->GetY1Axis()->RemoveDataView(views[i]);
	SAFE_DELETE(views[i]);
	for (int j = i; j < nbView - 1; j++) views[j] = views[j + 1];
	nbView--;

}

void ProfilePlotter::Reset() {

	chart->GetY1Axis()->ClearDataView();
	for (int i = 0; i < nbView; i++) SAFE_DELETE(views[i]);
	nbView = 0;

}

void ProfilePlotter::ProcessMessage(GLComponent *src, int message) {
	Geometry *geom = worker->GetGeometry();

	switch (message) {
	case MSG_BUTTON:
		if (src == dismissButton) {
			SetVisible(false);
		}
		else if (src == selButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) {
				geom->UnselectAll();
				int facetId = (int)((double)profCombo->GetUserValueAt(idx) / 6.0);
				geom->GetFacet(facetId)->selected = true;
				mApp->UpdateFacetParams(true);
				geom->UpdateSelection();
				mApp->facetList->SetSelectedRow(facetId);
				mApp->facetList->ScrollToVisible(facetId, 1, true);
			}
		}
		else if (src == addButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) {
				int facetId = (int)((double)profCombo->GetUserValueAt(idx) / 6.0);
				addView(facetId, profCombo->GetUserValueAt(idx) % 6);
				refreshViews();
			}
			
		}
		else if (src == removeButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) {
				int facetId = (int)((double)profCombo->GetUserValueAt(idx) / 6.0);
				remView(facetId, profCombo->GetUserValueAt(idx) % 6);
				refreshViews();
			}
		}
		else if (src == resetButton) {
			Reset();
		}
		else if (src == formulaBtn) {
			plot();
		}
	case MSG_TOGGLE:
		if (src == normToggle) {
			refreshViews();
		}
		else if (src == logToggle) {
			chart->GetY1Axis()->SetScale(logToggle->GetState() ? LOG_SCALE : LINEAR_SCALE);
		}
	}

	GLWindow::ProcessMessage(src, message);

}

void ProfilePlotter::SetViews(std::vector<int> views, int mode) {
	Reset();
	for (int view : views)
		addView(view, mode);
	Refresh();
}

std::vector<int> ProfilePlotter::GetViews() {
	std::vector<int>v;
	v.reserve(nbView);
	for (int i = 0; i < nbView; i++)
		v.push_back(views[i]->userData1);
	return v;
}

bool ProfilePlotter::IsLogScaled() {
	return chart->GetY1Axis()->GetScale();
}
void ProfilePlotter::SetLogScaled(bool logScale) {
	chart->GetY1Axis()->SetScale(logScale);
	logToggle->SetState(logScale);
}