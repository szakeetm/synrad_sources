/*
File:        TrajectoryDetails.cpp
Description: Trajectory details window
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

#include "TrajectoryDetails.h"
#include "GeneratePhoton.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Utils.h"
#include "SynRad.h"
#include "Random.h"

extern SynRad *mApp;

typedef struct {

	char *name;
	int   width;
	int   align;

} COLUMN;

static COLUMN allColumn[] = {
	{ "#", 40, ALIGN_CENTER },
	{ "L (cm)", 100, ALIGN_CENTER },
	{ "Orbit_pos_X (cm)", 100, ALIGN_CENTER },
	{ "Orbit_pos_Y (cm)", 100, ALIGN_CENTER },
	{ "Orbit_pos_Z (cm)", 100, ALIGN_CENTER },
	{ "Orbit_Alpha (rad)", 100, ALIGN_CENTER },
	{ "Orbit_Theta (rad)", 100, ALIGN_CENTER },
	{ "Orbit_B (T)", 100, ALIGN_CENTER },
	{ "Orbit_Bx (T)", 100, ALIGN_CENTER },
	{ "Orbit_By (T)", 100, ALIGN_CENTER },
	{ "Orbit_Bz (T)", 100, ALIGN_CENTER },
	{ "Orbit_Radius (cm)", 100, ALIGN_CENTER },
	{ "Orbit_Critical_E (eV)", 110, ALIGN_CENTER },
	{ "Orbit_Emittance_X (cm)", 110, ALIGN_CENTER },
	{ "Orbit_Emittance_Y (cm)", 110, ALIGN_CENTER },
	{ "Orbit_Beta_X (cm)", 100, ALIGN_CENTER },
	{ "Orbit_Beta_Y (cm)", 100, ALIGN_CENTER },
	{ "Orbit_Eta", 100, ALIGN_CENTER },
	{ "Orbit_Eta_Prime", 100, ALIGN_CENTER },
	{ "Orbit_Energy_Spread", 110, ALIGN_CENTER },
	{ "Orbit_Sigma_X (cm)", 100, ALIGN_CENTER },
	{ "Orbit_Sigma_Y (cm)", 100, ALIGN_CENTER },
	{ "Orbit_Sigma_X_Prime (rad)", 110, ALIGN_CENTER },
	{ "Orbit_Sigma_Y_Prime (rad)", 110, ALIGN_CENTER },

	{ "Photon_Natural_div_X (rad)", 110, ALIGN_CENTER },
	{ "Photon_Natural_div_Y (rad)", 110, ALIGN_CENTER },
	{ "Photon_Offset_X (cm)", 110, ALIGN_CENTER },
	{ "Photon_Offset_Y (cm)", 110, ALIGN_CENTER },
	{ "Photon_div_X (rad)", 100, ALIGN_CENTER },
	{ "Photon_div_Y (rad)", 100, ALIGN_CENTER },
	{ "Photon_B (T)", 100, ALIGN_CENTER },
	{ "Photon_Bx (T)", 100, ALIGN_CENTER },
	{ "Photon_By (T)", 100, ALIGN_CENTER },
	{ "Photon_Bz (T)", 100, ALIGN_CENTER },
	{ "Photon_B_ortho (T)", 100, ALIGN_CENTER },
	{ "Photon_B_par (T)", 100, ALIGN_CENTER },
	{ "Photon_Radius (cm)", 100, ALIGN_CENTER },
	{ "Photon_Critical_Energy (eV)", 110, ALIGN_CENTER },
	{ "Photon_Energy (eV)", 100, ALIGN_CENTER },
	{ "Photon_E/E_crit", 100, ALIGN_CENTER },
	{ "Photon_G1H2", 100, ALIGN_CENTER },
	{ "Photon_B_factor", 100, ALIGN_CENTER },
	{ "Photon_B_factor_P", 100, ALIGN_CENTER },
	{ "Photon_Flux (1/s)", 100, ALIGN_CENTER },
	{ "Photon_Power (W)", 100, ALIGN_CENTER },

	/*{ "Integral_Flux (1/s)", 100, ALIGN_CENTER },
	{ "Integral_Power (W)", 100, ALIGN_CENTER },
	{ "Integral_Bx (T)", 100, ALIGN_CENTER },
	{ "Integral_By (T)", 100, ALIGN_CENTER },
	{ "Integral_Bz (T)", 100, ALIGN_CENTER },*/

};

// -----------------------------------------------------------------

TrajectoryDetails::TrajectoryDetails() :GLWindow() {

	int wD = 1000;
	int hD = 700;
	worker = NULL;

	SetTitle("Trajectory Points");
	SetIconfiable(TRUE);
	SetResizable(TRUE);
	SetMinimumSize(600, 300);



	regionCombo = new GLCombo(0);
	Add(regionCombo);

	updateButton = new GLButton(0, "Update table!");
	Add(updateButton);
	dismissButton = new GLButton(0, "Dismiss");
	Add(dismissButton);

	pointList = new GLList(0);
	pointList->SetColumnLabelVisible(TRUE);
	pointList->SetGrid(TRUE);
	Add(pointList);

	everyLabel = new GLLabel("Display every");
	Add(everyLabel);

	freqText = new GLTextField(0, "");
	Add(freqText);

	nbPointLabel = new GLLabel("");
	Add(nbPointLabel);

	sPanel = new GLTitledPanel("Show column");
	sPanel->SetClosable(TRUE);
	sPanel->Close();
	Add(sPanel);

	for (int i = 1; i < NB_TPCOLUMN; i++) {
		show[i - 1] = new GLToggle(0, allColumn[i].name);
		show[i - 1]->SetState(TRUE);
		sPanel->Add(show[i - 1]);
	}

	checkAllButton = new GLButton(0, "Check All");
	sPanel->Add(checkAllButton);
	uncheckAllButton = new GLButton(0, "Uncheck All");
	sPanel->Add(uncheckAllButton);

	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	RestoreDeviceObjects();

}

// -----------------------------------------------------------------

void TrajectoryDetails::PlaceComponents() {

	// Show toggle panel
	int nbW = (width - 20) / 120;
	int nbL = (NB_TPCOLUMN - 1) / nbW + (((NB_TPCOLUMN - 1) % nbW) ? 1 : 0);
	int hp;
	if (!sPanel->IsClosed())
		hp = 20 * (nbL + 1) + 25;
	else
		hp = 20;
	sPanel->SetBounds(5, height - (hp + 52), width - 10, hp);
	for (int i = 0; i < (NB_TPCOLUMN - 1); i++)
		sPanel->SetCompBounds(show[i], 5 + 120 * ((i) % nbW), 18 + 20 * ((i) / nbW), 125, 19);

	pointList->SetBounds(5, 5, width - 10, height - (62 + hp));

	sPanel->SetCompBounds(uncheckAllButton, width / 2 - 90, hp - 25, 90, 19);
	sPanel->SetCompBounds(checkAllButton, width / 2 + 5, hp - 25, 90, 19);

	regionCombo->SetBounds(10, height - 45, 90, 19);
	everyLabel->SetBounds(150, height - 45, 70, 19);
	freqText->SetBounds(225, height - 45, 40, 19);
	nbPointLabel->SetBounds(268, height - 45, 80, 19);
	updateButton->SetBounds(width - 260, height - 45, 155, 19);
	dismissButton->SetBounds(width - 100, height - 45, 90, 19);

}

// -----------------------------------------------------------------

void TrajectoryDetails::SetBounds(int x, int y, int w, int h) {

	GLWindow::SetBounds(x, y, w, h);
	PlaceComponents();

}

char *TrajectoryDetails::FormatCell(int idx, int mode, GenPhoton* photon) {

	Worker *worker = &(mApp->worker);

	static char ret[256];
	strcpy(ret, "");

	Trajectory_Point* p = &worker->regions[displayedRegion].Points[idx];
	Region_full* reg = &worker->regions[displayedRegion];

	switch (mode) {
	case 0: //id
		sprintf(ret, "%d", idx + 1);
		break;
	case 1: //L
		sprintf(ret, "%g", (double)idx*reg->params.dL);
		break;
	case 2: //Orbit_posX
		sprintf(ret, "%g", p->position.x);
		break;
	case 3: //Orbit_posY
		sprintf(ret, "%g", p->position.y);
		break;
	case 4: //Orbit_posZ
		sprintf(ret, "%g", p->position.z);
		break;
	case 5: //Orbit_alpha
		sprintf(ret, "%g", -asin(p->direction.y));
		break;
	case 6:{ //Orbit_theta
		double theta;
		if (p->direction.z == 0) {
			if (p->direction.x >= 0) theta = -PI / 2;
			else theta = PI / 2;
		}
		else theta = -atan(p->direction.x / p->direction.z);
		sprintf(ret, "%g", theta);
		break; }
	case 7: //Orbit_B
		sprintf(ret, "%g", p->B.Norme());
		break;
	case 8: //Orbit_Bx
		sprintf(ret, "%g", p->B.x);
		break;
	case 9: //Orbit_By
		sprintf(ret, "%g", p->B.y);
		break;
	case 10: //Orbit_Bz
		sprintf(ret, "%g", p->B.z);
		break;
	case 11: //Orbit_Radius
		sprintf(ret, "%g", p->rho.Norme());
		break;
	case 12: //Orbit_critical_E
		sprintf(ret, "%g", p->critical_energy);
		break;
	case 13: //Orbit_Emittance_X
		sprintf(ret, "%g", p->emittance_X);
		break;
	case 14: //Orbit_Emittance_Y
		sprintf(ret, "%g", p->emittance_Y);
		break;
	case 15: //Orbit_Beta_X
		sprintf(ret, "%g", p->beta_X);
		break;
	case 16: //Orbit_Beta_Y
		sprintf(ret, "%g", p->beta_Y);
		break;
	case 17: //Orbit_Eta
		sprintf(ret, "%g", p->eta);
		break;
	case 18: //Orbit_Eta_Prime
		sprintf(ret, "%g", p->eta_prime);
		break;
	case 19: //Orbit_Energy_Spread
		sprintf(ret, "%g %%", p->energy_spread);
		break;
	case 20: //Orbit_Sigma_X
		sprintf(ret, "%g", p->sigma_x);
		break;
	case 21: //Orbit_Sigma_Y
		sprintf(ret, "%g", p->sigma_y);
		break;
	case 22: //Orbit_Sigma_X_Prime
		sprintf(ret, "%g", p->sigma_x_prime);
		break;
	case 23: //Orbit_Sigma_Y_Prime
		sprintf(ret, "%g", p->sigma_y_prime);
		break;

	case 24: //Photon_Natural_div_X (rad)
		sprintf(ret, "%g", photon->natural_divx);
		break;
	case 25: //Photon_Natural_div_Y (rad)
		sprintf(ret, "%g", photon->natural_divy);
		break;
	case 26: //Photon_Offset_X (cm)
		sprintf(ret, "%g", photon->offset_x);
		break;
	case 27: //Photon_Offset_Y (cm)
		sprintf(ret, "%g", photon->offset_y);
		break;
	case 28: //Photon_div_X (rad)
		sprintf(ret, "%g", photon->offset_divx);
		break;
	case 29: //Photon_div_Y (rad)
		sprintf(ret, "%g", photon->offset_divy);
		break;
	case 30: //Photon_B (T)
		sprintf(ret, "%g", photon->B.Norme());
		break;
	case 31: //Photon_Bx (T)
		sprintf(ret, "%g", photon->B.x);
		break;
	case 32: //Photon_By (T)
		sprintf(ret, "%g", photon->B.y);
		break;
	case 33: //Photon_Bz (T)
		sprintf(ret, "%g", photon->B.z);
		break;
	case 34: //Photon_B_ortho (T)
		sprintf(ret, "%g", photon->B_ort);
		break;
	case 35: //Photon_B_par (T)
		sprintf(ret, "%g", photon->B_par);
		break;
	case 36: //Photon_Radius (cm)
		sprintf(ret, "%g", photon->radius);
		break;
	case 37: //Photon_Critical_Energy (eV)
		sprintf(ret, "%g", photon->critical_energy);
		break;
	case 38: //Photon_Energy (eV)
		sprintf(ret, "%g", photon->energy);
		break;
	case 39: //Photon_E/E_crit
		sprintf(ret, "%g", photon->energy / photon->critical_energy);
		break;
	case 40: //Photon_G1H2
		sprintf(ret, "%s", "g1h2 obsolete");
		break;
	case 41: //Photon_B_factor
		sprintf(ret, "%g", photon->B_factor);
		break;
	case 42: //Photon_B_factor_P
		sprintf(ret, "%g", photon->B_factor_power);
		break;
	case 43: //Photon_Flux
		sprintf(ret, "%g", photon->SR_flux);
		break;
	case 44: //Photon_Power (W)
		sprintf(ret, "%g", photon->SR_power);
		break;


	default:
		sprintf(ret, "%s", "");
		break;
	}

	return ret;

}

// -----------------------------------------------------------------

void TrajectoryDetails::UpdateTable() {

	GLProgress *updatePrg = new GLProgress("Updating table...", "Please wait");
	updatePrg->SetProgress(0.0);
	updatePrg->SetVisible(TRUE);

	//creating distributions
	Distribution2D K_1_3_distribution = Generate_K_Distribution(1.0 / 3.0);
	Distribution2D K_2_3_distribution = Generate_K_Distribution(2.0 / 3.0);
	Distribution2D integral_N_photons = Generate_Integral(LOWER_LIMIT, UPPER_LIMIT, INTEGRAL_MODE_N_PHOTONS);
	Distribution2D integral_SR_power = Generate_Integral(LOWER_LIMIT, UPPER_LIMIT, INTEGRAL_MODE_SR_POWER);
	//Distribution2D polarization_distribution = Generate_Polarization_Distribution(true, true);
	//Distribution2D g1h2_distribution = Generate_G1_H2_Distribution();

	/*double integ_BX = 0.0;
	double integ_BY = 0.0;
	double integ_BZ = 0.0;
	double integ_flux = 0.0;
	double integ_power = 0.0;*/

	//initialize random number generator
	DWORD seed = (DWORD)((int)(mApp->GetTick()*1000.0)*_getpid());
	rseed(seed);


	static char ret[256];
	strcpy(ret, "");

	char *tmpName[NB_TPCOLUMN];
	int  tmpWidth[NB_TPCOLUMN];
	int  tmpAlign[NB_TPCOLUMN];


	int  nbCol = 0;

	for (int i = 0; i < NB_TPCOLUMN; i++) {
		if (i == 0 || show[i - 1]->GetState()) {
			tmpName[nbCol] = allColumn[i].name;
			tmpWidth[nbCol] = allColumn[i].width;
			tmpAlign[nbCol] = allColumn[i].align;
			shown[nbCol] = i;
			nbCol++;
		}
	}

	int nbPoints = (int)worker->regions[displayedRegion].Points.size();

	pointList->SetSize(nbCol, (int)((double)nbPoints / (double)freq)+1);
	pointList->SetColumnWidths(tmpWidth);
	pointList->SetColumnLabels(tmpName);
	pointList->SetColumnAligns(tmpAlign);

	for (int pointId = 0; pointId < nbPoints; pointId += freq) {
 		GenPhoton photon = GeneratePhoton(pointId, &worker->regions[displayedRegion], worker->generation_mode, worker->psi_distr, worker->chi_distr,pointId == 0);
		updatePrg->SetProgress((double)pointId / (double)nbPoints);
		for (int j = 0; j < nbCol; j++)
			pointList->SetValueAt(j, (int)((double)pointId / (double)freq), FormatCell(pointId, shown[j], &photon));
	}
	updatePrg->SetVisible(FALSE);
	SAFE_FREE(updatePrg);
}

// -----------------------------------------------------------------

void TrajectoryDetails::Update() {

	if (!worker) return;
	if (!IsVisible()) return;

	UpdateTable();

}

// -----------------------------------------------------------------

void TrajectoryDetails::Display(Worker *w, int regionId) {

	worker = w;
	displayedRegion = regionId;
	freq = MAX(1, (int)((double)worker->regions[displayedRegion].Points.size() / 1000.0)); //by default, display around 1000 trajectory points
	char tmp[1024];

	regionCombo->SetSize((int)worker->regions.size());
	for (size_t r = 0; r < worker->regions.size(); r++){

		sprintf(tmp, "%s %d", "Region", r + 1);
		regionCombo->SetValueAt(r, tmp);
	}
	sprintf(tmp, "%d", freq);
	freqText->SetText(tmp);
	sprintf(tmp, "th point (%d in total)", (int)((double)worker->regions[displayedRegion].Points.size() / (double)freq));
	nbPointLabel->SetText(tmp);
	regionCombo->SetSelectedIndex(regionId);
	SetVisible(TRUE);
	Update();

}

// -----------------------------------------------------------------

void TrajectoryDetails::ProcessMessage(GLComponent *src, int message) {

	switch (message) {

	case MSG_BUTTON:
		if (src == dismissButton) {
			SetVisible(FALSE);
		}
		else if (src == checkAllButton) {
			for (int i = 0; i < NB_TPCOLUMN - 1; i++) show[i]->SetState(TRUE);
			//UpdateTable();
		}
		else if (src == uncheckAllButton) {
			for (int i = 0; i < NB_TPCOLUMN - 1; i++) show[i]->SetState(FALSE);
			//UpdateTable();
		}
		else if (src == updateButton) {
			int freq_temp;
			if (!freqText->GetNumberInt(&freq_temp) || !(freq_temp >= 1)) {
				GLMessageBox::Display("Invalid point frequency", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			displayedRegion = regionCombo->GetSelectedIndex();
			freq = freq_temp;
			UpdateTable();
		}
		break;

	case MSG_TOGGLE:
		//UpdateTable();
		break;

		/*case MSG_LIST_COL:
			if( src==facetListD ) {
			// Save column width
			int c = facetListD->GetDraggedCol();
			allColumn[shown[c]].width = facetListD->GetColWidth(c);
			}
			break;*/

	case MSG_PANELR:
		PlaceComponents();
		break;

	case MSG_LIST:
		worker->regions[displayedRegion].selectedPoint = pointList->GetSelectedRow(TRUE);
		break;
	case MSG_TEXT_UPD:{
		int freq_temp;
		if (freqText->GetNumberInt(&freq_temp)) {
			char tmp[128];
			sprintf(tmp, "th point (%d in total)", (int)((double)worker->regions[displayedRegion].Points.size() / (double)freq_temp));
			nbPointLabel->SetText(tmp);
		}

		break;
	}
	}

	GLWindow::ProcessMessage(src, message);
}

int TrajectoryDetails::GetRegionId() {
	return displayedRegion;
}

void TrajectoryDetails::SelectPoint(int idx) {
	pointList->SetSelectedRow((int)((double)idx / (double)freq));
	pointList->ScrollToVisible((int)((double)idx / (double)freq), 0, FALSE);
}