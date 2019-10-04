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
#ifndef _REGIONEDITORH_
#define _REGIONEDITORH_

#include <GLApp/GLCombo.h>
#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLGradient.h"

class Worker;
class Region_full;

class RegionEditor : public GLWindow {

public:

	// Construction
	RegionEditor();

	// Implementation
	void ProcessMessage(GLComponent *src,int message);
	void Display(Worker *w,int regionId);
	int GetRegionId();
	

private:

	Worker      *worker;
	Region_full      *cr; //current region
	int regionId;
	
	GLTitledPanel	*beamPanel;
	GLLabel	*label6;
	GLLabel	*label4;
	GLTextField	*startPointYtext;
	GLTextField	*startPointZtext;
	GLLabel	*label7;
	GLLabel	*label5;
	GLLabel	*label3;
	GLTextField	*startPointXtext;
	GLLabel	*label2;
	GLLabel	*label1;
	GLCombo	*startDirDefinitionCombo;
	GLLabel	*label8;
	GLLabel	*label9;
	GLLabel	*label10;
	GLTextField	*theta0text;
	GLTextField	*alpha0text;
	GLLabel	*label11;
	GLLabel	*label12;
	GLButton	*startDirInfoButton;
	GLTextField	*startDirZtext;
	GLLabel	*label15;
	GLTextField	*startDirYtext;
	GLLabel	*label14;
	GLTextField	*startDirXtext;
	GLLabel	*label13;
	GLLabel	*label19;
	GLLabel	*label20;
	GLTextField	*limitsYtext;
	GLTextField	*limitsZtext;
	GLLabel	*label21;
	GLLabel	*label22;
	GLLabel	*label23;
	GLTextField	*limitsXtext;
	GLLabel	*label24;
	GLLabel	*label25;
	GLLabel	*label16;
	GLTextField	*dLtext;
	GLLabel	*label17;
	GLLabel	*label18;
	GLButton	*limitsInfoButton;
	GLButton	*dlInfoButton;
	GLTitledPanel	*particlePanel;
	GLLabel	*label27;
	GLLabel	*label26;
	GLTextField	*particleMassText;
	GLLabel	*label28;
	GLCombo	*particleChargeCombo;
	GLLabel	*label29;
	GLButton	*setParticleProtonButton;
	GLButton	*setParticlePositronButton;
	GLButton	*setParticleElectronButton;
	GLLabel	*label30;
	GLTitledPanel	*beamParamPanel;
	GLLabel	*label33;
	GLTextField	*beamCurrentText;
	GLLabel	*label34;
	GLLabel	*label31;
	GLTextField	*beamEnergyText;
	GLLabel	*label32;
	GLToggle	*idealBeamToggle;
	GLToggle	*emittanceCouplingToggle;
	GLTextField	*emittanceYtext;
	GLToggle	*emittanceXYtoggle;
	GLTextField	*couplingText;
	GLTextField	*emittanceXtext;
	GLTextField	*emittanceText;
	GLLabel	*label36;
	GLLabel	*label35;
	GLButton	*beamsizeInfoButton;
	GLLabel	*label40;
	GLLabel	*label39;
	GLLabel	*label38;
	GLLabel	*label37;
	GLLabel	*label49;
	GLLabel	*label47;
	GLLabel	*label48;
	GLLabel	*label50;
	GLTextField	*etaPrimeText;
	GLTextField	*energySpreadText;
	GLLabel	*label45;
	GLLabel	*label46;
	GLTextField	*etaText;
	GLLabel	*label43;
	GLLabel	*label44;
	GLTextField	*betaYtext;
	GLLabel	*label41;
	GLLabel	*label42;
	GLTextField	*betaXtext;
	GLButton	*gammaInfoButton;
	GLTextField	*BXYfileNameText;
	GLLabel	*label51;
	GLToggle	*constantBXYtoggle;
	GLButton	*bxyEditButton;
	GLButton	*bxyBrowseButton;
	GLTitledPanel	*photonGenPanel;
	GLToggle	*enableParPolarizationToggle;
	GLLabel	*label59;
	GLLabel	*label57;
	GLTextField	*EmaxText;
	GLLabel	*label58;
	GLLabel	*label54;
	GLLabel	*label56;
	GLTextField	*EminText;
	GLLabel	*label55;
	GLToggle	*enableOrtPolarizationToggle;
	GLLabel	*label61;
	GLLabel	*label60;
	GLTextField	*psiMaxXtext;
	GLTextField	*psiMaxYtext;
	GLToggle	*limitAngleToggle;
	GLLabel	*label63;
	GLLabel	*label62;
	GLTitledPanel	*magPanel;
	GLCombo	*BxtypeCombo;
	GLLabel	*label66;
	GLLabel	*label65;
	GLButton	*magxEditButton;
	GLTextField	*constBXtext;
	GLButton	*magxBrowseButton;
	GLTextField	*MAGfileXtext;
	GLLabel	*label64;
	GLCombo	*BztypeCombo;
	GLLabel	*label70;
	GLLabel	*label71;
	GLButton	*magzEditButton;
	GLTextField	*constBZtext;
	GLButton	*magzBrowseButton;
	GLTextField	*MAGfileZtext;
	GLLabel	*label72;
	GLCombo	*BytypeCombo;
	GLLabel	*label67;
	GLLabel	*label68;
	GLButton	*magyEditButton;
	GLTextField	*constBYtext;
	GLButton	*magyBrowseButton;
	GLTextField	*MAGfileYtext;
	GLLabel	*label69;
	GLButton	*applyButton;
	GLTextField	*copyFromRegionTextfield;
	GLButton	*copyFromRegionButton;
	GLLabel	*label52;
	GLTextField	*startingStructureTextfield;

	void EnableDisableControls(GLComponent* src);
	void FillValues();
	void ApplyChanges();
	void CopyFromRegion(const size_t& sourceRegionId);
};

#endif /* _REGIONEDITORH_ */
