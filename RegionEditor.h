/*
File:        RegionEditor.h
Description: Region Editor
Program:     SynRad

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*/

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLGradient.h"
#include "GLApp\GLFileBox.h"

#ifndef _REGIONEDITORH_
#define _REGIONEDITORH_

class RegionEditor : public GLWindow {

public:

	// Construction
	RegionEditor();

	// Implementation
	void ProcessMessage(GLComponent *src,int message);
	void Display(Worker *w,int regionId);
	

private:

	Worker      *worker;
	Region_full      *cr; //current region
	int regionId;
	
	GLTitledPanel *trajectoryPanel,*beamPanel,*generationPanel,*magneticPanel;

	//trajectory properties
	GLTextField *startPointXtext,*startPointYtext,*startPointZtext;
	GLCombo    *startDirDefinitionCombo;
	GLTextField *theta0text,*alpha0text,*startDirXtext,*startDirYtext,*startDirZtext;
	GLTextField *dLtext,*limitsXtext,*limitsYtext,*limitsZtext;
	GLButton *startDirInfoButton,*dlInfoButton,*limitsInfoButton,*beamsizeInfoButton;
	GLButton *bxyBrowseButton,*magxBrowseButton,*magyBrowseButton,*magzBrowseButton;
	GLButton *bxyEditButton,*magxEditButton,*magyEditButton,*magzEditButton;

	//beam properties
	GLTextField *particleMassText,*beamEnergyText,*beamCurrentText;
	GLCombo *particleChargeCombo;
	GLButton *setParticleElectronButton,*setParticlePositronButton,*setParticleProtonButton;

	GLToggle *idealBeamToggle;
	GLCombo *betaOrBXYCombo;
	GLTextField *betaXtext,*betaYtext,*etaText,*etaPrimeText,*emittanceText,*couplingText,*energySpreadText;
	//GLButton *browseBXYbutton;
	GLTextField *BXYfileNameText;

	//generation properties
	GLTextField *EminText,*EmaxText,*psiMaxXtext,*psiMaxYtext;
	GLToggle *limitAngleToggle;
	GLToggle *enableParPolarizationToggle,*enableOrtPolarizationToggle;

	//magnetic fields
	GLCombo *BxtypeCombo,*BytypeCombo,*BztypeCombo;
	GLTextField *constBXtext,*constBYtext,*constBZtext;
	GLTextField *MAGfileXtext,*MAGfileYtext,*MAGfileZtext;
	//GLButton *browseMAGXbutton,*browseMAGYbutton,*browseMAGZbutton;

	//common stuff
	GLButton *applyButton,*dismissButton;

	void EnableDisableControls();
	void FillValues();
	void ApplyChanges();
};

#endif /* _REGIONEDITORH_ */
