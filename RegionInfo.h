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

/*
File:        RegionInfo.h
Description: Region Info
*/

/*

#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "RegionEditor.h"

#include "Geometry_shared.h"
#include "Worker.h"
*/

#ifndef _TRAJINFOH_
#define _TRAJINFOH_

#include "GLApp/GLWindow.h"
class Worker;
class Region_full;
class GLButton;
class GLLabel;
class GLCombo;
class GLTextField;
class GLToggle;

class RegionInfo : public GLWindow {

public:

	// Construction
	RegionInfo(Worker *work);
	// Implementation
	void ProcessMessage(GLComponent *src,int message);
	void Update();

private:

	Region_full   *selectedRegion;
	GLButton    *notepadButton,*cancelButton,*reloadButton,*editButton,*saveAsButton,*viewPointsButton;
	GLLabel     *t1,*t2,*t3,*t4,*t5,*t6/*,*freqLabel*/;
	GLCombo     *regionSelector;
	//GLToggle    *integrateToggle;
	GLTextField *pathLabel;/* , *freqField;*/
	//GLTitledPanel *exportPanel;
	Worker *work;
	//int exportFrequency;
	//int calcIntegrals;
};

#endif /* _TRAJINFOH_ */
