/*
File:        RegionInfo.h
Description: Region Info
*/

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "RegionEditor.h"

#include "Geometry.h"
#include "Worker.h"

#ifndef _TRAJINFOH_
#define _TRAJINFOH_

class RegionInfo : public GLWindow {

public:

	// Construction
	RegionInfo(Worker *work);
	// Implementation
	void ProcessMessage(GLComponent *src,int message);
	void Update();
	void ExportPoints(int regionId,int exportFrequency=1,BOOL doFullScan=FALSE);

private:

	RegionEditor     *regionEditor;
	Region   *selectedRegion;
	GLButton    *notepadButton,*cancelButton,*reloadButton,*exportButton,*editButton,*saveAsButton;
	GLLabel     *pathLabel,*t1,*t2,*t3,*t4,*t5,*t6,*freqLabel;
	GLCombo     *regionSelector;
	GLToggle    *integrateToggle;
	GLTextField *freqField;
	GLTitledPanel *exportPanel;
	Worker *work;
	int exportFrequency;
	int calcIntegrals;
};

#endif /* _TRAJINFOH_ */
