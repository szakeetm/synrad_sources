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

#include "Geometry.h"
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
