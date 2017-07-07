/*
  File:        GlobalSettings.h
  Description: Global settings dialog (units, antialiasing)
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
#ifndef _GLOBALSETTINGSH_
#define _GLOBALSETTINGSH_

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLGradient.h"
#include "Shared.h" //MAX_PROCESS macro

class Worker;
class GLList;

class GlobalSettings : public GLWindow {

public:

  // Construction
  GlobalSettings();
  void SMPUpdate();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void Display(Worker *w);

private:

	void RestartProc();
	  Worker      *worker;
  GLList      *processList;
  GLButton    *restartButton;
  GLButton    *maxButton;
  GLTextField *nbProcText;
  GLTextField *autoSaveText;
  GLTextField *cutoffText;
 
  float lastUpdate;
  float lastCPUTime[MAX_PROCESS];
  float lastCPULoad[MAX_PROCESS];

  GLToggle      *chkAntiAliasing;
  GLToggle      *chkWhiteBg;
  //GLToggle      *chkNonIsothermal;
  GLToggle      *chkSimuOnly;
  GLToggle      *chkCheckForUpdates;
  GLToggle      *chkAutoUpdateFormulas;
  GLToggle      *chkNewReflectionModel;
  GLToggle      *chkCompressSavedFiles;
  GLToggle      *lowFluxToggle;
  GLButton    *applyButton;
  GLButton    *cancelButton;
  GLButton    *lowFluxInfo;
  GLButton    *newReflectmodeInfo;

  /*GLTextField *outgassingText;
  GLTextField *gasmassText;*/

};

#endif /* _GLOBALSETTINGSH_ */
