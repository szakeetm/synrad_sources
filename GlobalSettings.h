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

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLGradient.h"

#ifndef _GLOBALSETTINGSH_
#define _GLOBALSETTINGSH_

class GlobalSettings : public GLWindow {

public:

  // Construction
  GlobalSettings();
  void SMPUpdate(float appTime,BOOL forceUpdate=FALSE);

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

  GLToggle      *chkAntiAliasing;
  GLToggle      *chkWhiteBg;
  //GLToggle      *chkNonIsothermal;
  GLToggle      *chkSimuOnly;
  GLToggle      *chkCheckForUpdates;
  GLToggle      *chkAutoUpdateFormulas;
  GLToggle      *chkCompressSavedFiles;
  GLToggle      *lowFluxToggle;
  GLButton    *applyButton;
  GLButton    *cancelButton;
  GLButton    *lowFluxInfo;
  


};

#endif /* _GLOBALSETTINGSH_ */
