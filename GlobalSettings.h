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
#ifndef _GLOBALSETTINGSH_
#define _GLOBALSETTINGSH_

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLGradient.h"
#include "Buffer_shared.h" //MAX_PROCESS macro

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
 
  int lastUpdate;
  //float lastCPUTime[MAX_PROCESS];
  //float lastCPULoad[MAX_PROCESS];

  GLToggle      *chkAntiAliasing;
  GLToggle      *chkWhiteBg;
  GLToggle		*leftHandedToggle;
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
