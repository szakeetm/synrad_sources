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
#include "GLApp/GLWindow.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLProgress.h"
#include "Worker.h"

#ifndef _FACETMESHH_
#define _FACETMESHH_

class FacetMesh : public GLWindow {

public:

  // Construction
  FacetMesh();

  // Component method
  void EditFacet(Worker *w);

  // Implementation
  void ProcessMessage(GLComponent *src,int message);

private:

  bool Apply();
  void UpdateSize();
  void UpdateSizeForRatio();
  void UpdateToggle(GLComponent *src);
  void QuickApply(); //Apply View Settings without stopping the simulation

  Worker   *worker;
  Geometry *geom;
  int       fIdx;

  GLTitledPanel *iPanel;
  GLTitledPanel *vPanel;
  GLTextField   *uLength;
  GLTextField   *vLength;

  GLToggle      *enableBtn;
  GLToggle      *boundaryBtn;
  GLToggle      *recordAbsBtn;
  GLToggle      *recordReflBtn;
  GLToggle      *recordTransBtn;
  GLToggle      *recordDirBtn;
  GLToggle      *showTexture;
  GLToggle      *showVolume;
  GLTextField   *resolutionText;

  GLTextField   *ramText;
  GLTextField   *cellText;
  GLButton      *updateButton;

  GLButton    *applyButton;
  GLButton	  *quickApply; //Apply View Settings without stopping the simulation
  GLButton    *cancelButton;

  GLProgress  *progressDlg;

};

#endif /* _FACETMESHH_ */
