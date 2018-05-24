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
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLTitledPanel.h"
#include "Geometry_shared.h"
#include "GeometryViewer.h"

#ifndef _VIEWER3DSETTINGSH_
#define _VIEWER3DSETTINGSH_

extern int whiteBg;

class Viewer3DSettings : public GLWindow {

public:

  // Construction
  Viewer3DSettings();

  // Component methods
  void Refresh(Geometry *s, GeometryViewer *v);
  void Reposition(int wD = 0, int hD = 0);

  // Implementation
  void ProcessMessage(GLComponent *src,int message);

private:

  Geometry       *geom;
  GeometryViewer *viewer;

  GLTitledPanel *panel;
  GLCombo       *showMode;
  GLTextField   *traStepText;
  GLTextField   *angStepText;
  GLTextField   *dispNumHits;
  GLTextField   *dispNumLeaks;
  GLTextField   *dispNumTraj;
  GLToggle      *hiddenEdge;
  GLToggle      *hiddenVertex;
  GLToggle      *showMesh;
  GLToggle      *showDirection;
  GLToggle      *showTPtoggle;
  GLToggle      *shadeLines;
  GLToggle      *autoNorme;
  GLToggle      *centerNorme;
  GLToggle		*antiAliasing;
  GLToggle		*bigDots;
  GLToggle		*hideLotselected;
  GLTextField   *normeText;
  GLTextField   *hideLotText;

  GLButton    *applyButton;
  GLButton    *cancelButton;

};

#endif /* _VIEWER3DSETTINGSH_ */

