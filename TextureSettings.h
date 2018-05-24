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
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLGradient.h"
#include "GeometryViewer.h"

#ifndef _TEXTURESETTINGSH_
#define _TEXTURESETTINGSH_

class SynradGeometry;

class TextureSettings : public GLWindow {

public:

  // Construction
  TextureSettings();

  // Component methods
  void Display(Worker *w,GeometryViewer **v);
  void Update();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);

private:

  void UpdateSize();

  Worker         *worker;
  SynradGeometry *geom;
  GeometryViewer **viewers;

  GLToggle      *texAutoScale;
  GLTextField   *texMinText;
  GLTextField   *texMaxText;
  GLLabel       *texCMinText;
  GLLabel       *texCMaxText;
  GLToggle      *colormapBtn;
  GLTextField   *swapText;
  GLGradient    *gradient;
  GLToggle      *logBtn;
  GLButton      *setCurrentButton;
  GLCombo       *modeCombo;

  GLButton    *updateButton;

};

#endif /* _TEXTURESETTINGSH_ */
