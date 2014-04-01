/*
  File:        ProfilePlotter.h
  Description: Profile plotter window
  Program:     SynRad
  Author:      R. KERSEVAN / M ADY
  Copyright:   E.S.R.F / CERN

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
#include "GLApp/GLChart.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLParser.h"
#include "GLApp/GLTextField.h"
#include "Worker.h"
#include "Geometry.h"

#ifndef _FACETPLOTTERH_
#define _FACETPLOTTERH_

class ProfilePlotter : public GLWindow {

public:

  // Construction
  ProfilePlotter();

  // Component method
  void Display(Worker *w);
  void Refresh();
  void Update(float appTime,BOOL force=FALSE);
  void Reset();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void SetBounds(int x,int y,int w,int h);

private:

  void addView(int facet,int mode);
  void remView(int facet,int mode);
  void refreshViews();
  void plot();

  Worker      *worker;
  GLButton    *dismissButton;
  GLChart     *chart;
  GLCombo     *profCombo;
  GLToggle    *logToggle;
  GLToggle    *normToggle;
  GLButton    *selButton;
  GLButton    *addButton;
  GLButton    *removeButton;
  GLButton    *resetButton;
  GLTextField *formulaText;
  GLButton    *formulaBtn;

  GLDataView  *views[32];
  GLCColor    *colors[8];

  int          nbColors;
  int          nbView;
  float        lastUpdate;

};

#endif /* _FACETPLOTTERH_ */
