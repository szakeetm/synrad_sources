/*
  File:        SpectrumPlotter.h
  Description: Spectrum plotter window
  Program:     SynRad
  Author:      R. KERSEVAN / M SZAKACS
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

#ifndef _SPECTRUMPLOTTERH_
#define _SPECTRUMPLOTTERH_

class SpectrumPlotter : public GLWindow {

public:

  // Construction
  SpectrumPlotter();

  // Component method
  void Display(Worker *w);
  void SetScale();
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

  Worker      *worker;
  GLButton    *dismissButton;
  GLChart     *chart;
  GLCombo     *specCombo;
  GLToggle    *logToggle;
  GLToggle    *normToggle;
  GLButton    *selButton;
  GLButton    *addButton;
  GLButton    *removeButton;
  GLButton    *resetButton;

  GLDataView  *views[32];

  int          nbView;
  double       delta;
  float        lastUpdate;

};

#endif /* _FACETPLOTTERH_ */
