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
#include "GLApp/GLChart/GLChart.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLParser.h"
#include "GLApp/GLTextField.h"
#include "Worker.h"
#include "Geometry_shared.h"

#ifndef _PROFILEPLOTTERH_
#define _PROFILEPLOTTERH_

class ProfilePlotter : public GLWindow {

public:

  // Construction
  ProfilePlotter();

  // Component method
  void Display(Worker *w);
  void Refresh();
  void Update(float appTime,bool force=false);
  void Reset();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void SetBounds(int x,int y,int w,int h);
  void addView(int facet);
  std::vector<int> GetViews();
  void SetViews(std::vector<int> views,int mode);
  bool IsLogScaled();
  void SetLogScaled(bool logScale);
  void SetWorker(Worker *w);

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

#endif /* _PROFILEPLOTTERH_ */
