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
#include "GLApp/GLList.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLCombo.h"
#include "Worker.h"

#ifndef _TRAJECTORYDETAILSH_
#define _TRAJECTORYDETAILSH_

#define NB_TPCOLUMN 53

class Facet;

class TrajectoryDetails : public GLWindow {

public:

  // Construction
  TrajectoryDetails();

  // Component method
  void Display(Worker *w,int regionId);
  void Update();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void SetBounds(int x,int y,int w,int h);
  void SelectPoint(int idx);
  size_t GetRegionId();

private:

  char *GetCountStr(Facet *f);
  void UpdateTable();
  char *FormatCell(int idx,int mode,GenPhoton* photon);
  void PlaceComponents();

  size_t displayedRegion;
  int freq;
  Worker      *worker;
  GLList      *pointList;
  GLCombo     *regionCombo;
  GLTextField *freqText;
  GLLabel     *everyLabel;
  GLLabel     *nbPointLabel;

  GLTitledPanel *sPanel;          // Toggle panel
  GLToggle      *show[NB_TPCOLUMN-1];
  int            shown[NB_TPCOLUMN];

  GLButton    *checkAllButton;
  GLButton    *uncheckAllButton;
  GLButton    *dismissButton;
  GLButton	  *updateButton;

};

#endif /* _TRAJECTORYDETAILSH_ */
