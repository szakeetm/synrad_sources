/*
  File:        FacetCoordinates.h
  Description: Facet coordinates window
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
#include "GLApp/GLTextField.h"
#include "GLApp/GLList.h"
#include "GLApp/GLButton.h"
#include "Worker.h"

#ifndef _FACETCOORDINATESH_
#define _FACETCOORDINATESH_

class FacetCoordinates : public GLWindow {

public:

  // Construction
  FacetCoordinates();

  // Component method
  void Display(Worker *w);
  void Update();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);

private:

  void GetSelected();
  void setBounds(GLComponent *org,GLComponent *src,int x,int y,int w,int h);

  Worker       *worker;
  Facet        *selFacet;
  GLList       *facetListC;
  GLButton     *dismissButton;
  GLButton     *updateButton;
  GLButton     *insert1Button;
  GLButton     *insert2Button;
  GLButton     *removeButton;
  GLTextField  *insertPosText;

};

#endif /* _FACETCOORDINATESH_ */
