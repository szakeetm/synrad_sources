/*
  File:        SelectDialog.cpp
  Description: Facet selection by number dialog

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "SelectDialog.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "SynRad.h"
extern SynRad *mApp;

// Construct a message dialog box
SelectDialog::SelectDialog(Worker *w):GLWindow() {

  int xD,yD,wD,hD;
  wD=340;hD=80;
  SetTitle("Select facet by number");

  // Label
  GLLabel  *label = new GLLabel("Facet number:");
  label->SetBounds(10,3,100,18);
  Add(label);

  numText = new GLTextField(0,NULL);
  numText->SetBounds(90,3,100,18);
  Add(numText);

  // Buttons
  int startX = 5;
  int startY = 30;
    
     
    GLButton *btn = new GLButton(GLDLG_SELECT,"Select");
    btn->SetBounds(startX,startY,95,20);
    Add(btn);
    startX+=100;
  

	 
    GLButton *btn2 = new GLButton(GLDLG_SELECT_ADD,"Add to selection");
    btn2->SetBounds(startX,startY,95,20);
    Add(btn2);
    startX+=100;
  

  
    GLButton *btn3 = new GLButton(GLDLG_SELECT_REM,"Remove from selection");
    btn3->SetBounds(startX,startY,125,20);
    Add(btn3);
    startX+=130;  

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  if( wD > wS ) wD = wS;
  xD = (wS-wD)/2;
  yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  // Create objects
  RestoreDeviceObjects();

  rCode = GLDLG_CANCEL;

  work = w;
  geom = w->GetGeometry();
}

// -------------------------------------------------

void SelectDialog::ProcessMessage(GLComponent *src,int message) {
  if(message==MSG_BUTTON) {
		rCode = src->GetId();
		  int facetnumber;
		  numText->GetNumberInt(&facetnumber);
		  if( !facetnumber ) return;	
		  if((facetnumber<=0 )||(facetnumber>geom->GetNbFacet())) {
			  GLMessageBox::Display("Invalid number","Error",GLDLG_OK,GLDLG_ICONERROR);
			  return;
		  }
		  if (rCode == GLDLG_SELECT) geom->UnSelectAll();
		  geom->GetFacet(facetnumber-1)->selected = (rCode == GLDLG_SELECT || rCode == GLDLG_SELECT_ADD);
		  geom->UpdateSelection();
		  mApp->UpdateFacetParams(TRUE);
		  mApp->UpdateFacetlistSelected();
		  mApp->facetList->ScrollToVisible(facetnumber-1,1,TRUE);
  }
  GLWindow::ProcessMessage(src,message);
}