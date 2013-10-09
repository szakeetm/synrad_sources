/*
  File:        FacetCoordinates.cpp
  Description: Facet coordinates window
  Program:     SynRad
  Author:      R. Kerservan / M SZAKACS
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

#include "FacetCoordinates.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "Utils.h"
#include "SynRad.h"

//extern GLApplication *theApp;
extern SynRad *theApp;

static const int   flWidth[] = {35,40,100,100,100};
static const char *flName[] = {"#","Vertex","X","Y","Z"};
static const int   flAligns[] = { ALIGN_CENTER,ALIGN_CENTER,ALIGN_LEFT,ALIGN_LEFT,ALIGN_LEFT };
static const int   fEdits[] = { 0,0,EDIT_NUMBER,EDIT_NUMBER,EDIT_NUMBER };

FacetCoordinates::FacetCoordinates():GLWindow() {

  int wD = 405;
  int hD = 350;
  SetIconfiable(TRUE);
  SetResizable(FALSE);

  GLTitledPanel *p = new GLTitledPanel("Facet index");
  p->SetBounds(5,hD-145,wD-10,95);
  Add(p);

  insert1Button = new GLButton(0,"Insert (New vertex)");
  setBounds(p,insert1Button,5,20,120,19);
  Add(insert1Button);

  insert2Button = new GLButton(0,"Insert (Existing vertex)");
  setBounds(p,insert2Button,5,45,120,19);
  Add(insert2Button);

  removeButton = new GLButton(0,"Remove");;
  setBounds(p,removeButton,5,70,120,19);
  Add(removeButton);

  GLLabel *l1 = new GLLabel("Insertion position");
  setBounds(p,l1,130,20,90,19);
  Add(l1);

  insertPosText = new GLTextField(0,"");
  setBounds(p,insertPosText,220,20,30,19);
  Add(insertPosText);

  // TODO
  insert1Button->SetEnabled(FALSE);
  insert2Button->SetEnabled(FALSE);
  removeButton->SetEnabled(FALSE);
  insertPosText->SetEnabled(FALSE);

  updateButton = new GLButton(0,"Apply");
  updateButton->SetBounds(wD-195,hD-43,90,19);
  Add(updateButton);

  dismissButton = new GLButton(0,"Dismiss");
  dismissButton->SetBounds(wD-95,hD-43,90,19);
  Add(dismissButton);

  facetListC = new GLList(0);
  facetListC->SetBounds(5,5,wD-10,hD-60);
  facetListC->SetColumnLabelVisible(TRUE);
  facetListC->SetGrid(TRUE);
  //facetListC->SetSelectionMode(SINGLE_CELL);
  //facetListC->Sortable=TRUE;
  Add(facetListC);

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  RestoreDeviceObjects();

}

void FacetCoordinates::GetSelected() {

  selFacet = NULL;
  int i = 0;
  Geometry *g = worker->GetGeometry();
  int nb = g->GetNbFacet();
  while(!selFacet && i<nb) {
    if( g->GetFacet(i)->selected ) selFacet = g->GetFacet(i);
    if(!selFacet) i++;
  }

  char tmp[32];
  sprintf(tmp,"Facets coordinates #%d",i+1);
  SetTitle(tmp);

}

void FacetCoordinates::Update() {

  char tmp[256];
  if(!IsVisible()) return;
  GetSelected();
  if(!selFacet) return;

  Geometry *s = worker->GetGeometry();

  int nbIndex = selFacet->sh.nbIndex;

  facetListC->SetSize(5,nbIndex);
  facetListC->SetColumnWidths((int*)flWidth);
  facetListC->SetColumnLabels((char **)flName);
  facetListC->SetColumnAligns((int *)flAligns);
  facetListC->SetColumnEditable((int *)fEdits);
  for(int i=0;i<nbIndex;i++) {
    int idx = selFacet->indices[i];
    sprintf(tmp,"%d",i+1);
    facetListC->SetValueAt(0,i,tmp);
    sprintf(tmp,"%d",idx+1);
    facetListC->SetValueAt(1,i,tmp);
    VERTEX3D *v = s->GetVertex(idx);
    sprintf(tmp,"%.10g",v->x);
    facetListC->SetValueAt(2,i,tmp);
    sprintf(tmp,"%.10g",v->y);
    facetListC->SetValueAt(3,i,tmp);
    sprintf(tmp,"%.10g",v->z);
    facetListC->SetValueAt(4,i,tmp);
  }

}

void FacetCoordinates::Display(Worker *w) {

  worker = w;
  SetVisible(TRUE);
  Update();

}

void FacetCoordinates::ProcessMessage(GLComponent *src,int message) {

  Geometry *geom = worker->GetGeometry();
  SynRad *mApp = (SynRad *)theApp;
  switch(message) {
    case MSG_BUTTON:
      if(src==dismissButton) {
        SetVisible(FALSE);
      } else if(src==updateButton) {
        int rep = GLMessageBox::Display("Apply geometry changes ?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);
        if( rep == GLDLG_OK ) {
			if (mApp->AskToReset(worker)) {
			//if (worker->running) worker->Stop_Public();
			changedSinceSave=TRUE;
          for(int i=0;i<selFacet->sh.nbIndex;i++) {
            int idx = selFacet->indices[i];
            double x,y,z;
            sscanf(facetListC->GetValueAt(2,i),"%lf",&x);
            sscanf(facetListC->GetValueAt(3,i),"%lf",&y);
            sscanf(facetListC->GetValueAt(4,i),"%lf",&z);
            geom->SetVertex(idx,x,y,z);
          }
          geom->Rebuild();
          // Send to sub process
          try {
            worker->Reload();
          } catch(Error &e) {
            GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
            return;
          }
          GLWindowManager::FullRepaint();
		 }
        }
      }
    break;
    case MSG_LIST:
      if(src==facetListC) {
        char tmp[32];
        sprintf(tmp,"%d",facetListC->GetSelectedRow()+1);
        insertPosText->SetText(tmp);
      }
    break;

  }

  GLWindow::ProcessMessage(src,message);

}

void FacetCoordinates::setBounds(GLComponent *org,GLComponent *src,int x,int y,int w,int h) {
  int xc,yc,wc,hc;
  org->GetBounds(&xc,&yc,&wc,&hc);
  src->SetBounds(xc+x,yc+y,w,h);
}
