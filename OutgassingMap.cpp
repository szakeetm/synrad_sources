/*
  File:        OutgassingMap.cpp
  Description: Outgassing Map dialog
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

#include "OutgassingMap.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLFileBox.h"
#include "SynRad.h"
extern SynRad *theApp;

//static const char *fileFilters = "Text files\0*.txt";
//static const int   nbFilter = sizeof(fileFilters) / (2*sizeof(char *));

// --------------------------------------------------------------------

OutgassingMap::OutgassingMap():GLWindow() {

  int wD = 600;
  int hD = 300;
  lastUpdate = 0.0f;
  strcpy(currentDir,".");

  SetTitle("Outgassing map");
  SetResizable(TRUE);
  SetIconfiable(TRUE);
  SetMinimumSize(wD,hD);

  mapList = new GLList(0);
  mapList->SetColumnLabelVisible(TRUE);
  mapList->SetRowLabelVisible(TRUE);
  mapList->SetAutoColumnLabel(TRUE);
  mapList->SetAutoRowLabel(TRUE);
  //mapList->SetRowLabelMargin(20);
  mapList->SetGrid(TRUE);
  mapList->SetSelectionMode(SINGLE_CELL);
  mapList->SetCornerLabel("\202\\\201");
  Add(mapList);

  
  desLabel = new GLLabel("Desorption type:");
  Add(desLabel);
  desCombo = new GLCombo(0);
  desCombo->SetSize(3);
  desCombo->SetValueAt(0,"Uniform");
  desCombo->SetValueAt(1,"Cosine");
  desCombo->SetValueAt(2,"Cosine^N");
  
  desCombo->SetSelectedIndex(1);
  Add(desCombo);
  

  /*
  loadButton = new GLButton(0,"Load");
  Add(loadButton);
  */

  sizeButton = new GLButton(0,"Auto size");
  Add(sizeButton);

  explodeButton = new GLButton(0,"Explode");
  Add(explodeButton);

  pasteButton = new GLButton(0,"Paste");
  Add(pasteButton);

  /*
  maxButton = new GLButton(0,"Find Max.");
  Add(maxButton);
  */

  cancelButton = new GLButton(0,"Dismiss");
  Add(cancelButton);

  exponentText = new GLTextField(0,"");
  exponentText->SetEditable(FALSE);
  Add(exponentText);

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  RestoreDeviceObjects();

  worker = NULL;

}

// --------------------------------------------------------------------

void OutgassingMap::PlaceComponents() {

  mapList->SetBounds(5,5,width-15,height-55);
  //saveButton->SetBounds(10,height-45,70,19);
  explodeButton->SetBounds(10,height-45,70,19);
  sizeButton->SetBounds(90,height-45,70,19);
  pasteButton->SetBounds(170,height-45,70,19);
  //maxButton->SetBounds(170,height-45,70,19);
  desLabel->SetBounds(270,height-45,30,19);
  desCombo->SetBounds(370,height-45,80,19);
  cancelButton->SetBounds(width-90,height-45,80,19);
  exponentText->SetBounds(450,height-45,50,19);
}

// -----------------------------------------------------------------

void OutgassingMap::SetBounds(int x,int y,int w,int h) {

  GLWindow::SetBounds(x,y,w,h);
  PlaceComponents();

}

// --------------------------------------------------------------------

void OutgassingMap::GetSelected() {

  if(!worker) return;

  Geometry *geom = worker->GetGeometry();
  selFacet = NULL;
  int i = 0;
  int nb = geom->GetNbFacet();
  while(!selFacet && i<nb) {
    if( geom->GetFacet(i)->selected ) selFacet = geom->GetFacet(i);
    if(!selFacet) i++;
  }

  char tmp[32];
  sprintf(tmp,"Outgassing map for Facet #%d",i+1);
  SetTitle(tmp);

}

// --------------------------------------------------------------------

void OutgassingMap::Update(float appTime,BOOL force) {

  if(!IsVisible()) return;

  if(force) {
    UpdateTable();
    lastUpdate = appTime;
    return;
  }

  if( (appTime-lastUpdate>1.0f) ) {
    if(worker->running) UpdateTable();
    lastUpdate = appTime;
  }

}

// --------------------------------------------------------------------

void OutgassingMap::UpdateTable() {
	//maxValue=0.0f;
	//double scale;
  GetSelected();
  if( !selFacet ) {
    mapList->Clear();
    return;
  }

  SHELEM *mesh = selFacet->mesh;
  if( mesh ) {

    char tmp[256];
    int w = selFacet->sh.texWidth;
    int h = selFacet->sh.texHeight;
    mapList->SetSize(w,h);
    mapList->SetColumnAlign(ALIGN_CENTER);
	mapList->SetGrid(TRUE);
	mapList->SetColumnLabelVisible(TRUE);
	//mapList->SetSelectionMode(SINGLE_ROW);

    //int mode = viewCombo->GetSelectedIndex();

    
        for(int i=0;i<w;i++) { //width (cols)
			*(mapList->cEdits+i)=EDIT_NUMBER;
          for(int j=0;j<h;j++) { //height (rows)
			  sprintf(tmp,"0");
			  if (selFacet->mesh[i+j*w].area==0.0) sprintf(tmp,"Outside");
			  mapList->SetValueAt(i,j,tmp);
          }
        }

      
    }
  else {
	  mapList->Clear();
	  mapList->SetSize(0,0);
  }
  

}

// --------------------------------------------------------------------

void OutgassingMap::Display(Worker *w) {

  worker = w;
  UpdateTable();
  SetVisible(TRUE);

}

// --------------------------------------------------------------------

void OutgassingMap::Close() {
  worker = NULL;
  if(selFacet) selFacet->UnselectElem();
  mapList->Clear();
}

// --------------------------------------------------------------------
/*
void OutgassingMap::SaveFile() {

  if(!selFacet) return;

  FILENAME *fn = GLFileBox::SaveFile(currentDir,NULL,"Save File",fileFilters,nbFilter);

  if( fn ) {

    int u,v,wu,wv;
    if( !mapList->GetSelectionBox(&u,&v,&wu,&wv) ) {
      u=0;
      v=0;
      wu = mapList->GetNbRow();
      wv = mapList->GetNbColumn();
    }

    // Save tab separated text
    FILE *f = fopen(fn->fullName,"w");

    if( f==NULL ) {
      char errMsg[512];
      sprintf(errMsg,"Cannot open file\nFile:%s",fn->fullName);
      GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
      return;
    }

    for(int i=u;i<u+wu;i++) {
      for(int j=v;j<v+wv;j++) {
        char *str = mapList->GetValueAt(j,i);
        if( str ) fprintf(f,"%s",str);
        if( j<v+wv-1 ) 
          fprintf(f,"\t");
      }
      fprintf(f,"\r\n");
    }
    fclose(f);

  }

}
*/
// --------------------------------------------------------------------

void OutgassingMap::ProcessMessage(GLComponent *src,int message) {

  switch(message) {

    case MSG_CLOSE:
      Close();
      break;

    case MSG_BUTTON:
      if(src==cancelButton) {
        Close();
        GLWindow::ProcessMessage(NULL,MSG_CLOSE);
      } else if(src==pasteButton) {
		  mapList->PasteClipboardText();
      } else if (src==sizeButton) {
        mapList->AutoSizeColumn();
      } else if (src==explodeButton) {
		  if (worker->GetGeometry()->GetNbSelected()!=1) {
			  GLMessageBox::Display("Exactly one facet has to be selected","Error",GLDLG_OK,GLDLG_ICONERROR);
			  return;
		  }
		  if (!selFacet->hasMesh) {
			  GLMessageBox::Display("Selected facet must have a mesh","Error",GLDLG_OK,GLDLG_ICONERROR);
			  return;
		  }
		  
		  int w = mapList->GetNbColumn(); //width
		  int h = mapList->GetNbRow(); //height

		  double *values = (double *)malloc(w*h*sizeof(double));
		  int count=0;
		  for(int j=0;j<h;j++) { //height (rows)
			  for(int i=0;i<w;i++) { //width (cols)
				  if (selFacet->mesh[i+j*w].area>0.0) {
					  char *str = mapList->GetValueAt(i,j);
					  int conv = sscanf(str,"%lf",values+count);
					  if(!conv || !(*(values+count++)>=0.0)) {
						  mapList->SetSelectedCell(i,j);
						  char tmp[256];
						  sprintf(tmp,"Invalid outgassing number at Cell(%d,%d)",i,j);
						  GLMessageBox::Display(tmp,"Error",GLDLG_OK,GLDLG_ICONERROR);
						  return;
					  }
				  }
			  }
		  }
		  double desorbTypeN;
		  if (desCombo->GetSelectedIndex()==2) {
			  exponentText->GetNumber(&desorbTypeN) ;
				  if( !(desorbTypeN>1.0) ) {
					  exponentText->SetFocus(TRUE);
					  GLMessageBox::Display("Desorption type exponent must be greater than 1.0","Error",GLDLG_OK,GLDLG_ICONERROR);
					  return;
				  }
			  
		  }
		  if( GLMessageBox::Display("Explode selected facet?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {

			  if (mApp->AskToReset()) {
				  mApp->changedSinceSave=TRUE;
				  try { 
					  worker->GetGeometry()->ExplodeSelected(TRUE,desCombo->GetSelectedIndex(),desorbTypeN,values);
					  SAFE_FREE(values);
					  mApp->UpdateModelParams();
					  mApp->UpdateFacetParams(TRUE);
					  worker->GetGeometry()->CalcTotalOutGassing();
					  // Send to sub process
					  worker->Reload(); } catch(Error &e) {
						  GLMessageBox::Display((char *)e.GetMsg(),"Error exploding facet (not enough memory?)",GLDLG_OK,GLDLG_ICONERROR);
					}
				}
		  }
      }/*else if (src==saveButton) {
        SaveFile();
      } else if (src==maxButton) {
		         int u,v,wu,wv;
		  mapList->SetSelectedCell(maxX,maxY);
		  if( mapList->GetSelectionBox(&v,&u,&wv,&wu) )
          selFacet->SelectElem(u,v,wu,wv);
      }
	  */
      break;

    case MSG_LIST:
      if(src==mapList) {
        int u,v,wu,wv;
        if( mapList->GetSelectionBox(&v,&u,&wv,&wu) )
          selFacet->SelectElem(u,v,wu,wv);
      }
      break;

    case MSG_COMBO:
      if(src==desCombo) {
		  BOOL cosineN = desCombo->GetSelectedIndex()==2;
		  exponentText->SetEditable(cosineN);
		  if (!cosineN) exponentText->SetText("");
      }
      break;

  }

  GLWindow::ProcessMessage(src,message);
}
