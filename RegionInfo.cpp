/*
  File:        RegionInfo.cpp
  Description: Region info dialog
  Program:     SynRad

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "RegionInfo.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "Region_full.h"
#include "SynRad.h"
#include "GLApp/GLFileBox.h"

extern GLApplication *theApp;
static const char *fileTrajFilters = "CSV files\0*.csv\tTXT files\0*.txt\0All files\0*.*\0";
static const int   nbTrajFilter = 3;
bool EndsWithParam(const char* s);

RegionInfo::RegionInfo(Worker *w) {

  int wD = 550;
  int hD = 240;

  SetTitle("Region info");

  GLLabel *l1 = new GLLabel("Number of points:");
  l1->SetBounds(10,35,170,18);
  Add(l1);

  t1 = new GLLabel("No trajectory loaded.");
  t1->SetBounds(180,35,170,18);
  Add(t1);

  GLLabel *l2 = new GLLabel("Start:");
  l2->SetBounds(10,60,170,18);
  Add(l2);
  
  t2 = new GLLabel("No trajectory loaded.");
  t2->SetBounds(180,60,170,18);
  Add(t2);

  GLLabel *l3 = new GLLabel("End:");
  l3->SetBounds(10,85,170,18);
  Add(l3);

  t3 = new GLLabel("No trajectory loaded.");
  t3->SetBounds(180,85,170,18);
  Add(t3);

   GLLabel *l4 = new GLLabel("Selected Point:");
  l4->SetBounds(10,110,170,18);
  Add(l4);

  t4 = new GLLabel("None");
  t4->SetBounds(180,110,170,18);
  Add(t4);

  GLLabel *l5 = new GLLabel("Magnetic field:");
  l5->SetBounds(10,135,170,18);
  Add(l5);

  t5 = new GLLabel("");
  t5->SetBounds(180,135,170,18);
  Add(t5);

  GLLabel *l6 = new GLLabel("Curvature:");
  l6->SetBounds(10,160,170,18);
  Add(l6);

  t6 = new GLLabel("");
  t6->SetBounds(180,160,170,18);
  Add(t6);

  /*exportPanel = new GLTitledPanel("Export region points");
  exportPanel->SetBounds(5,180,wD-10,65);
  Add(exportPanel);

  GLLabel *freqTxt = new GLLabel("Export every          th point");
  freqTxt->SetBounds(10,195,150,18);
  Add(freqTxt);

  freqField = new GLTextField(0,"1");
  freqField->SetBounds(78,195,30,18);
  Add(freqField);

  freqLabel = new GLLabel("");
  freqLabel->SetBounds(180,195,150,18);
  Add(freqLabel);

  integrateToggle = new GLToggle(0,"Calculate integrals (go through all trajectory points)");
  integrateToggle->SetBounds(10,220,200,18);
  integrateToggle->SetState(true);
  Add(integrateToggle);*/

  viewPointsButton = new GLButton(0,"View Points");
  viewPointsButton->SetBounds(wD-540,hD-44,85,21);
  Add(viewPointsButton);

  saveAsButton = new GLButton(0,"Save as...");
  saveAsButton->SetBounds(wD-450,hD-44,85,21);
  Add(saveAsButton);
  
  editButton = new GLButton(0,"Edit");
  editButton->SetBounds(wD-360,hD-44,85,21);
  Add(editButton);

  notepadButton = new GLButton(0,"to Notepad");
  notepadButton->SetBounds(wD-270,hD-44,85,21);
  Add(notepadButton);
  
  reloadButton = new GLButton(0,"Reload file");
  reloadButton->SetBounds(wD-180,hD-44,85,21);
  Add(reloadButton);
  
  cancelButton = new GLButton(0,"Dismiss");
  cancelButton->SetBounds(wD-90,hD-44,85,21);
  Add(cancelButton);

  

  pathLabel=new GLTextField(0,"No PAR file loaded");
  pathLabel->SetEditable(false);
  pathLabel->SetBounds(100,10,430,18);
  Add(pathLabel);

  regionSelector=new GLCombo(0);
  regionSelector->SetBounds(10,10,80,20);
  regionSelector->SetSelectedIndex(0);
  Add(regionSelector);
  work=w;
  Update();

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);
  RestoreDeviceObjects();
}

void RegionInfo::ProcessMessage(GLComponent *src,int message) {
	SynRad *mApp = (SynRad *)theApp;
  //double dX,dY,dZ;

  switch(message) {
  /* 
  case MSG_TEXT_UPD:
		if( freqField->GetNumberInt(&exportFrequency) ) {
			char tmp[256];
			sprintf(tmp,"%d points will be exported with an interval of %g cm.",
				(int)((double)selectedRegion->Points.size()/(double)exportFrequency),selectedRegion->params.dL_cm*(double)exportFrequency);
			freqLabel->SetText(tmp);
			integrateToggle->SetState(exportFrequency==1);
		}
		*/
	case MSG_BUTTON:

    if(src==cancelButton) 
      GLWindow::ProcessMessage(NULL,MSG_CLOSE);
	else if (src==notepadButton) {
		char tmp[512];
		sprintf(tmp,"notepad.exe \"%s\"",work->regions[regionSelector->GetSelectedIndex()].fileName.c_str());
		StartProc(tmp,STARTPROC_FOREGROUND);
	} else if (src==reloadButton) {
		char tmp[512];
		sprintf(tmp,"%s",work->regions[regionSelector->GetSelectedIndex()].fileName.c_str());
		try {
			work->AddRegion(tmp,regionSelector->GetSelectedIndex());
		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),tmp);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
		Update();
		if (mApp->trajectoryDetails && mApp->trajectoryDetails->IsVisible() && mApp->trajectoryDetails->GetRegionId() == regionSelector->GetSelectedIndex()) mApp->trajectoryDetails->Update();
		if (mApp->spectrumPlotter) mApp->spectrumPlotter->SetScale();
	} else if (src==editButton) {
		if( mApp->regionEditor==NULL ) mApp->regionEditor = new RegionEditor();
		mApp->regionEditor->Display(work,regionSelector->GetSelectedIndex());
		//regionEditor->DoModal();
		//SAFE_DELETE(regionEditor);
	} else if (src==viewPointsButton) {
		if ( mApp->trajectoryDetails==NULL) mApp->trajectoryDetails = new TrajectoryDetails();
		mApp->trajectoryDetails->Display(work,regionSelector->GetSelectedIndex());
	} else if (src==saveAsButton) {
		FILENAME *fn = GLFileBox::SaveFile(NULL,NULL,"Save Region","param files\0*.param\0All files\0*.*\0",2);
		if (!fn || !fn->fullName) return;
		if (!EndsWithParam(fn->fullName))
			sprintf(fn->fullName,"%s.param",fn->fullName); //append .param extension
		work->SaveRegion(fn->fullName,regionSelector->GetSelectedIndex());
	}
	break;
	case MSG_COMBO:
    if(src==regionSelector) {

      Update();

    }
    break;
	case MSG_TOGGLE:
	/*if (src==integrateToggle) {
		calcIntegrals=integrateToggle->GetState();
	}*/
	break;
  }
  GLWindow::ProcessMessage(src,message);
}

void RegionInfo::Update() {
	
	int nbRegion=(int)work->regions.size();
	regionSelector->SetSize(nbRegion);
	if (nbRegion==0) return;
	char tmp2[16];
	for (int i=0;i<nbRegion;i++) {
		sprintf(tmp2,"Region %d",i+1);
		regionSelector->SetValueAt(i,tmp2);
	}
	if (regionSelector->GetSelectedIndex()>(nbRegion-1)) regionSelector->SetSelectedIndex(nbRegion-1);
	regionSelector->SetSelectedIndex(regionSelector->GetSelectedIndex()); //update text
	selectedRegion=&(work->regions[regionSelector->GetSelectedIndex()]);
	if (!selectedRegion->isLoaded) return;

	char tmp[256];
	/*sprintf(tmp,"%d points will be exported with an interval of %g cm.",
		(int)((double)selectedRegion->Points.size()/(double)exportFrequency),selectedRegion->dL*(double)exportFrequency);
	freqLabel->SetText(tmp);*/

	SynRad *mApp = (SynRad *)theApp;
	
	
	sprintf(tmp,"%s",work->regions[regionSelector->GetSelectedIndex()].fileName.c_str());
	pathLabel->SetText(tmp);
		
	sprintf(tmp,"%d",(int)selectedRegion->Points.size());
	t1->SetText(tmp);
	sprintf(tmp,"(%g , %g , %g)",selectedRegion->Points[0].position.x,selectedRegion->Points[0].position.y,selectedRegion->Points[0].position.z);
	t2->SetText(tmp);
	sprintf(tmp,"(%g , %g , %g)",selectedRegion->Points[(int)selectedRegion->Points.size()-1].position.x,selectedRegion->Points[(int)selectedRegion->Points.size()-1].position.y,selectedRegion->Points[(int)selectedRegion->Points.size()-1].position.z);
	t3->SetText(tmp);
	if (selectedRegion->selectedPointId!=-1) {
		sprintf(tmp,"#%d (%g , %g , %g)",selectedRegion->selectedPointId,selectedRegion->Points[selectedRegion->selectedPointId].position.x,
			selectedRegion->Points[selectedRegion->selectedPointId].position.y,selectedRegion->Points[selectedRegion->selectedPointId].position.z);
	} else sprintf(tmp,"None.");
	t4->SetText(tmp);
	
	if (selectedRegion->selectedPointId!=-1) {
		Vector3d B=selectedRegion->B(selectedRegion->selectedPointId,Vector3d(0,0,0));
		sprintf(tmp,"B=%g T (%g , %g , %g)",B.Norme(),B.x,B.y,B.z);
	} else sprintf(tmp,"");
	t5->SetText(tmp);
	
	if (selectedRegion->selectedPointId!=-1) {
		Vector3d Rho=selectedRegion->Points[selectedRegion->selectedPointId].rho;
		sprintf(tmp,"Rho=%g cm (%g , %g , %g)",Rho.Norme(),Rho.x,Rho.y,Rho.z);
	} else sprintf(tmp,"");
	t6->SetText(tmp);
	return;
}

bool EndsWithParam(const char* s)
{
  int ret = 0;

  if (s != NULL)
  {
    size_t size = strlen(s);

    if (size >= 6 &&
        s[size-6] == '.' &&
        s[size-5] == 'p' &&
        s[size-4] == 'a' &&
        s[size-3] == 'r' &&
        s[size-2] == 'a' &&
        s[size-1] == 'm')
    {
      ret = 1;
    }
  }

  return ret;
}