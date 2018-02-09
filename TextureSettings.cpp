/*
  File:        TextureSettings.cpp
  Description: Texture settings dialog (min,max,autoscale,gradient)
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

#include "TextureSettings.h"
#include "SynradGeometry.h"
#include "Worker.h"
#include "Facet_shared.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/MathTools.h" //FormatMemory
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLCombo.h"

TextureSettings::TextureSettings():GLWindow() {

  int wD = 450;
  int hD = 230;

  SetTitle("Texture Scaling");
  SetIconfiable(true);

  GLTitledPanel *panel = new GLTitledPanel("Texture Range");
  panel->SetBounds(5,2,315,98);
  Add(panel);

  GLLabel *l1 = new GLLabel("Min");
  l1->SetBounds(10,20,50,18);
  Add(l1);

  texMinText = new GLTextField(0,"");
  texMinText->SetBounds(60,20,85,19);
  texMinText->SetEditable(true);
  Add(texMinText);

  GLLabel *l2 = new GLLabel("Max");
  l2->SetBounds(10,45,50,18);
  Add(l2);

  texMaxText = new GLTextField(0,"");
  texMaxText->SetBounds(60,45,85,19);
  texMaxText->SetEditable(true);
  Add(texMaxText);

  setCurrentButton = new GLButton(0,"Set to current");
  setCurrentButton->SetBounds(105,70,90,19);
  Add(setCurrentButton);

  updateButton = new GLButton(0,"Apply");
  updateButton->SetBounds(10,70,90,19);
  Add(updateButton);
  texAutoScale = new GLToggle(0,"Autoscale");
  texAutoScale->SetBounds(150,20,80,19);
  Add(texAutoScale);

  colormapBtn = new GLToggle(0,"Use colors");
  colormapBtn->SetBounds(230,20,85,19);
  Add(colormapBtn);

  logBtn = new GLToggle(0,"Log scale");
  logBtn->SetBounds(150,45,80,19);
  Add(logBtn);

  GLLabel *l3 = new GLLabel("Swap");
  l3->SetBounds(230,45,30,18);
  Add(l3);

  swapText = new GLTextField(0,"");
  swapText->SetEditable(false);
  swapText->SetBounds(260,45,50,18);
  Add(swapText);

  

  GLTitledPanel *panel2 = new GLTitledPanel("Current");
  panel2->SetBounds(325,2,120,98);
  Add(panel2);

  GLLabel *l4 = new GLLabel("Min:");
  l4->SetBounds(340,30,20,19);
  Add(l4);

  texCMinText = new GLLabel("");
  texCMinText->SetBounds(370,30,70,19);
  Add(texCMinText);

  GLLabel *l5 = new GLLabel("Max:");
  l5->SetBounds(340,65,20,19);
  Add(l5);

  texCMaxText = new GLLabel("");
  texCMaxText->SetBounds(370,65,70,19);
  Add(texCMaxText);

  

  GLTitledPanel *panel3 = new GLTitledPanel("Gradient");
  panel3->SetBounds(5,102,440,65);
  Add(panel3);

  gradient = new GLGradient(0);
  gradient->SetMouseCursor(true);
  gradient->SetBounds(10,117,420,40);
  Add(gradient);

  modeCombo = new GLCombo(0);
  modeCombo->SetSize(3);
  modeCombo->SetValueAt(0,"MC Hits");
  modeCombo->SetValueAt(1,"Flux (ph/sec/cm\262");
  modeCombo->SetValueAt(2,"Power (W/cm\262");
  modeCombo->SetBounds(10,180,100,25);
  modeCombo->SetSelectedIndex(1);
  Add(modeCombo);

  SetBounds(8,30,wD,hD);

  /*// Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);*/

  RestoreDeviceObjects();

  geom = NULL;

}

void TextureSettings::UpdateSize() {

  size_t swap = 0;
  size_t nbFacet = geom->GetNbFacet();
  for(size_t i=0;i<nbFacet;i++) {
    Facet *f = geom->GetFacet(i);
    if(f->sh.isTextured) {
      swap += f->GetTexSwapSize(colormapBtn->GetState());
    }
  }
  swapText->SetText(FormatMemory(swap));

}

void TextureSettings::Update() {

  if(!IsVisible() || IsIconic()) return;  

  char tmp[128];
  
  /*
  //Manual minimum label
  if (geom->textureMode == TEXTURE_MODE_MCHITS) sprintf(tmp, "%.3E", (double)geom->texMin_MC);
  else if (geom->textureMode == TEXTURE_MODE_FLUX) sprintf(tmp, "%.3E", geom->texMin_flux);
  else if (geom->textureMode == TEXTURE_MODE_POWER) sprintf(tmp, "%.3E", geom->texMin_power);
  texMinText->SetText(tmp);
  //Manual maximum label
  if (geom->textureMode == TEXTURE_MODE_MCHITS) sprintf(tmp, "%.3E", (double)geom->texMax_MC);
  else if (geom->textureMode == TEXTURE_MODE_FLUX) sprintf(tmp, "%.3E", geom->texMax_flux);
  else if (geom->textureMode == TEXTURE_MODE_POWER) sprintf(tmp, "%.3E", geom->texMax_power);
  texMaxText->SetText(tmp);
  */

  //Autoscale minimum label
  if (geom->textureMode==TEXTURE_MODE_MCHITS) sprintf(tmp,"%.3E",(double)geom->textureMin_auto.count);
  else if (geom->textureMode==TEXTURE_MODE_FLUX) sprintf(tmp,"%.3E",geom->textureMin_auto.flux);
  else if (geom->textureMode==TEXTURE_MODE_POWER) sprintf(tmp,"%.3E",geom->textureMin_auto.power);
  texCMinText->SetText(tmp);
  //Autoscale maximum label
   if (geom->textureMode==TEXTURE_MODE_MCHITS) sprintf(tmp,"%.3E",(double)geom->textureMax_auto.count);
  else if (geom->textureMode==TEXTURE_MODE_FLUX) sprintf(tmp,"%.3E",geom->textureMax_auto.flux);
  else if (geom->textureMode==TEXTURE_MODE_POWER) sprintf(tmp,"%.3E",geom->textureMax_auto.power);
  texCMaxText->SetText(tmp);

  texAutoScale->SetState(geom->texAutoScale);
  logBtn->SetState(geom->texLogScale);
  gradient->SetScale(geom->texLogScale?LOG_SCALE:LINEAR_SCALE);
  if( !geom->texAutoScale ) {
    if (geom->textureMode==TEXTURE_MODE_MCHITS) gradient->SetMinMax((double)geom->textureMin_manual.count,(double)geom->textureMax_manual.count);
	else if (geom->textureMode==TEXTURE_MODE_FLUX) gradient->SetMinMax(geom->textureMin_manual.flux,geom->textureMax_manual.flux);
	else if (geom->textureMode==TEXTURE_MODE_POWER) gradient->SetMinMax(geom->textureMin_manual.power,geom->textureMax_manual.power);
  } else {
	  if (geom->textureMode == TEXTURE_MODE_MCHITS) gradient->SetMinMax((double)geom->textureMin_auto.count, (double)geom->textureMax_auto.count);
	  else if (geom->textureMode == TEXTURE_MODE_FLUX) gradient->SetMinMax(geom->textureMin_auto.flux, geom->textureMax_auto.flux);
	  else if (geom->textureMode == TEXTURE_MODE_POWER) gradient->SetMinMax(geom->textureMin_auto.power, geom->textureMax_auto.power);
  }
  colormapBtn->SetState(geom->texColormap);
  gradient->SetType( geom->texColormap?GRADIENT_COLOR:GRADIENT_BW );
  modeCombo->SetSelectedIndex(geom->textureMode);
  UpdateSize();

}

void TextureSettings::Display(Worker *w,GeometryViewer **v) {

  worker = w;
  geom = w->GetSynradGeometry();
  viewers = v;
  if(!geom->IsLoaded()) {
	  GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
	  return;
  }
  SetVisible(true);
  Update();

  char tmp[64];
  //Manual minimum label
  if (geom->textureMode == TEXTURE_MODE_MCHITS) sprintf(tmp, "%.3E", (double)geom->textureMin_manual.count);
  else if (geom->textureMode == TEXTURE_MODE_FLUX) sprintf(tmp, "%.3E", geom->textureMin_manual.flux);
  else if (geom->textureMode == TEXTURE_MODE_POWER) sprintf(tmp, "%.3E", geom->textureMin_manual.power);
  texMinText->SetText(tmp);
  //Manual maximum label
  if (geom->textureMode == TEXTURE_MODE_MCHITS) sprintf(tmp, "%.3E", (double)geom->textureMax_manual.count);
  else if (geom->textureMode == TEXTURE_MODE_FLUX) sprintf(tmp, "%.3E", geom->textureMax_manual.flux);
  else if (geom->textureMode == TEXTURE_MODE_POWER) sprintf(tmp, "%.3E", geom->textureMax_manual.power);
  texMaxText->SetText(tmp);

}

void TextureSettings::ProcessMessage(GLComponent *src,int message) {

  switch(message) {
    case MSG_BUTTON:

    if (src==updateButton) {

      double min,max;

      if( !texMinText->GetNumber(&min) ) {
        GLMessageBox::Display("Invalid minimum value","Error",GLDLG_OK,GLDLG_ICONERROR);
        Update();
        return;
      }
      if( !texMaxText->GetNumber(&max) ) {
        GLMessageBox::Display("Invalid maximum value","Error",GLDLG_OK,GLDLG_ICONERROR);
        Update();
        return;
      }
      if( min>=max ) {
        GLMessageBox::Display("min must be lower than max","Error",GLDLG_OK,GLDLG_ICONERROR);
        Update();
        return;
      }
	  if (geom->textureMode==TEXTURE_MODE_MCHITS) {
		  geom->textureMin_manual.count = (llong)min;
		  geom->textureMax_manual.count = (llong)max;
	  }
	  else if (geom->textureMode == TEXTURE_MODE_FLUX) {
		  geom->textureMin_manual.flux = min;
		  geom->textureMax_manual.flux = max;
	  }
	  else if (geom->textureMode==TEXTURE_MODE_POWER) {
		  geom->textureMin_manual.power = min;
		  geom->textureMax_manual.power = max;
	  }
      geom->texAutoScale = texAutoScale->GetState();
      worker->Update(0.0f);
      Update();

    } else if (src==setCurrentButton) {
		if (geom->textureMode==TEXTURE_MODE_MCHITS) {
			geom->textureMin_manual.count = geom->textureMin_auto.count;
			geom->textureMax_manual.count = geom->textureMax_auto.count;
		}
		else if (geom->textureMode==TEXTURE_MODE_FLUX) {
			geom->textureMin_manual.flux = geom->textureMin_auto.flux;
			geom->textureMax_manual.flux = geom->textureMax_auto.flux;
		}
		else if (geom->textureMode==TEXTURE_MODE_POWER) {
			geom->textureMin_manual.power = geom->textureMin_auto.power;
			geom->textureMax_manual.power = geom->textureMax_auto.power;
		}
		texMinText->SetText(texCMinText->GetText());
		texMaxText->SetText(texCMaxText->GetText());
		texAutoScale->SetState(false);
		geom->texAutoScale=false;
		try {
				worker->Update(0.0f);
			} catch(Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
			}
		Update();
	}
    break;

    case MSG_TOGGLE:
    if (src==colormapBtn) {
      //for(int i=0;i<MAX_VIEWER;i++) viewers[i]->showColormap = colormapBtn->GetState();
      geom->texColormap = colormapBtn->GetState();
      worker->Update(0.0f);
      Update();
    } else if (src==texAutoScale) {
      geom->texAutoScale = texAutoScale->GetState();
      worker->Update(0.0f);
      Update();
    } else if (src==logBtn) {
      geom->texLogScale = logBtn->GetState();
      gradient->SetScale(geom->texLogScale?LOG_SCALE:LINEAR_SCALE);
      worker->Update(0.0f);
      Update();
    }
    break;

    case MSG_TEXT:
		ProcessMessage(updateButton,MSG_BUTTON);
		break;

	case MSG_COMBO:
		if(src==modeCombo) {
			geom->textureMode=modeCombo->GetSelectedIndex();
			worker->Update(0.0f);
			char tmp[256];
			if (geom->textureMode==TEXTURE_MODE_MCHITS) {
				sprintf(tmp,"%g",(double)geom->textureMin_manual.count);
				texMinText->SetText(tmp);
				sprintf(tmp,"%g",(double)geom->textureMax_manual.count);
				texMaxText->SetText(tmp);
			} else if (geom->textureMode==TEXTURE_MODE_FLUX) {
				sprintf(tmp,"%g",geom->textureMin_manual.flux);
				texMinText->SetText(tmp);
				sprintf(tmp,"%g",geom->textureMax_manual.flux);
				texMaxText->SetText(tmp);
			} else if (geom->textureMode==TEXTURE_MODE_POWER) {
				sprintf(tmp,"%g",geom->textureMin_manual.power);
				texMinText->SetText(tmp);
				sprintf(tmp,"%g",geom->textureMax_manual.power);
				texMaxText->SetText(tmp);
			}
			Update();
		}
		break;
  }

  GLWindow::ProcessMessage(src,message);
}
