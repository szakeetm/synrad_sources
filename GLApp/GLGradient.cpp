/*
  File:        GLGradient.cpp
  Description: Button class (SDL/OpenGL OpenGL application framework)
  Author:      J-L PONS (2007)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#include "GLWindow.h"
#include "GLGradient.h"
#include "GLToolkit.h"
#include <math.h>

// Rainbow
COLORREF rainbowCol[] = { 0x000000,   //Black
	                        0xFF0000,   //Red
	                        0xFF8030,   //Orange
                          0xF8F800,   //Yellow
                          0x28FF28,   //Green
										    	0x0040FF,   //Light Blue
										    	0x0000B0,   //Dark Blue
                          0x800090,   //Purple 1
                          0xF00080};  //Purple 2



// ---------------------------------------------------------------------

GLGradient::GLGradient(int compId):GLComponent(compId) {

  colorTex = 0;
  gType = GRADIENT_BW;
  mouseCursor = FALSE;
  axis = new GLAxis(this,HORIZONTAL_DOWN);
  axis->SetAnnotation(VALUE_ANNO);
  axis->SetAutoScale(FALSE);
  axis->SetMinimum(0.0);
  axis->SetMaximum(100.0);
  mouseValue = new GLLabel("---");
  mouseValue->SetParent(this);
  cursorPos=1.0;
  gWidth = 1;
  gHeight = 1;
  gPosX = 0;
  gPosY = 0;

}

// ----------------------------------------------------------

void GLGradient::SetScale(int scale) {
  axis->SetScale(scale);
  if(scale==LOG_SCALE)
    axis->SetLabelFormat(SCIENTIFIC_FORMAT);
  else
    axis->SetLabelFormat(AUTO_FORMAT);
  UpdateValue();
}

// ----------------------------------------------------------

int GLGradient::GetScale() {
  return axis->GetScale();
}

// ----------------------------------------------------------

GLGradient::~GLGradient() {
  SAFE_DELETE(axis);
}

// ----------------------------------------------------------

void GLGradient::SetMouseCursor(BOOL enable) {
  mouseCursor = enable;
}

// ----------------------------------------------------------

void GLGradient::SetMinMax(double min,double max) {

  if( min<max ) {
    if( axis->GetScale()==LOG_SCALE ) {
      if(min<1e-20) min=1e-20;
      if(max<1e-20) max=1e-20;
    }
    axis->SetMinimum(min);
    axis->SetMaximum(max);
  }

}

void GLGradient::SetType(int type) {
  gType = type;
}

// ----------------------------------------------------------

void GLGradient::InvalidateDeviceObjects() {
  DELETE_TEX(colorTex);
  DELETE_TEX(bwTex);
}

// ----------------------------------------------------------

void GLGradient::RestoreDeviceObjects() {

  // Rainbow texture
  DWORD *buff32 = (DWORD *)malloc( 512*4*sizeof(DWORD) );

	for(int i=0;i<512;i++) {
 
    double r1,g1,b1;
    double r2,g2,b2;
		int colId = i/64;

    r1 = (double) ((rainbowCol[colId] >> 16) & 0xFF);
    g1 = (double) ((rainbowCol[colId] >>  8) & 0xFF);
    b1 = (double) ((rainbowCol[colId] >>  0) & 0xFF);
   
    r2 = (double) ((rainbowCol[colId+1] >> 16) & 0xFF);
    g2 = (double) ((rainbowCol[colId+1] >>  8) & 0xFF);
    b2 = (double) ((rainbowCol[colId+1] >>  0) & 0xFF);
  
    double rr = (double)(i-colId*64) / 64.0;
    SATURATE(rr,0.0,1.0);
    COLORREF c = (COLORREF)((int)( r1 + (r2-r1)*rr )+
                            (int)( g1 + (g2-g1)*rr )* 256 +
                            (int)( b1 + (b2-b1)*rr )* 65536 );

    buff32[i+0*512] = c;
    buff32[i+1*512] = c;
    buff32[i+2*512] = c;
    buff32[i+3*512] = c;

	}

  // Saturation color
  buff32[511+0*512] = 0xFFFFFF;
  buff32[511+1*512] = 0xFFFFFF;
  buff32[511+2*512] = 0xFFFFFF;
  buff32[511+3*512] = 0xFFFFFF;
  buff32[510+0*512] = 0xFFFFFF;
  buff32[510+1*512] = 0xFFFFFF;
  buff32[510+2*512] = 0xFFFFFF;
  buff32[510+3*512] = 0xFFFFFF;

  glEnable(GL_TEXTURE_2D);

  glGenTextures(1,&colorTex);
  glBindTexture(GL_TEXTURE_2D,colorTex);
  glTexImage2D (
    GL_TEXTURE_2D,      // Type
    0,                  // No Mipmap
    4,                  // Format RGBA
    512,                // Width
    4,                  // Height
    0,                  // Border
    GL_RGBA,            // Format RGBA
    GL_UNSIGNED_BYTE,   // 8 Bit/color
    buff32              // Data
  );   

	for(int i=0;i<512;i++) {
   
    double rr = (double)(i) / 2.0;
    COLORREF c = (COLORREF)((int)( rr  )+
                            (int)( rr  )* 256 +
                            (int)( rr  )* 65536 );

    buff32[i+0*512] = c;
    buff32[i+1*512] = c;
    buff32[i+2*512] = c;
    buff32[i+3*512] = c;

	}

  glGenTextures(1,&bwTex);
  glBindTexture(GL_TEXTURE_2D,bwTex);
  glTexImage2D (
    GL_TEXTURE_2D,      // Type
    0,                  // No Mipmap
    4,                  // Format RGBA
    512,                // Width
    4,                  // Height
    0,                  // Border
    GL_RGBA,            // Format RGBA
    GL_UNSIGNED_BYTE,   // 8 Bit/color
    buff32              // Data
  );   

  free(buff32);

}

void GLGradient::UpdateValue() {

  double min = axis->GetMin();
  double max = axis->GetMax();
  double v;

  if( axis->GetScale()==LOG_SCALE ) {
    v = pow( 10.0 , cursorPos*(max-min)+min );
  } else {
    v = cursorPos*(max-min)+min;
  }
        
  mouseValue->SetText(axis->FormatValue(v,0.0));

}

void GLGradient::ManageEvent(SDL_Event *evt) {

  GLContainer::ManageEvent(evt);
  GLContainer::RelayEvent(evt);

  if( !evtProcessed ) {

    switch(evt->type) {
      case SDL_MOUSEBUTTONUP:
      break;
    case SDL_MOUSEBUTTONDOWN:
      break;
    case SDL_MOUSEBUTTONDBLCLICK:
      break;
    case SDL_MOUSEMOTION:
      if( mouseCursor ) {
        int ex = GetWindow()->GetX(this,evt) + posX;
        int ey = GetWindow()->GetY(this,evt) + posY;
        cursorPos = (double)(ex-gPosX) / (double)(gWidth);
        SATURATE(cursorPos,0.0,1.0);
        UpdateValue();
      } else {
        cursorPos = -1.0;
      }

      break;
    }

  }

}

void GLGradient::Paint() {

  if(!parent) return;
  GLComponent::Paint();

  if (height <= 20 || width <= 20)
    return;

  int offX = 0;

  if( mouseCursor ) {
    mouseValue->SetBounds(posX+2,posY+8,60,19);
    mouseValue->Paint();
    offX = 35;
  }

  // Gradient texture
  glDisable(GL_CULL_FACE);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_2D);
  glDisable(GL_BLEND);
  switch(gType) {
    case GRADIENT_BW:
    glBindTexture(GL_TEXTURE_2D,bwTex);
    break;
    case GRADIENT_COLOR:
    glBindTexture(GL_TEXTURE_2D,colorTex);
    break;
  }

  gWidth  = width-(30+offX);
  gHeight = 20;
  gPosX   = posX+offX+22;
  gPosY   = posY+5;

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glColor3f(1.0f,1.0f,1.0f);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0f,0.0f);glVertex2i(gPosX       ,gPosY);
  glTexCoord2f(1.0f,0.0f);glVertex2i(gPosX+gWidth,gPosY);
  glTexCoord2f(1.0f,1.0f);glVertex2i(gPosX+gWidth,gPosY+gHeight);
  glTexCoord2f(0.0f,1.0f);glVertex2i(gPosX       ,gPosY+gHeight);
  glEnd();

  GLCColor black;
  black.r = 0;black.g = 0;black.b = 0;
  axis->MeasureAxis(gWidth,0);
  axis->PaintAxisDirect(gPosX,gPosY+gHeight, black, 0, 0);

  if( mouseCursor ) {

    if( cursorPos>0.0 ) {
      glColor3f(0.0f,0.0f,0.0f);
      glBegin(GL_LINES);
      int xpos = (int)( cursorPos*(double)gWidth + gPosX );
      glVertex2i(xpos,gPosY);
      glVertex2i(xpos,gPosY+gHeight);

      glVertex2i(xpos-3,gPosY-2);
      glVertex2i(xpos-1,gPosY);
      glVertex2i(xpos,gPosY);
      glVertex2i(xpos+2,gPosY-2);

      glVertex2i(xpos-3,gPosY+gHeight+3);
      glVertex2i(xpos-1,gPosY+gHeight+1);
      glVertex2i(xpos,gPosY+gHeight+1);
      glVertex2i(xpos+2,gPosY+gHeight+3);

      glEnd();
    }

  }

}