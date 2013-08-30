/*
  File:        GLToolkit.cpp
  Description: SDL/OpenGL OpenGL application framework
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

#include "GLToolkit.h"
#include "GLApp.h"
#include "GLWindow.h"
#include "GLMatrix.h"
#include "GLSprite.h"
#include <math.h>
#include <malloc.h>
#include <CImage.h>
#include "../SynRad.h"

static GLFont2D *dlgFont  = NULL;
static GLFont2D *dlgFontB = NULL;
static int scrWidth;
static int scrHeight;

// dashed line
static GLushort dashDotPattern  = 0xE4E4;
static GLushort dotPattern      = 0x8888;
static GLushort dashPattern     = 0xF0F0;
static GLushort longDashPattern = 0xFF00;

// Cursors
static int currentCursor = -1;
static SDL_Cursor *defCursor = NULL;
static SDL_Cursor *busyCursor = NULL;
static SDL_Cursor *sizeCursor = NULL;
static SDL_Cursor *sizeHCursor = NULL;
static SDL_Cursor *sizeHSCursor = NULL;
static SDL_Cursor *sizeVCursor = NULL;
static SDL_Cursor *zoomCursor = NULL;
static SDL_Cursor *zoom2Cursor = NULL;
static SDL_Cursor *textCursor = NULL;
static SDL_Cursor *selAddCursor = NULL;
static SDL_Cursor *selDelCursor = NULL;
static SDL_Cursor *handCursor = NULL;
static SDL_Cursor *rotateCursor = NULL;
static SDL_Cursor *vertexCursor = NULL;
static SDL_Cursor *trajCursor = NULL;
static SDL_Cursor *vertexAddCursor = NULL;
static SDL_Cursor *vertexClrCursor = NULL;

// Texture for GUI components
static Sprite2D *compTex = NULL;

// Error log
#define MAX_LOG 256

static int nbLog=0;
static int logLength=0;
static char *logs[MAX_LOG];

extern GLApplication *theApp;

// ---------------------------------------------------------------------

GLFont2D *GLToolkit::GetDialogFont() {
  return dlgFont;
}

// ---------------------------------------------------------------------

GLFont2D *GLToolkit::GetDialogFontBold() {
  return dlgFontB;
}

// ---------------------------------------------------------------------

void GLToolkit::InvalidateDeviceObjects() {

  if( dlgFont ) {
    dlgFont->InvalidateDeviceObjects();
    SAFE_DELETE(dlgFont);
  }
  if( dlgFontB ) {
    dlgFontB->InvalidateDeviceObjects();
    SAFE_DELETE(dlgFontB);
  }
  if( compTex ) {
    compTex->InvalidateDeviceObjects();
    SAFE_DELETE(compTex);
  }

}

// ---------------------------------------------------------------

int GLToolkit::GetCursor() {
  return currentCursor;
}

// ---------------------------------------------------------------

void GLToolkit::SetCursor(int cursor) {

  if( cursor==currentCursor ) 
    return;

  currentCursor = cursor;
  switch(currentCursor) {
    case CURSOR_DEFAULT:
      SDL_SetCursor(defCursor);
      break;
	case CURSOR_BUSY:
	  SDL_SetCursor(busyCursor);
	  break;
    case CURSOR_SIZE:
      SDL_SetCursor(sizeCursor);
      break;
    case CURSOR_SIZEH:
      SDL_SetCursor(sizeHCursor);
      break;
    case CURSOR_SIZEHS:
      SDL_SetCursor(sizeHSCursor);
      break;
    case CURSOR_SIZEV:
      SDL_SetCursor(sizeVCursor);
      break;
    case CURSOR_ZOOM:
      SDL_SetCursor(zoomCursor);
      break;
    case CURSOR_ZOOM2:
      SDL_SetCursor(zoom2Cursor);
      break;
    case CURSOR_TEXT:
      SDL_SetCursor(textCursor);
      break;
    case CURSOR_SELADD:
      SDL_SetCursor(selAddCursor);
      break;
    case CURSOR_SELDEL:
      SDL_SetCursor(selDelCursor);
      break;
    case CURSOR_HAND:
      SDL_SetCursor(handCursor);
      break;
    case CURSOR_ROTATE:
      SDL_SetCursor(rotateCursor);
      break;
	case CURSOR_VERTEX:
      SDL_SetCursor(vertexCursor);
      break;
	case CURSOR_VERTEX_ADD:
      SDL_SetCursor(vertexAddCursor);
      break;
	case CURSOR_VERTEX_CLR:
      SDL_SetCursor(vertexClrCursor);
      break;
	case CURSOR_TRAJ:
      SDL_SetCursor(trajCursor);
      break;
  }

}

// ---------------------------------------------------------------------

SDL_Cursor *InitCursor(char *pngName,int tx,int ty) {

  int i, row, col;
  Uint8 data[4*32];
  Uint8 mask[4*32];
  memset(data,0,4*32);
  memset(mask,0,4*32);

  CImage img;
  if( img.LoadImage(pngName) ) {

    i = -1;
    BYTE *image = img.GetData();
    for ( row=0; row<32; ++row ) {
      for ( col=0; col<32; ++col ) {
        if ( col % 8 ) {
          data[i] <<= 1;
          mask[i] <<= 1;
        } else {
          ++i;
        }
        switch (image[(col+row*32)*3]) {
          case 0:
            data[i] |= 0x01;
            mask[i] |= 0x01;
            break;
          case 0xFF:
            mask[i] |= 0x01;
            break;
        }
      }
    }
    img.Release();
    return SDL_CreateCursor(data,mask,32,32,tx,ty);

  } else {
    return NULL;
  }


}

// ---------------------------------------------------------------------

void GLToolkit::SetIcon32x32(char *pngName) {

  CImage img;
  if( img.LoadImage(pngName) ) {
    SDL_Surface *s = SDL_CreateRGBSurface(SDL_SWSURFACE,32,32,24,0,0,0,0);
    SDL_LockSurface(s);
    memcpy(s->pixels , img.GetData() , 3*32*32);
    SDL_UnlockSurface(s);
    SDL_WM_SetIcon(s, NULL);
    img.Release();
  }

}

// ---------------------------------------------------------------------

BOOL GLToolkit::RestoreDeviceObjects(int width,int height) {

  scrWidth = width;
  scrHeight = height;

  dlgFont = new GLFont2D("images/fnt_sansserif_8.png");
  dlgFont->SetTextSize(12,14);
  dlgFont->SetVariableWidth(TRUE);
  if( !dlgFont->RestoreDeviceObjects(scrWidth,scrHeight) )
    return FALSE;

  dlgFontB = new GLFont2D("images/fnt_sansserif_8b.png");
  dlgFontB->SetTextSize(12,14);
  dlgFontB->SetVariableWidth(TRUE);
  if( !dlgFontB->RestoreDeviceObjects(scrWidth,scrHeight) )
    return FALSE;

  if(!defCursor) defCursor = InitCursor("images/cursor_default.png",0,0);
  if(!busyCursor) busyCursor = InitCursor("images/cursor_busy.png",0,0);
  if(!sizeCursor) sizeCursor = InitCursor("images/cursor_resize.png",5,5);
  if(!sizeHCursor) sizeHCursor = InitCursor("images/cursor_resizeh.png",10,4);
  if(!sizeHSCursor) sizeHSCursor = InitCursor("images/cursor_resizehs.png",6,5);
  if(!sizeVCursor) sizeVCursor = InitCursor("images/cursor_resizev.png",4,10);
  if(!zoomCursor) zoomCursor = InitCursor("images/cursor_zoom.png",9,9);
  if(!zoom2Cursor) zoom2Cursor = InitCursor("images/cursor_zoom2.png",9,9);
  if(!textCursor) textCursor = InitCursor("images/cursor_text.png",3,7);
  if(!selAddCursor) selAddCursor = InitCursor("images/cursor_selp.png",0,0);
  if(!selDelCursor) selDelCursor = InitCursor("images/cursor_selm.png",0,0);
  if(!handCursor) handCursor = InitCursor("images/cursor_hand.png",9,9);
  if(!rotateCursor) rotateCursor = InitCursor("images/cursor_rotate.png",9,9);
  if(!vertexCursor) vertexCursor = InitCursor("images/cursor_vertex.png",0,0);
  if(!trajCursor) trajCursor = InitCursor("images/cursor_traj.png",0,0);
  if(!vertexAddCursor) vertexAddCursor = InitCursor("images/cursor_vertexp.png",0,0);
  if(!vertexClrCursor) vertexClrCursor = InitCursor("images/cursor_vertexm.png",0,0);

  compTex = new Sprite2D();
  if( !compTex->RestoreDeviceObjects("images/gui_background.png","images/gui_backgrounda.png",256,256) )
    return FALSE;

  SetCursor(CURSOR_DEFAULT);
  return TRUE;

}

// ---------------------------------------------------------------------

void GLToolkit::Log(char *message) {

  if( nbLog < MAX_LOG ) {
    logs[nbLog] = _strdup(message);
    logLength += (int)strlen(logs[nbLog])+2;
    nbLog++;
  }

}

// ---------------------------------------------------------------------

char *GLToolkit::GetLogs() {

  if( logLength ) {
    char *ret = (char *)malloc(logLength+512);
    strcpy(ret,"");
    for(int i=0;i<nbLog;i++) {
      strcat(ret,logs[i]);
      strcat(ret,"\n");
    }
    return ret;
  }

  return NULL;

}

void GLToolkit::ClearLogs() {

  for(int i=0;i<nbLog;i++)
    SAFE_FREE(logs[i]);
  logLength = 0;
  nbLog=0;

}

// ---------------------------------------------------------------------

void GLToolkit::GetScreenSize(int *width,int *height) {
  *width = scrWidth;
  *height = scrHeight;
}

// ---------------------------------------------------------------

void GLToolkit::SetViewport(GLVIEWPORT v) {
  GLToolkit::SetViewport(v.x,v.y,v.width,v.height);
}

// ---------------------------------------------------------------

void GLToolkit::SetViewport(int x,int y,int width,int height) {

  int vy = scrHeight - (y+height);
  int vx = x;
  glViewport(vx,vy,width,height);

}

// ---------------------------------------------------------------

#define TW 256.0f

void GLToolkit::DrawButtonBack(int x,int y,int width,int height,int state) {

  int w1 = width / 2;
  int w2 = width - w1;

  float fw1 = (float)w1/TW;
  float fw2 = (float)w2/TW;
  float ft  = 21.0f / TW;

  if( !state ) compTex->SetColor(1.0f,1.0f,1.0f);
  else         compTex->SetColor(0.8f,0.8f,0.8f);

  // Left part
  compTex->UpdateSprite(x,y,x+w1,y+height);
  compTex->SetSpriteMapping(0.0f,0.0f,fw1,ft);
  compTex->Render(FALSE);

  // Right part
  compTex->UpdateSprite(x+w1,y,x+width,y+height);
  compTex->SetSpriteMapping(1.0f-fw2,0.0f,1.0f,ft);
  compTex->Render(FALSE);

}

// ---------------------------------------------------------------

void GLToolkit::DrawTinyButton(int x,int y,int state) {

  if( !state ) compTex->SetColor(1.0f,1.0f,1.0f);
  else         compTex->SetColor(0.8f,0.8f,0.8f);

  compTex->UpdateSprite(x,y,x+12,y+12);
  compTex->SetSpriteMapping(0.0f,95.0f/TW,12.0f/TW,107.0f/TW);
  compTex->Render(FALSE);

}

// ---------------------------------------------------------------


void GLToolkit::DrawSmallButton(int x,int y,int state) {

  if( !state ) compTex->SetColor(1.0f,1.0f,1.0f);
  else         compTex->SetColor(0.8f,0.8f,0.8f);

  compTex->UpdateSprite(x,y,x+15,y+17);
  compTex->SetSpriteMapping(0.0f,24.0f/TW,15.0f/TW,41.0f/TW);
  compTex->Render(FALSE);

}

// ---------------------------------------------------------------

void GLToolkit::DrawToggle(int x,int y) {

  compTex->SetColor(1.0f,1.0f,1.0f);
  compTex->UpdateSprite(x,y,x+12,y+12);
  compTex->SetSpriteMapping(17.0f/TW,24.0f/TW,29.0f/TW,36.0f/TW);
  compTex->Render(FALSE);

}

// ---------------------------------------------------------------

void GLToolkit::DrawBar(int x,int y,int width,int height) {

  compTex->SetColor(1.0f,1.0f,1.0f);
  compTex->UpdateSprite(x,y,x+width,y+height);
  compTex->SetSpriteMapping(0.0f,44.0f/TW,1.0f,51.0f/TW);
  compTex->Render(FALSE);

}

// ---------------------------------------------------------------

void GLToolkit::DrawTextBack(int x,int y,int width,int height,int rBack, int gBack, int bBack) {

  int w1 = width / 2;
  int w2 = width - w1;

  float fw1 = (float)w1/TW;
  float fw2 = (float)w2/TW;
  float ft  = 21.0f / TW;

  compTex->SetColor(rBack/240.0f,gBack/240.0f,bBack/240.0f);

  // Left part
  compTex->UpdateSprite(x,y,x+w1,y+height);
  compTex->SetSpriteMapping(0.0f,54.0f/TW,fw1,75.0f/TW);
  compTex->Render(FALSE);

  // Right part
  compTex->UpdateSprite(x+w1,y,x+width,y+height);
  compTex->SetSpriteMapping(1.0f-fw2,54.0f/TW,1.0f,75.0f/TW);
  compTex->Render(FALSE);

}

// ---------------------------------------------------------------

void GLToolkit::DrawVGradientBox(int x,int y,int width,int height,BOOL shadow,BOOL iShadow,BOOL isEtched) {

  compTex->SetColor(1.0f,1.0f,1.0f);
  compTex->UpdateSprite(x,y,x+width,y+height);
  compTex->SetSpriteMapping(0.1f,3.0f/TW,0.5f,17.0f/TW);
  compTex->Render(TRUE);

  glDisable(GL_TEXTURE_2D);
  glDisable(GL_BLEND);
  glLineWidth(1.0f);
  DrawBorder(x,y,width,height,shadow,iShadow,isEtched);

}

// ---------------------------------------------------------------

void GLToolkit::DrawHGradientBox(int x,int y,int width,int height,BOOL shadow,BOOL iShadow,BOOL isEtched) {

  compTex->SetColor(1.0f,1.0f,1.0f);
  compTex->UpdateSprite(x,y,x+width,y+height);
  compTex->SetSpriteMapping(32.0f/TW,25.0f/TW,45.0f/TW,35.0f/TW);
  compTex->Render(TRUE);

  glDisable(GL_TEXTURE_2D);
  glDisable(GL_BLEND);
  glLineWidth(1.0f);
  DrawBorder(x,y,width,height,shadow,iShadow,isEtched);

}

void GLToolkit::DrawHScroll(int x,int y,int width,int height,int state) {

  int w1 = width / 2;
  int w2 = width - w1;
  if( w1>200 ) {
    w1 = 200;
    w2 = 200;
  }
  int left = width-(w1+w2);
  int nbStep = left/200 + ((left%200)?1:0);

  float fw1 = (float)w1/TW;
  float fw2 = (float)w2/TW;
  float ftu  = 109.0f / TW;
  float ftd  = 121.0f / TW;

  switch(state) {

    case 0:
     compTex->SetColor(1.0f,1.0f,1.0f);
     break;
    case 1:
     compTex->SetColor(0.9f,0.9f,0.9f);
     break;
    case 2:
     ftu  = 123.0f / TW;
     ftd  = 135.0f / TW;
     compTex->SetColor(1.0f,1.0f,1.0f);
     break;

  }

  // First part
  compTex->UpdateSprite(x,y,x+w1,y+height);
  compTex->SetSpriteMapping(0.0f,ftu,fw1,ftd);
  compTex->Render(FALSE);

  // Middle parts
  compTex->SetSpriteMapping(0.1f,ftu,0.9f,ftd);
  for(int i=0;i<left;i+=200) {
    int sleft = MIN(left-i,200);
    compTex->UpdateSprite(x+w1+i,y,x+w1+i+sleft,y+height);
    compTex->Render(FALSE);
  }

  // Last part
  compTex->UpdateSprite(x+w1+left,y,x+width,y+height);
  compTex->SetSpriteMapping(1.0f-fw2,ftu,1.0f,ftd);
  compTex->Render(FALSE);

}

void GLToolkit::DrawVScroll(int x,int y,int width,int height,int state) {

  int h1 = height / 2;
  int h2 = height - h1;
  if( h1>200 ) {
    h1 = 200;
    h2 = 200;
  }
  int left = height-(h1+h2);
  int nbStep = left/200 + ((left%200)?1:0);

  float fw1 = (float)h1/TW;
  float fw2 = (float)h2/TW;
  float ftu  = 109.0f / TW;
  float ftd  = 121.0f / TW;

  switch(state) {

    case 0:
     compTex->SetColor(1.0f,1.0f,1.0f);
     break;
    case 1:
     compTex->SetColor(0.9f,0.9f,0.9f);
     break;
    case 2:
     ftu  = 123.0f / TW;
     ftd  = 135.0f / TW;
     compTex->SetColor(1.0f,1.0f,1.0f);
     break;

  }

  // First part
  compTex->UpdateSprite(x,y,x+width,y+h1);
  compTex->SetSpriteMapping(0.0f,ftu,fw1,ftd);
  compTex->Render90(FALSE);

  // Middle parts
  compTex->SetSpriteMapping(0.1f,ftu,0.9f,ftd);
  for(int i=0;i<left;i+=200) {
    int sleft = MIN(left-i,200);
    compTex->UpdateSprite(x,y+h1+i,x+width,y+h1+i+sleft);
    compTex->Render90(FALSE);
  }

  // Last part
  compTex->UpdateSprite(x,y+h1+left,x+width,y+height);
  compTex->SetSpriteMapping(1.0f-fw2,ftu,1.0f,ftd);
  compTex->Render90(FALSE);

}

// ---------------------------------------------------------------

void GLToolkit::DrawHIGradientBox(int x,int y,int width,int height,BOOL shadow,BOOL iShadow,BOOL isEtched) {

  compTex->SetColor(1.0f,1.0f,1.0f);
  compTex->UpdateSprite(x,y,x+width,y+height);
  compTex->SetSpriteMapping(49.0f/TW,25.0f/TW,62.0f/TW,35.0f/TW);
  compTex->Render(TRUE);

  glDisable(GL_TEXTURE_2D);
  glDisable(GL_BLEND);
  glLineWidth(1.0f);
  DrawBorder(x,y,width,height,shadow,iShadow,isEtched);

}

// ---------------------------------------------------------------

void GLToolkit::Draw16x16(int x,int y,int xt,int yt) {

  compTex->SetColor(1.0f,1.0f,1.0f);
  compTex->UpdateSprite(x,y,x+16,y+16);
  compTex->SetSpriteMapping((float)xt/TW,(float)yt/TW,(float)(xt+16)/TW,(float)(yt+16)/TW);
  compTex->Render(FALSE);

}

// ---------------------------------------------------------------

void GLToolkit::DrawLumBitmap(int x,int y,int width,int height,BYTE *buffer) {

  glDisable(GL_TEXTURE_2D);
  glEnable(GL_BLEND);
  glRasterPos2i(x-1,y+1);
  glDrawPixels(width,height,GL_LUMINANCE_ALPHA,GL_UNSIGNED_BYTE,buffer);

}

// ---------------------------------------------------------------------

void GLToolkit::DrawBorder(int x,int y,int width,int height,BOOL shadow,BOOL iShadow,BOOL isEtched) {

  float rL = 1.0f;
  float gL = 1.0f;
  float bL = 1.0f;
  float rD = 0.5f;
  float gD = 0.5f;
  float bD = 0.5f;

  if( shadow ) {
  
    // Shadow colors
    /* 
    float rL = (float)r / 100.0f;
    float gL = (float)g / 100.0f;
    float bL = (float)b / 100.0f;
    SATURATE(rL,0.0f,1.0f);
    SATURATE(gL,0.0f,1.0f);
    SATURATE(bL,0.0f,1.0f);

    float rD = (float)r / 500.0f;
    float gD = (float)g / 500.0f;
    float bD = (float)b / 500.0f;
    SATURATE(rD,0.0f,1.0f);
    SATURATE(gD,0.0f,1.0f);
    SATURATE(bD,0.0f,1.0f);
    */

    if(!iShadow) glColor3f(rL,gL,bL);
    else         glColor3f(rD,gD,bD);
    glBegin(GL_LINES);
    _glVertex2i(x,y);
    _glVertex2i(x+width,y);
    _glVertex2i(x,y);
    _glVertex2i(x,y+height);
    glEnd();

    if(iShadow) glColor3f(rL,gL,bL);
    else        glColor3f(rD,gD,bD);
    glBegin(GL_LINES);
    _glVertex2i(x+width,y+height);
    _glVertex2i(x-1    ,y+height);
    _glVertex2i(x+width,y+height);
    _glVertex2i(x+width,y);
    glEnd();

    #ifdef _DEBUG
      theApp->nbLine+=4;
    #endif

  } else if( isEtched ) {

    // Etched colors
    /*
    float rL = (float)r / 100.0f;
    float gL = (float)g / 100.0f;
    float bL = (float)b / 100.0f;
    SATURATE(rL,0.0f,1.0f);
    SATURATE(gL,0.0f,1.0f);
    SATURATE(bL,0.0f,1.0f);

    float rD = (float)r / 500.0f;
    float gD = (float)g / 500.0f;
    float bD = (float)b / 500.0f;
    SATURATE(rD,0.0f,1.0f);
    SATURATE(gD,0.0f,1.0f);
    SATURATE(bD,0.0f,1.0f);
    */

    glColor3f(rD,gD,bD);

    glBegin(GL_LINES);
    _glVertex2i(x,y);
    _glVertex2i(x+width,y);
    _glVertex2i(x,y);
    _glVertex2i(x,y+height);
    _glVertex2i(x+width-1,y+height-1);
    _glVertex2i(x+width-1,y);
    _glVertex2i(x+width-1,y+height-1);
    _glVertex2i(x,y+height-1);
    glEnd();

    glColor3f(rL,gL,bL);

    glBegin(GL_LINES);
    _glVertex2i(x+1,y+1);
    _glVertex2i(x+width,y+1);
    _glVertex2i(x+1,y+1);
    _glVertex2i(x+1,y+height);
    _glVertex2i(x+width,y+height);
    _glVertex2i(x+width,y);
    _glVertex2i(x+width,y+height);
    _glVertex2i(x,y+height);
    glEnd();

    #ifdef _DEBUG
      theApp->nbLine+=4;
    #endif

  }


}

// ---------------------------------------------------------------------

void GLToolkit::DrawBox(int x,int y,int width,int height,int r,int g,int b,BOOL shadow,BOOL iShadow,BOOL isEtched) {

  float rN = (float)r / 255.0f;
  float gN = (float)g / 255.0f;
  float bN = (float)b / 255.0f;
  SATURATE(rN,0.0f,1.0f);
  SATURATE(gN,0.0f,1.0f);
  SATURATE(bN,0.0f,1.0f);

  glDisable(GL_TEXTURE_2D);
  glDisable(GL_BLEND);
  glLineWidth(1.0f);
  glColor3f(rN,gN,bN);

  glBegin(GL_QUADS);
  
  glVertex2i(x,y);
  glVertex2i(x+width,y);
  
  glVertex2i(x+width,y+height);
  glVertex2i(x,y+height);
  glEnd();

#ifdef _DEBUG
  theApp->nbPoly++;
#endif

  DrawBorder(x,y,width,height,shadow,iShadow,isEtched);

}

// -------------------------------------------

void GLToolkit::DrawPoly(int lineWidth,int dashStyle,int r,int g,int b,
                         int nbPoint,int *pointX,int *pointY) {

  // Draw unclosed polygon
  if( lineWidth==0 ) return;

  glDisable(GL_TEXTURE_2D);
  glDisable(GL_CULL_FACE);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
  glLineWidth((float)lineWidth);

  float rN = (float)r / 255.0f;
  float gN = (float)g / 255.0f;
  float bN = (float)b / 255.0f;
  SATURATE(rN,0.0f,1.0f);
  SATURATE(gN,0.0f,1.0f);
  SATURATE(bN,0.0f,1.0f);

  if( dashStyle==DASHSTYLE_NONE ) {
    glDisable(GL_LINE_STIPPLE);
  } else {
    glEnable(GL_LINE_STIPPLE);
    switch(dashStyle) {
      case DASHSTYLE_DOT:
        glLineStipple(1,dotPattern);
        break;
      case DASHSTYLE_DASH:
        glLineStipple(1,dashPattern);
        break;
      case DASHSTYLE_LONG_DASH:
        glLineStipple(1,longDashPattern);
        break;
      case DASHSTYLE_DASH_DOT:
        glLineStipple(1,dashDotPattern);
        break;
    }
  }
  glColor3f(rN,gN,bN);

  glBegin(GL_LINE_STRIP);
  for(int i=0;i<nbPoint;i++)
    _glVertex2i(pointX[i],pointY[i]);
  glEnd();

  glDisable(GL_LINE_STIPPLE);

}

// -------------------------------------------

BOOL GLToolkit::IsInsidePoly(int x,int y,int *PointX,int *PointY,int nbPts) {

   // Fast method to check if a point is inside a polygon or not.
   // Works with convex and concav polys, orientation independant
   int n_updown,n_found,j;
   double x1,x2,y1,y2,a,b;
   double xp = (double)x;
   double yp = (double)y;
   BOOL inside = FALSE;

   n_updown=0;
   n_found=0;

   for (j = 0; j < nbPts; j++) {

     x1 = (double) PointX[j];
     y1 = (double) PointY[j];

     if (j < nbPts - 1) {
       x2 = (double) PointX[j + 1];
       y2 = (double) PointY[j + 1];
     } else {
       x2 = (double) PointX[0];
       y2 = (double) PointY[0];
     }

     if (xp > MIN(x1, x2) && xp <= MAX(x1, x2)) {
       a = (y2 - y1) / (x2 - x1);
       b = y1 - a * x1;
       if ((a * xp + b) < yp) {
         n_updown = n_updown + 1;
       } else {
         n_updown = n_updown - 1;
       }
       n_found++;
     }
   }

   if (n_updown<0) n_updown=-n_updown;
   inside = (((n_found/2)&1) ^ ((n_updown/2)&1));
   return inside;

}


// -------------------------------------------

BOOL GLToolkit::Get2DScreenCoord(float x,float y,float z,int *xe,int *ye) {

  GLfloat mProj[16];
  GLfloat mView[16];
  GLVIEWPORT g;

  // Compute location on screen
  glGetFloatv(GL_PROJECTION_MATRIX , mProj);
  glGetFloatv(GL_MODELVIEW_MATRIX , mView);
  glGetIntegerv(GL_VIEWPORT,(GLint *)&g);

  GLMatrix proj; proj.LoadGL(mProj);
  GLMatrix view; view.LoadGL(mView);
  GLMatrix m; m.Multiply(&proj,&view);

  float rx,ry,rz,rw;
  m.TransfomVec(x,y,z,1.0f,&rx,&ry,&rz,&rw);
  if(rw<=0.0f) return FALSE;
  *xe = (int)(((rx / rw) + 1.0f)  * (float)g.width/2.0f);
  *ye = (int)(((-ry / rw) + 1.0f) * (float)g.height/2.0f);

  return TRUE;

}

// -------------------------------------------

float GLToolkit::GetVisibility(double x,double y,double z,double nx,double ny,double nz) {

  GLfloat mView[16];
  float rx,ry,rz,rw;
  float ntx,nty,ntz,ntw;

  glGetFloatv( GL_MODELVIEW_MATRIX , mView );
  GLMatrix view; view.LoadGL(mView);
  view.TransfomVec((float)x,(float)y,(float)z,1.0f,&rx,&ry,&rz,&rw);
  view.TransfomVec((float)nx,(float)ny,(float)nz,0.0f,&ntx,&nty,&ntz,&ntw);
  return rx*ntx + ry*nty + rz*ntz;

}

// -------------------------------------------

float GLToolkit::GetCamDistance(GLfloat *mView,double x,double y,double z) {

  float rx,ry,rz,rw;
  GLMatrix view; view.LoadGL(mView);
  view.TransfomVec((float)x,(float)y,(float)z,1.0f,&rx,&ry,&rz,&rw);
  return sqrtf(rx*rx + ry*ry + rz*rz);

}

// -------------------------------------------
static GLMatrix   dsm;
static GLfloat    dsmProj[16];
static GLfloat    dsmView[16];
static GLVIEWPORT dsg;

void GLToolkit::DrawStringInit() {

  glGetFloatv( GL_PROJECTION_MATRIX , dsmProj );
  glGetFloatv( GL_MODELVIEW_MATRIX , dsmView );
  glGetIntegerv(GL_VIEWPORT,(GLint *)&dsg);

  GLMatrix proj; proj.LoadGL(dsmProj);
  GLMatrix view; view.LoadGL(dsmView);
  dsm.Multiply(&proj,&view);

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glOrtho(dsg.x,dsg.x+dsg.width,dsg.y+dsg.height,dsg.y,-1.0,1.0);
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

  glEnable(GL_TEXTURE_2D);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

}

void GLToolkit::DrawString(float x,float y,float z,char *str,GLFont2D *fnt,int offx,int offy) {

  // Compute location on screen
  float rx,ry,rz,rw;
  dsm.TransfomVec(x,y,z,1.0f,&rx,&ry,&rz,&rw);
  if(rw<0.0f) return;
  int xe = dsg.x +(int)(((rx / rw) + 1.0f) * (float)dsg.width/2.0f);
  int ye = dsg.y +(int)(((-ry / rw) + 1.0f) * (float)dsg.height/2.0f);

  fnt->DrawTextFast(xe+offx,ye+offy,str);

}

void GLToolkit::DrawStringRestore() {

  // restore transform matrix
  glMatrixMode( GL_PROJECTION );
  glLoadMatrixf(dsmProj);
  glMatrixMode( GL_MODELVIEW );
  glLoadMatrixf(dsmView);

}

// -------------------------------------------

void GLToolkit::DrawRule(double length,BOOL invertX,BOOL invertY,BOOL invertZ,double n) {

  DrawVector(0.0,0.0,0.0,(invertX)?-length:length,0.0,0.0,n);
  DrawVector(0.0,0.0,0.0,0.0,(invertY)?-length:length,0.0,n);
  DrawVector(0.0,0.0,0.0,0.0,0.0,(invertZ)?-length:length,n);
  glPointSize(4.0f);
  glBegin(GL_POINTS);
  glVertex3d(0.0,0.0,0.0);
  glEnd();

}

// -------------------------------------------

void GLToolkit::DrawVector(double x1,double y1,double z1,double x2,double y2,double z2,double nr) {

  double n = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) )*nr;
  double nx;
  double ny;
  double nz;

  // Choose a normal vector
  if( abs(x2-x1) > 1e-3 ) {
    // Oy
    nx = (z2-z1)/n;
    ny = (y2-y1)/n;
    nz = (x1-x2)/n;
  } else if( abs(y2-y1) > 1e-3 ) {
    // Oz
    nx = (y2-y1)/n;
    ny = (x1-x2)/n;
    nz = (z2-z1)/n;
  } else {
    // Ox
    nx = (x2-x1)/n;
    ny = (z2-z1)/n;
    nz = (y1-y2)/n;
  }

  glDisable(GL_TEXTURE_2D);
  glDisable(GL_LIGHTING);
  glDisable(GL_BLEND);

  glDisable(GL_CULL_FACE);

 glBegin(GL_LINES);
  glVertex3d(x1,y1,z1);
  glVertex3d(x2,y2,z2);

  if( n>3.0 ) {
    glVertex3d(x2,y2,z2);
    glVertex3d(x2+nx-(x2-x1)/n, y2+ny-(y2-y1)/n, z2+nz-(z2-z1)/n);

    glVertex3d(x2,y2,z2);
    glVertex3d(x2-nx-(x2-x1)/n, y2-ny-(y2-y1)/n ,z2-nz-(z2-z1)/n);

    glVertex3d(x2+nx-(x2-x1)/n, y2+ny-(y2-y1)/n, z2+nz-(z2-z1)/n);
    glVertex3d(x2-nx-(x2-x1)/n, y2-ny-(y2-y1)/n ,z2-nz-(z2-z1)/n);
  }

  glEnd();

}

// -------------------------------------------

void GLToolkit::PerspectiveLH(double fovy,double aspect,double zNear,double zFar) {

  // Create and load a left handed proj matrix
  double fov = (fovy/360.0) * (2.0*PI);
  double f = cos(fov/2.0) / sin(fov/2.0);

  GLfloat mProj[16];

  mProj[0] = (float)(f/aspect); mProj[1] = 0.0f; mProj[2] = 0.0f; mProj[3] = 0.0f;
  mProj[4] = 0.0f; mProj[5] = (float)(f); mProj[6] = 0.0f; mProj[7] = 0.0f;
  mProj[8] = 0.0f; mProj[9] = 0.0f; mProj[10] = (float)((zFar+zNear)/(zFar-zNear)); mProj[11] = 1.0f;
  mProj[12] = 0.0f; mProj[13] = 0.0f; mProj[14] = (float)((2.0*zNear*zFar)/(zNear-zFar)) ; mProj[15] = 0.0f;

  glLoadMatrixf(mProj);

}

// -------------------------------------------

void GLToolkit::LookAtLH(double xEye,double yEye,double zEye,
                             double xAt,double yAt,double zAt,
                             double xUp,double yUp,double zUp) {

  // Create and load a left handed view matrix
  double Zx = xAt - xEye;
  double Zy = yAt - yEye;
  double Zz = zAt - zEye;
  double nZ = sqrt( Zx*Zx + Zy*Zy + Zz*Zz );
  Zx /= nZ;
  Zy /= nZ;
  Zz /= nZ;

  double Xx = (yUp*Zz - zUp*Zy);
  double Xy = (zUp*Zx - xUp*Zz);
  double Xz = (xUp*Zy - yUp*Zx);
  double nX = sqrt( Xx*Xx + Xy*Xy + Xz*Xz );
  Xx /= nX;
  Xy /= nX;
  Xz /= nX;

  double Yx = (Zy*Xz - Zz*Xy);
  double Yy = (Zz*Xx - Zx*Xz);
  double Yz = (Zx*Xy - Zy*Xx);

  double dotXE = xEye*Xx + yEye*Xy + zEye*Xz;
  double dotYE = xEye*Yx + yEye*Yy + zEye*Yz;
  double dotZE = xEye*Zx + yEye*Zy + zEye*Zz;

  GLfloat mView[16];

  mView[0] = (float)Xx; mView[1] = (float)Yx; mView[2] = (float)Zx; mView[3] = 0.0f;
  mView[4] = (float)Xy; mView[5] = (float)Yy; mView[6] = (float)Zy; mView[7] = 0.0f;
  mView[8] = (float)Xz; mView[9] = (float)Yz; mView[10] = (float)Zz; mView[11] = 0.0f;
  mView[12] = (float)-dotXE; mView[13] = (float)-dotYE; mView[14] = (float)-dotZE; mView[15] = 1.0f;

  glLoadMatrixf(mView);
                             
}

// -------------------------------------------

void GLToolkit::SetMaterial(GLMATERIAL *mat) {

  float acolor[] = { mat->Ambient.r, mat->Ambient.g, mat->Ambient.b, mat->Ambient.a };
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, acolor);
  float dcolor[] = { mat->Diffuse.r, mat->Diffuse.g, mat->Diffuse.b, mat->Diffuse.a };
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, dcolor);
  float scolor[] = { mat->Specular.r, mat->Specular.g, mat->Specular.b, mat->Specular.a };
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, scolor);
  float ecolor[] = { mat->Emissive.r, mat->Emissive.g, mat->Emissive.b, mat->Emissive.a };
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, ecolor);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, mat->Power);
  glColor4f(mat->Ambient.r, mat->Ambient.g, mat->Ambient.b, mat->Ambient.a);

}

// -------------------------------------------

void GLToolkit::printGlError(GLenum glError) {

  switch( glError ) {
    case GL_INVALID_ENUM:
      Log("OpenGL failure: An unacceptable value is specified for an enumerated argument. The offending function is ignored, having no side effect other than to set the error flag.\n");
      break;
    case GL_INVALID_VALUE:
      Log("OpenGL failure: A numeric argument is out of range. The offending function is ignored, having no side effect other than to set the error flag.\n");
      break;
    case GL_INVALID_OPERATION:
      Log("OpenGL failure: The specified operation is not allowed in the current state. The offending function is ignored, having no side effect other than to set the error flag.\n");
      break;
    case GL_STACK_OVERFLOW:
      Log("OpenGL failure: This function would cause a stack overflow. The offending function is ignored, having no side effect other than to set the error flag.\n");
      break;
    case GL_STACK_UNDERFLOW:
      Log("OpenGL failure: This function would cause a stack underflow. The offending function is ignored, having no side effect other than to set the error flag.\n");
      break;
    case GL_OUT_OF_MEMORY:
      Log("OpenGL failure: There is not enough memory left to execute the function. The state of OpenGL is undefined, except for the state of the error flags, after this error is recorded.\n");
      break;
  }
}

void GLToolkit::CheckGLErrors(char *compname) {
	static BOOL savedOnGLError=FALSE;   
	GLenum glError = glGetError();
       if( glError != GL_NO_ERROR ) { 
		   char tmp[256];
		   sprintf(tmp,"OpenGL: Error painting %s.",compname);
         GLToolkit::Log(tmp);
         GLToolkit::printGlError(glError); 
         //Exit();
#ifdef _DEBUG //stop 
		 __debugbreak();
#endif
		 /*if (!savedOnGLError) {
			 SynRad *mApp = (SynRad *)theApp;
			 mApp->AutoSave();
			 savedOnGLError = TRUE;
		 }*/
		 throw Error(tmp);
       }
}