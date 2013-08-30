/*
   File:        GLSprite.cpp
  Description: 2D Sprites (SDL/OpenGL OpenGL application framework)
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

#include "GLSprite.h"
#include "GLToolkit.h"
#include "GLApp.h"
#include <CImage.h>
#include <malloc.h>
#include <stdio.h>

extern GLApplication *theApp;

// -------------------------------------------

Sprite2D::Sprite2D() {
  rC = 1.0f;
  gC = 1.0f;
  bC = 1.0f;
}

// -------------------------------------------

void Sprite2D::SetSpriteMapping(float mx1,float my1,float mx2,float my2) {

  this->mx1 = mx1;
  this->my1 = my1;
  this->mx2 = mx2;
  this->my2 = my2;

}

// -------------------------------------------

void Sprite2D::UpdateSprite(int x1,int y1,int x2,int y2,float mx1,float my1,float mx2,float my2) {

  this->x1 = x1;
  this->y1 = y1;
  this->x2 = x2;
  this->y2 = y2;
  this->mx1 = mx1;
  this->my1 = my1;
  this->mx2 = mx2;
  this->my2 = my2;

}

// -------------------------------------------

void Sprite2D::UpdateSprite(int x1,int y1,int x2,int y2) {

  this->x1 = x1;
  this->y1 = y1;
  this->x2 = x2;
  this->y2 = y2;

}

// -------------------------------------------

void Sprite2D::SetColor(float r,float g,float b) {
  rC = r;
  gC = g;
  bC = b;
}

// -------------------------------------------

int Sprite2D::RestoreDeviceObjects(char *diffName,char *alphaName,int scrWidth,int scrHeight) {

  GLint  bpp;
  GLenum format;
  BYTE *buff32;
  char tmp[512];

  // Load the image
  CImage img;
  CImage imga;

  
  
  if( !img.LoadImage(diffName) ) {
    sprintf(tmp,"Failed to load \"%s\"",diffName);
    GLToolkit::Log(tmp); //Commented out for File open dialog
    return 0;
  }
  

  hasAlpha = strcmp(alphaName,"none")!=0;
  if( hasAlpha ) {
    if( !imga.LoadImage(alphaName) ) {
      sprintf(tmp,"Failed to load \"%s\"",alphaName);
      GLToolkit::Log(tmp);
      return 0;
    }
  }

  // Make 32 Bit RGB/RGBA buffer
  int fWidth  = img.Width();
  int fHeight = img.Height();

  if( hasAlpha ) {

    format = GL_RGBA;
    bpp    = 4;
    buff32 = (BYTE *)malloc(fWidth*fHeight*4);
    BYTE *data   = img.GetData();
    BYTE *adata   = imga.GetData();

    for(int y=0;y<fHeight;y++) {
      for(int x=0;x<fWidth;x++) {
        buff32[x*4 + 0 + y*4*fWidth] = data[x*3+2 + y*3*fWidth];
        buff32[x*4 + 1 + y*4*fWidth] = data[x*3+1 + y*3*fWidth];
        buff32[x*4 + 2 + y*4*fWidth] = data[x*3+0 + y*3*fWidth];
        buff32[x*4 + 3 + y*4*fWidth] = adata[x*3+1 + y*3*fWidth];
      }
    }

  } else {

    format = GL_RGB;
    bpp    = 3;
    buff32 = (BYTE *)malloc(fWidth*fHeight*3);
    BYTE *data   = img.GetData();

    for(int y=0;y<fHeight;y++) {
      for(int x=0;x<fWidth;x++) {
        buff32[x*3 + 0 + y*3*fWidth] = data[x*3+2 + y*3*fWidth];
        buff32[x*3 + 1 + y*3*fWidth] = data[x*3+1 + y*3*fWidth];
        buff32[x*3 + 2 + y*3*fWidth] = data[x*3+0 + y*3*fWidth];
      }
    }

  }

  glGenTextures(1,&texId);
  glBindTexture(GL_TEXTURE_2D,texId);

  glTexImage2D (
    GL_TEXTURE_2D,       // Type
    0,                   // No Mipmap
    bpp,                // Byte per pixel
    fWidth,             // Width
    fHeight,            // Height
    0,                   // Border
    format,             // Format RGB/RGBA
    GL_UNSIGNED_BYTE,   // 8 Bit/color
    buff32              // Data
  );   

  free(buff32);
  img.Release();
  imga.Release();

  GLenum glError = glGetError();
  if( glError != GL_NO_ERROR )
  {
    sprintf(tmp,"Failed to create Sprite2D \"%s\"",diffName);
    GLToolkit::Log(tmp);
    GLToolkit::printGlError(glError);
    return 0;    
  }

  // Compute othographic matrix (for Transfomed Lit vertex)
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glOrtho( 0, scrWidth, scrHeight, 0, -1, 1 );
  glGetFloatv( GL_PROJECTION_MATRIX , pMatrix );

  return 1;

}

// -------------------------------------------

void Sprite2D::InvalidateDeviceObjects() {

  if(texId) glDeleteTextures(1, &texId);
  texId = 0;

}

// -------------------------------------------

void Sprite2D::Render(BOOL doLinear) {

  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D,texId);
  if( doLinear ) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  } else {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  }
  if( hasAlpha ) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  } else {
    glDisable(GL_BLEND);
  }
  glColor3f(rC,gC,bC);
  /*
  if(loadMatrix) {
    glMatrixMode( GL_PROJECTION );
    glLoadMatrixf(pMatrix);
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
  }
  */
  glBegin(GL_QUADS);
  glTexCoord2f(mx1,my1);glVertex2i(x1,y1);
  glTexCoord2f(mx2,my1);glVertex2i(x2,y1);
  glTexCoord2f(mx2,my2);glVertex2i(x2,y2);
  glTexCoord2f(mx1,my2);glVertex2i(x1,y2);
  glEnd();

#ifdef _DEBUG
  theApp->nbPoly++;
#endif

}

// -------------------------------------------

void Sprite2D::Render90(BOOL doLinear) {

  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D,texId);
  if( doLinear ) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  } else {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  }
  if( hasAlpha ) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  } else {
    glDisable(GL_BLEND);
  }
  glColor3f(rC,gC,bC);
  /*
  if(loadMatrix) {
    glMatrixMode( GL_PROJECTION );
    glLoadMatrixf(pMatrix);
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
  }
  */
  glBegin(GL_QUADS);
  glTexCoord2f(mx1,my1);glVertex2i(x1,y1);
  glTexCoord2f(mx1,my2);glVertex2i(x2,y1);
  glTexCoord2f(mx2,my2);glVertex2i(x2,y2);
  glTexCoord2f(mx2,my1);glVertex2i(x1,y2);
  glEnd();

#ifdef _DEBUG
  theApp->nbPoly++;
#endif

}
