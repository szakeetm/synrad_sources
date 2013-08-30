/*
   File:        GLMatrix.h
  Description: 4x4 Matrix (SDL/OpenGL OpenGL application framework)
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

#include <SDL_opengl.h>

#ifndef _GLMATRIXH_
#define _GLMATRIXH_

class GLMatrix {

public:

  GLMatrix();

  void Init33(float _11,float _12,float _13,
              float _21,float _22,float _23,
              float _31,float _32,float _33);

  void Identity();
  void RotateX(float angle);
  void RotateY(float angle);
  void RotateZ(float angle);
  void Translate(float x,float y,float z);
  void Transpose();

  void Multiply(GLMatrix *m1,GLMatrix *m2);
  void Multiply(GLMatrix *m2);

  void TransfomVec(float x,float y,float z,float w,float *rx,float *ry,float *rz,float *rw);

  GLfloat *GetGL();
  void LoadGL(GLfloat *m);
  void LoadGL(GLenum which);

  float _11;float _12;float _13;float _14;
  float _21;float _22;float _23;float _24;
  float _31;float _32;float _33;float _34;
  float _41;float _42;float _43;float _44;

};

#endif /* _GLMATRIXH_ */
