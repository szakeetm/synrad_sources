/*
   File:        GLMatrix.cpp
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

#include "GLMatrix.h"
#include <math.h>

// -------------------------------------

GLMatrix::GLMatrix() {
  Identity();
}

// -------------------------------------

void GLMatrix::Init33(float _11,float _12,float _13,float _21,float _22,float _23,float _31,float _32,float _33) {

  this->_11 = _11; this->_12 = _12; this->_13 = _13;
  this->_21 = _21; this->_22 = _22; this->_23 = _23;
  this->_31 = _31; this->_32 = _32; this->_33 = _33;

}

// -------------------------------------
  
void GLMatrix::Identity() {

  _11=1.0f;  _12=0.0f;  _13=0.0f; _14=0.0f;
  _21=0.0f;  _22=1.0f;  _23=0.0f; _24=0.0f;
  _31=0.0f;  _32=0.0f;  _33=1.0f; _34=0.0f;
  _41=0.0f;  _42=0.0f;  _43=0.0f; _44=1.0f;

}

// -------------------------------------

void GLMatrix::Translate(float x,float y,float z) {

  Identity();
  _14 = x;
  _24 = y;
  _34 = z;

}
    
// -------------------------------------

void GLMatrix::Multiply(GLMatrix *m2) {

  GLMatrix m1 = *this;
  
  this->_11 = m1._11*m2->_11 + m1._12*m2->_21 + m1._13*m2->_31 + m1._14*m2->_41;
  this->_12 = m1._11*m2->_12 + m1._12*m2->_22 + m1._13*m2->_32 + m1._14*m2->_42;
  this->_13 = m1._11*m2->_13 + m1._12*m2->_23 + m1._13*m2->_33 + m1._14*m2->_43;
  this->_14 = m1._11*m2->_14 + m1._12*m2->_24 + m1._13*m2->_34 + m1._14*m2->_44;

  this->_21 = m1._21*m2->_11 + m1._22*m2->_21 + m1._23*m2->_31 + m1._24*m2->_41;
  this->_22 = m1._21*m2->_12 + m1._22*m2->_22 + m1._23*m2->_32 + m1._24*m2->_42;
  this->_23 = m1._21*m2->_13 + m1._22*m2->_23 + m1._23*m2->_33 + m1._24*m2->_43;
  this->_24 = m1._21*m2->_14 + m1._22*m2->_24 + m1._23*m2->_34 + m1._24*m2->_44;

  this->_31 = m1._31*m2->_11 + m1._32*m2->_21 + m1._33*m2->_31 + m1._34*m2->_41;
  this->_32 = m1._31*m2->_12 + m1._32*m2->_22 + m1._33*m2->_32 + m1._34*m2->_42;
  this->_33 = m1._31*m2->_13 + m1._32*m2->_23 + m1._33*m2->_33 + m1._34*m2->_43;
  this->_34 = m1._31*m2->_14 + m1._32*m2->_24 + m1._33*m2->_34 + m1._34*m2->_44;

  this->_41 = m1._41*m2->_11 + m1._42*m2->_21 + m1._43*m2->_31 + m1._44*m2->_41;
  this->_42 = m1._41*m2->_12 + m1._42*m2->_22 + m1._43*m2->_32 + m1._44*m2->_42;
  this->_43 = m1._41*m2->_13 + m1._42*m2->_23 + m1._43*m2->_33 + m1._44*m2->_43;
  this->_44 = m1._41*m2->_14 + m1._42*m2->_24 + m1._43*m2->_34 + m1._44*m2->_44;

}

void GLMatrix::Multiply(GLMatrix *m1,GLMatrix *m2) {

  this->_11 = m1->_11*m2->_11 + m1->_12*m2->_21 + m1->_13*m2->_31 + m1->_14*m2->_41;
  this->_12 = m1->_11*m2->_12 + m1->_12*m2->_22 + m1->_13*m2->_32 + m1->_14*m2->_42;
  this->_13 = m1->_11*m2->_13 + m1->_12*m2->_23 + m1->_13*m2->_33 + m1->_14*m2->_43;
  this->_14 = m1->_11*m2->_14 + m1->_12*m2->_24 + m1->_13*m2->_34 + m1->_14*m2->_44;

  this->_21 = m1->_21*m2->_11 + m1->_22*m2->_21 + m1->_23*m2->_31 + m1->_24*m2->_41;
  this->_22 = m1->_21*m2->_12 + m1->_22*m2->_22 + m1->_23*m2->_32 + m1->_24*m2->_42;
  this->_23 = m1->_21*m2->_13 + m1->_22*m2->_23 + m1->_23*m2->_33 + m1->_24*m2->_43;
  this->_24 = m1->_21*m2->_14 + m1->_22*m2->_24 + m1->_23*m2->_34 + m1->_24*m2->_44;

  this->_31 = m1->_31*m2->_11 + m1->_32*m2->_21 + m1->_33*m2->_31 + m1->_34*m2->_41;
  this->_32 = m1->_31*m2->_12 + m1->_32*m2->_22 + m1->_33*m2->_32 + m1->_34*m2->_42;
  this->_33 = m1->_31*m2->_13 + m1->_32*m2->_23 + m1->_33*m2->_33 + m1->_34*m2->_43;
  this->_34 = m1->_31*m2->_14 + m1->_32*m2->_24 + m1->_33*m2->_34 + m1->_34*m2->_44;

  this->_41 = m1->_41*m2->_11 + m1->_42*m2->_21 + m1->_43*m2->_31 + m1->_44*m2->_41;
  this->_42 = m1->_41*m2->_12 + m1->_42*m2->_22 + m1->_43*m2->_32 + m1->_44*m2->_42;
  this->_43 = m1->_41*m2->_13 + m1->_42*m2->_23 + m1->_43*m2->_33 + m1->_44*m2->_43;
  this->_44 = m1->_41*m2->_14 + m1->_42*m2->_24 + m1->_43*m2->_34 + m1->_44*m2->_44;
  
}

// -------------------------------------

void GLMatrix::TransfomVec(float x,float y,float z,float w,float *rx,float *ry,float *rz,float *rw) {

  *rx = x * _11 + y * _12 + z * _13 + w * _14;
  *ry = x * _21 + y * _22 + z * _23 + w * _24;
  *rz = x * _31 + y * _32 + z * _33 + w * _34;
  *rw = x * _41 + y * _42 + z * _43 + w * _44;

}

// -------------------------------------

GLfloat *GLMatrix::GetGL() {
  
  static GLfloat ret[16];
  ret[0] =  _11; ret[4] =  _12; ret[8]  = _13;  ret[12] = _14; 
  ret[1] =  _21; ret[5] =  _22; ret[9]  = _23;  ret[13] = _24; 
  ret[2] =  _31; ret[6] =  _32; ret[10] = _33;  ret[14] = _34; 
  ret[3] =  _41; ret[7] =  _42; ret[11] = _43;  ret[15] = _44; 
  return ret;

}

// -------------------------------------

void GLMatrix::LoadGL(GLfloat *m) {
  
  _11 = m[0]; _12 = m[4]; _13 = m[8];  _14 = m[12]; 
  _21 = m[1]; _22 = m[5]; _23 = m[9];  _24 = m[13]; 
  _31 = m[2]; _32 = m[6]; _33 = m[10]; _34 = m[14]; 
  _41 = m[3]; _42 = m[7]; _43 = m[11]; _44 = m[15]; 

}

// -------------------------------------

void GLMatrix::LoadGL(GLenum which) {

  GLfloat m[16];
  glGetFloatv(which,m);
  _11 = m[0]; _12 = m[4]; _13 = m[8];  _14 = m[12]; 
  _21 = m[1]; _22 = m[5]; _23 = m[9];  _24 = m[13]; 
  _31 = m[2]; _32 = m[6]; _33 = m[10]; _34 = m[14]; 
  _41 = m[3]; _42 = m[7]; _43 = m[11]; _44 = m[15]; 

}

// -------------------------------------

void GLMatrix::Transpose() {

  GLfloat m[16];
  m[0] =  _11; m[4] =  _12; m[8]  = _13;  m[12] = _14; 
  m[1] =  _21; m[5] =  _22; m[9]  = _23;  m[13] = _24; 
  m[2] =  _31; m[6] =  _32; m[10] = _33;  m[14] = _34; 
  m[3] =  _41; m[7] =  _42; m[11] = _43;  m[15] = _44; 

  _11 = m[0];  _12 = m[1];  _13 = m[2];  _14 = m[3]; 
  _21 = m[4];  _22 = m[5];  _23 = m[6];  _24 = m[7]; 
  _31 = m[8];  _32 = m[9];  _33 = m[10]; _34 = m[11]; 
  _41 = m[12]; _42 = m[13]; _43 = m[14]; _44 = m[15]; 

}

// -------------------------------------

void GLMatrix::RotateX(float angle) {

  float cs = cosf(angle);
  float sn = sinf(angle);

  Identity();
  _22=cs;   _23=-sn;
  _32=sn;   _33=cs;
  
}

// -------------------------------------

void GLMatrix::RotateY(float angle) {

  float cs = cosf(angle);
  float sn = sinf(angle);

  Identity();
  _11=cs;  _13=sn;
  _31=-sn; _33=cs;
  
}

// -------------------------------------

void GLMatrix::RotateZ(float angle) {

  float cs = cosf(angle);
  float sn = sinf(angle);

  Identity();
  _11=cs;   _12=-sn;
  _21=sn;   _22=cs;
  
}
