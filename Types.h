/*
  File:        Types.h
  Description: Various type/macro definitions and util functions
  Program:     SynRad
  Author:      R. KERSEVAN / M ADY / M ADY
  Copyright:   E.S.R.F / CERN / CERN

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef TYPESH
#define TYPESH

// 64 bit integer declaration

#ifdef WIN
  typedef unsigned __int64 llong;
#else
  typedef long long llong;
#endif

// AC matrix floating type

typedef float ACFLOAT;

// Desorption type

#define DES_NONE    0   // No desorption
#define DES_UNIFORM 1   // Uniform
#define DES_COSINE  2   // cos(theta)
/*#define DES_COSINE2 3   // cos(theta)^2
#define DES_COSINE3 4   // cos(theta)^3
#define DES_COSINE4 5   // cos(theta)^4*/
#define DES_COSINE_N 3 // cos(theta)^N

// Reflection type
#define REF_DIFFUSE 0   // Diffuse (cosine law)
#define REF_MIRROR  1   // Mirror
#define REF_MATERIAL 10 //Real rough surface reflection

  // Profile type

#define REC_NONE       0  // No recording
#define REC_PRESSUREU  1  // Pressure profile (U direction)
#define REC_PRESSUREV  2  // Pressure profile (V direction)
#define REC_ANGULAR    3  // Angular profile

#define PROFILE_SIZE  100 // Size of profile
#define SPECTRUM_SIZE 100 //number of histogram bins

// Hit type

#define HIT_DES   1
#define HIT_ABS   2
#define HIT_REF   3
#define HIT_TRANS 4
#define HIT_TELEPORT 5
#define LASTHIT 6

// Geometry structure definitions
typedef struct {

  double x;
  double y;
  double z;
  int selected;

} VERTEX3D;

/*typedef struct {
	VERTEX3D v1;
	VERTEX3D v2;
	VERTEX3D v3;
} MATRIX;*/

typedef struct {

  double u;
  double v;

} VERTEX2D;

typedef struct {

  VERTEX3D min;
  VERTEX3D max;

} AABB;

typedef struct {

  VERTEX3D pos;
  int      type;
  double dF;
  double dP;
} HIT;

typedef struct {

  VERTEX3D pos;
  VERTEX3D dir;

} LEAK;

typedef struct {
   VERTEX2D* points;
   size_t nbPoints;
   float   area;     // Area of element
   float   uCenter;  // Center coordinates
   float   vCenter;  // Center coordinates
   //int     elemId;   // Element index (MESH array)
   //int full;
 } CELLPROPERTIES;

// Density/Hit field stuff
#define HITMAX_INT64 18446744073709551615
#define HITMAX_DOUBLE 1E308


// Velocity field
typedef struct {
  VERTEX3D dir;
  llong count;
} VHIT;

#define IS_ZERO(x) (fabs((x))<1e-10)

#define DOT2(x1,y1,x2,y2) ((x1)*(x2) + (y1)*(y2))
#define DOT3(x1,y1,z1,x2,y2,z2) ((x1)*(x2) + (y1)*(y2) + (z1)*(z2))

#define DET22(_11,_12,_21,_22) ( (_11)*(_22) - (_21)*(_12) )
#define DET33(_11,_12,_13,_21,_22,_23,_31,_32,_33)  \
  ((_11)*( (_22)*(_33) - (_32)*(_23) ) +            \
   (_12)*( (_23)*(_31) - (_33)*(_21) ) +            \
   (_13)*( (_21)*(_32) - (_31)*(_22) ))


#endif /* TYPESH */
