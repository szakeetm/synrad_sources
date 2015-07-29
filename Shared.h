/*
  File:        Shared.h
  Description: Shared memory structure (inter process communication)
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

#include "Types.h"
#include <Windows.h>
#include <vector>
#include "Tools.h"

#ifndef SHAREDH
#define SHAREDH

#define NBHLEAK     2048  // Leak history max length
#define NBHHIT      2048  // Max. displayed number of lines and Porto (OPO)hits.
#define MAX_PROCESS 16    // Maximum number of process

typedef  struct {
    // Counts
    double fluxAbs,powerAbs;
    llong nbHit;           // Number of hits
    llong nbAbsorbed;      // Number of absorbed molec
	llong nbDesorbed;
  } SHHITS;

typedef struct {

  SHHITS total;               // Global counts
  llong  nbLeakTotal;         // Total leaks
  int    nbLastLeaks;         // Last leaks
  int    nbHHit;              // Last hits
  HIT    pHits[NBHHIT];       // Hit history
  LEAK   pLeak[NBHLEAK];      // Leak history
  llong  minHit_MC,maxHit_MC;
  double   minHit_flux;              // Minimum on texture
  double   maxHit_flux;              // Maximum on texture
  double   minHit_power;              // Minimum on texture
  double   maxHit_power;              // Maximum on texture
  double distTraveledTotal;
} SHGHITS;

// -----------------------------------------------------------------
// Master control shared memory block  (name: SRDCTRL[masterPID])
// 
// -----------------------------------------------------------------
#define PROCESS_STARTING 0   // Loading state
#define PROCESS_RUN      1   // Running state
#define PROCESS_READY    2   // Waiting state
#define PROCESS_KILLED   3   // Process killed
#define PROCESS_ERROR    4   // Process in error
#define PROCESS_DONE     5   // Simulation ended

#define COMMAND_NONE     10  // No change
#define COMMAND_LOAD     11  // Load geometry
#define COMMAND_START    12  // Start simu
#define COMMAND_PAUSE     13  // Pause simu
#define COMMAND_RESET    14  // Reset simu
#define COMMAND_EXIT     15  // Exit
#define COMMAND_CLOSE    16  // Release handles

static const char *prStates[] = {

"Not started",
"Running",
"Waiting",
"Killed",
"Error",
"Done",
"",
"",
"",
"No command",
"Load",
"Start",
"Stop",
"Reset",
"Exit",
"Close",

};

typedef struct {
  
  // Process control
  int    states[MAX_PROCESS];        // Process states/commands
  int    cmdParam[MAX_PROCESS];      // Command param 1
  llong  cmdParam2[MAX_PROCESS];     // Command param 2
  char   statusStr[MAX_PROCESS][64]; // Status message
} SHMASTER;

typedef struct {

  int        nbFacet;   // Number of facets (total)
  int        nbVertex;  // Number of 3D vertices
  int        nbSuper;   // Number of superstructures
  int        nbRegion;  //number of magnetic regions
  int        nbTrajPoints; //total number of trajectory points (calculated at CopyGeometryBuffer)
  int        generation_mode; //fluxwise or powerwise
  BOOL		 lowFluxMode;
  double	 lowFluxCutoff;
  char       name[64];  // (Short file name)
} SHGEOM;

typedef struct {

  // Facet parameters
  double sticking;       // Sticking        (0=>reflection  , 1=>absorption)
  double opacity;        // opacity         (0=>transparent , 1=>opaque)
  double area;           // Facet area (m^2)
  double roughness;      // sigmaX/sigmaH, only used if reflection type is "material scattering"
  int    reflectType;    // Reflection type
  int    profileType;    // Profile type
  int    superIdx;       // Super structure index (Indexed from 0)
  int    superDest;      // Super structure destination index (Indexed from 1, 0=>current)
  int	 teleportDest;   // Teleport destination facet id (for periodic boundary condition) (Indexed from 1, 0=>none)
  BOOL   countAbs;       // Count absoprtion (MC texture)
  BOOL   countRefl;      // Count reflection (MC texture)
  BOOL   countTrans;     // Count transparent (MC texture)
  BOOL   countDirection; // Record avergare direction (MC texture)
  BOOL   hasSpectrum;    // Calculate energy spectrum (histogram)

  // Flags
  BOOL   is2sided;     // 2 sided
  BOOL   isProfile;    // Profile facet
  BOOL   isOpaque;     // Opacity != 0
  BOOL   isTextured;   // texture
  BOOL   isVolatile;   // Volatile facet (absorbtion facet which does not affect particule trajectory)

  // Global hit counters
  SHHITS counter;

  // Normal vector
  VERTEX3D    N;    // normalized
  VERTEX3D    Nuv;  // normal to (u,v) not normlized

  // Axis Aligned Bounding Box (AABB)
  AABB       bb;
  VERTEX3D   center;

  // Geometry
  int    nbIndex;   // Number of index/vertex
  double sign;      // Facet vertex rotation (see Facet::DetectOrientation())

  // Plane basis (O,U,V) (See Geometry::InitializeGeometry() for info)
  VERTEX3D   O;  // Origin
  VERTEX3D   U;  // U vector
  VERTEX3D   V;  // V vector
  VERTEX3D   nU; // Normalized U
  VERTEX3D   nV; // Normalized V

  // Hit/Abs/Des/Density recording on 2D texture map
  int    texWidth;    // Rounded texture resolution (U)
  int    texHeight;   // Rounded texture resolution (V)
  double texWidthD;   // Actual texture resolution (U)
  double texHeightD;  // Actual texture resolution (V)

  int   hitOffset;      // Hit address offset for this facet

} SHFACET;

// -----------------------------------------------------------------
// Mesh shared structure  (name: SRDLOAD[masterPID])
//
//  SHELEM

typedef struct {

  float   area;     // Area of element
  float   uCenter;  // Center coordinates
  float   vCenter;  // Center coordinates
  int     elemId;   // Element index (MESH array)
  BOOL    full;     // Element is full

} SHELEM;

#endif /* SHAREDH */
