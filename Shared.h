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

#include <Windows.h>
#include <vector>
#include "GLApp/GLTypes.h"
#include "SynradTypes.h"
//#include "Worker.h"

#ifndef SHAREDH
#define SHAREDH

#define SPECTRUM_SIZE 100 //number of histogram bins
#define PROFILE_SIZE  100 // Size of profile
#define LEAKCACHESIZE     2048  // Leak history max length
#define HITCACHESIZE      2048  // Max. displayed number of lines and Porto (OPO)hits.
#define MAX_PROCESS 32    // Maximum number of process

typedef struct {

	Vector3d pos;
	int    type;
	double dF;
	double dP;
} HIT;

// Velocity field
typedef struct {
	Vector3d dir;
	llong count;
} VHIT;

typedef struct {

	Vector3d pos;
	Vector3d dir;

} LEAK;

typedef  struct {
    // Counts
    double fluxAbs,powerAbs;
    llong nbHit;           // Number of hits
    llong nbAbsorbed;      // Number of absorbed molec
	llong nbDesorbed;
  } SHHITS;

typedef struct {

  SHHITS total;               // Global counts
  size_t hitCacheSize;              // Number of valid hits in cache
  size_t lastHitIndex;					//Index of last recorded hit in gHits (turns over when reaches HITCACHESIZE)
  HIT    hitCache[HITCACHESIZE];       // Hit history
  
  size_t  lastLeakIndex;		  //Index of last recorded leak in gHits (turns over when reaches LEAKCACHESIZE)
  size_t  leakCacheSize;        //Number of valid leaks in the cache
  size_t  nbLeakTotal;         // Total leaks
  LEAK   leakCache[LEAKCACHESIZE];      // Leak history

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
#define COMMAND_PAUSE    13  // Pause simu
#define COMMAND_RESET    14  // Reset simu
#define COMMAND_EXIT     15  // Exit
#define COMMAND_CLOSE    16  // Release handles
#define COMMAND_UPDATEPARAMS 17 //Update simulation mode (low flux, fluxwise/powerwise, displayed regions)

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
"",
"No command",
"Loading",
"Starting",
"Stopping",
"Resetting",
"Exiting",
"Closing",
"Update params"
};

typedef struct {
  
  // Process control
  int    states[MAX_PROCESS];        // Process states/commands
  int    cmdParam[MAX_PROCESS];      // Command param 1
  llong  cmdParam2[MAX_PROCESS];     // Command param 2
  char   statusStr[MAX_PROCESS][64]; // Status message
} SHCONTROL;

typedef struct {

  size_t        nbFacet;   // Number of facets (total)
  size_t        nbVertex;  // Number of 3D vertices
  int        nbSuper;   // Number of superstructures
  size_t        nbRegion;  //number of magnetic regions
  size_t        nbTrajPoints; //total number of trajectory points (calculated at CopyGeometryBuffer)
  BOOL       newReflectionModel;
  char       name[64];  // (Short file name)
} SHGEOM; //Shared memory interface with main program

typedef struct {
	int      generation_mode; //fluxwise or powerwise
	BOOL	 lowFluxMode;
	double	 lowFluxCutoff;
} SHMODE;

typedef struct {

  // Facet parameters
  double sticking;       // Sticking        (0=>reflection  , 1=>absorption)
  double opacity;        // opacity         (0=>transparent , 1=>opaque)
  double area;           // Facet area (m^2)
  int    doScattering;   // Do rough surface scattering
  double rmsRoughness;   // RMS height roughness, in meters
  double autoCorrLength; // Autocorrelation length, in meters
  int    reflectType;    // Reflection type. 0=Diffuse, 1=Mirror, 10,11,12... : Material 0, Material 1, Material 2...., 9:invalid 
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
  //BOOL   isOpaque;     // Opacity != 0
  BOOL   isTextured;   // texture
  BOOL   isVolatile;   // Volatile facet (absorbtion facet which does not affect particule trajectory)

  // Normal vector
  Vector3d    N;    // normalized
  Vector3d    Nuv;  // normal to (u,v) not normlized

  // Axis Aligned Bounding Box (AABB)
  AABB       bb;
  Vector3d   center;

  // Geometry
  int    nbIndex;   // Number of index/vertex
  double sign;      // Facet vertex rotation (see Facet::DetectOrientation())

  // Plane basis (O,U,V) (See Geometry::InitializeGeometry() for info)
  Vector3d   O;  // Origin
  Vector3d   U;  // U vector
  Vector3d   V;  // V vector
  Vector3d   nU; // Normalized U
  Vector3d   nV; // Normalized V

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

// Hit type

#define HIT_DES   1
#define HIT_ABS   2
#define HIT_REF   3
#define HIT_TRANS 4
#define HIT_TELEPORT 5
#define LASTHIT 6

//Reflection type
#define REFL_ABSORB 0
#define REFL_FORWARD 1
#define REFL_DIFFUSE 2
#define REFL_BACK 3
#define REFL_TRANS 4

// Reflection type
#define REF_DIFFUSE 0   // Diffuse (cosine law)
#define REF_MIRROR  1   // Mirror
#define REF_MATERIAL 10 //Real rough surface reflection

// Profile type

#define REC_NONE       0  // No recording
#define REC_PRESSUREU  1  // Pressure profile (U direction)
#define REC_PRESSUREV  2  // Pressure profile (V direction)
#define REC_ANGULAR    3  // Angular profile

// Density/Hit field stuff
#define HITMAX_INT64 18446744073709551615
#define HITMAX_DOUBLE 1E308

#endif /* SHAREDH */
