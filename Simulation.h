/*
  File:        Simulation.h
  Description: Monte-Carlo Simulation for UHV
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

#define PI 3.14159265358979323846
typedef int BOOL;
#define MAX(x,y) (((x)<(y))?(y):(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define TRUE  1
#define FALSE 0
#define SAFE_FREE(x) if(x) { free(x);x=NULL; }
#define SATURATE(x,min,max) {if(x<(min)) x=(min); if(x>(max)) x=(max);}
#define MAX_STRUCT 512

#include "Shared.h"
#include "smp/SMP.h"
#include "Tools.h"
#include <vector>
#include "Tools.h"
#include "Region_mathonly.h"

#ifndef _SIMULATIONH_
#define _SIMULATIONH_

// Local facet structure

typedef struct {

  SHFACET sh;

  int      *indices;   // Indices (Reference to geometry vertex)
  VERTEX2D *vertices2; // Vertices (2D plane space, UV coordinates)
  llong     *hits_MC;      // Number of MC hits
  double    *hits_flux;  //absorbed SR flux
  double    *hits_power; //absorbed SR power
  double	 fullSizeInc; // 1/Texture FULL element area
  double   *inc;        // reciprocial of element area
  BOOL     *largeEnough; //cells that are NOT too small for autoscaling
  VHIT     *direction; // Direction field recording (average)
  char     *fullElem;  // Direction field recording (only on full element) (WHY???)
  llong    *profile_hits;   // MC hits
  double   *profile_flux;   // SR Flux
  double   *profile_power;  // SR power
  Histogram    *spectrum_fluxwise; //Energy spectrum fluxwise
  Histogram    *spectrum_powerwise;

  // Temporary var (used in Intersect for collision)
  double colDist;
  double colU;
  double colV;
  double rw;
  double rh;
  double iw;
  double ih;

  // Temporary var (used in FillHit for hit recording)
  BOOL   hitted;
  BOOL   ready;         // Volatile state
  int    textureSize;   // Texture size (in bytes)
  int    profileSize;   // profile size (in bytes)
  int    directionSize; // direction field size (in bytes)
  int    spectrumSize;  // spectrum size in bytes

  int globalId; //Global index (to identify when superstructures are present)

} FACET;

// Temporary transparent hit
#define MAX_THIT    16384
extern  FACET     **THits;

// AABBTree node

struct AABBNODE {

  AABB             bb;
  struct AABBNODE *left;
  struct AABBNODE *right;
  FACET          **list;
  int              nbFacet;

};

// Local simulation structure

typedef struct {

  int              nbFacet;  // Number of facet
  FACET          **facets;   // Facet handles
  struct AABBNODE *aabbTree; // Structure AABB tree

} SUPERSTRUCT;

typedef struct {

  SHHITS tmpCount;            // Temporary number of hits (between 2 updates)
  llong nbLeakTotal;          // Total number of unexpected leak (simulation error)
  int    nbLastLeak;          // Last leaks
  int    nbHHit;              // Last hits
  llong  maxDesorption;       // Maximum number of desorption
  HIT    pHits[NBHHIT];       // Last hit history
  LEAK   pLeak[NBHLEAK];      // Leak history
  //llong  wallHits[BOUNCEMAX]; // 'Wall collision count before absoprtion' density histogram
  llong totalDesorbed;        //total number of generated photons, for process state reporting and simulation end check
  int   generation_mode;      //fluxwise or powerwise

  // Geometry
  char        name[64];         // Global name
  int         nbVertex;         // Number of vertex
  int         totalFacet;       // Total number of facet
  VERTEX3D   *vertices3;        // Vertices
  int         nbSuper;          // Number of super structure
  size_t         nbRegion;
  size_t         nbMaterials;
  size_t         nbTrajPoints;
  double      sourceArea;       //number of trajectory points weighed by 1/dL
  size_t         nbDistrPoints_MAG;
  size_t         nbDistrPoints_BXY;
  int         curStruct;        // Current structure
  int			teleportedFrom;
  SUPERSTRUCT str[MAX_STRUCT];

  std::vector<Region_mathonly> regions;// Regions
  std::vector<Material> materials;//materials
  std::vector<std::vector<double>> psi_distr;
  std::vector<std::vector<double>> chi_distr;

  FACET *lastHit;     // Last hitted facet
  //int sourceArea;  // Number of trajectory points
  double stepPerSec;  // Avg number of step per sec
  size_t textTotalSize;  // Texture total size
  size_t profTotalSize;  // Profile total size
  size_t dirTotalSize;   // Direction field total size
  size_t spectrumTotalSize; //Spectrums total size
  BOOL loadOK;        // Load OK flag
  BOOL lastUpdateOK;  // Last hit update timeout
  BOOL hasVolatile;   // Contains volatile facet
  BOOL hasDirection;  // Contains direction field

  // Particle coordinates (MC)
  VERTEX3D pPos;    // Position
  VERTEX3D pDir;    // Direction
  //int      nbPHit;  // Number of hit (current particle) //Uncommented as it had no function
  double   dF;  //Flux carried by photon
  double   dP;  //Power carried by photon
  double   energy; //energy of the generated photon
  double   distTraveledCurrentParticle; //Distance traveled by particle before absorption
  double   distTraveledSinceUpdate;

  BOOL lowFluxMode;
  double lowFluxCutoff;
  double oriRatio;

} SIMULATION;

// Handle to simulation object
extern SIMULATION *sHandle;

// -- Macros ---------------------------------------------------




// -- Methods ---------------------------------------------------

void RecordHitOnTexture(FACET *f,double dF,double dP);
void InitSimulation();
void ClearSimulation();
BOOL LoadSimulation(Dataport *loader);
BOOL StartSimulation();
void ResetSimulation();
BOOL SimulationRun();
BOOL SimulationMCStep(int nbStep);
void RecordHit(const int &type,const double &dF,const double &dP);
void RecordLeakPos();
BOOL StartFromSource();
void ComputeSourceArea();
//int PerformBounce(FACET *iFacet,double sigmaRatio=0.0,double theta=0.0,double phi=0.0,
//	Vector N_rotated=Vector(0,0,0),Vector nU_rotated=Vector(0,0,0),Vector nV_rotated=Vector(0,0,0),double thetaOffset=0.0,double phiOffset=0.0, double randomAngle=0.0, int reflType=1);
void PerformBounce(FACET *iFacet, const double &theta, const double &phi, const int &reflType);
int Stick(FACET *collidedFacet);
void PerformTeleport(FACET *iFacet);
void PolarToCartesian(FACET *iFacet,double theta,double phi,BOOL reverse,double rotateUV=0.0);

void CartesianToPolar(FACET *iFacet,double *theta,double *phi);
void UpdateHits(Dataport *dpHit,int prIdx,DWORD timeout);
void UpdateMCHits(Dataport *dpHit,int prIdx,DWORD timeout);
void ResetCounter();
struct AABBNODE *BuildAABBTree(FACET **list,int nbFacet,int depth);
int FindBestCuttingPlane(struct AABBNODE *node,int *left,int *right);
void ComputeBB(struct AABBNODE *node);
void DestroyAABB(struct AABBNODE *node);
void IntersectTree(struct AABBNODE *node);
BOOL Intersect(VERTEX3D *rayPos,VERTEX3D *rayDir,double *dist,FACET **iFact,FACET *last);
BOOL Visible(VERTEX3D *c1,VERTEX3D *c2,FACET *f1,FACET *f2);
void ProfileFacet(FACET *f,const double &dF,const double &dP,const double &E);
BOOL IsInFacet(FACET *f,const double &u,const double &v);
double GetTick();
size_t   GetHitsSize();

#endif /* _SIMULATIONH_ */
