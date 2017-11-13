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

#ifndef _SIMULATIONH_
#define _SIMULATIONH_

#define MAX_STRUCT 512
#define MAX_THIT 16384

#include "SynradTypes.h"
#include "SMP.h"
#include <vector>
#include "Vector.h"
#include "Region_mathonly.h"
#include "TruncatedGaussian\rtnorm.hpp"
#include "SynradDistributions.h"
#include <tuple>

// Local facet structure

typedef struct {

  FacetProperties sh;
  FacetHitBuffer counter;

  int      *indices;   // Indices (Reference to geometry vertex)
  Vector2d *vertices2; // Vertices (2D plane space, UV coordinates)
  llong     *hits_MC;      // Number of MC hits
  double    *hits_flux;  //absorbed SR flux
  double    *hits_power; //absorbed SR power
  double	 fullSizeInc; // 1/Texture FULL element area
  double   *inc;        // reciprocial of element area
  bool     *largeEnough; //cells that are NOT too small for autoscaling
  VHIT     *direction; // Direction field recording (average)
  //char     *fullElem;  // Direction field recording (only on full element) (WHY???)
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
  bool   hitted;
  bool   ready;         // Volatile state
  size_t    textureSize;   // Texture size (in bytes)
  size_t    profileSize;   // profile size (in bytes)
  size_t    directionSize; // direction field size (in bytes)
  size_t    spectrumSize;  // spectrum size in bytes

  size_t globalId; //Global index (to identify when superstructures are present)

  void RegisterTransparentPass(); //Allows one shared Intersect routine between MolFlow and Synrad

} SubprocessFacet;

// Temporary transparent hit
extern  SubprocessFacet     **THitCache; //Global variable

// Local simulation structure

typedef struct {

  int              nbFacet;  // Number of facet
  SubprocessFacet          **facets;   // Facet handles
  struct AABBNODE *aabbTree; // Structure AABB tree

} SUPERSTRUCT;

typedef struct {

  FacetHitBuffer tmpCount;            // Temporary number of hits (between 2 calls of UpdateMC)
  llong  desorptionLimit;       // Maximum number of desorption

  int    hitCacheSize;              // Last hits  
  HIT    hitCache[HITCACHESIZE];       // Last hit history

  size_t    nbLeakSinceUpdate;   // Leaks since last UpdateMC
  size_t	leakCacheSize;		// Leaks from regions with displayed photons since last UpdateMC (registered in cache)
  LEAK		leakCache[LEAKCACHESIZE];      // Leak cache since last UpdateMC

  llong totalDesorbed;        //total number of generated photons, for process state reporting and simulation end check

  // Geometry
  char        name[64];         // Global name
  size_t         nbVertex;         // Number of vertex
  size_t         totalFacet;       // Total number of facet
  Vector3d   *vertices3;        // Vertices
  size_t         nbSuper;          // Number of super structure
  size_t         nbRegion;
  size_t         nbMaterials;
  size_t         nbTrajPoints;
  size_t      sourceArea;       //number of trajectory points weighed by 1/dL
  //size_t         nbDistrPoints_MAG;
  size_t         nbDistrPoints_BXY;
  size_t         curStruct;        // Current structure
  int			teleportedFrom;
  size_t sourceRegionId;
  SUPERSTRUCT str[MAX_STRUCT];

  std::vector<Region_mathonly> regions;// Regions
  std::vector<Material> materials;//materials
  std::vector<std::vector<double>> psi_distro;
  std::vector<std::vector<std::vector<double>>> chi_distros;
  std::vector<std::vector<double>> parallel_polarization;

  SubprocessFacet* lastHitFacet;     // Last hitted facet. Pointer, not a reference, that way it can be NULL (example: desorption from the middle of the volume)
  //int sourceArea;  // Number of trajectory points
  double stepPerSec;  // Avg number of step per sec
  size_t textTotalSize;  // Texture total size
  size_t profTotalSize;  // Profile total size
  size_t dirTotalSize;   // Direction field total size
  size_t spectrumTotalSize; //Spectrums total size
  bool loadOK;        // Load OK flag
  bool lastUpdateOK;  // Last hit update timeout
  bool hasVolatile;   // Contains volatile facet
  //bool hasDirection;  // Contains direction field

  gsl_rng *gen; //rnd gen stuff

  // Particle coordinates (MC)
  Vector3d pPos;    // Position
  Vector3d pDir;    // Direction
  //int      nbPHit;  // Number of hit (current particle) //Uncommented as it had no function
  double   dF;  //Flux carried by photon
  double   dP;  //Power carried by photon
  double   energy; //energy of the generated photon
  double   distTraveledCurrentParticle; //Distance traveled by particle before absorption
  double   distTraveledSinceUpdate;

  double oriRatio;
  bool newReflectionModel;

  SHMODE mode; //Low flux, generation mode, photon cache display

} SIMULATION;

// Handle to simulation object
extern SIMULATION *sHandle;

// -- Macros ---------------------------------------------------

// -- Methods ---------------------------------------------------

void RecordHitOnTexture(SubprocessFacet& f,double dF,double dP);
void RecordDirectionVector(SubprocessFacet& f);
void ProfileFacet(SubprocessFacet *f, const double &dF, const double &dP, const double &E);
void InitSimulation();
void ClearSimulation();
void SetState(int state, const char *status, bool changeState = true, bool changeStatus = true);
void SetErrorSub(const char *msg);
bool LoadSimulation(Dataport *loader);
bool UpdateSimuParams(Dataport *loader);
bool StartSimulation();
void ResetSimulation();
bool SimulationRun();
bool SimulationMCStep(const size_t& nbStep);
std::tuple<double, std::vector<double>, bool> GetStickingProbability(const SubprocessFacet& collidedFacet, const double& theta);

bool DoOldRegularReflection(SubprocessFacet& collidedFacet, const double& stickingProbability, const int& reflType, const double& inTheta, const double& inPhi,
	const Vector3d& N_rotated, const Vector3d& nU_rotated, const Vector3d& nV_rotated);
bool DoLowFluxReflection(SubprocessFacet& collidedFacet, const double& stickingProbability, const int& reflType, const double& inTheta, const double& inPhi,
	const Vector3d& N_rotated = Vector3d(0, 0, 0), const Vector3d& nU_rotated = Vector3d(0, 0, 0), const Vector3d& nV_rotated = Vector3d(0, 0, 0)); //old or new model
int GetHardHitType(const double& stickingProbability, const std::vector<double>& materialReflProbabilities, const bool& complexScattering);
std::tuple<Vector3d, Vector3d, Vector3d> PerturbateSurface( const SubprocessFacet& collidedFacet, const double& sigmaRatio);
//std::tuple<double,double,double> GetDirComponents(const Vector3d & nU_rotated, const Vector3d & nV_rotated, const Vector3d & N_rotated);
void RecordHit(const int &type,const double &dF,const double &dP);
void RecordLeakPos();
bool StartFromSource();
void ComputeSourceArea();
void PerformBounce_new(SubprocessFacet& collidedFacet, const int &reflType, const double &inTheta, const double &inPhi);
bool PerformBounce_old(SubprocessFacet& collidedFacet, const int& reflType, const double& inTheta, const double& inPhi, const Vector3d& N_rotated, const Vector3d& nU_rotated, const Vector3d& nV_rotated);
bool VerifiedSpecularReflection(const SubprocessFacet& collidedFacet, const int& reflType, const double& inTheta, const double& inPhi, const Vector3d& nU_rotated, const Vector3d& nV_rotated, const Vector3d& N_rotated);
void Stick(SubprocessFacet& collidedFacet);
void PerformTeleport(const SubprocessFacet& collidedFacet);

void UpdateHits(Dataport *dpHit,int prIdx,DWORD timeout);
void UpdateMCHits(Dataport *dpHit,int prIdx,DWORD timeout);
void ResetTmpCounters();
//bool RaySphereIntersect(Vector3d *center, double radius, Vector3d *rPos, Vector3d *rDir, double *dist);
double GetTick();
size_t   GetHitsSize();

#endif /* _SIMULATIONH_ */