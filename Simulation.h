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

#include "SynradTypes.h"
#include "Buffer_shared.h"
#include "SMP.h"
#include <vector>
#include "Vector.h"
#include "Region_mathonly.h"
#include "TruncatedGaussian\rtnorm.hpp"
#include "SynradDistributions.h"
#include <tuple>

// Local facet structure

class SubprocessFacet {
public:
	FacetProperties sh;
	FacetHitBuffer tmpCounter;

    std::vector<size_t> indices;   // Indices (Reference to geometry vertex)
    std::vector<Vector2d> vertices2; // Vertices (2D plane space, UV coordinates)
    std::vector<TextureCell> texture;

    double	 fullSizeInc; // 1/Texture FULL element area
    std::vector<double> textureCellIncrements;        // reciprocial of element area
    std::vector<bool> largeEnough; //cells that are NOT too small for autoscaling
    std::vector<DirectionCell>direction; // Direction field recording (average)
	//char     *fullElem;  // Direction field recording (only on full element) (WHY???)
    std::vector<ProfileSlice> profile;
    Histogram spectrum;

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
	void ResetCounter();
    void	ResizeCounter(size_t nbMoments);
    bool  InitializeOnLoad(const size_t& globalId);

    void InitializeHistogram();

    bool InitializeDirectionTexture();

    bool InitializeProfile();

    bool InitializeTexture();

    bool InitializeSpectrum();

    bool InitializeLinkAndVolatile(const size_t & id);
};


class AABBNODE;
// Local simulation structure
class SuperStructure {
public:
    SuperStructure();
    ~SuperStructure();
    std::vector<SubprocessFacet> facets;   // Facet handles
	AABBNODE *aabbTree; // Structure AABB tree
} ;

class CurrentParticleStatus {
public:



    Vector3d position;    // Position
    Vector3d direction;    // Direction
    double oriRatio; //Represented ratio of desorbed, used for low flux mode

    //Recordings for histogram
    size_t   nbBounces; // Number of hit (current particle) since desorption
    double   distanceTraveled;
    double   distTraveledSinceUpdate;

    double   dF;  //Flux carried by photon
    double   dP;  //Power carried by photon
    double   energy; //energy of the generated photon

    size_t   structureId;        // Current structure
    int      teleportedFrom;   // We memorize where the particle came from: we can teleport back
    SubprocessFacet *lastHitFacet;     // Last hitted facet
    std::vector<SubprocessFacet*> transparentHitBuffer; //Storing this buffer simulation-wide is cheaper than recreating it at every Intersect() call
};

class Simulation {
public:
    Simulation();
    GlobalHitBuffer tmpGlobalResult;            // Temporary number of hits (between 2 calls of UpdateMC)

	size_t    nbLeakSinceUpdate;   // Leaks since last UpdateMC

	size_t totalDesorbed;        //total number of generated photons, for process state reporting and simulation end check

	// Geometry
	/*char name[64];         // Global name
	size_t nbVertex;         // Number of vertex
	size_t totalFacet;       // Total number of facet*/
    //size_t nbSuper;          // Number of super structure

    GeomProperties sh;
    WorkerParams wp;
    OntheflySimulationParams ontheflyParams; //Low flux, generation mode, photon cache display

    // to worker
    //size_t nbRegion;
    //size_t nbTrajPoints;
    //bool newReflectionModel;

    std::vector<Vector3d> vertices3;        // Vertices
	size_t nbMaterials;
	size_t sourceArea;       //number of trajectory points weighed by 1/dL
	//size_t nbDistrPoints_MAG;
	size_t nbDistrPoints_BXY;
	size_t sourceRegionId;
    std::vector<SuperStructure> structures; //They contain the facets

	std::vector<Region_mathonly> regions;// Regions
	std::vector<Material> materials;//materials
	std::vector<std::vector<double>> psi_distro;
	std::vector<std::vector<std::vector<double>>> chi_distros;
	std::vector<std::vector<double>> parallel_polarization;

	//SubprocessFacet* lastHitFacet;     // Last hitted facet. Pointer, not a reference, that way it can be NULL (example: desorption from the middle of the volume)
	//int sourceArea;  // Number of trajectory points
	double stepPerSec;  // Avg number of step per sec
	size_t textTotalSize;  // Texture total size
	size_t profTotalSize;  // Profile total size
	size_t dirTotalSize;   // Direction field total size
	size_t spectrumTotalSize; //Spectrums total size
	bool loadOK;        // Load OK flag
	bool lastHitUpdateOK;  // Last hit update timeout
	bool lastLogUpdateOK;  // Last log update timeout
	bool hasVolatile;   // Contains volatile facet
	//bool hasDirection;  // Contains direction field

	gsl_rng *gen; //rnd gen stuff

    // Particle coordinates (MC)
    CurrentParticleStatus currentParticle;


    std::vector<ParticleLoggerItem> tmpParticleLog;

};

// -- Macros ---------------------------------------------------

// -- Methods ---------------------------------------------------

void RecordHitOnTexture(SubprocessFacet& f, double dF, double dP);
void RecordDirectionVector(SubprocessFacet& f);
void LogHit(SubprocessFacet *f);
void ProfileFacet(SubprocessFacet &f, const double& energy, const ProfileSlice &increment);
void InitSimulation();
void ClearSimulation();
void SetState(size_t state, const char *status, bool changeState = true, bool changeStatus = true);
void SetErrorSub(const char *msg);
bool LoadSimulation(Dataport *loader);
bool UpdateOntheflySimuParams(Dataport *loader);
bool StartSimulation();
void ResetSimulation();
bool SimulationRun();
bool SimulationMCStep(const size_t& nbStep);
std::tuple<double, std::vector<double>, bool> GetStickingProbability(const SubprocessFacet& collidedFacet, const double& theta);

bool DoOldRegularReflection(SubprocessFacet& collidedFacet, const int& reflType, const double& inTheta, const double& inPhi,
	const Vector3d& N_rotated, const Vector3d& nU_rotated, const Vector3d& nV_rotated);
bool DoLowFluxReflection(SubprocessFacet& collidedFacet, const double& stickingProbability, const bool& complexScattering, const std::vector<double>& materialReflProbabilities,
	const double& inTheta, const double& inPhi,
	const Vector3d& N_rotated = Vector3d(0, 0, 0), const Vector3d& nU_rotated = Vector3d(0, 0, 0), const Vector3d& nV_rotated = Vector3d(0, 0, 0)); //old or new model
int GetHardHitType(const double& stickingProbability, const std::vector<double>& materialReflProbabilities, const bool& complexScattering);
std::tuple<Vector3d, Vector3d, Vector3d> PerturbateSurface(const SubprocessFacet& collidedFacet, const double& sigmaRatio);
//std::tuple<double,double,double> GetDirComponents(const Vector3d & nU_rotated, const Vector3d & nV_rotated, const Vector3d & N_rotated);
void RecordHit(const int &type, const double &dF, const double &dP);
void RecordLeakPos();
bool StartFromSource();
void ComputeSourceArea();
void PerformBounce_new(SubprocessFacet& collidedFacet, const int &reflType, const double &inTheta, const double &inPhi);
bool PerformBounce_old(SubprocessFacet& collidedFacet, const int& reflType, const double& inTheta, const double& inPhi, const Vector3d& N_rotated, const Vector3d& nU_rotated, const Vector3d& nV_rotated);
bool VerifiedSpecularReflection(const SubprocessFacet& collidedFacet, const int& reflType, const double& inTheta, const double& inPhi, const Vector3d& nU_rotated, const Vector3d& nV_rotated, const Vector3d& N_rotated);
void Stick(SubprocessFacet& collidedFacet);
void PerformTeleport(const SubprocessFacet& collidedFacet);
void UpdateHits(Dataport *dpHit, Dataport *dpLog, int prIdx, DWORD timeout);
void UpdateLog(Dataport *dpLog, DWORD timeout);
void UpdateMCHits(Dataport *dpHit, int prIdx, DWORD timeout);
void ResetTmpCounters();
//bool RaySphereIntersect(Vector3d *center, double radius, Vector3d *rPos, Vector3d *rDir, double *dist);
double GetTick();
size_t   GetHitsSize();

#endif /* _SIMULATIONH_ */