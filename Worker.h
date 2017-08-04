/*
  File:        Worker.h
  Description: Sub processes handling
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

#ifndef _WORKERH_
#define _WORKERH_

#include "Region_full.h"
//#include "LoadStatus.h"
//#include "Distributions.h" //for materials
//#include "PugiXML/pugixml.hpp"
#include <string>
#include "GLApp/GLTypes.h"
#include "Shared.h"
#include "Smp.h"

typedef float ACFLOAT;

class Geometry;
class SynradGeometry;
class GLProgress;
class LoadStatus;

class Worker
{

public:

  // Constructor
  Worker();
  ~Worker();

  // Return a handle to the currently loaded geometry
  Geometry *GetGeometry();
  SynradGeometry* GetSynradGeometry();

  void AddMaterial(std::string *fileName);
  void LoadGeometry(char *fileName, bool insert=false, bool newStr=false);// Loads or inserts a geometry (throws Error)
  //void InsertGeometry(bool newStr,char *fileName); // Inserts a new geometry (throws Error)
  void LoadTexturesSYN(FileReader* f,int version);  // Load a textures(throws Error)
  void RebuildTextures();

  std::vector<std::vector<double>> ImportCSV(FileReader *file);
  
  // Save a geometry (throws Error)
  void SaveGeometry(char *fileName,GLProgress *prg,bool askConfirm=true,bool saveSelected=false,bool autoSave=false,bool crashSave=false);

  // Export textures (throws Error)
  bool IsDpInitialized();
  void ExportTextures(char *fileName,int grouping,int mode,bool askConfirm=true,bool saveSelected=false);
  //void ExportRegionPoints(char *fileName,GLProgress *prg,int regionId,int exportFrequency,bool doFullScan);
  void ExportDesorption(char *fileName,bool selectedOnly,int mode,double eta0,double alpha,Distribution2D *distr);

  // Return/Set the current filename
  char *GetFileName();
  char *GetShortFileName();
  char *GetShortFileName(char* longFileName);
  void  SetFileName(char *fileName);

  void SetProcNumber(size_t n);// Set number of processes [1..32] (throws Error)
  size_t GetProcNumber();  // Get number of processes
  void SetMaxDesorption(llong max);// Set the number of maximum desorption
  DWORD GetPID(size_t prIdx);// Get PID
  void ResetStatsAndHits(float appTime);
  void Reload();    // Reload simulation (throws Error)
  void RealReload();
  void ChangeSimuParams();
  void StartStop(float appTime,int mode);    // Switch running/stopped
  void Stop_Public();// Switch running/stopped
  //void Exit(); // Free all allocated resource
  void KillAll();// Kill all sub processes
  void Update(float appTime);// Get hit counts for sub process
  void SendHits();// Send total and facet hit counts to subprocesses
  void SetLeakCache(LEAK *buffer,size_t *nb,Dataport *dpHit);// Set Leak
  void SetHitCache(HIT *buffer,size_t *nb, Dataport *dpHit);  // Set HHit
  void GetProcStatus(int *states,std::vector<std::string>& statusStrings);// Get process status
  BYTE *GetHits(); // Access to dataport (HIT)
  void ReleaseHits();
  void ClearRegions();
  void RemoveRegion(int index);
  void AddRegion(const char *fileName,int position=-1); //load region (position==-1: add as new region)
  void RecalcRegion(int regionId);
  void SaveRegion(char *fileName,int position,bool overwrite=false);

  // Global simulation parameters
  size_t  nbAbsorption;      // Total number of molecules absorbed (64 bit integer)
  size_t  nbDesorption;      // Total number of molecules generated (64 bit integer)
  size_t  nbHit;             // Total number of hit (64 bit integer)
  size_t  nbLeakTotal;            // Total number of leak
  double totalFlux;         // Total desorbed Flux
  double totalPower;        // Total desorbed power
  size_t  desorptionLimit;     // Number of desoprtion before halting
  double distTraveledTotal; // Total distance traveled by particles (for mean free path calc.)
  bool   running;           // Started/Stopped state
  float  startTime;         // Start time
  float  stopTime;          // Stop time
  float  simuTime;          // Total simulation time
  int    mode;              // Simulation mode
  size_t    nbTrajPoints;       // number of all points in trajectory
  double no_scans;           // = nbDesorption/nbTrajPoints. Stored separately for saving/loading
  int    generation_mode;   //fluxwise or powerwise
  bool   lowFluxMode;
  double lowFluxCutoff;
  bool   newReflectionModel;
  std::vector<Region_full> regions;
  std::vector<Material> materials;
  std::vector<std::vector<double>> psi_distro; //psi-energy map for full (par+ort) polarization
  std::vector<std::vector<double>> parallel_polarization; //ratio of parallel polarization for a given E/E_crit ratio and psi vertical angle
  std::vector<std::vector<std::vector<double>>> chi_distros; //3 psi-chi    maps for full/parallel/orthogonal polarizations
  
  char fullFileName[512]; // Current loaded file

  bool needsReload;
  bool abortRequested;

  bool calcAC; //Not used in Synrad, kept for ResetStatsAndHits function shared with Molflow

  // Caches
  HIT  hitCache[HITCACHESIZE];
  LEAK leakCache[LEAKCACHESIZE];
  size_t    hitCacheSize;            // Total number of hhit
  size_t leakCacheSize;

private:

  // Process management
  size_t    nbProcess;
  DWORD  pID[MAX_PROCESS];
  DWORD  pid;
  bool   allDone;

  // Geometry handle
  SynradGeometry *geom;

  // Dataport handles and names
  Dataport *dpControl;
  Dataport *dpHit,*dpMat;
  char      ctrlDpName[32];
  char      loadDpName[32];
  char      hitsDpName[32];
  char      materialsDpName[32];

  // Methods
  bool ExecuteAndWait(int command, int waitState, size_t param=0);
  bool Wait(int waitState, LoadStatus *statusWindow);
  void ResetWorkerStats();
  void ClearHits(bool noReload);
  char *GetErrorDetails();
  void ThrowSubProcError(std::string message);
  void ThrowSubProcError(char *message=NULL);
  void Start();
  void Stop();
  //void OneStep();
  void InnerStop(float appTime);

};

#endif /* _WORKERH_ */
