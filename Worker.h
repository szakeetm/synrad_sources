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

#include "Geometry.h"
#include "Region_full.h"
#include "Distributions.h" //for materials

extern float m_fTime;
class Worker
{

public:

  // Constructor
  Worker();
  ~Worker();

  // Return a handle to the currently loaded geometry
  Geometry *GetGeometry();
  void AddMaterial(string *fileName);
  void LoadGeometry(char *fileName);// Load a geometry (throws Error)
  void InsertGeometry(BOOL newStr,char *fileName); // Inserts a new geometry (throws Error)
  void loadTextures(char *fileName,int version);  // Load a textures(throws Error)
    void RebuildTextures();
  
  // Save a geometry (throws Error)
  void SaveGeometry(char *fileName,GLProgress *prg,BOOL askConfirm=TRUE,BOOL saveSelected=FALSE,BOOL autoSave=FALSE,BOOL crashSave=FALSE);

  // Export textures (throws Error)
  BOOL IsDpInitialized();
  void ExportTextures(char *fileName,int grouping,int mode,BOOL askConfirm=TRUE,BOOL saveSelected=FALSE);
  void ExportRegionPoints(char *fileName,GLProgress *prg,int regionId,int exportFrequency,BOOL doFullScan);
  void ExportDesorption(char *fileName,bool selectedOnly,int mode,double eta0,double alpha,Distribution2D *distr);

  // Save a geometry using the current file name (throws Error)
  void SaveGeometry(GLProgress *prg);

  // Return/Set the current filename
  char *GetFileName();
  char *GetShortFileName();
  void  SetFileName(char *fileName);
  void SetProcNumber(int n);// Set number of processes [1..16] (throws Error)
  int GetProcNumber();  // Get number of processes
  void SetMaxDesorption(llong max);// Set the number of maximum desorption
  DWORD GetPID(int prIdx);// Get PID
  void Reset(float appTime);// Reset simulation
  void Reload();    // Reload simulation (throws Error)
  void RealReload();
  void StartStop(float appTime,int mode);    // Switch running/stopped
  void Stop_Public();// Switch running/stopped
  void Exit(); // Free all allocated resource
  void KillAll();// Kill all sub processes
  void Update(float appTime);// Get hit counts for sub process
  void SendHits();// Send total and facet hit counts to subprocesses
  void GetLeak(LEAK *buffer,int *nb);  // Get Leak
  void SetLeak(LEAK *buffer,int *nb,SHGHITS *gHits);// Set Leak
  void GetHHit(HIT *buffer,int *nb); // Get HHit
  void SetHHit(HIT *buffer,int *nb,SHGHITS *gHits);  // Set HHit
  void GetProcStatus(int *states,char **status);// Get process status
  BYTE *GetHits(); // Access to dataport (HIT)
  void  ReleaseHits();
  void ClearRegions();
  void AddRegion(char *fileName,int position=-1); //load region (position==-1: add as new region)
  void RecalcRegion(int regionId);
  void SaveRegion(char *fileName,int position,BOOL overwrite=FALSE);


  // Global simulation parameters
  llong  nbAbsorption;      // Total number of molecules absorbed (64 bit integer)
  llong  nbDesorption;      // Total number of molecules generated (64 bit integer)
  llong  nbHit;             // Total number of hit (64 bit integer)
  double totalFlux;         // Total desorbed Flux
  double totalPower;        // Total desorbed power
  llong  maxDesorption;     // Number of desoprtion before halting
  double distTraveledTotal; // Total distance traveled by particles (for mean free path calc.)
  llong  nbLeakTotal;            // Total number of leak
  int    nbLastLeaks;
  int    nbHHit;            // Total number of hhit
  BOOL   running;           // Started/Stopped state
  float  startTime;         // Start time
  float  stopTime;          // Stop time
  float  simuTime;          // Total simulation time
  int    mode;              // Simulation mode
  int    nbTrajPoints;       // number of all points in trajectory
  double no_scans;           // = nbDesorption/nbTrajPoints. Stored separately for saving/loading
  int    generation_mode;   //fluxwise or powerwise
  BOOL   lowFluxMode;
  double lowFluxCutoff;
  std::vector<Region_full> regions;
  std::vector<Material> materials;
  char fullFileName[512]; // Current loaded file

  

private:

  // Process management
  int    nbProcess;
  DWORD  pID[MAX_PROCESS];
  DWORD  pid;
  BOOL   allDone;

  // Geometry handle
  Geometry *geom;

  // Dataport handles and names
  Dataport *dpControl;
  Dataport *dpHit,*dpMat;
  char      ctrlDpName[32];
  char      loadDpName[32];
  char      hitsDpName[32];
  char      materialsDpName[32];

  // Caches
  HIT  hhitCache[NBHHIT];
  LEAK leakCache[NBHHIT];

  // Methods
  BOOL ExecuteAndWait(int command,int waitState,int param=0,GLProgress *prg=NULL);
  BOOL Wait(int waitState,int timeout,GLProgress *prg=NULL);
  void ResetWorkerStats();
  void ClearHits();
  char *GetErrorDetails();
  void ThrowSubProcError(char *message=NULL);
  void Start();
  void Stop();
  void OneStep();
  void InnerStop(float appTime);

};

#endif /* _WORKERH_ */
