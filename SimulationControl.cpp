/*
File:        SimulationControl.c
Description: Simulation control routines
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

#pragma once

#ifdef WIN
#include <windows.h> // For GetTickCount()
#include <Process.h> // For _getpid()
#else
#include <time.h>
#include <sys/time.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "Simulation.h"
#include "Random.h"
#include "SynradTypes.h" //Histogram
//#include "Tools.h"

extern void SetErrorSub(const char *message);

// -------------------------------------------------------
// Global handles
// -------------------------------------------------------

FACET     **THits;
SIMULATION *sHandle;

// -------------------------------------------------------
// Timing stuff
// -------------------------------------------------------

#ifdef WIN
BOOL usePerfCounter;         // Performance counter usage
LARGE_INTEGER perfTickStart; // First tick
double perfTicksPerSec;      // Performance counter (number of tick per second)
#endif
DWORD tickStart;

// -------------------------------------------------------

void InitSimulation() {

	// Global handle allocation
	sHandle = (SIMULATION *)malloc(sizeof(SIMULATION));
	memset(sHandle,0,sizeof(SIMULATION));
	THits = (FACET **)malloc(MAX_THIT*sizeof(FACET *)); // Transparent hit cache

#ifdef WIN
	{
		LARGE_INTEGER qwTicksPerSec;
		usePerfCounter = QueryPerformanceFrequency( &qwTicksPerSec );
		if( usePerfCounter ) {
			QueryPerformanceCounter( &perfTickStart );
			perfTicksPerSec = (double)qwTicksPerSec.QuadPart;
		}
		tickStart = GetTickCount();
	}
#else
	tickStart = (DWORD)time(NULL);
#endif

}

// -------------------------------------------------------

void ClearSimulation() {

	int i,j;

	// Free old stuff
	/*for (i=0;i<(int)sHandle->regions.size();i++)
	delete &(sHandle->regions[i]);*/
	sHandle->regions.clear();
	//sHandle->regions=std::vector<Region_mathonly>();
	SAFE_FREE(sHandle->vertices3);
	for(j=0;j<sHandle->nbSuper;j++) {
		for(i=0;i<sHandle->str[j].nbFacet;i++) {
			FACET *f = sHandle->str[j].facets[i];
			if( f ) {
				SAFE_FREE(f->indices);
				SAFE_FREE(f->vertices2);
				SAFE_FREE(f->inc);
				SAFE_FREE(f->largeEnough);
				SAFE_FREE(f->hits_MC);
				SAFE_FREE(f->hits_flux);
				SAFE_FREE(f->hits_power);
				//SAFE_FREE(f->area);
				SAFE_FREE(f->profile_hits);
				SAFE_FREE(f->profile_flux);
				SAFE_FREE(f->profile_power);
				delete f->spectrum_fluxwise;
				delete f->spectrum_powerwise;
				SAFE_FREE(f->direction);
				//SAFE_FREE(f->fullElem);
			}
			SAFE_FREE(f);
		}
		SAFE_FREE(sHandle->str[j].facets);
		if( sHandle->str[j].aabbTree ) {
			DestroyAABB(sHandle->str[j].aabbTree->left);
			DestroyAABB(sHandle->str[j].aabbTree->right);
			free(sHandle->str[j].aabbTree);
			sHandle->str[j].aabbTree=NULL;
		}
	}
	memset(sHandle,0,sizeof(SIMULATION));

}
// -------------------------------------------------------

DWORD RevertBit(DWORD dw) {
	DWORD dwIn   = dw;
	DWORD dwOut  = 0;
	DWORD dWMask = 1;
	int i;
	for(i=0;i<32;i++) {
		if( dwIn & 0x80000000UL ) dwOut |= dWMask;
		dwIn = dwIn << 1;
		dWMask = dWMask << 1;
	}
	return dwOut;
}

DWORD GetSeed() {

	/*#ifdef WIN
	DWORD r;
	_asm {
	rdtsc
	mov r,eax
	}
	return RevertBit(r ^ (DWORD)(_getpid()*65519));
	#else*/
	return (DWORD)((int)(GetTick()*1000.0)*_getpid());
	//#endif

}



double Norme(Vector3d *v) {
	return sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
}



BOOL LoadSimulation(Dataport *loader) {

	int i,j,idx;
	BYTE *buffer,*bufferAfterMaterials;
	BYTE *areaBuff;
	BYTE *bufferStart;
	Vector3d *shVert;
	double t1,t0;
	DWORD seed;
	char err[128];

	t0 = GetTick();

	sHandle->loadOK = FALSE;
	SetState(PROCESS_STARTING, "Clearing previous simulation");
	try {
		ClearSimulation();
	}
	catch (...) {
		SetErrorSub("Error clearing geometry");
		return FALSE;
	}
	/* //Mutex not necessary: by the time the COMMAND_LOAD is issued the interface releases the handle, concurrent reading is safe and it's only destroyed by the interface when all processes are ready loading
	//Result: faster, parallel loading
	// Connect the dataport
	SetState(PROCESS_STARTING, "Waiting for 'loader' dataport access...");
	if( !AccessDataportTimed(loader,40000) ) {
		SetErrorSub("Failed to connect to DP");
		return FALSE;
	}
	*/

	bufferStart = (BYTE *)loader->buff;
	buffer = bufferStart;

	// Load new geom from the dataport

	SHGEOM* shGeom = (SHGEOM *)buffer;
	if(shGeom->nbSuper>MAX_STRUCT) {
		ReleaseDataport(loader);
		SetErrorSub("Too many structures");
		return FALSE;
	}
	if(shGeom->nbSuper<=0) {
		ReleaseDataport(loader);
		SetErrorSub("No structures");
		return FALSE;
	}

	sHandle->newReflectionModel = shGeom->newReflectionModel;
	sHandle->nbVertex = shGeom->nbVertex;
	sHandle->nbSuper = shGeom->nbSuper;
	sHandle->nbRegion = shGeom->nbRegion;
	sHandle->totalFacet = shGeom->nbFacet;
	
	buffer+=sizeof(SHGEOM);

	SHMODE* shMode = (SHMODE *)buffer;
	sHandle->mode.generation_mode = shMode->generation_mode;
	sHandle->mode.lowFluxCutoff = shMode->lowFluxCutoff;
	sHandle->mode.lowFluxMode = shMode->lowFluxMode;
	buffer += sizeof(SHMODE);

	sHandle->nbTrajPoints=0;
	
	//Regions
	try{
	for (size_t r = 0; r<sHandle->nbRegion; r++) {
		RegionParams *regparam = (RegionParams *) buffer;
		Region_mathonly newreg;
		newreg.params = *regparam;//memcpy(&newreg.params, regparam, sizeof(RegionParams));
		buffer += sizeof(RegionParams);
		sHandle->regions.push_back(newreg);
	}
	//copy trajectory points	
	for (size_t r = 0; r<sHandle->nbRegion; r++) {
		/*sHandle->regions[r].Points.reserve(sHandle->regions[r].params.nbPointsToCopy);
		for (int k=0;k<sHandle->regions[r].params.nbPointsToCopy;k++) {
			sHandle->regions[r].Points.push_back(*(Trajectory_Point*)(buffer));
			buffer+=sizeof(Trajectory_Point);
		}*/

		sHandle->regions[r].Points.resize(sHandle->regions[r].params.nbPointsToCopy);
		memcpy(&sHandle->regions[r].Points[0], buffer, sizeof(Trajectory_Point)*sHandle->regions[r].params.nbPointsToCopy);
		buffer += sizeof(Trajectory_Point)*sHandle->regions[r].params.nbPointsToCopy;
		sHandle->nbTrajPoints+=sHandle->regions[r].params.nbPointsToCopy;
	}

	//copy distribution points
	//sHandle->nbDistrPoints_MAG=0;
	sHandle->nbDistrPoints_BXY=0;
	for (size_t r = 0; r<sHandle->nbRegion; r++) {

		sHandle->regions[r].Bx_distr.Resize(sHandle->regions[r].params.nbDistr_MAG.x);
		for (size_t j = 0; j<sHandle->regions[r].params.nbDistr_MAG.x; j++) {
			sHandle->regions[r].Bx_distr.valuesX[j]=READBUFFER(double);
			sHandle->regions[r].Bx_distr.valuesY[j]=READBUFFER(double);
		}

		sHandle->regions[r].By_distr.Resize(sHandle->regions[r].params.nbDistr_MAG.y);
		for (size_t j = 0; j<sHandle->regions[r].params.nbDistr_MAG.y; j++) {
			sHandle->regions[r].By_distr.valuesX[j] = READBUFFER(double);
			sHandle->regions[r].By_distr.valuesY[j] = READBUFFER(double);
		}

		sHandle->regions[r].Bz_distr.Resize(sHandle->regions[r].params.nbDistr_MAG.z);
		for (size_t j = 0; j<sHandle->regions[r].params.nbDistr_MAG.z; j++) {
			sHandle->regions[r].Bz_distr.valuesX[j] = READBUFFER(double);
			sHandle->regions[r].Bz_distr.valuesY[j] = READBUFFER(double);
		}

		sHandle->regions[r].betaFunctions.Clear();

		for (size_t j = 0; j<sHandle->regions[r].params.nbDistr_BXY; j++) {
			double x = READBUFFER(double);
			std::vector<double> betaValues;
			betaValues.resize(6);
			memcpy(betaValues.data(), buffer, 6 * sizeof(double)); //Vector data is contigious. The 6 values are: BetaX, BetaY, EtaX, EtaX', AlphaX, AlphaY
			buffer += 6 * sizeof(double);
		}

		sHandle->nbDistrPoints_BXY+=sHandle->regions[r].params.nbDistr_BXY;
	}
	}
	catch (...) {
		SetErrorSub("Error loading regions");
		return FALSE;
	}

	//Load materials
	try{
	sHandle->nbMaterials = READBUFFER(size_t); //copying number of materials
	for (size_t i = 0; i<sHandle->nbMaterials; i++) { //going through all materials
		Material newMaterial;
		newMaterial.hasBackscattering = READBUFFER(BOOL);
		size_t angleValsSize = READBUFFER(size_t); //copying number of angles (columns)
		size_t energyValsSize = READBUFFER(size_t); //copying number of energies (rows)
		newMaterial.angleVals.reserve(angleValsSize);
		for (size_t j = 0; j<angleValsSize; j++) {
			newMaterial.angleVals.push_back(*((double*)buffer)); //copying angles (header)
			buffer+=sizeof(double);
		}
		newMaterial.energyVals.reserve(energyValsSize);
		for (size_t j = 0; j<energyValsSize; j++) {
			newMaterial.energyVals.push_back(*((double*)buffer)); //copying energies (column1)
			buffer+=sizeof(double);
		}
		//copying reflectivity probabilities (cells)
		for (size_t j = 0; j<newMaterial.energyVals.size(); j++) {
			std::vector<std::vector<double>> currentEnergy;
			for (size_t k = 0; k<newMaterial.angleVals.size(); k++) {
				std::vector<double> reflProbability;
				reflProbability.push_back(*((double*)buffer)); //forward scattering
				buffer+=sizeof(double);
				if (newMaterial.hasBackscattering) {
					reflProbability.push_back(*((double*)buffer)); //diffuse scattering
					buffer += sizeof(double);
					reflProbability.push_back(*((double*)buffer)); //backscattering
					buffer += sizeof(double);
					reflProbability.push_back(*((double*)buffer)); //transparent pass
					buffer += sizeof(double);
				}
				currentEnergy.push_back(reflProbability);
			}
			newMaterial.reflVals.push_back(currentEnergy);
		}
		sHandle->materials.push_back(newMaterial);
	}
	}
	catch (...) {
		SetErrorSub("Error loading materials");
		return FALSE;
	}

	try {
		//Load psi_distr
		size_t psi_size = READBUFFER(size_t);
		sHandle->psi_distr.reserve(psi_size);
		for (size_t j = 0; j < psi_size; j++) {
			std::vector<double> row;
			size_t row_size = READBUFFER(size_t);
			sHandle->psi_distr.reserve(row_size);
			for (size_t k = 0; k < row_size; k++) {
				row.push_back(*((double*)buffer)); //cum. distr
				buffer += sizeof(double);
			}
			sHandle->psi_distr.push_back(row);
		}

		//Load chi_distr
		size_t chi_size = READBUFFER(size_t);
		sHandle->chi_distr.reserve(chi_size);
		for (size_t j = 0; j < chi_size; j++) {
			std::vector<double> row;
			int row_size = READBUFFER(size_t);
			sHandle->chi_distr.reserve(row_size);
			for (size_t k = 0; k < row_size; k++) {
				row.push_back(*((double*)buffer)); //cum. distr
				buffer += sizeof(double);
			}
			sHandle->chi_distr.push_back(row);
		}
	} 
	catch (...) {
		SetErrorSub("Error loading distributions");
		return FALSE;
	}

	bufferAfterMaterials=buffer;

	// Prepare super structure
	buffer +=  sizeof(Vector3d)*sHandle->nbVertex;
	for(i=0;i<sHandle->totalFacet;i++) {
		SHFACET *shFacet = (SHFACET *)buffer;
		sHandle->str[shFacet->superIdx].nbFacet++;
		buffer+=sizeof(SHFACET) + shFacet->nbIndex*(sizeof(int) + sizeof(Vector2d));
	}
	for(i=0;i<sHandle->nbSuper;i++) {
		int nbF = sHandle->str[i].nbFacet;
		if( nbF==0 ) {
		} else {
			sHandle->str[i].facets=(FACET **)malloc(nbF*sizeof(FACET *));
			memset(sHandle->str[i].facets,0,nbF*sizeof(FACET *));
		}
		sHandle->str[i].nbFacet = 0;
	}
	areaBuff = buffer;
	//buffer = (BYTE *)loader->buff;

	// Name
	memcpy(sHandle->name,shGeom->name,64);

	// Vertex
	sHandle->vertices3 = (Vector3d *)malloc(sHandle->nbVertex*sizeof(Vector3d));
	buffer=bufferAfterMaterials;
	shVert = (Vector3d *)(buffer);
	memcpy(sHandle->vertices3,shVert,sHandle->nbVertex*sizeof(Vector3d));
	buffer+=sizeof(Vector3d)*sHandle->nbVertex;

	// Facets
	for(i=0;i<sHandle->totalFacet;i++) {

		SHFACET *shFacet = (SHFACET *)buffer;
		FACET *f = (FACET *)malloc(sizeof(FACET));
		if (!f) {
			SetErrorSub("Not enough memory to load facets");
			return FALSE;
		}
		memset(f,0,sizeof(FACET));
		memcpy(&(f->sh),shFacet,sizeof(SHFACET));

		sHandle->hasVolatile |= f->sh.isVolatile;
		//sHandle->hasDirection |= f->sh.countDirection;

		idx = f->sh.superIdx;
		sHandle->str[idx].facets[sHandle->str[idx].nbFacet] = f;
		sHandle->str[idx].facets[sHandle->str[idx].nbFacet]->globalId = i;
		sHandle->str[idx].nbFacet++;


		if( f->sh.superDest || f->sh.isVolatile ) {
			// Link or volatile facet, overides facet settings
			// Must be full opaque and 0 sticking
			// (see SimulationMC.c::PerformBounce)
			//f->sh.isOpaque = TRUE;
			f->sh.opacity = 1.0;
			f->sh.sticking = 0.0;
			if( (f->sh.superDest-1) >= sHandle->nbSuper || f->sh.superDest<0 ) {
				// Geometry error
				ClearSimulation();
				ReleaseDataport(loader);
				sprintf(err,"Invalid structure (wrong link on F#%d)",i+1);
				SetErrorSub(err);
				return FALSE;
			}
		}

		// Reset counter in local memory
		memset(&(f->counter),0,sizeof(SHHITS));
		f->indices = (int *)malloc(f->sh.nbIndex*sizeof(int));
		buffer+=sizeof(SHFACET);
		memcpy(f->indices,buffer,f->sh.nbIndex*sizeof(int));
		buffer+=f->sh.nbIndex*sizeof(int);
		f->vertices2 = (Vector2d *)malloc(f->sh.nbIndex * sizeof(Vector2d));
		if (!f->vertices2) {
			SetErrorSub("Not enough memory to load vertices");
			return FALSE;
		}
		memcpy(f->vertices2,buffer,f->sh.nbIndex * sizeof(Vector2d));
		buffer+=f->sh.nbIndex*sizeof(Vector2d);

		//Textures
		if(f->sh.isTextured) {
			int nbE = f->sh.texWidth*f->sh.texHeight;
			f->textureSize = nbE*(2*sizeof(double)+sizeof(llong));
			f->hits_MC = (llong *)malloc(nbE*sizeof(llong));
			memset(f->hits_MC,0,nbE*sizeof(llong));
			f->hits_flux = (double *)malloc(nbE*sizeof(double));
			memset(f->hits_flux,0,nbE*sizeof(double));
			f->hits_power = (double *)malloc(nbE*sizeof(double));
			memset(f->hits_power,0,nbE*sizeof(double));

			f->inc = (double*)malloc(f->textureSize);
			f->largeEnough = (BOOL *)malloc(sizeof(BOOL)*nbE);
			//f->fullElem = (char *)malloc(nbE);
			if (!(f->inc && f->largeEnough/* && f->fullElem*/)) {
				SetErrorSub("Not enough memory to load");
				return FALSE;
			}
			f->fullSizeInc=(float)1E30;
			for(j=0;j<nbE;j++) {
				double incVal = ((double *)areaBuff)[j];
				/*if( incVal<0 ) {
					f->fullElem[j] = 1;
					f->inc[j] = -incVal;
				} else {
					f->fullElem[j] = 0;*/
					f->inc[j] = incVal;
				/*}*/
				if ((f->inc[j]>0.0)&&(f->inc[j]<f->fullSizeInc)) f->fullSizeInc = f-> inc[j];
			}
			for(j=0;j<nbE;j++) { //second pass, filter out very small cells
				f->largeEnough[j]=(f->inc[j]<((5.0)*f->fullSizeInc));
			}
			sHandle->textTotalSize += f->textureSize;
			areaBuff += nbE*sizeof(double);
			f->iw = 1.0 / (double)f->sh.texWidthD;
			f->ih = 1.0 / (double)f->sh.texHeightD;
			f->rw = Norme(&(f->sh.U)) * f->iw;
			f->rh = Norme(&(f->sh.V)) * f->ih;
		}
		if(f->sh.isProfile) {
			f->profileSize = PROFILE_SIZE*(sizeof(llong)+2*sizeof(double));

			f->profile_hits = (llong *)malloc(PROFILE_SIZE*sizeof(llong));
			memset(f->profile_hits,0,PROFILE_SIZE*sizeof(llong));

			f->profile_flux = (double *)malloc(PROFILE_SIZE*sizeof(double));
			memset(f->profile_flux,0,PROFILE_SIZE*sizeof(double));

			f->profile_power = (double *)malloc(PROFILE_SIZE*sizeof(double));
			memset(f->profile_power,0,PROFILE_SIZE*sizeof(double));

			sHandle->profTotalSize += f->profileSize; 
		}
		if(f->sh.countDirection) {
			f->directionSize = f->sh.texWidth*f->sh.texHeight*sizeof(VHIT);
			f->direction = (VHIT *)malloc(f->directionSize);
			memset(f->direction,0,f->directionSize);
			sHandle->dirTotalSize += f->directionSize;
		}
		if(f->sh.hasSpectrum) {
			f->spectrumSize = 2*sizeof(double)*SPECTRUM_SIZE;
			double min_energy,max_energy;
			if (sHandle->nbRegion>0) {
				min_energy=sHandle->regions[0].params.energy_low_eV;
				max_energy=sHandle->regions[0].params.energy_hi_eV;
			} else {
				min_energy=10.0;
				max_energy=1000000.0;
			}
			f->spectrum_fluxwise = new Histogram(min_energy,max_energy,SPECTRUM_SIZE,true);
			f->spectrum_powerwise = new Histogram(min_energy,max_energy,SPECTRUM_SIZE,true);
			sHandle->spectrumTotalSize += f->spectrumSize;
		}

	}

	//ReleaseDataport(loader); //Commented out as AccessDataport removed

	// Build all AABBTrees
	for(i=0;i<sHandle->nbSuper;i++)
		sHandle->str[i].aabbTree = BuildAABBTree(sHandle->str[i].facets,sHandle->str[i].nbFacet,0);

	// Initialise simulation
	ComputeSourceArea();
	seed = GetSeed();
	rseed( seed );

	//--- GSL random init ---
	gsl_rng_env_setup();                          // Read variable environnement
	const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
	sHandle->gen = gsl_rng_alloc(type);          // Rand generator allocation

	sHandle->loadOK = TRUE;
	t1 = GetTick();
	printf("  Load %s successful\n",sHandle->name);
	printf("  Geometry: %d vertex %d facets\n",sHandle->nbVertex,sHandle->totalFacet);
	printf("  Region: %zd regions\n",sHandle->nbRegion);
	printf("  Trajectory points: %zd points\n", sHandle->nbTrajPoints);
	printf("  Geom size: %d bytes\n",(int)(buffer-bufferStart));
	printf("  Number of stucture: %d\n",sHandle->nbSuper);
	printf("  Global Hit: %d bytes\n",sizeof(SHGHITS));
	printf("  Facet Hit : %d bytes\n",sHandle->totalFacet*(int)sizeof(SHHITS));
	printf("  Texture   : %zd bytes\n",sHandle->textTotalSize);
	printf("  Profile   : %zd bytes\n",sHandle->profTotalSize);
	printf("  Direction : %zd bytes\n",sHandle->dirTotalSize);
	printf("  Spectrum  : %zd bytes\n",sHandle->spectrumTotalSize);
	printf("  Total     : %zd bytes\n",GetHitsSize());
	printf("  Seed: %u\n",seed);
	printf("  Loading time: %.3f ms\n",(t1-t0)*1000.0);
	return TRUE;

}

BOOL UpdateSimuParams(Dataport *loader) {
	// Connect the dataport
	if (!AccessDataportTimed(loader, 2000)) {
		SetErrorSub("Failed to connect to DP");
		return FALSE;
	}
	BYTE* buffer = (BYTE *)loader->buff;
	
	SHMODE* shMode = (SHMODE *)buffer;
	sHandle->mode.generation_mode = shMode->generation_mode;
	sHandle->mode.lowFluxCutoff = shMode->lowFluxCutoff;
	sHandle->mode.lowFluxMode = shMode->lowFluxMode;
	buffer += sizeof(SHMODE);

	for (size_t i = 0; i < sHandle->nbRegion; i++) {
		sHandle->regions[i].params.showPhotons = READBUFFER(BOOL);
	}

	ReleaseDataport(loader);

	//Reset hit cache
	sHandle->hitCacheSize = 0;
	//memset(sHandle->hitCache, 0, sizeof(HIT)*HITCACHESIZE);
	
	sHandle->leakCacheSize = 0;
	//memset(sHandle->leakCache, 0, sizeof(LEAK)*LEAKCACHESIZE); //No need to reset, will gradually overwrite

	return TRUE;
}

void UpdateHits(Dataport *dpHit,int prIdx,DWORD timeout) {

		UpdateMCHits(dpHit,prIdx,timeout);

}

size_t GetHitsSize() {
	return sHandle->textTotalSize + sHandle->profTotalSize + sHandle->dirTotalSize +
		sHandle->spectrumTotalSize + sHandle->totalFacet*sizeof(SHHITS) + sizeof(SHGHITS);
}

void ResetTmpCounters() {
	SetState(NULL, "Resetting local cache...", FALSE, TRUE);

	memset(&sHandle->tmpCount, 0, sizeof(SHHITS));

	sHandle->distTraveledSinceUpdate = 0.0;
	sHandle->nbLeakSinceUpdate = 0;
	sHandle->hitCacheSize = 0;
	sHandle->leakCacheSize = 0;

	for(int j=0;j<sHandle->nbSuper;j++) {
		for(int i=0;i<sHandle->str[j].nbFacet;i++) {
			FACET *f = sHandle->str[j].facets[i];
			f->counter.fluxAbs=0.0;
			f->counter.powerAbs=0.0;
			f->counter.nbHit=0;
			f->counter.nbAbsorbed=0;
			f->hitted = FALSE;
			int textureElemNb=f->sh.texHeight*f->sh.texWidth;
			if( f->hits_MC ) memset(f->hits_MC,0,textureElemNb*sizeof(llong));
			if( f->hits_flux ) memset(f->hits_flux,0,textureElemNb*sizeof(double));
			if( f->hits_power ) memset(f->hits_power,0,textureElemNb*sizeof(double));

			if( f->profile_hits ) memset(f->profile_hits,0,PROFILE_SIZE*sizeof(llong));
			if( f->profile_flux ) memset(f->profile_flux,0,PROFILE_SIZE*sizeof(double));
			if( f->profile_power ) memset(f->profile_power,0,PROFILE_SIZE*sizeof(double));
			if( f->spectrum_fluxwise ) f->spectrum_fluxwise->ResetCounts();
			if( f->spectrum_powerwise ) f->spectrum_powerwise->ResetCounts();

			if( f->direction ) memset(f->direction,0,f->directionSize);
		}
	}

}

void ResetSimulation() {
	sHandle->lastHit = NULL;
	sHandle->totalDesorbed = 0;
	ResetTmpCounters();

}

BOOL StartSimulation() {
	if (sHandle->regions.size()==0) {
		SetErrorSub("No regions");
		return FALSE;
	}
	
	//Check for invalid material (ID==9) passed from the interface
	//It is valid to load invalid materials (to extract textures, etc.) but not to launch the simulation
	for (int s = 0; s<sHandle->nbSuper; s++) {
		for (int f = 0; f < sHandle->str[s].nbFacet; f++) {
			if (sHandle->str[s].facets[f]->sh.reflectType == 9) {
				char tmp[32];
				sprintf(tmp, "Invalid material on Facet %d.", f + 1);
				SetErrorSub(tmp);
				return FALSE;
			}
		}
	}

	if(!sHandle->lastHit) StartFromSource();
	return true;
}

// -------------------------------------------------------

void RecordHit(const int &type,const double &dF,const double &dP) {
	if (sHandle->regions[sHandle->sourceRegionId].params.showPhotons) {
		if (sHandle->hitCacheSize < HITCACHESIZE) {
			sHandle->hitCache[sHandle->hitCacheSize].pos = sHandle->pPos;
			sHandle->hitCache[sHandle->hitCacheSize].type = type;
			sHandle->hitCache[sHandle->hitCacheSize].dF = sHandle->dF;
			sHandle->hitCache[sHandle->hitCacheSize].dP = sHandle->dP;
			sHandle->hitCacheSize++;
		}
	}
}

void RecordLeakPos() {
	// Source region check performed when calling this routine 
	// Record leak for debugging
	RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
	RecordHit(HIT_LAST, sHandle->dF, sHandle->dP);
	if (sHandle->leakCacheSize < LEAKCACHESIZE) {
		sHandle->leakCache[sHandle->leakCacheSize].pos = sHandle->pPos;
		sHandle->leakCache[sHandle->leakCacheSize].dir = sHandle->pDir;
		sHandle->leakCacheSize++;
	}
}

// -------------------------------------------------------

BOOL SimulationRun() {

	// 1s step
	double t0,t1;
	int    nbStep=1;
	BOOL   goOn;

	if( sHandle->stepPerSec==0.0 ) {
		
			nbStep=250;
			
	}

	if( sHandle->stepPerSec!=0.0 )
		nbStep = (int)(sHandle->stepPerSec+0.5);
	if(nbStep<1) nbStep = 1;
	t0 = GetTick();

		goOn = SimulationMCStep(nbStep);

	
	t1 = GetTick();
	sHandle->stepPerSec = (double)(nbStep) / (t1-t0);
#ifdef _DEBUG
	printf("Running: stepPerSec = %f\n",sHandle->stepPerSec);
#endif

	return !goOn;

}

// -------------------------------------------------------

double GetTick() {

	// Number of sec since the application startup

#ifdef WIN

	if( usePerfCounter ) {

		LARGE_INTEGER t,dt;
		QueryPerformanceCounter( &t );
		dt.QuadPart = t.QuadPart - perfTickStart.QuadPart;
		return (double)(dt.QuadPart)/perfTicksPerSec;

	} else {

		return (double)((GetTickCount() - tickStart)/1000.0);

	}

#else

	if(tickStart < 0 )
		tickStart = time(NULL);

	struct timeval tv;
	gettimeofday(&tv,NULL);
	return ( (double)(tv.tv_sec-tickStart)*1000.0 + (double)tv.tv_usec/1000.0 );

#endif

}
