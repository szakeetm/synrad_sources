/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/
#pragma once

#ifdef WIN
#define NOMINMAX
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
#include <sstream>
#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include "Random.h"
#include "SynradTypes.h" //Histogram
//#include "Tools.h"
#include <cereal/types/utility.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/json.hpp>

extern void SetErrorSub(const char *message);

// Global handles

extern Simulation *sHandle;

// Timing stuff

#ifdef WIN
bool usePerfCounter;         // Performance counter usage
LARGE_INTEGER perfTickStart; // First tick
double perfTicksPerSec;      // Performance counter (number of tick per second)
#endif
DWORD tickStart;

void InitSimulation() {

	// Global handle allocation
	sHandle = new Simulation();

#ifdef WIN
	{
		LARGE_INTEGER qwTicksPerSec;
		usePerfCounter = QueryPerformanceFrequency(&qwTicksPerSec);
		if (usePerfCounter) {
			QueryPerformanceCounter(&perfTickStart);
			perfTicksPerSec = (double)qwTicksPerSec.QuadPart;
		}
		tickStart = GetTickCount();
	}
#else
	tickStart = (DWORD)time(NULL);
#endif

}

void ClearSimulation() {

    delete sHandle;
    sHandle = new Simulation;

}

DWORD RevertBit(DWORD dw) {
	DWORD dwIn = dw;
	DWORD dwOut = 0;
	DWORD dWMask = 1;
	int i;
	for (i = 0; i < 32; i++) {
		if (dwIn & 0x80000000UL) dwOut |= dWMask;
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

bool LoadSimulation(Dataport *loader) {

	double t1, t0;
	DWORD seed;
	t0 = GetTick();

	sHandle->loadOK = false;
	SetState(PROCESS_STARTING, "Clearing previous simulation");
	try {
		ClearSimulation();
	}
	catch (...) {
		SetErrorSub("Error clearing geometry");
		return false;
	}
	/* //Mutex not necessary: by the time the COMMAND_LOAD is issued the interface releases the handle, concurrent reading is safe and it's only destroyed by the interface when all processes are ready loading
	//Result: faster, parallel loading
	// Connect the dataport
	SetState(PROCESS_STARTING, "Waiting for 'loader' dataport access...");
	if( !AccessDataportTimed(loader,40000) ) {
		SetErrorSub("Failed to connect to DP");
		return false;
	}
	*/

    SetState(PROCESS_STARTING, "Loading simulation");
    sHandle->textTotalSize =
    sHandle->profTotalSize =
    sHandle->dirTotalSize =
            sHandle->spectrumTotalSize = 0;

    {

        std::string inputString(loader->size,NULL);
        BYTE* buffer = (BYTE*)loader->buff;
        std::copy(buffer, buffer + loader->size, inputString.begin());
        std::stringstream inputStream;
        inputStream << inputString;
        cereal::BinaryInputArchive inputarchive(inputStream);

        //std::ifstream is("data.json");
        //cereal::JSONInputArchive inputarchive(is);

        //Worker params
        inputarchive(sHandle->wp);
        sHandle->regions.resize(sHandle->wp.nbRegion); //Create structures
        inputarchive(sHandle->ontheflyParams);
        inputarchive(sHandle->regions); // mathonly
        inputarchive(sHandle->materials);
        inputarchive(sHandle->psi_distro);
        inputarchive(sHandle->chi_distros);
        inputarchive(sHandle->parallel_polarization);

        //Geometry
        inputarchive(sHandle->sh);
        inputarchive(sHandle->vertices3);

        // Prepare super structure
        sHandle->structures.resize(sHandle->sh.nbSuper); //Create structures
        /*for (size_t i = 0; i < sHandle->sh.nbSuper; i++) {
            int nbF = sHandle->structures[i].facets.size();
            if (nbF == 0) {

            }
            else {
                std::vector<SubprocessFacet>(nbF).swap(sHandle->structures[i].facets);
            }
        }*/

        SetState(PROCESS_STARTING, ("Loading facets"));

        //Facets
        for (size_t i = 0; i < sHandle->sh.nbFacet; i++) { //Necessary because facets is not (yet) a vector in the interface
            SubprocessFacet f;
            inputarchive(f.sh);
            inputarchive(f.indices);
            inputarchive(f.vertices2);
            inputarchive(f.textureCellIncrements);

            //Some initialization
            if (!f.InitializeOnLoad(i)) return false;
            if (f.sh.superIdx == -1) { //Facet in all structures
                for (auto& s : sHandle->structures) {
                    s.facets.push_back(f);
                }
            }
            else {
                sHandle->structures[f.sh.superIdx].facets.push_back(f); //Assign to structure
            }
        }
    }//inputarchive goes out of scope, file released

    sHandle->wp.nbTrajPoints = 0;

    if (sHandle->sh.nbSuper <= 0) {
        //ReleaseDataport(loader);
        SetErrorSub("No structures");
        return false;
    }

    try {
        //copy trajectory points
        for (auto& reg : sHandle->regions) {
            reg.Points.resize(reg.params.nbPointsToCopy);
            sHandle->wp.nbTrajPoints += reg.params.nbPointsToCopy;
        }
        //copy distribution points

        sHandle->nbDistrPoints_BXY = 0;
        for (auto& reg : sHandle->regions) {
            reg.Bx_distr.Resize((size_t)reg.params.nbDistr_MAG.x);
            reg.By_distr.Resize((size_t)reg.params.nbDistr_MAG.y);
            reg.Bz_distr.Resize((size_t)reg.params.nbDistr_MAG.z);
            reg.latticeFunctions.Resize(0);

            sHandle->nbDistrPoints_BXY += reg.params.nbDistr_BXY;
        }
    }
    catch (...) {
        SetErrorSub("Error loading regions");
        return false;
    }



    //Reserve particle log
    if (sHandle->ontheflyParams.enableLogging)
        sHandle->tmpParticleLog.reserve(sHandle->ontheflyParams.logLimit / sHandle->ontheflyParams.nbProcess);


    /*BYTE *bufferStart = (BYTE *)loader->buff;
    BYTE *buffer = bufferStart;
    BYTE *areaBuff;

	// Load new geom from the dataport

    WorkerParams* wPms = (WorkerParams *)buffer;
    sHandle->wp.newReflectionModel = wPms->newReflectionModel;
    sHandle->wp.nbRegion = wPms->nbRegion;
    buffer += sizeof(WorkerParams);

	GeomProperties* shGeom = (GeomProperties *)buffer;
	if (shGeom->nbSuper > MAX_STRUCT) {
		//ReleaseDataport(loader);
		SetErrorSub("Too many structures");
		return false;
	}
	if (shGeom->nbSuper <= 0) {
		//ReleaseDataport(loader);
		SetErrorSub("No structures");
		return false;
	}*/

    /*sHandle->totalFacet = shGeom->nbFacet;
    sHandle->nbVertex = shGeom->nbVertex;
    sHandle->nbSuper = shGeom->nbSuper;*/


/*
	buffer += sizeof(GeomProperties);

    sHandle->sh = READBUFFER(GeomProperties); //Copy all geometry properties
    sHandle->structures.resize(sHandle->sh.nbSuper); //Create structures

    sHandle->ontheflyParams = READBUFFER(OntheflySimulationParams);

	sHandle->wp.nbTrajPoints = 0;

	// RegionParams / Region_mathonly + Trajectory points
	// Distribution Points Bx By Bz LatticeFunctions
	// Material
	// chi_distros
	// psi_distro
	// parallel_polarization

	// FacetProperties
	// SubprocessFacet


	//Regions
	try {
		for (size_t r = 0; r < sHandle->wp.nbRegion; r++) {
			RegionParams *regparam = (RegionParams *)buffer;
			Region_mathonly newreg;
			newreg.params = *regparam;//memcpy(&newreg.params, regparam, sizeof(RegionParams));
			buffer += sizeof(RegionParams);
			sHandle->regions.push_back(newreg);
		}
		//copy trajectory points	
		for (size_t r = 0; r < sHandle->wp.nbRegion; r++) {

			sHandle->regions[r].Points.resize(sHandle->regions[r].params.nbPointsToCopy);
			memcpy(&sHandle->regions[r].Points[0], buffer, sizeof(Trajectory_Point)*sHandle->regions[r].params.nbPointsToCopy);
			buffer += sizeof(Trajectory_Point)*sHandle->regions[r].params.nbPointsToCopy;
			sHandle->wp.nbTrajPoints += sHandle->regions[r].params.nbPointsToCopy;
		}

		//copy distribution points
		//sHandle->nbDistrPoints_MAG=0;
		sHandle->nbDistrPoints_BXY = 0;
		for (size_t r = 0; r < sHandle->wp.nbRegion; r++) {

			sHandle->regions[r].Bx_distr.Resize((size_t)sHandle->regions[r].params.nbDistr_MAG.x);
			for (size_t j = 0; j < sHandle->regions[r].params.nbDistr_MAG.x; j++) {
				double x = READBUFFER(double);
				double y = READBUFFER(double);
				sHandle->regions[r].Bx_distr.SetPair(j, x, y);
			}

			sHandle->regions[r].By_distr.Resize((size_t)sHandle->regions[r].params.nbDistr_MAG.y);
			for (size_t j = 0; j < sHandle->regions[r].params.nbDistr_MAG.y; j++) {
				double x = READBUFFER(double);
				double y = READBUFFER(double);
				sHandle->regions[r].By_distr.SetPair(j, x, y);
			}

			sHandle->regions[r].Bz_distr.Resize((size_t)sHandle->regions[r].params.nbDistr_MAG.z);
			for (size_t j = 0; j < sHandle->regions[r].params.nbDistr_MAG.z; j++) {
				double x = READBUFFER(double);
				double y = READBUFFER(double);
				sHandle->regions[r].Bz_distr.SetPair(j, x, y);
			}

			sHandle->regions[r].latticeFunctions.Resize(0);

			for (size_t j = 0; j < sHandle->regions[r].params.nbDistr_BXY; j++) {
				double x = READBUFFER(double);
				std::vector<double> betaValues;
				betaValues.resize(6);
				memcpy(betaValues.data(), buffer, 6 * sizeof(double)); //Vector data is contigious. The 6 values are: BetaX, BetaY, EtaX, EtaX', AlphaX, AlphaY
				buffer += 6 * sizeof(double);
				//sHandle->regions[r].latticeFunctions.AddPair(x, betaValues);
			}

			sHandle->nbDistrPoints_BXY += sHandle->regions[r].params.nbDistr_BXY;
		}
	}
	catch (...) {
		SetErrorSub("Error loading regions");
		return false;
	}

	//Load materials
	try {
		sHandle->nbMaterials = READBUFFER(size_t); //copying number of materials
		for (size_t i = 0; i < sHandle->nbMaterials; i++) { //going through all materials
			Material newMaterial;
			newMaterial.hasBackscattering = READBUFFER(bool);
			size_t angleValsSize = READBUFFER(size_t); //copying number of angles (columns)
			size_t energyValsSize = READBUFFER(size_t); //copying number of energies (rows)
			newMaterial.angleVals.reserve(angleValsSize);
			for (size_t j = 0; j < angleValsSize; j++) {
				newMaterial.angleVals.push_back(*((double*)buffer)); //copying angles (header)
				buffer += sizeof(double);
			}
			newMaterial.energyVals.reserve(energyValsSize);
			for (size_t j = 0; j < energyValsSize; j++) {
				newMaterial.energyVals.push_back(*((double*)buffer)); //copying energies (column1)
				buffer += sizeof(double);
			}
			//copying reflectivity probabilities (cells)
			for (size_t j = 0; j < newMaterial.energyVals.size(); j++) {
				std::vector<std::vector<double>> currentEnergy;
				for (size_t k = 0; k < newMaterial.angleVals.size(); k++) {
					std::vector<double> reflProbability;
					reflProbability.push_back(*((double*)buffer)); //forward scattering
					buffer += sizeof(double);
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
		return false;
	}


	size_t totalDistroSize = READBUFFER(size_t);

	if (sHandle->chi_distros.size() == 3) buffer += totalDistroSize; //already loaded, skip
	else {
		sHandle->chi_distros.resize(3);
		try {
			//Load psi_distr
			size_t psi_size = READBUFFER(size_t);
			sHandle->psi_distro.reserve(psi_size);
			for (size_t j = 0; j < psi_size; j++) {
				std::vector<double> row;
				size_t row_size = READBUFFER(size_t);
				row.reserve(row_size);
				for (size_t k = 0; k < row_size; k++) {
					row.push_back(*((double*)buffer)); //cum. distr
					buffer += sizeof(double);
				}
				sHandle->psi_distro.push_back(row);
			}
			for (size_t c = 0; c < 3; c++) { //3 polarization components

				//Load chi_distr
				size_t chi_size = READBUFFER(size_t);
				sHandle->chi_distros[c].reserve(chi_size);
				for (size_t j = 0; j < chi_size; j++) {
					std::vector<double> row;
					size_t row_size = READBUFFER(size_t);
					row.reserve(row_size);
					for (size_t k = 0; k < row_size; k++) {
						row.push_back(*((double*)buffer)); //cum. distr
						buffer += sizeof(double);
					}
					sHandle->chi_distros[c].push_back(row);
				}
			}
			//Load polarization distribution
			psi_size = READBUFFER(size_t);
			sHandle->parallel_polarization.reserve(psi_size);
			for (size_t j = 0; j < psi_size; j++) {
				std::vector<double> row;
				size_t row_size = READBUFFER(size_t);
				row.reserve(row_size);
				for (size_t k = 0; k < row_size; k++) {
					row.push_back(*((double*)buffer)); //cum. distr
					buffer += sizeof(double);
				}
				sHandle->parallel_polarization.push_back(row);
			}
		}
		catch (...) {
			SetErrorSub("Error loading distributions");
			return false;
		}
	}


    BYTE *bufferAfterMaterials = buffer;

	// Prepare super structure

    sHandle->structures.resize(sHandle->sh.nbSuper); //Create structures

    buffer += sizeof(Vector3d)*sHandle->sh.nbVertex;
	for (size_t i = 0; i < sHandle->sh.nbFacet; i++) {
		FacetProperties *shFacet = (FacetProperties *)buffer;

		buffer += sizeof(FacetProperties) + shFacet->nbIndex*(sizeof(size_t) + sizeof(Vector2d));
	}
	for (size_t i = 0; i < sHandle->sh.nbSuper; i++) {
		int nbF = sHandle->structures[i].facets.size();
		if (nbF == 0) {

		}
		else {
            std::vector<SubprocessFacet>(nbF).swap(sHandle->structures[i].facets);
		}
	}

    buffer = bufferAfterMaterials;
    areaBuff = buffer;
	//buffer = (BYTE *)loader->buff;

	// Name
	sHandle->sh.name = shGeom->name;

	// Vertex
    try {
        sHandle->vertices3.resize(sHandle->sh.nbVertex);
    }
    catch (...) {
        SetErrorSub("Not enough memory to load vertices");
        return false;
    }


    memcpy(sHandle->vertices3.data(), buffer, sHandle->sh.nbVertex * sizeof(Vector3d));
    buffer += sizeof(Vector3d)*sHandle->sh.nbVertex;

	// Facets
	for (size_t i = 0; i < sHandle->sh.nbFacet; i++) {

		SubprocessFacet f;
		f.sh = READBUFFER(FacetProperties);

        f.globalId = i;

		if (f.sh.superIdx == -1) {
			for (size_t s = 0; s < sHandle->sh.nbSuper; s++) {
				SubprocessFacet f_copy = f;
				if (s > 0) { //Create copy
					f_copy = SubprocessFacet(f);
				}
				sHandle->structures[s].facets[sHandle->structures[s].facets.size()] = f_copy;
			}
		}
		else {
			sHandle->structures[f.sh.superIdx].facets[sHandle->structures[f.sh.superIdx].facets.size()] = f;
		}
		//sHandle->str[idx].facets[sHandle->str[idx].nbFacet]->globalId = i;

		// Reset counter in local memory
		//memset(&(f.tmpCounter), 0, sizeof(FacetHitBuffer));
		//buffer += sizeof(FacetProperties);

		// Indicies
        f.indices.resize(f.sh.nbIndex);
        memcpy(f.indices.data(), buffer, f.sh.nbIndex * sizeof(size_t));
        buffer += f.sh.nbIndex * sizeof(size_t);

        // Vertices2
        try {
            f.vertices2.resize(f.sh.nbIndex);
        } catch (...) {
            SetErrorSub("Not enough memory to load vertices");
            return false;
        }
        memcpy(f.vertices2.data(), buffer, f.sh.nbIndex * sizeof(Vector2d));
        buffer += f.sh.nbIndex * sizeof(Vector2d);

        //Some initialization
        if (!f.InitializeOnLoad(i)) return false;
        if (f.sh.superIdx == -1) { //Facet in all structures
            for (auto& s : sHandle->structures) {
                s.facets.push_back(f);
            }
        }
        else {
            sHandle->structures[f.sh.superIdx].facets.push_back(f); //Assign to structure
        }
	}

    //Reserve particle log
    if (sHandle->ontheflyParams.enableLogging)
        sHandle->tmpParticleLog.reserve(sHandle->ontheflyParams.logLimit / sHandle->ontheflyParams.nbProcess);
*/
    //ReleaseDataport(loader); //Commented out as AccessDataport removed

	// Build all AABBTrees
	size_t maxDepth = 0;

    for (auto& s : sHandle->structures) {
        std::vector<SubprocessFacet*> facetPointers; facetPointers.reserve(s.facets.size());
        for (auto& f : s.facets) {
            facetPointers.push_back(&f);
        }
        s.aabbTree = BuildAABBTree(facetPointers, 0, maxDepth);
    }

	// Initialise simulation
	ComputeSourceArea();
	seed = GetSeed();
	rseed(seed);

	//--- GSL random init ---
    gsl_rng_env_setup();                          // Read variable environnement
    const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
    sHandle->gen = gsl_rng_alloc(type);          // Rand generator allocation

	sHandle->loadOK = true;
	t1 = GetTick();
	printf("  Load %s successful\n", sHandle->sh.name.c_str());
	printf("  Geometry: %zd vertex %zd facets\n", sHandle->sh.nbVertex, sHandle->sh.nbFacet);
	printf("  Region: %zd regions\n", sHandle->wp.nbRegion);
	printf("  Trajectory points: %zd points\n", sHandle->wp.nbTrajPoints);
	printf("  Geom size: %d bytes\n", /*(size_t)(buffer - bufferStart)*/0);
	printf("  Number of stucture: %zd\n", sHandle->sh.nbSuper);
	printf("  Global Hit: %zd bytes\n", sizeof(GlobalHitBuffer));
	printf("  Facet Hit : %zd bytes\n", sHandle->sh.nbFacet*(int)sizeof(FacetHitBuffer));
	printf("  Texture   : %zd bytes\n", sHandle->textTotalSize);
	printf("  Profile   : %zd bytes\n", sHandle->profTotalSize);
	printf("  Direction : %zd bytes\n", sHandle->dirTotalSize);
	printf("  Spectrum  : %zd bytes\n", sHandle->spectrumTotalSize);
	printf("  Total     : %zd bytes\n", GetHitsSize());
	printf("  Seed: %u\n", seed);
	printf("  Loading time: %.3f ms\n", (t1 - t0)*1000.0);
	return true;

}

bool UpdateOntheflySimuParams(Dataport *loader) {
	// Connect the dataport
	if (!AccessDataportTimed(loader, 2000)) {
		SetErrorSub("Failed to connect to loader DP");
		return false;
	}
	BYTE* buffer = (BYTE *)loader->buff;

	sHandle->ontheflyParams = READBUFFER(OntheflySimulationParams);

	for (size_t i = 0; i < sHandle->wp.nbRegion; i++) {
		sHandle->regions[i].params.showPhotons = READBUFFER(bool);
	}

	ReleaseDataport(loader);

	//Reset hit cache
    sHandle->tmpGlobalResult.hitCacheSize = 0;
	//memset(sHandle->hitCache, 0, sizeof(HIT)*HITCACHESIZE);
    sHandle->tmpGlobalResult.leakCacheSize = 0;
	//memset(sHandle->leakCache, 0, sizeof(LEAK)*LEAKCACHESIZE); //No need to reset, will gradually overwrite

	return true;
}

void UpdateHits(Dataport *dpHit, Dataport* dpLog, int prIdx, DWORD timeout) {

	UpdateMCHits(dpHit, prIdx, timeout);
	if (dpLog) UpdateLog(dpLog, timeout);
}

size_t GetHitsSize() {
	return sHandle->textTotalSize + sHandle->profTotalSize + sHandle->dirTotalSize +
		sHandle->spectrumTotalSize + sHandle->sh.nbFacet * sizeof(FacetHitBuffer) + sizeof(GlobalHitBuffer);
}

void ResetTmpCounters() {
	SetState(NULL, "Resetting local cache...", false, true);

	memset(&sHandle->tmpGlobalResult, 0, sizeof(GlobalHitBuffer));

	sHandle->currentParticle.distTraveledSinceUpdate = 0.0;
	sHandle->nbLeakSinceUpdate = 0;
	sHandle->tmpGlobalResult.hitCacheSize = 0;
	sHandle->tmpGlobalResult.leakCacheSize = 0;

    for (int j = 0; j < sHandle->sh.nbSuper; j++) {
        for (auto& f : sHandle->structures[j].facets) {

            f.ResetCounter();
            f.hitted = false;

            std::vector<TextureCell>(f.texture.size()).swap(f.texture);
            std::vector<ProfileSlice>(f.profile.size()).swap(f.profile);
            //if (f.sh.recordSpectrum)
                f.spectrum.ResetCounts();
            std::vector<DirectionCell>(f.direction.size()).swap(f.direction);



        }
    }
}

void ResetSimulation() {
    sHandle->currentParticle.lastHitFacet = NULL;
	sHandle->totalDesorbed = 0;
	ResetTmpCounters();
	sHandle->tmpParticleLog.clear();
}

bool StartSimulation() {
	if (sHandle->regions.size() == 0) {
		SetErrorSub("No regions");
		return false;
	}

	//Check for invalid material (ID==9) passed from the interface
	//It is valid to load invalid materials (to extract textures, etc.) but not to launch the simulation
    for(auto& s : sHandle->structures) {
		for (auto& f : s.facets) {
			if (f.sh.reflectType == 9) {
				char tmp[32];
				sprintf(tmp, "Invalid material on Facet %zd.", f.globalId + 1);
				SetErrorSub(tmp);
				return false;
			}
		}
	}

	//if (!sHandle->lastHitFacet) StartFromSource();
    if (!sHandle->currentParticle.lastHitFacet) StartFromSource();
    //return (sHandle->currentParticle.lastHitFacet != NULL);
	return true;
}

void RecordHit(const int &type, const double &dF, const double &dP) {
	if (sHandle->regions[sHandle->sourceRegionId].params.showPhotons) {
        if (sHandle->tmpGlobalResult.hitCacheSize < HITCACHESIZE) {
            sHandle->tmpGlobalResult.hitCache[sHandle->tmpGlobalResult.hitCacheSize].pos = sHandle->currentParticle.position;
            sHandle->tmpGlobalResult.hitCache[sHandle->tmpGlobalResult.hitCacheSize].type = type;
            sHandle->tmpGlobalResult.hitCache[sHandle->tmpGlobalResult.hitCacheSize].dF = sHandle->currentParticle.dF;
            sHandle->tmpGlobalResult.hitCache[sHandle->tmpGlobalResult.hitCacheSize].dP = sHandle->currentParticle.dP;

            sHandle->tmpGlobalResult.hitCacheSize++;
        }
	}
}

void RecordLeakPos() {
	// Source region check performed when calling this routine 
	// Record leak for debugging
	RecordHit(HIT_REF, sHandle->currentParticle.dF, sHandle->currentParticle.dP);
	RecordHit(HIT_LAST, sHandle->currentParticle.dF, sHandle->currentParticle.dP);
    if (sHandle->tmpGlobalResult.leakCacheSize < LEAKCACHESIZE) {
        sHandle->tmpGlobalResult.leakCache[sHandle->tmpGlobalResult.leakCacheSize].pos = sHandle->currentParticle.position;
        sHandle->tmpGlobalResult.leakCache[sHandle->tmpGlobalResult.leakCacheSize].dir = sHandle->currentParticle.direction;
        sHandle->tmpGlobalResult.leakCacheSize++;
    }
}

bool SimulationRun() {

	// 1s step
	double t0, t1;
	size_t    nbStep = 1;
	bool   goOn;

	if (sHandle->stepPerSec == 0.0) {

		nbStep = 250;

	}

	if (sHandle->stepPerSec != 0.0)
		nbStep = (size_t)(sHandle->stepPerSec + 0.5);
	if (nbStep < 1) nbStep = 1;
	t0 = GetTick();

	goOn = SimulationMCStep(nbStep);

	t1 = GetTick();
	sHandle->stepPerSec = (double)(nbStep) / (t1 - t0);
#ifdef _DEBUG
	printf("Running: stepPerSec = %f\n", sHandle->stepPerSec);
#endif

	return !goOn;

}

double GetTick() {

	// Number of sec since the application startup

#ifdef WIN

	if (usePerfCounter) {

		LARGE_INTEGER t, dt;
		QueryPerformanceCounter(&t);
		dt.QuadPart = t.QuadPart - perfTickStart.QuadPart;
		return (double)(dt.QuadPart) / perfTicksPerSec;

	}
	else {

		return (double)((GetTickCount() - tickStart) / 1000.0);

	}

#else

	if (tickStart < 0)
		tickStart = time(NULL);

	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((double)(tv.tv_sec - tickStart)*1000.0 + (double)tv.tv_usec / 1000.0);

#endif

}

bool SubprocessFacet::InitializeOnLoad(const size_t& id) {
    globalId = id;
    if (!InitializeLinkAndVolatile(id)) return false;
    if (!InitializeTexture()) return false;
    if (!InitializeProfile()) return false;
    if (!InitializeDirectionTexture()) return false;
    if (!InitializeSpectrum()) return false;

    return true;
}

bool SubprocessFacet::InitializeProfile() {
    //Profiles
    if (sh.isProfile) {
        profileSize = PROFILE_SIZE * sizeof(ProfileSlice);
        try {
            profile = std::vector<ProfileSlice>(PROFILE_SIZE);
        }
        catch (...) {
            SetErrorSub("Not enough memory to load profiles");
            return false;
        }
        sHandle->profTotalSize += profileSize * (1);
    }
    else profileSize = 0;
    return true;
}

bool SubprocessFacet::InitializeTexture(){
    //Textures
    if (sh.isTextured) {
        size_t nbE = sh.texWidth*sh.texHeight;
        textureSize = nbE*sizeof(TextureCell);
        try {
            texture.resize(nbE);
            textureCellIncrements.resize(nbE);
            largeEnough.resize(nbE);
        }
        catch (...) {
            SetErrorSub("Not enough memory to load textures");
            return false;
        }


        fullSizeInc = 1E30;

        for (size_t j = 0; j < nbE; j++) {
            if ((textureCellIncrements[j] > 0.0) && (textureCellIncrements[j] < fullSizeInc)) fullSizeInc = textureCellIncrements[j];
        }
        for (size_t j = 0; j < nbE; j++) { //second pass, filter out very small cells
            largeEnough[j] = (textureCellIncrements[j] < ((5.0f)*fullSizeInc));
        }
        sHandle->textTotalSize += textureSize;
        iw = 1.0 / (double)sh.texWidthD;
        ih = 1.0 / (double)sh.texHeightD;
        rw = sh.U.Norme() * iw;
        rh = sh.V.Norme() * ih;
    }
    else textureSize = 0;
    return true;
}

bool SubprocessFacet::InitializeDirectionTexture(){
    //Direction
    if (sh.countDirection) {
        directionSize = sh.texWidth*sh.texHeight * sizeof(DirectionCell);
        try {
            direction = std::vector<DirectionCell>(directionSize);
        }
        catch (...) {
            SetErrorSub("Not enough memory to load direction textures");
            return false;
        }
        sHandle->dirTotalSize += directionSize;
    }
    else directionSize = 0;
    return true;
}

bool SubprocessFacet::InitializeSpectrum(){
    //Spectrum
    if (sh.recordSpectrum) {
        spectrumSize = sizeof(ProfileSlice)*SPECTRUM_SIZE;
        double min_energy, max_energy;
        if (sHandle->wp.nbRegion > 0) {
            min_energy = sHandle->regions[0].params.energy_low_eV;
            max_energy = sHandle->regions[0].params.energy_hi_eV;
        }
        else {
            min_energy = 10.0;
            max_energy = 1000000.0;
        }

        try {
            spectrum = Histogram(min_energy, max_energy, SPECTRUM_SIZE, true);
            //spectrum = std::vector<ProfileSlice>(SPECTRUM_SIZE);
        }
        catch (...) {
            SetErrorSub("Not enough memory to load spectrum");
            return false;
        }
        sHandle->spectrumTotalSize += spectrumSize;
    }
    else spectrumSize = 0;

    return true;
}

bool SubprocessFacet::InitializeLinkAndVolatile(const size_t & id){
    sHandle->hasVolatile |= sh.isVolatile;

    if (sh.superDest || sh.isVolatile) {
        // Link or volatile facet, overides facet settings
        // Must be full opaque and 0 sticking
        // (see SimulationMC.c::PerformBounce)
        //sh.isOpaque = true;
        sh.opacity = 1.0;
        sh.sticking = 0.0;
        if (((sh.superDest - 1) >= sHandle->sh.nbSuper || sh.superDest < 0)) {
            // Geometry error
            ClearSimulation();
            //ReleaseDataport(loader);
            std::ostringstream err;
            err << "Invalid structure (wrong link on F#" << id + 1 << ")";
            SetErrorSub(err.str().c_str());
            return false;
        }
    }
    return true;
}