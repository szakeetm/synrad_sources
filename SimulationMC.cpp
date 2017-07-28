/*
File:        SimulationMC.c
Description: Monte-Carlo Simulation for UHV (Physics related routines)
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Simulation.h"
#include "Random.h"
#include "Distributions.h"
#include <vector>
#include "Region_mathonly.h"
#include "GeneratePhoton.h"
#include "GLApp/MathTools.h"
#include "SynradTypes.h" //Histogram
#include <tuple>

extern SIMULATION *sHandle;
extern void SetErrorSub(const char *message);

extern Distribution2D K_1_3_distribution;
extern Distribution2D K_2_3_distribution;
extern Distribution2D integral_N_photons;
extern Distribution2D integral_SR_power;
extern Distribution2D polarization_distribution;
//extern Distribution2D g1h2_distribution;


void ComputeSourceArea() {
	sHandle->sourceArea = sHandle->nbTrajPoints;
}

// -------------------------------------------------------

void PolarToCartesian(FACET *iFacet, double theta, double phi, bool reverse, double rotateUV) {

	Vector3d U, V, N;
	double u, v, n;

	// Polar in (nU,nV,N) to Cartesian(x,y,z) transformation  ( nU = U/|U| , nV = V/|V| )
	// tetha is the angle to the normal of the facet N, phi to U
	// ! See Geometry::InitializeGeometry() for further informations on the (U,V,N) basis !
	// (nU,nV,N) and (x,y,z) are both left handed

	/*#ifdef WIN
	_asm {                    // FPU stack
	fld qword ptr [theta]
	fsincos                 // cos(t)        sin(t)
	fld qword ptr [phi]
	fsincos                 // cos(p)        sin(p) cos(t) sin(t)
	fmul st(0),st(3)        // cos(p)*sin(t) sin(p) cos(t) sin(t)
	fstp qword ptr [u]      // sin(p)        cos(t) sin(t)
	fmul st(0),st(2)        // sin(p)*sin(t) cos(t) sin(t)
	fstp qword ptr [v]      // cos(t) sin(t)
	fstp qword ptr [n]      // sin(t)
	fstp qword ptr [dummy]  // Flush the sin(t)
	}
	#else*/
	u = sin(theta)*cos(phi);
	v = sin(theta)*sin(phi);
	n = cos(theta);
	//#endif

	// Get the (nU,nV,N) orthonormal basis of the facet
	U = iFacet->sh.nU;
	V = iFacet->sh.nV;
	N = iFacet->sh.N;
	if (reverse) {
		N.x = N.x*(-1.0);
		N.y = N.y*(-1.0);
		N.z = N.z*(-1.0);
	}

	/*if (rotateUV>0.0) {
		Rotate(&U, iFacet->sh.O, N, rotateUV);
		Rotate(&V, iFacet->sh.O, N, rotateUV);
		}*/

	// Basis change (nU,nV,N) -> (x,y,z)
	sHandle->pDir.x = u*U.x + v*V.x + n*N.x;
	sHandle->pDir.y = u*U.y + v*V.y + n*N.y;
	sHandle->pDir.z = u*U.z + v*V.z + n*N.z;
	//_ASSERTE(Norme(&sHandle->pDir)<=1.0);
}

// -------------------------------------------------------

void CartesianToPolar(FACET *iFacet, double *theta, double *phi) {

	// Get polar coordinates of the incoming particule direction in the (U,V,N) facet space.
	// Note: The facet is parallel to (U,V), we use its (nU,nV,N) orthonormal basis here.
	// (nU,nV,N) and (x,y,z) are both left handed

	// Cartesian(x,y,z) to polar in (nU,nV,N) transformation

	// Basis change (x,y,z) -> (nU,nV,N)
	// We use the fact that (nU,nV,N) belongs to SO(3)
	double u = Dot(sHandle->pDir, iFacet->sh.nU);
	double v = Dot(sHandle->pDir, iFacet->sh.nV);
	double n = Dot(sHandle->pDir, iFacet->sh.N);
	Saturate(n, -1.0, 1.0); //sometimes rounding errors do occur, 'acos' function would return no value for theta

	// (u,v,n) -> (theta,phi)
	double rho = sqrt(v*v + u*u);
	*theta = acos(n);              // Angle to normal (PI/2 => PI)
	*phi = asin(v / rho);
	if (u < 0.0) *phi = PI - *phi;  // Angle to U

}

void UpdateMCHits(Dataport *dpHit, int prIdx, DWORD timeout) {

	BYTE *buffer;
	SHGHITS *gHits;
	llong minHitsOld_MC;
	llong maxHitsOld_MC;
	double minHitsOld_flux, minHitsOld_power;
	double maxHitsOld_flux, maxHitsOld_power;
	int i, j, s, x, y, nb;
#ifdef _DEBUG
	double t0, t1;
	t0 = GetTick();
#endif
	SetState(NULL, "Waiting for 'hits' dataport access...", false, true);
	sHandle->lastUpdateOK = AccessDataportTimed(dpHit, timeout);
	SetState(NULL, "Updating MC hits...", false, true);
	if (!sHandle->lastUpdateOK) return;

	buffer = (BYTE*)dpHit->buff;
	gHits = (SHGHITS *)buffer;

	// Global hits and leaks
	gHits->total.nbHit += sHandle->tmpCount.nbHit;
	gHits->total.nbAbsorbed += sHandle->tmpCount.nbAbsorbed;
	gHits->total.nbDesorbed += sHandle->tmpCount.nbDesorbed;
	gHits->distTraveledTotal += sHandle->distTraveledSinceUpdate;
	gHits->total.fluxAbs += sHandle->tmpCount.fluxAbs;
	gHits->total.powerAbs += sHandle->tmpCount.powerAbs;

	minHitsOld_MC = gHits->minHit_MC;
	maxHitsOld_MC = gHits->maxHit_MC;
	minHitsOld_flux = gHits->minHit_flux;
	maxHitsOld_flux = gHits->maxHit_flux;
	minHitsOld_power = gHits->minHit_power;
	maxHitsOld_power = gHits->maxHit_power;

	gHits->minHit_MC = HITMAX_INT64;
	gHits->maxHit_MC = 0;
	gHits->minHit_flux = HITMAX_DOUBLE;
	gHits->maxHit_flux = 0;
	gHits->minHit_power = HITMAX_DOUBLE;
	gHits->maxHit_power = 0;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak
	for (size_t leakIndex = 0; leakIndex < sHandle->leakCacheSize; leakIndex++)
		gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = sHandle->leakCache[leakIndex];
	gHits->nbLeakTotal += sHandle->nbLeakSinceUpdate;
	gHits->lastLeakIndex = (gHits->lastLeakIndex + sHandle->leakCacheSize) % LEAKCACHESIZE;
	gHits->leakCacheSize = Min(LEAKCACHESIZE, gHits->leakCacheSize + sHandle->leakCacheSize);

	// HHit (Only prIdx 0)
	if (prIdx == 0) {
		for (size_t hitIndex = 0; hitIndex < sHandle->hitCacheSize; hitIndex++)
			gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = sHandle->hitCache[hitIndex];

		if (sHandle->hitCacheSize > 0) {
			gHits->lastHitIndex = (gHits->lastHitIndex + sHandle->hitCacheSize) % HITCACHESIZE;

			if (gHits->lastHitIndex < (HITCACHESIZE - 1)) {
				gHits->lastHitIndex++;
				gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
			}

			gHits->hitCacheSize = Min(HITCACHESIZE, gHits->hitCacheSize + sHandle->hitCacheSize);
		}
	}

	// Facets
	for (s = 0; s < sHandle->nbSuper; s++) {
		for (i = 0; i < sHandle->str[s].nbFacet; i++) {

			FACET *f = sHandle->str[s].facets[i];
			if (f->hitted) {

				SHHITS *fFit = (SHHITS *)(buffer + f->sh.hitOffset);
				fFit->nbAbsorbed += f->counter.nbAbsorbed;
				fFit->nbDesorbed += f->counter.nbDesorbed;
				fFit->fluxAbs += f->counter.fluxAbs;
				fFit->powerAbs += f->counter.powerAbs;
				fFit->nbHit += f->counter.nbHit;

				if (f->sh.isProfile) {
					llong *shProfile = (llong *)(buffer + (f->sh.hitOffset + sizeof(SHHITS)));
					for (j = 0; j < PROFILE_SIZE; j++)
						shProfile[j] += f->profile_hits[j];

					double *shProfile2 = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS)) + PROFILE_SIZE*sizeof(llong));
					for (j = 0; j < PROFILE_SIZE; j++)
						shProfile2[j] += f->profile_flux[j];

					double *shProfile3 = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS)) + PROFILE_SIZE*(sizeof(double) + sizeof(llong)));
					for (j = 0; j < PROFILE_SIZE; j++)
						shProfile3[j] += f->profile_power[j];
				}


				int profileSize = (f->sh.isProfile) ? PROFILE_SIZE*(2 * sizeof(double) + sizeof(llong)) : 0;
				if (f->sh.isTextured) {
					int textureElements = f->sh.texHeight*f->sh.texWidth;
					llong *shTexture = (llong *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize));
					for (y = 0; y < f->sh.texHeight; y++) {
						for (x = 0; x < f->sh.texWidth; x++) {
							int add = x + y*f->sh.texWidth;
							llong val = shTexture[add] + f->hits_MC[add];
							if (val > gHits->maxHit_MC)	//no cell size check for MC						
								gHits->maxHit_MC = val;
							if (val > 0 && val < gHits->minHit_MC)
								gHits->minHit_MC = val;
							shTexture[add] = val;
						}
					}

					double *shTexture2 = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize + textureElements*sizeof(llong)));
					for (y = 0; y < f->sh.texHeight; y++) {
						for (x = 0; x < f->sh.texWidth; x++) {
							int add = x + y*f->sh.texWidth;
							double val2 = shTexture2[add] + f->hits_flux[add];
							if (val2 > gHits->maxHit_flux&& f->largeEnough[add]) {
								gHits->maxHit_flux = val2;
							}
							if (val2 > 0.0 && val2 < gHits->minHit_flux&& f->largeEnough[add]) gHits->minHit_flux = val2;
							shTexture2[add] = val2;
						}
					}

					double *shTexture3 = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize + textureElements*(sizeof(llong) + sizeof(double))));
					for (y = 0; y < f->sh.texHeight; y++) {
						for (x = 0; x < f->sh.texWidth; x++) {
							int add = x + y*f->sh.texWidth;
							double val3 = shTexture3[add] + f->hits_power[add];
							if (val3 > gHits->maxHit_power&& f->largeEnough[add])
								gHits->maxHit_power = val3;
							if (val3 > 0.0 && val3 < gHits->minHit_power&& f->largeEnough[add]) gHits->minHit_power = val3; //disregard very small elements
							shTexture3[add] = val3;
						}
					}
				}

				if (f->sh.countDirection) {
					VHIT *shDir = (VHIT *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize + f->textureSize));
					for (y = 0; y < f->sh.texHeight; y++) {
						for (x = 0; x < f->sh.texWidth; x++) {
							int add = x + y*f->sh.texWidth;
							shDir[add].dir.x += f->direction[add].dir.x;
							shDir[add].dir.y += f->direction[add].dir.y;
							shDir[add].dir.z += f->direction[add].dir.z;
							shDir[add].count += f->direction[add].count;
						}
					}
				}

				if (f->sh.hasSpectrum) {
					double *shSpectrum_fluxwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize
						+ f->textureSize + f->directionSize));
					for (j = 0; j < SPECTRUM_SIZE; j++)
						shSpectrum_fluxwise[j] += f->spectrum_fluxwise->GetCount(j);

					double *shSpectrum_powerwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize
						+ f->textureSize + f->directionSize + SPECTRUM_SIZE*sizeof(double)));
					for (j = 0; j < SPECTRUM_SIZE; j++)
						shSpectrum_powerwise[j] += f->spectrum_powerwise->GetCount(j);
				}

			} // End if(hitted)
		} // End nbFacet
	} // End nbSuper

	//if there were no textures:
	if (gHits->minHit_MC == HITMAX_INT64) gHits->minHit_MC = minHitsOld_MC;
	if (gHits->maxHit_MC == 0) gHits->maxHit_MC = maxHitsOld_MC;
	if (gHits->minHit_flux == HITMAX_DOUBLE) gHits->minHit_flux = minHitsOld_flux;
	if (gHits->maxHit_flux == 0.0) gHits->maxHit_flux = maxHitsOld_flux;
	if (gHits->minHit_power == HITMAX_DOUBLE) gHits->minHit_power = minHitsOld_power;
	if (gHits->maxHit_power == 0.0) gHits->maxHit_power = maxHitsOld_power;

	ReleaseDataport(dpHit);
	ResetTmpCounters();
	extern char* GetSimuStatus();
	SetState(NULL, GetSimuStatus(), false, true);

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif

}






// -------------------------------------------------------------
// Compute particle teleport
// -------------------------------------------------------------

void PerformTeleport(FACET *iFacet) {

	double inPhi, inTheta;
	//Search destination
	FACET *destination;
	bool found = false;
	int destIndex;

	if (iFacet->sh.teleportDest == -1) {
		destIndex = sHandle->teleportedFrom;
		if (destIndex == -1) {
			/*char err[128];
			sprintf(err, "Facet %d tried to teleport to the facet where the particle came from, but there is no such facet.", iFacet->globalId + 1);
			SetErrorSub(err);*/
			RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
			sHandle->lastHit = iFacet;
			return; //LEAK
		}
	}
	else destIndex = iFacet->sh.teleportDest - 1;

	//Look in which superstructure is the destination facet:
	for (int i = 0; i < sHandle->nbSuper && (!found); i++) {
		for (int j = 0; j < sHandle->str[i].nbFacet && (!found); j++) {
			if (destIndex == sHandle->str[i].facets[j]->globalId) {
				destination = sHandle->str[i].facets[j];
				sHandle->curStruct = destination->sh.superIdx; //change current superstructure
				sHandle->teleportedFrom = iFacet->globalId; //memorize where the particle came from
				found = true;
			}
		}
	}
	if (!found) {
		/*char err[128];
		sprintf(err, "Teleport destination of facet %d not found (facet %d does not exist)", iFacet->globalId + 1, iFacet->sh.teleportDest);
		SetErrorSub(err);*/
		RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
		sHandle->lastHit = iFacet;
		return; //LEAK
	}

	// Count this hit as a transparent pass
	RecordHit(HIT_TELEPORT, sHandle->dF, sHandle->dP);
	if (iFacet->hits_MC && iFacet->sh.countTrans) RecordHitOnTexture(iFacet, sHandle->dF, sHandle->dP);


	// Relaunch particle from new facet
	CartesianToPolar(iFacet, &inTheta, &inPhi);
	PolarToCartesian(destination, inTheta, inPhi, false);
	// Move particle to teleport destination point
	sHandle->pPos.x = destination->sh.O.x + iFacet->colU*destination->sh.U.x + iFacet->colV*destination->sh.V.x;
	sHandle->pPos.y = destination->sh.O.y + iFacet->colU*destination->sh.U.y + iFacet->colV*destination->sh.V.y;
	sHandle->pPos.z = destination->sh.O.z + iFacet->colU*destination->sh.U.z + iFacet->colV*destination->sh.V.z;
	RecordHit(HIT_TELEPORT, sHandle->dF, sHandle->dP);
	sHandle->lastHit = destination;

	//Count hits on teleport facets (only TP source)
	//iFacet->counter.nbAbsorbed++;
	//destination->counter.nbDesorbed++;

	iFacet->counter.nbHit++;//destination->counter.nbHit++;
	iFacet->counter.fluxAbs += sHandle->dF;//destination->counter.fluxAbs+=sHandle->dF;
	iFacet->counter.powerAbs += sHandle->dP;//destination->counter.powerAbs+=sHandle->dP;
}

// -------------------------------------------------------------
// Perform nbStep simulation steps (a step is a bounce)
// -------------------------------------------------------------

bool SimulationMCStep(int nbStep) {

	FACET   *collidedFacet;
	double   d;
	bool     found;
	int      i;

	// Perform simulation steps
	for (i = 0; i < nbStep; i++) {

		found = Intersect(&(sHandle->pPos), &(sHandle->pDir), &d, &collidedFacet, sHandle->lastHit); //May decide reflection type

		if (found) {

			// Move particule to intersection point
			sHandle->pPos = sHandle->pPos + d*sHandle->pDir;
			sHandle->distTraveledSinceUpdate += d;

			if (collidedFacet->sh.teleportDest) {
				PerformTeleport(collidedFacet);
			}
			else {
				// Hit count
				sHandle->tmpCount.nbHit++;
				collidedFacet->counter.nbHit++;
				if (collidedFacet->sh.superDest) {	// Handle super structure link facet
					sHandle->curStruct = collidedFacet->sh.superDest - 1;
					// Count this hit as a transparent pass
					RecordHit(HIT_TRANS, sHandle->dF, sHandle->dP);
					if (collidedFacet->hits_MC && collidedFacet->sh.countTrans) RecordHitOnTexture(collidedFacet, sHandle->dF, sHandle->dP);
				}
				else { //Not superDest or Teleport
					if (sHandle->newReflectionModel) {
						//New reflection model (Synrad 1.4)
						double theta, phi;
						CartesianToPolar(collidedFacet, &theta, &phi);

						double stickingProbability;

						bool complexScattering; //forward/diffuse/back/transparent/stick
						std::vector<double> materialReflProbabilities; //0: forward reflection, 1: diffuse, 2: backscattering, 3: transparent, 100%-(0+1+2+3): absorption
						if (sHandle->dF == 0.0 || sHandle->dP == 0.0 || sHandle->energy < 1E-3) {
							//stick non real photons (from beam beginning)
							stickingProbability = 1.0; 
							complexScattering = false;
						}
						else if (collidedFacet->sh.reflectType < 2) {
							//Diffuse or Mirror
							stickingProbability = collidedFacet->sh.sticking;
							complexScattering = false;
						}
						else { //Material
							stickingProbability=GetStickingProbability(collidedFacet, theta, materialReflProbabilities, complexScattering);
						}
						if (!sHandle->mode.lowFluxMode || complexScattering) {
							int reflType = GetHardHitType(stickingProbability, materialReflProbabilities, complexScattering);
							if (reflType == REFL_ABSORB) {
								Stick(collidedFacet);
								if (!StartFromSource()) return false;
							}
							else PerformBounce_new(collidedFacet, theta, phi, reflType);
						}
						else {
							//Low flux mode and simple scattering:
							DoLowFluxReflection(collidedFacet, stickingProbability, theta, phi);
						} //end low flux mode
					}
					else {
						//Old reflection model (Synrad <=1.3)						
						if (sHandle->dF == 0.0 || sHandle->dP == 0.0 || sHandle->energy < 1E-3) { //stick non real photons (from beam beginning)
							Stick(collidedFacet);
							if (!StartFromSource()) return false;
						}
						else {
							double sigmaRatio = collidedFacet->sh.doScattering?collidedFacet->sh.rmsRoughness/collidedFacet->sh.autoCorrLength:0.0;
							if ((collidedFacet->sh.reflectType - 10) < (int)sHandle->materials.size()) { //found material type
																										//Generate incident angle
								Vector3d nU_rotated, N_rotated, nV_rotated;
								double n, u, v;
								bool reflected = false;
								do { //generate angles until reflecting away from surface
										
										double n_ori = Dot(sHandle->pDir, collidedFacet->sh.N); //original incident angle (negative if front collision, positive if back collision);
										do { //generate angles until incidence is from front
												
											PerturbateSurface(sigmaRatio, collidedFacet, nU_rotated, nV_rotated, N_rotated);
	
											//Calculate incident angles to bent surface (Cartesian to Polar routine modified)
											GetDirComponents(nU_rotated, nV_rotated, N_rotated, u, v, n);

										} while (n*n_ori <= 0.0); //regenerate random numbers if grazing angle would be over 90deg

									double rho = sqrt(v*v + u*u);
									double theta = acos(n);       // Angle to normal (PI/2 => PI if front collision, 0..PI/2 if back)
									double phi = asin(v / rho);
									if (u < 0.0) phi = PI - phi;  // Angle to U
										
									double stickingProbability;
									std::vector<double> materialReflProbabilities;
									bool complexScattering;
									if (collidedFacet->sh.reflectType < 2) {
										stickingProbability = collidedFacet->sh.sticking;
										complexScattering = false;
									}
									else {
										stickingProbability = GetStickingProbability(collidedFacet, theta, materialReflProbabilities, complexScattering);
									}
									if (!sHandle->mode.lowFluxMode) {
										reflected = DoOldRegularReflection(collidedFacet, stickingProbability, materialReflProbabilities,
											complexScattering, theta, phi, N_rotated, nU_rotated, nV_rotated);
									}
									else {
										reflected = DoLowFluxReflection(collidedFacet, stickingProbability, theta, phi, N_rotated, nU_rotated, nV_rotated);
									}
								} while (!reflected); //do it again if reflection wasn't successful (reflected against the surface due to roughness)
							}
							else {
								std::string err = "Facet " + (collidedFacet->globalId + 1);
								err += "reflection material type not found";
								SetErrorSub(err.c_str());
							}
						}
						
					}
				}
			}
		}
		else { // Leak (simulation error)
			sHandle->nbLeakSinceUpdate++;
			if (sHandle->regions[sHandle->sourceRegionId].params.showPhotons){
				RecordLeakPos();
			}
			if (!StartFromSource())
				// desorptionLimit reached
				return false;
		} //end intersect or leak
	} //end step
	return true;
}

double GetStickingProbability(FACET* collidedFacet, const double& theta, std::vector<double>& materialReflProbabilities, bool& complexScattering) {
	//Returns sticking probability, but also fills materialReflProbabilities and complexScattering (pass by reference)
	Material *mat = &(sHandle->materials[collidedFacet->sh.reflectType - 10]);
	complexScattering = mat->hasBackscattering;
	materialReflProbabilities = mat->Interpolate(sHandle->energy, abs(theta - PI / 2));
	double stickingProbability = complexScattering
		? 1.0 - materialReflProbabilities[0] - materialReflProbabilities[1] - materialReflProbabilities[2] //100% - forward - diffuse - back (transparent already excluded in Intersect() routine)
		: 1.0 - materialReflProbabilities[0]; //100% - forward (transparent already excluded in Intersect() routine)
	return stickingProbability;
}

int GetHardHitType(const double& stickingProbability,const std::vector<double>& materialReflProbabilities, const bool& complexScattering) {
	//Similar to Material::GetReflectionType, but transparent pass is excluded and treats Mirror/Diffuse surfaces too with single sticking factor
	double nonTransparentProbability = (complexScattering) ? 1.0 - materialReflProbabilities[3] : 1.0;
	double random = rnd()*nonTransparentProbability;
	if (random > (nonTransparentProbability - stickingProbability)) return REFL_ABSORB;
	else if (!complexScattering) {
		return REFL_FORWARD; //Not absorbed, so reflected
	}
	else { //complex scattering
		if (random < materialReflProbabilities[0]) return REFL_FORWARD; //forward reflection
		else if (random < (materialReflProbabilities[0] + materialReflProbabilities[1])) return REFL_DIFFUSE;
		else if (random < (materialReflProbabilities[0] + materialReflProbabilities[1] + materialReflProbabilities[2])) return REFL_BACK;
		else return REFL_BACK; //should never be the case. Trasparent already excluded in Intersect() routine, and in this "else" branch we know that it doesn't stick
	}
}

void PerturbateSurface(double& sigmaRatio, FACET* collidedFacet, Vector3d& nU_rotated, Vector3d& nV_rotated, Vector3d& N_rotated) {

	double rnd1 = rnd();
	double rnd2 = rnd(); //for debug
	double thetaOffset = atan(sigmaRatio*tan(PI*(rnd1 - 0.5)));
	double phiOffset = atan(sigmaRatio*tan(PI*(rnd2 - 0.5)));

	Vector3d nU_facet = Vector3d(collidedFacet->sh.nU.x, collidedFacet->sh.nU.y, collidedFacet->sh.nU.z);
	Vector3d nV_facet = Vector3d(collidedFacet->sh.nV.x, collidedFacet->sh.nV.y, collidedFacet->sh.nV.z);
	Vector3d N_facet = Vector3d(collidedFacet->sh.N.x, collidedFacet->sh.N.y, collidedFacet->sh.N.z);

	//Random rotation around N (to discard U orientation thus make scattering isotropic)
	double randomAngle = rnd() * 2 * PI;
	nU_facet = Rotate(nU_facet,Vector3d(0,0,0),N_facet, randomAngle);
	nV_facet = Rotate(nV_facet,Vector3d(0,0,0),N_facet, randomAngle);

	//Bending surface base vectors
	nU_rotated = Rotate(nU_facet,Vector3d(0,0,0),nV_facet, thetaOffset);
	N_rotated = Rotate(N_facet,Vector3d(0,0,0),nV_facet, thetaOffset);
	nU_rotated = Rotate(nU_rotated,Vector3d(0,0,0),nU_facet, phiOffset); //check if correct
	nV_rotated = Rotate(nV_facet,Vector3d(0,0,0),nU_facet, phiOffset);
	N_rotated = Rotate(N_rotated,Vector3d(0,0,0),nU_facet, phiOffset);
}

void GetDirComponents(Vector3d& nU_rotated, Vector3d& nV_rotated, Vector3d& N_rotated, double& u, double& v, double& n) {
	u = Dot(sHandle->pDir, nU_rotated);
	v = Dot(sHandle->pDir, nV_rotated);
	n = Dot(sHandle->pDir, N_rotated);
	Saturate(n, -1.0, 1.0); //sometimes rounding errors do occur, 'acos' function would return no value for theta
}



bool DoLowFluxReflection(FACET* collidedFacet, double stickingProbability, double theta, double phi, Vector3d N_rotated, Vector3d nU_rotated, Vector3d nV_rotated) {
	collidedFacet->counter.fluxAbs += sHandle->dF*stickingProbability;
	collidedFacet->counter.powerAbs += sHandle->dP*stickingProbability;
	//sHandle->distTraveledSinceUpdate += sHandle->distTraveledCurrentParticle;
	if (collidedFacet->hits_MC && collidedFacet->sh.countAbs) RecordHitOnTexture(collidedFacet,
		sHandle->dF*stickingProbability, sHandle->dP*stickingProbability);
	//okay, absorbed part recorded, let's see how much is left
	sHandle->oriRatio *= (1.0 - stickingProbability);
	if (sHandle->oriRatio < sHandle->mode.lowFluxCutoff) {//reflected part not important, throw it away
		RecordHit(HIT_ABS, sHandle->dF, sHandle->dP); //for hits and lines display
		return StartFromSource(); //false if maxdesorption reached
	}
	else { //reflect remainder
		sHandle->dF *= (1.0 - stickingProbability);
		sHandle->dP *= (1.0 - stickingProbability);
		if (sHandle->newReflectionModel) {
			PerformBounce_new(collidedFacet, theta, phi, REFL_FORWARD);
			return true;
		}
		else {
			return PerformBounce_old(collidedFacet, REFL_FORWARD, theta, phi, N_rotated, nU_rotated, nV_rotated);
		}
	}
}

bool DoOldRegularReflection(FACET* collidedFacet, double stickingProbability, const std::vector<double>& materialReflProbabilities, 
	const bool& complexScattering, double theta, double phi, Vector3d N_rotated, Vector3d nU_rotated, Vector3d nV_rotated) {
		int reflType;
		if (complexScattering) {
			reflType = GetHardHitType(stickingProbability, materialReflProbabilities, true);
			if (reflType == REFL_ABSORB) {
				Stick(collidedFacet);
				return StartFromSource(); //false if maxdesorption reached
			} else return PerformBounce_old(collidedFacet, reflType, theta, phi, N_rotated, nU_rotated, nV_rotated);
		}
		else { //reflect or absorb
			if (rnd() < stickingProbability) {
				Stick(collidedFacet);
				return StartFromSource(); //false if maxdesorption reached
			}
			else {
				reflType = REFL_FORWARD;
				return PerformBounce_old(collidedFacet, reflType, theta, phi, N_rotated, nU_rotated, nV_rotated);
			}
		}
		
}

// -------------------------------------------------------------
// Launch photon from a trajectory point
// -------------------------------------------------------------

bool StartFromSource() {

	// Check end of simulation
	if (sHandle->desorptionLimit > 0) {
		if (sHandle->totalDesorbed >= sHandle->desorptionLimit) {
			sHandle->lastHit = NULL;
			return false;
		}
	}

	//find source point
	int pointIdGlobal = (int)(rnd()*sHandle->sourceArea);
	bool found = false;
	int regionId;
	int pointIdLocal;
	int sum = 0;
	for (regionId = 0; regionId<sHandle->nbRegion&&!found; regionId++) {
		if ((pointIdGlobal >= sum) && (pointIdGlobal < (sum + (int)sHandle->regions[regionId].Points.size()))) {
			pointIdLocal = pointIdGlobal - sum;
			found = true;
		}
		else sum += (int)sHandle->regions[regionId].Points.size();
	}
	if (!found) SetErrorSub("No start point found");
	regionId--;
	//Trajectory_Point *source=&(sHandle->regions[regionId].Points[pointIdLocal]);
	Region_mathonly *sourceRegion = &(sHandle->regions[regionId]);
	if (!(sourceRegion->params.psimaxX_rad > 0.0 && sourceRegion->params.psimaxY_rad>0.0)) SetErrorSub("psiMaxX or psiMaxY not positive. No photon can be generated");
	
	size_t retries = 0;bool validEnergy;GenPhoton photon;
	do {
		photon = GeneratePhoton(pointIdLocal, sourceRegion, sHandle->mode.generation_mode,
			sHandle->psi_distro, sHandle->chi_distros[sourceRegion->params.polarizationCompIndex],
			sHandle->parallel_polarization, sHandle->tmpCount.nbDesorbed == 0);
		validEnergy = (photon.energy >= sourceRegion->params.energy_low_eV && photon.energy <= sourceRegion->params.energy_hi_eV);
		if (!validEnergy && photon.energy>0.0) {
			retries++;
			pointIdLocal = (int)(rnd()*sourceRegion->Points.size());
		}
	} while (!validEnergy && photon.energy>0.0 && retries < 5);

	if (!validEnergy && photon.energy>0.0) {
		char tmp[1024];
		sprintf(tmp, "Region %d point %d: can't generate within energy limits (%geV .. %geV)", regionId + 1, pointIdLocal + 1,
			sourceRegion->params.energy_low_eV , sourceRegion->params.energy_hi_eV);
		SetErrorSub(tmp);
		return false;
	}

	//sHandle->distTraveledCurrentParticle=0.0;
	sHandle->dF = photon.SR_flux;
	sHandle->dP = photon.SR_power;
	sHandle->energy = photon.energy;

	sHandle->oriRatio = 1.0;

	//starting position
	sHandle->pPos.x = photon.start_pos.x;
	sHandle->pPos.y = photon.start_pos.y;
	sHandle->pPos.z = photon.start_pos.z;

	sHandle->sourceRegionId = regionId;

	RecordHit(HIT_DES, sHandle->dF, sHandle->dP);

	//angle
	sHandle->pDir.x = photon.start_dir.x;
	sHandle->pDir.y = photon.start_dir.y;
	sHandle->pDir.z = photon.start_dir.z;

	//_ASSERTE(Norme(&sHandle->pDir)<=1.0);
	// Current structure = 0
	//sHandle->curStruct = src->sh.superIdx;
	//sHandle->curStruct=sourceRegion->i_struct1;
	sHandle->curStruct = 0;
	sHandle->teleportedFrom = -1;

	// Count
	sHandle->totalDesorbed++;
	sHandle->tmpCount.fluxAbs += sHandle->dF;
	sHandle->tmpCount.powerAbs += sHandle->dP;
	sHandle->tmpCount.nbDesorbed++;

	sHandle->lastHit = NULL;

	return true;

}

// -------------------------------------------------------------
// Compute bounce against a facet
// -------------------------------------------------------------
double TruncatedGaussian(gsl_rng *gen, const double &mean, const double &sigma, const double &lowerBound, const double &upperBound);

void PerformBounce_new(FACET *iFacet, const double &inTheta, const double &inPhi, const int &reflType) {

	double outTheta, outPhi; //perform bounce without scattering, will perturbate these angles later if it's a rough surface
	if (iFacet->sh.reflectType == REF_DIFFUSE) {
		outTheta = acos(sqrt(rnd()));
		outPhi = rnd()*2.0*PI;
	}
	else if (iFacet->sh.reflectType == REF_MIRROR) {
		outTheta = PI - inTheta;
		outPhi = inPhi;
	}
	else { //material reflection, might have backscattering
		switch (reflType) {
		case 1: //forward scattering
			outTheta = PI - inTheta;
			outPhi = inPhi;
			break;
		case 2: //diffuse scattering
			outTheta = acos(sqrt(rnd()));
			outPhi = rnd()*2.0*PI;
			break;
		case 3: //back scattering
			outTheta = PI - inTheta;
			outPhi = PI + inPhi;
			break;
		} //end switch (transparent pass treated at the Intersect() routine
	} //end material reflection

	if (iFacet->sh.doScattering) {
		double incidentAngle = abs(inTheta);
		if (incidentAngle > PI / 2) incidentAngle = PI - incidentAngle; //coming from the normal side
		double y = cos(incidentAngle);
		double wavelength = 3E8*6.626E-34 / (sHandle->energy*1.6E-19); //energy[eV] to wavelength[m]
		double specularReflProbability = exp(-Sqr(4 * PI*iFacet->sh.rmsRoughness*y / wavelength)); //Debye-Wallers factor, See "Measurements of x-ray scattering..." by Dugan, Sonnad, Cimino, Ishibashi, Scafers, eq.2
		bool specularReflection = rnd() < specularReflProbability;
		if (!specularReflection) {
			//Smooth surface reflection performed, now let's perturbate the angles
			//Using Gaussian approximated distributions of eq.14. of the above article
			double onePerTau = iFacet->sh.rmsRoughness / iFacet->sh.autoCorrLength;
			
			//Old acceptance-rejection algorithm
			/*
			size_t nbTries = 0; double outThetaPerturbated;
			do {
				double dTheta = Gaussian(2.9264*onePerTau); //Grazing angle perturbation, depends only on roughness, must assure it doesn't go against the surface
				outThetaPerturbated = outTheta + dTheta;
				nbTries++;
			} while (((outTheta < PI / 2) != (outThetaPerturbated < PI / 2)) && nbTries < 10);*/
			

			//New truncated Gaussian algorithm, see N. Chopin: Fast simulation of truncated Gaussian distributions, DOI: 10.1007/s11222-009-9168-1
			
			double lowerBound = 0.0;
			double upperBound = PI/2;
			if (outTheta > (PI / 2)) { //Limits: PI/2 .. PI instead of 0..PI/2
				lowerBound += PI / 2;
				upperBound += PI / 2;
			}
			double outThetaPerturbated = TruncatedGaussian(sHandle->gen, outTheta, 2.9264*onePerTau, lowerBound, upperBound);

			double dPhi = Gaussian((2.80657*pow(incidentAngle, -1.00238) - 1.00293*pow(incidentAngle, 1.22425))*onePerTau); //Out-of-plane angle perturbation, depends on roughness and incident angle
			outTheta = outThetaPerturbated;
			outPhi += dPhi;
		}
	}

	PolarToCartesian(iFacet, outTheta, outPhi, false);

	RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
	sHandle->lastHit = iFacet;
	if (iFacet->hits_MC && iFacet->sh.countRefl) RecordHitOnTexture(iFacet, sHandle->dF, sHandle->dP);
}

bool PerformBounce_old(FACET *iFacet, int reflType, double theta, double phi,
	Vector3d N_rotated, Vector3d nU_rotated, Vector3d nV_rotated) {

	double inPhi, inTheta;

	// Relaunch particle
	if (iFacet->sh.reflectType == REF_DIFFUSE) {
		//See docs/theta_gen.png for further details on angular distribution generation
		PolarToCartesian(iFacet, acos(sqrt(rnd())), rnd()*2.0*PI, false);
	} else { //Mirror reflection, optionally with surface perturbation
		if (iFacet->sh.reflectType == REF_MIRROR) reflType = REFL_FORWARD;
		if (!VerifiedSpecularReflection(iFacet, reflType, theta, phi,
			N_rotated, nU_rotated, nV_rotated)) {
			return false;
		}
	}
	RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
	sHandle->lastHit = iFacet;


	if (iFacet->hits_MC && iFacet->sh.countRefl) RecordHitOnTexture(iFacet, sHandle->dF, sHandle->dP);
	return true;
}

bool VerifiedSpecularReflection(FACET *iFacet, int reflType, double inTheta, double inPhi,
	Vector3d N_rotated, Vector3d nU_rotated, Vector3d nV_rotated) {
	//Specular reflection that returns false if going against surface
	double outTheta, outPhi;

	Vector3d N_facet = Vector3d(iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z);

	Vector3d newDir;
	double u, v, n;

	switch (reflType) {
	case 1: //forward scattering
		outTheta = PI - inTheta;
		outPhi = inPhi;
		break;
	case 2: //diffuse scattering
		outTheta = acos(sqrt(rnd()));
		outPhi = rnd()*2.0*PI;
		break;
	case 3: //back scattering
		outTheta = PI - inTheta;
		outPhi = PI + inPhi;
		break;
	} //end switch (transparent pass treated at the Intersect() routine

	u = sin(outTheta)*cos(outPhi);
	v = sin(outTheta)*sin(outPhi);
	n = cos(outTheta);

	newDir = Vector3d(u*nU_rotated.x + v*nV_rotated.x + n*N_rotated.x,
		u*nU_rotated.y + v*nV_rotated.y + n*N_rotated.y,
		u*nU_rotated.z + v*nV_rotated.z + n*N_rotated.z);

	if ((Dot(newDir,N_facet) > 0) != (inTheta > PI / 2)) {
		//inTheta > PI/2: ray coming from normal side
		//Dot(newDir,N_facet)>0: ray leaving towards normal side
		return false; //if reflection would go against the surface, generate new angles
	}

	sHandle->pDir.x = newDir.x;	
	sHandle->pDir.y = newDir.y;
	sHandle->pDir.z = newDir.z;
	return true;
}

void RecordHitOnTexture(FACET *f, double dF, double dP) {
	int tu = (int)(f->colU * f->sh.texWidthD);
	int tv = (int)(f->colV * f->sh.texHeightD);
	int add = tu + tv*f->sh.texWidth;
	f->hits_MC[add]++;
	f->hits_flux[add] += dF*f->inc[add]; //normalized by area
	f->hits_power[add] += dP*f->inc[add]; //normalized by area
}

void RecordDirectionVector(FACET *f) {
	int tu = (int)(f->colU * f->sh.texWidthD);
	int tv = (int)(f->colV * f->sh.texHeightD);
	int add = tu + tv*(f->sh.texWidth);

	f->direction[add].dir.x += sHandle->pDir.x;
	f->direction[add].dir.y += sHandle->pDir.y;
	f->direction[add].dir.z += sHandle->pDir.z;
	f->direction[add].count++;
}

void Stick(FACET* collidedFacet) {
	collidedFacet->counter.nbAbsorbed++;
	collidedFacet->counter.fluxAbs += sHandle->dF;
	collidedFacet->counter.powerAbs += sHandle->dP;
	sHandle->tmpCount.nbAbsorbed++;
	//sHandle->distTraveledSinceUpdate+=sHandle->distTraveledCurrentParticle;
	//sHandle->counter.nbAbsorbed++;
	//sHandle->counter.fluxAbs+=sHandle->dF;
	//sHandle->counter.powerAbs+=sHandle->dP;
	RecordHit(HIT_ABS, sHandle->dF, sHandle->dP); //for hits and lines display
	if (collidedFacet->hits_MC && collidedFacet->sh.countAbs) RecordHitOnTexture(collidedFacet, sHandle->dF, sHandle->dP);
}