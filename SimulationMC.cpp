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
#include "Tools.h"
#include "Utils.h" //debug only
#include <vector>
#include "Region_mathonly.h"
#include "GeneratePhoton.h"
#include "Utils.h" //rotation

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

void PolarToCartesian(FACET *iFacet, double theta, double phi, BOOL reverse, double rotateUV) {

	VERTEX3D U, V, N;
	double u, v, n;

	// Polar in (nU,nV,N) to Cartesian(x,y,z) transformation  ( nU = U/|U| , nV = V/|V| )
	// tetha is the angle to the normal of the facet N, phi to U
	// ! See Geometry::InitializeGeometry() for further informations on the (U,V,N) basis !
	// (nU,nV,N) and (x,y,z) are both left handed

	/*#ifdef WIN32
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
	double u = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.nU.x, iFacet->sh.nU.y, iFacet->sh.nU.z);
	double v = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.nV.x, iFacet->sh.nV.y, iFacet->sh.nV.z);
	double n = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z);
	SATURATE(n, -1.0, 1.0); //sometimes rounding errors do occur, 'acos' function would return no value for theta

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
	sHandle->lastUpdateOK = AccessDataportTimed(dpHit, timeout);
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
	nb = gHits->nbLastLeaks;
	for (i = 0; i < sHandle->nbLastLeak && i < NBHLEAK; i++)
		gHits->pLeak[(i + nb) % NBHLEAK] = sHandle->pLeak[i];
	gHits->nbLeakTotal += sHandle->nbLeakTotal;
	gHits->nbLastLeaks += sHandle->nbLastLeak;

	// HHit (Only prIdx 0)
	if (prIdx == 0) {
		gHits->nbHHit = sHandle->nbHHit;
		memcpy(gHits->pHits, sHandle->pHits, NBHHIT*sizeof(HIT));
	}

	// Facets
	for (s = 0; s < sHandle->nbSuper; s++) {
		for (i = 0; i < sHandle->str[s].nbFacet; i++) {

			FACET *f = sHandle->str[s].facets[i];
			if (f->hitted) {

				SHHITS *fFit = (SHHITS *)(buffer + f->sh.hitOffset);
				fFit->nbAbsorbed += f->sh.counter.nbAbsorbed;
				fFit->nbDesorbed += f->sh.counter.nbDesorbed;
				fFit->fluxAbs += f->sh.counter.fluxAbs;
				fFit->powerAbs += f->sh.counter.powerAbs;
				fFit->nbHit += f->sh.counter.nbHit;

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
					VHIT *shDir = (VHIT *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + f->profileSize + f->textureSize));
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
					double *shSpectrum_fluxwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + f->profileSize
						+ f->textureSize + f->directionSize));
					for (j = 0; j < SPECTRUM_SIZE; j++)
						shSpectrum_fluxwise[j] += f->spectrum_fluxwise->GetCount(j);

					double *shSpectrum_powerwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + f->profileSize
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

	//printf("\nResetCounter called from UpdateMCHits");
	ResetCounter();

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
	BOOL found = FALSE;
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
				found = TRUE;
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
	PolarToCartesian(destination, inTheta, inPhi, FALSE);
	// Move particle to teleport destination point
	sHandle->pPos.x = destination->sh.O.x + iFacet->colU*destination->sh.U.x + iFacet->colV*destination->sh.V.x;
	sHandle->pPos.y = destination->sh.O.y + iFacet->colU*destination->sh.U.y + iFacet->colV*destination->sh.V.y;
	sHandle->pPos.z = destination->sh.O.z + iFacet->colU*destination->sh.U.z + iFacet->colV*destination->sh.V.z;
	RecordHit(HIT_TELEPORT, sHandle->dF, sHandle->dP);
	sHandle->lastHit = destination;

	//Count hits on teleport facets (only TP source)
	//iFacet->sh.counter.nbAbsorbed++;
	//destination->sh.counter.nbDesorbed++;

	iFacet->sh.counter.nbHit++;//destination->sh.counter.nbHit++;
	iFacet->sh.counter.fluxAbs += sHandle->dF;//destination->sh.counter.fluxAbs+=sHandle->dF;
	iFacet->sh.counter.powerAbs += sHandle->dP;//destination->sh.counter.powerAbs+=sHandle->dP;
}

// -------------------------------------------------------------
// Perform nbStep simulation steps (a step is a bounce)
// -------------------------------------------------------------

BOOL SimulationMCStep(int nbStep) {

	FACET   *collidedFacet;
	double   d;
	BOOL     found;
	int      i;

	// Perform simulation steps
	for (i = 0; i < nbStep; i++) {

		found = Intersect(&(sHandle->pPos), &(sHandle->pDir), &d, &collidedFacet, sHandle->lastHit); //May decide reflection type

		if (found) {

			// Move particule to intersection point
			sHandle->pPos.x += d*sHandle->pDir.x;
			sHandle->pPos.y += d*sHandle->pDir.y;
			sHandle->pPos.z += d*sHandle->pDir.z;
			//sHandle->distTraveledCurrentParticle += d;
			sHandle->distTraveledSinceUpdate += d;

			if (!collidedFacet->sh.teleportDest){
				// Hit count
				sHandle->tmpCount.nbHit++;
				collidedFacet->sh.counter.nbHit++;
			}
			else {
				PerformTeleport(collidedFacet);
			}

			if (collidedFacet->sh.superDest) {	// Handle super structure link facet
				sHandle->curStruct = collidedFacet->sh.superDest - 1;
				// Count this hit as a transparent pass
				RecordHit(HIT_TRANS, sHandle->dF, sHandle->dP);
				if (collidedFacet->hits_MC && collidedFacet->sh.countTrans) RecordHitOnTexture(collidedFacet, sHandle->dF, sHandle->dP);
			}
			else if (collidedFacet->sh.isVolatile) { // Handle volatile facet
				if (collidedFacet->ready) {
					collidedFacet->sh.counter.nbAbsorbed++;
					collidedFacet->sh.counter.fluxAbs += sHandle->dF;
					collidedFacet->sh.counter.powerAbs += sHandle->dP;
					collidedFacet->ready = FALSE;
					if (collidedFacet->hits_MC && collidedFacet->sh.countAbs) RecordHitOnTexture(collidedFacet, sHandle->dF, sHandle->dP);
				}
			}
			else { //Just a regular facet, stick or bounce
				//Get incident angles
				double theta, phi;
				CartesianToPolar(collidedFacet, &theta, &phi);

				double stickingProbability, random;
				double nonTransparentProbability = 1.0;
				BOOL complexScattering = FALSE; //forward/diffuse/back/transparent/stick
				std::vector<double> materialReflProbabilities;
				if (sHandle->dF == 0.0 || sHandle->dP == 0.0 || sHandle->energy < 1E-3)
					stickingProbability = 1.0; //stick non real photons (from beam beginning)
				else if (collidedFacet->sh.reflectType < 2) //Diffuse or Mirror
					stickingProbability = collidedFacet->sh.sticking;
				else { //Material
					Material *mat = &(sHandle->materials[collidedFacet->sh.reflectType - 10]);
					complexScattering = mat->hasBackscattering;
					materialReflProbabilities = mat->Interpolate(sHandle->energy, abs(theta - PI / 2));

					//At this point we already know that it's not a transparent pass (Intersect() routine has ran)

					if (complexScattering) nonTransparentProbability -= materialReflProbabilities[3];

					stickingProbability = 1.0;
					size_t nbComponents = complexScattering ? 3 : 1;
					for (size_t i = 0; i < nbComponents; i++) //exclude transparent pass probability
						stickingProbability -= materialReflProbabilities[i];
				}
				if (!sHandle->lowFluxMode || complexScattering) {
					int reflType;
					//Do the original bounce or stick algorithm
					random = rnd()*nonTransparentProbability;
					if (random >(nonTransparentProbability - stickingProbability)) reflType = REFL_ABSORB;
					else if (!complexScattering) {
						reflType = REFL_FORWARD; //Not absorbed, so reflected
					}
					else { //complex scattering

						if (random < materialReflProbabilities[0]) reflType = REFL_FORWARD; //forward reflection
						else if (random < (materialReflProbabilities[0] + materialReflProbabilities[1])) reflType = REFL_DIFFUSE;
						else if (random < (materialReflProbabilities[0] + materialReflProbabilities[1] + materialReflProbabilities[2])) reflType = REFL_BACK;
						else reflType = REFL_BACK; //trasparent already excluded in Intersect() routine, and in this "else" branch we know that it doesn't stick
					}
					if (reflType == REFL_ABSORB) { if (!Stick(collidedFacet)) return FALSE; }
					else PerformBounce(collidedFacet, theta, phi, reflType);
				}
				else {
					//Low flux mode and simple scattering:
					//modified Stick() routine to record only the absorbed part	
					//No transparent pass in this case
					collidedFacet->sh.counter.fluxAbs += sHandle->dF*stickingProbability;
					collidedFacet->sh.counter.powerAbs += sHandle->dP*stickingProbability;
					if (collidedFacet->hits_MC && collidedFacet->sh.countAbs) RecordHitOnTexture(collidedFacet,
						sHandle->dF*stickingProbability, sHandle->dP*stickingProbability);
					//okay, absorbed part recorded, let's see how much is left
					sHandle->oriRatio *= 1.0 - stickingProbability;
					if (sHandle->oriRatio < sHandle->lowFluxCutoff) {//reflected part not important, throw it away
						RecordHit(HIT_ABS, sHandle->dF, sHandle->dP); //for hits and lines display
						if (!StartFromSource())
							// maxDesorption reached
							return FALSE; //FALSE
					}
					else { //reflect remainder
						sHandle->dF *= 1.0 - stickingProbability;
						sHandle->dP *= 1.0 - stickingProbability;
						PerformBounce(collidedFacet, theta, phi, REFL_FORWARD);
					}
				}
			} //end treating regular facet

			//Below is the old routine with the old scattering model
			//------------------------------------------------------

			//else if (collidedFacet->sh.reflectType >= 2) { //rough surface reflection
			//	if (sHandle->dF == 0.0 || sHandle->dP == 0.0 || sHandle->energy < 1E-3) { //stick non real photons (from beam beginning)
			//		if (!Stick(collidedFacet)) return FALSE;
			//	}
			//	else {
			//		double sigmaRatio = collidedFacet->sh.rmsRoughness;
			//		if ((collidedFacet->sh.reflectType - 2) < (int)sHandle->materials.size()) { //found material type
			//			//Generate incident angle
			//			Vector nU_facet = Vector(collidedFacet->sh.nU.x, collidedFacet->sh.nU.y, collidedFacet->sh.nU.z);
			//			Vector nV_facet = Vector(collidedFacet->sh.nV.x, collidedFacet->sh.nV.y, collidedFacet->sh.nV.z);
			//			Vector N_facet = Vector(collidedFacet->sh.N.x, collidedFacet->sh.N.y, collidedFacet->sh.N.z);
			//			Vector nU_rotated, N_rotated, nV_rotated;
			//			double n, u, v, thetaOffset, phiOffset, randomAngle;
			//			BOOL reflected = FALSE;
			//			do { //generate angles until reflecting away from surface
			//				double n_ori = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
			//					N_facet.x, N_facet.y, N_facet.z); //original incident angle (negative if front collision, positive if back collision);
			//				do { //generate angles until incidence is from front
			//					//Generating bending angle due to local roughness
			//					double rnd1 = rnd();
			//					double rnd2 = rnd(); //for debug
			//					thetaOffset = atan(sigmaRatio*tan(PI*(rnd1 - 0.5)));
			//					phiOffset = atan(sigmaRatio*tan(PI*(rnd2 - 0.5)));

			//					//Random rotation around N (to discard U orientation thus make scattering isotropic)
			//					randomAngle = rnd() * 2 * PI;
			//					nU_facet = nU_facet.Rotate(N_facet, randomAngle);
			//					nV_facet = nV_facet.Rotate(N_facet, randomAngle);

			//					//Bending surface base vectors
			//					nU_rotated = nU_facet.Rotate(nV_facet, thetaOffset);
			//					N_rotated = N_facet.Rotate(nV_facet, thetaOffset);
			//					nU_rotated = nU_rotated.Rotate(nU_facet, phiOffset); //check if correct
			//					nV_rotated = nV_facet.Rotate(nU_facet, phiOffset);
			//					N_rotated = N_rotated.Rotate(nU_facet, phiOffset);

			//					//Calculate incident angles to bent surface (Cartesian to Polar routine modified)
			//					u = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
			//						nU_rotated.x, nU_rotated.y, nU_rotated.z);
			//					v = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
			//						nV_rotated.x, nV_rotated.y, nV_rotated.z);
			//					n = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
			//						N_rotated.x, N_rotated.y, N_rotated.z);
			//					SATURATE(n, -1.0, 1.0); //sometimes rounding errors do occur, 'acos' function would return no value for theta

			//				} while (n*n_ori <= 0.0); //regenerate random numbers if grazing angle would be over 90deg

			//				double rho = sqrt(v*v + u*u);
			//				double theta = acos(n);              // Angle to normal (PI/2 => PI if front collision, 0..PI/2 if back)
			//				double phi = asin(v / rho);
			//				if (u < 0.0) phi = PI - phi;  // Angle to U

			//				if (!sHandle->lowFluxMode || sHandle->materials[collidedFacet->sh.reflectType - 2].hasBackscattering) { //original algorithm, bounce or stick
			//					int reflType = sHandle->materials[collidedFacet->sh.reflectType - 2].GetReflectionType(sHandle->energy, abs(theta - PI / 2),rnd());
			//					if (reflType) { //not absorbed
			//						reflected = PerformBounce(collidedFacet, sigmaRatio, theta, phi, N_rotated, nU_rotated, nV_rotated,thetaOffset,phiOffset,randomAngle, reflType);
			//					}
			//					else { //sticking
			//						if (!Stick(collidedFacet)) return FALSE;
			//						reflected = TRUE;
			//					}
			//				}
			//				else { //low flux mode
			//					//modified Stick() routine to record only the absorbed part
			//					double reflProbability = sHandle->materials[collidedFacet->sh.reflectType - 2].Interpolate(sHandle->energy, abs(theta - PI / 2))[0];	//first component							
			//					collidedFacet->sh.counter.fluxAbs += sHandle->dF*(1.0 - reflProbability);
			//					collidedFacet->sh.counter.powerAbs += sHandle->dP*(1.0 - reflProbability);
			//					//sHandle->distTraveledSinceUpdate += sHandle->distTraveledCurrentParticle;
			//					if (collidedFacet->hits_MC && collidedFacet->sh.countAbs) RecordHitOnTexture(collidedFacet,
			//						sHandle->dF*(1.0 - reflProbability), sHandle->dP*(1.0 - reflProbability));
			//					//okay, absorbed part recorded, let's see how much is left
			//					sHandle->oriRatio *= reflProbability;
			//					if (sHandle->oriRatio < sHandle->lowFluxCutoff) {//reflected part not important, throw it away
			//						RecordHit(HIT_ABS, sHandle->dF, sHandle->dP); //for hits and lines display
			//						if (!StartFromSource())
			//							// maxDesorption reached
			//							return FALSE; //FALSE
			//						reflected = TRUE;
			//					}
			//					else { //reflect remainder
			//						sHandle->dF *= reflProbability;
			//						sHandle->dP *= reflProbability;
			//						reflected = PerformBounce(collidedFacet, sigmaRatio, theta, phi, N_rotated, nU_rotated, nV_rotated, thetaOffset, phiOffset, randomAngle);
			//					}
			//				}
			//			} while (!reflected); //do it again if reflection wasn't successful (reflected against the surface due to roughness)
			//		}
			//		else {
			//			std::string err = "Facet " + (collidedFacet->globalId + 1);
			//			err += "reflection material type not found\0";
			//			SetErrorSub(err.c_str());
			//		}

			//	}
			//}
			//else if (collidedFacet->sh.sticking > 0.0) {
			//	if (!sHandle->lowFluxMode) {//original algorithm, bounce or stick
			//		if (collidedFacet->sh.sticking == 1.0 || rnd() < collidedFacet->sh.sticking) { //sticking
			//			if (!Stick(collidedFacet)) return FALSE;
			//		}
			//		else {
			//			PerformBounce(collidedFacet);
			//		}
			//	}
			//	else { //low flux mode
			//		{ //low flux mode
			//			//modified Stick() routine to record only the absorbed part
			//			collidedFacet->sh.counter.fluxAbs += sHandle->dF*collidedFacet->sh.sticking;
			//			collidedFacet->sh.counter.powerAbs += sHandle->dP*collidedFacet->sh.sticking;
			//			//sHandle->distTraveledSinceUpdate += sHandle->distTraveledCurrentParticle;
			//			if (collidedFacet->hits_MC && collidedFacet->sh.countAbs) RecordHitOnTexture(collidedFacet,
			//				sHandle->dF*collidedFacet->sh.sticking, sHandle->dP*collidedFacet->sh.sticking);
			//			//okay, absorbed part recorded, let's see how much is left
			//			sHandle->oriRatio *= (1.0 - collidedFacet->sh.sticking);
			//			if (sHandle->oriRatio < sHandle->lowFluxCutoff) {//reflected part not important, throw it away
			//				RecordHit(HIT_ABS, sHandle->dF, sHandle->dP); //for hits and lines display
			//				if (!StartFromSource())
			//					// maxDesorption reached
			//					return FALSE; //FALSE
			//			}
			//			else { //reflect remainder
			//				sHandle->dF *= 1.0 - collidedFacet->sh.sticking;
			//				sHandle->dP *= 1.0 - collidedFacet->sh.sticking;
			//				PerformBounce(collidedFacet);
			//			}
			//		}
			//	}
			//}
			//else {
			//	PerformBounce(collidedFacet); //if sticking==0, bounce without generating random number
			//}

		}
		else { // Leak (simulation error)
			RecordLeakPos();
			sHandle->nbLeakTotal++;
			if (!StartFromSource())
				// maxDesorption reached
				return FALSE;
		} //end intersect or leak
	} //end step
	return TRUE;
}

// -------------------------------------------------------------
// Launch photon from a trajectory point
// -------------------------------------------------------------

BOOL StartFromSource() {

	// Check end of simulation
	if (sHandle->maxDesorption > 0) {
		if (sHandle->totalDesorbed >= sHandle->maxDesorption) {
			sHandle->lastHit = NULL;
			return FALSE;
		}
	}

	//find source point
	int pointIdGlobal = (int)(rnd()*sHandle->sourceArea);
	BOOL found = false;
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
	if (!(sourceRegion->psimaxX > 0.0 && sourceRegion->psimaxY>0.0)) SetErrorSub("psiMaxX or psiMaxY not positive. No photon can be generated");
	GenPhoton photon = GeneratePhoton(pointIdLocal, sourceRegion, sHandle->generation_mode, sHandle->psi_distr, sHandle->chi_distr, sHandle->tmpCount.nbDesorbed == 0);

	//sHandle->distTraveledCurrentParticle=0.0;
	sHandle->dF = photon.SR_flux;
	sHandle->dP = photon.SR_power;
	sHandle->energy = photon.energy;

	sHandle->oriRatio = 1.0;

	/*Vector Z_local=source->direction.Normalize(); //Z' base vector
	Vector X_local=source->rho.Normalize(); //X' base vector
	Vector Y_local=Crossproduct(Z_local,X_local);*/

	/*Vector start_pos=source->position;
	start_pos=Add(start_pos,ScalarMult(X_local,photon.dX)); //apply dX offset
	start_pos=Add(start_pos,ScalarMult(Y_local,photon.dY)); //apply dY offset*/

	//starting position
	sHandle->pPos.x = photon.start_pos.x;
	sHandle->pPos.y = photon.start_pos.y;
	sHandle->pPos.z = photon.start_pos.z;

	RecordHit(HIT_DES, sHandle->dF, sHandle->dP);

	//angle
	/*Vector start_dir=Z_local;
	start_dir=start_dir.Rotate(Y_local,photon.divX);
	start_dir=start_dir.Rotate(X_local,photon.divY);*/

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

	return TRUE;

}

// -------------------------------------------------------------
// Compute bounce against a facet
// -------------------------------------------------------------

void PerformBounce(FACET *iFacet, const double &inTheta, const double &inPhi, const int &reflType) {

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
			outTheta = PI + inTheta;
			outPhi = inPhi;
			break;
		} //end switch (transparent pass treated at the Intersect() routine
	} //end material reflection

	if (iFacet->sh.doScattering) {
		double incidentAngle = abs(inTheta);
		if (incidentAngle > PI / 2) incidentAngle = PI - incidentAngle; //coming from the normal side
		double y = cos(incidentAngle);
		double wavelength = 3E8*6.626E-34 / (sHandle->energy*1.6E-19); //energy[eV] to wavelength[m]
		double specularReflProbability = exp(-Sqr(4 * PI*iFacet->sh.rmsRoughness*y / wavelength)); //Debye-Wallers factor, See "Measurements of x-ray scattering..." by Dugan, Sonnad, Cimino, Ishibashi, Scafers, eq.2
		BOOL specularReflection = rnd() < specularReflProbability;
		if (!specularReflection) {
			//Smooth surface reflection performed, now let's perturbate the angles
			//Using Gaussian approximated distributions of eq.14. of the above article
			double onePerTau = iFacet->sh.rmsRoughness / iFacet->sh.autoCorrLength;

			size_t nbTries = 0; double outThetaPerturbated;
			do {
				double dTheta = Gaussian(2.9264*onePerTau); //Grazing angle perturbation, depends only on roughness, must assure it doesn't go against the surface
				outThetaPerturbated = outTheta + dTheta;
				nbTries++;
			} while (((outTheta < PI / 2) != (outThetaPerturbated < PI / 2)) && nbTries < 10);
			double dPhi = Gaussian((2.80657*pow(incidentAngle, -1.00238) - 1.00293*pow(incidentAngle, 1.22425))*onePerTau); //Out-of-plane angle perturbation, depends on roughness and incident angle
			outTheta = outThetaPerturbated;
			outPhi += dPhi;
		}
	}

	PolarToCartesian(iFacet, outTheta, outPhi, FALSE);

	RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
	sHandle->lastHit = iFacet;
	if (iFacet->hits_MC && iFacet->sh.countRefl) RecordHitOnTexture(iFacet, sHandle->dF, sHandle->dP);
}

//Old PerformBounce Routine
//int PerformBounce(FACET *iFacet, double sigmaRatio, double theta, double phi,
//	Vector N_rotated, Vector nU_rotated, Vector nV_rotated,  double thetaOffset, double phiOffset, double randomAngle,int reflType) {
//
//	double inPhi, inTheta;
//	//BOOL revert=FALSE;
//
//	/*if( iFacet->sh.is2sided ) {
//		// We may need to revert normal in case of 2 sided hit
//		revert = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
//		iFacet->sh.N.x,iFacet->sh.N.y,iFacet->sh.N.z)>0.0;
//		} //on the other hand, theta alrady contains this information*/
//
//	// 0<theta<PI/2:  back reflection
//	// PI/2<theta<PI: front reflection
//
//	// Relaunch particle
//	switch (iFacet->sh.reflectType) {
//	case REF_MIRROR:
//		CartesianToPolar(iFacet, &inTheta, &inPhi);
//		PolarToCartesian(iFacet, PI - inTheta, inPhi, FALSE);
//		break;
//	case REF_DIFFUSE:
//		//See docs/theta_gen.png for further details on angular distribution generation
//		PolarToCartesian(iFacet, acos(sqrt(rnd())), rnd()*2.0*PI, FALSE);
//		//To do: check if against facet
//		break;
//	default: //Material reflection - theta, etc. are already calculated
//		switch (reflType) {
//		case 1: //forward reflection
//		{
//			Vector N_facet = Vector(iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z);
//
//			double u, v, n;
//
//			theta = PI - theta; //perform reflection
//
//			u = sin(theta)*cos(phi);
//			v = sin(theta)*sin(phi);
//			n = cos(theta);
//
//			Vector newDir = Vector(u*nU_rotated.x + v*nV_rotated.x + n*N_rotated.x,
//				u*nU_rotated.y + v*nV_rotated.y + n*N_rotated.y,
//				u*nU_rotated.z + v*nV_rotated.z + n*N_rotated.z);
//
//			if ((DOT3(newDir.x, newDir.y, newDir.z,
//				N_facet.x, N_facet.y, N_facet.z) > 0) != (theta < PI / 2)) {
//				return 0; //if reflection would go against the surface, generate new angles
//			}
//
//			sHandle->pDir.x = newDir.x;
//			sHandle->pDir.y = newDir.y;
//			sHandle->pDir.z = newDir.z;
//
//			//CartesianToPolar(iFacet, &inTheta, &inPhi);
//			//PolarToCartesian(iFacet, PI - inTheta+2.0*thetaOffset, inPhi+2.0*phiOffset, FALSE,TRUE);
//
//			break;
//		}
//		case 2: //diffuse reflection
//			PolarToCartesian(iFacet, acos(sqrt(rnd())), rnd()*2.0*PI, FALSE);
//			break;
//
//		case 3: //back reflection
//		{
//			/*Vector N_facet = Vector(iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z);
//
//			Vector newDir;
//			double u, v, n;
//
//			theta = PI + theta; //turn back ray
//
//			u = sin(theta)*cos(phi);
//			v = sin(theta)*sin(phi);
//			n = cos(theta);
//
//			newDir = Vector(u*nU_rotated.x + v*nV_rotated.x + n*N_facet.x,
//				u*nU_rotated.y + v*nV_rotated.y + n*N_facet.y,
//				u*nU_rotated.z + v*nV_rotated.z + n*N_facet.z);
//
//			if ((DOT3(newDir.x, newDir.y, newDir.z,
//				N_facet.x, N_facet.y, N_facet.z) > 0) != (theta < 2*PI)) {
//				return 0; //if reflection would go against the surface, generate new angles
//			}*/
//
//			Vector backDir = Vector(-sHandle->pDir.x, -sHandle->pDir.y, -sHandle->pDir.z); //Turn back ray
//			Vector newDir = backDir.Rotate(Crossproduct(backDir,Vector(iFacet->sh.N.x,iFacet->sh.N.y,iFacet->sh.N.z)), 5.0*thetaOffset); //Random perturbation, coeff. 5.0 empirically found
//			newDir=newDir.Rotate(backDir, randomAngle);
//
//			sHandle->pDir.x = newDir.x;
//			sHandle->pDir.y = newDir.y;
//			sHandle->pDir.z = newDir.z;
//			break;
//		}
//		}
//		break;
//	}
//	RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
//	sHandle->lastHit = iFacet;
//
//
//	if (iFacet->hits_MC && iFacet->sh.countRefl) RecordHitOnTexture(iFacet, sHandle->dF, sHandle->dP);
//	return 1;
//}

void RecordHitOnTexture(FACET *f, double dF, double dP) {
	int tu = (int)(f->colU * f->sh.texWidthD);
	int tv = (int)(f->colV * f->sh.texHeightD);
	int add = tu + tv*f->sh.texWidth;
	f->hits_MC[add]++;
	f->hits_flux[add] += dF*f->inc[add]; //normalized by area
	f->hits_power[add] += dP*f->inc[add]; //normalized by area
}

double Gaussian(const double &sigma) {
	//Box-Muller transform
	//return sigma*sqrt(-2 * log(rnd()))*cos(2 * PI*rnd());

	//Generates a random number following the Gaussian distribution around 0 with 'sigma' standard deviation
	double v1, v2, r, fac;
	do {
		v1 = 2.0*rnd() - 1.0;
		v2 = 2.0*rnd() - 1.0;
		r = Sqr(v1) + Sqr(v2);
	} while (r >= 1.0);
	fac = sqrt(-2.0*log(r) / r);
	return v2*fac*sigma;
}

int Stick(FACET* collidedFacet) {
	collidedFacet->sh.counter.nbAbsorbed++;
	collidedFacet->sh.counter.fluxAbs += sHandle->dF;
	collidedFacet->sh.counter.powerAbs += sHandle->dP;
	sHandle->tmpCount.nbAbsorbed++;
	//sHandle->distTraveledSinceUpdate+=sHandle->distTraveledCurrentParticle;
	//sHandle->counter.nbAbsorbed++;
	//sHandle->counter.fluxAbs+=sHandle->dF;
	//sHandle->counter.powerAbs+=sHandle->dP;
	RecordHit(HIT_ABS, sHandle->dF, sHandle->dP); //for hits and lines display
	if (collidedFacet->hits_MC && collidedFacet->sh.countAbs) RecordHitOnTexture(collidedFacet, sHandle->dF, sHandle->dP);
	if (!StartFromSource())
		// maxDesorption reached
		return 0; //FALSE
	return 1;
}