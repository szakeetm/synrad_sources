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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include "Random.h"
#include "SynradDistributions.h"
#include <vector>
#include "Region_mathonly.h"
#include "GeneratePhoton.h"
#include "GLApp/MathTools.h"
#include "SynradTypes.h" //Histogram
#include <tuple>

extern Simulation *sHandle;
extern void SetErrorSub(const char *message);

//extern Distribution2D K_1_3_distribution;
//extern Distribution2D K_2_3_distribution;
//extern DistributionND SR_spectrum_CDF;
extern Distribution2D integral_N_photons;
extern Distribution2D integral_SR_power;
//extern Distribution2D polarization_distribution;
//extern Distribution2D g1h2_distribution;

void ComputeSourceArea() {
	sHandle->sourceArea = sHandle->nbTrajPoints;
}

void UpdateMCHits(Dataport *dpHit, int prIdx, DWORD timeout) {

	BYTE *buffer;
	GlobalHitBuffer *gHits;
	TextureCell oldMin, oldMax;
	size_t i, j, s, x, y;
#ifdef _DEBUG
	double t0, t1;
	t0 = GetTick();
#endif
	SetState(NULL, "Waiting for 'hits' dataport access...", false, true);
	sHandle->lastHitUpdateOK = AccessDataportTimed(dpHit, timeout);
	SetState(NULL, "Updating MC hits...", false, true);
	if (!sHandle->lastHitUpdateOK) return;

	buffer = (BYTE*)dpHit->buff;
	gHits = (GlobalHitBuffer *)buffer;

	// Global hits and leaks
	gHits->total.nbMCHit += sHandle->tmpGlobalCount.nbMCHit;
	gHits->total.nbHitEquiv += sHandle->tmpGlobalCount.nbHitEquiv;
	gHits->total.nbAbsEquiv += sHandle->tmpGlobalCount.nbAbsEquiv;
	gHits->total.nbDesorbed += sHandle->tmpGlobalCount.nbDesorbed;
	gHits->distTraveledTotal += sHandle->distTraveledSinceUpdate;
	gHits->total.fluxAbs += sHandle->tmpGlobalCount.fluxAbs;
	gHits->total.powerAbs += sHandle->tmpGlobalCount.powerAbs;

	oldMin = gHits->hitMin;
	oldMax = gHits->hitMax;

	gHits->hitMin.count = HITMAX_INT64;
	gHits->hitMax.count = 0;
	gHits->hitMin.flux = HITMAX_DOUBLE;
	gHits->hitMax.flux = 0;
	gHits->hitMin.power = HITMAX_DOUBLE;
	gHits->hitMax.power = 0;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak
	for (size_t leakIndex = 0; leakIndex < sHandle->leakCacheSize; leakIndex++)
		gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = sHandle->leakCache[leakIndex];
	gHits->nbLeakTotal += sHandle->nbLeakSinceUpdate;
	gHits->lastLeakIndex = (gHits->lastLeakIndex + sHandle->leakCacheSize) % LEAKCACHESIZE;
	gHits->leakCacheSize = Min(LEAKCACHESIZE, gHits->leakCacheSize + sHandle->leakCacheSize);

	// Hit cache (Only prIdx 0)
	if (prIdx == 0) {
		for (size_t hitIndex = 0; hitIndex < sHandle->hitCacheSize; hitIndex++)
			gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = sHandle->hitCache[hitIndex];

		if (sHandle->hitCacheSize > 0) {
			gHits->lastHitIndex = (gHits->lastHitIndex + sHandle->hitCacheSize) % HITCACHESIZE;

			//if (gHits->lastHitIndex < (HITCACHESIZE - 1)) {
			//	gHits->lastHitIndex++;
				gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
			//}

			gHits->hitCacheSize = Min(HITCACHESIZE, gHits->hitCacheSize + sHandle->hitCacheSize);
		}
	}

	// Facets
	for (s = 0; s < sHandle->nbSuper; s++) {
		for (i = 0; i < sHandle->str[s].nbFacet; i++) {

			SubprocessFacet *f = sHandle->str[s].facets[i];
			if (f->hitted) {

				FacetHitBuffer *fFit = (FacetHitBuffer *)(buffer + f->sh.hitOffset);
				*fFit += f->counter;

				if (f->sh.isProfile) {
					ProfileSlice *shProfile = (ProfileSlice *)(buffer + (f->sh.hitOffset + sizeof(FacetHitBuffer)));
					for (j = 0; j < PROFILE_SIZE; j++) {
						shProfile[j] += f->profile[j];
					}
				}

				size_t profileSize = (f->sh.isProfile) ? PROFILE_SIZE*sizeof(ProfileSlice) : 0;

				if (f->sh.isTextured) {
					size_t nbTextureCells = f->sh.texHeight*f->sh.texWidth;
					TextureCell *shTexture = (TextureCell *)(buffer + (f->sh.hitOffset + sizeof(FacetHitBuffer) + profileSize));
					for (y = 0; y < f->sh.texHeight; y++) {
						for (x = 0; x < f->sh.texWidth; x++) {
							size_t index = x + y*f->sh.texWidth;
							//Increase value
							shTexture[index] += f->texture[index];
							//Adjust min/max
							if (shTexture[index].count > gHits->hitMax.count)	gHits->hitMax.count = shTexture[index].count;
							if (shTexture[index].count > 0 && shTexture[index].count < gHits->hitMin.count) gHits->hitMin.count = shTexture[index].count;
							if (shTexture[index].flux > gHits->hitMax.flux && f->largeEnough[index]) gHits->hitMax.flux = shTexture[index].flux;
							if (shTexture[index].flux > 0.0 && shTexture[index].flux < gHits->hitMin.flux && f->largeEnough[index]) gHits->hitMin.flux = shTexture[index].flux;
							if (shTexture[index].power > gHits->hitMax.power && f->largeEnough[index]) gHits->hitMax.power = shTexture[index].power;
							if (shTexture[index].power > 0.0 && shTexture[index].power < gHits->hitMin.power && f->largeEnough[index]) gHits->hitMin.power = shTexture[index].power;
						}
					}
				}

				if (f->sh.countDirection) {
					DirectionCell *shDir = (DirectionCell *)(buffer + (f->sh.hitOffset + sizeof(FacetHitBuffer) + profileSize + f->textureSize));
					for (y = 0; y < f->sh.texHeight; y++) {
						for (x = 0; x < f->sh.texWidth; x++) {
							size_t add = x + y*f->sh.texWidth;
							shDir[add].dir += f->direction[add].dir;
							shDir[add].count += f->direction[add].count;
						}
					}
				}

				if (f->sh.recordSpectrum) {
					ProfileSlice *shSpectrum = (ProfileSlice *)(buffer + (f->sh.hitOffset + sizeof(FacetHitBuffer) + profileSize
						+ f->textureSize + f->directionSize));
					for (j = 0; j < SPECTRUM_SIZE; j++) {
						shSpectrum[j] += f->spectrum->GetCounts(j);
					}
				}

			} // End if(hitted)
		} // End nbFacet
	} // End nbSuper

	//if there were no textures:
	if (gHits->hitMin.count == HITMAX_INT64) gHits->hitMin.count = oldMin.count;
	if (gHits->hitMax.count == 0) gHits->hitMax.count = oldMax.count;
	if (gHits->hitMin.flux == HITMAX_DOUBLE) gHits->hitMin.flux = oldMin.flux;
	if (gHits->hitMax.flux == 0.0) gHits->hitMax.flux = oldMax.flux;
	if (gHits->hitMin.power == HITMAX_DOUBLE) gHits->hitMin.power = oldMin.power;
	if (gHits->hitMax.power == 0.0) gHits->hitMax.power = oldMax.power;

	ReleaseDataport(dpHit);
	ResetTmpCounters();
	extern char* GetSimuStatus();
	SetState(NULL, GetSimuStatus(), false, true);

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif

}

// Compute particle teleport

void PerformTeleport(SubprocessFacet& collidedFacet) {

	//Search destination
	SubprocessFacet *destination;
	bool found = false;
	int destIndex;

	if (collidedFacet.sh.teleportDest == -1) {
		destIndex = sHandle->teleportedFrom;
		if (destIndex == -1) {
			/*char err[128];
			sprintf(err, "Facet %d tried to teleport to the facet where the particle came from, but there is no such facet.", collidedFacet.globalId + 1);
			SetErrorSub(err);*/
			RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
			sHandle->lastHitFacet = &collidedFacet;
			return; //LEAK
		}
	}
	else destIndex = collidedFacet.sh.teleportDest - 1;

	//Look in which superstructure is the destination facet:
	for (size_t i = 0; i < sHandle->nbSuper && (!found); i++) {
		for (size_t j = 0; j < sHandle->str[i].nbFacet && (!found); j++) {
			if (destIndex == sHandle->str[i].facets[j]->globalId) {
				destination = sHandle->str[i].facets[j];
				sHandle->curStruct = (int)destination->sh.superIdx; //change current superstructure
				sHandle->teleportedFrom = (int)collidedFacet.globalId; //memorize where the particle came from
				found = true;
			}
		}
	}
	if (!found) {
		/*char err[128];
		sprintf(err, "Teleport destination of facet %d not found (facet %d does not exist)", collidedFacet.globalId + 1, collidedFacet.sh.teleportDest);
		SetErrorSub(err);*/
		RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
		sHandle->lastHitFacet = &collidedFacet;
		return; //LEAK
	}

	// Count this hit as a transparent pass
	RecordHit(HIT_TELEPORTSOURCE, sHandle->dF, sHandle->dP);
	if (collidedFacet.texture && collidedFacet.sh.countTrans) RecordHitOnTexture(collidedFacet, sHandle->dF, sHandle->dP);

	// Relaunch particle from new facet
	double inPhi, inTheta;
	std::tie(inTheta,inPhi) = CartesianToPolar(collidedFacet.sh.nU,collidedFacet.sh.nV, collidedFacet.sh.N);
	PolarToCartesian(destination, inTheta, inPhi, false);
	// Move particle to teleport destination point
	sHandle->pPos = destination->sh.O + collidedFacet.colU*destination->sh.U + collidedFacet.colV*destination->sh.V;

	RecordHit(HIT_TELEPORTDEST, sHandle->dF, sHandle->dP);
	sHandle->lastHitFacet = destination;

	//Count hits on teleport facets (only TP source)
	//collidedFacet.counter.nbAbsorbed++;
	//destination->counter.nbDesorbed++;

	ProfileSlice increment;
	increment.count_absorbed = 0;
	increment.count_incident = 1;
	increment.flux_absorbed = 0.0;
	increment.flux_incident = sHandle->dF;
	increment.power_absorbed = 0.0;
	increment.power_incident = sHandle->dP;
	ProfileFacet(collidedFacet, sHandle->energy, increment); //Put here since removed from Intersect() routine

	collidedFacet.counter.nbMCHit++;//destination->counter.nbMCHit++;
	collidedFacet.counter.nbHitEquiv += sHandle->oriRatio;
	collidedFacet.counter.fluxAbs += sHandle->dF;//destination->counter.fluxAbs+=sHandle->dF;
	collidedFacet.counter.powerAbs += sHandle->dP;//destination->counter.powerAbs+=sHandle->dP;
}

// Perform nbStep simulation steps (a step is a bounce)

bool SimulationMCStep(const size_t& nbStep) {

	// Perform simulation steps
	for (size_t i = 0; i < nbStep; i++) {

		SubprocessFacet*   collidedFacetPtr=NULL;
		double   d;
		bool     found;

		std::tie(found,collidedFacetPtr,d) = Intersect(sHandle->pPos, sHandle->pDir,THitCache); //May decide reflection type

		if (found) {
			
			SubprocessFacet& collidedFacet = *collidedFacetPtr; //better readability and passing as reference argument during this function

			// Move particle to intersection point
			sHandle->pPos = sHandle->pPos + d*sHandle->pDir;
			sHandle->distTraveledSinceUpdate += d;
			LogHit(&collidedFacet);

			if (collidedFacet.sh.teleportDest) {
				PerformTeleport(collidedFacet); //increases hit, flux, power counters
			}
			else {
				// Record incident
				sHandle->tmpGlobalCount.nbMCHit++;
				sHandle->tmpGlobalCount.nbHitEquiv += sHandle->oriRatio;
				collidedFacet.counter.nbMCHit++;
				collidedFacet.counter.nbHitEquiv += sHandle->oriRatio;
				if (collidedFacet.sh.superDest) {	// Handle super structure link facet
					sHandle->curStruct = collidedFacet.sh.superDest - 1;
					// Count this hit as a transparent pass
					RecordHit(HIT_TRANS, sHandle->dF, sHandle->dP);
					ProfileSlice increment;
					increment.count_absorbed = 0;
					increment.count_incident = 1;
					increment.flux_absorbed = 0.0;
					increment.flux_incident = sHandle->dF;
					increment.power_absorbed = 0.0;
					increment.power_incident = sHandle->dP;
					ProfileFacet(collidedFacet, sHandle->energy, increment);
					
					if (collidedFacet.texture && collidedFacet.sh.countTrans) RecordHitOnTexture(collidedFacet, sHandle->dF, sHandle->dP);
				}
				else { //Not superDest or Teleport
					if (sHandle->newReflectionModel) {
						//New reflection model (Synrad 1.4)
						double inTheta, inPhi;
						std::tie(inTheta,inPhi) = CartesianToPolar(collidedFacet.sh.nU,collidedFacet.sh.nV,collidedFacet.sh.N);

						double stickingProbability;
						bool complexScattering; //forward/diffuse/back/transparent/stick
						std::vector<double> materialReflProbabilities; //0: forward reflection, 1: diffuse, 2: backscattering, 3: transparent, 100%-(0+1+2+3): absorption

						if (sHandle->dF == 0.0 || sHandle->dP == 0.0 || sHandle->energy < 1E-3) {
							//stick non real photons (from beam beginning)
							stickingProbability = 1.0; 
							complexScattering = false;
						}
						else if (collidedFacet.sh.reflectType < 2) {
							//Diffuse or Mirror
							stickingProbability = collidedFacet.sh.sticking;
							complexScattering = false;
						}
						else { //Material
							std::tie(stickingProbability,materialReflProbabilities, complexScattering) = GetStickingProbability(collidedFacet, inTheta);
							//In the new reflection model, reflection probabilities don't take into account surface roughness
						}
						
						if (!sHandle->ontheflyParams.lowFluxMode) {
							//Regular mode, stick or bounce
							int reflType = GetHardHitType(stickingProbability, materialReflProbabilities, complexScattering);
							if (reflType == REFL_ABSORB) {
								Stick(collidedFacet);
								if (!StartFromSource()) return false;
							}
							else PerformBounce_new(collidedFacet, reflType, inTheta, inPhi);
						}
						else {
							//Low flux mode and simple scattering:
							Vector3d dummyNullVector(0.0, 0.0, 0.0); //DoLowFluxReflection will not use it since sHandle->newReflectionModel == true
							DoLowFluxReflection(collidedFacet, stickingProbability, complexScattering, materialReflProbabilities,
								inTheta, inPhi, dummyNullVector, dummyNullVector, dummyNullVector);
						} //end low flux mode
					}
					else {
						//Old reflection model (Synrad <=1.3 or if user deselects new model in Global Settings)						
						if (sHandle->dF == 0.0 || sHandle->dP == 0.0 || sHandle->energy < 1E-3) { //stick non real photons (from beam beginning)
							Stick(collidedFacet);
							if (!StartFromSource()) return false;
						}
						else {
							//We first generate a perturbated surface, and from that point on we treat the whole reflection process as if it were locally flat
							

							if ((collidedFacet.sh.reflectType - 10) < (int)sHandle->materials.size()) { //found material type
																										//Generate incident angle
								Vector3d nU_rotated, N_rotated, nV_rotated;
								bool reflected = false;
								do { //generate surfaces until reflected ray goes away from facet (and not through it)
									
									//First step: generate a random surface and determine incident angles
									double inTheta, inPhi;
									if (collidedFacet.sh.doScattering) {
										double n_ori = Dot(sHandle->pDir, collidedFacet.sh.N); //incident angle with original facet surface (negative if front collision, positive if back collision);
										double n_new;
										do { //generate angles until incidence is from front
											double sigmaRatio = collidedFacet.sh.rmsRoughness/collidedFacet.sh.autoCorrLength;
											std::tie(nU_rotated, nV_rotated, N_rotated) = PerturbateSurface(collidedFacet, sigmaRatio);
											n_new = Dot(sHandle->pDir, N_rotated);
										} while (n_new*n_ori <= 0.0); //generate new random surface if grazing angle would be over 90deg (shadowing)
										std::tie(inTheta, inPhi) = CartesianToPolar(nU_rotated, nV_rotated, N_rotated); //Finally, get incident angles
									}
									else { //no scattering, use original surface
										nU_rotated = collidedFacet.sh.nU;
										nV_rotated = collidedFacet.sh.nV;
										N_rotated = collidedFacet.sh.N;
										std::tie(inTheta, inPhi) = CartesianToPolar(collidedFacet.sh.nU, collidedFacet.sh.nV, collidedFacet.sh.N); //Incident angles with original facet surface
									}
										
									//Second step: determine sticking/scattering probabilities
									double stickingProbability; 
									std::vector<double> materialReflProbabilities; //vector of fwd/diffuse/back/transparent scattering probs
									bool complexScattering; //does the material support multiple (diffuse, back, transparent) reflection modes?

									if (collidedFacet.sh.reflectType < 2) {
										//diffuse or mirror with fixed probability
										stickingProbability = collidedFacet.sh.sticking;
										complexScattering = false;
									}
									else {
										//material reflection, depends on incident angle and energy
										std::tie(stickingProbability,materialReflProbabilities, complexScattering) = GetStickingProbability(collidedFacet, inTheta);
									}
									if (!sHandle->ontheflyParams.lowFluxMode) {
										//Regular Monte-Carlo, stick or reflect fwd/diff/back/through
										int reflType = GetHardHitType(stickingProbability, materialReflProbabilities, complexScattering);
										reflected = DoOldRegularReflection(collidedFacet, reflType, inTheta, inPhi, N_rotated, nU_rotated, nV_rotated);
									}
									else {
										//Low flux mode
										reflected = DoLowFluxReflection(collidedFacet, stickingProbability, complexScattering, materialReflProbabilities,
											inTheta, inPhi, 
											N_rotated, nU_rotated, nV_rotated);
									}
								} while (!reflected); //do it again if reflection wasn't successful (reflected against the surface due to roughness)
							}
							else {
								std::string err = "Facet " + (collidedFacet.globalId + 1);
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

std::tuple<double,std::vector<double>,bool> GetStickingProbability(const SubprocessFacet& collidedFacet, const double& theta) {
	//Returns sticking probability, materialReflProbabilities and complexScattering
	Material *mat = &(sHandle->materials[collidedFacet.sh.reflectType - 10]);
	bool complexScattering = mat->hasBackscattering;
	std::vector<double> materialReflProbabilities = mat->BilinearInterpolate(sHandle->energy, abs(theta - PI / 2));
	double stickingProbability = complexScattering
		? 1.0 - materialReflProbabilities[0] - materialReflProbabilities[1] - materialReflProbabilities[2] //100% - forward - diffuse - back (transparent already excluded in Intersect() routine)
		: 1.0 - materialReflProbabilities[0]; //100% - forward (transparent already excluded in Intersect() routine)
	return std::make_tuple(stickingProbability,materialReflProbabilities,complexScattering);
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

std::tuple<Vector3d,Vector3d,Vector3d> PerturbateSurface(const SubprocessFacet& collidedFacet, const double& sigmaRatio) {

	double rnd1 = rnd();
	double rnd2 = rnd(); //for debug
	//Saturate(rnd1, 0.01, 0.99); //Otherwise thetaOffset would go to +/- infinity
	//Saturate(rnd2, 0.01, 0.99); //Otherwise phiOffset would go to +/- infinity
	double thetaOffset = atan(sigmaRatio*tan(PI*(rnd1 - 0.5)));
	double phiOffset = atan(sigmaRatio*tan(PI*(rnd2 - 0.5)));

	//Make a local copy
	Vector3d nU_facet = collidedFacet.sh.nU;
	Vector3d nV_facet = collidedFacet.sh.nV;
	Vector3d N_facet = collidedFacet.sh.N;

	//Random rotation around N (to discard U orientation thus make scattering isotropic)
	double randomAngle = rnd() * 2 * PI;
	nU_facet = Rotate(nU_facet,Vector3d(0,0,0),N_facet, randomAngle);
	nV_facet = Rotate(nV_facet,Vector3d(0,0,0),N_facet, randomAngle);

	//Bending surface base vectors
	Vector3d nU_rotated = Rotate(nU_facet,Vector3d(0,0,0),nV_facet, thetaOffset);
	Vector3d N_rotated = Rotate(N_facet,Vector3d(0,0,0),nV_facet, thetaOffset);
	nU_rotated = Rotate(nU_rotated,Vector3d(0,0,0),nU_facet, phiOffset); //check if correct
	Vector3d nV_rotated = Rotate(nV_facet,Vector3d(0,0,0),nU_facet, phiOffset);
	N_rotated = Rotate(N_rotated,Vector3d(0,0,0),nU_facet, phiOffset);

	return std::make_tuple(nU_rotated, nV_rotated, N_rotated);
}

/*
std::tuple<double,double,double> GetDirComponents(const Vector3d& nU_rotated, const Vector3d& nV_rotated, const Vector3d& N_rotated) {
	double u = Dot(sHandle->pDir, nU_rotated);
	double v = Dot(sHandle->pDir, nV_rotated);
	double n = Dot(sHandle->pDir, N_rotated);
	Saturate(n, -1.0, 1.0); //sometimes rounding errors do occur, 'acos' function would return no value for theta
	return std::make_tuple(u, v, n);
}*/

bool DoLowFluxReflection(SubprocessFacet& collidedFacet, const double& stickingProbability, const bool& complexScattering, const std::vector<double>& materialReflProbabilities,
	const double& inTheta, const double& inPhi,
	const Vector3d& N_rotated, const Vector3d& nU_rotated, const Vector3d& nV_rotated) {

	//First register sticking part:
	collidedFacet.counter.fluxAbs += sHandle->dF*stickingProbability;
	collidedFacet.counter.powerAbs += sHandle->dP*stickingProbability;
	collidedFacet.counter.nbAbsEquiv += stickingProbability;
	if (collidedFacet.texture && collidedFacet.sh.countAbs) RecordHitOnTexture(collidedFacet,
		sHandle->dF*stickingProbability, sHandle->dP*stickingProbability);
	ProfileSlice increment;
	increment.count_absorbed = 0;
	increment.count_incident = 1;
	increment.flux_absorbed = sHandle->dF*stickingProbability;
	increment.flux_incident = sHandle->dF;
	increment.power_absorbed = sHandle->dP*stickingProbability;
	increment.power_incident = sHandle->dP;
	ProfileFacet(collidedFacet, sHandle->energy, increment);
	//Absorbed part recorded, let's see how much is left
	double survivalProbability = 1.0 - stickingProbability;
	sHandle->oriRatio *= survivalProbability;
	if (sHandle->oriRatio < sHandle->ontheflyParams.lowFluxCutoff) {//reflected part not important, throw it away
		RecordHit(HIT_ABS, sHandle->dF, sHandle->dP); //for hits and lines display
		return StartFromSource(); //false if maxdesorption reached
	}
	else { //reflect remainder
		sHandle->dF *= survivalProbability;
		sHandle->dP *= survivalProbability;

		//Decide reflection type (fwd/diff/back, transparent excluded by intersect)
		int reflType;
		if (!complexScattering) reflType = REFL_FORWARD;
		else {
			//Like GetHardHitType() but already excluding absorption
			double anyReflectionProbability = survivalProbability - materialReflProbabilities[3]; //Not absorbed, not transparent
			double random = rnd() * anyReflectionProbability;
			if (random < materialReflProbabilities[0]) reflType = REFL_FORWARD;
			else if (random < (materialReflProbabilities[0] + materialReflProbabilities[1])) reflType = REFL_DIFFUSE;
			else reflType = REFL_BACK;
		}

		if (sHandle->newReflectionModel) {
			PerformBounce_new(collidedFacet, reflType, inTheta, inPhi);
			return true;
		}
		else {
			return PerformBounce_old(collidedFacet, reflType, inTheta, inPhi, N_rotated, nU_rotated, nV_rotated);
		}
	}
}

bool DoOldRegularReflection(SubprocessFacet& collidedFacet, const int& reflType, const double& theta, const double& phi,
	const Vector3d& N_rotated, const Vector3d& nU_rotated, const Vector3d& nV_rotated) {
	
	sHandle->lastHitFacet = &collidedFacet;
	ProfileSlice increment;
	increment.count_absorbed = 0;
	increment.count_incident = 1;
	increment.flux_absorbed = 0.0;
	increment.flux_incident = sHandle->dF;
	increment.power_absorbed = 0.0;
	increment.power_incident = sHandle->dP;
	ProfileFacet(collidedFacet, sHandle->energy, increment);
	if (collidedFacet.texture && collidedFacet.sh.countRefl) RecordHitOnTexture(collidedFacet, sHandle->dF, sHandle->dP);
	if (reflType == REFL_ABSORB) {
				Stick(collidedFacet);
				return StartFromSource(); //false if maxdesorption reached
		} else return PerformBounce_old(collidedFacet, reflType, theta, phi, N_rotated, nU_rotated, nV_rotated);		
}

// Launch photon from a trajectory point

bool StartFromSource() {

	// Check end of simulation
	if (sHandle->ontheflyParams.desorptionLimit > 0) {
		if (sHandle->totalDesorbed >= sHandle->ontheflyParams.desorptionLimit / sHandle->ontheflyParams.nbProcess) {
			sHandle->lastHitFacet = NULL;
			return false;
		}
	}

	//find source point
	size_t pointIdGlobal = (size_t)(rnd()*(double)sHandle->sourceArea);
	bool found = false;
	size_t regionId;
	size_t pointIdLocal;
	size_t sum = 0;
	for (regionId = 0; regionId<sHandle->nbRegion&&!found; regionId++) {
		if ((pointIdGlobal >= sum) && (pointIdGlobal < (sum + sHandle->regions[regionId].Points.size()))) {
			pointIdLocal = pointIdGlobal - sum;
			found = true;
		}
		else sum += sHandle->regions[regionId].Points.size();
	}
	if (!found) SetErrorSub("No start point found");
	regionId--;
	//Trajectory_Point *source=&(sHandle->regions[regionId].Points[pointIdLocal]);
	Region_mathonly *sourceRegion = &(sHandle->regions[regionId]);
	if (!(sourceRegion->params.psimaxX_rad > 0.0 && sourceRegion->params.psimaxY_rad>0.0)) SetErrorSub("psiMaxX or psiMaxY not positive. No photon can be generated");
	
	size_t retries = 0;bool validEnergy;GenPhoton photon;
	do {
		photon = GeneratePhoton(pointIdLocal, sourceRegion, sHandle->ontheflyParams.generation_mode,
			sHandle->psi_distro, sHandle->chi_distros[sourceRegion->params.polarizationCompIndex],
			sHandle->parallel_polarization, sHandle->tmpGlobalCount.nbDesorbed == 0);
		validEnergy = (photon.energy >= sourceRegion->params.energy_low_eV && photon.energy <= sourceRegion->params.energy_hi_eV);
		if (!validEnergy && photon.energy>0.0) {
			retries++;
			pointIdLocal = (size_t)(rnd()*(double)sourceRegion->Points.size());
		}
	} while (!validEnergy && photon.energy>0.0 && retries < 5);

	if (!validEnergy && photon.energy>0.0) {
		char tmp[1024];
		sprintf(tmp, "Region %zd point %zd: can't generate within energy limits (%geV .. %geV)", regionId + 1, pointIdLocal + 1,
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
	sHandle->pPos = photon.start_pos;

	sHandle->sourceRegionId = regionId;

	RecordHit(HIT_DES, sHandle->dF, sHandle->dP);

	//angle
	sHandle->pDir = photon.start_dir;

	if (sourceRegion->params.structureId > sHandle->nbSuper) {
		char tmp[1024];
		sprintf(tmp, "Region %zd is in structure %zd which doesn't exist\n(This geometry has only %zd structures)",
			regionId + 1, sourceRegion->params.structureId + 1, sHandle->nbSuper);
		SetErrorSub(tmp);
		return false;
	}
	sHandle->curStruct = sourceRegion->params.structureId;
	sHandle->teleportedFrom = -1;

	// Count
	sHandle->totalDesorbed++;
	sHandle->tmpGlobalCount.fluxAbs += sHandle->dF;
	sHandle->tmpGlobalCount.powerAbs += sHandle->dP;
	sHandle->tmpGlobalCount.nbDesorbed++;

	sHandle->lastHitFacet = NULL; //Photon originates from the volume, not from a facet

	return true;

}

// Compute bounce against a facet

double TruncatedGaussian(gsl_rng *gen, const double &mean, const double &sigma, const double &lowerBound, const double &upperBound);

void PerformBounce_new(SubprocessFacet& collidedFacet,  const int &reflType, const double &inTheta, const double &inPhi) {

	double outTheta, outPhi; //perform bounce without scattering, will perturbate these angles later if it's a rough surface
	if (collidedFacet.sh.reflectType == REFLECTION_DIFFUSE) {
		outTheta = acos(sqrt(rnd()));
		outPhi = rnd()*2.0*PI;
	}
	else if (collidedFacet.sh.reflectType == REFLECTION_SPECULAR) {
		outTheta = PI - inTheta;
		outPhi = inPhi;
	}
	else { //material reflection, might have backscattering
		switch (reflType) {
		case REFL_FORWARD: //forward scattering
			outTheta = PI - inTheta;
			outPhi = inPhi;
			break;
		case REFL_DIFFUSE: //diffuse scattering
			outTheta = acos(sqrt(rnd()));
			outPhi = rnd()*2.0*PI;
			break;
		case REFL_BACK: //back scattering
			outTheta = PI - inTheta;
			outPhi = PI + inPhi;
			break;
		} //end switch (transparent pass treated at the Intersect() routine
	} //end material reflection

	if (collidedFacet.sh.doScattering) {
		double incidentAngle = abs(inTheta);
		if (incidentAngle > PI / 2) incidentAngle = PI - incidentAngle; //coming from the normal side
		double y = cos(incidentAngle);
		double wavelength = 3E8*6.626E-34 / (sHandle->energy*1.6E-19); //energy[eV] to wavelength[m]
		double specularReflProbability = exp(-Sqr(4 * PI*collidedFacet.sh.rmsRoughness*y / wavelength)); //Debye-Wallers factor, See "Measurements of x-ray scattering..." by Dugan, Sonnad, Cimino, Ishibashi, Scafers, eq.2
		bool specularReflection = rnd() < specularReflProbability;
		if (!specularReflection) {
			//Smooth surface reflection performed, now let's perturbate the angles
			//Using Gaussian approximated distributions of eq.14. of the above article
			double onePerTau = collidedFacet.sh.rmsRoughness / collidedFacet.sh.autoCorrLength;
			
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

	PolarToCartesian(&collidedFacet, outTheta, outPhi, false);

	RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
	sHandle->lastHitFacet = &collidedFacet;
	if (collidedFacet.texture && collidedFacet.sh.countRefl) RecordHitOnTexture(collidedFacet, sHandle->dF, sHandle->dP);
}

bool PerformBounce_old(SubprocessFacet& collidedFacet, const int& reflType, const double& inTheta, const double& inPhi,
	const Vector3d& N_rotated, const Vector3d& nU_rotated, const Vector3d& nV_rotated) {
	RecordHit(HIT_REF, sHandle->dF, sHandle->dP);
	// Relaunch particle, regular monte-carlo
	if (collidedFacet.sh.reflectType == REFLECTION_DIFFUSE) {
		//See docs/theta_gen.png for further details on angular distribution generation
		PolarToCartesian(&collidedFacet, acos(sqrt(rnd())), rnd()*2.0*PI, false);
	} else { //Fwd/diff/back reflection, optionally with surface perturbation
		if (!VerifiedSpecularReflection(collidedFacet, (collidedFacet.sh.reflectType == REFLECTION_SPECULAR)?REFL_FORWARD:reflType, inTheta, inPhi,
			nU_rotated, nV_rotated, N_rotated)) {
			return false;
		}
	}
	return true;
}

bool VerifiedSpecularReflection(const SubprocessFacet& collidedFacet, const int& reflType, const double& inTheta, const double& inPhi,
	const Vector3d& nU_rotated, const Vector3d& nV_rotated,const Vector3d& N_rotated) {
	
	//Specular reflection that returns false if going against surface
	//Changes sHandle->pDir

	double outTheta, outPhi;

	//Vector3d N_facet = Vector3d(collidedFacet->sh.N.x, collidedFacet->sh.N.y, collidedFacet->sh.N.z);

	Vector3d newDir;
	double u, v, n;
	bool calcNewDir = true;

	switch (reflType) {
	case REFL_FORWARD: //forward scattering
		outTheta = PI - inTheta;
		outPhi = inPhi;
		break;
	case REFL_DIFFUSE: //diffuse scattering
		outTheta = acos(sqrt(rnd()));
		outPhi = rnd()*2.0*PI;
		break;
	case REFL_BACK: //back scattering

		//we need to perturbate the backscattered ray with the angle difference between the original and the rotated surface
		
		//Get rotation matrix that transforms collidedFacet->sh.N into N_rotated
		//See https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
		
		Vector3d v = CrossProduct(collidedFacet.sh.N, N_rotated);
		double c = Dot(collidedFacet.sh.N, N_rotated);
		
		double factor = 1.0 / (1.0 + c);

		std::vector<std::vector<double>> nullMatrix = { {0.0,0.0,0.0} , {0.0,0.0,0.0} , {0.0,0.0,0.0} };
		std::vector<std::vector<double>> identityMatrix = { {1.0,0.0,0.0} , {0.0,1.0,0.0} , {0.0,0.0,1.0} }; //rows
		std::vector<std::vector<double>> v_skew = { {0.0,-v.z,v.y} , {v.z,0.0,-v.x} , {-v.y,v.x,0.0} };
		
		std::vector<std::vector<double>> v_skew_square = nullMatrix;
		//3x3 matrix multiplication
		for (size_t row = 0; row < 3; row++) {
			for (size_t col = 0; col < 3; col++) {
				for (size_t comp = 0; comp < 3; comp++) {
					v_skew_square[row][col] += v_skew[row][comp] * v_skew[comp][col];
				}
			}
		}

		std::vector<std::vector<double>> rotationMatrix = nullMatrix;
		for (size_t row = 0; row < 3; row++) {
			for (size_t col = 0; col < 3; col++) {
				rotationMatrix[row][col] = identityMatrix[row][col] + v_skew[row][col] + factor*v_skew_square[row][col];
			}
		}

		newDir = Vector3d(-sHandle->pDir.x * rotationMatrix[0][0] + -sHandle->pDir.y * rotationMatrix[0][1] + -sHandle->pDir.z * rotationMatrix[0][2],
			-sHandle->pDir.x * rotationMatrix[1][0] + -sHandle->pDir.y * rotationMatrix[1][1] + -sHandle->pDir.z * rotationMatrix[1][2],
			-sHandle->pDir.x * rotationMatrix[2][0] + -sHandle->pDir.y * rotationMatrix[2][1] + -sHandle->pDir.z * rotationMatrix[2][2]);

		calcNewDir = false;
		break;
	} //end switch (transparent pass treated at the Intersect() routine

	if (calcNewDir) {

		u = sin(outTheta)*cos(outPhi);
		v = sin(outTheta)*sin(outPhi);
		n = cos(outTheta);

		newDir = u*nU_rotated + v*nV_rotated + n*N_rotated;

	}

	if ((Dot(newDir,collidedFacet.sh.N) > 0) != (inTheta > PI / 2)) {
		//inTheta > PI/2: ray coming from normal side
		//Dot(newDir,N_facet)>0: ray leaving towards normal side
		return false; //if reflection would go against the surface, generate new angles
	}

	sHandle->pDir = newDir;
	return true;
}

void RecordHitOnTexture(SubprocessFacet& f, double dF, double dP) {
	size_t tu = (size_t)(f.colU * f.sh.texWidthD);
	size_t tv = (size_t)(f.colV * f.sh.texHeightD);
	size_t index = tu + tv*f.sh.texWidth;
	f.texture[index].count++;
	f.texture[index].flux += dF*f.inc[index]; //normalized by area
	f.texture[index].power += dP*f.inc[index]; //normalized by area
}

void RecordDirectionVector(SubprocessFacet& f) {
	size_t tu = (size_t)(f.colU * f.sh.texWidthD);
	size_t tv = (size_t)(f.colV * f.sh.texHeightD);
	size_t add = tu + tv*(f.sh.texWidth);

	f.direction[add].dir.x += sHandle->pDir.x;
	f.direction[add].dir.y += sHandle->pDir.y;
	f.direction[add].dir.z += sHandle->pDir.z;
	f.direction[add].count++;
}

void Stick(SubprocessFacet& collidedFacet) {
	collidedFacet.counter.nbAbsEquiv+=1.0;
	collidedFacet.counter.fluxAbs += sHandle->dF;
	collidedFacet.counter.powerAbs += sHandle->dP;
	sHandle->tmpGlobalCount.nbAbsEquiv+=1.0;
	//sHandle->distTraveledSinceUpdate+=sHandle->distTraveledCurrentParticle;
	//sHandle->counter.nbAbsorbed++;
	//sHandle->counter.fluxAbs+=sHandle->dF;
	//sHandle->counter.powerAbs+=sHandle->dP;
	RecordHit(HIT_ABS, sHandle->dF, sHandle->dP); //for hits and lines display
	ProfileSlice increment;
	increment.count_absorbed = 1;
	increment.count_incident = 1;
	increment.flux_absorbed = sHandle->dF;
	increment.flux_incident = sHandle->dF;
	increment.power_absorbed = sHandle->dP;
	increment.power_incident = sHandle->dP;
	ProfileFacet(collidedFacet, sHandle->energy, increment);
	if (collidedFacet.texture && collidedFacet.sh.countAbs) RecordHitOnTexture(collidedFacet, sHandle->dF, sHandle->dP);
}

void SubprocessFacet::RegisterTransparentPass()
{
	//Low flux mode not supported (ray properties can't change on transparent pass since it's inside the Intersect() routine)
	this->hitted = true;
	this->counter.nbMCHit++; //count MC hits
	this->counter.nbHitEquiv += sHandle->oriRatio;
	this->counter.fluxAbs += sHandle->dF;
	this->counter.powerAbs += sHandle->dP;
	ProfileSlice increment;
	increment.count_absorbed = 0;
	increment.count_incident = 1;
	increment.flux_absorbed = 0.0;
	increment.flux_incident = sHandle->dF;
	increment.power_absorbed = 0.0;
	increment.power_incident = sHandle->dP;
	ProfileFacet(*this, sHandle->energy, increment); //count profile
	LogHit(this);
	if (this->texture && this->sh.countTrans) RecordHitOnTexture(*this, sHandle->dF, sHandle->dP); //count texture
	if (this->direction && this->sh.countDirection) RecordDirectionVector(*this);
}

void ProfileFacet(SubprocessFacet &f, const double &energy, const ProfileSlice& increment) {

	size_t pos;

	switch (f.sh.profileType) {

	case PROFILE_ANGULAR: {
		double dot = abs(Dot(f.sh.N, sHandle->pDir));
		double theta = acos(dot);              // Angle to normal (0 to PI/2)
		size_t grad = (size_t)(((double)PROFILE_SIZE)*(theta) / (PI / 2)); // To Grad
		Saturate(grad, 0, PROFILE_SIZE - 1);
		f.profile[grad] += increment;

	} break;

	case PROFILE_U:
		pos = (size_t)((f.colU)*(double)PROFILE_SIZE);
		Saturate(pos, 0, PROFILE_SIZE - 1);
		f.profile[pos]+= increment;
		break;

	case PROFILE_V:
		pos = (size_t)((f.colV)*(double)PROFILE_SIZE);
		Saturate(pos, 0, PROFILE_SIZE - 1);
		f.profile[pos] += increment;
		break;

	}

	if (f.sh.recordSpectrum) {
		f.spectrum->Add(energy, increment);
	}
}

void SubprocessFacet::ResetCounter() {
	memset(&counter, 0, sizeof(counter));
}

void UpdateLog(Dataport * dpLog, DWORD timeout)
{
	if (sHandle->tmpParticleLog.size()) {
		double t0, t1;
		t0 = GetTick();
		SetState(NULL, "Waiting for 'dpLog' dataport access...", false, true);
		sHandle->lastLogUpdateOK = AccessDataportTimed(dpLog, timeout);
		SetState(NULL, "Updating Log...", false, true);
		if (!sHandle->lastLogUpdateOK) return;

		size_t* logBuff = (size_t*)dpLog->buff;
		size_t recordedLogSize = *logBuff;
		ParticleLoggerItem* logBuff2 = (ParticleLoggerItem*)(logBuff + 1);

		size_t writeNb;
		if (recordedLogSize > sHandle->ontheflyParams.logLimit) writeNb = 0;
		else writeNb = Min(sHandle->tmpParticleLog.size(), sHandle->ontheflyParams.logLimit - recordedLogSize);
		memcpy(&logBuff2[recordedLogSize], &sHandle->tmpParticleLog[0], writeNb * sizeof(ParticleLoggerItem)); //Knowing that vector memories are contigious
		(*logBuff) += writeNb;
		ReleaseDataport(dpLog);
		sHandle->tmpParticleLog.clear();
		extern char* GetSimuStatus();
		SetState(NULL, GetSimuStatus(), false, true);

#ifdef _DEBUG
		t1 = GetTick();
		printf("Update log: %f us\n", (t1 - t0)*1000000.0);
#endif
	}
}

void LogHit(SubprocessFacet * f)
{
	if (sHandle->ontheflyParams.enableLogging &&
		sHandle->ontheflyParams.logFacetId == f->globalId &&
		sHandle->tmpParticleLog.size() < (sHandle->ontheflyParams.logLimit / sHandle->ontheflyParams.nbProcess)) {
		ParticleLoggerItem log;
		log.facetHitPosition = Vector2d(f->colU, f->colV);
		std::tie(log.hitTheta, log.hitPhi) = CartesianToPolar(f->sh.nU, f->sh.nV, f->sh.N);
		log.oriRatio = sHandle->oriRatio;
		log.energy = sHandle->energy;
		log.dF = sHandle->dF;
		log.dP = sHandle->dP;
		sHandle->tmpParticleLog.push_back(log);
	}
}
