/*
File:        SimulationMC.c
Description: Monte-Carlo Simulation for UHV (Physics related routines) 
Program:     SynRad
Author:      R. KERSEVAN / M SZAKACS
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

extern SIMULATION *sHandle;
extern void SetErrorSub(char *message);

extern Distribution2D K_1_3_distribution;
extern Distribution2D K_2_3_distribution;
extern Distribution2D integral_N_photons;
extern Distribution2D integral_SR_power;
extern Distribution2D polarization_distribution;
extern Distribution2D g1h2_distribution;


void ComputeSourceArea() {
	sHandle->sourceArea=sHandle->nbTrajPoints;
}




// -------------------------------------------------------

void PolarToCartesian(FACET *iFacet,double theta,double phi,BOOL reverse) {

	VERTEX3D U,V,N;
	double u,v,n;

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
		N.x=N.x*(-1.0);
		N.y=N.y*(-1.0);
		N.z=N.z*(-1.0);
	}

	// Basis change (nU,nV,N) -> (x,y,z)
	sHandle->pDir.x = u*U.x + v*V.x + n*N.x;
	sHandle->pDir.y = u*U.y + v*V.y + n*N.y;
	sHandle->pDir.z = u*U.z + v*V.z + n*N.z;
	//_ASSERTE(Norme(&sHandle->pDir)<=1.0);
}

// -------------------------------------------------------

void CartesianToPolar(FACET *iFacet,double *theta,double *phi) {

	// Get polar coordinates of the incoming particule direction in the (U,V,N) facet space.
	// Note: The facet is parallel to (U,V), we use its (nU,nV,N) orthonormal basis here.
	// (nU,nV,N) and (x,y,z) are both left handed

	// Cartesian(x,y,z) to polar in (nU,nV,N) transformation

	// Basis change (x,y,z) -> (nU,nV,N)
	// We use the fact that (nU,nV,N) belongs to SO(3)
	double u = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
		iFacet->sh.nU.x,iFacet->sh.nU.y,iFacet->sh.nU.z);
	double v = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
		iFacet->sh.nV.x,iFacet->sh.nV.y,iFacet->sh.nV.z);
	double n = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
		iFacet->sh.N.x,iFacet->sh.N.y,iFacet->sh.N.z);
	SATURATE(n,-1.0,1.0); //sometimes rounding errors do occur, 'acos' function would return no value for theta

	// (u,v,n) -> (theta,phi)
	double rho = sqrt( v*v + u*u );
	*theta = acos(n);              // Angle to normal (PI/2 => PI)
	*phi = asin(v/rho);
	if( u<0.0 ) *phi = PI - *phi;  // Angle to U

}

int RoughReflection(FACET *iFacet,double sigmaRatio,double theta,double phi,
	Vector N_rotated,Vector nU_rotated,Vector nV_rotated) {
		//_ASSERTE(theta>(PI/2.0));

	Vector nU_facet=Vector(iFacet->sh.nU.x,iFacet->sh.nU.y,iFacet->sh.nU.z);
	Vector nV_facet=Vector(iFacet->sh.nV.x,iFacet->sh.nV.y,iFacet->sh.nV.z);
	Vector N_facet=Vector(iFacet->sh.N.x,iFacet->sh.N.y,iFacet->sh.N.z);

	Vector newDir;
	double u,v,n;

	//Polar to cartesian routineF-_
	theta=PI-theta; //perform reflection

	u = sin(theta)*cos(phi);
	v = sin(theta)*sin(phi);
	n = cos(theta);

	newDir=Vector(u*nU_rotated.x + v*nV_rotated.x + n*N_rotated.x,
		u*nU_rotated.y + v*nV_rotated.y + n*N_rotated.y,
		u*nU_rotated.z + v*nV_rotated.z + n*N_rotated.z);

	if ((DOT3(newDir.x,newDir.y,newDir.z,
		N_facet.x,N_facet.y,N_facet.z)>0) != (theta<PI/2)) {
			return 0; //if reflection would go against the surface, generate new angles
	}

	sHandle->pDir.x = newDir.x;
	sHandle->pDir.y = newDir.y;
	sHandle->pDir.z = newDir.z;
	//_ASSERTE(Norme(&sHandle->pDir)<=1.0);
	return 1;
}

// -------------------------------------------------------

void UpdateMCHits(Dataport *dpHit,int prIdx,DWORD timeout) {

	BYTE *buffer;
	SHGHITS *gHits;
	llong minHitsOld_MC;
	llong maxHitsOld_MC;
	double minHitsOld_flux,minHitsOld_power;
	double maxHitsOld_flux,maxHitsOld_power;
	int i,j,s,x,y,nb;
#ifdef _DEBUG
	double t0,t1;
	t0 = GetTick();
#endif
	sHandle->lastUpdateOK = AccessDataportTimed(dpHit,timeout);
	if( !sHandle->lastUpdateOK ) return;

	buffer = (BYTE*)dpHit->buff;
	gHits = (SHGHITS *)buffer;

	// Global hits and leaks
	gHits->total.nbHit      += sHandle->tmpCount.nbHit;
	gHits->total.nbAbsorbed += sHandle->tmpCount.nbAbsorbed;
	gHits->total.nbDesorbed += sHandle->tmpCount.nbDesorbed;
	gHits->distTraveledTotal    += sHandle->distTraveledSinceUpdate;
	gHits->total.fluxAbs += sHandle->tmpCount.fluxAbs;
	gHits->total.powerAbs += sHandle->tmpCount.powerAbs;

	minHitsOld_MC=gHits->minHit_MC;
	maxHitsOld_MC=gHits->maxHit_MC;
	minHitsOld_flux=gHits->minHit_flux;
	maxHitsOld_flux=gHits->maxHit_flux;
	minHitsOld_power=gHits->minHit_power;
	maxHitsOld_power=gHits->maxHit_power;

	gHits->minHit_MC=HITMAX_INT64;
	gHits->maxHit_MC=0;
	gHits->minHit_flux=HITMAX_DOUBLE;
	gHits->maxHit_flux=0;
	gHits->minHit_power=HITMAX_DOUBLE;
	gHits->maxHit_power=0;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak
	nb = gHits->nbLastLeaks;
	for(i=0;i<sHandle->nbLastLeak && i<NBHLEAK;i++)
		gHits->pLeak[(i+nb) % NBHLEAK] = sHandle->pLeak[i];
	gHits->nbLeakTotal += sHandle->nbLeakTotal;
	gHits->nbLastLeaks += sHandle->nbLastLeak;

	// HHit (Only prIdx 0)
	if( prIdx==0 ) {
		gHits->nbHHit = sHandle->nbHHit;
		memcpy(gHits->pHits,sHandle->pHits,NBHHIT*sizeof(HIT));
	}

	// Facets
	for(s=0;s<sHandle->nbSuper;s++) {
		for(i=0;i<sHandle->str[s].nbFacet;i++) {

			FACET *f = sHandle->str[s].facets[i];
			if( f->hitted ) {

				SHHITS *fFit = (SHHITS *)(buffer + f->sh.hitOffset);
				fFit->nbAbsorbed += f->sh.counter.nbAbsorbed;
				fFit->nbDesorbed += f->sh.counter.nbDesorbed;
				fFit->fluxAbs += f->sh.counter.fluxAbs;
				fFit->powerAbs += f->sh.counter.powerAbs;
				fFit->nbHit += f->sh.counter.nbHit;

				if( f->sh.isProfile ) {
					llong *shProfile = (llong *)(buffer + (f->sh.hitOffset + sizeof(SHHITS)));
					for(j=0;j<PROFILE_SIZE;j++)
						shProfile[j] += f->profile_hits[j];

					double *shProfile2 = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS))+PROFILE_SIZE*sizeof(llong));
					for(j=0;j<PROFILE_SIZE;j++)
						shProfile2[j] += f->profile_flux[j];

					double *shProfile3 = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS))+PROFILE_SIZE*(sizeof(double)+sizeof(llong)));
					for(j=0;j<PROFILE_SIZE;j++)
						shProfile3[j] += f->profile_power[j];
				}


				int profileSize=(f->sh.isProfile)?PROFILE_SIZE*(2*sizeof(double)+sizeof(llong)):0;
				if( f->sh.isTextured ) {
					int textureElements=f->sh.texHeight*f->sh.texWidth;
					llong *shTexture = (llong *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize));
					for(y=0;y<f->sh.texHeight;y++) {
						for(x=0;x<f->sh.texWidth;x++) {
							int add = x + y*f->sh.texWidth;
							llong val = shTexture[add] + f->hits_MC[add];
							if(val>gHits->maxHit_MC)	//no cell size check for MC						
								gHits->maxHit_MC=val;
							if (val>0 && val<gHits->minHit_MC)
								gHits->minHit_MC=val;
							shTexture[add] = val;
						}
					}

					double *shTexture2 = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize+textureElements*sizeof(llong)));
					for(y=0;y<f->sh.texHeight;y++) {
						for(x=0;x<f->sh.texWidth;x++) {
							int add = x + y*f->sh.texWidth;
							double val2 = shTexture2[add] + f->hits_flux[add];
							if(val2>gHits->maxHit_flux&& f->largeEnough[add]) {
								gHits->maxHit_flux=val2;
							}
							if (val2>0.0 && val2<gHits->minHit_flux&& f->largeEnough[add]) gHits->minHit_flux=val2;
							shTexture2[add] = val2;
						}
					}

					double *shTexture3 = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize+textureElements*(sizeof(llong)+sizeof(double))));
					for(y=0;y<f->sh.texHeight;y++) {
						for(x=0;x<f->sh.texWidth;x++) {
							int add = x + y*f->sh.texWidth;
							double val3 = shTexture3[add] + f->hits_power[add];
							if(val3>gHits->maxHit_power&& f->largeEnough[add]) 
								gHits->maxHit_power=val3;
							if (val3>0.0 && val3<gHits->minHit_power&& f->largeEnough[add]) gHits->minHit_power=val3; //disregard very small elements
							shTexture3[add] = val3;
						}
					}
				}

				if( f->sh.countDirection ) {
					VHIT *shDir = (VHIT *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + f->profileSize + f->textureSize));
					for(y=0;y<f->sh.texHeight;y++) {
						for(x=0;x<f->sh.texWidth;x++) {
							int add = x + y*f->sh.texWidth;
							shDir[add].dir.x += f->direction[add].dir.x;
							shDir[add].dir.y += f->direction[add].dir.y;
							shDir[add].dir.z += f->direction[add].dir.z;
							shDir[add].count += f->direction[add].count;
						}
					}
				}

				if ( f->sh.hasSpectrum ) {
					double *shSpectrum_fluxwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + f->profileSize 
						+ f->textureSize + f->directionSize));
					for(j=0;j<SPECTRUM_SIZE;j++)
						shSpectrum_fluxwise[j] += f->spectrum_fluxwise->GetCount(j);

					double *shSpectrum_powerwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + f->profileSize 
						+ f->textureSize + f->directionSize + SPECTRUM_SIZE*sizeof(double)));
					for(j=0;j<SPECTRUM_SIZE;j++)
						shSpectrum_powerwise[j] += f->spectrum_powerwise->GetCount(j);
				}

			} // End if(hitted)
		} // End nbFacet
	} // End nbSuper

	//if there were no textures:
	if (gHits->minHit_MC==HITMAX_INT64) gHits->minHit_MC=minHitsOld_MC;
	if (gHits->maxHit_MC==0) gHits->maxHit_MC=maxHitsOld_MC;
	if (gHits->minHit_flux==HITMAX_DOUBLE) gHits->minHit_flux=minHitsOld_flux;
	if (gHits->maxHit_flux==0.0) gHits->maxHit_flux=maxHitsOld_flux;
	if (gHits->minHit_power==HITMAX_DOUBLE) gHits->minHit_power=minHitsOld_power;
	if (gHits->maxHit_power==0.0) gHits->maxHit_power=maxHitsOld_power;

	ReleaseDataport(dpHit);

	//printf("\nResetCounter called from UpdateMCHits");
	ResetCounter();

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n",(t1-t0)*1000000.0);
#endif

}






// -------------------------------------------------------------
// Compute particle teleport
// -------------------------------------------------------------

void PerformTeleport(FACET *iFacet) {

	double inPhi,inTheta;
	//Search destination
	FACET *destination;
	BOOL found=FALSE;
	BOOL revert=FALSE;
	int i,j;
	//Look in which superstructure is the destination facet:
	if (sHandle->nbSuper==1) { //speedup for mono-structure systems
		found=iFacet->sh.teleportDest<=sHandle->str[0].nbFacet;
		if (found) destination=sHandle->str[0].facets[iFacet->sh.teleportDest-1];
	} else {
		for (i=0;i<sHandle->nbSuper&&(!found);i++) {
			for (j=0;j<sHandle->str[i].nbFacet&&(!found);j++) {
				if ((iFacet->sh.teleportDest-1)==sHandle->str[i].facets[j]->globalId) {
					destination=sHandle->str[i].facets[j];
					sHandle->curStruct = destination->sh.superIdx; //change current superstructure
					found=TRUE;
				}
			}
		}
	}
	if (!found) {
		char err[128];
		sprintf(err,"Teleport destination of facet %d not found (facet %d does not exist)",iFacet->globalId+1,iFacet->sh.teleportDest);
		SetErrorSub(err);
		return;
	}

	// Count this hit as a transparent pass
	RecordHit(HIT_TELEPORT,sHandle->dF,sHandle->dP);
	if( iFacet->hits_MC && iFacet->sh.countTrans ) RecordHitOnTexture(iFacet,sHandle->dF,sHandle->dP);


	// Relaunch particle from new facet
	CartesianToPolar(iFacet,&inTheta,&inPhi);
	PolarToCartesian(destination,inTheta,inPhi,FALSE);
	// Move particle to teleport destination point
	sHandle->pPos.x = destination->sh.O.x+iFacet->colU*destination->sh.U.x+iFacet->colV*destination->sh.V.x;
	sHandle->pPos.y = destination->sh.O.y+iFacet->colU*destination->sh.U.y+iFacet->colV*destination->sh.V.y;
	sHandle->pPos.z = destination->sh.O.z+iFacet->colU*destination->sh.U.z+iFacet->colV*destination->sh.V.z;
	RecordHit(HIT_TELEPORT,sHandle->dF,sHandle->dP);
	sHandle->lastHit = destination;

	//Count hits on teleport facets (only TP source)
	iFacet->sh.counter.nbAbsorbed++;
	//destination->sh.counter.nbDesorbed++;

	iFacet->sh.counter.nbHit++;//destination->sh.counter.nbHit++;
	iFacet->sh.counter.fluxAbs+=sHandle->dF;//destination->sh.counter.fluxAbs+=sHandle->dF;
	iFacet->sh.counter.powerAbs+=sHandle->dP;//destination->sh.counter.powerAbs+=sHandle->dP;
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
	for(i=0;i<nbStep;i++) {

		found = Intersect(&(sHandle->pPos),&(sHandle->pDir),&d,&collidedFacet,sHandle->lastHit);

		if( found ) {

			// Move particule to intersection point
			sHandle->pPos.x += d*sHandle->pDir.x;
			sHandle->pPos.y += d*sHandle->pDir.y;
			sHandle->pPos.z += d*sHandle->pDir.z;
			sHandle->distTraveledCurrentParticle += d;

			if (collidedFacet->sh.teleportDest) {
				PerformTeleport(collidedFacet); 
			} else if(collidedFacet->sh.reflectType>=2) { //rough surface reflection
				if (sHandle->dF==0.0 || sHandle->dP==0.0 || sHandle->energy<1E-3) { //stick non real photons (from beam beginning)
					if (!Stick(collidedFacet)) return FALSE;
				} else {
				double sigmaRatio=collidedFacet->sh.roughness; 
				for (int m=0;m<(int)sHandle->materials.size();m++) {
					if (collidedFacet->sh.reflectType==m+2) { //found material type
						//Generate incident angle
							Vector nU_facet=Vector(collidedFacet->sh.nU.x,collidedFacet->sh.nU.y,collidedFacet->sh.nU.z);
							Vector nV_facet=Vector(collidedFacet->sh.nV.x,collidedFacet->sh.nV.y,collidedFacet->sh.nV.z);
							Vector N_facet=Vector(collidedFacet->sh.N.x,collidedFacet->sh.N.y,collidedFacet->sh.N.z);
						Vector nU_rotated,N_rotated,nV_rotated;
						double n,u,v,thetaOffset,phiOffset;
						BOOL reflected;
						do {
							double n_ori= DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
									N_facet.x,N_facet.y,N_facet.z); //original incident angle (negative if front collision, positive if back collision);
							do {
								//Generating bending angle due to local roughness
								double rnd1=rnd();
								double rnd2=rnd(); //for debug
								thetaOffset=atan(sigmaRatio*tan(PI*(rnd1-0.5)));
								phiOffset=atan(sigmaRatio*tan(PI*(rnd2-0.5)));

								//Random rotation around N (to discard U orientation thus make scattering isotropic)
								double randomAngle=rnd()*2*PI;
								nU_facet=nU_facet.Rotate(N_facet,randomAngle);
								nV_facet=nV_facet.Rotate(N_facet,randomAngle);
								
								//Bending surface base vectors
								nU_rotated=nU_facet.Rotate(nV_facet,thetaOffset);
								N_rotated=N_facet.Rotate(nV_facet,thetaOffset);
								nU_rotated=nU_rotated.Rotate(nU_facet,phiOffset); //check if correct
								nV_rotated=nV_facet.Rotate(nU_facet,phiOffset);
								N_rotated=N_rotated.Rotate(nU_facet,phiOffset);

								//Calculate incident angles to bent surface (Cartesian to Polar routine modified)
								u = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
									nU_rotated.x,nU_rotated.y,nU_rotated.z);
								v = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
									nV_rotated.x,nV_rotated.y,nV_rotated.z);
								n = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
									N_rotated.x,N_rotated.y,N_rotated.z);
								SATURATE(n,-1.0,1.0); //sometimes rounding errors do occur, 'acos' function would return no value for theta

							} while (n*n_ori<=0.0 ); //regenerate random numbers if grazing angle would be over 90deg

							double rho = sqrt( v*v + u*u );
							double theta = acos(n);              // Angle to normal (PI/2 => PI)
							double phi = asin(v/rho);
							if( u<0.0 ) phi = PI - phi;  // Angle to U
							double reflProbability=sHandle->materials[m].Interpolate(sHandle->energy,(theta-PI/2)*1000.0)/100.0;
							if (rnd()<reflProbability) {
								reflected=PerformBounce(collidedFacet,sigmaRatio,theta,phi,N_rotated,nU_rotated,nV_rotated);
							} else { //sticking
								if (!Stick(collidedFacet)) return FALSE;
								reflected=TRUE;
							}
						} while (!reflected); //do it again if reflection wasn't successful (ie. new angle was over 90deg)
						//}while (FALSE);
						break;
					}
				}
				}
			} else if( collidedFacet->sh.sticking>0.0 ) {

				if( collidedFacet->sh.sticking==1.0 || rnd()<collidedFacet->sh.sticking ) { //sticking
					if (!Stick(collidedFacet)) return FALSE;
				} else {
					PerformBounce(collidedFacet);
				}

			} else {
				PerformBounce(collidedFacet); //if sticking==0, bounce without generating random number
			}

			if (!collidedFacet->sh.teleportDest){
				// Hit count
				sHandle->tmpCount.nbHit++;
				//sHandle->counter.nbHit++;
				collidedFacet->sh.counter.nbHit++;
				//dF,dP comes here?
			}

		} else {

			// Leak (simulation error)
			RecordLeakPos();
			sHandle->nbLeakTotal++;
			if( !StartFromSource() )
				// maxDesorption reached
				return FALSE;

		}

	}

	return TRUE;

}

// -------------------------------------------------------------
// Launch photon from a trajectory point
// -------------------------------------------------------------

BOOL StartFromSource() {

	// Check end of simulation
	if( sHandle->maxDesorption>0 ) {
		if( sHandle->totalDesorbed>=sHandle->maxDesorption ) {
			sHandle->lastHit=NULL;
			return FALSE;
		}
	}

	//find source point
	double pointIdGlobal=rnd()*sHandle->sourceArea;
	BOOL found=false;
	int regionId;
	double pointIdLocal;
	double sum=0.0;
	for (regionId=0;regionId<sHandle->nbRegion&&!found;regionId++) {
		if ((pointIdGlobal>=sum) && (pointIdGlobal<(sum+(double)sHandle->regions[regionId].Points.size()))) {
			pointIdLocal=pointIdGlobal-sum;
			found=true;
		} else sum+=(int)sHandle->regions[regionId].Points.size();
	}
	if (!found) SetErrorSub("No start point found");
	regionId--;
	//Trajectory_Point *source=&(sHandle->regions[regionId].Points[pointIdLocal]);
	Region_mathonly *sourceRegion=&(sHandle->regions[regionId]);
	if (!(sourceRegion->psimaxX>0.0 && sourceRegion->psimaxY>0.0)) SetErrorSub("psiMaxX or psiMaxY not positive. No photon can be generated");
	GenPhoton photon=Radiate(pointIdLocal,sourceRegion);

	sHandle->distTraveledCurrentParticle=0.0;
	sHandle->dF=photon.dF;
	sHandle->dP=photon.dP;
	sHandle->energy=photon.energy;

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
	
	RecordHit(HIT_DES,sHandle->dF,sHandle->dP);

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
	sHandle->curStruct=0;

	// Count
	sHandle->totalDesorbed++;
	//sHandle->counter.nbDesorbed++;
	sHandle->tmpCount.nbDesorbed++;
	sHandle->nbPHit = 0;

	sHandle->lastHit=NULL;

	return TRUE;

}

// -------------------------------------------------------------
// Compute bounce against a facet
// -------------------------------------------------------------

int PerformBounce(FACET *iFacet,double sigmaRatio,double theta,double phi,
	Vector N_rotated,Vector nU_rotated,Vector nV_rotated) {

	double inPhi,inTheta;
	BOOL revert=FALSE;

	// Handle super structure link facet
	if( iFacet->sh.superDest ) {

		sHandle->curStruct = iFacet->sh.superDest - 1;
		// Count this hit as a transparent pass
		RecordHit(HIT_TRANS,sHandle->dF,sHandle->dP);
		if( iFacet->hits_MC && iFacet->sh.countTrans ) RecordHitOnTexture(iFacet,sHandle->dF,sHandle->dP);
		return 1;

	}

	// Handle volatile facet
	if( iFacet->sh.isVolatile ) {

		if( iFacet->ready ) {
			iFacet->sh.counter.nbAbsorbed++;
			iFacet->sh.counter.fluxAbs+=sHandle->dF;
			iFacet->sh.counter.powerAbs+=sHandle->dP;
			iFacet->ready = FALSE;
			if( iFacet->hits_MC && iFacet->sh.countAbs ) RecordHitOnTexture(iFacet,sHandle->dF,sHandle->dP);
		}
		return 1;
	}

	if( iFacet->sh.is2sided ) {
		// We may need to revert normal in case of 2 sided hit
		revert = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
			iFacet->sh.N.x,iFacet->sh.N.y,iFacet->sh.N.z)>0.0;
	}

	// Relaunch particle
	switch( iFacet->sh.reflectType ) {
	case REF_MIRROR:
		//here comes a super complicated photon reflection routine?
		CartesianToPolar(iFacet,&inTheta,&inPhi);
		PolarToCartesian(iFacet,PI-inTheta,inPhi,FALSE);
		break;
	case REF_DIFFUSE:
		//See docs/theta_gen.png for further details on angular distribution generation
		PolarToCartesian(iFacet,acos(sqrt(rnd())),rnd()*2.0*PI,FALSE);
		break;
	default: //Material reflection - theta, etc. are already calculated
		for (int i=0;i<(sHandle->nbMaterials);i++) {
			if (iFacet->sh.reflectType==(2+i)) {
				if (!RoughReflection(iFacet,sigmaRatio,theta,phi,
					N_rotated,nU_rotated,nV_rotated)) {
						return 0;
				}
				break;
			}
		}
	}
	if( revert ) {
		sHandle->pDir.x = -sHandle->pDir.x;
		sHandle->pDir.y = -sHandle->pDir.y;
		sHandle->pDir.z = -sHandle->pDir.z;
	}

	RecordHit(HIT_REF,sHandle->dF,sHandle->dP);
	sHandle->lastHit = iFacet;
	sHandle->nbPHit++;


	if( iFacet->hits_MC && iFacet->sh.countRefl ) RecordHitOnTexture(iFacet,sHandle->dF,sHandle->dP);
	//sHandle->temperature = iFacet->sh.temperature; //Thermalize particle
	return 1;
}

void RecordHitOnTexture(FACET *f,double dF,double dP) {                            
	int tu = (int)(f->colU * f->sh.texWidthD);  
	int tv = (int)(f->colV * f->sh.texHeightD); 
	int add = tu+tv*f->sh.texWidth;  
	f->hits_MC[add]++;
	f->hits_flux[add]+=dF*f->inc[add]; //normalized by area
	f->hits_power[add]+=dP*f->inc[add]; //normalized by area
}

GenPhoton Radiate(double sourceId,Region_mathonly *current_region) { //Generates a photon from point number 'sourceId'

	//Interpolate source point
	SATURATE(sourceId,0,(double)current_region->Points.size()-1.0000001); //make sure we stay within limits
	Trajectory_Point previousPoint=current_region->Points[(int)sourceId];
	Trajectory_Point nextPoint=current_region->Points[(int)sourceId+1];
	double overshoot=sourceId-(int)sourceId;
	
	Trajectory_Point source;
	source.position=Weigh(previousPoint.position,nextPoint.position,overshoot);
	source.direction=Weigh(previousPoint.direction,nextPoint.direction,overshoot);
	source.rho=Weigh(previousPoint.rho,nextPoint.rho,overshoot);
	
	static double last_critical_energy,last_Bfactor,last_Bfactor_power; //to speed up calculation if critical energy didn't change
	static double last_average_ans;
	double average_;
	double emittance_x_,emittance_y_;
	double B_factor,B_factor_power;
	double x_offset,y_offset,divx_offset,divy_offset;
	double radius=1E30;
	GenPhoton result;
	Vector X_local,Y_local,Z_local;

	//Calculate local base vectors
	Z_local=source.direction.Normalize(); //Z' base vector
	/*if (source->rho.Norme()<1E29) {
		X_local=source->rho.Normalize(); //X' base vector
	} else {
		X_local=RandomPerpendicularVector(Z_local,1.0);
		}*/
	Y_local=Vector(0.0,1.0,0.0); //same as absolute Y - assuming that machine's orbit is in the XZ plane
	X_local=Crossproduct(Z_local,Y_local);

	double critical_energy;
	if (current_region->emittance>0.0) { //if beam is not ideal
		//calculate non-ideal beam's offset (the four sigmas)
		double betax_,betay_,eta_,etaprime_,e_spread_;
		if (current_region->betax<0.0) { //negative betax value: load BXY file
			double coordinate; //interpolation X value (first column of BXY file)
			if (current_region->beta_kind==0) coordinate=sourceId*current_region->dL;
			else if (current_region->beta_kind==1) coordinate=source.position.x;
			else if (current_region->beta_kind==2) coordinate=source.position.y;
			else if (current_region->beta_kind==3) coordinate=source.position.z;

			betax_=current_region->beta_x_distr->InterpolateY(coordinate);
			betay_=current_region->beta_y_distr->InterpolateY(coordinate);
			eta_=current_region->eta_distr->InterpolateY(coordinate);
			etaprime_=current_region->etaprime_distr->InterpolateY(coordinate);
			e_spread_=current_region->e_spread_distr->InterpolateY(coordinate);
			//the six above distributions are the ones that are read from the BXY file
			//interpolateY finds the Y value corresponding to X input

		} else { //if no BXY file, use average values
			betax_=current_region->betax;
			betay_=current_region->betay;
			eta_=current_region->eta;
			etaprime_=current_region->etaprime;
			e_spread_=current_region->energy_spread;
		}

		emittance_x_=current_region->emittance/(1.0+current_region->coupling);
		emittance_y_=emittance_x_*current_region->coupling;

		double sigmaxprime=sqrt(emittance_x_/betax_+Sqr(etaprime_*e_spread_));
		//{ hor lattice-dependent divergence, radians }
		double sigmayprime=sqrt(emittance_y_/betay_);
		//{ same for vertical divergence, radians }
		double sigmax=sqrt(emittance_x_*betax_+Sqr(eta_*e_spread_));
		//{ hor lattice-dependent beam dimension, cm }
		double sigmay=sqrt(emittance_y_*betay_);
		//{ same for vertical, cm }

		x_offset=Gaussian(sigmax);
		y_offset=Gaussian(sigmay);
		divx_offset=Gaussian(sigmaxprime);
		divy_offset=Gaussian(sigmayprime);

		Vector offset=Vector(0,0,0); //choose ideal beam as origin
		offset=Add(offset,ScalarMult(X_local,x_offset)); //apply dX offset
		offset=Add(offset,ScalarMult(Y_local,y_offset)); //apply dY offset
		result.start_pos=Add(source.position,offset);

		Vector B_local=current_region->B(sourceId,offset); //recalculate B at offset position
		double B_parallel=DotProduct(source.direction,B_local);
		double B_orthogonal=sqrt(Sqr(B_local.Norme())-Sqr(B_parallel));
		
		if (B_orthogonal>VERY_SMALL)
			radius=current_region->E/0.00299792458/B_orthogonal; //Energy in GeV divided by speed of light/1E9 converted to centimeters
		else radius=1.0E30;

		critical_energy=2.959E-5*pow(current_region->gamma,3)/radius; //becomes ~1E-30 if radius is 1E30
	} else {
		//0 emittance, ideal beam
		critical_energy=source.Critical_Energy(current_region->gamma);
		result.start_pos=source.position;
		radius=source.rho.Norme();
	}

	double generated_energy=SYNGEN1(current_region->energy_low/critical_energy,current_region->energy_hi/critical_energy,
		sHandle->generation_mode);

	if (critical_energy==last_critical_energy)
		B_factor=last_Bfactor;
	else B_factor=(exp(integral_N_photons.InterpolateY(log(current_region->energy_hi/critical_energy)))
		-exp(integral_N_photons.InterpolateY(log(current_region->energy_low/critical_energy))))
		/exp(integral_N_photons.valuesY[integral_N_photons.size-1]);

	double SR_flux,SR_power;
	SR_flux=current_region->dL/(radius*2*PI)*current_region->gamma*4.1289E14*B_factor*current_region->current;
	//Total flux per revolution for electrons: 8.08E17*E[GeV]*I[mA] photons/sec
	//8.08E17 * 0.000511GeV = 4.1289E14

	//if (generation_mode==SYNGEN_MODE_POWERWISE) {
	if (critical_energy==last_critical_energy)
		B_factor_power=last_Bfactor_power;
	else B_factor_power=(exp(integral_SR_power.InterpolateY(log(current_region->energy_hi/critical_energy)))
		-exp(integral_SR_power.InterpolateY(log(current_region->energy_low/critical_energy))))
		/exp(integral_SR_power.valuesY[integral_SR_power.size-1]);

	if (critical_energy==last_critical_energy)
		average_=last_average_ans;
	else average_=integral_N_photons.Interval_Mean(current_region->energy_low/critical_energy,current_region->energy_hi/critical_energy);

	double average=integral_N_photons.average;
	//double average1=integral_N_photons.average1; //unused

	if (sHandle->generation_mode==SYNGEN_MODE_POWERWISE)
		SR_flux=SR_flux/generated_energy*average_;

	//}
	double f=polarization_distribution.InterpolateY(generated_energy);
	double g1h2=exp(g1h2_distribution.InterpolateY(generated_energy));
	double f_times_g1h2=f*g1h2;
	double natural_divx,natural_divy;
	
	do {
		natural_divy=find_psi(generated_energy,Sqr(current_region->gamma),f_times_g1h2,
			current_region->enable_par_polarization,current_region->enable_ort_polarization)/current_region->gamma;
	} while (natural_divy>current_region->psimaxY);
	do {
		natural_divx=find_chi(natural_divy,current_region->gamma,f_times_g1h2,
			current_region->enable_par_polarization,current_region->enable_ort_polarization); //divided by sHandle->gamma inside the function
	} while (natural_divx>current_region->psimaxX);

	if (rnd()<0.5) natural_divx=-natural_divx;
	if (rnd()<0.5) natural_divy=-natural_divy;

	if (B_factor>0.0 && average_>VERY_SMALL) {
		SR_power=SR_flux*generated_energy*critical_energy*1.602189E-19*average/average_*B_factor_power/B_factor; //flux already multiplied by current
	} else SR_power=0.0;

	last_critical_energy=critical_energy;
	last_Bfactor=B_factor;
	last_Bfactor_power=B_factor_power;
	last_average_ans=average_;

	//return values
	
	if (current_region->emittance>0.0) {
		//do the transf
		result.start_dir=Z_local;
		result.start_dir=result.start_dir.Rotate(Y_local,natural_divx+divx_offset);
		result.start_dir=result.start_dir.Rotate(X_local,natural_divy+divy_offset);
	} else {
		result.start_dir=Z_local;
		result.start_dir=result.start_dir.Rotate(Y_local,natural_divx);
		result.start_dir=result.start_dir.Rotate(X_local,natural_divy);
	}
	result.dF=SR_flux;
	//if (!(SR_flux==SR_flux)) __debugbreak();
	result.dP=SR_power;
	result.energy=generated_energy*critical_energy;

	return result;
}

double Gaussian(const double &sigma) {
	//Generates a random number following the Gaussian distribution around 0 with 'sigma' standard deviation
	double v1,v2,r,fac;
	do {
		v1=2.0*rnd()-1.0;
		v2=2.0*rnd()-1.0;
		r=Sqr(v1)+Sqr(v2);
	} while (r>=1.0);
	fac=sqrt(-2.0*log(r)/r);
	return v2*fac*sigma;
}

int Stick(FACET* collidedFacet) {
	collidedFacet->sh.counter.nbAbsorbed++;
	collidedFacet->sh.counter.fluxAbs+=sHandle->dF;
	collidedFacet->sh.counter.powerAbs+=sHandle->dP;
	sHandle->tmpCount.nbAbsorbed++;
	sHandle->tmpCount.fluxAbs+=sHandle->dF;
	sHandle->tmpCount.powerAbs+=sHandle->dP;
	sHandle->distTraveledSinceUpdate+=sHandle->distTraveledCurrentParticle;
	//sHandle->counter.nbAbsorbed++;
	//sHandle->counter.fluxAbs+=sHandle->dF;
	//sHandle->counter.powerAbs+=sHandle->dP;
	RecordHit(HIT_ABS,sHandle->dF,sHandle->dP); //for hits and lines display
	if( collidedFacet->hits_MC && collidedFacet->sh.countAbs ) RecordHitOnTexture(collidedFacet,sHandle->dF,sHandle->dP);
	if( !StartFromSource() )
		// maxDesorption reached
		return 0; //FALSE
	return 1;
}