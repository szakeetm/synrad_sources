/*
  File:        IntersectAABB.cpp
  Description: Ray geometry intersection (Using AABB tree optimisation)
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
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "Simulation.h"
#include "Random.h"
#include "GLApp/MathTools.h"
#include "SynradTypes.h" //Histogram

// -----------------------------------------------------------
// AABB tree stuff
// -----------------------------------------------------------

// Temporary for intersection
extern  double    intMinLgth;
extern  BOOL      intFound;
extern  Vector3d  intD;
extern  Vector3d  intZ;
extern  int       intNbTHits;
extern  double    iRx;
extern  double    iRy;
extern  double    iRz;
extern  BOOL      nullRx;
extern  BOOL      nullRy;
extern  BOOL      nullRz;
extern  Vector3d *rayPos;
extern  Vector3d *rayDir;
extern  FACET   **iFacet;
extern  FACET    *fLast;
extern  double    tNear;
extern  double    tFar;
extern  double    it1, it2;
extern  BOOL      AABBHit;



void ProfileFacet(FACET *f, const double &dF, const double &dP, const double &E) {

	int pos;

	f->hitted = TRUE;

	switch (f->sh.profileType) {

	case REC_ANGULAR: {
		double dot = abs(Dot(f->sh.N, *rayDir));
		double theta = acos(dot);              // Angle to normal (0 to PI/2)
		int grad = (int)(((double)PROFILE_SIZE)*(theta) / (PI / 2)); // To Grad
		SATURATE(grad, 0, PROFILE_SIZE - 1);
		f->profile_hits[grad]++;
		f->profile_flux[grad] += dF;
		f->profile_power[grad] += dP;

	} break;

	case REC_PRESSUREU:
		pos = (int)((f->colU)*(double)PROFILE_SIZE);
		SATURATE(pos, 0, PROFILE_SIZE - 1);
		f->profile_hits[pos]++;
		f->profile_flux[pos] += dF;
		f->profile_power[pos] += dP;
		break;

	case REC_PRESSUREV:
		pos = (int)((f->colV)*(double)PROFILE_SIZE);
		SATURATE(pos, 0, PROFILE_SIZE - 1);
		f->profile_hits[pos]++;
		f->profile_flux[pos] += dF;
		f->profile_power[pos] += dP;
		break;

	}

	if (f->sh.hasSpectrum) {
		f->spectrum_fluxwise->Add(E, dF, 0.001);
		f->spectrum_powerwise->Add(E, dP, 0.001);
	}
}

BOOL Intersect(Vector3d *rPos, Vector3d *rDir,  // Source ray (rayDir vector must be normalized)
	double *dist,                   // Distance to collision point
	FACET **iFact, FACET *last) {    // Collided facet, previous collision

	int i, j;
	intMinLgth = 1e100;
	intFound = FALSE;
	intNbTHits = 0;
	rayPos = rPos;
	rayDir = rDir;
	intD.x = -rayDir->x;
	intD.y = -rayDir->y;
	intD.z = -rayDir->z;
	nullRx = (rayDir->x == 0.0);
	nullRy = (rayDir->y == 0.0);
	nullRz = (rayDir->z == 0.0);
	if (!nullRx) iRx = 1.0 / rayDir->x;
	if (!nullRy) iRy = 1.0 / rayDir->y;
	if (!nullRz) iRz = 1.0 / rayDir->z;
	iFacet = iFact;
	fLast = last;

	IntersectTree(sHandle->str[sHandle->curStruct].aabbTree);

	if (intFound) {

		FACET *f = *iFacet;
		*dist = intMinLgth;

		ProfileFacet(f, sHandle->dF, sHandle->dP, sHandle->energy);

		// Second pass for transparent hits
		for (i = 0; i < intNbTHits; i++) {

			f = THits[i];
			if (f->colDist < intMinLgth) {
				f->counter.nbHit++; //count MC hits
				f->counter.fluxAbs += sHandle->dF;
				f->counter.powerAbs += sHandle->dP;
				ProfileFacet(f, sHandle->dF, sHandle->dP, sHandle->energy); //count profile
				if (f->hits_MC && f->sh.countTrans) RecordHitOnTexture(f, sHandle->dF, sHandle->dP); //count texture
				if (f->direction && f->sh.countDirection) RecordDirectionVector(f);
			}

		}

		/*
		// Compute intersection with spheric volume element
		if (sHandle->hasDirection) {

			for (j = 0; j < sHandle->nbSuper; j++) {
				for (i = 0; i < sHandle->str[j].nbFacet; i++) {
					f = sHandle->str[j].facets[i];
					if (f->direction && f->sh.countDirection) {

						int      x, y;
						Vector3d center;
						double   d;
						double   r = f->rw*0.45; // rw/2 - 10% (avoid side FX)

						for (x = 0; x < f->sh.texWidth; x++) {
							for (y = 0; y < f->sh.texHeight; y++) {
								int add = x + y*f->sh.texWidth;
								if (isFull) {

									double uC = ((double)x + 0.5) * f->iw;
									double vC = ((double)y + 0.5) * f->ih;
									center.x = f->sh.O.x + f->sh.U.x*uC + f->sh.V.x*vC;
									center.y = f->sh.O.y + f->sh.U.y*uC + f->sh.V.y*vC;
									center.z = f->sh.O.z + f->sh.U.z*uC + f->sh.V.z*vC;
									if (RaySphereIntersect(&center, r, rPos, rDir, &d)) {
										if (d < intMinLgth) {
											f->direction[add].dir.x += sHandle->pDir.x;
											f->direction[add].dir.y += sHandle->pDir.y;
											f->direction[add].dir.z += sHandle->pDir.z;
											f->direction[add].count++;
										}
									}

								}
							}
						}
					}
				}
			}
		}*/

	}

	return intFound;

}
