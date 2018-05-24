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
#include "Region_mathonly.h"
#include "Random.h"
#include <ctime>
#include <vector>
#include <string>
using namespace std;

Vector3d Quadrupole::B(const Vector3d &position) {
	
	double dX,dY,dZ; //offsets from quadrupole center
	double Bx_,By_,Bz_; //magnetic field in the transformed coordinates
	
	//Calculate offset
	dX=position.x-center.x;
	dY=position.y-center.y;
	dZ=position.z-center.z;

	//Transform 
	double xp=dX*cosbeta_q+dZ*sinbeta_q;
	double yp=-dX*sinalfa_q*sinbeta_q+dY*cosalfa_q+dZ*sinalfa_q*cosbeta_q;
	double zp=-dX*cosalfa_q*sinbeta_q-dY*sinalfa_q+dZ*cosalfa_q*cosbeta_q;

	dX=xp*cosrot_q+yp*sinrot_q;
	dY=-xp*sinrot_q+yp*cosrot_q;
	xp=dX;
	yp=dY;

	if (zp>=0.0 && zp<=L_q) { //if inside the quadrupole
		Bx_=-K_q*yp;
		By_=-K_q*xp;
	} else {
		Bx_=0.0;
		By_=0.0;
	}
	Bz_=0.0;

	//Inverse transformation
	xp=Bx_*cosrot_q-By_*sinrot_q;
	yp=Bx_*sinrot_q+By_*cosrot_q;
	zp=Bz_;

	return Vector3d(xp*cosbeta_q-yp*sinalfa_q*sinbeta_q -zp*cosalfa_q*sinbeta_q,
		yp*cosalfa_q -zp*sinalfa_q,
		xp*sinbeta_q+yp*sinalfa_q*cosbeta_q +zp*cosalfa_q*cosbeta_q);
}