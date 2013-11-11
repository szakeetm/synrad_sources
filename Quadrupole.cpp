#include "Region.h"
#include "Random.h"
#include "Tools.h"
#include <ctime>
#include <vector>
#include <string>
using namespace std;

Vector Quadrupole::B(const Vector &position) {
	
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

	return Vector(xp*cosbeta_q-yp*sinalfa_q*sinbeta_q -zp*cosalfa_q*sinbeta_q,
		yp*cosalfa_q -zp*sinalfa_q,
		xp*sinbeta_q+yp*sinalfa_q*cosbeta_q +zp*cosalfa_q*cosbeta_q);
}