#ifndef _QUADRUPOLE_
#define _QUADRUPOLE_

#include <math.h>
#include "Tools.h"

class Quadrupole {
public:
	Vector center,direction;
	double alfa_q,beta_q,rot_q;
	double sinalfa_q,cosalfa_q,sinbeta_q,cosbeta_q,sinrot_q,cosrot_q;
	double K_q,L_q,period_q;
	Vector B(const Vector &position);
};

#endif