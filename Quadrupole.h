#ifndef _QUADRUPOLE_
#define _QUADRUPOLE_

class Quadrupole {
public:
	Vector3d center,direction;
	double alfa_q,beta_q,rot_q;
	double sinalfa_q,cosalfa_q,sinbeta_q,cosbeta_q,sinrot_q,cosrot_q;
	double K_q,L_q,period_q;
	Vector3d B(const Vector3d &position);
};

#endif