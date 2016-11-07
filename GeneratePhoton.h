#ifndef _GENPHOTONH_
#define _GENPHOTONH_
#include "Tools.h"
#include "Region_mathonly.h"
GenPhoton GeneratePhoton(int pointId, Region_mathonly *current_region, int generation_mode,
	std::vector<std::vector<double>> &psi_distr,std::vector<std::vector<double>> &chi_distr,BOOL recalc = 0); //Generates a photon from point number 'pointId'
double Interval_Mean(const double &min, const double &max);
#endif