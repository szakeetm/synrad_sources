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
#pragma once

#include "Region_mathonly.h"
GenPhoton GeneratePhoton(size_t pointId, Region_mathonly *current_region, int generation_mode,
	std::vector<std::vector<double>> &psi_distro, std::vector<std::vector<double>> &chi_distro,
	std::vector<std::vector<double>> &parallel_polarization, bool recalc = 0); //Generates a photon from point number 'pointId'
double Interval_Mean(const double &min, const double &max);