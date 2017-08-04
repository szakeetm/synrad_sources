/*
  File:        Geometry.h
  Description: Main geometry class (Handles sets of facets)
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

#pragma once

#include "Geometry.h"
#include "Region_full.h"

#define TEXTURE_MODE_MCHITS 0
#define TEXTURE_MODE_FLUX 1
#define TEXTURE_MODE_POWER 2
#define SYNVERSION   9

#define PARAMVERSION 3
class Worker;

class SynradGeometry: public Geometry {

#pragma region Geometry.cpp
private:
	//void CalculateFacetParam(int facetId); // Facet parameters
	std::vector<std::string> InsertSYNGeom(FileReader *file, size_t *nbV, size_t *nbF, InterfaceVertex **V, Facet ***F, size_t strIdx = 0, bool newStruct = false);
	void SaveProfileSYN(FileWriter *file, Dataport *dpHit, int super = -1, bool saveSelected = false, bool crashSave = false);
	void SaveSpectrumSYN(FileWriter *file, Dataport *dpHit, int super = -1, bool saveSelected = false, bool crashSave = false);
	void SaveProfileGEO(FileWriter *file, int super = -1, bool saveSelected = false);
public:
	SynradGeometry();
	void LoadGEO(FileReader *file, GLProgress *prg, int *version);
	std::vector<std::string> LoadSYN(FileReader *file, GLProgress *prg, LEAK *pleak, size_t *nbleakLoad, HIT *hitCache, size_t *nbHHitLoad, int *version);
	bool LoadTextures(FileReader *file, GLProgress *prg, Dataport *dpHit, int version);
	std::vector<std::string> InsertSYN(FileReader *file, GLProgress *prg, bool newStr);
	void SaveTXT(FileWriter *file, Dataport *dhHit, bool saveSelected);
	void ExportTextures(FILE *file, int grouping, int mode, double no_scans, Dataport *dhHit, bool saveSelected);
	void SaveDesorption(FILE *file, Dataport *dhHit, bool selectedOnly, int mode, double eta0, double alpha, Distribution2D *distr);

	//void SaveGEO(FileWriter *file,GLProgress *prg,Dataport *dpHit,bool saveSelected,LEAK *pleak,int *nbleakSave,HIT *hitCache,int *nbHHitSave,bool crashSave=false);
	void SaveSYN(FileWriter *file, GLProgress *prg, Dataport *dpHit, bool saveSelected, LEAK *leakCache, size_t *nbLeakTotal, HIT *hitCache, size_t *nbHitSave, bool crashSave = false);
	void SaveXML_geometry(pugi::xml_node saveDoc, Worker *work, GLProgress *prg, bool saveSelected);
	bool SaveXML_simustate(pugi::xml_node saveDoc, Worker *work, BYTE *buffer, SHGHITS *gHits, int nbLeakSave, int nbHHitSave,
		LEAK *leakCache, HIT *hitCache, GLProgress *prg, bool saveSelected);
	void LoadXML_geom(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg);
	void InsertXML(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg, bool newStr);
	bool LoadXML_simustate(pugi::xml_node loadXML, Dataport *dpHit, Worker *work, GLProgress *progressDlg);
	void BuildPipe(double L, double R, double s, int step);
	void LoadProfileSYN(FileReader *file, Dataport *dpHit);
	void LoadSpectrumSYN(FileReader *file, Dataport *dpHit);
	size_t GetGeometrySize(std::vector<Region_full> *regions, std::vector<Material> *materials, 
		std::vector<std::vector<double>> &psi_distro, std::vector<std::vector<std::vector<double>>> &chi_distros,
		const std::vector<std::vector<double>>& parallel_polarization);
	size_t GetHitsSize();
	void CopyGeometryBuffer(BYTE *buffer, std::vector<Region_full> *regions, std::vector<Material> *materials,
		std::vector<std::vector<double>> &psi_distro, const std::vector<std::vector<std::vector<double>>> &chi_distros,
		const std::vector<std::vector<double>> &parallel_polarization, int generation_mode, bool lowFluxMode, double lowFluxCutoff, bool newReflectionModel);
#pragma endregion

#pragma region GeometryRender.cpp
	void BuildFacetTextures(BYTE *hits, bool renderRegularTexture, bool renderDirectionTexture);
#pragma endregion

	// Temporary variable (used by LoadXXX)
	double loaded_totalFlux;
	double loaded_totalPower;
	double loaded_no_scans;

};

