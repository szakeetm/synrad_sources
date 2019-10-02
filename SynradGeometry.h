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

#include "Geometry_shared.h"
#include "Region_full.h"
#include <cereal/archives/json.hpp>

#define SYNVERSION   10

#define PARAMVERSION 4
class Worker;
class Material;

class SynradGeometry: public Geometry {

#pragma region Geometry.cpp
private:
	//void CalculateFacetParam(int facetId); // Facet parameters
	std::vector<std::string> InsertSYNGeom(FileReader *file, size_t strIdx, bool newStruct = false);
	void SaveProfileSYN(FileWriter *file, Dataport *dpHit, int super = -1, bool saveSelected = false, bool crashSave = false);
	void SaveSpectrumSYN(FileWriter *file, Dataport *dpHit, int super = -1, bool saveSelected = false, bool crashSave = false);
	//void SaveProfileGEO(FileWriter *file, int super = -1, bool saveSelected = false);
public:
	SynradGeometry();
	void LoadGEO(FileReader *file, GLProgress *prg, int *version, Worker *worker);
	std::vector<std::string> LoadSYN(FileReader *file, GLProgress *prg, int *version, Worker *worker);
	bool LoadTextures(FileReader *file, GLProgress *prg, Dataport *dpHit, int version);
	std::vector<std::string> InsertSYN(FileReader *file, GLProgress *prg, bool newStr);
	void SaveTXT(FileWriter *file, Dataport *dhHit, bool saveSelected);
	void ExportTextures(FILE *file, int grouping, int mode, double no_scans, Dataport *dhHit, bool saveSelected);
	void SaveDesorption(FILE *file, Dataport *dhHit, bool selectedOnly, int mode, double eta0, double alpha, const Distribution2D &distr); //Deprecated, not used anymore

	//void SaveGEO(FileWriter *file,GLProgress *prg,Dataport *dpHit,bool saveSelected,LEAK *pleak,int *nbleakSave,HIT *hitCache,int *nbHHitSave,bool crashSave=false);
	void SaveSYN(FileWriter *file, GLProgress *prg, Dataport *dpHit, bool saveSelected, LEAK *leakCache, size_t *nbLeakTotal, HIT *hitCache, size_t *nbHitSave, bool crashSave = false);
	void SaveXML_geometry(pugi::xml_node saveDoc, Worker *work, GLProgress *prg, bool saveSelected);
	bool SaveXML_simustate(pugi::xml_node saveDoc, Worker *work, BYTE *buffer, GlobalHitBuffer *gHits, int nbLeakSave, int nbHHitSave,
		LEAK *leakCache, HIT *hitCache, GLProgress *prg, bool saveSelected);
	void LoadXML_geom(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg);
	void InsertXML(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg, bool newStr);
	bool LoadXML_simustate(pugi::xml_node loadXML, Dataport *dpHit, Worker *work, GLProgress *progressDlg);
	void BuildPipe(double L, double R, double s, int step);
	void LoadProfileSYN(FileReader *file, Dataport *dpHit, const int& version);
	void LoadSpectrumSYN(FileReader *file, Dataport *dpHit, const int& version);
	size_t GetGeometrySize(std::vector<Region_full> &regions, std::vector<Material> &materials, 
		std::vector<std::vector<double>> &psi_distro, std::vector<std::vector<std::vector<double>>> &chi_distros,
		const std::vector<std::vector<double>>& parallel_polarization);
	size_t GetHitsSize();
	void CopyGeometryBuffer(BYTE *buffer, std::vector<Region_full> &regions, std::vector<Material> &materials,
		std::vector<std::vector<double>> &psi_distro, const std::vector<std::vector<std::vector<double>>> &chi_distros,
		const std::vector<std::vector<double>> &parallel_polarization, const bool& newReflectionModel, const OntheflySimulationParams& ontheflyParams);
#pragma endregion

#pragma region GeometryRender.cpp
	void BuildFacetTextures(BYTE *hits, bool renderRegularTexture, bool renderDirectionTexture);
#pragma endregion

    void SerializeForLoader(cereal::BinaryOutputArchive &outputArchive);

	// Temporary variable (used by LoadXXX)
	double loaded_totalFlux;
	double loaded_totalPower;
	double loaded_no_scans;

};

