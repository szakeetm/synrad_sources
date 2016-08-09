/*
File:        Geometry.cpp
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

#include "Geometry.h"
#include "ASELoader.h"
#include "Utils.h"
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "GLApp/GLMatrix.h"
#include "GLApp\GLMessageBox.h"
#include "Synrad.h"
#include "GLApp\GLWindowManager.h"
#include "Region_full.h"

#define WRITEBUFFER(_value,_type) *((_type *)buffer)=_value;buffer += sizeof(_type)

using namespace pugi;
extern SynRad *mApp;

Geometry::Geometry() {

	facets = NULL;
	vertices3 = NULL;
	polyList = 0;
	selectList = 0;
	selectList2 = 0;
	selectList3 = 0;
	arrowList = 0;
	sphereList = 0;
	textureMode = TEXTURE_MODE_FLUX;
	texMin_MC = 0;
	texMax_MC = 1;
	texMin_flux = 0.0;
	texMax_flux = 1.0;
	texMin_power = 0.0;
	texMax_power = 1.0;
	autoNorme = TRUE;
	centerNorme = TRUE;
	normeRatio = 1.0f;
	texAutoScale = TRUE;
	texLogScale = FALSE;

	sh.nbSuper = 0;
	viewStruct = -1;
	strcpy(strPath, "");
	Clear();

}

Geometry::~Geometry() {
	Clear();
}

void Geometry::Clear() {
	viewStruct = -1; //otherwise a nonexistent structure could stay selected
	// Free memory
	if (facets) {
		for (int i = 0; i < sh.nbFacet; i++)
			SAFE_DELETE(facets[i]);
		free(facets);
	}
	if (vertices3) free(vertices3);
	for (int i = 0; i < sh.nbSuper; i++) {
		SAFE_FREE(strName[i]);
		SAFE_FREE(strFileName[i]);
	}
	memset(strName, 0, MAX_SUPERSTR * sizeof(char *));
	memset(strFileName, 0, MAX_SUPERSTR * sizeof(char *));
	DeleteGLLists(TRUE, TRUE);

	if (mApp && mApp->splitFacet) mApp->splitFacet->ClearUndoFacets();


	// Init default
	facets = NULL;         // Facets array
	vertices3 = NULL;      // Facets vertices in (x,y,z) space
	sh.nbFacet = 0;        // Number of facets
	sh.nbVertex = 0;       // Number of vertex
	isLoaded = FALSE;      // isLoaded flag
	texCMin_MC = 0;             // Current minimum
	texCMax_MC = 1;             // Current maximum
	texCMin_flux = texCMin_power = 0.0;
	texCMax_flux = texCMax_power = 1.0;
	sh.nbSuper = 0;          // Structure number
	nbSelected = 0;          // Number of selected facet

	memset(lineList, 0, sizeof(lineList));
	memset(strName, 0, sizeof(strName));
	memset(strFileName, 0, sizeof(strFileName));

	// Init OpenGL material
	memset(&whiteMaterial, 0, sizeof(GLMATERIAL));
	whiteMaterial.Diffuse.r = 0.9f;
	whiteMaterial.Diffuse.g = 0.9f;
	whiteMaterial.Diffuse.b = 0.9f;
	whiteMaterial.Ambient.r = 0.9f;
	whiteMaterial.Ambient.g = 0.9f;
	whiteMaterial.Ambient.b = 0.9f;

	memset(&fillMaterial, 0, sizeof(GLMATERIAL));


	fillMaterial.Diffuse.r = 0.6f;
	fillMaterial.Diffuse.g = 0.65f;
	fillMaterial.Diffuse.b = 0.65f;
	fillMaterial.Ambient.r = 0.45f;
	fillMaterial.Ambient.g = 0.41f;
	fillMaterial.Ambient.b = 0.41f;

	memset(&arrowMaterial, 0, sizeof(GLMATERIAL));
	arrowMaterial.Diffuse.r = 0.4f;
	arrowMaterial.Diffuse.g = 0.2f;
	arrowMaterial.Diffuse.b = 0.0f;
	arrowMaterial.Ambient.r = 0.6f;
	arrowMaterial.Ambient.g = 0.3f;
	arrowMaterial.Ambient.b = 0.0f;

	nbSelectedHist = 0;    // Selection history
	nbSelectedHistVertex = 0;

}

void Geometry::CalculateFacetParam(int facetId) {

	CalculateFacetParam_geometry(facets[facetId]);

}

void Geometry::InitializeGeometry(int facet_number, BOOL BBOnly) {

	// Perform various precalculation here for a faster simulation.

	//GLProgress *initGeoPrg = new GLProgress("Initializing geometry...","Please wait");
	//initGeoPrg->SetProgress(0.0);
	//initGeoPrg->SetVisible(TRUE);
	if (facet_number == -1) { //bounding box for all vertices
		bb.min.x = 1e100;
		bb.min.y = 1e100;
		bb.min.z = 1e100;
		bb.max.x = -1e100;
		bb.max.y = -1e100;
		bb.max.z = -1e100;

		// Axis Aligned Bounding Box
		for (int i = 0; i < sh.nbVertex; i++) {
			VERTEX3D p = vertices3[i];
			if (!(vertices3[i].selected == FALSE || vertices3[i].selected == TRUE)) vertices3[i].selected = FALSE; //initialize selection
			if (p.x < bb.min.x) bb.min.x = p.x;
			if (p.y < bb.min.y) bb.min.y = p.y;
			if (p.z < bb.min.z) bb.min.z = p.z;
			if (p.x > bb.max.x) bb.max.x = p.x;
			if (p.y > bb.max.y) bb.max.y = p.y;
			if (p.z > bb.max.z) bb.max.z = p.z;
		}

		Worker *worker = &(mApp->worker);
		for (int i = 0; i < (int)worker->regions.size(); i++) {
			if (worker->regions[i].AABBmin.x < bb.min.x) bb.min.x = worker->regions[i].AABBmin.x;
			if (worker->regions[i].AABBmin.y < bb.min.y) bb.min.y = worker->regions[i].AABBmin.y;
			if (worker->regions[i].AABBmin.z < bb.min.z) bb.min.z = worker->regions[i].AABBmin.z;
			if (worker->regions[i].AABBmax.x > bb.max.x) bb.max.x = worker->regions[i].AABBmax.x;
			if (worker->regions[i].AABBmax.y > bb.max.y) bb.max.y = worker->regions[i].AABBmax.y;
			if (worker->regions[i].AABBmax.z > bb.max.z) bb.max.z = worker->regions[i].AABBmax.z;
		}
	}
	else { //bounding box only for the changed facet
		for (int i = 0; i < facets[facet_number]->sh.nbIndex; i++) {
			VERTEX3D p = vertices3[facets[facet_number]->indices[i]];
			//if(!(vertices3[i].selected==FALSE || vertices3[i].selected==TRUE)) vertices3[i].selected=FALSE; //initialize selection
			if (p.x < bb.min.x) bb.min.x = p.x;
			if (p.y < bb.min.y) bb.min.y = p.y;
			if (p.z < bb.min.z) bb.min.z = p.z;
			if (p.x > bb.max.x) bb.max.x = p.x;
			if (p.y > bb.max.y) bb.max.y = p.y;
			if (p.z > bb.max.z) bb.max.z = p.z;
		}
	}

	center.x = (bb.max.x + bb.min.x) / 2.0;
	center.y = (bb.max.y + bb.min.y) / 2.0;
	center.z = (bb.max.z + bb.min.z) / 2.0;

	if (BBOnly) return;


	// Choose an orthogonal (U,V) 2D basis for each facet. (This method can be 
	// improved. stub). The algorithm chooses the longest vedge for the U vector.
	// then it computes V (orthogonal to U and N). Afterwards, U and V are rescaled 
	// so each facet vertex are included in the rectangle defined by uU + vV (0<=u<=1 
	// and 0<=v<=1) of origin f->sh.O, U and V are always orthogonal and (U,V,N) 
	// form a 3D left handed orthogonal basis (not necessary orthonormal).
	// This coordinates system allows to prevent from possible "almost degenerated"
	// basis on fine geometry. It also greatly eases the facet/ray instersection routine 
	// and ref/abs/des hit recording and visualization. In order to ease some calculations, 
	// nU et nV (normalized U et V) are also stored in the Facet structure.
	// The local coordinates of facet vertex are stored in (U,V) coordinates (vertices2).



	int fOffset = sizeof(SHGHITS);
	for (int i = 0; i < sh.nbFacet; i++) {
		//initGeoPrg->SetProgress((double)i/(double)sh.nbFacet);
		if ((facet_number == -1) || (i == facet_number)) { //permits to initialize only one facet
			// Main facet params
			CalculateFacetParam(i);

			// Current facet
			Facet *f = facets[i];

			// Detect non visible edge
			f->InitVisibleEdge();

			// Detect orientation
			f->DetectOrientation();

			if (facet_number == -1) {
				// Hit address
				f->sh.hitOffset = fOffset;
				fOffset += f->GetHitsSize();
			}
		}
	}

	isLoaded = TRUE;
	if (facet_number == -1) {
		BuildGLList();
		mApp->UpdateModelParams();
		mApp->UpdateFacetParams();
	}
	//initGeoPrg->SetVisible(FALSE);
	//SAFE_DELETE(initGeoPrg);
	_ASSERTE(_CrtCheckMemory());
}

void Geometry::RebuildLists() {
	BuildGLList();
}

size_t Geometry::GetGeometrySize(std::vector<Region_full> *regions, std::vector<Material> *materials, std::vector<std::vector<double>> &psi_distr, std::vector<std::vector<double>> &chi_distr) {

	// Compute number of bytes allocated
	size_t memoryUsage = 0;
	memoryUsage += sizeof(SHGEOM);
	
	//Regions
	memoryUsage += (*regions).size()*sizeof(RegionParams);
	for (int i = 0; i < (int)(*regions).size(); i++) {
		memoryUsage += sizeof(Trajectory_Point)*(*regions)[i].Points.size();
		memoryUsage += 2 * sizeof(double)*(*regions)[i].Bx_distr.size;
		memoryUsage += 2 * sizeof(double)*(*regions)[i].By_distr.size;
		memoryUsage += 2 * sizeof(double)*(*regions)[i].Bz_distr.size;
		memoryUsage += 6 * sizeof(double)*(*regions)[i].beta_x_distr.size; //Coord, BetaX, BetaY, EtaX, EtaX' , E_spread
	}
	//Material library
	memoryUsage += sizeof(size_t); //number of (*materials)
	for (size_t i = 0; i < (*materials).size(); i++) { //going through all (*materials)
		memoryUsage += sizeof(BOOL);//hasBackscattering
		memoryUsage += sizeof(size_t);//copying number of angles (columns)
		memoryUsage += sizeof(size_t);//copying number of energies (rows)
		memoryUsage += ((*materials)[i].angleVals.size())*sizeof(double);//copying angles (header)
		memoryUsage += ((*materials)[i].energyVals.size())*sizeof(double);//copying energies (column1)
		size_t nbComponents = (*materials)[i].hasBackscattering ? 4 : 1;
		memoryUsage += ((*materials)[i].angleVals.size())*(*materials)[i].energyVals.size()*nbComponents*sizeof(double);//copying reflectivity probabilities (cells)
	}

	memoryUsage += sizeof(size_t); //psi number of rows
	if (psi_distr.size() > 0) memoryUsage += psi_distr.size()*(sizeof(size_t) + psi_distr[0].size()*sizeof(double));

	memoryUsage += sizeof(size_t); //chi number of rows
	if (chi_distr.size() > 0) memoryUsage += chi_distr.size()*(sizeof(size_t) + chi_distr[0].size()*sizeof(double));

	memoryUsage += sh.nbVertex * sizeof(VERTEX3D);
	
	for (int i = 0; i < sh.nbFacet; i++) {
		memoryUsage += facets[i]->GetGeometrySize();
	}

	return memoryUsage;
}

void Geometry::CopyGeometryBuffer(BYTE *buffer, std::vector<Region_full> *regions, std::vector<Material> *materials,
	std::vector<std::vector<double>> &psi_distr, std::vector<std::vector<double>> &chi_distr, int generation_mode, BOOL lowFluxMode, double lowFluxCutoff) {

	// Build shared buffer for geometry (see Shared.h)
	int fOffset = sizeof(SHGHITS);
	SHGEOM *shGeom = (SHGEOM *)buffer;
	sh.nbRegion = (*regions).size();
	sh.generation_mode = generation_mode;
	sh.lowFluxMode = lowFluxMode;
	sh.lowFluxCutoff = lowFluxCutoff;
	memcpy(shGeom, &(this->sh), sizeof(SHGEOM));
	buffer += sizeof(SHGEOM);

	// Build shared buffer for trajectory (see Shared.h)
	for (size_t i = 0; i < sh.nbRegion; i++) {
		RegionParams* regparam = (RegionParams*)buffer;
		(*regions)[i].params.nbDistr_MAG = Vector((*regions)[i].Bx_distr.size, (*regions)[i].By_distr.size, (*regions)[i].Bz_distr.size);
		(*regions)[i].params.nbPointsToCopy = (*regions)[i].Points.size();
		memcpy(regparam, &((*regions)[i].params),sizeof(RegionParams));
		//reg = regions[i];
		buffer += sizeof(RegionParams);
	}
	//copy trajectory points
	for (size_t i = 0; i < sh.nbRegion; i++) {
		/*for (size_t j = 0; j < (*regions)[i].params.nbPointsToCopy; j++){
			WRITEBUFFER((*regions)[i].Points[j], Trajectory_Point);
		}*/
		
		memcpy(buffer, &(*regions)[i].Points[0], (*regions)[i].Points.size()  * sizeof(Trajectory_Point));
		buffer += (*regions)[i].Points.size()  * sizeof(Trajectory_Point);
	}

	//copy distribution points
	for (size_t i = 0; i < sh.nbRegion; i++) {
		for (size_t j = 0; j < (size_t)(*regions)[i].params.nbDistr_MAG.x; j++) {
			WRITEBUFFER((*regions)[i].Bx_distr.valuesX[j], double);
			WRITEBUFFER((*regions)[i].Bx_distr.valuesY[j], double);
		}

		for (size_t j = 0; j < (size_t)(*regions)[i].params.nbDistr_MAG.y; j++) {
			WRITEBUFFER((*regions)[i].By_distr.valuesX[j], double);
			WRITEBUFFER((*regions)[i].By_distr.valuesY[j], double);
		}

		for (size_t j = 0; j < (size_t)(*regions)[i].params.nbDistr_MAG.z; j++) {
			WRITEBUFFER((*regions)[i].Bz_distr.valuesX[j], double);
			WRITEBUFFER((*regions)[i].Bz_distr.valuesY[j], double);
		}

		for (size_t j = 0; j < (*regions)[i].params.nbDistr_BXY; j++) {
			WRITEBUFFER((*regions)[i].beta_x_distr.valuesX[j], double);
			WRITEBUFFER((*regions)[i].beta_x_distr.valuesY[j], double);
			WRITEBUFFER((*regions)[i].beta_y_distr.valuesY[j], double);
			WRITEBUFFER((*regions)[i].eta_distr.valuesY[j], double);
			WRITEBUFFER((*regions)[i].etaprime_distr.valuesY[j], double);
			WRITEBUFFER((*regions)[i].e_spread_distr.valuesY[j], double);
		}
	}

	WRITEBUFFER((*materials).size(), size_t); //copying number of materials

	for (auto material : *materials) { //going through all materials
		WRITEBUFFER(material.hasBackscattering, BOOL);
		WRITEBUFFER(material.angleVals.size(), size_t); //copying number of angles (columns)
		WRITEBUFFER(material.energyVals.size(), size_t); //copying number of energies (rows)

		for (auto angleval : material.angleVals){
			WRITEBUFFER(angleval, double);//copying angles (header)
		}
		for (auto energyval : material.energyVals){
			WRITEBUFFER(energyval, double); //copying energies (column1)
		}

		for (size_t j = 0; j < material.energyVals.size(); j++) {
			for (size_t k = 0; k < material.angleVals.size(); k++) {
				WRITEBUFFER(material.reflVals[j][k][0], double); //forward reflection probability
				if (material.hasBackscattering) {
					WRITEBUFFER(material.reflVals[j][k][1], double); //diffuse reflection probability
					WRITEBUFFER(material.reflVals[j][k][2], double); //back reflection probability
					WRITEBUFFER(material.reflVals[j][k][3], double); //transparent pass probability
				}
			}
		}
	}

	//psi_distr
	WRITEBUFFER(psi_distr.size(), size_t); //copying number of rows
	for (auto row : psi_distr) {
		WRITEBUFFER(row.size(), size_t); //copying number of values
		for (auto value : row){
			WRITEBUFFER(value, double);
		}
	}

	//chi_distr
	WRITEBUFFER(chi_distr.size(), size_t); //copying number of rows
	for (auto row : chi_distr) {
		WRITEBUFFER(row.size(), size_t); //copying number of values
		for (auto value : row) {
			WRITEBUFFER(value, double);
		}
	}


	memcpy(buffer, vertices3, sizeof(VERTEX3D)*sh.nbVertex);
	buffer += sizeof(VERTEX3D)*sh.nbVertex;
	for (size_t i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		f->sh.hitOffset = fOffset;
		fOffset += f->GetHitsSize();
		memcpy(buffer, &(f->sh), sizeof(SHFACET));
		buffer += sizeof(SHFACET);
		memcpy(buffer, f->indices, sizeof(int)*f->sh.nbIndex);
		buffer += sizeof(int)*f->sh.nbIndex;
		memcpy(buffer, f->vertices2, sizeof(VERTEX2D)*f->sh.nbIndex);
		buffer += sizeof(VERTEX2D)*f->sh.nbIndex;
	}

	// Add surface elements area (reciprocal)
	for (int k = 0; k < sh.nbFacet; k++) {
		Facet *f = facets[k];
		DWORD add = 0;
		if (f->sh.isTextured) {

			if (f->cellPropertiesIds) {

				for (int j = 0; j < f->sh.texHeight; j++) {
					for (int i = 0; i < f->sh.texWidth; i++) {
						float area = f->GetMeshArea(add);
						if (area > 0.0f) {
							// Use the sign bit to store isFull flag
							/*if (f->mesh[add].full) {
								WRITEBUFFER(-1.0 / area, double);
							}
							else {*/
								WRITEBUFFER(1.0 / area, double);
							/*}*/
						}
						else {
							WRITEBUFFER(0.0, double);
						}
						add++;
					}
				}

			}
			else {

				double rw = Norme(f->sh.U) / (double)(f->sh.texWidthD);
				double rh = Norme(f->sh.V) / (double)(f->sh.texHeightD);
				float area = (float)(rw*rh);

				for (int j = 0; j < f->sh.texHeight; j++) {
					for (int i = 0; i < f->sh.texWidth; i++) {
						if (area > 0.0f) {
							WRITEBUFFER(1.0 / area, double);
						}
						else {
							WRITEBUFFER(0.0, double);
						}
					}
				}

			}
		}
	}

}


void Geometry::SetAutoNorme(BOOL enable) {
	autoNorme = enable;
}

BOOL Geometry::GetAutoNorme() {
	return autoNorme;
}

void Geometry::SetCenterNorme(BOOL enable) {
	centerNorme = enable;
}

BOOL Geometry::GetCenterNorme() {
	return centerNorme;
}

void Geometry::SetNormeRatio(float r) {
	normeRatio = r;
}

float Geometry::GetNormeRatio() {
	return normeRatio;
}

DWORD Geometry::GetHitsSize() {

	// Compute number of bytes allocated
	DWORD memoryUsage = 0;
	memoryUsage += sizeof(SHGHITS);
	for (int i = 0; i < sh.nbFacet; i++) {
		memoryUsage += facets[i]->GetHitsSize();
	}

	return memoryUsage;
}

/*DWORD Geometry::GetMaxElemNumber() {

	int nbElem = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		if (f->mesh) nbElem += f->sh.texWidth*f->sh.texHeight;
		else          return 0;
	}
	return nbElem;

}*/

/*void Geometry::CopyElemBuffer(BYTE *buffer) {

	int idx = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		int sz = f->sh.texWidth * f->sh.texHeight * sizeof(SHELEM);
		memcpy(buffer + idx, f->mesh, sz);
		idx += sz;
	}

}*/

void Geometry::BuildShapeList() {

	// Shapes used for direction field rendering

	// 3D arrow (direction field)
	int nbDiv = 10;
	double alpha = 2.0*PI / (double)nbDiv;

	arrowList = glGenLists(1);
	glNewList(arrowList, GL_COMPILE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	glBegin(GL_TRIANGLES);

	// Arrow
	for (int i = 0; i < nbDiv; i++) {

		double y1 = sin(alpha*(double)i);
		double z1 = cos(alpha*(double)i);
		double y2 = sin(alpha*(double)((i + 1) % nbDiv));
		double z2 = cos(alpha*(double)((i + 1) % nbDiv));

		glNormal3d(0.0, y1, z1);
		glVertex3d(-0.5, 0.5*y1, 0.5*z1);
		glNormal3d(1.0, 0.0, 0.0);
		glVertex3d(0.5, 0.0, 0.0);
		glNormal3d(0.0, y2, z2);
		glVertex3d(-0.5, 0.5*y2, 0.5*z2);

	}

	// Cap facets
	for (int i = 0; i < nbDiv; i++) {

		double y1 = sin(alpha*(double)i);
		double z1 = cos(alpha*(double)i);
		double y2 = sin(alpha*(double)((i + 1) % nbDiv));
		double z2 = cos(alpha*(double)((i + 1) % nbDiv));

		glNormal3d(-1.0, 0.0, 0.0);
		glVertex3d(-0.5, 0.5*y1, 0.5*z1);
		glNormal3d(-1.0, 0.0, 0.0);
		glVertex3d(-0.5, 0.5*y2, 0.5*z2);
		glNormal3d(-1.0, 0.0, 0.0);
		glVertex3d(-0.5, 0.0, 0.0);

	}

	glEnd();
	glEndList();

	// Shpere list (isotropic case)
	int nbPhi = 16;
	int nbTetha = 7;
	double dphi = 2.0*PI / (double)(nbPhi);
	double dtetha = PI / (double)(nbTetha + 1);

	sphereList = glGenLists(1);
	glNewList(sphereList, GL_COMPILE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glDisable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	glBegin(GL_TRIANGLES);

	for (int i = 0; i <= nbTetha; i++) {
		for (int j = 0; j < nbPhi; j++) {

			VERTEX3D v1, v2, v3, v4;

			v1.x = sin(dtetha*(double)i)*cos(dphi*(double)j);
			v1.y = sin(dtetha*(double)i)*sin(dphi*(double)j);
			v1.z = cos(dtetha*(double)i);

			v2.x = sin(dtetha*(double)(i + 1))*cos(dphi*(double)j);
			v2.y = sin(dtetha*(double)(i + 1))*sin(dphi*(double)j);
			v2.z = cos(dtetha*(double)(i + 1));

			v3.x = sin(dtetha*(double)(i + 1))*cos(dphi*(double)(j + 1));
			v3.y = sin(dtetha*(double)(i + 1))*sin(dphi*(double)(j + 1));
			v3.z = cos(dtetha*(double)(i + 1));

			v4.x = sin(dtetha*(double)i)*cos(dphi*(double)(j + 1));
			v4.y = sin(dtetha*(double)i)*sin(dphi*(double)(j + 1));
			v4.z = cos(dtetha*(double)i);

			if (i < nbTetha) {
				glNormal3d(v1.x, v1.y, v1.z);
				glVertex3d(v1.x, v1.y, v1.z);
				glNormal3d(v2.x, v2.y, v2.z);
				glVertex3d(v2.x, v2.y, v2.z);
				glNormal3d(v3.x, v3.y, v3.z);
				glVertex3d(v3.x, v3.y, v3.z);
			}

			if (i > 0) {
				glNormal3d(v1.x, v1.y, v1.z);
				glVertex3d(v1.x, v1.y, v1.z);
				glNormal3d(v3.x, v3.y, v3.z);
				glVertex3d(v3.x, v3.y, v3.z);
				glNormal3d(v4.x, v4.y, v4.z);
				glVertex3d(v4.x, v4.y, v4.z);
			}

		}
	}

	glEnd();
	glEndList();


}

void Geometry::BuildSelectList() {

	nbSelected = 0;

	selectList = glGenLists(1);
	glNewList(selectList, GL_COMPILE);
	/*
	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	if (antiAliasing){
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	//glBlendFunc(GL_ONE,GL_ZERO);
	}
	glLineWidth(2.0f);


	for(int i=0;i<sh.nbFacet;i++ ) {
	Facet *f = facets[i];
	if( f->selected ) {
	//DrawFacet(f,FALSE);
	DrawFacet(f,1,1,1);
	nbSelected++;
	}
	}
	glLineWidth(1.0f);
	if (antiAliasing) {
	glDisable(GL_BLEND);
	glDisable(GL_LINE_SMOOTH);
	}*/
	glEndList();



	// Second list for usage with POLYGON_OFFSET
	selectList2 = glGenLists(1);
	glNewList(selectList2, GL_COMPILE);
	/*
	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	if (antiAliasing){
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	}
	glLineWidth(2.0f);

	for(int i=0;i<sh.nbFacet;i++ ) {
	Facet *f = facets[i];
	if( f->selected )
	{
	//DrawFacet(f,TRUE,FALSE,TRUE);
	DrawFacet(f,1,1,1);
	}
	}
	glLineWidth(1.0f);
	if (antiAliasing) {
	glDisable(GL_BLEND);
	glDisable(GL_LINE_SMOOTH);
	}*/
	glEndList();


	// Third list with hidden (hole join) edge visible
	selectList3 = glGenLists(1);
	glNewList(selectList3, GL_COMPILE);

	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	if (mApp->antiAliasing){
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	}
	glLineWidth(2.0f);

	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		if (f->selected) {
			//DrawFacet(f,FALSE,TRUE,TRUE);
			DrawFacet(f, 1, 1, 1);
			nbSelected++;
		}
	}
	glLineWidth(1.0f);
	if (mApp->antiAliasing) {
		glDisable(GL_BLEND);
		glDisable(GL_LINE_SMOOTH);
	}
	glEndList();


}

void Geometry::UpdateSelection() {

	DeleteGLLists();
	BuildSelectList();

}

void Geometry::BuildGLList() {

	// Compile geometry for OpenGL
	for (int j = 0; j < sh.nbSuper; j++) {
		lineList[j] = glGenLists(1);
		glNewList(lineList[j], GL_COMPILE);
		for (int i = 0; i < sh.nbFacet; i++) {
			if (facets[i]->sh.superIdx == j)
				DrawFacet(facets[i]);
		}
		glEndList();
	}

	polyList = glGenLists(1);
	glNewList(polyList, GL_COMPILE);
	DrawPolys();
	glEndList();

	BuildSelectList();

}

void Geometry::Rebuild() {

	// Rebuild internal structure on geometry change

	// Remove texture (improvement TODO)
	for (int i = 0; i < sh.nbFacet; i++)
		facets[i]->SetTexture(0.0, 0.0, FALSE);
	// Delete old resources
	DeleteGLLists(TRUE, TRUE);

	// Reinitialise geom
	InitializeGeometry();

}

int Geometry::InvalidateDeviceObjects() {

	DeleteGLLists(TRUE, TRUE);
	DELETE_LIST(arrowList);
	DELETE_LIST(sphereList);
	for (int i = 0; i < sh.nbFacet; i++)
		facets[i]->InvalidateDeviceObjects();

	return GL_OK;

}

void Geometry::ClearFacetMeshLists()
{
	GLProgress *prg = new GLProgress("Please wait...", "Clearing facet meshes...");
	prg->SetVisible(TRUE);
	int nbFacet = mApp->worker.GetGeometry()->GetNbFacet();
	for (int i = 0;i < nbFacet;i++) {
		prg->SetProgress((double)i / (double)nbFacet);
		DELETE_LIST(mApp->worker.GetGeometry()->GetFacet(i)->glElem);
	}
	prg->SetVisible(FALSE);
	SAFE_DELETE(prg);
}

void Geometry::BuildFacetMeshLists()
{
	GLProgress *prg = new GLProgress("Please wait...", "Building facet meshes...");
	prg->SetVisible(TRUE);
	int nbFacet = mApp->worker.GetGeometry()->GetNbFacet();
	for (int i = 0;i < nbFacet;i++) {
		prg->SetProgress((double)i / (double)nbFacet);
		mApp->worker.GetGeometry()->GetFacet(i)->BuildMeshList();
	}
	prg->SetVisible(FALSE);
	SAFE_DELETE(prg);
}

int Geometry::RestoreDeviceObjects() {

	if (!IsLoaded()) return GL_OK;

	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		f->RestoreDeviceObjects();
		BuildFacetList(f);
	}

	BuildGLList();

	return GL_OK;

}

void Geometry::BuildFacetList(Facet *f) {

	// Rebuild OpenGL geomtetry with texture
	//if( f->sh.isTextured && f->textureVisible ) {
	if (f->sh.isTextured) {
		// Facet geometry
		glNewList(f->glList, GL_COMPILE);
		if (f->sh.nbIndex == 3) {
			glBegin(GL_TRIANGLES);
			FillFacet(f, TRUE);
			glEnd();
		}
		else if (f->sh.nbIndex == 4) {
			glBegin(GL_QUADS);
			FillFacet(f, TRUE);
			glEnd();
		}
		else {
			glBegin(GL_TRIANGLES);
			Triangulate(f, TRUE);
			glEnd();
		}
		glEndList();

	}

}

void Geometry::SetFacetTexture(int facet, double ratio, BOOL mesh) {

	Facet *f = facets[facet];
	double nU = Norme(f->sh.U);
	double nV = Norme(f->sh.V);

	if (!f->SetTexture(nU*ratio, nV*ratio, mesh)) {
		char errMsg[512];
		sprintf(errMsg, "Not enough memory to build mesh on Facet %d. ", facet + 1);
		throw Error(errMsg);
	}
	f->tRatio = ratio;
	BuildFacetList(f);

}

void  Geometry::BuildPipe(double L, double R, double s, int step) {
	Clear();

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	sprintf(sh.name, "PIPE%g", L / R);

	int nbDecade = 0;
	int nbTF = 9 * nbDecade;
	int nbTV = 4 * nbTF;

	sh.nbVertex = 2 * step + nbTV;
	if (!(vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D))))
		throw Error("Couldn't allocate memory for vertices");
	memset(vertices3, 0, sh.nbVertex * sizeof(VERTEX3D));

	sh.nbFacet = step + 2 + nbTF;
	sh.nbSuper = 1;
	if (!(facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *))))
		throw Error("Couldn't allocate memory for facets");
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));

	// Vertices
	for (int i = 0; i < step; i++) {
		double angle = (double)i / (double)step * 2 * PI;
		vertices3[2 * i + nbTV].x = R*cos(angle);
		vertices3[2 * i + nbTV].y = R*sin(angle);
		vertices3[2 * i + nbTV].z = 0.0;
		vertices3[2 * i + 1 + nbTV].x = R*cos(angle);
		vertices3[2 * i + 1 + nbTV].y = R*sin(angle);
		vertices3[2 * i + 1 + nbTV].z = L;
	}

	try{
		// Cap facet
		facets[0 + nbTF] = new Facet(step);
		facets[0 + nbTF]->sh.sticking = 1.0;
		for (int i = 0; i < step; i++)
			facets[0 + nbTF]->indices[i] = 2 * i + nbTV;

		facets[1 + nbTF] = new Facet(step);
		facets[1 + nbTF]->sh.sticking = 1.0;
		for (int i = 0; i < step; i++)
			facets[1 + nbTF]->indices[step - i - 1] = 2 * i + 1 + nbTV;

		// Wall facet
		for (int i = 0; i < step; i++) {
			facets[i + 2 + nbTF] = new Facet(4);
			facets[i + 2 + nbTF]->sh.reflectType = REF_MIRROR;
			facets[i + 2 + nbTF]->sh.sticking = s;
			facets[i + 2 + nbTF]->indices[0] = 2 * i + nbTV;
			facets[i + 2 + nbTF]->indices[1] = 2 * i + 1 + nbTV;
			if (i < step - 1) {
				facets[i + 2 + nbTF]->indices[2] = 2 * (i + 1) + 1 + nbTV;
				facets[i + 2 + nbTF]->indices[3] = 2 * (i + 1) + nbTV;
			}
			else {
				facets[i + 2 + nbTF]->indices[2] = 1 + nbTV;
				facets[i + 2 + nbTF]->indices[3] = 0 + nbTV;
			}
		}

		// Volatile facet
		for (int d = 0; d < nbDecade; d++) {
			for (int i = 0; i < 9; i++) {

				double z = (double)(i + 1) * pow(10, (double)d);
				int idx = d * 36 + i * 4;

				vertices3[idx + 0].x = -R;
				vertices3[idx + 0].y = R;
				vertices3[idx + 0].z = z;
				vertices3[idx + 1].x = R;
				vertices3[idx + 1].y = R;
				vertices3[idx + 1].z = z;
				vertices3[idx + 2].x = R;
				vertices3[idx + 2].y = -R;
				vertices3[idx + 2].z = z;
				vertices3[idx + 3].x = -R;
				vertices3[idx + 3].y = -R;
				vertices3[idx + 3].z = z;

				facets[9 * d + i] = new Facet(4);
				facets[9 * d + i]->sh.sticking = 0.0;
				facets[9 * d + i]->sh.opacity = 0.0;
				facets[9 * d + i]->sh.isVolatile = TRUE;
				facets[9 * d + i]->indices[0] = idx + 0;
				facets[9 * d + i]->indices[1] = idx + 1;
				facets[9 * d + i]->indices[2] = idx + 2;
				facets[9 * d + i]->indices[3] = idx + 3;

			}
		}
	}
	catch (std::bad_alloc) {
		Clear();
		throw Error("Couldn't reserve memory for the facets");
	}
	catch (...) {
		throw Error("Unspecified Error while building pipe");
	}
	InitializeGeometry();
	//isLoaded = TRUE;
	strName[0] = _strdup("Pipe");
	strFileName[0] = _strdup("pipe.txt");

}

void Geometry::UpdateName(FileReader *file) {
	UpdateName(file->GetName());
}

void Geometry::UpdateName(const char *fileName) {
	strcpy(sh.name, FileUtils::GetFilename(fileName).c_str());
}

void Geometry::AdjustProfile() {

	// Backward compatibily with TXT profile (To be improved)
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		if (f->sh.profileType == REC_PRESSUREU) {
			VERTEX3D v0;
			Sub(&v0, vertices3 + (f->indices[1]), vertices3 + (f->indices[0])); // v0 = P0P1
			double n0 = Norme(v0);
			double nU = Norme(f->sh.U);
			if (IS_ZERO(n0 - nU)) f->sh.profileType = REC_PRESSUREU; // Select U
			else               f->sh.profileType = REC_PRESSUREV; // Select V
		}
	}

}

void Geometry::LoadASE(FileReader *file, GLProgress *prg) {

	Clear();

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	ASELoader ase(file);
	ase.Load();

	// Compute total of facet
	sh.nbFacet = 0;
	for (int i = 0; i < ase.nbObj; i++) sh.nbFacet += ase.OBJ[i].nb_face;

	// Allocate mem
	sh.nbVertex = 3 * sh.nbFacet;
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memset(vertices3, 0, sh.nbVertex * sizeof(VERTEX3D));

	// Fill 
	int nb = 0;
	for (int i = 0; i < ase.nbObj; i++) {

		for (int j = 0; j < ase.OBJ[i].nb_face; j++) {
			vertices3[3 * nb + 0] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v1];
			vertices3[3 * nb + 1] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v2];
			vertices3[3 * nb + 2] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v3];
			facets[nb] = new Facet(3);
			facets[nb]->indices[0] = 3 * nb + 0;
			facets[nb]->indices[1] = 3 * nb + 1;
			facets[nb]->indices[2] = 3 * nb + 2;
			nb++;
		}

	}

	UpdateName(file);
	sh.nbSuper = 1;
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	isLoaded = TRUE;

}

void Geometry::LoadSTR(FileReader *file, GLProgress *prg) {

	char nPath[512];
	char fPath[512];
	char fName[512];
	char sName[512];
	int nF, nV;
	Facet **F;
	VERTEX3D *V;
	FileReader *fr;

	Clear();

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	// Load multiple structure file
	sh.nbSuper = file->ReadInt();

	strcpy(fPath, file->ReadLine());
	strcpy(nPath, FileUtils::GetPath(file->GetName()).c_str());

	for (int n = 0; n < sh.nbSuper; n++) {

		int i1 = file->ReadInt();
		int i2 = file->ReadInt();
		fr = NULL;
		strcpy(sName, file->ReadWord());
		strName[n] = _strdup(sName);
		char *e = strrchr(strName[n], '.');
		if (e) *e = 0;

		sprintf(fName, "%s%s", nPath, sName);
		if (FileUtils::Exist(fName)) {
			fr = new FileReader(fName);
			strcpy(strPath, nPath);
		}
		else {
			sprintf(fName, "%s\\%s", fPath, sName);
			if (FileUtils::Exist(fName)) {
				fr = new FileReader(fName);
				strcpy(strPath, fPath);
			}
		}

		if (!fr) {
			char errMsg[512];
			sprintf(errMsg, "Cannot find %s", sName);
			throw Error(errMsg);
		}

		strFileName[n] = _strdup(sName);
		LoadTXTGeom(fr, &nV, &nF, &V, &F, n);
		Merge(nV, nF, V, F);
		SAFE_FREE(V);
		SAFE_FREE(F);
		delete fr;

	}

	UpdateName(file);
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::LoadSTL(FileReader *file, GLProgress *prg, double scaleFactor) {

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	char *w;

	prg->SetMessage("Clearing current geometry...");
	Clear();

	// First pass
	prg->SetMessage("Counting facets in STL file...");
	//file->ReadKeyword("solid");
	file->ReadLine(); // solid name
	w = file->ReadWord();
	while (strcmp(w, "facet") == 0) {
		sh.nbFacet++;
		file->JumpSection("endfacet");
		w = file->ReadWord();
	}
	if (strcmp(w, "endsolid") != 0) throw Error("Unexpected or not supported STL keyword, 'endsolid' required\nMaybe the STL file was saved in binary instead of ASCII format?");

	// Allocate mem
	sh.nbVertex = 3 * sh.nbFacet;
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	if (!facets) throw Error("Out of memory: LoadSTL");
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	if (!vertices3) throw Error("Out of memory: LoadSTL");
	memset(vertices3, 0, sh.nbVertex * sizeof(VERTEX3D));

	// Second pass
	prg->SetMessage("Reading facets...");
	file->SeekStart();
	//file->ReadKeyword("solid");
	file->ReadLine();
	for (int i = 0; i < sh.nbFacet; i++) {

		double p = (double)i / (double)(sh.nbFacet);
		prg->SetProgress(p);

		file->ReadKeyword("facet");
		file->ReadKeyword("normal");
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
		file->ReadKeyword("outer");
		file->ReadKeyword("loop");

		file->ReadKeyword("vertex");
		vertices3[3 * i + 0].x = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 0].y = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 0].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		vertices3[3 * i + 1].x = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 1].y = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 1].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		vertices3[3 * i + 2].x = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 2].y = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 2].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("endloop");
		file->ReadKeyword("endfacet");

		try {
			facets[i] = new Facet(3);
		}
		catch (std::bad_alloc& badalloc) {

			throw Error("Out of memory");

		}
		facets[i]->indices[0] = 3 * i + 0;
		facets[i]->indices[1] = 3 * i + 2;
		facets[i]->indices[2] = 3 * i + 1;

	}

	sh.nbSuper = 1;
	UpdateName(file);
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	prg->SetMessage("Initializing geometry...");
	InitializeGeometry();
	isLoaded = TRUE;

}

void Geometry::LoadTXT(FileReader *file, GLProgress *prg) {

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	Clear();
	LoadTXTGeom(file, &(sh.nbVertex), &(sh.nbFacet), &vertices3, &facets);
	UpdateName(file);
	sh.nbSuper = 1;
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;
}

void Geometry::InsertTXT(FileReader *file, GLProgress *prg, BOOL newStr) {

	//Clear();
	int structId = viewStruct;
	if (structId == -1) structId = 0;
	InsertTXTGeom(file, &(sh.nbVertex), &(sh.nbFacet), &vertices3, &facets, structId, newStr);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::InsertSTL(FileReader *file, GLProgress *prg, double scaleFactor, BOOL newStr) {

	//Clear();
	int structId = viewStruct;
	if (structId == -1) structId = 0;
	InsertSTLGeom(file, &(sh.nbVertex), &(sh.nbFacet), &vertices3, &facets, structId, newStr, scaleFactor);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::InsertGEO(FileReader *file, GLProgress *prg, BOOL newStr) {

	//Clear();
	int structId = viewStruct;
	if (structId == -1) structId = 0;
	InsertGEOGeom(file, &(sh.nbVertex), &(sh.nbFacet), &vertices3, &facets, structId, newStr);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;

}

std::vector<std::string> Geometry::InsertSYN(FileReader *file, GLProgress *prg, BOOL newStr) {

	int structId = viewStruct;
	if (structId == -1) structId = 0;
	std::vector<std::string> result = InsertSYNGeom(file, &(sh.nbVertex), &(sh.nbFacet), &vertices3, &facets, structId, newStr);
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	//AdjustProfile();
	return result;
}

void Geometry::LoadTXTGeom(FileReader *file, int *nbV, int *nbF, VERTEX3D **V, Facet ***F, int strIdx) {

	file->ReadInt(); // Unused
	tNbHit = file->ReadLLong();
	tNbLeak = file->ReadLLong();
	tNbDesorption = file->ReadLLong();
	tNbDesorptionMax = file->ReadLLong();

	int nV = file->ReadInt();
	int nF = file->ReadInt();

	// Allocate memory
	Facet   **f = (Facet **)malloc(nF * sizeof(Facet *));
	memset(f, 0, nF * sizeof(Facet *));
	VERTEX3D *v = (VERTEX3D *)malloc(nV * sizeof(VERTEX3D));
	memset(v, 0, nV * sizeof(VERTEX3D)); //avoid selected flag

	// Read geometry vertices
	for (int i = 0; i < nV; i++) {
		v[i].x = file->ReadDouble();
		v[i].y = file->ReadDouble();
		v[i].z = file->ReadDouble();
	}

	// Read geometry facets (indexed from 1)
	for (int i = 0; i < nF; i++) {
		int nb = file->ReadInt();
		f[i] = new Facet(nb);
		for (int j = 0; j < nb; j++)
			f[i]->indices[j] = file->ReadInt() - 1;
	}

	// Read facets params
	for (int i = 0; i<nF; i++) {
		f[i]->LoadTXT(file);
		while ((f[i]->sh.superDest) > sh.nbSuper) { //If facet refers to a structure that doesn't exist, create it
			AddStruct("TXT linked");
		}
		f[i]->sh.superIdx = strIdx;
	}

	SAFE_FREE(*V);
	SAFE_FREE(*F);

	*nbV = nV;
	*nbF = nF;
	*V = v;
	*F = f;

}

void Geometry::InsertTXTGeom(FileReader *file, int *nbVertex, int *nbFacet, VERTEX3D **vertices3, Facet ***facets, int strIdx, BOOL newStruct) {

	UnSelectAll();

	//tNbHit = file->ReadLLong();
	//tNbLeak = file->ReadInt();
	//tNbDesorption = file->ReadLLong();
	//tNbDesorptionMax = file->ReadLLong(); 
	for (int i = 0; i < 5; i++) file->ReadInt(); //leading lines

	int nbNewVertex = file->ReadInt();
	int nbNewFacets = file->ReadInt();

	// Allocate memory
	*facets = (Facet **)realloc(*facets, (nbNewFacets + *nbFacet) * sizeof(Facet **));
	memset(*facets + *nbFacet, 0, nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex + *nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3, *vertices3, (*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3 + *nbVertex, 0, nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3 = tmp_vertices3;

	// Read geometry vertices
	for (int i = *nbVertex; i < (*nbVertex + nbNewVertex); i++) {
		(*vertices3 + i)->x = file->ReadDouble();
		(*vertices3 + i)->y = file->ReadDouble();
		(*vertices3 + i)->z = file->ReadDouble();
		(*vertices3 + i)->selected = FALSE;
	}

	// Read geometry facets (indexed from 1)
	for (int i = *nbFacet; i < (*nbFacet + nbNewFacets); i++) {
		int nb = file->ReadInt();
		*(*facets + i) = new Facet(nb);
		(*facets)[i]->selected = TRUE;
		for (int j = 0; j < nb; j++)
			(*facets)[i]->indices[j] = file->ReadInt() - 1 + *nbVertex;
	}

	// Read facets params
	for (int i = *nbFacet; i<(*nbFacet + nbNewFacets); i++) {
		(*facets)[i]->LoadTXT(file);
		while (((*facets)[i]->sh.superDest) > sh.nbSuper) { //If facet refers to a structure that doesn't exist, create it
			AddStruct("TXT linked");
		}
		if (newStruct) {
			(*facets)[i]->sh.superIdx = sh.nbSuper;
		}
		else {

			(*facets)[i]->sh.superIdx = strIdx;
		}
	}

	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) AddStruct("Inserted TXT file");

}

void Geometry::InsertGEOGeom(FileReader *file, int *nbVertex, int *nbFacet, VERTEX3D **vertices3, Facet ***facets, int strIdx, BOOL newStruct) {

	UnSelectAll();
	//char tmp[512];

	file->ReadKeyword("version"); file->ReadKeyword(":");
	int version2;
	version2 = file->ReadInt();
	if (version2 > GEOVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported GEO version V%d", version2);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalDes"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalLeak"); file->ReadKeyword(":");
	file->ReadLLong();
	if (version2 >= 12) {
		file->ReadKeyword("totalAbs"); file->ReadKeyword(":");
		file->ReadLLong();
		if (version2 >= 15) {
			file->ReadKeyword("totalDist_total");
		}
		else { //between versions 12 and 15
			file->ReadKeyword("totalDist");
		}
		file->ReadKeyword(":");
		file->ReadDouble();
		if (version2 >= 15) {
			file->ReadKeyword("totalDist_fullHitsOnly"); file->ReadKeyword(":");
			file->ReadDouble();
		}
	}
	file->ReadKeyword("maxDes"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("nbVertex"); file->ReadKeyword(":");
	int nbNewVertex = file->ReadInt();
	file->ReadKeyword("nbFacet"); file->ReadKeyword(":");
	int nbNewFacets = file->ReadInt();
	file->ReadKeyword("nbSuper"); file->ReadKeyword(":");
	int nbNewSuper = file->ReadInt();
	int nbF = 0;
	int nbV = 0;
	if (version2 >= 2) {
		file->ReadKeyword("nbFormula"); file->ReadKeyword(":");
		nbF = file->ReadInt();
		file->ReadKeyword("nbView"); file->ReadKeyword(":");
		nbV = file->ReadInt();
	}
	int nbS = 0;
	if (version2 >= 8) {
		file->ReadKeyword("nbSelection"); file->ReadKeyword(":");
		nbS = file->ReadInt();
	}
	if (version2 >= 7) {
		file->ReadKeyword("gasMass"); file->ReadKeyword(":");
		file->ReadDouble();
	}
	if (version2 >= 10) { //time-dependent version
		file->ReadKeyword("userMoments"); file->ReadKeyword("{");
		file->ReadKeyword("nb"); file->ReadKeyword(":");
		int nb = file->ReadInt();

		for (int i = 0; i < nb; i++)
			file->ReadString();
		file->ReadKeyword("}");
	}
	if (version2 >= 11) { //gas pulse parameters
		file->ReadKeyword("desorptionStart"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("desorptionStop"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("timeWindow"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("useMaxwellian"); file->ReadKeyword(":");
		file->ReadInt();
	}

	if (version2 >= 12) { //2013.aug.22
		file->ReadKeyword("calcConstantFlow"); file->ReadKeyword(":");
		file->ReadInt();
	}
	if (version2 >= 2) {
		file->ReadKeyword("formulas"); file->ReadKeyword("{");
		for (int i = 0; i < nbF; i++) {
			char tmpName[256];
			char tmpExpr[512];
			strcpy(tmpName, file->ReadString());
			strcpy(tmpExpr, file->ReadString());
			//mApp->OffsetFormula(tmpExpr, sh.nbFacet-1);
			//mApp->AddFormula(tmpName, tmpExpr);
		}
		file->ReadKeyword("}");

		file->ReadKeyword("views"); file->ReadKeyword("{");
		for (int i = 0; i < nbV; i++) {
			char tmpName[256];
			AVIEW v;
			strcpy(tmpName, file->ReadString());
			v.projMode = file->ReadInt();
			v.camAngleOx = file->ReadDouble();
			v.camAngleOy = file->ReadDouble();
			v.camDist = file->ReadDouble();
			v.camOffset.x = file->ReadDouble();
			v.camOffset.y = file->ReadDouble();
			v.camOffset.z = file->ReadDouble();
			v.performXY = file->ReadInt();

			v.vLeft = file->ReadDouble();
			v.vRight = file->ReadDouble();
			v.vTop = file->ReadDouble();
			v.vBottom = file->ReadDouble();
			mApp->AddView(tmpName, v);
		}
		file->ReadKeyword("}");
	}

	if (version2 >= 8) {
		file->ReadKeyword("selections"); file->ReadKeyword("{");
		for (int i = 0; i < nbS; i++) {
			ASELECTION s;
			char tmpName[256];
			strcpy(tmpName, file->ReadString());
			s.name = _strdup(tmpName);
			s.nbSel = file->ReadInt();
			s.selection = (int *)malloc((s.nbSel)*sizeof(int));

			for (int j = 0; j < s.nbSel; j++) {
				s.selection[j] = file->ReadInt() + sh.nbFacet; //offset facet number by current number of facets
			}
			mApp->AddSelection(s.name, s);
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("structures"); file->ReadKeyword("{");
	for (int i = 0; i < nbNewSuper; i++) {
		strName[sh.nbSuper + i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		/* //Commented out for GEO
		sprintf(tmp,"%s.txt",strName[i]);
		strFileName[i] = _strdup(tmp);
		*/
	}
	file->ReadKeyword("}");

	// Reallocate memory
	*facets = (Facet **)realloc(*facets, (nbNewFacets + *nbFacet) * sizeof(Facet **));
	memset(*facets + *nbFacet, 0, nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex + *nbVertex) * sizeof(VERTEX3D));
	if (!tmp_vertices3) throw Error("Out of memory: InsertGEOGeom");
	memmove(tmp_vertices3, *vertices3, (*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3 + *nbVertex, 0, nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3 = tmp_vertices3;

	// Read geometry vertices
	file->ReadKeyword("vertices"); file->ReadKeyword("{");
	for (int i = *nbVertex; i < (*nbVertex + nbNewVertex); i++) {
		// Check idx
		int idx = file->ReadInt();
		if (idx != i - *nbVertex + 1) throw Error(file->MakeError("Wrong vertex index !"));
		(*vertices3 + i)->x = file->ReadDouble();
		(*vertices3 + i)->y = file->ReadDouble();
		(*vertices3 + i)->z = file->ReadDouble();
		(*vertices3 + i)->selected = FALSE;
	}
	file->ReadKeyword("}");

	if (version2 >= 6) {
		// Read leaks
		file->ReadKeyword("leaks"); file->ReadKeyword("{");
		file->ReadKeyword("nbLeak"); file->ReadKeyword(":");
		int nbleak2 = file->ReadInt();
		for (int i = 0; i < nbleak2; i++) {
			int idx = file->ReadInt();
			//if( idx != i ) throw Error(file->MakeError("Wrong leak index !"));
			file->ReadDouble();
			file->ReadDouble();
			file->ReadDouble();

			file->ReadDouble();
			file->ReadDouble();
			file->ReadDouble();
		}
		file->ReadKeyword("}");

		// Read hit cache
		file->ReadKeyword("hits"); file->ReadKeyword("{");
		file->ReadKeyword("nbHHit"); file->ReadKeyword(":");
		int nbHHit2 = file->ReadInt();
		for (int i = 0; i < nbHHit2; i++) {
			int idx = file->ReadInt();
			//if( idx != i ) throw Error(file->MakeError("Wrong hit cache index !"));
			file->ReadDouble();
			file->ReadDouble();
			file->ReadDouble();

			file->ReadInt();
		}
		file->ReadKeyword("}");
	}

	// Read geometry facets (indexed from 1)
	for (int i = *nbFacet; i < (*nbFacet + nbNewFacets); i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1 - *nbFacet) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nb = file->ReadInt();

		if (nb < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %d has only %d vertices. ", i, nb);
			throw Error(errMsg);
		}

		*(*facets + i) = new Facet(nb);
		(*facets)[i]->LoadGEO(file, version2, nbNewVertex);
		(*facets)[i]->selected = TRUE;
		for (int j = 0; j < nb; j++)
			(*facets)[i]->indices[j] += *nbVertex;
		file->ReadKeyword("}");
		if (newStruct) {
			(*facets)[i]->sh.superIdx += sh.nbSuper;
			if ((*facets)[i]->sh.superDest > 0) (*facets)[i]->sh.superDest += sh.nbSuper;
		}
		else {
			(*facets)[i]->sh.superIdx += strIdx;
			if ((*facets)[i]->sh.superDest > 0) (*facets)[i]->sh.superDest += strIdx;
		}
	}


	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) sh.nbSuper += nbNewSuper;
	else if (sh.nbSuper < strIdx + nbNewSuper) sh.nbSuper = strIdx + nbNewSuper;

}

std::vector<std::string> Geometry::InsertSYNGeom(FileReader *file, int *nbVertex, int *nbFacet, VERTEX3D **vertices3, Facet ***facets, int strIdx, BOOL newStruct) {

	std::vector<std::string> parFileList;
	UnSelectAll();
	//char tmp[512];

	file->ReadKeyword("version"); file->ReadKeyword(":");
	int version2;
	version2 = file->ReadInt();
	if (version2 > SYNVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported SYN version V%d", version2);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalDes"); file->ReadKeyword(":");
	file->ReadLLong();
	if (version2 >= 6) {
		file->ReadKeyword("no_scans"); file->ReadKeyword(":");
		loaded_no_scans = file->ReadDouble();
	}
	else loaded_no_scans = 0;
	file->ReadKeyword("totalLeak"); file->ReadKeyword(":");
	file->ReadLLong();
	if (version2 > 2) {
		file->ReadKeyword("totalFlux"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("totalPower"); file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("maxDes"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("nbVertex"); file->ReadKeyword(":");
	int nbNewVertex = file->ReadInt();
	file->ReadKeyword("nbFacet"); file->ReadKeyword(":");
	int nbNewFacets = file->ReadInt();
	file->ReadKeyword("nbSuper"); file->ReadKeyword(":");
	int nbNewSuper = file->ReadInt();
	int nbF = 0;  std::vector<std::vector<string>> loadFormulas;
	int nbV = 0;
	file->ReadKeyword("nbFormula"); file->ReadKeyword(":");
	nbF = file->ReadInt();
	file->ReadKeyword("nbView"); file->ReadKeyword(":");
	nbV = file->ReadInt();
	int nbS = 0;

	file->ReadKeyword("nbSelection"); file->ReadKeyword(":");
	nbS = file->ReadInt();
	if (version2 > 1) {
		file->ReadKeyword("nbRegions"); file->ReadKeyword(":");
		int nbR = file->ReadInt();
		parFileList.reserve(nbR);

		file->ReadKeyword("PARfiles"); file->ReadKeyword("{");
		for (int i = 0; i < nbR; i++) {
			parFileList.push_back(file->ReadString());
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("formulas"); file->ReadKeyword("{");
	for (int i = 0; i < nbF; i++) {
		char tmpName[256];
		char tmpExpr[512];
		strcpy(tmpName, file->ReadString());
		strcpy(tmpExpr, file->ReadString());
		//mApp->AddFormula(tmpName, tmpExpr); //parse after selection groups are loaded
		std::vector<string> newFormula;
		newFormula.push_back(tmpName);
		mApp->OffsetFormula(tmpExpr, sh.nbFacet); //offset formula
		newFormula.push_back(tmpExpr);
		loadFormulas.push_back(newFormula);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("views"); file->ReadKeyword("{");
	for (int i = 0; i < nbV; i++) {
		char tmpName[256];
		AVIEW v;
		strcpy(tmpName, file->ReadString());
		v.projMode = file->ReadInt();
		v.camAngleOx = file->ReadDouble();
		v.camAngleOy = file->ReadDouble();
		v.camDist = file->ReadDouble();
		v.camOffset.x = file->ReadDouble();
		v.camOffset.y = file->ReadDouble();
		v.camOffset.z = file->ReadDouble();
		v.performXY = file->ReadInt();

		v.vLeft = file->ReadDouble();
		v.vRight = file->ReadDouble();
		v.vTop = file->ReadDouble();
		v.vBottom = file->ReadDouble();
		mApp->AddView(tmpName, v);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("selections"); file->ReadKeyword("{");
	for (int i = 0; i < nbS; i++) {
		ASELECTION s;
		char tmpName[256];
		strcpy(tmpName, file->ReadString());
		s.name = _strdup(tmpName);
		s.nbSel = file->ReadInt();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));

		for (int j = 0; j < s.nbSel; j++) {
			s.selection[j] = file->ReadInt() + sh.nbFacet;
		}
		mApp->AddSelection(s.name, s);
	}
	file->ReadKeyword("}");

	for (int i = 0; i < nbF; i++) { //parse formulas now that selection groups are loaded
		mApp->AddFormula(loadFormulas[i][0].c_str(), loadFormulas[i][1].c_str());
	}

	file->ReadKeyword("structures"); file->ReadKeyword("{");
	for (int i = 0; i < nbNewSuper; i++) {
		strName[sh.nbSuper + i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		/*sprintf(tmp,"%s.txt",strName[i]);
		strFileName[i] = _strdup(tmp);*/
	}
	file->ReadKeyword("}");

	// Reallocate memory
	*facets = (Facet **)realloc(*facets, (nbNewFacets + *nbFacet) * sizeof(Facet **));
	memset(*facets + *nbFacet, 0, nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex + *nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3, *vertices3, (*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3 + *nbVertex, 0, nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3 = tmp_vertices3;

	// Read geometry vertices
	file->ReadKeyword("vertices"); file->ReadKeyword("{");
	for (int i = *nbVertex; i < (*nbVertex + nbNewVertex); i++) {
		// Check idx
		int idx = file->ReadInt();
		if (idx != i - *nbVertex + 1) throw Error(file->MakeError("Wrong vertex index !"));
		(*vertices3 + i)->x = file->ReadDouble();
		(*vertices3 + i)->y = file->ReadDouble();
		(*vertices3 + i)->z = file->ReadDouble();
		(*vertices3 + i)->selected = FALSE;
	}
	file->ReadKeyword("}");


	// Read leaks
	file->ReadKeyword("leaks"); file->ReadKeyword("{");
	file->ReadKeyword("nbLeak"); file->ReadKeyword(":");
	int nbleak_local = file->ReadInt();
	for (int i = 0; i < nbleak_local; i++) {
		int idx = file->ReadInt();
		//if( idx != i ) throw Error(file->MakeError("Wrong leak index !"));
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();

		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
	}
	file->ReadKeyword("}");

	// Read hit cache
	file->ReadKeyword("hits"); file->ReadKeyword("{");
	file->ReadKeyword("nbHHit"); file->ReadKeyword(":");
	int nbHHit_local = file->ReadInt();
	for (int i = 0; i < nbHHit_local; i++) {
		int idx = file->ReadInt();
		//if( idx != i ) throw Error(file->MakeError("Wrong hit cache index !"));
		file->ReadDouble(); //x
		file->ReadDouble(); //y
		file->ReadDouble(); //z
		file->ReadDouble(); //dF
		file->ReadDouble(); //dP
		file->ReadInt();    //type
	}
	file->ReadKeyword("}");


	// Read geometry facets (indexed from 1)
	for (int i = *nbFacet; i < (*nbFacet + nbNewFacets); i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1 - *nbFacet) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nb = file->ReadInt();

		if (nb < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %d has only %d vertices. ", i, nb);
			throw Error(errMsg);
		}

		*(*facets + i) = new Facet(nb);
		(*facets)[i]->LoadSYN(file, mApp->worker.materials, version2, nbNewVertex);
		(*facets)[i]->selected = TRUE;
		for (int j = 0; j < nb; j++)
			(*facets)[i]->indices[j] += *nbVertex;
		file->ReadKeyword("}");
		if (newStruct) {
			(*facets)[i]->sh.superIdx += sh.nbSuper;
			if ((*facets)[i]->sh.superDest > 0) (*facets)[i]->sh.superDest += sh.nbSuper;
		}
		else {

			(*facets)[i]->sh.superIdx += strIdx;
			if ((*facets)[i]->sh.superDest > 0) (*facets)[i]->sh.superDest += strIdx;
		}
	}


	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) sh.nbSuper += nbNewSuper;
	else if (sh.nbSuper < strIdx + nbNewSuper) sh.nbSuper = strIdx + nbNewSuper;
	return parFileList;
}

void Geometry::InsertSTLGeom(FileReader *file, int *nbVertex, int *nbFacet, VERTEX3D **vertices3, Facet ***facets, int strIdx, BOOL newStruct, double scaleFactor) {

	UnSelectAll();
	char *w;

	int nbNewFacets = 0;
	// First pass
	file->ReadKeyword("solid");
	file->ReadLine(); // solid name
	w = file->ReadWord();
	while (strcmp(w, "facet") == 0) {
		nbNewFacets++;
		file->JumpSection("endfacet");
		w = file->ReadWord();
	}
	if (strcmp(w, "endsolid") != 0) throw Error("Unexpected or not supported STL keyword, 'endsolid' required");

	// Allocate memory
	int nbNewVertex = 3 * nbNewFacets;
	*facets = (Facet **)realloc(*facets, (nbNewFacets + *nbFacet) * sizeof(Facet **));
	memset(*facets + *nbFacet, 0, nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex + *nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3, *vertices3, (*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3 + *nbVertex, 0, nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3 = tmp_vertices3;

	// Second pass
	file->SeekStart();
	file->ReadKeyword("solid");
	file->ReadLine();
	for (int i = 0; i < nbNewFacets; i++) {

		file->ReadKeyword("facet");
		file->ReadKeyword("normal");
		file->ReadDouble(); //ignoring normal vector, will be calculated from triangle orientation
		file->ReadDouble();
		file->ReadDouble();
		file->ReadKeyword("outer");
		file->ReadKeyword("loop");

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex + 3 * i + 0].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 0].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 0].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex + 3 * i + 1].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 1].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 1].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex + 3 * i + 2].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 2].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 2].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("endloop");
		file->ReadKeyword("endfacet");

		*(*facets + i + *nbFacet) = new Facet(3);
		(*facets)[i + *nbFacet]->selected = TRUE;
		(*facets)[i + *nbFacet]->indices[0] = *nbVertex + 3 * i + 0;
		(*facets)[i + *nbFacet]->indices[1] = *nbVertex + 3 * i + 1;
		(*facets)[i + *nbFacet]->indices[2] = *nbVertex + 3 * i + 2;

		if (newStruct) {
			(*facets)[i + *nbFacet]->sh.superIdx = sh.nbSuper;
		}
		else {
			(*facets)[i + *nbFacet]->sh.superIdx = strIdx;
		}
	}

	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) AddStruct("Inserted STL file");
}

void Geometry::SaveProfileTXT(FileWriter *file, int super, BOOL saveSelected) {
	// Profiles
	for (int j = 0; j < PROFILE_SIZE; j++)
		file->Write("\n");
}

void Geometry::SaveProfileGEO(FileWriter *file, int super, BOOL saveSelected) {
	file->Write("profiles {\n");
	// Profiles
	int nbProfile = 0;
	file->Write(" number: "); file->WriteInt(nbProfile, "\n");
	file->Write(" facets: ");
	file->Write("\n");
	for (int j = 0; j < PROFILE_SIZE; j++)
		file->Write("\n");
	file->Write("}\n");
}

void Geometry::LoadGEO(FileReader *file, GLProgress *prg, LEAK *pleak, int *nbleak, HIT *pHits, int *nbHHit, int *version) {

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	prg->SetMessage("Clearing current geometry...");
	Clear();
	//mApp->ClearFormula();


	// Globals
	char tmp[512];
	prg->SetMessage("Reading GEO file header...");
	file->ReadKeyword("version"); file->ReadKeyword(":");
	*version = file->ReadInt();
	if (*version > GEOVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported GEO version V%d", *version);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit"); file->ReadKeyword(":");
	tNbHit = 0; file->ReadLLong();
	file->ReadKeyword("totalDes"); file->ReadKeyword(":");
	tNbDesorption = 0; file->ReadLLong();
	file->ReadKeyword("totalLeak"); file->ReadKeyword(":");
	tNbLeak = 0; file->ReadLLong();
	if (*version >= 12) {
		file->ReadKeyword("totalAbs"); file->ReadKeyword(":");
		tNbAbsorption = 0; file->ReadLLong();
		if (*version >= 15) {
			file->ReadKeyword("totalDist_total");
		}
		else { //between versions 12 and 15
			file->ReadKeyword("totalDist");
		}
		file->ReadKeyword(":");
		distTraveledTotal = 0.0; file->ReadDouble();
		if (*version >= 15) {
			file->ReadKeyword("totalDist_fullHitsOnly"); file->ReadKeyword(":");
			file->ReadDouble();
		}
	}
	else {
		tNbAbsorption = 0;
		distTraveledTotal = 0.0;
	}
	file->ReadKeyword("maxDes"); file->ReadKeyword(":");
	tNbDesorptionMax = file->ReadLLong();
	file->ReadKeyword("nbVertex"); file->ReadKeyword(":");
	sh.nbVertex = file->ReadInt();
	file->ReadKeyword("nbFacet"); file->ReadKeyword(":");
	sh.nbFacet = file->ReadInt();
	file->ReadKeyword("nbSuper"); file->ReadKeyword(":");
	sh.nbSuper = file->ReadInt();
	int nbF = 0; std::vector<std::vector<string>> loadFormulas;
	int nbV = 0;
	if (*version >= 2) {
		file->ReadKeyword("nbFormula"); file->ReadKeyword(":");
		nbF = file->ReadInt(); loadFormulas.reserve(nbF);
		file->ReadKeyword("nbView"); file->ReadKeyword(":");
		nbV = file->ReadInt();
	}
	int nbS = 0;
	if (*version >= 8) {
		file->ReadKeyword("nbSelection"); file->ReadKeyword(":");
		nbS = file->ReadInt();
	}
	if (*version >= 7) {
		file->ReadKeyword("gasMass"); file->ReadKeyword(":");
		file->ReadDouble(); //gas mass
	}
	if (*version >= 10) { //time-dependent version
		file->ReadKeyword("userMoments"); file->ReadKeyword("{");
		file->ReadKeyword("nb"); file->ReadKeyword(":");
		int nb = file->ReadInt();

		for (int i = 0; i < nb; i++) {
			char tmpExpr[512];
			strcpy(tmpExpr, file->ReadString());
			//mApp->worker.userMoments.push_back(tmpExpr);
			//mApp->worker.AddMoment(mApp->worker.ParseMoment(tmpExpr));
		}
		file->ReadKeyword("}");
	}
	if (*version >= 11) { //pulse version
		file->ReadKeyword("desorptionStart"); file->ReadKeyword(":");
		//worker->desorptionStartTime=
		file->ReadDouble();
		file->ReadKeyword("desorptionStop"); file->ReadKeyword(":");
		//worker->desorptionStopTime=
		file->ReadDouble();
		file->ReadKeyword("timeWindow"); file->ReadKeyword(":");
		//worker->timeWindowSize=
		file->ReadDouble();
		file->ReadKeyword("useMaxwellian"); file->ReadKeyword(":");
		//worker->useMaxwellDistribution=
		file->ReadInt();
	}
	if (*version >= 12) { //2013.aug.22
		file->ReadKeyword("calcConstantFlow"); file->ReadKeyword(":");
		//worker->calcConstantFlow=
		file->ReadInt();
	}
	if (*version >= 2) {
		file->ReadKeyword("formulas"); file->ReadKeyword("{");
		for (int i = 0; i < nbF; i++) {
			char tmpName[256];
			char tmpExpr[512];
			strcpy(tmpName, file->ReadString());
			strcpy(tmpExpr, file->ReadString());
			//mApp->AddFormula(tmpName,tmpExpr);
			std::vector<string> newFormula;
			newFormula.push_back(tmpName);
			newFormula.push_back(tmpExpr);
			loadFormulas.push_back(newFormula);
		}
		file->ReadKeyword("}");

		file->ReadKeyword("views"); file->ReadKeyword("{");
		for (int i = 0; i < nbV; i++) {
			char tmpName[256];
			AVIEW v;
			strcpy(tmpName, file->ReadString());
			v.projMode = file->ReadInt();
			v.camAngleOx = file->ReadDouble();
			v.camAngleOy = file->ReadDouble();
			v.camDist = file->ReadDouble();
			v.camOffset.x = file->ReadDouble();
			v.camOffset.y = file->ReadDouble();
			v.camOffset.z = file->ReadDouble();
			v.performXY = file->ReadInt();

			v.vLeft = file->ReadDouble();
			v.vRight = file->ReadDouble();
			v.vTop = file->ReadDouble();
			v.vBottom = file->ReadDouble();
			mApp->AddView(tmpName, v);
		}
		file->ReadKeyword("}");
	}

	if (*version >= 8) {
		file->ReadKeyword("selections"); file->ReadKeyword("{");
		for (int i = 0; i < nbS; i++) {
			ASELECTION s;
			char tmpName[256];
			strcpy(tmpName, file->ReadString());
			s.name = _strdup(tmpName);
			s.nbSel = file->ReadInt();
			s.selection = (int *)malloc((s.nbSel)*sizeof(int));

			for (int j = 0; j < s.nbSel; j++) {
				s.selection[j] = file->ReadInt();
			}
			mApp->AddSelection(s.name, s);
		}
		file->ReadKeyword("}");
	}

	for (int i = 0; i < nbF; i++) { //parse formulas now that selection groups are loaded
		mApp->AddFormula(loadFormulas[i][0].c_str(), loadFormulas[i][1].c_str());
	}

	file->ReadKeyword("structures"); file->ReadKeyword("{");
	for (int i = 0; i < sh.nbSuper; i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		sprintf(tmp, "%s.txt", strName[i]);
		strFileName[i] = _strdup(tmp);
	}
	file->ReadKeyword("}");

	// Allocate memory
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memset(vertices3, 0, sh.nbVertex * sizeof(VERTEX3D));

	// Read vertices
	prg->SetMessage("Reading vertices...");
	file->ReadKeyword("vertices"); file->ReadKeyword("{");
	for (int i = 0; i < sh.nbVertex; i++) {
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1) throw Error(file->MakeError("Wrong vertex index !"));
		vertices3[i].x = file->ReadDouble();
		vertices3[i].y = file->ReadDouble();
		vertices3[i].z = file->ReadDouble();
		vertices3[i].selected = FALSE;
	}
	file->ReadKeyword("}");

	if (*version >= 6) {
		prg->SetMessage("Reading leaks and hits...");
		// Read leaks
		file->ReadKeyword("leaks"); file->ReadKeyword("{");
		file->ReadKeyword("nbLeak"); file->ReadKeyword(":");
		*nbleak = file->ReadInt();
		for (int i = 0; i < *nbleak; i++) {
			int idx = file->ReadInt();
			if (idx != i) throw Error(file->MakeError("Wrong leak index !"));
			(pleak + i)->pos.x = file->ReadDouble();
			(pleak + i)->pos.y = file->ReadDouble();
			(pleak + i)->pos.z = file->ReadDouble();

			(pleak + i)->dir.x = file->ReadDouble();
			(pleak + i)->dir.y = file->ReadDouble();
			(pleak + i)->dir.z = file->ReadDouble();
		}
		file->ReadKeyword("}");

		// Read hit cache
		file->ReadKeyword("hits"); file->ReadKeyword("{");
		file->ReadKeyword("nbHHit"); file->ReadKeyword(":");
		*nbHHit = file->ReadInt();
		for (int i = 0; i < *nbHHit; i++) {
			int idx = file->ReadInt();
			if (idx != i) throw Error(file->MakeError("Wrong hit cache index !"));
			(pHits + i)->pos.x = file->ReadDouble();
			(pHits + i)->pos.y = file->ReadDouble();
			(pHits + i)->pos.z = file->ReadDouble();

			(pHits + i)->type = file->ReadInt();
		}
		file->ReadKeyword("}");
	}

	// Read facets
	prg->SetMessage("Reading facets...");
	for (int i = 0; i < sh.nbFacet; i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nbI = file->ReadInt();
		if (nbI < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %d has only %d vertices. ", i, nbI);
			throw Error(errMsg);
		}
		prg->SetProgress((float)i / sh.nbFacet);
		facets[i] = new Facet(nbI);
		facets[i]->LoadGEO(file, *version, sh.nbVertex);
		file->ReadKeyword("}");
	}

	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;
	UpdateName(file);

	// Update mesh
	/*prg->SetMessage("Building mesh..."); GEO file textures not loaded in Synrad
	for(int i=0;i<sh.nbFacet;i++) {
	double p = (double)i/(double)sh.nbFacet;
	prg->SetProgress(p);
	Facet *f = facets[i];
	if (!f->SetTexture(f->sh.texWidthD, f->sh.texHeightD, f->hasMesh)) {
	char errMsg[512];
	sprintf(errMsg, "Not enough memory to build mesh on Facet %d. ", i + 1);
	throw Error(errMsg);
	}
	BuildFacetList(f);
	double nU = Norme(&(f->sh.U));
	f->tRatio = f->sh.texWidthD / nU;
	}*/
}

bool Geometry::LoadTextures(FileReader *file, GLProgress *prg, Dataport *dpHit, int version) {

	if (file->SeekFor("{textures}")) {
		char tmp[256];
		//versions 3+
		// Block dpHit during the whole disc reading

		AccessDataport(dpHit);

		// Globals
		BYTE *buffer = (BYTE *)dpHit->buff;
		SHGHITS *gHits = (SHGHITS *)buffer;

		gHits->total.nbHit = tNbHit;
		gHits->total.nbDesorbed = tNbDesorption;
		gHits->total.nbAbsorbed = tNbAbsorption;
		gHits->nbLeakTotal = tNbLeak;
		gHits->total.fluxAbs = tFlux;
		gHits->total.powerAbs = tPower;
		gHits->distTraveledTotal = distTraveledTotal;

		// Read facets
		file->ReadKeyword("minHit_MC"); file->ReadKeyword(":");
		gHits->minHit_MC = file->ReadLLong();
		file->ReadKeyword("maxHit_MC"); file->ReadKeyword(":");
		gHits->maxHit_MC = file->ReadLLong();
		file->ReadKeyword("minHit_flux"); file->ReadKeyword(":");
		gHits->minHit_flux = file->ReadDouble();
		file->ReadKeyword("maxHit_flux"); file->ReadKeyword(":");
		gHits->maxHit_flux = file->ReadDouble();
		file->ReadKeyword("minHit_power"); file->ReadKeyword(":");
		gHits->minHit_power = file->ReadDouble();
		file->ReadKeyword("maxHit_power"); file->ReadKeyword(":");
		gHits->maxHit_power = file->ReadDouble();

		for (int i = 0; i < sh.nbFacet; i++) {
			Facet *f = facets[i];
			if (f->hasMesh) {
				prg->SetProgress((double)i / (double)sh.nbFacet);
				file->ReadKeyword("texture_facet");
				// Check idx
				int idx = file->ReadInt();

				if (idx != i + 1) {
					sprintf(tmp, "Wrong facet index. Expected %d, read %d.", i + 1, idx);
					throw Error(file->MakeError(tmp));
				}

				//Now load values
				file->ReadKeyword("{");

				int ix, iy;

				int profSize = (f->sh.isProfile) ? (PROFILE_SIZE*(sizeof(llong) + 2 * sizeof(double))) : 0;
				int nbE = f->sh.texHeight*f->sh.texWidth;
				llong *hits_MC = (llong *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize));
				double *hits_flux = (double *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + nbE*sizeof(llong)));
				double *hits_power = (double *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + nbE*(sizeof(llong) + sizeof(double))));

				int texWidth_file, texHeight_file;
				//In case of rounding errors, the file might contain different texture dimensions than expected.
				if (version >= 8) {
					file->ReadKeyword("width"); file->ReadKeyword(":"); texWidth_file = file->ReadInt();
					file->ReadKeyword("height"); file->ReadKeyword(":"); texHeight_file = file->ReadInt();
				}
				else {
					texWidth_file = f->sh.texWidth;
					texHeight_file = f->sh.texHeight;
				}

				for (iy = 0; iy < (MIN(f->sh.texHeight, texHeight_file)); iy++) { //MIN: If stored texture is larger, don't read extra cells
					for (ix = 0; ix<(MIN(f->sh.texWidth, texWidth_file)); ix++) { //MIN: If stored texture is larger, don't read extra cells
						int index = iy*(f->sh.texWidth) + ix;
						*(hits_MC + index) = file->ReadLLong();
						if (version >= 7) file->ReadDouble(); //cell area
						*(hits_flux + index) = file->ReadDouble();
						*(hits_power + index) = file->ReadDouble();

						//Normalize by area
						if (f->GetMeshArea(index)>0.0) {
							*(hits_flux + index) /= f->GetMeshArea(index);
							*(hits_power + index) /= f->GetMeshArea(index);
						}
					}
					for (int ie = 0; ie < texWidth_file - f->sh.texWidth; ie++) {//Executed if file texture is bigger than expected texture
						//Read extra cells from file without doing anything
						file->ReadLLong();
						if (version >= 7) file->ReadDouble(); //cell area
						file->ReadDouble();
						file->ReadDouble();
					}
				}
				for (int ie = 0; ie < texHeight_file - f->sh.texHeight; ie++) {//Executed if file texture is bigger than expected texture
					//Read extra cells ffrom file without doing anything
					for (int iw = 0; iw < texWidth_file; iw++) {
						file->ReadLLong();
						if (version >= 7) file->ReadDouble(); //cell area
						file->ReadDouble();
						file->ReadDouble();
					}
				}
				file->ReadKeyword("}");

			}
		}

		ReleaseDataport(dpHit);
		return true;

	}
	else
	{
		//old versions
		return false;
	}

}

// -----------------------------------------------------------
/*
void Geometry::SaveGEO(FileWriter *file,GLProgress *prg,Dataport *dpHit,BOOL saveSelected,LEAK *pleak,int *nbleakSave,HIT *pHits,int *nbHHitSave,BOOL crashSave) {

SynRad *mApp = (SynRad *)theApp;
prg->SetMessage("Counting hits...");
if(!IsLoaded()) throw Error("Nothing to save !");


// Block dpHit during the whole disc writing
AccessDataport(dpHit);

// Globals
//BYTE *buffer = (BYTE *)dpHit->buff;
//SHGHITS *gHits = (SHGHITS *)buffer;


float dCoef = 1.0f;


prg->SetMessage("Writing geometry details...");
file->Write("version:");file->WriteInt(GEOVERSION,"\n");
file->Write("totalHit:");file->WriteLLong(0,"\n");
file->Write("totalDes:");file->WriteLLong(0,"\n");
file->Write("totalLeak:");file->WriteInt(0,"\n");
file->Write("maxDes:");file->WriteLLong(tNbDesorptionMax,"\n");
file->Write("nbVertex:");file->WriteInt(sh.nbVertex,"\n");
file->Write("nbFacet:");file->WriteInt(saveSelected?nbSelected:sh.nbFacet,"\n");
file->Write("nbSuper:");file->WriteInt(sh.nbSuper,"\n");
file->Write("nbFormula:");file->WriteInt(mApp->nbFormula,"\n");
file->Write("nbView:");file->WriteInt(mApp->nbView,"\n");
file->Write("nbSelection:");file->WriteInt(mApp->nbSelection,"\n");
file->Write("gasMass:");file->WriteDouble(gasMass,"\n");

file->Write("formulas {\n");
for(int i=0;i<mApp->nbFormula;i++) {
file->Write("  \"");
file->Write(mApp->formulas[i].parser->GetName());
file->Write("\" \"");
file->Write(mApp->formulas[i].parser->GetExpression());
file->Write("\"\n");
}
file->Write("}\n");

file->Write("views {\n");
for(int i=0;i<mApp->nbView;i++) {
file->Write("  \"");
file->Write(mApp->views[i].name);
file->Write("\"\n");
file->WriteInt(mApp->views[i].projMode," ");
file->WriteDouble(mApp->views[i].camAngleOx," ");
file->WriteDouble(mApp->views[i].camAngleOy," ");
file->WriteDouble(mApp->views[i].camDist," ");
file->WriteDouble(mApp->views[i].camOffset.x," ");
file->WriteDouble(mApp->views[i].camOffset.y," ");
file->WriteDouble(mApp->views[i].camOffset.z," ");
file->WriteInt(mApp->views[i].performXY," ");
file->WriteDouble(mApp->views[i].vLeft," ");
file->WriteDouble(mApp->views[i].vRight," ");
file->WriteDouble(mApp->views[i].vTop," ");
file->WriteDouble(mApp->views[i].vBottom,"\n");
}
file->Write("}\n");

file->Write("selections {\n");
for(int i=0;i<mApp->nbSelection;i++) {
file->Write("  \"");
file->Write(mApp->selections[i].name);
file->Write("\"\n ");
file->WriteInt(mApp->selections[i].nbSel,"\n");
for (int j=0;j<mApp->selections[i].nbSel;j++) {
file->Write("  ");
file->WriteInt(mApp->selections[i].selection[j],"\n");
}
}
file->Write("}\n");

file->Write("structures {\n");
for(int i=0;i<sh.nbSuper;i++) {
file->Write("  \"");
file->Write(strName[i]);
file->Write("\"\n");
}
file->Write("}\n");
//vertices
prg->SetMessage("Writing vertices...");
file->Write("vertices {\n");
for(int i=0;i<sh.nbVertex;i++) {
prg->SetProgress(0.33*((double)i/(double)sh.nbVertex));
file->Write("  ");
file->WriteInt(i+1," ");
file->WriteDouble(vertices3[i].x," ");
file->WriteDouble(vertices3[i].y," ");
file->WriteDouble(vertices3[i].z,"\n");
}
file->Write("}\n");

//leaks
prg->SetMessage("Writing leaks...");
file->Write("leaks {\n");
file->Write("  nbLeak:");file->WriteInt(0,"\n");
file->Write("}\n");

//hit cache (lines and dots)
prg->SetMessage("Writing hits...");
file->Write("hits {\n");
file->Write("  nbHHit:");file->WriteInt(0,"\n");
file->Write("}\n");

//facets

prg->SetMessage("Writing facets...");

for(int i=0,k=0;i<sh.nbFacet;i++) {
prg->SetProgress(0.33+((double)i/(double)sh.nbFacet) *0.33);
if( saveSelected ) {
if( facets[i]->selected ) { facets[i]->SaveGEO(file,k);k++; }
} else {
facets[i]->SaveGEO(file,i);
}
}

prg->SetMessage("Writing profiles...");
SaveProfileGEO(file,-1,saveSelected);

///Save textures, for GEO file version 3
char tmp[256];
sprintf(tmp,"{textures}\n");
file->Write(tmp);
file->Write("minHit:");file->WriteDouble(0,"\n");
file->Write("maxHit:");file->WriteDouble(1,"\n");

ReleaseDataport(dpHit);

}
*/
// -----------------------------------------------------------

void Geometry::SaveTXT(FileWriter *file, Dataport *dpHit, BOOL saveSelected) {

	if (!IsLoaded()) throw Error("Nothing to save !");

	// Unused
	file->WriteInt(0, "\n");

	// Block dpHit during the whole disc writing
	AccessDataport(dpHit);

	// Globals
	//BYTE *buffer = (BYTE *)dpHit->buff;
	//SHGHITS *gHits = (SHGHITS *)buffer;

	// Unused
	file->WriteLLong(0, "\n");
	file->WriteInt(0, "\n");
	file->WriteLLong(0, "\n");
	file->WriteLLong(tNbDesorptionMax, "\n");

	file->WriteInt(sh.nbVertex, "\n");
	file->WriteInt(saveSelected ? nbSelected : sh.nbFacet, "\n");

	// Read geometry vertices
	for (int i = 0; i < sh.nbVertex; i++) {
		file->WriteDouble(vertices3[i].x, " ");
		file->WriteDouble(vertices3[i].y, " ");
		file->WriteDouble(vertices3[i].z, "\n");
	}

	// Facets
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		int j;
		if (saveSelected) {
			if (f->selected) {
				file->WriteInt(f->sh.nbIndex, " ");
				for (j = 0; j < f->sh.nbIndex - 1; j++)
					file->WriteInt(f->indices[j] + 1, " ");
				file->WriteInt(f->indices[j] + 1, "\n");
			}
		}
		else {
			file->WriteInt(f->sh.nbIndex, " ");
			for (j = 0; j < f->sh.nbIndex - 1; j++)
				file->WriteInt(f->indices[j] + 1, " ");
			file->WriteInt(f->indices[j] + 1, "\n");
		}
	}

	// Params
	for (int i = 0; i < sh.nbFacet; i++) {

		// Update facet hits from shared mem
		Facet *f = facets[i];
		//SHHITS *shF = (SHHITS *)(buffer + f->sh.hitOffset);
		//memcpy(&(f->sh.counter),shF,sizeof(SHHITS));
		if (saveSelected) {
			if (f->selected) f->SaveTXT(file);
		}
		else {
			f->SaveTXT(file);
		}

	}

	SaveProfileTXT(file, -1, saveSelected);

	ReleaseDataport(dpHit);

}

void Geometry::ExportTextures(FILE *file, int grouping, int mode, double no_scans, Dataport *dpHit, BOOL saveSelected) {

	//if(!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing
	BYTE *buffer = NULL;
	if (dpHit)
		if (AccessDataport(dpHit))
			buffer = (BYTE *)dpHit->buff;

	// Globals
	//BYTE *buffer = (BYTE *)dpHit->buff;
	//SHGHITS *gHits = (SHGHITS *)buffer;

	if (grouping == 1) fprintf(file, "X_coord_cm\tY_coord_cm\tZ_coord_cm\tValue\t\n"); //mode 10: special ANSYS export

	// Facets
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];

		if (f->selected) {
			if (grouping == 0) fprintf(file, "FACET%d\n", i + 1); //mode 10: special ANSYS export
			//SHELEM *mesh = f->mesh;
			if (f->cellPropertiesIds) {
				char tmp[256];
				char out[256];
				float dCoef = 1.0f;
				if (!buffer) return;
				SHGHITS *shGHit = (SHGHITS *)buffer;
				int w = f->sh.texWidth;
				int h = f->sh.texHeight;
				int nbE = w*h;
				int profSize = (f->sh.isProfile) ? (PROFILE_SIZE*(sizeof(llong) + 2 * sizeof(double))) : 0;

				llong *hits_MC = (llong *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize));
				double *hits_flux = (double *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize + nbE*sizeof(llong)));
				double *hits_power = (double *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize + nbE*(sizeof(llong) + sizeof(double))));
				double norm = 1.0 / no_scans; //normalize values by number of scans (and don't normalize by area...)

				for (int i = 0; i < w; i++) {
					for (int j = 0; j < h; j++) {
						int index = i + j*w;
						tmp[0] = out[0] = 0;
						switch (mode) {

						case 0: // Element area
							sprintf(tmp, "%g", f->GetMeshArea(index));
							break;

						case 1: // MC_hits
							if (!grouping || hits_MC[index]) sprintf(tmp, "%I64d", hits_MC[index]);
							break;

						case 2: // Flux
							if (!grouping || hits_flux[index]) sprintf(tmp, "%g", hits_flux[index] * f->GetMeshArea(i+j*w)*norm);
							break;

						case 3: // Power
							if (!grouping || hits_power[index]) sprintf(tmp, "%g", hits_power[index] * f->GetMeshArea(i+j*w)*norm);
							break;

						case 4: // Flux/area
							if (!grouping || hits_flux[index]) sprintf(tmp, "%g", hits_flux[index] * norm);
							break;

						case 5: // Power/area
							if (!grouping || hits_power[index]) sprintf(tmp, "%g", hits_power[index] * 0.01*norm); //Don't write 0 powers
							break;
						}

						if (grouping == 1 && tmp && tmp[0]) {
							VERTEX2D center = f->GetMeshCenter(index);
							sprintf(out, "%g\t%g\t%g\t%s\t\n",
								f->sh.O.x + center.u*f->sh.U.x + center.v*f->sh.V.x,
								f->sh.O.y + center.u*f->sh.U.y + center.v*f->sh.V.y,
								f->sh.O.z + center.u*f->sh.U.z + center.v*f->sh.V.z,
								tmp);
						}
						else sprintf(out, "%s", tmp);

						if (out) fprintf(file, "%s", out);
						if (j < w - 1 && grouping == 0)
							fprintf(file, "\t");
					}
					if (grouping == 0) fprintf(file, "\n");
				}
			}
			else {
				fprintf(file, "No mesh.\n");
			}
			if (grouping == 0) fprintf(file, "\n"); //Current facet exported. 
		}

	}

	ReleaseDataport(dpHit);

}

void Geometry::SaveDesorption(FILE *file, Dataport *dpHit, BOOL selectedOnly, int mode, double eta0, double alpha, Distribution2D *distr) {

	if (!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing
	BYTE *buffer = NULL;
	if (dpHit)
		if (AccessDataport(dpHit))
			buffer = (BYTE *)dpHit->buff;

	Worker *worker = &(mApp->worker);

	// Globals
	//BYTE *buffer = (BYTE *)dpHit->buff;
	//SHGHITS *gHits = (SHGHITS *)buffer;

	// Facets
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];

		if (f->selected || !selectedOnly) {

			//SHELEM *mesh = f->mesh;
			if (f->cellPropertiesIds) {
				fprintf(file, "facet %d {\n", i + 1);
				char tmp[256];
				fprintf(file, "cell_size_cm: %g\n", 1.0 / f->tRatio);
				float dCoef = 1.0f;
				if (!buffer) return;

				int w = f->sh.texWidth;
				int h = f->sh.texHeight;
				int textureSize_double = w*h*sizeof(double);
				int textureSize_llong = w*h*sizeof(llong);
				SHGHITS *shGHit = (SHGHITS *)buffer;
				int profile_memory = PROFILE_SIZE*(2 * sizeof(double) + sizeof(llong));
				int profSize = (f->sh.isProfile) ? profile_memory : 0;
				double *hits_flux = (double *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize + textureSize_llong));


				for (int i = 0; i < w; i++) {
					for (int j = 0; j < h; j++) {
						double dose = hits_flux[i + j*w] * f->GetMeshArea(i+j*w) / worker->no_scans;
						double val;
						if (mode == 1) { //no conversion
							val = dose;
						}
						else if (mode == 2) { //use equation
							double eta;
							if (dose < VERY_SMALL) eta = 0.0;
							else eta = eta0*pow(dose, alpha);
							val = dose*eta;
						}
						else if (mode == 3) {  //use file
							double eta = distr->InterpolateY(dose);
							val = dose*eta;
						}
						sprintf(tmp, "%g", val);
						if (tmp) fprintf(file, "%s", tmp);
						//fprintf(file,"%s",tmp);
						if (j < h - 1)
							fprintf(file, "\t");
					}
					fprintf(file, "\n");
				}
				fprintf(file, "}\n"); //Current facet exported. 
			}

		}

	}

	ReleaseDataport(dpHit);

}

void Geometry::SaveSTR(Dataport *dpHit, BOOL saveSelected) {

	if (!IsLoaded()) throw Error("Nothing to save !");
	if (sh.nbSuper < 1) throw Error("Cannot save single structure in STR format");

	// Block dpHit during the whole disc writting
	AccessDataport(dpHit);
	for (int i = 0; i < sh.nbSuper; i++)
		SaveSuper(dpHit, i);
	ReleaseDataport(dpHit);

}

void Geometry::SaveSuper(Dataport *dpHit, int s) {

	char fName[512];
	sprintf(fName, "%s/%s", strPath, strFileName[s]);
	FileWriter *file = new FileWriter(fName);

	// Unused
	file->WriteInt(0, "\n");

	// Globals
	BYTE *buffer = (BYTE *)dpHit->buff;
	SHGHITS *gHits = (SHGHITS *)buffer;

	//Extract data of the specified super structure
	llong totHit = 0;
	llong totAbs = 0;
	llong totDes = 0;
	int *refIdx = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(refIdx, 0xFF, sh.nbVertex*sizeof(int));
	int nbV = 0;
	int nbF = 0;

	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		if (f->sh.superIdx == s) {
			totHit += f->sh.counter.nbHit;
			totAbs += f->sh.counter.nbAbsorbed;
			totDes += f->sh.counter.nbDesorbed;
			for (int j = 0; j < f->sh.nbIndex; j++)
				refIdx[f->indices[j]] = 1;
			nbF++;
		}
	}

	for (int i = 0; i < sh.nbVertex; i++) {
		if (refIdx[i] >= 0) {
			refIdx[i] = nbV;
			nbV++;
		}
	}

	file->WriteLLong(0, "\n");
	file->WriteInt(0, "\n");
	file->WriteLLong(0, "\n");
	file->WriteLLong(0, "\n");

	file->WriteInt(nbV, "\n");
	file->WriteInt(nbF, "\n");

	// Read geometry vertices
	for (int i = 0; i < sh.nbVertex; i++) {
		if (refIdx[i] >= 0) {
			file->WriteDouble(vertices3[i].x, " ");
			file->WriteDouble(vertices3[i].y, " ");
			file->WriteDouble(vertices3[i].z, "\n");
		}
	}

	// Facets
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		int j;
		if (f->sh.superIdx == s) {
			file->WriteInt(f->sh.nbIndex, " ");
			for (j = 0; j < f->sh.nbIndex - 1; j++)
				file->WriteInt(refIdx[f->indices[j]] + 1, " ");
			file->WriteInt(refIdx[f->indices[j]] + 1, "\n");
		}
	}

	// Params
	for (int i = 0; i < sh.nbFacet; i++) {

		// Update facet hits from shared mem
		Facet *f = facets[i];
		if (f->sh.superIdx == s) {
			SHHITS *shF = (SHHITS *)(buffer + f->sh.hitOffset);
			memcpy(&(f->sh.counter), shF, sizeof(SHHITS));
			f->SaveTXT(file);
		}

	}

	SaveProfileTXT(file, s);

	SAFE_DELETE(file);
	free(refIdx);

}



BOOL Geometry::IsLoaded() {
	return isLoaded;
}

void Geometry::SaveSYN(FileWriter *file, GLProgress *prg, Dataport *dpHit, BOOL saveSelected, LEAK *pleak,
	int *nbleakSave, HIT *pHits, int *nbHHitSave, BOOL crashSave) {

	prg->SetMessage("Counting hits...");
	if (!IsLoaded()) throw Error("Nothing to save !");


	// Block dpHit during the whole disc writing
	if (!crashSave && !saveSelected) AccessDataport(dpHit);

	// Globals
	BYTE *buffer;
	if (!crashSave && !saveSelected) buffer = (BYTE *)dpHit->buff;
	SHGHITS *gHits;
	if (!crashSave && !saveSelected) gHits = (SHGHITS *)buffer;


	double dCoef = 1.0;
	int ix, iy;

	Worker *worker = &(mApp->worker);
	texCMin_MC = (!crashSave && !saveSelected) ? gHits->minHit_MC : 0;
	texCMax_MC = (!crashSave && !saveSelected) ? gHits->maxHit_MC : 1;
	texCMin_flux = (!crashSave && !saveSelected) ? gHits->minHit_flux : 0;
	texCMax_flux = (!crashSave && !saveSelected) ? gHits->maxHit_flux : 1;
	texCMin_power = (!crashSave && !saveSelected) ? gHits->minHit_power : 0;
	texCMax_power = (!crashSave && !saveSelected) ? gHits->maxHit_power : 1;

	prg->SetMessage("Writing geometry details...");
	file->Write("version:"); file->WriteInt(SYNVERSION, "\n");
	file->Write("totalHit:"); file->WriteLLong((!crashSave && !saveSelected) ? gHits->total.nbHit : 0, "\n");
	file->Write("totalDes:"); file->WriteLLong((!crashSave && !saveSelected) ? gHits->total.nbDesorbed : 0, "\n");
	file->Write("no_scans:"); file->WriteDouble((!crashSave && !saveSelected) ? worker->no_scans : 0, "\n");
	file->Write("totalLeak:"); file->WriteLLong((!crashSave && !saveSelected) ? gHits->nbLeakTotal : 0, "\n");
	file->Write("totalFlux:"); file->WriteDouble((!crashSave && !saveSelected) ? gHits->total.fluxAbs : 0, "\n");
	file->Write("totalPower:"); file->WriteDouble((!crashSave && !saveSelected) ? gHits->total.powerAbs : 0, "\n");
	file->Write("maxDes:"); file->WriteLLong((!crashSave && !saveSelected) ? tNbDesorptionMax : 0, "\n");
	file->Write("nbVertex:"); file->WriteInt(sh.nbVertex, "\n");
	file->Write("nbFacet:"); file->WriteInt(saveSelected ? nbSelected : sh.nbFacet, "\n");
	file->Write("nbSuper:"); file->WriteInt(sh.nbSuper, "\n");
	file->Write("nbFormula:"); file->WriteInt((!saveSelected) ? mApp->nbFormula : 0, "\n");
	file->Write("nbView:"); file->WriteInt(mApp->nbView, "\n");
	file->Write("nbSelection:"); file->WriteInt((!saveSelected) ? mApp->nbSelection : 0, "\n");

	file->Write("nbRegions:"); file->WriteInt((!saveSelected) ? (int)worker->regions.size() : 0, "\n");
	file->Write("PARfiles {\n");
	for (int i = 0; (!saveSelected) && (i < (int)worker->regions.size()); i++) { //write only the filenames without the path
		file->Write("  \"");
		char tmp[512];
		size_t foundDir = worker->regions[i].fileName.rfind('\\');
		if (foundDir != string::npos) {//found
			strcpy(tmp, &(worker->regions[i].fileName.c_str()[(int)foundDir + 1]));
		}
		else {//no "\" sign
			strcpy(tmp, worker->regions[i].fileName.c_str());
		}
		file->Write(tmp);
		file->Write("\"\n");
	}
	file->Write("}\n");

	file->Write("formulas {\n");
	for (int i = 0; (!saveSelected) && (i < mApp->nbFormula); i++) {
		file->Write("  \"");
		file->Write(mApp->formulas[i].parser->GetName());
		file->Write("\" \"");
		file->Write(mApp->formulas[i].parser->GetExpression());
		file->Write("\"\n");
	}
	file->Write("}\n");

	file->Write("views {\n");
	for (int i = 0; i < mApp->nbView; i++) {
		file->Write("  \"");
		file->Write(mApp->views[i].name);
		file->Write("\"\n");
		file->WriteInt(mApp->views[i].projMode, " ");
		file->WriteDouble(mApp->views[i].camAngleOx, " ");
		file->WriteDouble(mApp->views[i].camAngleOy, " ");
		file->WriteDouble(mApp->views[i].camDist, " ");
		file->WriteDouble(mApp->views[i].camOffset.x, " ");
		file->WriteDouble(mApp->views[i].camOffset.y, " ");
		file->WriteDouble(mApp->views[i].camOffset.z, " ");
		file->WriteInt(mApp->views[i].performXY, " ");
		file->WriteDouble(mApp->views[i].vLeft, " ");
		file->WriteDouble(mApp->views[i].vRight, " ");
		file->WriteDouble(mApp->views[i].vTop, " ");
		file->WriteDouble(mApp->views[i].vBottom, "\n");
	}
	file->Write("}\n");

	file->Write("selections {\n");
	for (int i = 0; (!saveSelected) && (i < mApp->nbSelection); i++) {
		file->Write("  \"");
		file->Write(mApp->selections[i].name);
		file->Write("\"\n ");
		file->WriteInt(mApp->selections[i].nbSel, "\n");
		for (int j = 0; j < mApp->selections[i].nbSel; j++) {
			file->Write("  ");
			file->WriteInt(mApp->selections[i].selection[j], "\n");
		}
		//file->Write("\n");
	}
	file->Write("}\n");

	file->Write("structures {\n");
	for (int i = 0; i < sh.nbSuper; i++) {
		file->Write("  \"");
		file->Write(strName[i]);
		file->Write("\"\n");
	}
	file->Write("}\n");
	//vertices
	prg->SetMessage("Writing vertices...");
	file->Write("vertices {\n");
	for (int i = 0; i < sh.nbVertex; i++) {
		prg->SetProgress(0.33*((double)i / (double)sh.nbVertex));
		file->Write("  ");
		file->WriteInt(i + 1, " ");
		file->WriteDouble(vertices3[i].x, " ");
		file->WriteDouble(vertices3[i].y, " ");
		file->WriteDouble(vertices3[i].z, "\n");
	}
	file->Write("}\n");

	//leaks
	prg->SetMessage("Writing leaks...");
	file->Write("leaks {\n");
	file->Write("  nbLeak:"); file->WriteInt((!crashSave && !saveSelected) ? *nbleakSave : 0, "\n");
	for (int i = 0; (i < *nbleakSave) && (!crashSave && !saveSelected); i++) {

		file->Write("  ");
		file->WriteInt(i, " ");
		file->WriteDouble((pleak + i)->pos.x, " ");
		file->WriteDouble((pleak + i)->pos.y, " ");
		file->WriteDouble((pleak + i)->pos.z, " ");

		file->WriteDouble((pleak + i)->dir.x, " ");
		file->WriteDouble((pleak + i)->dir.y, " ");
		file->WriteDouble((pleak + i)->dir.z, "\n");
	}
	file->Write("}\n");

	//hit cache (lines and dots)
	prg->SetMessage("Writing hit cache...");
	file->Write("hits {\n");
	file->Write("  nbHHit:"); file->WriteInt((!crashSave && !saveSelected) ? *nbHHitSave : 0, "\n");
	for (int i = 0; (i < *nbHHitSave) && (!crashSave && !saveSelected); i++) {

		file->Write("  ");
		file->WriteInt(i, " ");
		file->WriteDouble((pHits + i)->pos.x, " ");
		file->WriteDouble((pHits + i)->pos.y, " ");
		file->WriteDouble((pHits + i)->pos.z, " ");
		file->WriteDouble((pHits + i)->dF, " ");
		file->WriteDouble((pHits + i)->dP, " ");
		file->WriteInt((pHits + i)->type, "\n");
	}
	file->Write("}\n");

	//facets

	prg->SetMessage("Writing facets...");

	for (int i = 0, k = 0; i < sh.nbFacet; i++) {
		prg->SetProgress(0.33 + ((double)i / (double)sh.nbFacet) *0.33);
		if (!saveSelected || facets[i]->selected) { facets[i]->SaveSYN(file, worker->materials, k, crashSave); k++; }
	}

	prg->SetMessage("Writing profiles and spectrum...");
	SaveProfileSYN(file, dpHit, -1, saveSelected, crashSave);
	SaveSpectrumSYN(file, dpHit, -1, saveSelected, crashSave);

	///Save textures, for GEO file version 3+
	char tmp[256];
	file->Write("{textures}\n");

	file->Write("minHit_MC:"); file->WriteLLong(texCMin_MC, "\n");
	file->Write("maxHit_MC:"); file->WriteLLong(texCMax_MC, "\n");
	file->Write("minHit_flux:"); file->WriteDouble(texCMin_flux, "\n");
	file->Write("maxHit_flux:"); file->WriteDouble(texCMax_flux, "\n");
	file->Write("minHit_power:"); file->WriteDouble(texCMin_power, "\n");
	file->Write("maxHit_power:"); file->WriteDouble(texCMax_power, "\n");

	//Selections
	//SaveSelections();

	prg->SetMessage("Writing textures...");
	for (int i = 0; i < sh.nbFacet; i++) {
		prg->SetProgress(((double)i / (double)sh.nbFacet) *0.33 + 0.66);
		Facet *f = facets[i];
		if (f->hasMesh) {
			int profSize = (f->sh.isProfile) ? (PROFILE_SIZE*(sizeof(llong) + 2 * sizeof(double))) : 0;
			int nbE = f->sh.texHeight*f->sh.texWidth;
			llong *hits_MC;
			if (!crashSave && !saveSelected) hits_MC = (llong *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize));
			double *hits_flux;
			if (!crashSave && !saveSelected) hits_flux = (double *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + nbE*sizeof(llong)));
			double *hits_power;
			if (!crashSave && !saveSelected) hits_power = (double *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + nbE*(sizeof(llong) + sizeof(double))));

			//char tmp[256];
			sprintf(tmp, "texture_facet %d {\n", i + 1);
			file->Write(tmp);
			file->Write("width:"); file->WriteInt(f->sh.texWidth); file->Write(" height:"); file->WriteInt(f->sh.texHeight); file->Write("\n");
			for (iy = 0; iy < (f->sh.texHeight); iy++) {
				for (ix = 0; ix < (f->sh.texWidth); ix++) {
					int index = iy*(f->sh.texWidth) + ix;
					file->WriteLLong((!crashSave && !saveSelected) ? *(hits_MC + index) : 0, "\t");
					file->WriteDouble(f->GetMeshArea(index), "\t");
					file->WriteDouble((!crashSave && !saveSelected) ? *(hits_flux + index)*f->GetMeshArea(index) : 0, "\t");
					file->WriteDouble((!crashSave && !saveSelected) ? *(hits_power + index)*f->GetMeshArea(index) : 0, "\t");
				}
				file->Write("\n");
			}
			file->Write("}\n");
		}
	}

	if (!crashSave && !saveSelected) ReleaseDataport(dpHit);

}

std::vector<std::string> Geometry::LoadSYN(FileReader *file, GLProgress *prg, LEAK *pleak, int *nbleak, HIT *pHits, int *nbHHit, int *version) {

	prg->SetMessage("Clearing current geometry...");
	Clear();

	std::vector<std::string> result;

	// Globals
	char tmp[512];
	prg->SetMessage("Reading SYN file header...");
	file->ReadKeyword("version"); file->ReadKeyword(":");
	*version = file->ReadInt();
	if (*version > SYNVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported SYN version V%d", *version);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit"); file->ReadKeyword(":");
	tNbHit = file->ReadLLong();
	file->ReadKeyword("totalDes"); file->ReadKeyword(":");
	tNbDesorption = file->ReadLLong();
	if (*version >= 6) {
		file->ReadKeyword("no_scans"); file->ReadKeyword(":");
		loaded_no_scans = file->ReadDouble();
	}
	else loaded_no_scans = 0;
	file->ReadKeyword("totalLeak"); file->ReadKeyword(":");
	tNbLeak = file->ReadInt();
	if (*version > 2) {
		file->ReadKeyword("totalFlux"); file->ReadKeyword(":");
		tFlux = file->ReadDouble();
		file->ReadKeyword("totalPower"); file->ReadKeyword(":");
		tPower = file->ReadDouble();
	}
	tNbAbsorption = 0;
	distTraveledTotal = 0.0;
	file->ReadKeyword("maxDes"); file->ReadKeyword(":");
	tNbDesorptionMax = file->ReadLLong();
	file->ReadKeyword("nbVertex"); file->ReadKeyword(":");
	sh.nbVertex = file->ReadInt();
	file->ReadKeyword("nbFacet"); file->ReadKeyword(":");
	sh.nbFacet = file->ReadInt();
	file->ReadKeyword("nbSuper"); file->ReadKeyword(":");
	sh.nbSuper = file->ReadInt();
	int nbF = 0; std::vector<std::vector<string>> loadFormulas;
	int nbV = 0;

	file->ReadKeyword("nbFormula"); file->ReadKeyword(":");
	nbF = file->ReadInt(); loadFormulas.reserve(nbF);
	file->ReadKeyword("nbView"); file->ReadKeyword(":");
	nbV = file->ReadInt();
	int nbS = 0;
	file->ReadKeyword("nbSelection"); file->ReadKeyword(":");
	nbS = file->ReadInt();

	if (*version > 1) {
		file->ReadKeyword("nbRegions"); file->ReadKeyword(":");
		int nbR = file->ReadInt();
		result.reserve(nbR);

		file->ReadKeyword("PARfiles"); file->ReadKeyword("{");
		for (int i = 0; i < nbR; i++) {
			result.push_back(file->ReadString());
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("formulas"); file->ReadKeyword("{");
	for (int i = 0; i < nbF; i++) {
		char tmpName[256];
		char tmpExpr[512];
		strcpy(tmpName, file->ReadString());
		strcpy(tmpExpr, file->ReadString());
		//mApp->AddFormula(tmpName, tmpExpr); //parse after selection groups are loaded
		std::vector<string> newFormula;
		newFormula.push_back(tmpName);
		newFormula.push_back(tmpExpr);
		loadFormulas.push_back(newFormula);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("views"); file->ReadKeyword("{");
	for (int i = 0; i < nbV; i++) {
		char tmpName[256];
		AVIEW v;
		strcpy(tmpName, file->ReadString());
		v.projMode = file->ReadInt();
		v.camAngleOx = file->ReadDouble();
		v.camAngleOy = file->ReadDouble();
		v.camDist = file->ReadDouble();
		v.camOffset.x = file->ReadDouble();
		v.camOffset.y = file->ReadDouble();
		v.camOffset.z = file->ReadDouble();
		v.performXY = file->ReadInt();
		v.vLeft = file->ReadDouble();
		v.vRight = file->ReadDouble();
		v.vTop = file->ReadDouble();
		v.vBottom = file->ReadDouble();
		mApp->AddView(tmpName, v);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("selections"); file->ReadKeyword("{");
	for (int i = 0; i < nbS; i++) {
		ASELECTION s;
		char tmpName[256];
		strcpy(tmpName, file->ReadString());
		s.name = _strdup(tmpName);
		s.nbSel = file->ReadInt();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));

		for (int j = 0; j < s.nbSel; j++) {
			s.selection[j] = file->ReadInt();
		}
		mApp->AddSelection(s.name, s);
	}
	file->ReadKeyword("}");

	for (int i = 0; i < nbF; i++) { //parse formulas now that selection groups are loaded
		mApp->AddFormula(loadFormulas[i][0].c_str(), loadFormulas[i][1].c_str());
	}

	file->ReadKeyword("structures"); file->ReadKeyword("{");
	for (int i = 0; i < sh.nbSuper; i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		sprintf(tmp, "%s.txt", strName[i]);
		strFileName[i] = _strdup(tmp);
	}
	file->ReadKeyword("}");

	// Allocate memory
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memset(vertices3, 0, sh.nbVertex * sizeof(VERTEX3D));

	// Read vertices
	prg->SetMessage("Reading vertices...");
	file->ReadKeyword("vertices"); file->ReadKeyword("{");
	for (int i = 0; i < sh.nbVertex; i++) {
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1) throw Error(file->MakeError("Wrong vertex index !"));
		vertices3[i].x = file->ReadDouble();
		vertices3[i].y = file->ReadDouble();
		vertices3[i].z = file->ReadDouble();
		vertices3[i].selected = FALSE;
	}
	file->ReadKeyword("}");
	prg->SetMessage("Reading leaks and hits...");
	// Read leaks
	file->ReadKeyword("leaks"); file->ReadKeyword("{");
	file->ReadKeyword("nbLeak"); file->ReadKeyword(":");
	*nbleak = file->ReadInt();
	for (int i = 0; i < *nbleak; i++) {
		int idx = file->ReadInt();
		if (idx != i) throw Error(file->MakeError("Wrong leak index !"));
		(pleak + i)->pos.x = file->ReadDouble();
		(pleak + i)->pos.y = file->ReadDouble();
		(pleak + i)->pos.z = file->ReadDouble();

		(pleak + i)->dir.x = file->ReadDouble();
		(pleak + i)->dir.y = file->ReadDouble();
		(pleak + i)->dir.z = file->ReadDouble();
	}
	file->ReadKeyword("}");

	// Read hit cache
	file->ReadKeyword("hits"); file->ReadKeyword("{");
	file->ReadKeyword("nbHHit"); file->ReadKeyword(":");
	*nbHHit = file->ReadInt();
	for (int i = 0; i < *nbHHit; i++) {
		int idx = file->ReadInt();
		if (idx != i) throw Error(file->MakeError("Wrong hit cache index !"));
		(pHits + i)->pos.x = file->ReadDouble();
		(pHits + i)->pos.y = file->ReadDouble();
		(pHits + i)->pos.z = file->ReadDouble();
		(pHits + i)->dF = file->ReadDouble();
		(pHits + i)->dP = file->ReadDouble();
		(pHits + i)->type = file->ReadInt();
	}
	file->ReadKeyword("}");
	// Read facets
	prg->SetMessage("Reading facets...");
	for (int i = 0; i < sh.nbFacet; i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nbI = file->ReadInt();
		if (nbI < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %d has only %d vertices. ", i, nbI);
			throw Error(errMsg);
		}
		prg->SetProgress((float)i / sh.nbFacet);
		facets[i] = new Facet(nbI);
		facets[i]->LoadSYN(file, mApp->worker.materials, *version, sh.nbVertex);
		file->ReadKeyword("}");
	}

	prg->SetMessage("Initalizing geometry and building mesh...");
	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;
	UpdateName(file);

	// Update mesh
	prg->SetMessage("Building mesh...");
	for (int i = 0; i < sh.nbFacet; i++) {
		double p = (double)i / (double)sh.nbFacet;
		prg->SetProgress(p);
		Facet *f = facets[i];
		if (!f->SetTexture(f->sh.texWidthD, f->sh.texHeightD, f->hasMesh)) {
			char errMsg[512];
			sprintf(errMsg, "Not enough memory to build mesh on Facet %d. ", i + 1);
			throw Error(errMsg);
		}
		BuildFacetList(f);
		double nU = Norme(f->sh.U);
		f->tRatio = f->sh.texWidthD / nU;
	}
	return result;
}

void Geometry::LoadProfileSYN(FileReader *file, Dataport *dpHit) { //profiles and spectrums
	AccessDataport(dpHit);
	BYTE *buffer = (BYTE *)dpHit->buff;
	file->ReadKeyword("profiles"); file->ReadKeyword("{");
	// Profiles
	int nbProfile;
	file->ReadKeyword("number"); file->ReadKeyword(":"); nbProfile = file->ReadInt();
	int *profileFacet = (int *)malloc((nbProfile)*sizeof(int));
	file->ReadKeyword("facets"); file->ReadKeyword(":");
	for (int i = 0; i < nbProfile; i++)
		profileFacet[i] = file->ReadInt();
	for (int j = 0; j < PROFILE_SIZE; j++) {
		for (int i = 0; i < nbProfile; i++) {
			Facet *f = GetFacet(profileFacet[i]);
			llong *profilePtr_MC = (llong *)(buffer + f->sh.hitOffset + sizeof(SHHITS));
			double *profilePtr_flux = (double *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + PROFILE_SIZE*sizeof(llong));
			double *profilePtr_power = (double *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + PROFILE_SIZE*(sizeof(llong) + sizeof(double)));

			profilePtr_MC[j] = file->ReadLLong();
			profilePtr_flux[j] = file->ReadDouble();
			profilePtr_power[j] = file->ReadDouble();
		}
	}
	file->ReadKeyword("}");
	ReleaseDataport(dpHit);
	SAFE_FREE(profileFacet);
}
void Geometry::LoadSpectrumSYN(FileReader *file, Dataport *dpHit) { //spectrums and spectrums
	AccessDataport(dpHit);
	BYTE *buffer = (BYTE *)dpHit->buff;
	file->ReadKeyword("spectrums"); file->ReadKeyword("{");
	// Spectrums
	int nbSpectrum;
	file->ReadKeyword("number"); file->ReadKeyword(":"); nbSpectrum = file->ReadInt();
	int *spectrumFacet = (int *)malloc((nbSpectrum)*sizeof(int));
	file->ReadKeyword("facets"); file->ReadKeyword(":");
	for (int i = 0; i < nbSpectrum; i++)
		spectrumFacet[i] = file->ReadInt();

	for (int j = 0; j < PROFILE_SIZE; j++) {
		for (int i = 0; i < nbSpectrum; i++) {
			Facet *f = GetFacet(spectrumFacet[i]);
			int profileSize = (f->sh.isProfile) ? PROFILE_SIZE*(2 * sizeof(double) + sizeof(llong)) : 0;
			int textureSize = f->sh.texWidth*f->sh.texHeight*(2 * sizeof(double) + sizeof(llong));
			int directionSize = f->sh.countDirection*f->sh.texWidth*f->sh.texHeight*sizeof(VHIT);
			double *shSpectrum_fluxwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize
				+ textureSize + directionSize));
			double *shSpectrum_powerwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize
				+ textureSize + directionSize + SPECTRUM_SIZE*sizeof(double)));

			shSpectrum_fluxwise[j] = file->ReadDouble();
			shSpectrum_powerwise[j] = file->ReadDouble();
		}
	}
	file->ReadKeyword("}");
	ReleaseDataport(dpHit);
	SAFE_FREE(spectrumFacet);
}
void Geometry::SaveProfileSYN(FileWriter *file, Dataport *dpHit, int super, BOOL saveSelected, BOOL crashSave) {
	//Profiles

	BYTE *buffer;
	if (!crashSave && !saveSelected) buffer = (BYTE *)dpHit->buff;
	file->Write("profiles {\n");
	// Profiles
	int nbProfile = 0;
	int *profileFacet = (int *)malloc((sh.nbFacet)*sizeof(int));
	for (int i = 0; i < sh.nbFacet; i++)
		if ((!saveSelected && !crashSave) && facets[i]->sh.isProfile)
			profileFacet[nbProfile++] = i;
	file->Write(" number: "); file->WriteInt(nbProfile, "\n");
	file->Write(" facets: ");
	for (int i = 0; i < nbProfile; i++)  //doesn't execute when crashSave or saveSelected...
		file->WriteInt(profileFacet[i], "\t");
	file->Write("\n");
	for (int j = 0; j < PROFILE_SIZE; j++) {
		for (int i = 0; i < nbProfile; i++) { //doesn't execute when crashSave or saveSelected...
			Facet *f = GetFacet(profileFacet[i]);
			llong *profilePtr_MC;
			if (!crashSave) profilePtr_MC = (llong *)(buffer + f->sh.hitOffset + sizeof(SHHITS));
			double *profilePtr_flux;
			if (!crashSave) profilePtr_flux = (double *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + PROFILE_SIZE*sizeof(llong));
			double *profilePtr_power;
			if (!crashSave) profilePtr_power = (double *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + PROFILE_SIZE*(sizeof(llong) + sizeof(double)));

			file->WriteLLong((!crashSave) ? profilePtr_MC[j] : 0); file->Write("\t");
			file->WriteDouble((!crashSave) ? profilePtr_flux[j] : 0); file->Write("\t");
			file->WriteDouble((!crashSave) ? profilePtr_power[j] : 0); file->Write("\t");
		}
		file->Write("\n");
	}
	file->Write("}\n");
	SAFE_FREE(profileFacet);
}

void Geometry::SaveSpectrumSYN(FileWriter *file, Dataport *dpHit, int super, BOOL saveSelected, BOOL crashSave) {
	//Spectrums

	BYTE *buffer;
	if (!crashSave && !saveSelected) buffer = (BYTE *)dpHit->buff;
	file->Write("spectrums {\n");
	// Profiles
	int nbSpectrum = 0;
	int *spectrumFacet = (int *)malloc((sh.nbFacet)*sizeof(int));
	for (int i = 0; i < sh.nbFacet; i++)
		if ((!saveSelected && !crashSave) && facets[i]->sh.hasSpectrum)
			spectrumFacet[nbSpectrum++] = i;
	file->Write(" number: "); file->WriteInt(nbSpectrum, "\n");
	file->Write(" facets: ");
	for (int i = 0; i < nbSpectrum; i++)   //doesn't execute when crashSave or saveSelected...
		file->WriteInt(spectrumFacet[i], "\t");
	file->Write("\n");
	for (int j = 0; j < SPECTRUM_SIZE; j++) {
		for (int i = 0; i < nbSpectrum; i++) {  //doesn't execute when crashSave or saveSelected...
			Facet *f = GetFacet(spectrumFacet[i]);
			int profileSize = (f->sh.isProfile) ? PROFILE_SIZE*(2 * sizeof(double) + sizeof(llong)) : 0;
			int textureSize = f->sh.texWidth*f->sh.texHeight*(2 * sizeof(double) + sizeof(llong));
			int directionSize = f->sh.countDirection*f->sh.texWidth*f->sh.texHeight*sizeof(VHIT);
			double *shSpectrum_fluxwise;
			if (!crashSave) shSpectrum_fluxwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize
				+ textureSize + directionSize));
			double *shSpectrum_powerwise;
			if (!crashSave) shSpectrum_powerwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize
				+ textureSize + directionSize + SPECTRUM_SIZE*sizeof(double)));

			file->WriteDouble((!crashSave) ? shSpectrum_fluxwise[j] : 0); file->Write("\t");
			file->WriteDouble((!crashSave) ? shSpectrum_powerwise[j] : 0); file->Write("\t");
		}
		file->Write("\n");
	}
	file->Write("}\n");
	SAFE_FREE(spectrumFacet);
}

//Temporary placeholders
void Geometry::SaveXML_geometry(pugi::xml_node saveDoc, Worker *work, GLProgress *prg, BOOL saveSelected){}
BOOL Geometry::SaveXML_simustate(pugi::xml_node saveDoc, Worker *work, BYTE *buffer, SHGHITS *gHits, int nbLeakSave, int nbHHitSave,
	LEAK *pLeak, HIT *pHits, GLProgress *prg, BOOL saveSelected){
	return FALSE;
}
void Geometry::LoadXML_geom(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg){
	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	//mApp->ClearFormula();
	Clear();
	xml_node geomNode = loadXML.child("Geometry");

	//Vertices
	sh.nbVertex = geomNode.child("Vertices").select_nodes("Vertex").size();

	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	int idx = 0;
	for (xml_node vertex : geomNode.child("Vertices").children("Vertex")) {
		vertices3[idx].x = vertex.attribute("x").as_double();
		vertices3[idx].y = vertex.attribute("y").as_double();
		vertices3[idx].z = vertex.attribute("z").as_double();
		vertices3[idx].selected = FALSE;
		idx++;

	}

	//Structures
	sh.nbSuper = geomNode.child("Structures").select_nodes("Structure").size();
	idx = 0;
	for (xml_node structure : geomNode.child("Structures").children("Structure")) {
		strName[idx] = _strdup(structure.attribute("name").value());
		// For backward compatibilty with STR
		char tmp[256];
		sprintf(tmp, "%s.txt", strName[idx]);
		strFileName[sh.nbSuper + idx] = _strdup(tmp);
		idx++;
	}

	//Parameters (needs to precede facets)
	xml_node simuParamNode = loadXML.child("SynradSimuSettings");
	BOOL isSynradFile = (simuParamNode != NULL); //if no "SynradSimuSettings" node, it's a Molflow XML file

	if (isSynradFile) {
		//Nothing yet...
	}

	//Facets
	sh.nbFacet = geomNode.child("Facets").select_nodes("Facet").size();
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	idx = 0;
	for (xml_node facetNode : geomNode.child("Facets").children("Facet")) {
		size_t nbIndex = facetNode.child("Indices").select_nodes("Indice").size();
		if (nbIndex < 3) {
			char errMsg[128];
			sprintf(errMsg, "Facet %d has only %d vertices. ", idx + 1, nbIndex);
			throw Error(errMsg);
		}

		facets[idx] = new Facet(nbIndex);
		facets[idx]->LoadXML(facetNode, sh.nbVertex, isSynradFile, isSynradFile);

		if (isSynradFile) {
			//Nothing yet...
		}
		idx++;
	}

	xml_node interfNode = loadXML.child("Interface");

	xml_node selNode = interfNode.child("Selections");
	//int nbS = selNode.select_nodes("Selection").size();

	for (xml_node sNode : selNode.children("Selection")) {
		ASELECTION s;
		s.name = _strdup(sNode.attribute("name").as_string());
		s.nbSel = sNode.select_nodes("selItem").size();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));
		idx = 0;
		for (xml_node iNode : sNode.children("selItem"))
			s.selection[idx++] = iNode.attribute("facet").as_int();
		mApp->AddSelection(s.name, s);
	}

	xml_node viewNode = interfNode.child("Views");
	for (xml_node newView : selNode.children("View")) {
		AVIEW v;
		v.name = _strdup(newView.attribute("name").as_string());
		v.projMode = newView.attribute("projMode").as_int();
		v.camAngleOx = newView.attribute("camAngleOx").as_double();
		v.camAngleOy = newView.attribute("camAngleOy").as_double();
		v.camDist = newView.attribute("camDist").as_double();
		v.camOffset.x = newView.attribute("camOffset.x").as_double();
		v.camOffset.y = newView.attribute("camOffset.y").as_double();
		v.camOffset.z = newView.attribute("camOffset.z").as_double();
		v.performXY = newView.attribute("performXY").as_int();
		v.vLeft = newView.attribute("vLeft").as_double();
		v.vRight = newView.attribute("vRight").as_double();
		v.vTop = newView.attribute("vTop").as_double();
		v.vBottom = newView.attribute("vBottom").as_double();
		mApp->AddView(v.name, v);
	}


	if (isSynradFile) {
		xml_node formulaNode = interfNode.child("Formulas");
		for (xml_node newFormula : formulaNode.children("Formula")) {
			mApp->AddFormula(newFormula.attribute("name").as_string(),
				newFormula.attribute("expression").as_string());
		}
	}

	InitializeGeometry(); //Includes Buildgllist

	isLoaded = TRUE;

	// Update mesh
	progressDlg->SetMessage("Building mesh...");
	for (int i = 0; i < sh.nbFacet; i++) {
		double p = (double)i / (double)sh.nbFacet;

		progressDlg->SetProgress(p);
		Facet *f = facets[i];
		if (!f->SetTexture(f->sh.texWidthD, f->sh.texHeightD, f->hasMesh)) {
			char errMsg[512];
			sprintf(errMsg, "Not enough memory to build mesh on Facet %d. ", i + 1);
			throw Error(errMsg);
		}
		BuildFacetList(f);
		double nU = Norme(f->sh.U);
		f->tRatio = f->sh.texWidthD / nU;
	}
}
void Geometry::InsertXML(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg, BOOL newStr){
	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	//mApp->ClearFormula();
	//Clear();
	int structId = viewStruct;
	if (structId == -1) structId = 0;
	UnSelectAll();

	xml_node geomNode = loadXML.child("Geometry");
	//Vertices
	int nbNewVertex = geomNode.child("Vertices").select_nodes("Vertex").size();
	int nbNewFacets = geomNode.child("Facets").select_nodes("Facet").size();

	// reallocate memory
	facets = (Facet **)realloc(facets, (nbNewFacets + sh.nbFacet) * sizeof(Facet **));
	memset(facets + sh.nbFacet, 0, nbNewFacets * sizeof(Facet *));

	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex + sh.nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3, vertices3, (sh.nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3 + sh.nbVertex, 0, nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(vertices3);
	vertices3 = tmp_vertices3;

	// Read geometry vertices
	int idx = sh.nbVertex;
	for (xml_node vertex : geomNode.child("Vertices").children("Vertex")) {
		vertices3[idx].x = vertex.attribute("x").as_double();
		vertices3[idx].y = vertex.attribute("y").as_double();
		vertices3[idx].z = vertex.attribute("z").as_double();
		vertices3[idx].selected = FALSE;
		idx++;
	}

	//Structures
	int nbNewSuper = geomNode.child("Structures").select_nodes("Structure").size();
	idx = 0;
	for (xml_node structure : geomNode.child("Structures").children("Structure")) {
		strName[sh.nbSuper + idx] = _strdup(structure.attribute("name").value());
		// For backward compatibilty with STR
		char tmp[256];
		sprintf(tmp, "%s.txt", strName[idx]);
		strFileName[sh.nbSuper + idx] = _strdup(tmp);
		idx++;
	}

	//Parameters (needs to precede facets)
	xml_node simuParamNode = loadXML.child("SynradSimuSettings");
	BOOL isSynradFile = (simuParamNode != NULL); //if no "SynradSimuSettings" node, it's a Molflow XML file

	if (isSynradFile) {
		//Nothing yet...
	}

	//Facets
	idx = sh.nbFacet;
	for (xml_node facetNode : geomNode.child("Facets").children("Facet")) {
		int nbIndex = facetNode.child("Indices").select_nodes("Indice").size();
		if (nbIndex < 3) {
			char errMsg[128];
			sprintf(errMsg, "Facet %d has only %d vertices. ", idx + 1, nbIndex);
			throw Error(errMsg);
		}
		facets[idx] = new Facet(nbIndex);
		facets[idx]->LoadXML(facetNode, sh.nbVertex + nbNewVertex, isSynradFile, sh.nbVertex);
		facets[idx]->selected = TRUE;
		facets[idx]->sh.superIdx += structId; //offset structure
		if (facets[idx]->sh.superDest>0) facets[idx]->sh.superDest += structId;
		if (isSynradFile) {
			//Nothing yet...
		}
		idx++;
	}

	xml_node interfNode = loadXML.child("Interface");
	xml_node selNode = interfNode.child("Selections");
	//int nbS = selNode.select_nodes("Selection").size();

	for (xml_node sNode : selNode.children("Selection")) {
		ASELECTION s;
		s.name = _strdup(sNode.attribute("name").as_string());
		s.nbSel = sNode.select_nodes("selItem").size();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));
		idx = 0;
		for (xml_node iNode : sNode.children("selItem"))
			s.selection[idx++] = iNode.attribute("facet").as_int() + sh.nbFacet; //offset selection numbers
		mApp->AddSelection(s.name, s);
	}

	xml_node viewNode = interfNode.child("Views");
	for (xml_node newView : selNode.children("View")) {
		AVIEW v;
		v.name = _strdup(newView.attribute("name").as_string());
		v.projMode = newView.attribute("projMode").as_int();
		v.camAngleOx = newView.attribute("camAngleOx").as_double();
		v.camAngleOy = newView.attribute("camAngleOy").as_double();
		v.camDist = newView.attribute("camDist").as_double();
		v.camOffset.x = newView.attribute("camOffset.x").as_double();
		v.camOffset.y = newView.attribute("camOffset.y").as_double();
		v.camOffset.z = newView.attribute("camOffset.z").as_double();
		v.performXY = newView.attribute("performXY").as_int();
		v.vLeft = newView.attribute("vLeft").as_double();
		v.vRight = newView.attribute("vRight").as_double();
		v.vTop = newView.attribute("vTop").as_double();
		v.vBottom = newView.attribute("vBottom").as_double();
		mApp->AddView(v.name, v);
	}


	sh.nbVertex += nbNewVertex;
	sh.nbFacet += nbNewFacets; //formulas can refer to newly inserted facets

	if (isSynradFile) {
		xml_node formulaNode = interfNode.child("Formulas");
		for (xml_node newFormula : formulaNode.children("Formula")) {
			char tmpExpr[512];
			strcpy(tmpExpr, newFormula.attribute("expression").as_string());
			mApp->OffsetFormula(tmpExpr, sh.nbFacet);
			mApp->AddFormula(newFormula.attribute("name").as_string(),
				tmpExpr);
		}
	}

	if (newStr) sh.nbSuper += nbNewSuper;
	else if (sh.nbSuper < structId + nbNewSuper) sh.nbSuper = structId + nbNewSuper;
	InitializeGeometry(); //Includes Buildgllist
	isLoaded = TRUE;

}
BOOL Geometry::LoadXML_simustate(pugi::xml_node loadXML, Dataport *dpHit, Worker *work, GLProgress *progressDlg){ return FALSE; }