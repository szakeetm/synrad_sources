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

using namespace pugi;
extern SynRad *mApp;

// -----------------------------------------------------------

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
	texMin_MC=0;
	texMax_MC=1;
	texMin_flux=0.0;
	texMax_flux=1.0;
	texMin_power=0.0;
	texMax_power=1.0;
	autoNorme = TRUE;
	centerNorme = TRUE;
	normeRatio = 1.0f;
	texAutoScale=TRUE;
	texLogScale=FALSE;

	sh.nbSuper=0;
	viewStruct = -1;
	strcpy(strPath,"");
	Clear();

}

// -----------------------------------------------------------

Geometry::~Geometry() {
	Clear();
}

// -----------------------------------------------------------

void Geometry::Clear() {

	// Free memory
	if(facets) {
		for(int i=0;i<sh.nbFacet;i++)
			SAFE_DELETE(facets[i]);
		free(facets);
	}
	if(vertices3) free(vertices3);
	for(int i=0;i<sh.nbSuper;i++) {
		SAFE_FREE(strName[i]);
		SAFE_FREE(strFileName[i]);
	}
	memset(strName,0,MAX_SUPERSTR * sizeof(char *));
	memset(strFileName,0,MAX_SUPERSTR * sizeof(char *));
	DeleteGLLists(TRUE,TRUE);

	// Init default
	facets = NULL;         // Facets array
	vertices3 = NULL;      // Facets vertices in (x,y,z) space
	sh.nbFacet = 0;        // Number of facets
	sh.nbVertex = 0;       // Number of vertex
	isLoaded = FALSE;      // isLoaded flag
	texCMin_MC=0;             // Current minimum
	texCMax_MC=1;             // Current maximum
	texCMin_flux=texCMin_power=0.0;
	texCMax_flux=texCMax_power=1.0;
	sh.nbSuper=0;          // Structure number
	nbSelected=0;          // Number of selected facet

	memset(lineList,0,sizeof(lineList));
	memset(strName,0,sizeof(strName));
	memset(strFileName,0,sizeof(strFileName));

	// Init OpenGL material
	memset (&whiteMaterial, 0, sizeof (GLMATERIAL));
	whiteMaterial.Diffuse.r = 0.9f;
	whiteMaterial.Diffuse.g = 0.9f;
	whiteMaterial.Diffuse.b = 0.9f;
	whiteMaterial.Ambient.r = 0.9f;
	whiteMaterial.Ambient.g = 0.9f;
	whiteMaterial.Ambient.b = 0.9f;

	memset (&fillMaterial, 0, sizeof (GLMATERIAL));

	
	fillMaterial.Diffuse.r = 0.6f;
	fillMaterial.Diffuse.g = 0.65f;
	fillMaterial.Diffuse.b = 0.65f;
	fillMaterial.Ambient.r = 0.45f;
	fillMaterial.Ambient.g = 0.41f;
	fillMaterial.Ambient.b = 0.41f;

	memset (&arrowMaterial, 0, sizeof (GLMATERIAL));
	arrowMaterial.Diffuse.r = 0.4f;
	arrowMaterial.Diffuse.g = 0.2f;
	arrowMaterial.Diffuse.b = 0.0f;
	arrowMaterial.Ambient.r = 0.6f;
	arrowMaterial.Ambient.g = 0.3f;
	arrowMaterial.Ambient.b = 0.0f;

	nbSelectedHist = 0;    // Selection history
	nbSelectedHistVertex = 0;

}

void Geometry::CalculateFacetParam(int facet) {

	Facet *f = facets[facet];

	// Calculate facet normal
	VERTEX3D p0 = vertices3[facets[facet]->indices[0]];
	VERTEX3D v1;
	VERTEX3D v2;
	BOOL consecutive = TRUE;
	int ind=2;

	// TODO: Handle possible collinear consequtive vectors
	int i0 = facets[facet]->indices[0];
	int i1 = facets[facet]->indices[1];
	while (ind<f->sh.nbIndex && consecutive) {
		int i2 = facets[facet]->indices[ind++];

		Sub(&v1 , vertices3+i1 , vertices3+i0); // v1 = P0P1
		Sub(&v2 , vertices3+i2 , vertices3+i1); // v2 = P1P2
		Cross(&(f->sh.N),&v1,&v2);              // Cross product
		consecutive=(Norme(&(f->sh.N))<1e-11);
	}
	f->collinear=consecutive; //mark for later that this facet was on a line
	Normalize(&(f->sh.N));                  // Normalize

	// Calculate Axis Aligned Bounding Box
	f->sh.bb.min.x =  1e100;
	f->sh.bb.min.y =  1e100;
	f->sh.bb.min.z =  1e100;
	f->sh.bb.max.x = -1e100;
	f->sh.bb.max.y = -1e100;
	f->sh.bb.max.z = -1e100;

	for(int i=0;i<f->sh.nbIndex;i++) {
		VERTEX3D p = vertices3[f->indices[i]];
		if( p.x < f->sh.bb.min.x ) f->sh.bb.min.x = p.x;
		if( p.y < f->sh.bb.min.y ) f->sh.bb.min.y = p.y;
		if( p.z < f->sh.bb.min.z ) f->sh.bb.min.z = p.z;
		if( p.x > f->sh.bb.max.x ) f->sh.bb.max.x = p.x;
		if( p.y > f->sh.bb.max.y ) f->sh.bb.max.y = p.y;
		if( p.z > f->sh.bb.max.z ) f->sh.bb.max.z = p.z;
	}

	// Facet center (AABB center)
	f->sh.center.x = (f->sh.bb.max.x + f->sh.bb.min.x)/2.0;
	f->sh.center.y = (f->sh.bb.max.y + f->sh.bb.min.y)/2.0;
	f->sh.center.z = (f->sh.bb.max.z + f->sh.bb.min.z)/2.0;

	// Plane equation
	double A = f->sh.N.x;
	double B = f->sh.N.y;
	double C = f->sh.N.z;
	double D = -Dot(&(f->sh.N),&p0);

	// Facet planeity
	int nb = f->sh.nbIndex;
	double max = 0.0;
	for(int i=1;i<nb;i++) {
		VERTEX3D p = vertices3[f->indices[i]];
		double d = A * p.x + B * p.y + C * p.z + D;
		if(fabs(d)>fabs(max)) max = d;
	}

	// Plane coef
	f->a = A;
	f->b = B;
	f->c = C;
	f->d = D;
	f->err = max;

}






void Geometry::InitializeGeometry(int facet_number,BOOL BBOnly) {

	// Perform various precalculation here for a faster simulation.

	//GLProgress *initGeoPrg = new GLProgress("Initializing geometry...","Please wait");
	//initGeoPrg->SetProgress(0.0);
	//initGeoPrg->SetVisible(TRUE);
	if (facet_number==-1) { //bounding box for all vertices
		bb.min.x =  1e100;
		bb.min.y =  1e100;
		bb.min.z =  1e100;
		bb.max.x = -1e100;
		bb.max.y = -1e100;
		bb.max.z = -1e100;

		// Axis Aligned Bounding Box
		for(int i=0;i<sh.nbVertex;i++) {
			VERTEX3D p = vertices3[i];
			if(!(vertices3[i].selected==FALSE || vertices3[i].selected==TRUE)) vertices3[i].selected=FALSE; //initialize selection
			if( p.x < bb.min.x ) bb.min.x = p.x;
			if( p.y < bb.min.y ) bb.min.y = p.y;
			if( p.z < bb.min.z ) bb.min.z = p.z;
			if( p.x > bb.max.x ) bb.max.x = p.x;
			if( p.y > bb.max.y ) bb.max.y = p.y;
			if( p.z > bb.max.z ) bb.max.z = p.z;
		}

		Worker *worker=&(mApp->worker);
		for (int i=0;i<(int)worker->regions.size();i++) {
			if (worker->regions[i].AABBmin.x<bb.min.x) bb.min.x=worker->regions[i].AABBmin.x;
			if (worker->regions[i].AABBmin.y<bb.min.y) bb.min.y=worker->regions[i].AABBmin.y;
			if (worker->regions[i].AABBmin.z<bb.min.z) bb.min.z=worker->regions[i].AABBmin.z;
			if (worker->regions[i].AABBmax.x>bb.max.x) bb.max.x=worker->regions[i].AABBmax.x;
			if (worker->regions[i].AABBmax.y>bb.max.y) bb.max.y=worker->regions[i].AABBmax.y;
			if (worker->regions[i].AABBmax.z>bb.max.z) bb.max.z=worker->regions[i].AABBmax.z;
		}
	} else { //bounding box only for the changed facet
		for(int i=0;i<facets[facet_number]->sh.nbIndex;i++) {
			VERTEX3D p = vertices3[facets[facet_number]->indices[i]];
			//if(!(vertices3[i].selected==FALSE || vertices3[i].selected==TRUE)) vertices3[i].selected=FALSE; //initialize selection
			if( p.x < bb.min.x ) bb.min.x = p.x;
			if( p.y < bb.min.y ) bb.min.y = p.y;
			if( p.z < bb.min.z ) bb.min.z = p.z;
			if( p.x > bb.max.x ) bb.max.x = p.x;
			if( p.y > bb.max.y ) bb.max.y = p.y;
			if( p.z > bb.max.z ) bb.max.z = p.z;
		}
	}

	center.x = (bb.max.x + bb.min.x)/2.0;
	center.y = (bb.max.y + bb.min.y)/2.0;
	center.z = (bb.max.z + bb.min.z)/2.0;

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
	for(int i=0;i<sh.nbFacet;i++) {
		//initGeoPrg->SetProgress((double)i/(double)sh.nbFacet);
		if((facet_number==-1)||(i==facet_number)) { //permits to initialize only one facet
			// Main facet params
			CalculateFacetParam(i);

			// Current facet
			Facet *f = facets[i];

			/*
			// Search longest edge
			double dMax = 0.0;
			int    i1Max,i2Max;
			for(int j=0;j<f->sh.nbIndex;j++) {
				int j1 = IDX(j+1,f->sh.nbIndex);
				int i1 = f->indices[j];
				int i2 = f->indices[j1];
				double xl = (vertices3[i1].x - vertices3[i2].x);
				double yl = (vertices3[i1].y - vertices3[i2].y);
				double zl = (vertices3[i1].z - vertices3[i2].z);
				double l = DOT3( xl,yl,zl,xl,yl,zl ); //distance of vertices at i1 and i2
				if( l>dMax + 1e-8 ) { //disregard differences below double precision
					dMax = l;
					i1Max = i1;
					i2Max = i2;
				}
			}

			// First vertex
			int i0 = f->indices[0];
			VERTEX3D p0 = vertices3[i0];
			*/

			VERTEX3D p0 = vertices3[f->indices[0]];
			VERTEX3D p1 = vertices3[f->indices[1]];

			VERTEX3D U,V;

			/* OLD fashion (no longueur used)
			// Intersection line (directed by U and including p0) 
			U.x = f->c;
			U.y = 0.0;
			U.z = -f->a;
			double nU = Norme(&U);

			if( IS_ZERO(nU) ) {
			// Plane parallel to (O,x,z)
			U.x = 1.0;
			U.y = 0.0;
			U.z = 0.0;
			} else {
			ScalarMult(&U,1.0/nU); // Normalize U
			}
			*/
			/*
			U.x = vertices3[i2Max].x - vertices3[i1Max].x;
			U.y = vertices3[i2Max].y - vertices3[i1Max].y;
			U.z = vertices3[i2Max].z - vertices3[i1Max].z;
			*/

			U.x = p1.x - p0.x;
			U.y = p1.y - p0.y;
			U.z = p1.z - p0.z;
			
			double nU = Norme(&U);
			ScalarMult(&U,1.0/nU); // Normalize U

			// Construct a normal vector V
			Cross(&V,&(f->sh.N),&U); // |U|=1 and |N|=1 => |V|=1

			// u,v vertices (we start with p0 at 0,0)
			f->vertices2[0].u = 0.0;
			f->vertices2[0].v = 0.0;
			VERTEX2D min; min.u=0.0; min.v=0.0;
			VERTEX2D max; max.u=0.0; max.v=0.0;

			for(int j=1;j<f->sh.nbIndex;j++) {

				VERTEX3D p = vertices3[f->indices[j]];
				VERTEX3D v;
				Sub(&v,&p,&p0); // v = p0p
				f->vertices2[j].u = Dot(&U,&v);  // Project p on U along the V direction
				f->vertices2[j].v = Dot(&V,&v);  // Project p on V along the U direction

				// Bounds
				if( f->vertices2[j].u > max.u ) max.u = f->vertices2[j].u;
				if( f->vertices2[j].v > max.v ) max.v = f->vertices2[j].v;
				if( f->vertices2[j].u < min.u ) min.u = f->vertices2[j].u;
				if( f->vertices2[j].v < min.v ) min.v = f->vertices2[j].v;

			}

			// Calculate facet area (Meister/Gauss formula)
			double A=0.0;
			for(int j=0;j<f->sh.nbIndex;j++) {
				int j1 = IDX(j+1,f->sh.nbIndex);
				A += f->vertices2[j].u*f->vertices2[j1].v - f->vertices2[j1].u*f->vertices2[j].v;
			}
			f->sh.area = fabs(0.5 * A);

			// Compute the 2D basis (O,U,V)
			double uD = (max.u - min.u);
			double vD = (max.v - min.v);

			// Origin
			f->sh.O.x = min.u * U.x + min.v * V.x + p0.x;
			f->sh.O.y = min.u * U.y + min.v * V.y + p0.y;
			f->sh.O.z = min.u * U.z + min.v * V.z + p0.z;

			// Rescale U and V vector
			f->sh.nU = U;
			ScalarMult(&U,uD);
			f->sh.U = U;

			f->sh.nV = V;
			ScalarMult(&V,vD);
			f->sh.V = V;

			Cross(&(f->sh.Nuv),&U,&V);

			// Rescale u,v coordinates
			for(int j=0;j<f->sh.nbIndex;j++) {

				VERTEX2D p = f->vertices2[j];
				f->vertices2[j].u = (p.u - min.u) / uD;
				f->vertices2[j].v = (p.v - min.v) / vD;

			}

			// Detect non visible edge
			f->InitVisibleEdge();

			// Detect orientation
			f->DetectOrientation();

			// Hit address
			f->sh.hitOffset = fOffset;
			fOffset += f->GetHitsSize();

			// Texture ratio
			f->tRatio = f->sh.texWidthD / uD;
		}
	}

	/*
	//Update mesh on all facets
	if (facet_number==-1){
		for(int i=0;i<sh.nbFacet;i++) {
			if(facets[i]->hasMesh) {
				SetFacetTexture(i,facets[i]->tRatio,TRUE);
			}
		}
	}*/

	isLoaded = TRUE;
	if (facet_number==-1) {
		BuildGLList();
		mApp->UpdateModelParams();
		mApp->UpdateFacetParams();
	}
	//initGeoPrg->SetVisible(FALSE);
	//SAFE_DELETE(initGeoPrg);
}

// -----------------------------------------------------------
void Geometry::RebuildLists() {
	BuildGLList();
}

// -----------------------------------------------------------

DWORD Geometry::GetGeometrySize(std::vector<Region_full> *regions,std::vector<Material> *materials) {

	// Compute number of bytes allocated
	DWORD memoryUsage = 0;
	memoryUsage += sizeof(SHGEOM);
	memoryUsage += sh.nbVertex * sizeof(VERTEX3D);
	for(int i=0;i<sh.nbFacet;i++) {
		memoryUsage += facets[i]->GetGeometrySize();
	}
	//Regions
	memoryUsage += (int)(*regions).size()*sizeof(Region_mathonly);
	for (int i=0;i<(int)(*regions).size();i++) {
		memoryUsage += sizeof(Trajectory_Point)*(*regions)[i].Points.size();
		memoryUsage += 2*sizeof(double)*(*regions)[i].Bx_distr->size;
		memoryUsage += 2*sizeof(double)*(*regions)[i].By_distr->size;
		memoryUsage += 2*sizeof(double)*(*regions)[i].Bz_distr->size;
		memoryUsage += 6*sizeof(double)*(*regions)[i].nbDistr_BXY;
	}
	//Material library
	memoryUsage+=sizeof(int); //number of (*materials)
	for (int i=0;i<(int)(*materials).size();i++) { //going through all (*materials)
		memoryUsage+=sizeof(int);//copying number of angles (columns)
		memoryUsage+=sizeof(int);//copying number of energies (rows)
		memoryUsage+=(int)((*materials)[i].angleVals.size())*sizeof(double);//copying angles (header)
		memoryUsage+=(int)((*materials)[i].energyVals.size())*sizeof(double);//copying energies (column1)
		memoryUsage+=(int)((*materials)[i].angleVals.size())*(int)(*materials)[i].energyVals.size()*sizeof(double);//copying reflectivity probabilities (cells)
	}
	return memoryUsage;
}

// -----------------------------------------------------------

void Geometry::CopyGeometryBuffer(BYTE *buffer, std::vector<Region_full> *regions, std::vector<Material> *materials, int generation_mode, BOOL lowFluxMode, double lowFluxCutoff) {

	// Build shared buffer for geometry (see Shared.h)
	int fOffset = sizeof(SHGHITS);
	SHGEOM *shGeom = (SHGEOM *)buffer;
	sh.nbRegion=(*regions).size();
	sh.generation_mode=generation_mode;
	sh.lowFluxMode = lowFluxMode;
	sh.lowFluxCutoff = lowFluxCutoff;
	memcpy(shGeom,&(this->sh),sizeof(SHGEOM));
	buffer += sizeof(SHGEOM);

	// Build shared buffer for trajectory (see Shared.h)
	for (int i=0;i<sh.nbRegion;i++) {
		Region_mathonly *reg = (Region_mathonly *)buffer;
		(*regions)[i].nbDistr_MAG=Vector((*regions)[i].Bx_distr->size,(*regions)[i].By_distr->size,(*regions)[i].Bz_distr->size);
		(*regions)[i].nbPointsToCopy=(int)(*regions)[i].Points.size();
		*reg=(*regions)[i];
		buffer +=sizeof(Region_mathonly);
	}	
		//copy trajectory points
	for (int i=0;i<sh.nbRegion;i++) {
		for (int j=0;j<(*regions)[i].nbPointsToCopy;j++) {
			(*((Trajectory_Point*)(buffer)))=(*regions)[i].Points[j];
			buffer+=sizeof(Trajectory_Point);
		}
	}
		//copy distribution points
	for (int i=0;i<sh.nbRegion;i++) {
		for (int j=0;j<(*regions)[i].nbDistr_MAG.x;j++) {
			(*((double*)(buffer)))=(*regions)[i].Bx_distr->valuesX[j];
			buffer+=sizeof(double);
			(*((double*)(buffer)))=(*regions)[i].Bx_distr->valuesY[j];
			buffer+=sizeof(double);
		}

		for (int j=0;j<(*regions)[i].nbDistr_MAG.y;j++) {
			(*((double*)(buffer)))=(*regions)[i].By_distr->valuesX[j];
			buffer+=sizeof(double);
			(*((double*)(buffer)))=(*regions)[i].By_distr->valuesY[j];
			buffer+=sizeof(double);
		}

		for (int j=0;j<(*regions)[i].nbDistr_MAG.z;j++) {
			(*((double*)(buffer)))=(*regions)[i].Bz_distr->valuesX[j];
			buffer+=sizeof(double);
			(*((double*)(buffer)))=(*regions)[i].Bz_distr->valuesY[j];
			buffer+=sizeof(double);
		}

		for (int j=0;j<(*regions)[i].nbDistr_BXY;j++) {
			(*((double*)(buffer)))=(*regions)[i].beta_x_distr->valuesX[j];
			buffer+=sizeof(double);
			(*((double*)(buffer)))=(*regions)[i].beta_x_distr->valuesY[j];
			buffer+=sizeof(double);
			(*((double*)(buffer)))=(*regions)[i].beta_y_distr->valuesY[j];
			buffer+=sizeof(double);
			(*((double*)(buffer)))=(*regions)[i].eta_distr->valuesY[j];
			buffer+=sizeof(double);
			(*((double*)(buffer)))=(*regions)[i].etaprime_distr->valuesY[j];
			buffer+=sizeof(double);
			(*((double*)(buffer)))=(*regions)[i].e_spread_distr->valuesY[j];
			buffer+=sizeof(double);
		}
	}

	*((int*)buffer)=(int)(*materials).size(); //copying number of materials
			buffer+=sizeof(int);
			for (int i=0;i<(int)(*materials).size();i++) { //going through all materials
				*((int*)buffer)=(int)(*materials)[i].angleVals.size(); //copying number of angles (columns)
				buffer+=sizeof(int);
				*((int*)buffer)=(int)(*materials)[i].energyVals.size(); //copying number of energies (rows)
				buffer+=sizeof(int);
				for (int j=0;j<(int)(*materials)[i].angleVals.size();j++) {
					*((double*)buffer)=(*materials)[i].angleVals[j]; //copying angles (header)
					buffer+=sizeof(double);
				}
				for (int j=0;j<(int)(*materials)[i].energyVals.size();j++) {
					*((double*)buffer)=(*materials)[i].energyVals[j]; //copying energies (column1)
					buffer+=sizeof(double);
				}
				for (int j=0;j<(int)(*materials)[i].energyVals.size();j++) {
					for (int k=0;k<(int)(*materials)[i].angleVals.size();k++) {
						*((double*)buffer)=(*materials)[i].reflVals[j][k]; //copying reflectivity probabilities (cells)
						buffer+=sizeof(double);
					}
				}
			}

	memcpy(buffer,vertices3,sizeof(VERTEX3D)*sh.nbVertex);
	buffer += sizeof(VERTEX3D)*sh.nbVertex;
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		f->sh.hitOffset = fOffset;
		fOffset += f->GetHitsSize();
		memcpy(buffer,&(f->sh),sizeof(SHFACET));
		buffer += sizeof(SHFACET);
		memcpy(buffer,f->indices,sizeof(int)*f->sh.nbIndex);
		buffer += sizeof(int)*f->sh.nbIndex;
		memcpy(buffer,f->vertices2,sizeof(VERTEX2D)*f->sh.nbIndex);
		buffer += sizeof(VERTEX2D)*f->sh.nbIndex;
	}

	// Add surface elements area (reciprocal)
	for(int k=0;k<sh.nbFacet;k++) {
		Facet *f = facets[k];
		DWORD add = 0;
		if(f->sh.isTextured) {

			if( f->mesh ) {

				for(int j=0;j<f->sh.texHeight;j++) {
					for(int i=0;i<f->sh.texWidth;i++) {
						float area = f->mesh[add].area;
						if( area>0.0f ) {
							// Use the sign bit to store isFull flag
							if( f->mesh[add].full )
								*((double *)buffer) = -1.0f / area;
							else
								*((double *)buffer) = 1.0f / area;
						} else {
							*((double *)buffer) = 0.0f;
						}
						add++;
						buffer+=sizeof(double);
					}
				}

			} else {

				double rw = Norme(&(f->sh.U)) / (double)(f->sh.texWidthD);
				double rh = Norme(&(f->sh.V)) / (double)(f->sh.texHeightD);
				float area = (float)(rw*rh);

				for(int j=0;j<f->sh.texHeight;j++) {
					for(int i=0;i<f->sh.texWidth;i++) {
						if( area>0.0f ) {
							*((double *)buffer) = 1.0f / area;
						} else {
							*((double *)buffer) = 0.0f;
						}
						buffer+=sizeof(double);
					}
				}

			}
		}
	}

}

// -----------------------------------------------------------

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

// -----------------------------------------------------------

DWORD Geometry::GetHitsSize() {

	// Compute number of bytes allocated
	DWORD memoryUsage = 0;
	memoryUsage += sizeof(SHGHITS);
	for(int i=0;i<sh.nbFacet;i++) {
		memoryUsage += facets[i]->GetHitsSize();
	}

	return memoryUsage;
}

// -----------------------------------------------------------

DWORD Geometry::GetMaxElemNumber() {

	int nbElem = 0;
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		if( f->mesh ) nbElem += f->sh.texWidth*f->sh.texHeight;
		else          return 0;
	}
	return nbElem;

}

// -----------------------------------------------------------

void Geometry::CopyElemBuffer(BYTE *buffer) {

	int idx = 0;
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		int sz = f->sh.texWidth * f->sh.texHeight * sizeof(SHELEM);
		memcpy(buffer+idx,f->mesh,sz);
		idx+=sz;
	}

}

// -----------------------------------------------------------

void Geometry::BuildShapeList() {

	// Shapes used for direction field rendering

	// 3D arrow (direction field)
	int nbDiv = 10;
	double alpha = 2.0*PI/(double)nbDiv;

	arrowList = glGenLists(1);
	glNewList(arrowList,GL_COMPILE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	glBegin(GL_TRIANGLES);

	// Arrow
	for(int i=0;i<nbDiv;i++) {

		double y1 = sin(alpha*(double)i);
		double z1 = cos(alpha*(double)i);
		double y2 = sin(alpha*(double)((i+1)%nbDiv));
		double z2 = cos(alpha*(double)((i+1)%nbDiv));

		glNormal3d( 0.0 , y1 , z1 );
		glVertex3d(-0.5 , 0.5*y1 , 0.5*z1 );
		glNormal3d( 1.0 , 0.0 , 0.0 );
		glVertex3d( 0.5 , 0.0 , 0.0 );
		glNormal3d( 0.0 , y2 , z2 );
		glVertex3d(-0.5 , 0.5*y2 , 0.5*z2 );

	}

	// Cap facets
	for(int i=0;i<nbDiv;i++) {

		double y1 = sin(alpha*(double)i);
		double z1 = cos(alpha*(double)i);
		double y2 = sin(alpha*(double)((i+1)%nbDiv));
		double z2 = cos(alpha*(double)((i+1)%nbDiv));

		glNormal3d( -1.0 , 0.0 , 0.0 );
		glVertex3d( -0.5 , 0.5*y1 , 0.5*z1 );
		glNormal3d( -1.0 , 0.0 , 0.0 );
		glVertex3d( -0.5 , 0.5*y2 , 0.5*z2 );
		glNormal3d( -1.0 , 0.0 , 0.0 );
		glVertex3d( -0.5 , 0.0 , 0.0 );

	}

	glEnd();
	glEndList();

	// Shpere list (isotropic case)
	int nbPhi   = 16;
	int nbTetha = 7;
	double dphi   = 2.0*PI / (double)(nbPhi);
	double dtetha = PI / (double)(nbTetha+1);

	sphereList = glGenLists(1);
	glNewList(sphereList,GL_COMPILE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glDisable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	glBegin(GL_TRIANGLES);

	for(int i=0;i<=nbTetha;i++) {
		for(int j=0;j<nbPhi;j++) {

			VERTEX3D v1,v2,v3,v4;

			v1.x = sin(dtetha*(double)i)*cos(dphi*(double)j);
			v1.y = sin(dtetha*(double)i)*sin(dphi*(double)j);
			v1.z = cos(dtetha*(double)i);

			v2.x = sin(dtetha*(double)(i+1))*cos(dphi*(double)j);
			v2.y = sin(dtetha*(double)(i+1))*sin(dphi*(double)j);
			v2.z = cos(dtetha*(double)(i+1));

			v3.x = sin(dtetha*(double)(i+1))*cos(dphi*(double)(j+1));
			v3.y = sin(dtetha*(double)(i+1))*sin(dphi*(double)(j+1));
			v3.z = cos(dtetha*(double)(i+1));

			v4.x = sin(dtetha*(double)i)*cos(dphi*(double)(j+1));
			v4.y = sin(dtetha*(double)i)*sin(dphi*(double)(j+1));
			v4.z = cos(dtetha*(double)i);

			if( i<nbTetha ) {
				glNormal3d(v1.x , v1.y , v1.z);
				glVertex3d(v1.x , v1.y , v1.z);
				glNormal3d(v2.x , v2.y , v2.z);
				glVertex3d(v2.x , v2.y , v2.z);
				glNormal3d(v3.x , v3.y , v3.z);
				glVertex3d(v3.x , v3.y , v3.z);
			}

			if( i>0 ) {
				glNormal3d(v1.x , v1.y , v1.z);
				glVertex3d(v1.x , v1.y , v1.z);
				glNormal3d(v3.x , v3.y , v3.z);
				glVertex3d(v3.x , v3.y , v3.z);
				glNormal3d(v4.x , v4.y , v4.z);
				glVertex3d(v4.x , v4.y , v4.z);
			}

		}
	}

	glEnd();
	glEndList();


}

// -----------------------------------------------------------

void Geometry::BuildSelectList() {

	nbSelected = 0;

	selectList = glGenLists(1);
	glNewList(selectList,GL_COMPILE);
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
	glNewList(selectList2,GL_COMPILE);
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
	glNewList(selectList3,GL_COMPILE);

	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	if (mApp->antiAliasing){
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	}
	glLineWidth(2.0f);

	for(int i=0;i<sh.nbFacet;i++ ) {
		Facet *f = facets[i];
		if( f->selected ) {
			//DrawFacet(f,FALSE,TRUE,TRUE);
			DrawFacet(f,1,1,1);
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



// -----------------------------------------------------------

void Geometry::UpdateSelection() {

	DeleteGLLists();
	BuildSelectList();

}

// -----------------------------------------------------------

void Geometry::BuildGLList() {

	// Compile geometry for OpenGL
	for(int j=0;j<sh.nbSuper;j++) {
		lineList[j] = glGenLists(1);
		glNewList(lineList[j],GL_COMPILE);
		for(int i=0;i<sh.nbFacet;i++) {
			if(facets[i]->sh.superIdx==j)
				DrawFacet(facets[i]);
		}
		glEndList();
	}

	polyList = glGenLists(1);
	glNewList(polyList,GL_COMPILE);
	DrawPolys();
	glEndList();

	BuildSelectList();

}


// -----------------------------------------------------------
// Collapse stuff
// -----------------------------------------------------------


// -----------------------------------------------------------

Facet *Geometry::MergeFacet(Facet *f1,Facet *f2) {
	mApp->changedSinceSave = TRUE;
	// Merge 2 facets into 1 when possible and create a new facet
	// otherwise return NULL.
	int  c1;
	int  c2;
	int  l;
	BOOL end = FALSE;
	Facet *nF = NULL;

	if( GetCommonEdges(f1,f2,&c1,&c2,&l) ) {
		int commonNo = f1->sh.nbIndex + f2->sh.nbIndex - 2 * l;
		if (commonNo == 0) { //two identical facets, so return a copy of f1
			nF = new Facet(f1->sh.nbIndex);
			nF->Copy(f1);
			for (int i = 0; i < f1->sh.nbIndex; i++)
				nF->indices[i] = f1->GetIndex(i);
			return nF;
		}

		int nbI=0;
		nF = new Facet(commonNo);
		// Copy params from f1
		//nF->Copy(f1);
		nF->Copy(f1);

		if( l==f1->sh.nbIndex ) {

			// f1 absorbed, copy indices from f2
			for(int i=0;i<f2->sh.nbIndex - l;i++)
				nF->indices[nbI++] = f2->GetIndex(c2+2+i);

		} else if ( l==f2->sh.nbIndex ) {

			// f2 absorbed, copy indices from f1
			for(int i=0;i<f1->sh.nbIndex - l;i++)
				nF->indices[nbI++] = f1->GetIndex(c1+l+i);

		} else {    

			// Copy indices from f1
			for(int i=0;i<f1->sh.nbIndex - (l-1);i++)
				nF->indices[nbI++] = f1->GetIndex(c1+l+i);
			// Copy indices from f2
			for(int i=0;i<f2->sh.nbIndex - (l+1);i++)
				nF->indices[nbI++] = f2->GetIndex(c2+2+i);

		}

	}

	return nF;

}

// -----------------------------------------------------------
void Geometry::Collapse(double vT,double fT,double lT,BOOL doSelectedOnly,GLProgress *prg) {	
	mApp->changedSinceSave = TRUE;
	Facet *fi,*fj;
	Facet *merged;

	vThreshold = vT;
	double totalWork=(1.0 + (double)(fT>0.0) + (double)(lT>0.0)); //for progress indicator
	// Collapse vertex
	if (vT>0.0) {
		CollapseVertex(prg,totalWork);
		RemoveCollinear();
		RemoveNullFacet();
		InitializeGeometry();
	}


	if( fT>0.0 ) {

		// Collapse facets
		int i = 0;
		prg->SetMessage("Collapsing facets...");
		while(i<sh.nbFacet) {
			prg->SetProgress((1.0+((double)i/(double)sh.nbFacet)) / totalWork);
			fi = facets[i];
			// Search a coplanar facet
			int j = i+1;
			while((!doSelectedOnly || fi->selected) && j<sh.nbFacet) {
				fj = facets[j];
				merged = NULL;
				if( (!doSelectedOnly || fj->selected) && fi->IsCoplanar(fj,fT) ) {
					// Collapse
					merged = MergeFacet(fi,fj);
					if( merged ) {
						// Replace the old 2 facets by the new one
						SAFE_DELETE(fi);
						SAFE_DELETE(fj);
						for(int k=j;k<sh.nbFacet-1;k++)
							facets[k] = facets[k+1];
						sh.nbFacet--;
						facets[i] = merged;
						//InitializeGeometry(i);
						//SetFacetTexture(i,facets[i]->tRatio,facets[i]->hasMesh);  //rebuild mesh
						fi = facets[i];
						mApp->RenumberSelections(j);
						mApp->RenumberFormulas(j);
						j=i+1;

					}
				}
				if(!merged) j++;
			}
			i++;
		}
	}
	//Collapse collinear sides. Takes some time, so only if threshold>0
	prg->SetMessage("Collapsing collinear sides...");
	if(lT>0.0) {
		for(int i=0;i<sh.nbFacet;i++) {
			prg->SetProgress((1.0+(double)(fT>0.0)+((double)i/(double)sh.nbFacet)) / totalWork);
			if (!doSelectedOnly || facets[i]->selected)
				MergecollinearSides(facets[i],lT);
		}
	}
	prg->SetMessage("Rebuilding geometry...");
	for(int i=0;i<sh.nbFacet;i++) {

		Facet *f = facets[i];

		// Revert orientation if normal has been swapped
		// This happens when the second vertex is no longer convex
		VERTEX3D n,v1,v2;
		double   d;
		int i0 = facets[i]->indices[0];
		int i1 = facets[i]->indices[1];
		int i2 = facets[i]->indices[2];

		Sub(&v1 , vertices3+i1 , vertices3+i0); // v1 = P0P1
		Sub(&v2 , vertices3+i2 , vertices3+i1); // v2 = P1P2
		Cross(&n,&v1,&v2);                      // Cross product
		d = Dot(&n,&(f->sh.N));
		if( d<0.0 ) f->SwapNormal();

	}



	// Delete old resources
	for(int i=0;i<sh.nbSuper;i++)
		DeleteGLLists(TRUE,TRUE);

	// Reinitialise geom
	InitializeGeometry();

}

void Geometry::MergecollinearSides(Facet *f,double lT) {
	mApp->changedSinceSave = TRUE;
	BOOL collinear;
	double linTreshold=cos(lT*PI/180);
	// Merge collinear sides
	for(int k=0;(k<f->sh.nbIndex&&f->sh.nbIndex>3);k++){
		k=k;
		do {
			//collinear=FALSE;
			int p0=f->indices[k];
			int p1=f->indices[(k+1)%f->sh.nbIndex];
			int p2=f->indices[(k+2)%f->sh.nbIndex]; //to compare last side with first too
			VERTEX3D p0p1;
			VERTEX3D p0p2;
			Sub(&p0p1,&vertices3[p1],&vertices3[p0]);
			Sub(&p0p2,&vertices3[p2],&vertices3[p1]);
			Normalize(&p0p1);
			Normalize(&p0p2);
			collinear=(Dot(&p0p1,&p0p2)>=linTreshold);
			if (collinear&&f->sh.nbIndex>3) { //collinear
				for (int l=(k+1)%f->sh.nbIndex;l<f->sh.nbIndex-1;l++){
					f->indices[l]=f->indices[l+1];
				}
				f->sh.nbIndex--;
			}
		} while (collinear&&f->sh.nbIndex>3);
	}
}

// -----------------------------------------------------------

void Geometry::Rebuild() {

	// Rebuild internal structure on geometry change

	// Remove texture (improvement TODO)
	for(int i=0;i<sh.nbFacet;i++)
		facets[i]->SetTexture(0.0,0.0,FALSE);
	// Delete old resources
	DeleteGLLists(TRUE,TRUE);

	// Reinitialise geom
	InitializeGeometry();

}

int Geometry::InvalidateDeviceObjects() {

	DeleteGLLists(TRUE,TRUE);
	DELETE_LIST(arrowList);
	DELETE_LIST(sphereList);
	for(int i=0;i<sh.nbFacet;i++)
		facets[i]->InvalidateDeviceObjects();

	return GL_OK;

}

int Geometry::RestoreDeviceObjects() {

	if(!IsLoaded()) return GL_OK;

	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		f->RestoreDeviceObjects();
		BuildFacetList(f);
	}

	BuildGLList();

	return GL_OK;

}

// -----------------------------------------------------------

void Geometry::BuildFacetList(Facet *f) {

	// Rebuild OpenGL geomtetry with texture
	//if( f->sh.isTextured && f->textureVisible ) {
	if( f->sh.isTextured ) {
		// Facet geometry
		glNewList(f->glList,GL_COMPILE);
		if( f->sh.nbIndex == 3 ) {
			glBegin(GL_TRIANGLES);
			FillFacet(f,TRUE);
			glEnd();
		} else if( f->sh.nbIndex == 4 ) {
			glBegin(GL_QUADS);
			FillFacet(f,TRUE);
			glEnd();
		} else {
			glBegin(GL_TRIANGLES);
			Triangulate(f,TRUE);
			glEnd();
		}
		glEndList();

	}

}

// -----------------------------------------------------------

void Geometry::SetFacetTexture(int facet,double ratio,BOOL mesh) {

	Facet *f = facets[facet];
	double nU = Norme(&(f->sh.U));
	double nV = Norme(&(f->sh.V));
	
	if (!f->SetTexture(nU*ratio, nV*ratio, mesh)) {
		char errMsg[512];
		sprintf(errMsg, "Not enough memory to build mesh on Facet %d. ", facet + 1);
		throw Error(errMsg);
	}
	f->tRatio = ratio;
	BuildFacetList(f);

}

// -----------------------------------------------------------
// Testing purpose function, construct a PIPE
// -----------------------------------------------------------
void  Geometry::BuildPipe(double L,double R,double s,int step) {
	Clear();

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	sprintf(sh.name,"PIPE%g",L/R);

	int nbDecade = 0;
	int nbTF = 9*nbDecade;
	int nbTV = 4*nbTF;

	sh.nbVertex  = 2*step + nbTV;
	if (!(vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D))))
		throw Error("Couldn't allocate memory for vertices");
	memset(vertices3,0,sh.nbVertex * sizeof(VERTEX3D));

	sh.nbFacet   = step + 2 + nbTF;
	sh.nbSuper = 1;
	if (!(facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *))))
		throw Error("Couldn't allocate memory for facets");
	memset(facets,0,sh.nbFacet * sizeof(Facet *));

	// Vertices
	for(int i=0;i<step;i++) {
		double angle = (double)i/(double)step * 2*PI;
		vertices3[2*i+nbTV].x = R*cos(angle);
		vertices3[2*i+nbTV].y = R*sin(angle);
		vertices3[2*i+nbTV].z = 0.0;
		vertices3[2*i+1+nbTV].x = R*cos(angle);
		vertices3[2*i+1+nbTV].y = R*sin(angle);
		vertices3[2*i+1+nbTV].z = L;
	}

    try{
	// Cap facet
	facets[0+nbTF] = new Facet(step);
	facets[0+nbTF]->sh.sticking = 1.0;
	for(int i=0;i<step;i++)
		facets[0+nbTF]->indices[i] = 2*i+nbTV;

	facets[1+nbTF] = new Facet(step);
	facets[1+nbTF]->sh.sticking = 1.0;
	for(int i=0;i<step;i++)
		facets[1+nbTF]->indices[step-i-1] = 2*i+1+nbTV;

	// Wall facet
	for(int i=0;i<step;i++) {
		facets[i+2+nbTF] = new Facet(4);
		facets[i+2+nbTF]->sh.reflectType = REF_DIFFUSE;
		facets[i+2+nbTF]->sh.sticking = s;
		facets[i+2+nbTF]->indices[0]  = 2*i+nbTV;
		facets[i+2+nbTF]->indices[1]  = 2*i+1+nbTV;
		if( i<step-1 ) {
			facets[i+2+nbTF]->indices[2] = 2*(i+1)+1+nbTV;
			facets[i+2+nbTF]->indices[3] = 2*(i+1)+nbTV;
		} else {
			facets[i+2+nbTF]->indices[2] = 1+nbTV;
			facets[i+2+nbTF]->indices[3] = 0+nbTV;
		}
	}

	// Volatile facet
	for(int d=0;d<nbDecade;d++) {
		for(int i=0;i<9;i++) {

			double z = (double)(i+1) * pow(10,(double)d);
			int idx = d*36+i*4;

			vertices3[idx+0].x = -R;
			vertices3[idx+0].y =  R;
			vertices3[idx+0].z =  z;
			vertices3[idx+1].x =  R;
			vertices3[idx+1].y =  R;
			vertices3[idx+1].z =  z;
			vertices3[idx+2].x =  R;
			vertices3[idx+2].y = -R;
			vertices3[idx+2].z =  z;
			vertices3[idx+3].x = -R;
			vertices3[idx+3].y = -R;
			vertices3[idx+3].z =  z;

			facets[9*d+i] = new Facet(4);
			facets[9*d+i]->sh.sticking = 0.0;
			facets[9*d+i]->sh.opacity = 0.0;
			facets[9*d+i]->sh.isVolatile = TRUE;
			facets[9*d+i]->indices[0] = idx+0;
			facets[9*d+i]->indices[1] = idx+1;
			facets[9*d+i]->indices[2] = idx+2;
			facets[9*d+i]->indices[3] = idx+3;

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

// -----------------------------------------------------------
// File handling
// -----------------------------------------------------------

void Geometry::UpdateName(FileReader *file) {

	char *p = strrchr(file->GetName(),'\\');
	if(!p) p = strrchr(file->GetName(),'/');

	if(!p)
		strcpy(sh.name,file->GetName());
	else
		strcpy(sh.name,p+1);

}



void Geometry::AdjustProfile() {

	// Backward compatibily with TXT profile (To be improved)
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		if( f->sh.profileType==REC_PRESSUREU ) {
			VERTEX3D v0;
			Sub(&v0 , vertices3+(f->indices[1]) , vertices3+(f->indices[0])); // v0 = P0P1
			double n0 = Norme(&v0);
			double nU = Norme(&(f->sh.U));
			if(IS_ZERO(n0-nU)) f->sh.profileType = REC_PRESSUREU; // Select U
			else               f->sh.profileType = REC_PRESSUREV; // Select V
		}
	}

}



void Geometry::LoadASE(FileReader *file,GLProgress *prg) {

	Clear();

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	ASELoader ase(file);
	ase.Load();

	// Compute total of facet
	sh.nbFacet = 0;
	for(int i=0;i<ase.nbObj;i++) sh.nbFacet += ase.OBJ[i].nb_face;

	// Allocate mem
	sh.nbVertex = 3*sh.nbFacet;
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets,0,sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memset(vertices3,0,sh.nbVertex * sizeof(VERTEX3D));

	// Fill 
	int nb = 0;
	for(int i=0;i<ase.nbObj;i++) {

		for(int j=0;j<ase.OBJ[i].nb_face;j++ ) {
			vertices3[3*nb+0] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v1];
			vertices3[3*nb+1] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v2];
			vertices3[3*nb+2] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v3];
			facets[nb] = new Facet(3);
			facets[nb]->indices[0]=3*nb+0;
			facets[nb]->indices[1]=3*nb+1;
			facets[nb]->indices[2]=3*nb+2;
			nb++;
		}

	}

	UpdateName(file);
	sh.nbSuper = 1;
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	isLoaded = TRUE;

}



void Geometry::LoadSTR(FileReader *file,GLProgress *prg) {

	char nPath[512];
	char fPath[512];
	char fName[512];
	char sName[512];
	int nF,nV;
	Facet **F;
	VERTEX3D *V;
	FileReader *fr;

	Clear();

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	// Load multiple structure file
	sh.nbSuper = file->ReadInt();

	strcpy(fPath,file->ReadLine());
	strcpy(nPath,FileUtils::GetPath(file->GetName()));

	for(int n=0;n<sh.nbSuper;n++) {

		int i1 = file->ReadInt();
		int i2 = file->ReadInt();
		fr = NULL;
		strcpy(sName,file->ReadWord());
		strName[n] = _strdup(sName);
		char *e = strrchr(strName[n],'.');
		if(e) *e=0;

		sprintf(fName,"%s\\%s",nPath,sName);
		if( FileUtils::Exist(fName) ) {
			fr = new FileReader(fName);
			strcpy(strPath,nPath);
		} else {
			sprintf(fName,"%s\\%s",fPath,sName);
			if( FileUtils::Exist(fName) ) {
				fr = new FileReader(fName);
				strcpy(strPath,fPath);
			}
		}

		if( !fr ) {
			char errMsg[512];
			sprintf(errMsg,"Cannot find %s",sName);
			throw Error(errMsg);
		}

		strFileName[n] = _strdup(sName);
		LoadTXTGeom(fr,&nV,&nF,&V,&F,n);
		Merge(nV,nF,V,F);
		SAFE_FREE(V);
		SAFE_FREE(F);
		delete fr;

	}

	UpdateName(file);
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::LoadSTL(FileReader *file,GLProgress *prg,double scaleFactor) {

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
	while(strcmp(w,"facet")==0) {
		sh.nbFacet++;
		file->JumpSection("endfacet");
		w = file->ReadWord();
	}
	if( strcmp(w,"endsolid")!=0 ) throw Error("Unexpected or not supported STL keyword, 'endsolid' required");

	// Allocate mem
	sh.nbVertex = 3*sh.nbFacet;
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	if (!facets) throw Error("Out of memory: LoadSTL");
	memset(facets,0,sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	if (!vertices3) throw Error("Out of memory: LoadSTL");
	memset(vertices3,0,sh.nbVertex * sizeof(VERTEX3D));

	// Second pass
	prg->SetMessage("Reading facets...");
	file->SeekStart();
	//file->ReadKeyword("solid");
	file->ReadLine();
	for(int i=0;i<sh.nbFacet;i++) {

		double p = (double)i/(double)(sh.nbFacet);
		prg->SetProgress(p);

		file->ReadKeyword("facet");
		file->ReadKeyword("normal");
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
		file->ReadKeyword("outer");
		file->ReadKeyword("loop");

		file->ReadKeyword("vertex");
		vertices3[3*i+0].x = file->ReadDouble()*scaleFactor;
		vertices3[3*i+0].y = file->ReadDouble()*scaleFactor;
		vertices3[3*i+0].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		vertices3[3*i+1].x = file->ReadDouble()*scaleFactor;
		vertices3[3*i+1].y = file->ReadDouble()*scaleFactor;
		vertices3[3*i+1].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		vertices3[3*i+2].x = file->ReadDouble()*scaleFactor;
		vertices3[3*i+2].y = file->ReadDouble()*scaleFactor;
		vertices3[3*i+2].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("endloop");
		file->ReadKeyword("endfacet");

		try {
			facets[i] = new Facet(3);
		}
		catch (std::bad_alloc& badalloc) {

			throw Error("Out of memory");

		}
		facets[i]->indices[0]=3*i+0;
		facets[i]->indices[1]=3*i+2;
		facets[i]->indices[2]=3*i+1;

	}

	sh.nbSuper = 1;
	UpdateName(file);
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	prg->SetMessage("Initializing geometry...");
	InitializeGeometry();
	isLoaded = TRUE;

}
void Geometry::LoadTXT(FileReader *file,GLProgress *prg) {

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	Clear();
	LoadTXTGeom(file,&(sh.nbVertex),&(sh.nbFacet),&vertices3,&facets);
	UpdateName(file);
	sh.nbSuper = 1;
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;
}



void Geometry::InsertTXT(FileReader *file,GLProgress *prg,BOOL newStr) {

	//Clear();
	int structId=viewStruct;
	if (structId==-1) structId=0;
	InsertTXTGeom(file,&(sh.nbVertex),&(sh.nbFacet),&vertices3,&facets,structId,newStr);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::InsertSTL(FileReader *file,GLProgress *prg,double scaleFactor,BOOL newStr) {

	//Clear();
	int structId=viewStruct;
	if (structId==-1) structId=0;
	InsertSTLGeom(file,&(sh.nbVertex),&(sh.nbFacet),&vertices3,&facets,structId,newStr,scaleFactor);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;

}



void Geometry::InsertGEO(FileReader *file,GLProgress *prg,BOOL newStr) {

	//Clear();
	int structId=viewStruct;
	if (structId==-1) structId=0;
	InsertGEOGeom(file,&(sh.nbVertex),&(sh.nbFacet),&vertices3,&facets,structId,newStr);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;

}

PARfileList Geometry::InsertSYN(FileReader *file,GLProgress *prg,BOOL newStr) {

	int structId=viewStruct;
	if (structId==-1) structId=0;
	PARfileList result=InsertSYNGeom(file,&(sh.nbVertex),&(sh.nbFacet),&vertices3,&facets,structId,newStr);
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	//AdjustProfile();
	return result;
}



void Geometry::LoadTXTGeom(FileReader *file,int *nbV,int *nbF,VERTEX3D **V,Facet ***F,int strIdx) {

	file->ReadInt(); // Unused
	tNbHit = file->ReadLLong();
	tNbLeak = file->ReadLLong();
	tNbDesorption = file->ReadLLong();
	tNbDesorptionMax = file->ReadLLong(); 

	int nV = file->ReadInt();
	int nF = file->ReadInt();

	// Allocate memory
	Facet   **f = (Facet **)malloc(nF * sizeof(Facet *));
	memset(f,0,nF * sizeof(Facet *));
	VERTEX3D *v = (VERTEX3D *)malloc(nV * sizeof(VERTEX3D));
	memset(v,0,nV * sizeof(VERTEX3D)); //avoid selected flag

	// Read geometry vertices
	for(int i=0;i<nV;i++) {
		v[i].x = file->ReadDouble();
		v[i].y = file->ReadDouble();
		v[i].z = file->ReadDouble();
	}

	// Read geometry facets (indexed from 1)
	for(int i=0;i<nF;i++) {
		int nb = file->ReadInt();
		f[i] = new Facet(nb);
		for(int j=0;j<nb;j++)
			f[i]->indices[j]=file->ReadInt()-1;
	}

	// Read facets params
	for(int i=0;i<nF;i++) {
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

// -----------------------------------------------------------
// -----------------------------------------------------------

void Geometry::InsertTXTGeom(FileReader *file,int *nbVertex,int *nbFacet,VERTEX3D **vertices3,Facet ***facets,int strIdx,BOOL newStruct) {

	UnSelectAll();

	//tNbHit = file->ReadLLong();
	//tNbLeak = file->ReadInt();
	//tNbDesorption = file->ReadLLong();
	//tNbDesorptionMax = file->ReadLLong(); 
	for (int i=0;i<5;i++) file->ReadInt(); //leading lines

	int nbNewVertex = file->ReadInt();
	int nbNewFacets = file->ReadInt();

	// Allocate memory
	*facets = (Facet **)realloc(*facets,(nbNewFacets+*nbFacet) * sizeof(Facet **));
	memset(*facets+*nbFacet,0,nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3,*vertices3,(*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3+*nbVertex,0,nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3=tmp_vertices3;

	// Read geometry vertices
	for(int i=*nbVertex;i<(*nbVertex+nbNewVertex);i++) {
		(*vertices3+i)->x = file->ReadDouble();
		(*vertices3+i)->y = file->ReadDouble();
		(*vertices3+i)->z = file->ReadDouble();
		(*vertices3+i)->selected = FALSE;
	}

	// Read geometry facets (indexed from 1)
	for(int i=*nbFacet;i<(*nbFacet+nbNewFacets);i++) {
		int nb = file->ReadInt();
		*(*facets+i) = new Facet(nb);
		(*facets)[i]->selected = TRUE;
		for(int j=0;j<nb;j++)
			(*facets)[i]->indices[j]=file->ReadInt()-1+*nbVertex;
	}

	// Read facets params
	for(int i=*nbFacet;i<(*nbFacet+nbNewFacets);i++) {
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


void Geometry::InsertGEOGeom(FileReader *file,int *nbVertex,int *nbFacet,VERTEX3D **vertices3,Facet ***facets,int strIdx,BOOL newStruct) {

	UnSelectAll();
	//char tmp[512];

	file->ReadKeyword("version");file->ReadKeyword(":");
	int version2;
	version2 = file->ReadInt();
	if (version2 > GEOVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported GEO version V%d", version2);
		throw Error(errMsg);
	}
	
	file->ReadKeyword("totalHit");file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalDes");file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalLeak");file->ReadKeyword(":");
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
	file->ReadKeyword("maxDes");file->ReadKeyword(":");
	file->ReadLLong(); 
	file->ReadKeyword("nbVertex");file->ReadKeyword(":");
	int nbNewVertex = file->ReadInt();
	file->ReadKeyword("nbFacet");file->ReadKeyword(":");
	int nbNewFacets = file->ReadInt();
	file->ReadKeyword("nbSuper");file->ReadKeyword(":");
	int nbNewSuper = file->ReadInt();
	int nbF=0;
	int nbV=0;
	if (version2>=2) {
		file->ReadKeyword("nbFormula");file->ReadKeyword(":");
		nbF = file->ReadInt();
		file->ReadKeyword("nbView");file->ReadKeyword(":");
		nbV = file->ReadInt();
	}
	int nbS=0;
	if (version2>=8) {
		file->ReadKeyword("nbSelection");file->ReadKeyword(":");
		nbS = file->ReadInt();
	}
	if (version2>=7) {
		file->ReadKeyword("gasMass");file->ReadKeyword(":");
		file->ReadDouble();
	}
	if (version2>=10) { //time-dependent version
		file->ReadKeyword("userMoments");file->ReadKeyword("{");
		file->ReadKeyword("nb");file->ReadKeyword(":");
		int nb=file->ReadInt();

		for(int i=0;i<nb;i++) 
			file->ReadString();
		file->ReadKeyword("}");
	}
	if (version2>=11) { //gas pulse parameters
		file->ReadKeyword("desorptionStart");file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("desorptionStop");file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("timeWindow");file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("useMaxwellian");file->ReadKeyword(":");
		file->ReadInt();
	}

	if (version2>=12) { //2013.aug.22
		file->ReadKeyword("calcConstantFlow");file->ReadKeyword(":");
		file->ReadInt();
	}
	if (version2>=2) {
		file->ReadKeyword("formulas");file->ReadKeyword("{");
		for(int i=0;i<nbF;i++) {
			char tmpName[256];
			char tmpExpr[512];
			strcpy(tmpName,file->ReadString());
			strcpy(tmpExpr,file->ReadString());
			mApp->OffsetFormula(tmpExpr, sh.nbFacet);
			mApp->AddFormula(tmpName,tmpExpr);
		}
		file->ReadKeyword("}");

		file->ReadKeyword("views");file->ReadKeyword("{");
		for(int i=0;i<nbV;i++) {
			char tmpName[256];
			AVIEW v;
			strcpy(tmpName,file->ReadString());
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
			mApp->AddView(tmpName,v);
		}
		file->ReadKeyword("}");
	}

	if (version2>=8) {
		file->ReadKeyword("selections");file->ReadKeyword("{");
		for(int i=0;i<nbS;i++) {
			ASELECTION s;
			char tmpName[256];
			strcpy(tmpName,file->ReadString());
			s.name=_strdup(tmpName);
			s.nbSel = file->ReadInt();
			s.selection = (int *)malloc((s.nbSel)*sizeof(int));

			for (int j = 0; j < s.nbSel; j++) {
				s.selection[j] = file->ReadInt() + sh.nbFacet; //offset facet number by current number of facets
			}
			mApp->AddSelection(s.name,s);
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("structures");file->ReadKeyword("{");
	for(int i=0;i<nbNewSuper;i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		/* //Commented out for GEO
		sprintf(tmp,"%s.txt",strName[i]);
		strFileName[i] = _strdup(tmp);
		*/
	}
	file->ReadKeyword("}");

	// Reallocate memory
	*facets = (Facet **)realloc(*facets,(nbNewFacets+*nbFacet) * sizeof(Facet **));
	memset(*facets+*nbFacet,0,nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	if (!tmp_vertices3) throw Error("Out of memory: InsertGEOGeom");
	memmove(tmp_vertices3,*vertices3,(*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3+*nbVertex,0,nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3=tmp_vertices3;

	// Read geometry vertices
	file->ReadKeyword("vertices");file->ReadKeyword("{");
	for(int i=*nbVertex;i<(*nbVertex+nbNewVertex);i++) {
		// Check idx
		int idx = file->ReadInt();
		if( idx != i-*nbVertex+1 ) throw Error(file->MakeError("Wrong vertex index !"));
		(*vertices3+i)->x = file->ReadDouble();
		(*vertices3+i)->y = file->ReadDouble();
		(*vertices3+i)->z = file->ReadDouble();
		(*vertices3+i)->selected = FALSE;
	}
	file->ReadKeyword("}");

	if (version2>=6) {
		// Read leaks
		file->ReadKeyword("leaks");file->ReadKeyword("{");
		file->ReadKeyword("nbLeak");file->ReadKeyword(":");
		int nbleak2 = file->ReadInt();
		for (int i=0;i<nbleak2;i++) {
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
		file->ReadKeyword("hits");file->ReadKeyword("{");
		file->ReadKeyword("nbHHit");file->ReadKeyword(":");
		int nbHHit2 = file->ReadInt();
		for (int i=0;i<nbHHit2;i++) {
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
	for(int i=*nbFacet;i<(*nbFacet+nbNewFacets);i++) {
		file->ReadKeyword("facet");  
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1-*nbFacet ) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nb = file->ReadInt();

		if( nb<3 ) {
			char errMsg[512];
			sprintf(errMsg,"Facet %d has only %d vertices. ",i,nb);
			throw Error(errMsg);
		}

		*(*facets+i) = new Facet(nb);
		(*facets)[i]->LoadGEO(file,version2,nbNewVertex);
		(*facets)[i]->selected = TRUE;
		for(int j=0;j<nb;j++)
			(*facets)[i]->indices[j]+=*nbVertex;
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
	if (newStruct) sh.nbSuper+=nbNewSuper;
	else if (sh.nbSuper < strIdx + nbNewSuper) sh.nbSuper = strIdx + nbNewSuper;

}

PARfileList Geometry::InsertSYNGeom(FileReader *file,int *nbVertex,int *nbFacet,VERTEX3D **vertices3,Facet ***facets,int strIdx,BOOL newStruct) {

	PARfileList result(0);
	UnSelectAll();
	//char tmp[512];

	file->ReadKeyword("version");file->ReadKeyword(":");
	int version2;
	version2 = file->ReadInt();
		if( version2>SYNVERSION ) {
		char errMsg[512];
		sprintf(errMsg,"Unsupported SYN version V%d",version2);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit");file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalDes");file->ReadKeyword(":");
	file->ReadLLong();
	if (version2>=6) {
		file->ReadKeyword("no_scans");file->ReadKeyword(":");
		loaded_no_scans = file->ReadDouble();
	}
	else loaded_no_scans = 0;
	file->ReadKeyword("totalLeak");file->ReadKeyword(":");
	file->ReadLLong();
	if (version2>2) {
		file->ReadKeyword("totalFlux");file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("totalPower");file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("maxDes");file->ReadKeyword(":");
	file->ReadLLong(); 
	file->ReadKeyword("nbVertex");file->ReadKeyword(":");
	int nbNewVertex = file->ReadInt();
	file->ReadKeyword("nbFacet");file->ReadKeyword(":");
	int nbNewFacets = file->ReadInt();
	file->ReadKeyword("nbSuper");file->ReadKeyword(":");
	int nbNewSuper = file->ReadInt();
	int nbF = 0;  std::vector<std::vector<string>> loadFormulas;
	int nbV = 0; 
	file->ReadKeyword("nbFormula");file->ReadKeyword(":");
	nbF = file->ReadInt();
	file->ReadKeyword("nbView");file->ReadKeyword(":");
	nbV = file->ReadInt();
	int nbS=0;
	
	file->ReadKeyword("nbSelection");file->ReadKeyword(":");
	nbS = file->ReadInt();
	if (version2>1) {
		file->ReadKeyword("nbRegions");file->ReadKeyword(":");
		int nbR = file->ReadInt();
		result=PARfileList(nbR);

		file->ReadKeyword("PARfiles");file->ReadKeyword("{");
		for (int i=0;i<nbR;i++) {
			char tmp[512];
			strcpy(tmp,file->ReadString());
			result.fileNames[i]=_strdup(tmp);
		}
		file->ReadKeyword("}");
	}
	
	file->ReadKeyword("formulas");file->ReadKeyword("{");
	for(int i=0;i<nbF;i++) {
		char tmpName[256];
		char tmpExpr[512];
		strcpy(tmpName,file->ReadString());
		strcpy(tmpExpr,file->ReadString());
		//mApp->AddFormula(tmpName, tmpExpr); //parse after selection groups are loaded
		std::vector<string> newFormula;
		newFormula.push_back(tmpName);
		mApp->OffsetFormula(tmpExpr, sh.nbFacet); //offset formula
		newFormula.push_back(tmpExpr);
		loadFormulas.push_back(newFormula);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("views");file->ReadKeyword("{");
	for(int i=0;i<nbV;i++) {
		char tmpName[256];
		AVIEW v;
		strcpy(tmpName,file->ReadString());
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
		mApp->AddView(tmpName,v);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("selections");file->ReadKeyword("{");
	for(int i=0;i<nbS;i++) {
		ASELECTION s;
		char tmpName[256];
		strcpy(tmpName,file->ReadString());
		s.name=_strdup(tmpName);
		s.nbSel = file->ReadInt();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));

		for (int j=0;j<s.nbSel;j++) {
			s.selection[j]=file->ReadInt()+sh.nbFacet;
		}
		mApp->AddSelection(s.name,s);
	}
	file->ReadKeyword("}");

	for (int i = 0; i < nbF; i++) { //parse formulas now that selection groups are loaded
		mApp->AddFormula(loadFormulas[i][0].c_str(), loadFormulas[i][1].c_str());
	}

	file->ReadKeyword("structures");file->ReadKeyword("{");
	for(int i=0;i<nbNewSuper;i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		/*sprintf(tmp,"%s.txt",strName[i]);
		strFileName[i] = _strdup(tmp);*/
	}
	file->ReadKeyword("}");

	// Reallocate memory
	*facets = (Facet **)realloc(*facets,(nbNewFacets+*nbFacet) * sizeof(Facet **));
	memset(*facets+*nbFacet,0,nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3,*vertices3,(*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3+*nbVertex,0,nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3=tmp_vertices3;

	// Read geometry vertices
	file->ReadKeyword("vertices");file->ReadKeyword("{");
	for(int i=*nbVertex;i<(*nbVertex+nbNewVertex);i++) {
		// Check idx
		int idx = file->ReadInt();
		if( idx != i-*nbVertex+1 ) throw Error(file->MakeError("Wrong vertex index !"));
		(*vertices3+i)->x = file->ReadDouble();
		(*vertices3+i)->y = file->ReadDouble();
		(*vertices3+i)->z = file->ReadDouble();
		(*vertices3+i)->selected = FALSE;
	}
	file->ReadKeyword("}");

	
		// Read leaks
		file->ReadKeyword("leaks");file->ReadKeyword("{");
		file->ReadKeyword("nbLeak");file->ReadKeyword(":");
		int nbleak_local = file->ReadInt();
		for (int i=0;i<nbleak_local;i++) {
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
		file->ReadKeyword("hits");file->ReadKeyword("{");
		file->ReadKeyword("nbHHit");file->ReadKeyword(":");
		int nbHHit_local = file->ReadInt();
		for (int i=0;i<nbHHit_local;i++) {
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
	for(int i=*nbFacet;i<(*nbFacet+nbNewFacets);i++) {
		file->ReadKeyword("facet");  
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1-*nbFacet ) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nb = file->ReadInt();

		if( nb<3 ) {
			char errMsg[512];
			sprintf(errMsg,"Facet %d has only %d vertices. ",i,nb);
			throw Error(errMsg);
		}

		*(*facets+i) = new Facet(nb);
		(*facets)[i]->LoadSYN(file,version2,nbNewVertex);
		(*facets)[i]->selected = TRUE;
		for(int j=0;j<nb;j++)
			(*facets)[i]->indices[j]+=*nbVertex;
		file->ReadKeyword("}");
		if (newStruct) {
			(*facets)[i]->sh.superIdx += sh.nbSuper;
		} else {
			(*facets)[i]->sh.superIdx = strIdx;
		}
	}


	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) sh.nbSuper+=nbNewSuper;
	return result;
}

void Geometry::InsertSTLGeom(FileReader *file,int *nbVertex,int *nbFacet,VERTEX3D **vertices3,Facet ***facets,int strIdx,BOOL newStruct,double scaleFactor) {

	UnSelectAll();
	char *w;

	int nbNewFacets=0;
	// First pass
	file->ReadKeyword("solid");
	file->ReadLine(); // solid name
	w = file->ReadWord();
	while(strcmp(w,"facet")==0) {
		nbNewFacets++;
		file->JumpSection("endfacet");
		w = file->ReadWord();
	}
	if( strcmp(w,"endsolid")!=0 ) throw Error("Unexpected or not supported STL keyword, 'endsolid' required");

	// Allocate memory
	int nbNewVertex=3*nbNewFacets;
	*facets = (Facet **)realloc(*facets,(nbNewFacets+*nbFacet) * sizeof(Facet **));
	memset(*facets+*nbFacet,0,nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3,*vertices3,(*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3+*nbVertex,0,nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3=tmp_vertices3;

	// Second pass
	file->SeekStart();
	file->ReadKeyword("solid");
	file->ReadLine();
	for(int i=0;i<nbNewFacets;i++) {

		file->ReadKeyword("facet");
		file->ReadKeyword("normal");
		file->ReadDouble(); //ignoring normal vector, will be calculated from triangle orientation
		file->ReadDouble();
		file->ReadDouble();
		file->ReadKeyword("outer");
		file->ReadKeyword("loop");

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex+3*i+0].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+0].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+0].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex+3*i+1].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+1].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+1].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex+3*i+2].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+2].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+2].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("endloop");
		file->ReadKeyword("endfacet");

		*(*facets+i+*nbFacet) = new Facet(3);
		(*facets)[i+*nbFacet]->selected = TRUE;
		(*facets)[i+*nbFacet]->indices[0] = *nbVertex+3*i+0;
		(*facets)[i+*nbFacet]->indices[1] = *nbVertex+3*i+1;
		(*facets)[i+*nbFacet]->indices[2] = *nbVertex+3*i+2;

		if (newStruct) {
			(*facets)[i+*nbFacet]->sh.superIdx = sh.nbSuper;
		} else {
			(*facets)[i+*nbFacet]->sh.superIdx = strIdx;
		}
	}

	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) AddStruct("Inserted STL file");

}


void Geometry::SaveProfileTXT(FileWriter *file,int super,BOOL saveSelected) {
	// Profiles
	for(int j=0;j<PROFILE_SIZE;j++)
		file->Write("\n");
}

void Geometry::SaveProfileGEO(FileWriter *file,int super,BOOL saveSelected) {
	file->Write("profiles {\n");
	// Profiles
	int nbProfile=0;
	file->Write(" number: ");file->WriteInt(nbProfile,"\n");
	file->Write(" facets: ");
	file->Write("\n");
	for(int j=0;j<PROFILE_SIZE;j++)
		file->Write("\n");
	file->Write("}\n");
}





































void Geometry::LoadGEO(FileReader *file,GLProgress *prg,LEAK *pleak,int *nbleak,HIT *pHits,int *nbHHit,int *version) {

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	prg->SetMessage("Clearing current geometry...");
	Clear();
	//mApp->ClearFormula();


	// Globals
	char tmp[512];
	prg->SetMessage("Reading GEO file header...");
	file->ReadKeyword("version");file->ReadKeyword(":");
	*version = file->ReadInt();
	if( *version>GEOVERSION ) {
		char errMsg[512];
		sprintf(errMsg,"Unsupported GEO version V%d",*version);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit");file->ReadKeyword(":");
	tNbHit =0;file->ReadLLong();
	file->ReadKeyword("totalDes");file->ReadKeyword(":");
	tNbDesorption =0;file->ReadLLong();
	file->ReadKeyword("totalLeak");file->ReadKeyword(":");
	tNbLeak = 0;file->ReadLLong();
	if (*version>=12) {
		file->ReadKeyword("totalAbs");file->ReadKeyword(":");
		tNbAbsorption =0; file->ReadLLong();
		if (*version >= 15) {
			file->ReadKeyword("totalDist_total");
		}
		else { //between versions 12 and 15
			file->ReadKeyword("totalDist");
		}
		file->ReadKeyword(":");
		distTraveledTotal =0.0; file->ReadDouble();
		if (*version >= 15) {
			file->ReadKeyword("totalDist_fullHitsOnly"); file->ReadKeyword(":");
			file->ReadDouble();
		}
	} else {
		tNbAbsorption=0;
		distTraveledTotal=0.0;
	}
	file->ReadKeyword("maxDes");file->ReadKeyword(":");
	tNbDesorptionMax = file->ReadLLong(); 
	file->ReadKeyword("nbVertex");file->ReadKeyword(":");
	sh.nbVertex = file->ReadInt();
	file->ReadKeyword("nbFacet");file->ReadKeyword(":");
	sh.nbFacet = file->ReadInt();
	file->ReadKeyword("nbSuper");file->ReadKeyword(":");
	sh.nbSuper = file->ReadInt();
	int nbF=0; std::vector<std::vector<string>> loadFormulas;
	int nbV=0;
	if (*version>=2) {
		file->ReadKeyword("nbFormula");file->ReadKeyword(":");
		nbF = file->ReadInt(); loadFormulas.reserve(nbF);
		file->ReadKeyword("nbView");file->ReadKeyword(":");
		nbV = file->ReadInt();
	}
	int nbS=0;
	if (*version>=8) {
		file->ReadKeyword("nbSelection");file->ReadKeyword(":");
		nbS = file->ReadInt();
	}
	if (*version>=7) {
		file->ReadKeyword("gasMass");file->ReadKeyword(":");
		file->ReadDouble(); //gas mass
	}
		if (*version>=10) { //time-dependent version
		file->ReadKeyword("userMoments");file->ReadKeyword("{");
		file->ReadKeyword("nb");file->ReadKeyword(":");
		int nb=file->ReadInt();

		for(int i=0;i<nb;i++) {
			char tmpExpr[512];
			strcpy(tmpExpr,file->ReadString());
			//mApp->worker.userMoments.push_back(tmpExpr);
			//mApp->worker.AddMoment(mApp->worker.ParseMoment(tmpExpr));
		}
		file->ReadKeyword("}");
	}
	if (*version>=11) { //pulse version
		file->ReadKeyword("desorptionStart");file->ReadKeyword(":");
		//worker->desorptionStartTime=
		file->ReadDouble();
		file->ReadKeyword("desorptionStop");file->ReadKeyword(":");
		//worker->desorptionStopTime=
		file->ReadDouble();
		file->ReadKeyword("timeWindow");file->ReadKeyword(":");
		//worker->timeWindowSize=
		file->ReadDouble();
		file->ReadKeyword("useMaxwellian");file->ReadKeyword(":");
		//worker->useMaxwellDistribution=
		file->ReadInt();
	}
	if (*version>=12) { //2013.aug.22
		file->ReadKeyword("calcConstantFlow");file->ReadKeyword(":");
		//worker->calcConstantFlow=
		file->ReadInt();
	}
	if (*version>=2) {
		file->ReadKeyword("formulas");file->ReadKeyword("{");
		for(int i=0;i<nbF;i++) {
			char tmpName[256];
			char tmpExpr[512];
			strcpy(tmpName,file->ReadString());
			strcpy(tmpExpr,file->ReadString());
			//mApp->AddFormula(tmpName,tmpExpr);
			std::vector<string> newFormula;
			newFormula.push_back(tmpName);
			newFormula.push_back(tmpExpr);
			loadFormulas.push_back(newFormula);
		}
		file->ReadKeyword("}");

		file->ReadKeyword("views");file->ReadKeyword("{");
		for(int i=0;i<nbV;i++) {
			char tmpName[256];
			AVIEW v;
			strcpy(tmpName,file->ReadString());
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
			mApp->AddView(tmpName,v);
		}
		file->ReadKeyword("}");
	}

	if (*version>=8) {
		file->ReadKeyword("selections");file->ReadKeyword("{");
		for(int i=0;i<nbS;i++) {
			ASELECTION s;
			char tmpName[256];
			strcpy(tmpName,file->ReadString());
			s.name=_strdup(tmpName);
			s.nbSel = file->ReadInt();
			s.selection = (int *)malloc((s.nbSel)*sizeof(int));

			for (int j=0;j<s.nbSel;j++) {
				s.selection[j]=file->ReadInt();
			}
			mApp->AddSelection(s.name,s);
		}
		file->ReadKeyword("}");
	}

	for (int i = 0; i < nbF; i++) { //parse formulas now that selection groups are loaded
		mApp->AddFormula(loadFormulas[i][0].c_str(), loadFormulas[i][1].c_str());
	}
	
	file->ReadKeyword("structures");file->ReadKeyword("{");
	for(int i=0;i<sh.nbSuper;i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		sprintf(tmp,"%s.txt",strName[i]);
		strFileName[i] = _strdup(tmp);
	}
	file->ReadKeyword("}");

	// Allocate memory
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets,0,sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memset(vertices3,0,sh.nbVertex * sizeof(VERTEX3D));

	// Read vertices
	prg->SetMessage("Reading vertices...");
	file->ReadKeyword("vertices");file->ReadKeyword("{");
	for(int i=0;i<sh.nbVertex;i++) {
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1 ) throw Error(file->MakeError("Wrong vertex index !"));
		vertices3[i].x = file->ReadDouble();
		vertices3[i].y = file->ReadDouble();
		vertices3[i].z = file->ReadDouble();
		vertices3[i].selected = FALSE;
	}
	file->ReadKeyword("}");

	if (*version>=6) {
		prg->SetMessage("Reading leaks and hits...");
		// Read leaks
		file->ReadKeyword("leaks");file->ReadKeyword("{");
		file->ReadKeyword("nbLeak");file->ReadKeyword(":");
		*nbleak = file->ReadInt();
		for (int i=0;i<*nbleak;i++) {
			int idx = file->ReadInt();
			if( idx != i ) throw Error(file->MakeError("Wrong leak index !"));
			(pleak+i)->pos.x = file->ReadDouble();
			(pleak+i)->pos.y = file->ReadDouble();
			(pleak+i)->pos.z = file->ReadDouble();

			(pleak+i)->dir.x = file->ReadDouble();
			(pleak+i)->dir.y = file->ReadDouble();
			(pleak+i)->dir.z = file->ReadDouble();   
		}
		file->ReadKeyword("}");

		// Read hit cache
		file->ReadKeyword("hits");file->ReadKeyword("{");
		file->ReadKeyword("nbHHit");file->ReadKeyword(":");
		*nbHHit = file->ReadInt();
		for (int i=0;i<*nbHHit;i++) {
			int idx = file->ReadInt();
			if( idx != i ) throw Error(file->MakeError("Wrong hit cache index !"));
			(pHits+i)->pos.x = file->ReadDouble();
			(pHits+i)->pos.y = file->ReadDouble();
			(pHits+i)->pos.z = file->ReadDouble();

			(pHits+i)->type = file->ReadInt();
		}
		file->ReadKeyword("}");
	}

	// Read facets
	prg->SetMessage("Reading facets...");
	for(int i=0;i<sh.nbFacet;i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1 ) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nbI = file->ReadInt();
		if( nbI<3 ) {
			char errMsg[512];
			sprintf(errMsg,"Facet %d has only %d vertices. ",i,nbI);
			throw Error(errMsg);
		}
		prg->SetProgress((float)i/sh.nbFacet);
		facets[i] = new Facet(nbI);
		facets[i]->LoadGEO(file,*version,sh.nbVertex);
		file->ReadKeyword("}");
	}

	InitializeGeometry();
	//AdjustProfile();
	isLoaded=TRUE;
	UpdateName(file);

	// Update mesh //Unnecessary in SynRad
	/*prg->SetMessage("Building mesh...");
	for(int i=0;i<sh.nbFacet;i++) {
		double p = (double)i/(double)sh.nbFacet;
		prg->SetProgress(p);
		Facet *f = facets[i];
		f->SetTexture(f->sh.texWidthD,f->sh.texHeightD,f->hasMesh); //here texWidthD
		BuildFacetList(f);
		double nU = Norme(&(f->sh.U));
		f->tRatio = f->sh.texWidthD / nU;
	}*/
}

// -----------------------------------------------------------
bool Geometry::LoadTextures(FileReader *file,GLProgress *prg,Dataport *dpHit,int version) {

	if(file->SeekFor("{textures}")) {
		char tmp[256];
		//versions 3+
		// Block dpHit during the whole disc reading

		AccessDataport(dpHit);

		// Globals
		BYTE *buffer = (BYTE *)dpHit->buff;
		SHGHITS *gHits = (SHGHITS *)buffer;

		gHits->total.nbHit=tNbHit;
		gHits->total.nbDesorbed=tNbDesorption;
		gHits->total.nbAbsorbed=tNbAbsorption;
		gHits->nbLeakTotal = tNbLeak;
		gHits->total.fluxAbs=tFlux;
		gHits->total.powerAbs=tPower;
		gHits->distTraveledTotal = distTraveledTotal;

		// Read facets
		file->ReadKeyword("minHit_MC");file->ReadKeyword(":");
		gHits->minHit_MC = file->ReadLLong();
		file->ReadKeyword("maxHit_MC");file->ReadKeyword(":");
		gHits->maxHit_MC = file->ReadLLong();
		file->ReadKeyword("minHit_flux");file->ReadKeyword(":");
		gHits->minHit_flux = file->ReadDouble();
		file->ReadKeyword("maxHit_flux");file->ReadKeyword(":");
		gHits->maxHit_flux = file->ReadDouble();
		file->ReadKeyword("minHit_power");file->ReadKeyword(":");
		gHits->minHit_power = file->ReadDouble();
		file->ReadKeyword("maxHit_power");file->ReadKeyword(":");
		gHits->maxHit_power = file->ReadDouble();

		for(int i=0;i<sh.nbFacet;i++) {
			Facet *f = facets[i];
			if (f->hasMesh) {
				prg->SetProgress((double)i/(double)sh.nbFacet);
				file->ReadKeyword("texture_facet");
				// Check idx
				int idx = file->ReadInt();

				if( idx != i+1 ) {
					sprintf(tmp,"Wrong facet index. Expected %d, read %d.",i+1,idx);
					throw Error(file->MakeError(tmp));
				}

				//Now load values
				file->ReadKeyword("{");

				int ix,iy;

				int profSize = (f->sh.isProfile)?(PROFILE_SIZE*(sizeof(llong)+2*sizeof(double))):0;
				int nbE=f->sh.texHeight*f->sh.texWidth;
				llong *hits_MC = (llong *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize));
				double *hits_flux = (double *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + nbE*sizeof(llong)));
				double *hits_power = (double *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize+nbE*(sizeof(llong)+sizeof(double))));

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

				for (iy = 0; iy<(MIN(f->sh.texHeight, texHeight_file)); iy++) { //MIN: If stored texture is larger, don't read extra cells
					for (ix=0;ix<(MIN(f->sh.texWidth,texWidth_file));ix++) { //MIN: If stored texture is larger, don't read extra cells
						int index=iy*(f->sh.texWidth)+ix;
						*(hits_MC+index)=file->ReadLLong();
						if (version >= 7) file->ReadDouble(); //cell area
						*(hits_flux+index)=file->ReadDouble();
						*(hits_power+index)=file->ReadDouble();

						//Normalize by area
						if(f->mesh[index].area>0.0) {
							*(hits_flux+index)/=f->mesh[index].area;
							*(hits_power+index)/=f->mesh[index].area;
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

void Geometry::SaveTXT(FileWriter *file,Dataport *dpHit,BOOL saveSelected) {

	if(!IsLoaded()) throw Error("Nothing to save !");

	// Unused
	file->WriteInt(0,"\n");

	// Block dpHit during the whole disc writing
	AccessDataport(dpHit);

	// Globals
	//BYTE *buffer = (BYTE *)dpHit->buff;
	//SHGHITS *gHits = (SHGHITS *)buffer;

	// Unused
	file->WriteLLong(0,"\n");
	file->WriteInt(0,"\n");
	file->WriteLLong(0,"\n");
	file->WriteLLong(tNbDesorptionMax,"\n");

	file->WriteInt(sh.nbVertex,"\n");
	file->WriteInt(saveSelected?nbSelected:sh.nbFacet,"\n");

	// Read geometry vertices
	for(int i=0;i<sh.nbVertex;i++) {
		file->WriteDouble(vertices3[i].x," ");
		file->WriteDouble(vertices3[i].y," ");
		file->WriteDouble(vertices3[i].z,"\n");
	}

	// Facets
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		int j;
		if( saveSelected ) {
			if( f->selected ) {
				file->WriteInt(f->sh.nbIndex," ");
				for(j=0;j<f->sh.nbIndex-1;j++)
					file->WriteInt(f->indices[j]+1," ");
				file->WriteInt(f->indices[j]+1,"\n");
			}
		} else {
			file->WriteInt(f->sh.nbIndex," ");
			for(j=0;j<f->sh.nbIndex-1;j++)
				file->WriteInt(f->indices[j]+1," ");
			file->WriteInt(f->indices[j]+1,"\n");
		}
	}

	// Params
	for(int i=0;i<sh.nbFacet;i++) {

		// Update facet hits from shared mem
		Facet *f = facets[i];
		//SHHITS *shF = (SHHITS *)(buffer + f->sh.hitOffset);
		//memcpy(&(f->sh.counter),shF,sizeof(SHHITS));
		if( saveSelected ) {
			if(f->selected) f->SaveTXT(file);
		} else {
			f->SaveTXT(file);
		}

	}

	SaveProfileTXT(file,-1,saveSelected);

	ReleaseDataport(dpHit);

}

// -----------------------------------------------------------------------
void Geometry::ExportTextures(FILE *file,int grouping,int mode,double no_scans,Dataport *dpHit,BOOL saveSelected) {

	//if(!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing
	BYTE *buffer = NULL;
	if( dpHit )
		if( AccessDataport(dpHit) ) 
			buffer = (BYTE *)dpHit->buff;

	// Globals
	//BYTE *buffer = (BYTE *)dpHit->buff;
	//SHGHITS *gHits = (SHGHITS *)buffer;

	if (grouping==1) fprintf(file,"X_coord_cm\tY_coord_cm\tZ_coord_cm\tValue\t\n"); //mode 10: special ANSYS export

	// Facets
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];

		if( f->selected ) {
			if (grouping==0) fprintf(file,"FACET%d\n",i+1); //mode 10: special ANSYS export
			SHELEM *mesh = f->mesh;
			if( mesh ) {
				char tmp[256];
				char out[256];
				float dCoef = 1.0f;
				if(!buffer) return;
				SHGHITS *shGHit = (SHGHITS *)buffer;
				int w = f->sh.texWidth;
				int h = f->sh.texHeight;
				int nbE = w*h;
				int profSize = (f->sh.isProfile) ? (PROFILE_SIZE*(sizeof(llong)+2 * sizeof(double))) : 0;

				llong *hits_MC = (llong *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize));
				double *hits_flux = (double *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize + nbE*sizeof(llong)));
				double *hits_power = (double *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize+nbE*(sizeof(llong)+sizeof(double))));
				double norm = 1.0 / no_scans; //normalize values by number of scans (and don't normalize by area...)

				for(int i=0;i<w;i++) {
					for(int j=0;j<h;j++) {
						int index=i+j*w;
						tmp[0]=out[0]=0;
						switch(mode) {

						case 0: // Element area
							sprintf(tmp,"%g",f->mesh[index].area);
							break;

						case 1: // MC_hits
							if (!grouping || hits_MC[index]) sprintf(tmp, "%I64d", hits_MC[index]);
							break;

						case 2: // Flux
							if (!grouping || hits_flux[index]) sprintf(tmp, "%g", hits_flux[index] * f->mesh[i + j*w].area*norm);
							break;

						case 3: // Power
							if (!grouping || hits_power[index]) sprintf(tmp, "%g", hits_power[index] * f->mesh[i + j*w].area*norm);
							break;

						case 4: // Flux/area
							if (!grouping || hits_flux[index]) sprintf(tmp, "%g", hits_flux[index] * norm);
							break;

						case 5: // Power/area
							if (!grouping || hits_power[index]) sprintf(tmp, "%g", hits_power[index] * 0.01*norm); //Don't write 0 powers
							break;
						}

						if (grouping==1 && tmp && tmp[0])
							sprintf(out,"%g\t%g\t%g\t%s\t\n",
								f->sh.O.x+f->mesh[index].uCenter*f->sh.U.x+f->mesh[index].vCenter*f->sh.V.x,
								f->sh.O.y+f->mesh[index].uCenter*f->sh.U.y+f->mesh[index].vCenter*f->sh.V.y,
								f->sh.O.z+f->mesh[index].uCenter*f->sh.U.z+f->mesh[index].vCenter*f->sh.V.z,
								tmp);
						else sprintf(out,"%s",tmp);

						if( out ) fprintf(file,"%s",out);
						if( j<w-1 && grouping==0) 
							fprintf(file,"\t");
					}
					if (grouping==0) fprintf(file,"\n");
				}
			}
			else {
				fprintf(file,"No mesh.\n");
			}
			if (grouping==0) fprintf(file,"\n"); //Current facet exported. 
		}

	}

	ReleaseDataport(dpHit);

}

void Geometry::SaveDesorption(FILE *file,Dataport *dpHit,BOOL selectedOnly,int mode,double eta0,double alpha,Distribution2D *distr) {

	if(!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing
	BYTE *buffer = NULL;
	if( dpHit )
		if( AccessDataport(dpHit) ) 
			buffer = (BYTE *)dpHit->buff;
	
	Worker *worker=&(mApp->worker);
	
	// Globals
	//BYTE *buffer = (BYTE *)dpHit->buff;
	//SHGHITS *gHits = (SHGHITS *)buffer;

	// Facets
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];

		if( f->selected || !selectedOnly) {
			
			SHELEM *mesh = f->mesh;
			if( mesh ) {
				fprintf(file,"facet %d {\n",i+1);
				char tmp[256];
				fprintf(file,"cell_size_cm: %g\n",1.0/f->tRatio);
				float dCoef = 1.0f;
				if(!buffer) return;
				
				int w = f->sh.texWidth;
				int h = f->sh.texHeight;
				int textureSize_double=w*h*sizeof(double);
				int textureSize_llong=w*h*sizeof(llong);
				SHGHITS *shGHit = (SHGHITS *)buffer;
				int profile_memory=PROFILE_SIZE*(2*sizeof(double)+sizeof(llong));
				int profSize = (f->sh.isProfile)?profile_memory:0;
				double *hits_flux = (double *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize + textureSize_llong));


					for(int i=0;i<w;i++) {
						for(int j=0;j<h;j++) {
							double dose=hits_flux[i+j*w]*f->mesh[i+j*w].area/worker->no_scans;
							double val;
							if (mode==1) { //no conversion
								val=dose;
							} else if (mode==2) { //use equation
								double eta;
								if (dose < VERY_SMALL) eta=0.0;	
								else eta=eta0*pow(dose,alpha);
								val=dose*eta;
							} else if (mode==3) {  //use file
								double eta=distr->InterpolateY(dose);
								val=dose*eta;
							}
							sprintf(tmp,"%g",val);
							if( tmp ) fprintf(file,"%s",tmp);
							//fprintf(file,"%s",tmp);
							if( j<h-1 ) 
								fprintf(file,"\t");
						}
						fprintf(file,"\n");
					}
			fprintf(file,"}\n"); //Current facet exported. 
			}
			
		}

	}

	ReleaseDataport(dpHit);

}



void Geometry::SaveSTR(Dataport *dpHit,BOOL saveSelected) {

	if(!IsLoaded()) throw Error("Nothing to save !");
	if(sh.nbSuper<1) throw Error("Cannot save single structure in STR format");

	// Block dpHit during the whole disc writting
	AccessDataport(dpHit);
	for(int i=0;i<sh.nbSuper;i++)
		SaveSuper(dpHit,i);
	ReleaseDataport(dpHit);

}




void Geometry::SaveSuper(Dataport *dpHit,int s) {

	char fName[512];
	sprintf(fName,"%s/%s",strPath,strFileName[s]);
	FileWriter *file = new FileWriter(fName);

	// Unused
	file->WriteInt(0,"\n");

	// Globals
	BYTE *buffer = (BYTE *)dpHit->buff;
	SHGHITS *gHits = (SHGHITS *)buffer;

	//Extract data of the specified super structure
	llong totHit  = 0;
	llong totAbs  = 0;
	llong totDes  = 0;
	int *refIdx = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(refIdx,0xFF,sh.nbVertex*sizeof(int));
	int nbV=0;
	int nbF=0;

	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		if( f->sh.superIdx==s ) {
			totHit += f->sh.counter.nbHit;
			totAbs += f->sh.counter.nbAbsorbed;
			totDes += f->sh.counter.nbDesorbed;
			for(int j=0;j<f->sh.nbIndex;j++)
				refIdx[f->indices[j]] = 1;
			nbF++;
		}
	}

	for(int i=0;i<sh.nbVertex;i++) {
		if(refIdx[i]>=0) {
			refIdx[i] = nbV;
			nbV++;
		}
	}

	file->WriteLLong(0,"\n");
	file->WriteInt(0,"\n");
	file->WriteLLong(0,"\n");
	file->WriteLLong(0,"\n");

	file->WriteInt(nbV,"\n");
	file->WriteInt(nbF,"\n");

	// Read geometry vertices
	for(int i=0;i<sh.nbVertex;i++) {
		if( refIdx[i]>=0 ) {
			file->WriteDouble(vertices3[i].x," ");
			file->WriteDouble(vertices3[i].y," ");
			file->WriteDouble(vertices3[i].z,"\n");
		}
	}

	// Facets
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		int j;
		if( f->sh.superIdx==s ) {
			file->WriteInt(f->sh.nbIndex," ");
			for(j=0;j<f->sh.nbIndex-1;j++)
				file->WriteInt(refIdx[f->indices[j]]+1," ");
			file->WriteInt(refIdx[f->indices[j]]+1,"\n");
		}
	}

	// Params
	for(int i=0;i<sh.nbFacet;i++) {

		// Update facet hits from shared mem
		Facet *f = facets[i];
		if( f->sh.superIdx==s ) {
			SHHITS *shF = (SHHITS *)(buffer + f->sh.hitOffset);
			memcpy(&(f->sh.counter),shF,sizeof(SHHITS));
			f->SaveTXT(file);
		}

	}

	SaveProfileTXT(file,s);

	SAFE_DELETE(file);
	free(refIdx);

}

void Geometry::RemoveFromStruct(int numToDel) {
	mApp->changedSinceSave = TRUE;

	int nb = 0;
	for (int i = 0; i < sh.nbFacet; i++)
		if (facets[i]->sh.superIdx == numToDel) nb++;


	if (nb == 0) return;

	Facet   **f = (Facet **)malloc((sh.nbFacet - nb) * sizeof(Facet *));

	nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->sh.superIdx == numToDel) {

			delete facets[i];
			mApp->RenumberSelections(nb);
			mApp->RenumberFormulas(nb);
		}
		else {

			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

}

BOOL Geometry::IsLoaded() {
	return isLoaded;
}

void Geometry::SaveSYN(FileWriter *file,GLProgress *prg,Dataport *dpHit,BOOL saveSelected,LEAK *pleak,
					   int *nbleakSave,HIT *pHits,int *nbHHitSave,BOOL crashSave) {

	prg->SetMessage("Counting hits...");
	if(!IsLoaded()) throw Error("Nothing to save !");


	// Block dpHit during the whole disc writing
	if (!crashSave && !saveSelected) AccessDataport(dpHit);

	// Globals
	BYTE *buffer;
	if (!crashSave && !saveSelected) buffer = (BYTE *)dpHit->buff;
	SHGHITS *gHits; 
	if (!crashSave && !saveSelected) gHits = (SHGHITS *)buffer;


	double dCoef = 1.0;
	int ix,iy;

	Worker *worker=&(mApp->worker);
	texCMin_MC = (!crashSave && !saveSelected)?gHits->minHit_MC:0;
	texCMax_MC = (!crashSave && !saveSelected)?gHits->maxHit_MC:1;
	texCMin_flux = (!crashSave && !saveSelected)?gHits->minHit_flux:0;
	texCMax_flux = (!crashSave && !saveSelected)?gHits->maxHit_flux:1;
	texCMin_power = (!crashSave && !saveSelected)?gHits->minHit_power:0;
	texCMax_power = (!crashSave && !saveSelected)?gHits->maxHit_power:1;

	prg->SetMessage("Writing geometry details...");
	file->Write("version:");file->WriteInt(SYNVERSION,"\n");
	file->Write("totalHit:");file->WriteLLong((!crashSave && !saveSelected)?gHits->total.nbHit:0,"\n");
	file->Write("totalDes:");file->WriteLLong((!crashSave && !saveSelected)?gHits->total.nbDesorbed:0,"\n");
	file->Write("no_scans:");file->WriteDouble((!crashSave && !saveSelected)?worker->no_scans:0, "\n");
	file->Write("totalLeak:");file->WriteLLong((!crashSave && !saveSelected)?gHits->nbLeakTotal:0,"\n");
	file->Write("totalFlux:");file->WriteDouble((!crashSave && !saveSelected)?gHits->total.fluxAbs:0,"\n");
	file->Write("totalPower:");file->WriteDouble((!crashSave && !saveSelected)?gHits->total.powerAbs:0,"\n");
	file->Write("maxDes:");file->WriteLLong((!crashSave && !saveSelected)?tNbDesorptionMax:0,"\n");
	file->Write("nbVertex:");file->WriteInt(sh.nbVertex,"\n");
	file->Write("nbFacet:");file->WriteInt(saveSelected?nbSelected:sh.nbFacet,"\n");
	file->Write("nbSuper:");file->WriteInt(sh.nbSuper,"\n");
	file->Write("nbFormula:");file->WriteInt((!saveSelected)?mApp->nbFormula:0,"\n");
	file->Write("nbView:");file->WriteInt(mApp->nbView,"\n");
	file->Write("nbSelection:");file->WriteInt((!saveSelected)?mApp->nbSelection:0,"\n");

	file->Write("nbRegions:");file->WriteInt((!saveSelected)?(int)worker->regions.size():0,"\n");
	file->Write("PARfiles {\n");
	for(int i=0;(!saveSelected)&&(i<(int)worker->regions.size());i++) { //write only the filenames without the path
		file->Write("  \"");
		char tmp[512];
		size_t foundDir=worker->regions[i].fileName.rfind('\\');
		if (foundDir!=string::npos) {//found
			strcpy(tmp,&(worker->regions[i].fileName.c_str()[(int)foundDir+1]));
		} else {//no "\" sign
			strcpy(tmp,worker->regions[i].fileName.c_str());
		}
		file->Write(tmp);
		file->Write("\"\n");
	}
	file->Write("}\n");
	
	file->Write("formulas {\n");
	for(int i=0;(!saveSelected)&&(i<mApp->nbFormula);i++) {
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
	for(int i=0;(!saveSelected)&&(i<mApp->nbSelection);i++) {
		file->Write("  \"");
		file->Write(mApp->selections[i].name);
		file->Write("\"\n ");
		file->WriteInt(mApp->selections[i].nbSel,"\n");
		for (int j=0;j<mApp->selections[i].nbSel;j++) {
			file->Write("  ");
			file->WriteInt(mApp->selections[i].selection[j],"\n");
		}
		//file->Write("\n");
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
	file->Write("  nbLeak:"); file->WriteInt((!crashSave && !saveSelected) ? *nbleakSave : 0, "\n");
	for (int i = 0; (i < *nbleakSave) && (!crashSave && !saveSelected); i++) {

		file->Write("  ");
		file->WriteInt(i," ");
		file->WriteDouble((pleak+i)->pos.x," ");
		file->WriteDouble((pleak+i)->pos.y," ");
		file->WriteDouble((pleak+i)->pos.z," ");

		file->WriteDouble((pleak+i)->dir.x," ");
		file->WriteDouble((pleak+i)->dir.y," ");
		file->WriteDouble((pleak+i)->dir.z,"\n");
	}
	file->Write("}\n");

	//hit cache (lines and dots)
	prg->SetMessage("Writing hit cache...");
	file->Write("hits {\n");
	file->Write("  nbHHit:"); file->WriteInt((!crashSave && !saveSelected) ? *nbHHitSave : 0, "\n");
	for (int i = 0; (i < *nbHHitSave) && (!crashSave && !saveSelected); i++) {

		file->Write("  ");
		file->WriteInt(i," ");
		file->WriteDouble((pHits+i)->pos.x," ");
		file->WriteDouble((pHits+i)->pos.y," ");
		file->WriteDouble((pHits+i)->pos.z," ");
		file->WriteDouble((pHits+i)->dF," ");
		file->WriteDouble((pHits+i)->dP," ");
		file->WriteInt((pHits+i)->type,"\n");
	}
	file->Write("}\n");

	//facets

	prg->SetMessage("Writing facets..."); 

	for(int i=0,k=0;i<sh.nbFacet;i++) {
		prg->SetProgress(0.33+((double)i/(double)sh.nbFacet) *0.33);
		if( !saveSelected || facets[i]->selected ) { facets[i]->SaveSYN(file,k,crashSave);k++; }	
	}

	prg->SetMessage("Writing profiles and spectrum...");
	SaveProfileSYN(file,dpHit,-1,saveSelected,crashSave);
	SaveSpectrumSYN(file,dpHit,-1,saveSelected,crashSave);

	///Save textures, for GEO file version 3+
	char tmp[256];
	file->Write("{textures}\n");
	
	file->Write("minHit_MC:");file->WriteLLong(texCMin_MC,"\n");
	file->Write("maxHit_MC:");file->WriteLLong(texCMax_MC,"\n");
	file->Write("minHit_flux:");file->WriteDouble(texCMin_flux,"\n");
	file->Write("maxHit_flux:");file->WriteDouble(texCMax_flux,"\n");
	file->Write("minHit_power:");file->WriteDouble(texCMin_power,"\n");
	file->Write("maxHit_power:");file->WriteDouble(texCMax_power,"\n");

	//Selections
	//SaveSelections();

	prg->SetMessage("Writing textures..."); 
	for(int i=0;i<sh.nbFacet;i++) {
		prg->SetProgress(((double)i/(double)sh.nbFacet) *0.33+0.66);
		Facet *f = facets[i];
		if (f->hasMesh) {
			int profSize = (f->sh.isProfile)?(PROFILE_SIZE*(sizeof(llong)+2*sizeof(double))):0;
			int nbE=f->sh.texHeight*f->sh.texWidth;
			llong *hits_MC;
			if (!crashSave && !saveSelected) hits_MC = (llong *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize));
			double *hits_flux; 
			if (!crashSave && !saveSelected) hits_flux = (double *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize+nbE*sizeof(llong)));
			double *hits_power; 
			if (!crashSave && !saveSelected) hits_power = (double *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize+nbE*(sizeof(llong)+sizeof(double))));

			//char tmp[256];
			sprintf(tmp,"texture_facet %d {\n",i+1);
			file->Write(tmp);
			file->Write("width:"); file->WriteInt(f->sh.texWidth); file->Write(" height:"); file->WriteInt(f->sh.texHeight); file->Write("\n");
			for (iy=0;iy<(f->sh.texHeight);iy++) {
				for (ix=0;ix<(f->sh.texWidth);ix++) {
					int index=iy*(f->sh.texWidth)+ix;
					file->WriteLLong((!crashSave && !saveSelected)?*(hits_MC+index):0,"\t");
					file->WriteDouble(f->mesh[index].area, "\t");
					file->WriteDouble((!crashSave && !saveSelected)?*(hits_flux+index)*f->mesh[index].area:0,"\t");
					file->WriteDouble((!crashSave && !saveSelected)?*(hits_power+index)*f->mesh[index].area:0,"\t");
				}
				file->Write("\n");
			}
			file->Write("}\n");
		}
	}

	if (!crashSave && !saveSelected) ReleaseDataport(dpHit);

}


PARfileList Geometry::LoadSYN(FileReader *file,GLProgress *prg,LEAK *pleak,int *nbleak,HIT *pHits,int *nbHHit,int *version) {

	prg->SetMessage("Clearing current geometry...");
	Clear();

	PARfileList result(0);

	// Globals
	char tmp[512];
	prg->SetMessage("Reading SYN file header...");
	file->ReadKeyword("version");file->ReadKeyword(":");
	*version = file->ReadInt();
	if( *version>SYNVERSION ) {
		char errMsg[512];
		sprintf(errMsg,"Unsupported SYN version V%d",*version);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit");file->ReadKeyword(":");
	tNbHit = file->ReadLLong();
	file->ReadKeyword("totalDes");file->ReadKeyword(":");
	tNbDesorption = file->ReadLLong();
	if (*version>=6) {
		file->ReadKeyword("no_scans");file->ReadKeyword(":");
		loaded_no_scans = file->ReadDouble();
	}
	else loaded_no_scans = 0;
	file->ReadKeyword("totalLeak");file->ReadKeyword(":");
	tNbLeak = file->ReadInt();
	if (*version>2) {
		file->ReadKeyword("totalFlux");file->ReadKeyword(":");
		tFlux=file->ReadDouble();
		file->ReadKeyword("totalPower");file->ReadKeyword(":");
		tPower=file->ReadDouble();
	}
	tNbAbsorption=0;
	distTraveledTotal=0.0;
	file->ReadKeyword("maxDes");file->ReadKeyword(":");
	tNbDesorptionMax = file->ReadLLong(); 
	file->ReadKeyword("nbVertex");file->ReadKeyword(":");
	sh.nbVertex = file->ReadInt();
	file->ReadKeyword("nbFacet");file->ReadKeyword(":");
	sh.nbFacet = file->ReadInt();
	file->ReadKeyword("nbSuper");file->ReadKeyword(":");
	sh.nbSuper = file->ReadInt();
	int nbF = 0; std::vector<std::vector<string>> loadFormulas;
	int nbV = 0;

	file->ReadKeyword("nbFormula"); file->ReadKeyword(":");
	nbF = file->ReadInt(); loadFormulas.reserve(nbF);
	file->ReadKeyword("nbView"); file->ReadKeyword(":");
	nbV = file->ReadInt();
	int nbS=0;
	file->ReadKeyword("nbSelection");file->ReadKeyword(":");
	nbS = file->ReadInt();

	if (*version>1) {
		file->ReadKeyword("nbRegions");file->ReadKeyword(":");
		int nbR = file->ReadInt();
		result=PARfileList(nbR);

		file->ReadKeyword("PARfiles");file->ReadKeyword("{");
		for (int i=0;i<nbR;i++) {
			char tmp[512];
			strcpy(tmp,file->ReadString());
			result.fileNames[i]=_strdup(tmp);
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("formulas");file->ReadKeyword("{");
	for(int i=0;i<nbF;i++) {
		char tmpName[256];
		char tmpExpr[512];
		strcpy(tmpName,file->ReadString());
		strcpy(tmpExpr,file->ReadString());
		//mApp->AddFormula(tmpName, tmpExpr); //parse after selection groups are loaded
		std::vector<string> newFormula;
		newFormula.push_back(tmpName);
		newFormula.push_back(tmpExpr);
		loadFormulas.push_back(newFormula);
	}
	file->ReadKeyword("}");
	
	file->ReadKeyword("views");file->ReadKeyword("{");
	for(int i=0;i<nbV;i++) {
		char tmpName[256];
		AVIEW v;
		strcpy(tmpName,file->ReadString());
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
		mApp->AddView(tmpName,v);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("selections");file->ReadKeyword("{");
	for(int i=0;i<nbS;i++) {
		ASELECTION s;
		char tmpName[256];
		strcpy(tmpName,file->ReadString());
		s.name=_strdup(tmpName);
		s.nbSel = file->ReadInt();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));

		for (int j=0;j<s.nbSel;j++) {
			s.selection[j]=file->ReadInt();
		}
		mApp->AddSelection(s.name,s);
	}
	file->ReadKeyword("}");

	for (int i = 0; i < nbF; i++) { //parse formulas now that selection groups are loaded
		mApp->AddFormula(loadFormulas[i][0].c_str(), loadFormulas[i][1].c_str());
	}

	file->ReadKeyword("structures");file->ReadKeyword("{");
	for(int i=0;i<sh.nbSuper;i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		sprintf(tmp,"%s.txt",strName[i]);
		strFileName[i] = _strdup(tmp);
	}
	file->ReadKeyword("}");

	// Allocate memory
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets,0,sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memset(vertices3,0,sh.nbVertex * sizeof(VERTEX3D));

	// Read vertices
	prg->SetMessage("Reading vertices...");
	file->ReadKeyword("vertices");file->ReadKeyword("{");
	for(int i=0;i<sh.nbVertex;i++) {
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1 ) throw Error(file->MakeError("Wrong vertex index !"));
		vertices3[i].x = file->ReadDouble();
		vertices3[i].y = file->ReadDouble();
		vertices3[i].z = file->ReadDouble();
		vertices3[i].selected = FALSE;
	}
	file->ReadKeyword("}");
	prg->SetMessage("Reading leaks and hits...");
	// Read leaks
	file->ReadKeyword("leaks");file->ReadKeyword("{");
	file->ReadKeyword("nbLeak");file->ReadKeyword(":");
	*nbleak = file->ReadInt();
	for (int i=0;i<*nbleak;i++) {
		int idx = file->ReadInt();
		if( idx != i ) throw Error(file->MakeError("Wrong leak index !"));
		(pleak+i)->pos.x = file->ReadDouble();
		(pleak+i)->pos.y = file->ReadDouble();
		(pleak+i)->pos.z = file->ReadDouble();

		(pleak+i)->dir.x = file->ReadDouble();
		(pleak+i)->dir.y = file->ReadDouble();
		(pleak+i)->dir.z = file->ReadDouble();   
	}
	file->ReadKeyword("}");

	// Read hit cache
	file->ReadKeyword("hits");file->ReadKeyword("{");
	file->ReadKeyword("nbHHit");file->ReadKeyword(":");
	*nbHHit = file->ReadInt();
	for (int i=0;i<*nbHHit;i++) {
		int idx = file->ReadInt();
		if( idx != i ) throw Error(file->MakeError("Wrong hit cache index !"));
		(pHits+i)->pos.x = file->ReadDouble();
		(pHits+i)->pos.y = file->ReadDouble();
		(pHits+i)->pos.z = file->ReadDouble();
		(pHits+i)->dF = file->ReadDouble();
		(pHits+i)->dP = file->ReadDouble();
		(pHits+i)->type = file->ReadInt();
	}
	file->ReadKeyword("}");
	// Read facets
	prg->SetMessage("Reading facets...");
	for(int i=0;i<sh.nbFacet;i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1 ) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nbI = file->ReadInt();
		if( nbI<3 ) {
			char errMsg[512];
			sprintf(errMsg,"Facet %d has only %d vertices. ",i,nbI);
			throw Error(errMsg);
		}
		prg->SetProgress((float)i/sh.nbFacet);
		facets[i] = new Facet(nbI);
		facets[i]->LoadSYN(file,*version,sh.nbVertex);
		file->ReadKeyword("}");
	}

	prg->SetMessage("Initalizing geometry and building mesh...");
	InitializeGeometry(); //Contains SetFacetTexture
	//AdjustProfile();
	isLoaded = TRUE;
	UpdateName(file);

	// Update mesh //InitializeGeometry will do it for us
	/*prg->SetMessage("Drawing textures...");
	for(int i=0;i<sh.nbFacet;i++) {
		double p = (double)i/(double)sh.nbFacet;
		prg->SetProgress(p);
		Facet *f = facets[i];
		//f->sh.texWidthD,f->sh.texHeightD,f->hasMesh);
		BuildFacetList(f);
		//double nU = Norme(&(f->sh.U));
		//f->tRatio = f->sh.texWidthD / nU;
	}*/
	return result;
}

void Geometry::LoadProfileSYN(FileReader *file,Dataport *dpHit) { //profiles and spectrums
	AccessDataport(dpHit);
	BYTE *buffer = (BYTE *)dpHit->buff;
	file->ReadKeyword("profiles");file->ReadKeyword("{");
	// Profiles
	int nbProfile;
	file->ReadKeyword("number");file->ReadKeyword(":");nbProfile=file->ReadInt();
	int *profileFacet = (int *)malloc((nbProfile)*sizeof(int));
	file->ReadKeyword("facets");file->ReadKeyword(":");
	for (int i=0;i<nbProfile;i++) 
		profileFacet[i]=file->ReadInt();
	for(int j=0;j<PROFILE_SIZE;j++) {
		for(int i=0;i<nbProfile;i++) {
			Facet *f = GetFacet(profileFacet[i]);
			llong *profilePtr_MC = (llong *)(buffer + f->sh.hitOffset + sizeof(SHHITS));
			double *profilePtr_flux = (double *)(buffer + f->sh.hitOffset + sizeof(SHHITS)+PROFILE_SIZE*sizeof(llong));
			double *profilePtr_power = (double *)(buffer + f->sh.hitOffset + sizeof(SHHITS)+PROFILE_SIZE*(sizeof(llong)+sizeof(double)));
			
			profilePtr_MC[j]=file->ReadLLong();
			profilePtr_flux[j]=file->ReadDouble();
			profilePtr_power[j]=file->ReadDouble();
		}
	}
	file->ReadKeyword("}");
	ReleaseDataport(dpHit);
	SAFE_FREE(profileFacet);
}
void Geometry::LoadSpectrumSYN(FileReader *file,Dataport *dpHit) { //spectrums and spectrums
	AccessDataport(dpHit);
	BYTE *buffer = (BYTE *)dpHit->buff;
	file->ReadKeyword("spectrums");file->ReadKeyword("{");
	// Spectrums
	int nbSpectrum;
	file->ReadKeyword("number");file->ReadKeyword(":");nbSpectrum=file->ReadInt();
	int *spectrumFacet = (int *)malloc((nbSpectrum)*sizeof(int));
	file->ReadKeyword("facets");file->ReadKeyword(":");
	for (int i=0;i<nbSpectrum;i++) 
		spectrumFacet[i]=file->ReadInt();
	
	for(int j=0;j<PROFILE_SIZE;j++) {
		for(int i=0;i<nbSpectrum;i++) {
			Facet *f = GetFacet(spectrumFacet[i]);
			int profileSize=(f->sh.isProfile)?PROFILE_SIZE*(2*sizeof(double)+sizeof(llong)):0;
			int textureSize = f->sh.texWidth*f->sh.texHeight*(2*sizeof(double)+sizeof(llong));
			int directionSize= f->sh.countDirection*f->sh.texWidth*f->sh.texHeight*sizeof(VHIT);
			double *shSpectrum_fluxwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize 
						+ textureSize + directionSize));
			double *shSpectrum_powerwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize 
						+ textureSize + directionSize + SPECTRUM_SIZE*sizeof(double)));
			
			shSpectrum_fluxwise[j]=file->ReadDouble();
			shSpectrum_powerwise[j]=file->ReadDouble();
		}
	}
	file->ReadKeyword("}");
	ReleaseDataport(dpHit);
	SAFE_FREE(spectrumFacet);
}
void Geometry::SaveProfileSYN(FileWriter *file,Dataport *dpHit,int super,BOOL saveSelected,BOOL crashSave) {
	//Profiles

	BYTE *buffer;
	if (!crashSave && !saveSelected) buffer = (BYTE *)dpHit->buff;
	file->Write("profiles {\n");
	// Profiles
	int nbProfile=0;
	int *profileFacet = (int *)malloc((sh.nbFacet)*sizeof(int));
	for(int i=0;i<sh.nbFacet;i++) 
		if ((!saveSelected && !crashSave) && facets[i]->sh.isProfile)
			profileFacet[nbProfile++] = i;
	file->Write(" number: ");file->WriteInt(nbProfile,"\n");
	file->Write(" facets: ");
	for (int i=0;i<nbProfile;i++)  //doesn't execute when crashSave or saveSelected...
		file->WriteInt(profileFacet[i],"\t");
	file->Write("\n");
	for(int j=0;j<PROFILE_SIZE;j++) {
		for(int i=0;i<nbProfile;i++) { //doesn't execute when crashSave or saveSelected...
			Facet *f = GetFacet(profileFacet[i]);
			llong *profilePtr_MC;
			if (!crashSave) profilePtr_MC= (llong *)(buffer + f->sh.hitOffset + sizeof(SHHITS));
			double *profilePtr_flux;
			if (!crashSave) profilePtr_flux = (double *)(buffer + f->sh.hitOffset + sizeof(SHHITS)+PROFILE_SIZE*sizeof(llong));
			double *profilePtr_power; 
			if (!crashSave) profilePtr_power = (double *)(buffer + f->sh.hitOffset + sizeof(SHHITS)+PROFILE_SIZE*(sizeof(llong)+sizeof(double)));

			file->WriteLLong((!crashSave)?profilePtr_MC[j]:0);file->Write("\t");
			file->WriteDouble((!crashSave)?profilePtr_flux[j]:0);file->Write("\t");
			file->WriteDouble((!crashSave)?profilePtr_power[j]:0);file->Write("\t");
		}
		file->Write("\n");
	}
	file->Write("}\n");
	SAFE_FREE(profileFacet);
}

void Geometry::SaveSpectrumSYN(FileWriter *file,Dataport *dpHit,int super,BOOL saveSelected,BOOL crashSave) {
	//Spectrums

	BYTE *buffer;
	if (!crashSave && !saveSelected) buffer = (BYTE *)dpHit->buff;
	file->Write("spectrums {\n");
	// Profiles
	int nbSpectrum=0;
	int *spectrumFacet = (int *)malloc((sh.nbFacet)*sizeof(int));
	for(int i=0;i<sh.nbFacet;i++) 
		if ((!saveSelected && !crashSave) && facets[i]->sh.hasSpectrum)
			spectrumFacet[nbSpectrum++] = i;
	file->Write(" number: ");file->WriteInt(nbSpectrum,"\n");
	file->Write(" facets: ");
	for (int i=0;i<nbSpectrum;i++)   //doesn't execute when crashSave or saveSelected...
		file->WriteInt(spectrumFacet[i],"\t");
	file->Write("\n");
	for(int j=0;j<SPECTRUM_SIZE;j++) {
		for(int i=0;i<nbSpectrum;i++) {  //doesn't execute when crashSave or saveSelected...
			Facet *f = GetFacet(spectrumFacet[i]);
			int profileSize=(f->sh.isProfile)?PROFILE_SIZE*(2*sizeof(double)+sizeof(llong)):0;
			int textureSize = f->sh.texWidth*f->sh.texHeight*(2*sizeof(double)+sizeof(llong));
			int directionSize= f->sh.countDirection*f->sh.texWidth*f->sh.texHeight*sizeof(VHIT);
			double *shSpectrum_fluxwise; 
			if (!crashSave) shSpectrum_fluxwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize 
						+ textureSize + directionSize));
			double *shSpectrum_powerwise; 
			if (!crashSave) shSpectrum_powerwise = (double *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + profileSize 
						+ textureSize + directionSize + SPECTRUM_SIZE*sizeof(double)));
			
			file->WriteDouble((!crashSave)?shSpectrum_fluxwise[j]:0);file->Write("\t");
			file->WriteDouble((!crashSave)?shSpectrum_powerwise[j]:0);file->Write("\t");
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
void Geometry::LoadXML_geom(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg, BOOL isSynxml){}
void Geometry::InsertXML(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg, BOOL newStr, BOOL isSynxml){}
BOOL Geometry::LoadXML_simustate(pugi::xml_node loadXML, Dataport *dpHit, Worker *work, GLProgress *progressDlg){ return FALSE; }