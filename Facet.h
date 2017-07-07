/*
  File:        Facet.h
  Description: Facet strucure
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

#ifndef FACETH
#define FACETH

#include "GLApp/GLApp.h"
#include "Shared.h"
#include "File.h"
#include "PugiXML/pugixml.hpp"
#include "Distributions.h"
//#include "Geometry.h"

class CellProperties;

struct NeighborFacet {
	size_t id;
	double angleDiff;
};

class CellProperties {
public:
	Vector2d* points;
	size_t nbPoints;
	float   area;     // Area of element
	float   uCenter;  // Center coordinates
	float   vCenter;  // Center coordinates
					  //int     elemId;   // Element index (MESH array)
					  //int full;
};

class Facet {

public:

  typedef struct {
  
	size_t nbV;
	size_t nbF;
    Facet **facets;

  } FACETGROUP;

  typedef struct {
	  size_t u;
	  size_t v;
	  size_t width;
	  size_t height;
  } BOX;

  // Constructor/Desctructor/Initialisation
  Facet(int nbIndex);
  ~Facet();

  // Shared struct
  SHFACET sh;

  size_t      *indices;      // Indices (Reference to geometry vertex)
  Vector2d *vertices2;    // Vertices (2D plane space, UV coordinates)
  int     *cellPropertiesIds;      // -1 if full element, -2 if outside polygon, otherwise index in meshvector
  CellProperties* meshvector;
  size_t meshvectorsize;

  // Normalized plane equation (ax + by + cz + d = 0)
  double a;
  double b;
  double c;
  double d;
  double err;          // planeity error
  size_t texDimH;         // Texture dimension (a power of 2)
  size_t texDimW;         // Texture dimension (a power of 2)
  double tRatio;       // Texture sample per unit
  bool	textureVisible; //Draw the texture?
  bool  collinear;      //All vertices are on a line (non-simple)
  bool	volumeVisible;	//Draw volume?
  bool    hasMesh;     // Has texture
  VHIT   *dirCache;    // Direction field cache
  bool textureError;   // Disable rendering if the texture has an error

  // GUI stuff
  bool  *visible;         // Edge visible flag
  bool   selected;        // Selected flag
  BOX    selectedElem;    // Selected mesh element
  GLint  glElem;          // Surface elements boundaries
  GLint  glSelElem;       // Selected surface elements boundaries
  GLint  glList;          // Geometry with texture
  GLuint glTex;           // Handle to OpenGL texture

  //Smart selection
  std::vector<NeighborFacet> neighbors;

  // Global hit counters
  SHHITS counterCache;

  //Facet methods

  bool  IsTXTLinkFacet();
  Vector3d GetRealCenter();
  void  LoadTXT(FileReader *file);
  void  SaveTXT(FileWriter *file);
  void  LoadGEO(FileReader *file,int version,int nbVertex);
  //void  SaveGEO(FileWriter *file,int idx);
  void  LoadSYN(FileReader *file, const std::vector<Material> &materials, int version, int nbVertex);
  void  SaveSYN(FileWriter *file, const std::vector<Material> &materials, int idx,bool crashSave=false);
  void  LoadXML(pugi::xml_node f, int nbVertex, bool isMolflowFile, int vertexOffset);
  bool  IsCoplanarAndEqual(Facet *f,double threshold);
  size_t   GetIndex(int idx);
  size_t   GetIndex(size_t idx);
  void  CopyFacetProperties(Facet *f,bool copyMesh=false);
  void  SwapNormal();
  void  Explode(FACETGROUP *group);
  void  FillVertexArray(InterfaceVertex *v);
  void  BuildMeshList();
  void  InitVisibleEdge();
  bool  SetTexture(double width,double height,bool useMesh);
  size_t GetGeometrySize();
  size_t GetHitsSize();
  size_t GetTexSwapSize(bool useColormap);
  size_t GetTexRamSize();
  size_t GetTexSwapSizeForRatio(double ratio,bool useColor);
  size_t GetTexRamSizeForRatio(double ratio,bool useMesh,bool countDir);
  size_t GetNbCellForRatio(double ratio);
  size_t GetNbCell();
  void  UpdateFlags();
  void  BuildTexture(double *texBuffer,double min,double max,double no_scans,bool useColorMap,bool doLog,bool normalize=true);
  void  BuildTexture(llong *texBuffer,llong min,llong max,bool useColorMap,bool doLog);
  bool  BuildMesh();
  void  BuildSelElemList();
  int   RestoreDeviceObjects();
  int   InvalidateDeviceObjects();
  void  DetectOrientation();
  double GetSmooth(const int &i,const int &j,double *texBuffer,const float &scaleF);
  double GetSmooth(const int &i,const int &j,llong *texBuffer,const float &scaleF);
  void  glVertex2u(double u,double v);
  void  ShiftVertex();
  void  RenderSelectedElem();
  void  SelectElem(size_t u, size_t v, size_t width, size_t height);
  void  UnselectElem();
  float GetMeshArea(size_t index, bool correct2sides = false);
  size_t GetMeshNbPoint(size_t index);
  Vector2d GetMeshPoint(size_t index, size_t pointId);
  Vector2d GetMeshCenter(size_t index);
  double GetArea();
};

class DeletedFacet {
public:
	Facet *f;
	size_t ori_pos;
	bool replaceOri;
};

#endif /* FACETH */
