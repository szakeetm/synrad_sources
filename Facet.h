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

class Facet {

public:

  typedef struct {
  
    int nbV;
    int nbF;
    Facet **facets;

  } FACETGROUP;

  typedef struct {
    int u;
    int v;
    int width;
    int height;
  } BOX;

  // Constructor/Desctructor/Initialisation
  Facet(int nbIndex);
  ~Facet();

  // Shared struct
  SHFACET sh;

  int      *indices;      // Indices (Reference to geometry vertex)
  VERTEX2D *vertices2;    // Vertices (2D plane space, UV coordinates)
  MESH     *meshPts;      // Mesh poly
  size_t       nbElem;       // Number of mesh elem

  // Normalized plane equation (ax + by + cz + d = 0)
  double a;
  double b;
  double c;
  double d;
  double err;          // planeity error
  int texDimH;         // Texture dimension (a power of 2)
  int texDimW;         // Texture dimension (a power of 2)
  double tRatio;       // Texture sample per unit
  BOOL	textureVisible; //Draw the texture?
  BOOL  collinear;      //All vertices are on a line (non-simple)
  BOOL	volumeVisible;	//Draw volume?
  SHELEM *mesh;        // Element mesh
  BOOL    hasMesh;     // Has texture
  VHIT   *dirCache;    // Direction field cache
  BOOL textureError;   // Disable rendering if the texture has an error

  // GUI stuff
  BOOL  *visible;         // Edge visible flag
  BOOL   selected;        // Selected flag
  BOX    selectedElem;    // Selected mesh element
  GLint  glElem;          // Surface elements boundaries
  GLint  glSelElem;       // Selected surface elements boundaries
  GLint  glList;          // Geometry with texture
  GLuint glTex;           // Handle to OpenGL texture

  //Facet methods

  BOOL  IsLinkFacet();
  void  LoadTXT(FileReader *file);
  void  SaveTXT(FileWriter *file);
  void  LoadGEO(FileReader *file,int version,int nbVertex);
  //void  SaveGEO(FileWriter *file,int idx);
  void  LoadSYN(FileReader *file,int version,int nbVertex);
  void  SaveSYN(FileWriter *file,int idx,BOOL crashSave=FALSE);
  BOOL  IsCoplanar(Facet *f,double threshold);
  int   GetIndex(int idx);
  void  Copy(Facet *f,BOOL copyMesh=FALSE);
  void  SwapNormal();
  void  Explode(FACETGROUP *group);
  void  FillVertexArray(VERTEX3D *v);
  void  BuildMeshList();
  void  InitVisibleEdge();
  void  SetTexture(double width,double height,BOOL useMesh);
  DWORD GetGeometrySize();
  DWORD GetHitsSize();
  DWORD GetTexSwapSize(BOOL useColormap);
  DWORD GetTexRamSize();
  DWORD GetTexSwapSizeForRatio(double ratio,BOOL useColor);
  DWORD GetTexRamSizeForRatio(double ratio,BOOL useMesh,BOOL countDir);
  DWORD GetNbCellForRatio(double ratio);
  DWORD GetNbCell();
  void  UpdateFlags();
  void  BuildTexture(double *texBuffer,double min,double max,double no_scans,BOOL useColorMap,BOOL doLog,BOOL normalize=TRUE);
  void  BuildTexture(llong *texBuffer,llong min,llong max,BOOL useColorMap,BOOL doLog);
  void  BuildMesh();
  void  BuildSelElemList();
  int   RestoreDeviceObjects();
  int   InvalidateDeviceObjects();
  void  DetectOrientation();
  double GetSmooth(const int &i,const int &j,double *texBuffer,const float &scaleF);
  double GetSmooth(const int &i,const int &j,llong *texBuffer,const float &scaleF);
  void  glVertex2u(double u,double v);
  void  ShiftVertex();
  void  RenderSelectedElem();
  void  SelectElem(int u,int v,int width,int height);
  void  UnselectElem();

};

#endif /* FACETH */
