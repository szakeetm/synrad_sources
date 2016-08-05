/*
  File:        GeometryRender.cpp
  Description: Geometry class (OpenGL rendering/selection stuff)
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
#include "Worker.h"
#include "Utils.h"
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "GLApp/GLMatrix.h"
#include "SynRad.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"

extern SynRad *mApp;

// OpenGL stuff ----------------------------------------------

// -----------------------------------------------------------

void Geometry::Select(Facet *f) {
  f->selected = (viewStruct==-1) || (viewStruct==f->sh.superIdx);
  if( !f->selected ) f->UnselectElem();
}

// -----------------------------------------------------------

void Geometry::Select(int facet) {
  if(!isLoaded) return;
  Select(facets[facet]);
  nbSelectedHist = 0;
  AddToSelectionHist(facet);
}

// -----------------------------------------------------------

void Geometry::SelectArea(int x1,int y1,int x2,int y2,BOOL clear,BOOL unselect,BOOL vertexBound,BOOL circularSelection) {

  // Select a set of facet according to a 2D bounding rectangle
  // (x1,y1) and (x2,y2) are in viewport coordinates
  float rx,ry,rz,rw,r2;
  int _x1,_y1,_x2,_y2;

  _x1 = MIN(x1,x2);
  _x2 = MAX(x1,x2);
  _y1 = MIN(y1,y2);
  _y2 = MAX(y1,y2);

  if (circularSelection) {
	  r2=pow((float)(x1-x2),2)+pow((float)(y1-y2),2);
  }

  GLfloat mProj[16];
  GLfloat mView[16];
  GLVIEWPORT g;

  glGetFloatv(GL_PROJECTION_MATRIX , mProj);
  glGetFloatv(GL_MODELVIEW_MATRIX , mView);
  glGetIntegerv(GL_VIEWPORT,(GLint *)&g);

  GLMatrix proj; proj.LoadGL(mProj);
  GLMatrix view; view.LoadGL(mView);
  GLMatrix m; m.Multiply(&proj,&view);

  if(clear && !unselect) Unselect();
  nbSelectedHist = 0;
  int lastPaintedProgress=-1;
  char tmp[256];
  int paintStep=(int)((double)sh.nbFacet/10.0);

  for (int i=0;i<sh.nbFacet;i++) {
	  if (sh.nbFacet>5000) {
		  if((i-lastPaintedProgress)>paintStep) {
			  lastPaintedProgress=i;;
			  sprintf(tmp,"Facet search: %d%%",(int)(i*100.0/(double)sh.nbFacet));
			  mApp->SetFacetSearchPrg(TRUE,tmp);
		  }
	  }
    Facet *f = facets[i];
    if(viewStruct==-1 || f->sh.superIdx==viewStruct) {

      int nb = facets[i]->sh.nbIndex;
      BOOL isInside = TRUE;
      int j=0;
	  BOOL hasSelectedVertex = FALSE;
      while( j<nb && isInside ) {

        int idx = f->indices[j];
        m.TransfomVec((float)vertices3[idx].x,(float)vertices3[idx].y,(float)vertices3[idx].z,1.0f,
                      &rx,&ry,&rz,&rw);

        if( rw>0.0f ) {
          int xe = (int)(((rx / rw) + 1.0f) * (float)g.width/2.0f);
          int ye = (int)(((-ry / rw) + 1.0f) * (float)g.height/2.0f);
		  if (!circularSelection)
			  isInside = (xe>=_x1) && (xe<=_x2) && (ye>=_y1) && (ye<=_y2);
		  else //circular selection
			  isInside= (pow((float)(xe-x1),2)+pow((float)(ye-y1),2))<=r2;
		  if(vertices3[idx].selected) hasSelectedVertex=TRUE;
        } else {
          isInside = FALSE;
        }
        j++;

      }

      if( isInside &&(!vertexBound || hasSelectedVertex)) {
        if( !unselect ) {
          f->selected = !unselect;
        } else {
          f->selected = !unselect;
        }
      }

    }    
  }
  mApp->SetFacetSearchPrg(FALSE,NULL);
  UpdateSelection();
}

// -----------------------------------------------------------

void Geometry::Select(int x,int y,BOOL clear,BOOL unselect,BOOL vertexBound,int width,int height) {
  int i;
  if(!isLoaded) return;

  // Select a facet on a mouse click in 3D perspectivce view 
  // (x,y) are in screen coordinates
  // TODO: Handle clipped polygon

  // Check intersection of the facet and a "perspective ray"
  int *allXe = (int *)malloc(sh.nbVertex * sizeof(int));
  int *allYe = (int *)malloc(sh.nbVertex * sizeof(int));

  // Transform points to screen coordinates
  BOOL *ok = (BOOL *)malloc(sh.nbVertex*sizeof(BOOL));
  BOOL *onScreen = (BOOL *)malloc(sh.nbVertex*sizeof(BOOL));
  for(i=0;i<sh.nbVertex;i++) {//here we could speed up by choosing visible vertices only?
    ok[i] = GLToolkit::Get2DScreenCoord((float)vertices3[i].x,(float)vertices3[i].y,(float)vertices3[i].z,allXe+i,allYe+i);
    onScreen[i]=(ok[i]&& (*(allXe+i)>=0) && (*(allYe+i)>=0) && (*(allXe+i)<=width) && (*(allYe+i)<=height));
  }

  // Check facets
  BOOL found = FALSE;
  BOOL clipped;
  BOOL hasVertexOnScreen;
  BOOL hasSelectedVertex;
  i=0;
  char tmp[256];
  int lastFound=-1;
  int lastPaintedProgress=-1;
  int paintStep=(int)((double)sh.nbFacet/10.0);

  while(i<sh.nbFacet && !found) {
	  if (sh.nbFacet>5000) {
		  if((i-lastPaintedProgress)>paintStep) {
			  lastPaintedProgress=i;;
			  sprintf(tmp,"Facet search: %d%%",(int)(i*100.0/(double)sh.nbFacet));
			  mApp->SetFacetSearchPrg(TRUE,tmp);
		  }
	  }
    if(viewStruct==-1 || facets[i]->sh.superIdx==viewStruct) {

      clipped = FALSE;
	  hasVertexOnScreen = FALSE;
	  hasSelectedVertex = FALSE;
      int nb = facets[i]->sh.nbIndex;
      // Build array of 2D points
      int *xe = (int *)malloc(nb * sizeof(int));
      int *ye = (int *)malloc(nb * sizeof(int));
      for(int j=0;j<nb && !clipped;j++) {
        int idx = facets[i]->indices[j];
        if( ok[idx] ) {
          xe[j] = allXe[idx];
          ye[j] = allYe[idx];
		  if (onScreen[idx]) hasVertexOnScreen = TRUE;
		} else {
			clipped = TRUE;
		}
      }
	  if (vertexBound) { //CAPS LOCK on, select facets onyl with at least one seleted vertex
		  for(int j=0;j<nb&&(!hasSelectedVertex);j++) {
			  int idx = facets[i]->indices[j];
			  if (vertices3[idx].selected) hasSelectedVertex=TRUE;
		  }
	  }

      if( !clipped && hasVertexOnScreen&&(!vertexBound||hasSelectedVertex)) {

        found = GLToolkit::IsInsidePoly(x,y,xe,ye,nb);
        free(xe);
        free(ye);

        if(found) {
          if( unselect ) { 
            facets[i]->selected = FALSE;
            found=FALSE; 
          }
          if( AlreadySelected(i) ) {
            // Already selected, try to find an other facet
            lastFound=i;
            found = FALSE;
          }
        }
      }

    }

    if(!found) i++;

  }
  mApp->SetFacetSearchPrg(FALSE,NULL);
  if(clear && !unselect) Unselect();

  if(!found && lastFound>=0) {
    if( !unselect ) {
      // Restart
      nbSelectedHist = 0;
      AddToSelectionHist(lastFound);
    }
    facets[lastFound]->selected = !unselect;
	if (!unselect) mApp->facetList->ScrollToVisible(lastFound,0,TRUE); //scroll to selected facet
  } else {
    if(found) {
      if( !unselect ) AddToSelectionHist(i);
      facets[i]->selected = !unselect;
	if (!unselect) mApp->facetList->ScrollToVisible(i,0,TRUE); //scroll to selected facet
	} else {
      nbSelectedHist = 0;
    }
  }

  free(allXe);
  free(allYe);
  free(ok);
  free(onScreen);
  UpdateSelection();

}

// ----------------------------------------------------------

void Geometry::SelectVertex(int vertexId) {
  //isVertexSelected[vertexId] = (viewStruct==-1) || (viewStruct==f->sh.superIdx);
	//here we should look through facets if vertex is member of any
  //if( !f->selected ) f->UnselectElem();
	if(!isLoaded) return;
	vertices3[vertexId].selected = TRUE;
	//nbSelectedHistVertex = 0;
	//AddToSelectionHistVertex(vertexId);
}


// -----------------------------------------------------------

void Geometry::SelectVertex(int x1,int y1,int x2,int y2,BOOL shiftDown,BOOL ctrlDown,BOOL circularSelection) {

  // Select a set of vertices according to a 2D bounding rectangle
  // (x1,y1) and (x2,y2) are in viewport coordinates

  float rx,ry,rz,rw,r2;
  int _x1,_y1,_x2,_y2;
  
  _x1 = MIN(x1,x2);
  _x2 = MAX(x1,x2);
  _y1 = MIN(y1,y2);
  _y2 = MAX(y1,y2);

  if (circularSelection) {
	  r2=pow((float)(x1-x2),2)+pow((float)(y1-y2),2);
  }

  GLfloat mProj[16];
  GLfloat mView[16];
  GLVIEWPORT g;

  glGetFloatv(GL_PROJECTION_MATRIX , mProj);
  glGetFloatv(GL_MODELVIEW_MATRIX , mView);
  glGetIntegerv(GL_VIEWPORT,(GLint *)&g);

  GLMatrix proj; proj.LoadGL(mProj);
  GLMatrix view; view.LoadGL(mView);
  GLMatrix m; m.Multiply(&proj,&view);

  
  if(!ctrlDown && !shiftDown) {
	  UnselectAllVertex();EmptySelectedVertexList();
		//nbSelectedHistVertex = 0;
  }
  

  for(int i=0;i<sh.nbVertex;i++) {

    VERTEX3D *v = GetVertex(i);
    //if(viewStruct==-1 || f->sh.superIdx==viewStruct) {
	if (true) {

      
		BOOL isInside;
        int idx = i;
        m.TransfomVec((float)vertices3[idx].x,(float)vertices3[idx].y,(float)vertices3[idx].z,1.0f,
                      &rx,&ry,&rz,&rw);

        if( rw>0.0f ) {
          int xe = (int)(((rx / rw) + 1.0f) * (float)g.width/2.0f);
          int ye = (int)(((-ry / rw) + 1.0f) * (float)g.height/2.0f);
		  if (!circularSelection)
			  isInside = (xe>=_x1) && (xe<=_x2) && (ye>=_y1) && (ye<=_y2);
		  else //circular selection
			  isInside= (pow((float)(xe-x1),2)+pow((float)(ye-y1),2))<=r2;
        } else {
          isInside = FALSE;
        }
        

      

      if( isInside ) {
		  vertices3[i].selected = !ctrlDown;
		  if (ctrlDown) RemoveFromSelectedVertexList(i);
		  else {
			  AddToSelectedVertexList(i);
			  if (mApp->facetCoordinates) mApp->facetCoordinates->UpdateId(i);
		  }
	  }

    }    
  }
  
  //UpdateSelectionVertex();
  if (mApp->vertexCoordinates) mApp->vertexCoordinates->Update();
}

void Geometry::SelectVertex(int x,int y,BOOL shiftDown,BOOL ctrlDown) {
  int i;
  if(!isLoaded) return;

  // Select a vertex on a mouse click in 3D perspectivce view 
  // (x,y) are in screen coordinates
  // TODO: Handle clipped polygon

  // Check intersection of the facet and a "perspective ray"
  int *allXe = (int *)malloc(sh.nbVertex * sizeof(int));
  int *allYe = (int *)malloc(sh.nbVertex * sizeof(int));

  // Transform points to screen coordinates
  BOOL *ok = (BOOL *)malloc(sh.nbVertex*sizeof(BOOL));
  for(i=0;i<sh.nbVertex;i++)
    ok[i] = GLToolkit::Get2DScreenCoord((float)vertices3[i].x,(float)vertices3[i].y,(float)vertices3[i].z,allXe+i,allYe+i);

  //Get Closest Point to click
  double minDist=9999;
  double distance;
  int minId=-1;
  for(i=0;i<sh.nbVertex;i++) {
	  if (ok[i] && !(allXe[i]<0) && !(allYe[i]<0)) { //calculate only for points on screen
		  distance=pow((double)(allXe[i]-x),2)+pow((double)(allYe[i]-y),2);
		  if (distance<minDist) {
			  minDist=distance;
			  minId=i;
		  }
	  }
  }

    if(!ctrlDown && !shiftDown) {
	  UnselectAllVertex();EmptySelectedVertexList();
		//nbSelectedHistVertex = 0;
  }

  if (minDist<250.0) {
	  vertices3[minId].selected=!ctrlDown;
	  if (ctrlDown) RemoveFromSelectedVertexList(minId);
	  else {
		  AddToSelectedVertexList(minId);
		  if (mApp->facetCoordinates) mApp->facetCoordinates->UpdateId(minId);
	  //nbSelectedHistVertex++;
	  }
  }
  
  free(allXe);
  free(allYe);
  free(ok);
  //UpdateSelection();
  if (mApp->vertexCoordinates) mApp->vertexCoordinates->Update();
}

void Geometry::AddToSelectionHist(int f) {

  if( nbSelectedHist<SEL_HISTORY ) {
    selectHist[nbSelectedHist] = f;
    nbSelectedHist++;
  }

}


BOOL Geometry::AlreadySelected(int f) {
  
  // Check if the facet has already been selected
  BOOL found = FALSE;
  int i = 0;
  while(!found && i<nbSelectedHist) {
    found = (selectHist[i] == f);
    if(!found) i++;
  }
  return found;

}



void Geometry::SelectAll() {
  for(int i=0;i<sh.nbFacet;i++)
    Select(facets[i]);
  UpdateSelection();
}



void Geometry::UnSelectAll() {
  for(int i=0;i<sh.nbFacet;i++)
    facets[i]->selected = FALSE;
  UpdateSelection();
}



int Geometry::GetNbSelected() {
  return nbSelected;
}

void Geometry::AddToSelectionHistVertex(int f) {

  if( nbSelectedHistVertex<SEL_HISTORY ) {
    selectHistVertex[nbSelectedHistVertex] = f;
    nbSelectedHistVertex++;
  }

}



BOOL Geometry::AlreadySelectedVertex(int idx) {
  

  // Check if the vertex is in the selection history
  BOOL found = FALSE;
  int i = 0;
  while(!found && i<nbSelectedHistVertex) {
    found = (selectHistVertex[i] == idx);
    if(!found) i++;
  }
  return found;

}

void Geometry::EmptySelectedVertexList(){
	selectedVertexList = std::vector<int>();
}

void Geometry::RemoveFromSelectedVertexList(int vertexId){
	for (size_t j = 0; j < selectedVertexList.size(); j++)
		if (selectedVertexList[j] == vertexId)
			selectedVertexList.erase(selectedVertexList.begin() + j);
}

void Geometry::AddToSelectedVertexList(int vertexId){
	selectedVertexList.push_back(vertexId);
}


void Geometry::SelectAllVertex() {
	for (int i = 0; i < sh.nbVertex; i++)
		SelectVertex(i);
	//UpdateSelectionVertex();
}





int Geometry::GetNbSelectedVertex() {
  nbSelectedVertex=0;
  for (int i=0;i<sh.nbVertex;i++) {
	  if (vertices3[i].selected) nbSelectedVertex++;
  }
  return nbSelectedVertex;
}


void Geometry::Unselect() {
  for(int i=0;i<sh.nbFacet;i++) {
    facets[i]->selected = FALSE;
    facets[i]->UnselectElem();
  }
  UpdateSelection();
}



void Geometry::UnselectAllVertex() {
  for(int i=0;i<sh.nbVertex;i++) {
    vertices3[i].selected = FALSE;
    //facets[i]->UnselectElem(); //what is this?
  }
  //UpdateSelectionVertex();
}



void Geometry::DrawFacet(Facet *f,BOOL offset,BOOL showHidden,BOOL selOffset) {

  // Render a facet (wireframe)
  int nb = f->sh.nbIndex;
  int i1;

  if( offset ) {

    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glEnable(GL_POLYGON_OFFSET_LINE);
	if (selOffset) {
		glPolygonOffset(0.0f,1.0f);
	} else {
		glPolygonOffset(0.0f,5.0f);
	}
	glBegin(GL_POLYGON);
    for(int j=0;j<nb;j++) {
      i1 = f->indices[j];
      glEdgeFlag(f->visible[j] || showHidden);
      glVertex3d(vertices3[i1].x,vertices3[i1].y,vertices3[i1].z);
    }
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glDisable(GL_POLYGON_OFFSET_LINE);

  } else {

    if( nb<8 ) {
      // No hole possible
      glBegin(GL_LINE_LOOP);
      for(int j=0;j<nb;j++) {
        i1 = f->indices[j];
        glVertex3d( vertices3[i1].x , vertices3[i1].y , vertices3[i1].z );
      }
      glEnd();
    } else {
      glBegin(GL_LINES);
      int i1,i2,j;
      for(j=0;j<nb-1;j++) {
        if( f->visible[j] || showHidden ) {
          i1 = f->indices[j];
          i2 = f->indices[j+1];
          glVertex3d( vertices3[i1].x , vertices3[i1].y , vertices3[i1].z);
          glVertex3d( vertices3[i2].x , vertices3[i2].y , vertices3[i2].z);
        }
      }
      // Last segment
      if(f->visible[j] || showHidden ) {
        i1 = f->indices[j];
        i2 = f->indices[0];
        glVertex3d( vertices3[i1].x , vertices3[i1].y , vertices3[i1].z);
        glVertex3d( vertices3[i2].x , vertices3[i2].y , vertices3[i2].z);
      }
      glEnd();
    }

  }

}



void Geometry::PaintSelectedVertices(BOOL hiddenVertex) {

  nbSelectedVertex = 0;

  char tmp[256];
  int *vIdx = (int *)malloc(sh.nbVertex*sizeof(int));
  memset(vIdx,0xFF,sh.nbVertex*sizeof(int));
  for(int i=0;i<sh.nbVertex;i++ ) {
    //VERTEX3D *v = GetVertex(i);
    if( vertices3[i].selected ) {
		vIdx[nbSelectedVertex] = i;
		nbSelectedVertex++;
	}
   }

    // Draw dot
    if (!mApp->whiteBg) glPointSize(6.0f);
	else glPointSize(7.0f);
	glEnable(GL_POINT_SMOOTH);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
	if (hiddenVertex) glDisable(GL_DEPTH_TEST);
	else glEnable(GL_DEPTH_TEST);
	if (!mApp->whiteBg) glColor3f(1.0f, 0.9f, 0.2f);
	else glColor3f(1.0f,0.5f,0.2f);

    glBegin(GL_POINTS);
    for(int i=0;i<nbSelectedVertex;i++) {
        VERTEX3D *v = GetVertex(vIdx[i]);
        glVertex3d(v->x,v->y,v->z);
    }
    glEnd();

    // Save contexct
    GLToolkit::DrawStringInit();
	if (!mApp->whiteBg) GLToolkit::GetDialogFont()->SetTextColor(1.0f, 0.9f, 0.2f);
	else GLToolkit::GetDialogFont()->SetTextColor(1.0f,0.2f,0.0f);

    // Draw Labels
	glEnable(GL_BLEND);
    for(int i=0;i<nbSelectedVertex;i++) {
		sprintf(tmp,"%d ",vIdx[i]+1);
        VERTEX3D *v = GetVertex(vIdx[i]);
        GLToolkit::DrawString((float)v->x,(float)v->y,(float)v->z,tmp,GLToolkit::GetDialogFont(),2,2);
	  
      }
    glDisable(GL_BLEND);
    //Restore
    GLToolkit::DrawStringRestore();
    free(vIdx);
	if (hiddenVertex) glEnable(GL_DEPTH_TEST);

}



void Geometry::DrawPolys() {

  int *f3 = (int *)malloc(sh.nbFacet*sizeof(int));
  int *f4 = (int *)malloc(sh.nbFacet*sizeof(int));
  int *fp = (int *)malloc(sh.nbFacet*sizeof(int));
  int nbF3=0;
  int nbF4=0;
  int nbFP=0;

  // Group TRI,QUAD and POLY
  for(int i=0;i<sh.nbFacet;i++) {
    int nb = facets[i]->sh.nbIndex;
	if( facets[i]->volumeVisible) {
		if( nb==3 ) {
			f3[nbF3++]=i;
		} else if (nb==4) {
			f4[nbF4++]=i;
		} else {
			fp[nbFP++]=i;
	 }
	}
  }

  // Draw
  glBegin(GL_TRIANGLES);

  // Triangle
  for(int i=0;i<nbF3;i++)
    FillFacet(facets[f3[i]],FALSE);

  // Triangulate polygon
  for(int i=0;i<nbFP;i++)
    Triangulate(facets[fp[i]],FALSE);

  glEnd();

  // Quads
  glBegin(GL_QUADS);
  for(int i=0;i<nbF4;i++)
    FillFacet(facets[f4[i]],FALSE);
  glEnd();

  free(f3);
  free(f4);
  free(fp);

}



void Geometry::SetCullMode(int mode) {

  switch(mode) {
    case 1: // SHOW_FRONT
      glEnable(GL_CULL_FACE);
      glCullFace(GL_FRONT);
      break;
    case 2: // SHOW_BACK
      glEnable(GL_CULL_FACE);
      glCullFace(GL_BACK);
      break;
    default: //SHOW_FRONTANDBACK
      glDisable(GL_CULL_FACE);
  }

}



void Geometry::BuildTexture(BYTE *hits) {

  SHGHITS *shGHit = (SHGHITS *)hits;
  //float dCoef = 1.0f;

      //if( shGHit->total.nbDesorbed>0 ) {
      //  dCoef = (float)totalOutgassing / (float)shGHit->total.nbDesorbed;
	Worker *worker=&(mApp->worker);     
	texCMin_MC = shGHit->minHit_MC;// * dCoef;
	texCMax_MC = shGHit->maxHit_MC;// * dCoef;
	texCMin_flux = shGHit->minHit_flux/worker->no_scans;// * dCoef;
	texCMax_flux = shGHit->maxHit_flux/worker->no_scans;// * dCoef;
	texCMin_power = shGHit->minHit_power/worker->no_scans;// * dCoef;
	texCMax_power = shGHit->maxHit_power/worker->no_scans;// * dCoef;

      //} else {
      //  texCMin = shGHit->minHit;
      //  texCMax = shGHit->maxHit;
      //}

  

  double iDesorbed = 0.0;
  if( shGHit->total.nbDesorbed ) 
    iDesorbed = 1.0 / (double)shGHit->total.nbDesorbed;

  GLProgress *prg = new GLProgress("Building texture", "Frame update");
  prg->SetBounds(5, 28, 300, 90);
  prg->SetVisible(TRUE);
  for(int i=0;i<sh.nbFacet;i++) {
	  prg->SetProgress((double)i / (double)sh.nbFacet);
	  Facet *f = facets[i];
	  GLint max_t;
	  glGetIntegerv(GL_MAX_TEXTURE_SIZE,&max_t);
	  if (f->sh.texHeight>max_t || f->sh.texWidth>max_t) {
		  if (!f->textureError) {
			  char tmp[1024];
			  sprintf(tmp,"Facet #%d has a texture of %dx%d cells.\n"
				  "Your video card only supports texture dimensions (width or height) up to %d cells.\n"
				  "Texture rendering has been disabled on this facet, but you can still read texture values\n"
				  "using the Texture Plotter window. Consider using a smaller mesh resolution, or split the facet\n"
				  "into smaller parts. (Use Facet/Explode... command)",i+1,f->sh.texHeight,f->sh.texWidth,max_t);
			  GLMessageBox::Display(tmp,"OpenGL Error",GLDLG_OK,GLDLG_ICONWARNING);
		  }
		  f->textureError = TRUE;
		  return;
	  } else {
		  f->textureError=FALSE;
	  }

	  int profSize = (f->sh.isProfile)?(PROFILE_SIZE*(sizeof(llong)+2*sizeof(double))):0;
	  //int spectrumSize = (f->sh.hasSpectrum)?(SPECTRUM_SIZE*2*sizeof(double)):0;
    int nbElem = f->sh.texWidth*f->sh.texHeight;
    int tSize = nbElem*sizeof(double);
    int dSize = nbElem*sizeof(VHIT);
	
	if( f->sh.isTextured ) {
       // Retrieve texture from shared memory (every seconds)
       llong *hits_MC = (llong *)((BYTE *)shGHit + (f->sh.hitOffset + sizeof(SHHITS) + profSize));
	   double *hits_flux = (double *)((BYTE *)shGHit + (f->sh.hitOffset + sizeof(SHHITS) + profSize+nbElem*sizeof(llong)));
	   double *hits_power = (double *)((BYTE *)shGHit + (f->sh.hitOffset + sizeof(SHHITS) + profSize+nbElem*(sizeof(llong)+sizeof(double))));

       if( texAutoScale ) {
		   if (textureMode==TEXTURE_MODE_MCHITS) f->BuildTexture(hits_MC,texCMin_MC,texCMax_MC,texColormap,texLogScale);
		   else if (textureMode==TEXTURE_MODE_FLUX) f->BuildTexture(hits_flux,texCMin_flux,texCMax_flux,worker->no_scans,texColormap,texLogScale);
		   else if (textureMode==TEXTURE_MODE_POWER) f->BuildTexture(hits_power,texCMin_power,texCMax_power,worker->no_scans,texColormap,texLogScale);
	   } else {
		   if (textureMode==TEXTURE_MODE_MCHITS) f->BuildTexture(hits_MC,texMin_MC,texMax_MC,texColormap,texLogScale);
		   else if (textureMode==TEXTURE_MODE_FLUX) f->BuildTexture(hits_flux,texMin_flux,texMax_flux,worker->no_scans,texColormap,texLogScale);
		   else if (textureMode==TEXTURE_MODE_POWER) f->BuildTexture(hits_power,texMin_power,texMax_power,worker->no_scans,texColormap,texLogScale);
	   }
    }
    if( f->sh.countDirection && f->dirCache ) {
       VHIT *dirs = (VHIT *)((BYTE *)shGHit + (f->sh.hitOffset + sizeof(SHHITS) + profSize + tSize));
       for(int j=0;j<nbElem;j++) {
         f->dirCache[j].dir.x = dirs[j].dir.x * iDesorbed; 
         f->dirCache[j].dir.y = dirs[j].dir.y * iDesorbed; 
         f->dirCache[j].dir.z = dirs[j].dir.z * iDesorbed; 
         f->dirCache[j].count = dirs[j].count;
       }
    }

  }
  prg->SetVisible(FALSE);
  SAFE_DELETE(prg);
}



void Geometry::Render(GLfloat *matView,BOOL renderVolume,BOOL renderTexture,int showMode,BOOL filter,BOOL showHidden,BOOL showMesh,BOOL showDir) {

  if(!isLoaded) return;

  // Render the geometry
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_BLEND);


  // Render Volume
  if( renderVolume ) {
	  glPolygonOffset(1.0f,4.0f);
	  SetCullMode(showMode);
	  GLfloat global_ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	  glEnable(GL_LIGHT0);
	  glEnable(GL_LIGHT1);
	  glEnable(GL_LIGHTING);
	  GLToolkit::SetMaterial (&fillMaterial);
	  glEnable(GL_POLYGON_OFFSET_FILL);
	  glCallList(polyList);
	  glDisable(GL_POLYGON_OFFSET_FILL);
	  GLToolkit::SetMaterial (&whiteMaterial);    
	  glDisable(GL_LIGHTING);
  } else {

	  // Default material
    GLToolkit::SetMaterial (&whiteMaterial);    

    // Draw lines
    glDisable(GL_CULL_FACE);
    glDisable(GL_LIGHTING);
	
	float color=(mApp->whiteBg)?0.0f:1.0f; //whitebg here
    if(viewStruct==-1) {      
		glColor4f(color,color,color,0.5f);
		if (mApp->antiAliasing) {
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
		}
		glPolygonOffset(1.0f,5.0f);
		glEnable(GL_POLYGON_OFFSET_LINE);
		for(int i=0;i<sh.nbSuper;i++)
			glCallList(lineList[i]);
		glDisable(GL_POLYGON_OFFSET_LINE);
		glDisable(GL_BLEND);
		glDisable(GL_LINE_SMOOTH);
		glColor3f(1.0f,1.0f,1.0f);
    } else {

      // Draw non selectable facet in dark grey
      glColor3f(0.2f,0.2f,0.2f);
	  		if (mApp->antiAliasing) {
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
		}
      for(int i=0;i<sh.nbSuper;i++)
        if(viewStruct!=i) 
          glCallList(lineList[i]);
      // Selectable in white
      glColor3f(color,color,color);
	  		glPolygonOffset(1.0f,5.0f);
		glEnable(GL_POLYGON_OFFSET_LINE);
      glCallList(lineList[viewStruct]);
		glDisable(GL_POLYGON_OFFSET_LINE);
		if (mApp->antiAliasing) {
			 glDisable(GL_BLEND);
			 glDisable(GL_LINE_SMOOTH);
		 }
    }

  }

  // Paint texture
  if( renderTexture ) {
    glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE,GL_ZERO);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0f,3.0f);
    for(int i=0;i<sh.nbFacet && renderTexture;i++) {
      Facet *f = facets[i];
      
	  if( f->sh.isTextured && f->textureVisible ) {
        if( f->sh.is2sided )   glDisable(GL_CULL_FACE);
        else                   SetCullMode(showMode);
        glBindTexture(GL_TEXTURE_2D,f->glTex);
        if( filter ) {
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        } else {
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        }
        glCallList(f->glList);
      }
    }
	glDisable(GL_POLYGON_OFFSET_FILL);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_BLEND);
	glDisable(GL_TEXTURE_2D);
  }

  // Paint mesh
  if( showMesh ) {
	  glColor4f(0.7f,0.7f,0.7f,0.3f);
	  if (mApp->antiAliasing) {
		  glEnable(GL_BLEND);
		  glEnable(GL_LINE_SMOOTH);
	  }
	  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

    for( int i=0; i<sh.nbFacet ;i++) {
      
	  Facet *f = facets[i];
      //if( f->selected && f->mesh ) glCallList(f->glElem);
	  if( f->cellPropertiesIds && f->textureVisible ) {
		    glEnable(GL_POLYGON_OFFSET_LINE);
			glPolygonOffset(1.0f,2.0f);
			glCallList(f->glElem);
			glDisable(GL_POLYGON_OFFSET_LINE);
	  }
    }
	if (mApp->antiAliasing) {
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
	}
  }

  // Paint direction fields
  if( showDir ) {

    GLToolkit::SetMaterial (&arrowMaterial);
    for(int i=0;i<sh.nbFacet;i++) {
      Facet *f = facets[i];
      if( f->sh.countDirection && f->dirCache ) {
        double iw = 1.0 / (double)f->sh.texWidthD;
        double ih = 1.0 / (double)f->sh.texHeightD;
        double rw = Norme(f->sh.U) * iw;
        for(int x=0;x<f->sh.texWidth;x++) {
          for(int y=0;y<f->sh.texHeight;y++) {           
             int add = x + y*f->sh.texWidth;
             if( f->GetMeshArea(add) > 0.0 ) {
               double uC = ((double)x + 0.5) * iw;
               double vC = ((double)y + 0.5) * ih;
               float xc = (float)( f->sh.O.x + f->sh.U.x*uC + f->sh.V.x*vC );
               float yc = (float)( f->sh.O.y + f->sh.U.y*uC + f->sh.V.y*vC );
               float zc = (float)( f->sh.O.z + f->sh.U.z*uC + f->sh.V.z*vC );
               RenderArrow(matView,
                           (float)f->dirCache[add].dir.x,
                           (float)f->dirCache[add].dir.y,
                           (float)f->dirCache[add].dir.z,
                           xc,yc,zc,(float)rw);
             }
          }
        }
      }
    }

    // Restore default matrix
    glLoadMatrixf(matView);
  }

  // Paint selection
  if(nbSelected) {
	  if (mApp->antiAliasing) {
	  glEnable(GL_BLEND);
	  glEnable(GL_LINE_SMOOTH);
	  }
	  glBlendFunc(GL_ONE,GL_ZERO);
    glColor3f(1.0f,0.0f,0.0f);    //red
    if( showHidden ) {
		glDisable(GL_DEPTH_TEST);
      glCallList(selectList3);
	  glEnable(GL_DEPTH_TEST);
    } else {
      //if( renderVolume || renderTexture ) {
		   glCallList(selectList3);
	  //} else {
	  //   glCallList(selectList3);
	  //}
    }
	if (mApp->antiAliasing) {
	glDisable(GL_LINE_SMOOTH);
	glDisable(GL_BLEND);
	}
	
  }

  // Paint selected cell on mesh
  for( int i=0; i<sh.nbFacet ;i++) {
    Facet *f = facets[i];
    f->RenderSelectedElem();
  }

}

void Geometry::DeleteGLLists(BOOL deletePoly,BOOL deleteLine) {
	if (deleteLine) {
		for(int i=0;i<sh.nbSuper;i++)
			DELETE_LIST(lineList[i]);
	}
	if (deletePoly) DELETE_LIST(polyList);
	DELETE_LIST(selectList);
	DELETE_LIST(selectList2);
	DELETE_LIST(selectList3);
}



void Geometry::RenderArrow(GLfloat *matView,float dx,float dy,float dz,float px,float py,float pz,float d) {

  if(!arrowList) BuildShapeList();

  // Compute transformation matrix for the arrow
  GLMatrix mView;
  GLMatrix aView;
  GLMatrix mScale;
  GLMatrix mRot;
  GLMatrix mT;
  float v1x,v1y,v1z;

  mView.LoadGL(matView);

  // Direction
  float n = sqrtf(dx*dx + dy*dy + dz*dz);
  if(IS_ZERO(n)) {

    // Isotropic (render a sphere)
    mScale._11 = (d/4.0f);
    mScale._22 = (d/4.0f);
    mScale._33 = (d/4.0f);

    mT.Translate(px,py,pz);

    aView.Multiply(&mView);
    aView.Multiply(&mT);
    aView.Multiply(&mScale);

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(aView.GetGL());
    glCallList(sphereList);
    return;

  }

  dx /= n;
  dy /= n;
  dz /= n;
  mRot._11 = dx;
  mRot._21 = dy;
  mRot._31 = dz;

  // A point belonging to the plane
  // normal to the direction vector
  if( !IS_ZERO(dx) ) {
    v1x = -dz / dx;
    v1y = 0.0;
    v1z = 1.0;
  } else if( !IS_ZERO(dy) ) {
    v1x = 0.0;
    v1y = -dz / dy;
    v1z = 1.0;
  } else if( !IS_ZERO(dz) ) {
    // normal to z
    v1x = 1.0;
    v1y = 0.0;
    v1z = 0.0;
  } else {
    // Null vector -> isotropic
  }

  float n1 = sqrtf(v1x*v1x + v1y*v1y + v1z*v1z);
  v1x /= n1;
  v1y /= n1;
  v1z /= n1;
  mRot._12 = v1x;
  mRot._22 = v1y;
  mRot._32 = v1z;

  // Cross product
  mRot._13 = (dy)*(v1z) - (dz)*(v1y);
  mRot._23 = (dz)*(v1x) - (dx)*(v1z);
  mRot._33 = (dx)*(v1y) - (dy)*(v1x);

  // Scale
  if( !autoNorme ) {
    mScale._11 = (n*d*normeRatio);
    mScale._22 = (d/4.0f);
    mScale._33 = (d/4.0f);
  } else {
    // Show only direction
    mScale._11 = (d/1.1f);
    mScale._22 = (d/4.0f);
    mScale._33 = (d/4.0f);
  }

  mT.Translate(px,py,pz);

  aView.Multiply(&mView);
  aView.Multiply(&mT);
  aView.Multiply(&mRot);
  aView.Multiply(&mScale);
  if( !centerNorme ) {
    mT.Translate(0.5f,0.0,0.0);
    aView.Multiply(&mT);
  }

  glMatrixMode(GL_MODELVIEW);
  glLoadMatrixf(aView.GetGL());
  glCallList(arrowList);

}

// ---------------------------------------------------------------
// Triangulation stuff
// ---------------------------------------------------------------

int Geometry::FindEar(POLYGON *p) {

  int i = 0;
  BOOL earFound = FALSE;
  while( i<p->nbPts && !earFound ) { 
    if( IsConvex(p,i) )
      earFound = !ContainsConcave(p,i-1,i,i+1);
    if(!earFound) i++;
  }

  // REM: Theoritically, it should always find an ear (2-Ears theorem).
  // However on degenerated geometry (flat poly) it may not find one.
  // Returns first point in case of failure.
  if( earFound )
    return i;
  else
    return 0;

}



void Geometry::AddTextureCoord(Facet *f,VERTEX2D *p) {

  // Add texture coord with a 1 texel border (for bilinear filtering)
  double uStep = 1.0/(double)f->texDimW;
  double vStep = 1.0/(double)f->texDimH;

#if 1
  double fu = f->sh.texWidthD * uStep;
  double fv = f->sh.texHeightD * vStep;
  glTexCoord2f((float)(uStep + p->u*fu),(float)(vStep + p->v*fv));
#else
  // Show border (debugging purpose)
  double fu = (f->sh.texWidthD+2.0) * uStep;
  double fv = (f->sh.texHeightD+2.0) * vStep;
  glTexCoord2f((float)(p->u*fu),(float)(p->v*fv));
#endif

}



void Geometry::FillFacet(Facet *f,BOOL addTextureCoord) {
	
  for(int i=0;i<f->sh.nbIndex;i++ ) {
    int idx = f->indices[i];
    glNormal3d( -f->sh.N.x , -f->sh.N.y , -f->sh.N.z );
    if(addTextureCoord) AddTextureCoord(f,f->vertices2+i);
    glVertex3d( vertices3[idx].x , vertices3[idx].y , vertices3[idx].z);
  }
  
}



void Geometry::DrawEar(Facet *f,POLYGON *p,int ear,BOOL addTextureCoord) {

  VERTEX3D  p3D;
  VERTEX2D *p1;
  VERTEX2D *p2;
  VERTEX2D *p3;

  // Follow orientation
  if( p->sign > 0.0 ) {
    p1 = &(p->pts[IDX(ear-1,p->nbPts)]);
    p2 = &(p->pts[IDX(ear+1,p->nbPts)]);
    p3 = &(p->pts[IDX(ear  ,p->nbPts)]);
  } else {
    p1 = &(p->pts[IDX(ear-1,p->nbPts)]);
    p2 = &(p->pts[IDX(ear  ,p->nbPts)]);
    p3 = &(p->pts[IDX(ear+1,p->nbPts)]);
  }

  glNormal3d( -f->sh.N.x , -f->sh.N.y , -f->sh.N.z );
  if(addTextureCoord) AddTextureCoord(f,p1);
  // (U,V) -> (x,y,z)
  p3D.x = f->sh.O.x + p1->u*f->sh.U.x + p1->v*f->sh.V.x;
  p3D.y = f->sh.O.y + p1->u*f->sh.U.y + p1->v*f->sh.V.y;
  p3D.z = f->sh.O.z + p1->u*f->sh.U.z + p1->v*f->sh.V.z;
  glVertex3d( p3D.x , p3D.y , p3D.z);

  glNormal3d( -f->sh.N.x , -f->sh.N.y , -f->sh.N.z );
  if(addTextureCoord) AddTextureCoord(f,p2);
  // (U,V) -> (x,y,z)
  p3D.x = f->sh.O.x + p2->u*f->sh.U.x + p2->v*f->sh.V.x;
  p3D.y = f->sh.O.y + p2->u*f->sh.U.y + p2->v*f->sh.V.y;
  p3D.z = f->sh.O.z + p2->u*f->sh.U.z + p2->v*f->sh.V.z;
  glVertex3d( p3D.x , p3D.y , p3D.z);

  glNormal3d( -f->sh.N.x , -f->sh.N.y , -f->sh.N.z );
  if(addTextureCoord) AddTextureCoord(f,p3);
  // (U,V) -> (x,y,z)
  p3D.x = f->sh.O.x + p3->u*f->sh.U.x + p3->v*f->sh.V.x;
  p3D.y = f->sh.O.y + p3->u*f->sh.U.y + p3->v*f->sh.V.y;
  p3D.z = f->sh.O.z + p3->u*f->sh.U.z + p3->v*f->sh.V.z;
  glVertex3d( p3D.x , p3D.y , p3D.z);

}



void Geometry::Triangulate(Facet *f,BOOL addTextureCoord) {

  // Triangulate a facet (rendering purpose)
  // The facet must have at least 3 points
  // Use the very simple "Two-Ears" theorem. It computes in O(n^2).

  // Build a POLYGON
  POLYGON p;
  p.nbPts = f->sh.nbIndex;
  p.pts = (VERTEX2D *)malloc(p.nbPts*sizeof(VERTEX2D));
  memcpy(p.pts,f->vertices2,p.nbPts*sizeof(VERTEX2D));
  p.sign = f->sh.sign;

  if( f->sh.sign==0.0 ) {
    // Not a simple polygon
    // Abort triangulation
    free(p.pts);
    return;
  }

  // Perform triangulation
  while( p.nbPts > 3 ) {
    int e = FindEar(&p);
    DrawEar(f,&p,e,addTextureCoord);
    // Remove the ear
    for(int i=e;i<p.nbPts-1;i++)
      p.pts[i] = p.pts[i+1];
    p.nbPts--;
  }

  // Draw the last ear
  DrawEar(f,&p,0,addTextureCoord);
  free(p.pts);

}
