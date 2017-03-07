/*
File:        GeometryViewer.cpp
Description: Geometry 3D Viewer component
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
#include "GLApp/GLApp.h"
#include "GLApp/GLWindowManager.h"
#include "GeometryViewer.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMatrix.h"
#include "GLApp/MathTools.h"

#include <math.h>
#include <malloc.h>
#include "SynRad.h"

extern SynRad *mApp;

void GeometryViewer::SetBounds(int x,int y,int width,int height) {
	if( this->posX != x || this->posY != y || this->width != width || this->height != height ) {
		GLComponent::SetBounds(x,y,width,height);

		frontBtn->SetBounds(posX+5,posY+height-22,33,19);
		topBtn->SetBounds(posX+40,posY+height-22,33,19);
		sideBtn->SetBounds(posX+75,posY+height-22,33,19);
		projCombo->SetBounds(posX+110,posY+height-22,50,19);
		coordLab->SetBounds(posX+162,posY+height-20,100,18);
		facetSearchState->SetBounds(posX+10,posY+10,90,19);


		autoBtn->SetBounds(posX+width-142,posY+height-22,19,19);
		selBtn->SetBounds(posX+width-122,posY+height-22,19,19);
		selVxBtn->SetBounds(posX+width-102,posY+height-22,19,19);

		selTrajBtn->SetBounds(posX+width-82,posY+height-22,19,19);
		zoomBtn->SetBounds(posX+width-62,posY+height-22,19,19);
		handBtn->SetBounds(posX+width-42,posY+height-22,19,19);
		sysBtn->SetBounds(posX+width-22,posY+height-22,19,19);

		hideLotlabel->SetBounds(posX + 10, posY + height - 67, 0, 19);
		capsLockLabel->SetBounds(posX + 10, posY + height - 47, 0, 19);

		toolBack->SetBounds(posX+1,posY+height-DOWN_MARGIN,width-2,DOWN_MARGIN-1);
		SetCurrentView(view);
	}
}

void GeometryViewer::DrawLinesAndHits() {

	// Draw Lines and Hits
	if((showLine || showHit) && work->nbHit>0) {

		// Retrieve hit data
		HIT pHits[NBHHIT];
		int nbHHit;
		work->GetHHit(pHits,&nbHHit);

		if( nbHHit ) {

			// Lines
			if( showLine ) {

				glDisable(GL_TEXTURE_2D);
				glDisable(GL_LIGHTING);

				glDisable(GL_CULL_FACE);

				//Determine MAX dF value
				int count=0;

				double logOpacityMax=-99;
				double logOpacityMin = 99;
				if (shadeLines) {
					while(count<dispNumHits && pHits[count].type!=0) {
						double logVal;
						if (work->generation_mode == SYNGEN_MODE_FLUXWISE) {
							logVal = log10(pHits[count].dP);
						}
						else {
							logVal = log10(pHits[count].dF);
						}
						if (logVal > logOpacityMax) logOpacityMax = logVal;
						if (logVal < logOpacityMin) logOpacityMin = logVal;
						count++;
					}
				}
				count=0;
				logOpacityMin = MAX(logOpacityMin, logOpacityMax-6.0); //Span through max. 6 orders of magnitude
				double opacitySpan = logOpacityMax - logOpacityMin;

				while(count<dispNumHits && pHits[count].type!=0) {
					//if (count>0&&pHits[count].type==HIT_DES&&pHits[count-1].type!=HIT_ABS) __debugbreak(); //desorbed without being absorbed first
					float lineOpacity;
					
					if (mApp->whiteBg) { //whitebg
						lineOpacity=1.0f;
						glColor4f(0.2f,0.7f,0.2f,lineOpacity);
						
					}
					else {
						if (shadeLines) {
							double logVal;
							if (work->generation_mode == SYNGEN_MODE_FLUXWISE) {
								logVal = log10(pHits[count].dP);
							}
							else {
								logVal = log10(pHits[count].dF);
							}
							lineOpacity = (float)((logVal-logOpacityMin)/opacitySpan);
						}
						else lineOpacity=1.0;
						glColor4f(0.5f,1.0f,0.5f,lineOpacity);
					}
					if (mApp->antiAliasing) {
						glEnable(GL_BLEND);//,glEnable(GL_LINE_SMOOTH);
						//glLineWidth(	2.0f);
					}

					glBegin(GL_LINE_STRIP);
					while(count<dispNumHits && pHits[count].type!=HIT_ABS) {

						//teleport routine
						if (pHits[count].type==HIT_TELEPORT) {
							glVertex3d(pHits[count].pos.x , pHits[count].pos.y , pHits[count].pos.z);
							glEnd();
							if (showTP) {
								if (!mApp->whiteBg) {
									glColor3f(1.0f,0.7f,0.2f);
								} else {
									glColor3f(1.0f,0.0f,1.0f);
								}
								glPushAttrib(GL_ENABLE_BIT); 

								glLineStipple(1, 0x0101);
								glEnable(GL_LINE_STIPPLE);
								glBegin(GL_LINE_STRIP);
								glVertex3d(pHits[count].pos.x , pHits[count].pos.y , pHits[count].pos.z); //source point
								count++;            
								glVertex3d(pHits[count].pos.x , pHits[count].pos.y , pHits[count].pos.z);  //teleport dest.
								glEnd();
								glPopAttrib();

								if (mApp->whiteBg) { //whitebg
									glColor4f(0.2f,0.7f,0.2f,1.0);
						
								} else {
									if (shadeLines) {
										double logVal;
										if (work->generation_mode == SYNGEN_MODE_FLUXWISE) {
											logVal = log10(pHits[count].dP);
										}
										else {
											logVal = log10(pHits[count].dF);
										}
										lineOpacity = (float)((logVal - logOpacityMin) / opacitySpan);
									}
									else lineOpacity = 1.0;
									glColor4f(0.5f,1.0f,0.5f,lineOpacity);
								} 
								glBegin(GL_LINE_STRIP);
							} else {
								
								//glVertex3d(pHits[count].pos.x , pHits[count].pos.y , pHits[count].pos.z); //source point
								count++;   
								glBegin(GL_LINE_STRIP);
								glVertex3d(pHits[count].pos.x , pHits[count].pos.y , pHits[count].pos.z);  //teleport dest.

							}
						}
						/*if (pHits[count].type==HIT_DES) {
							glEnd();glBegin(GL_LINE_STRIP); //pen up pen down for leaks
						}*/



						glVertex3d(pHits[count].pos.x , pHits[count].pos.y , pHits[count].pos.z);
						count++;
						if (pHits[count].type==LASTHIT) { //pen up at cache refresh border
							glEnd();
							count++;
							glBegin(GL_LINE_STRIP);
						}
					}
					if(count<dispNumHits && pHits[count].type!=0) {
						glVertex3d(pHits[count].pos.x , pHits[count].pos.y , pHits[count].pos.z);
						count++;
					}
					glEnd();
					if (mApp->antiAliasing){
						glDisable(GL_LINE_SMOOTH);
					glDisable(GL_BLEND);}
				}

			}

			// Hit
			if(showHit) {

				glDisable(GL_TEXTURE_2D);
				glDisable(GL_LIGHTING);
				glDisable(GL_BLEND);
				glDisable(GL_CULL_FACE);

				// Refl

				float pointSize=(bigDots)?2.0f:1.0f;
				glPointSize(pointSize);
				if (mApp->whiteBg) { //whitebg
					glColor3f(0.2f,0.2f,0.2f);
				}
				else {
					glColor3f(0.0f,1.0f,0.0f);
				}
				glBegin(GL_POINTS);
				for(int i=0;i<dispNumHits;i++)
					if(pHits[i].type==HIT_REF)
						glVertex3d(pHits[i].pos.x , pHits[i].pos.y , pHits[i].pos.z);
				glEnd();

				// Trans
				pointSize=(bigDots)?3.0f:2.0f;
				glPointSize(pointSize);
				glColor3f(0.5f,1.0f,1.0f);
				glBegin(GL_POINTS);
				for(int i=0;i<dispNumHits;i++)
					if(pHits[i].type==HIT_TRANS)
						glVertex3d(pHits[i].pos.x , pHits[i].pos.y , pHits[i].pos.z);
				glEnd();

				// Teleport
				if (showTP) {
					//pointSize=(bigDots)?3.0f:2.0f;
					glPointSize(pointSize);
					if (!mApp->whiteBg) {
						glColor3f(1.0f,0.7f,0.2f);
					} else {
						glColor3f(1.0f,0.0f,1.0f);
					}
					glBegin(GL_POINTS);
					for(int i=0;i<dispNumHits;i++)
						if(pHits[i].type==HIT_TELEPORT)
							glVertex3d(pHits[i].pos.x , pHits[i].pos.y , pHits[i].pos.z);
					glEnd();
				}

				// Abs
				glPointSize(pointSize);
				glColor3f(1.0f,0.0f,0.0f);
				glBegin(GL_POINTS);
				for(int i=0;i<dispNumHits;i++)
					if(pHits[i].type==HIT_ABS)
						glVertex3d(pHits[i].pos.x , pHits[i].pos.y , pHits[i].pos.z);
				glEnd();

				// Des
				glColor3f(0.3f,0.3f,1.0f);
				glBegin(GL_POINTS);
				for(int i=0;i<dispNumHits;i++)
					if(pHits[i].type==HIT_DES)
						glVertex3d(pHits[i].pos.x , pHits[i].pos.y , pHits[i].pos.z);
				glEnd();
			}
		}
	}
}



/*

void GeometryViewer::DrawBB() {
if( showLeak ) {
glDisable(GL_TEXTURE_2D);
glDisable(GL_LIGHTING);
glDisable(GL_BLEND);
glDisable(GL_CULL_FACE);
glColor3f(1.0f,1.0f,0.0f);
glBegin(GL_LINES);
DrawBB(geom->aabbTree);
glEnd();
}
}

void GeometryViewer::DrawBB(AABBNODE *node) {

if( node ) {

if( node->left==NULL && node->right==NULL ) {

// Leaf
glVertex3d(node->bb.min.x, node->bb.min.y, node->bb.min.z);
glVertex3d(node->bb.max.x, node->bb.min.y, node->bb.min.z);

glVertex3d(node->bb.max.x, node->bb.min.y, node->bb.min.z);
glVertex3d(node->bb.max.x, node->bb.min.y, node->bb.max.z);

glVertex3d(node->bb.max.x, node->bb.min.y, node->bb.max.z);
glVertex3d(node->bb.min.x, node->bb.min.y, node->bb.max.z);

glVertex3d(node->bb.min.x, node->bb.min.y, node->bb.max.z);
glVertex3d(node->bb.min.x, node->bb.min.y, node->bb.min.z);

glVertex3d(node->bb.min.x, node->bb.min.y, node->bb.min.z);
glVertex3d(node->bb.min.x, node->bb.max.y, node->bb.min.z);

glVertex3d(node->bb.max.x, node->bb.min.y, node->bb.min.z);
glVertex3d(node->bb.max.x, node->bb.max.y, node->bb.min.z);

glVertex3d(node->bb.max.x, node->bb.min.y, node->bb.max.z);
glVertex3d(node->bb.max.x, node->bb.max.y, node->bb.max.z);

glVertex3d(node->bb.min.x, node->bb.min.y, node->bb.max.z);
glVertex3d(node->bb.min.x, node->bb.max.y, node->bb.max.z);

glVertex3d(node->bb.min.x, node->bb.max.y, node->bb.min.z);
glVertex3d(node->bb.max.x, node->bb.max.y, node->bb.min.z);

glVertex3d(node->bb.max.x, node->bb.max.y, node->bb.min.z);
glVertex3d(node->bb.max.x, node->bb.max.y, node->bb.max.z);

glVertex3d(node->bb.max.x, node->bb.max.y, node->bb.max.z);
glVertex3d(node->bb.min.x, node->bb.max.y, node->bb.max.z);

glVertex3d(node->bb.min.x, node->bb.max.y, node->bb.max.z);
glVertex3d(node->bb.min.x, node->bb.max.y, node->bb.min.z);

} else {
DrawBB(node->left);
DrawBB(node->right);
}

}

}
*/



#define TRANSFORMBB( X,Y,Z )                                                \
	mv.TransfomVec((float)bbO.X,(float)bbO.Y,(float)bbO.Z,1.0f,&rx,&ry,&rz,&rw);\
	dx = (double)rx;                                                            \
	dy = (double)ry;                                                            \
	dz = (double)rz;                                                            \
	if( dx < xMin ) xMin = dx;                                                  \
	if( dy < yMin ) yMin = dy;                                                  \
	if( dz < zNear) zNear = dz;                                                 \
	if( dx > xMax ) xMax = dx;                                                  \
	if( dy > yMax ) yMax = dy;                                                  \
	if( dz > zFar ) zFar = dz;

#define TRANSFORMVERTEX( X,Y,Z )                                  \
	mv.TransfomVec((float)X,(float)Y,(float)Z,1.0f,&rx,&ry,&rz,&rw);  \
	dx = (double)rx;                                                  \
	dy = (double)ry;                                                  \
	dz = (double)rz;                                                  \
	if( dx < xMin ) xMin = dx;                                        \
	if( dy < yMin ) yMin = dy;                                        \
	if( dz < zNear) zNear = dz;                                       \
	if( dx > xMax ) xMax = dx;                                        \
	if( dy > yMax ) yMax = dy;                                        \
	if( dz > zFar ) zFar = dz;

