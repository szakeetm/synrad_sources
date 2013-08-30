/*
File:        GeometryViewer.cpp
Description: Geometry 3D Viewer component
Program:     SynRad
Author:      R. KERSEVAN / M SZAKACS
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
#include "Utils.h"
#include <math.h>
#include <malloc.h>


#define DOWN_MARGIN 25

static const GLfloat position[]  = { -100000.0f, 100000.0f, -120000.0f, 1.0f }; //light1
static const GLfloat positionI[] = {  100000.0f,-100000.0f,  120000.0f, 1.0f }; //light2

extern GLApplication *theApp;

// ---------------------------------------------------------------------

GeometryViewer::GeometryViewer(int id):GLComponent(id) {

	work = NULL;

	// Material

	memset (&greenMaterial, 0, sizeof (GLMATERIAL));
	greenMaterial.Diffuse.r = 0.0f;
	greenMaterial.Diffuse.g = 1.0f;
	greenMaterial.Diffuse.b = 0.0f;
	greenMaterial.Ambient.r = 0.0f;
	greenMaterial.Ambient.g = 1.0f;
	greenMaterial.Ambient.b = 0.0f;

	memset (&blueMaterial, 0, sizeof (GLMATERIAL));
	blueMaterial.Diffuse.r = 0.5f;
	blueMaterial.Diffuse.g = 0.6f;
	blueMaterial.Diffuse.b = 1.0f;
	blueMaterial.Ambient.r = 0.5f;
	blueMaterial.Ambient.g = 0.6f;
	blueMaterial.Ambient.b = 1.0f;

	/// Default values
	draggMode = DRAGG_NONE;
	selected = FALSE;
	view.projMode = PERSPECTIVE_PROJ;
	view.camAngleOx = 0.0;
	view.camAngleOy = 0.0;
	view.camDist = 100.0;
	view.camOffset.x = 0.0;
	view.camOffset.y = 0.0;
	view.camOffset.z = 0.0;
	view.vLeft = 0.0;
	view.vRight = 0.0;
	view.vTop = 0.0;
	view.vBottom = 0.0;
	view.name = NULL;
	view.performXY = XYZ_NONE;
	showIndex = FALSE;
	showVertex = FALSE;
	showNormal = FALSE;
	showUV = FALSE;
	showRule = TRUE;
	showLeak = FALSE;
	showHit = FALSE;
	showLine = FALSE;
	showVolume = FALSE;
	showTexture = FALSE;
	showHidden = FALSE;
	showHiddenVertex = FALSE;
	showMesh = FALSE;
	showDir = TRUE;
	autoScaleOn = FALSE;
	mode = MODE_SELECT;
	showBack = SHOW_FRONTANDBACK;
	showFilter = FALSE;
	showColormap = FALSE;
	camDistInc = 1.0;
	transStep = 1.0;
	angleStep = 0.005;
	selX1 = 0;
	selY1 = 0;
	selX2 = 0;
	selY2 = 0;
	selectionChange = FALSE;
	vectorLength = 5.0f;
	arrowLength = 1.0;
	dispNumHits = 2048;
	dispNumLeaks = 2048;
	dispNumTraj = 100;

	// GL Component default
	SetBorder(BORDER_NONE);
	int bgCol=(false)?255:0; //not necessary?
	SetBackgroundColor(bgCol,bgCol,bgCol);

	// Components
	toolBack = new GLLabel("");
	toolBack->SetBackgroundColor(220,220,220);
	toolBack->SetOpaque(TRUE);
	toolBack->SetBorder(BORDER_BEVEL_IN);
	Add(toolBack);
	coordLab = new GLLabel("");
	coordLab->SetBackgroundColor(220,220,220);
	Add(coordLab);

	facetSearchState = new GLLabel("");
	//facetSearchState->SetBackgroundColor(220,220,220);
	facetSearchState->SetVisible(FALSE);
	Add(facetSearchState);

	frontBtn = new GLButton(0,"Front");
	frontBtn->SetToggle(TRUE);
	Add(frontBtn);
	topBtn = new GLButton(0,"Top");
	topBtn->SetToggle(TRUE);
	Add(topBtn);
	sideBtn = new GLButton(0,"Side");
	sideBtn->SetToggle(TRUE);
	Add(sideBtn);

	projCombo = new GLCombo(0);
	projCombo->SetSize(2);
	projCombo->SetValueAt(0,"Persp.");
	projCombo->SetValueAt(1,"Ortho.");
	projCombo->SetSelectedIndex(0);
	Add(projCombo);

	capsLockLabel = new GLLabel("CAPS LOCK On: select facets only with selected vertex");
	capsLockLabel->SetTextColor(255,0,0);
	Add(capsLockLabel);

	zoomBtn = new GLButton(0,"");
	zoomBtn->SetIcon("images/icon_zoom.png");
	zoomBtn->SetDisabledIcon("images/icon_zoomd.png");
	zoomBtn->SetEnabled(FALSE);
	Add(zoomBtn);
	sysBtn = new GLButton(0,"");
	sysBtn->SetIcon("images/icon_winup.png");
	Add(sysBtn);
	handBtn = new GLButton(0,"");
	handBtn->SetIcon("images/icon_hand.png");
	Add(handBtn);
	selBtn = new GLButton(0,"");
	selBtn->SetIcon("images/icon_arrow.png");
	Add(selBtn);
	selVxBtn = new GLButton(0,"");
	selVxBtn->SetIcon("images/icon_vertex_select.png");
	Add(selVxBtn);
	selTrajBtn = new GLButton(0,"");
	selTrajBtn->SetIcon("images/icon_traj_select.png");
	Add(selTrajBtn);
	autoBtn = new GLButton(0,"");
	autoBtn->SetIcon("images/icon_autoscale.png");
	autoBtn->SetToggle(TRUE);
	Add(autoBtn);

	// Light
	glShadeModel(GL_SMOOTH);

	GLfloat global_ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);

	GLfloat ambientLight[]  = { 0.0f, 0.0f, 0.0f, 1.0f   };
	GLfloat diffuseLight[]  = { 0.8f, 0.8f, 0.8f, 0.0f   };
	GLfloat specularLight[] = { 1.0f, 1.0f, 1.0f, 1.0f   };

	GLfloat ambientLight2[]  = { 0.0f, 0.0f, 0.0f, 1.0f   };
	GLfloat diffuseLight2[]  = { 0.1f, 0.1f, 0.1f, 0.0f   };
	GLfloat specularLight2[] = { 0.3f, 0.3f, 0.3f, 0.3f   };

	glLightfv(GL_LIGHT0, GL_AMBIENT,  ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuseLight);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
	glDisable(GL_LIGHT0);

	glLightfv(GL_LIGHT1, GL_AMBIENT,  ambientLight2);
	glLightfv(GL_LIGHT1, GL_DIFFUSE,  diffuseLight2);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specularLight2);
	glDisable(GL_LIGHT1);

}

//----------------------------------------------------------------------
void GeometryViewer::ToOrigo() {
	//view.projMode = PERSPECTIVE_PROJ;
	view.camAngleOx = 0.0;
	view.camAngleOy = (view.projMode == PERSPECTIVE_PROJ)?0.0:PI;
	view.camDist = 100.0;
	view.camOffset.x = 0.0;
	view.camOffset.y = 0.0;
	view.camOffset.z = 0.0;

	view.vLeft = 0.0;
	view.vRight = 0.0;
	view.vTop = 0.0;
	view.vBottom = 0.0;
}

// ---------------------------------------------------------------------

void GeometryViewer::SetBounds(int x,int y,int width,int height) {
	if( this->posX != x || this->posY != y || this->width != width || this->height != height ) {
		GLComponent::SetBounds(x,y,width,height);

		frontBtn->SetBounds(posX+5,posY+height-22,33,19);
		topBtn->SetBounds(posX+40,posY+height-22,33,19);
		sideBtn->SetBounds(posX+75,posY+height-22,33,19);
		projCombo->SetBounds(posX+110,posY+height-22,50,19);
		capsLockLabel->SetBounds(posX+165,posY+height-22,0,19);
		coordLab->SetBounds(posX+162,posY+height-20,100,18);
		facetSearchState->SetBounds(posX+10,posY+10,90,19);

		autoBtn->SetBounds(posX+width-142,posY+height-22,19,19);
		selBtn->SetBounds(posX+width-122,posY+height-22,19,19);
		selVxBtn->SetBounds(posX+width-102,posY+height-22,19,19);
		selTrajBtn->SetBounds(posX+width-82,posY+height-22,19,19);
		zoomBtn->SetBounds(posX+width-62,posY+height-22,19,19);
		handBtn->SetBounds(posX+width-42,posY+height-22,19,19);
		sysBtn->SetBounds(posX+width-22,posY+height-22,19,19);

		toolBack->SetBounds(posX+1,posY+height-DOWN_MARGIN,width-2,DOWN_MARGIN-1);
		SetCurrentView(view);
	}
}

// ---------------------------------------------------------------------

BOOL GeometryViewer::IsSelected() {
	return selected;
}

void GeometryViewer::SetSelected(BOOL s) {
	selected = s;
}

// ---------------------------------------------------------------------

void GeometryViewer::SetFocus(BOOL focus) {
	if(focus && parent)  parent->ProcessMessage(this,MSG_GEOMVIEWER_SELECT);
	GLComponent::SetFocus(focus);
}

// ---------------------------------------------------------------------

void GeometryViewer::UpdateMouseCursor(int mode) { //Sets mouse cursor to action

	this->mode = mode;

	if(parent) {

		if(!((draggMode==DRAGG_ROTATE)||(draggMode==DRAGG_MOVE))) {
			switch(mode) {



			case MODE_SELECT:
				if (GetWindow()->IsCtrlDown()) {
					SetCursor(CURSOR_SELDEL);
				}
				else if (GetWindow()->IsShiftDown()) {
					SetCursor(CURSOR_SELADD);
				}
				else if (GetWindow()->IsAltDown()) {
					SetCursor(CURSOR_HAND);
				}
				else {
					SetCursor(CURSOR_DEFAULT);
				}

				break;

			case MODE_SELECTVERTEX:
				if (GetWindow()->IsCtrlDown()) {
					SetCursor(CURSOR_VERTEX_CLR);
				}
				else if (GetWindow()->IsShiftDown()) {
					SetCursor(CURSOR_VERTEX_ADD);
				}
				else if (GetWindow()->IsAltDown()) {
					SetCursor(CURSOR_HAND);
				}
				else {
					SetCursor(CURSOR_VERTEX);
				}

				break;
							
			case MODE_SELECTTRAJ:
				if (GetWindow()->IsAltDown()) {
					SetCursor(CURSOR_HAND);
				}
				else {
					SetCursor(CURSOR_TRAJ);
				}

				break;

			case MODE_ZOOM:
				SetCursor(CURSOR_ZOOM);
				break;
			case MODE_MOVE:
				SetCursor(CURSOR_HAND);

				break;
			}	
		}
		if (draggMode==DRAGG_ROTATE) SetCursor(CURSOR_ROTATE);
		if (draggMode==DRAGG_MOVE) SetCursor(CURSOR_HAND);

	}
}

// ---------------------------------------------------------------------

BOOL GeometryViewer::IsDragging() {
	return draggMode != DRAGG_NONE;
}

// ---------------------------------------------------------------------

void GeometryViewer::ToTopView() {

	if(!work) return;
	view.camAngleOx = PI/2.0;
	view.camAngleOy = 0.0;
	view.performXY=(view.projMode==PERSPECTIVE_PROJ)?XYZ_NONE:XYZ_TOP;
	AutoScale();

}

void GeometryViewer::ToSideView() {

	if(!work) return;

	view.camAngleOx = 0.0;
	view.camAngleOy = PI/2.0;
	view.performXY=(view.projMode==PERSPECTIVE_PROJ)?XYZ_NONE:XYZ_SIDE;
	AutoScale();

}

void GeometryViewer::ToFrontView() {
	ToOrigo();
	if(!work) return;
	//view.camAngleOx = 0.0;
	//view.camAngleOy = 0.0;
	view.performXY=(view.projMode==PERSPECTIVE_PROJ)?XYZ_NONE:XYZ_FRONT;
	AutoScale();

}

//-----------------------------------------------------------------------------

void GeometryViewer::UpdateLight() {
	
	float ratio = 1.0f;

	// Update lights
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	
	if( view.projMode==PERSPECTIVE_PROJ ) {
		glRotated(ToDeg(-view.camAngleOx),1.0,0.0,0.0);
		glRotated(ToDeg(-view.camAngleOy),0.0,1.0,0.0);
		ratio = 1.0f;
	} else {
		glScaled(1.0,-1.0,1.0);
		glRotated(ToDeg(-view.camAngleOx),1.0,0.0,0.0);
		glRotated(ToDeg(-view.camAngleOy),0.0,1.0,0.0);
		ratio = (float)view.camDist;
	}
	


	GLfloat d0[4],d1[4];
	d0[0] = 0.8f * ratio;
	d0[1] = 0.8f * ratio;
	d0[2] = 0.8f * ratio;
	d0[3] = 0.0;

	d1[0] = 0.15f * ratio;
	d1[1] = 0.15f * ratio;
	d1[2] = 0.15f * ratio;
	d1[3] = 0.0f;

	if( showBack==1 ) {
		glLightfv(GL_LIGHT0, GL_DIFFUSE,   d0);
		glLightfv(GL_LIGHT0, GL_POSITION,  positionI);
		glLightfv(GL_LIGHT1, GL_DIFFUSE,   d1);
		glLightfv(GL_LIGHT1, GL_POSITION,  position);
	} else {
		glLightfv(GL_LIGHT0, GL_DIFFUSE,   d0);
		glLightfv(GL_LIGHT0, GL_POSITION,       position);
		glLightfv(GL_LIGHT1, GL_DIFFUSE,   d1);
		glLightfv(GL_LIGHT1, GL_POSITION,       positionI);
	}
	
}

//-----------------------------------------------------------------------------

void GeometryViewer::UpdateMatrix() {

	if(!work) return;
	Geometry *geom = work->GetGeometry();
	if(!geom) return;

	// Model view matrix ---------------------------------------------------

	// Scale angle in -PI,PI
	view.camAngleOx = RoundAngle(view.camAngleOx);
	view.camAngleOy = RoundAngle(view.camAngleOy);

	// Convert polar coordinates
	VERTEX3D org = geom->GetCenter();

	camDir.x  = -cos(view.camAngleOx) * sin(view.camAngleOy);
	camDir.y  =  sin(view.camAngleOx);
	camDir.z  = -cos(view.camAngleOx) * cos(view.camAngleOy);

	camLeft.x = -cos(view.camAngleOy);
	camLeft.y = 0.0;
	camLeft.z = sin(view.camAngleOy);

	Cross(&camUp,&camDir,&camLeft);

	glMatrixMode( GL_MODELVIEW );

	switch( view.projMode ) {
	case PERSPECTIVE_PROJ:
		GLToolkit::LookAtLH((camDir.x * view.camDist+org.x)+ view.camOffset.x,
			(camDir.y * view.camDist+org.y)+ view.camOffset.y,
			(camDir.z * view.camDist+org.z)+ view.camOffset.z,
			org.x + view.camOffset.x , org.y + view.camOffset.y , org.z + view.camOffset.z,
			camUp.x,camUp.y,camUp.z );
		break;
	case ORTHOGRAPHIC_PROJ:
		glLoadIdentity();
		glScaled(-view.camDist,-view.camDist,view.camDist);
		glRotated(ToDeg(-view.camAngleOx),1.0,0.0,0.0);
		glRotated(ToDeg(view.camAngleOy),0.0,1.0,0.0);
		glTranslated(-(org.x-view.camOffset.x),-(org.y+view.camOffset.y),-(org.z+view.camOffset.z));
		break;
	}

	glGetFloatv( GL_MODELVIEW_MATRIX , matView );

	// Projection matrix ---------------------------------------------------

	double aspect = (double)width/(double)(height-DOWN_MARGIN);
	ComputeBB(TRUE);

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();

	if( view.projMode==PERSPECTIVE_PROJ ) {

		double _zNear = MAX(zNear,0.1);
		double _zFar  = (_zNear<zFar)?zFar:_zNear+1.0;
		GLToolkit::PerspectiveLH(FOV_ANGLE,aspect,_zNear-0.05,_zFar+10.0);

	} else {

		// 30% margin for extra geometry
		double l = zFar - zNear;
		if( (l>0.0) && (view.vRight - view.vLeft)>0.0 && (view.vBottom - view.vTop)>0.0 )
			glOrtho(view.vLeft,view.vRight,view.vBottom,view.vTop,-zFar-l*0.3,-zNear+l*0.3);

	}

	glGetFloatv( GL_PROJECTION_MATRIX , matProj );

}

// ---------------------------------------------------------------------

double GeometryViewer::ToDeg(double radians) {
	return (radians/PI)*180.0f;
}

// ---------------------------------------------------------------------

BOOL GeometryViewer::SelectionChanged() {
	BOOL ret = selectionChange;
	selectionChange = FALSE;
	return ret;
}

// ---------------------------------------------------------------------

AVIEW GeometryViewer::GetCurrentView() {
	return view;
}

// ---------------------------------------------------------------------

void GeometryViewer::SetCurrentView(AVIEW v) {

	if(!work) return;
	view = v;

	if( view.projMode == ORTHOGRAPHIC_PROJ ) {

		// Rescale viewport (aspect ratio correction + recenter)
		double wA = (double)width/(double)(height-DOWN_MARGIN);
		double vA = (view.vRight - view.vLeft)/(view.vBottom - view.vTop);
		double AA = vA / wA;
		double c = (view.vTop + view.vBottom) / 2.0;
		double l = AA * fabs(view.vTop-view.vBottom) / 2.0;
		view.vTop = c - l;
		view.vBottom = c + l;

	}

	UpdateMatrix();
	projCombo->SetSelectedIndex(view.projMode);
	zoomBtn->SetEnabled(view.performXY!=0);
	//UpdateMouseCursor(MODE_SELECT);
	UpdateMouseCursor(mode);


}

// ---------------------------------------------------------------------

void GeometryViewer::SetProjection(int mode) {
	view.projMode = mode;
	projCombo->SetSelectedIndex(mode);
	ToFrontView();
}

// ---------------------------------------------------------------------

void GeometryViewer::SetWorker(Worker *w) {
	work = w;
	ToFrontView();
	// Auto size vector length (consider Front View)
	Geometry *geom = work->GetGeometry();
	AABB bb = geom->GetBB();
	vectorLength = MAX((bb.max.x-bb.min.x) , (bb.max.y-bb.min.y)) / 3.0;
	arrowLength = 10.0 / vectorLength;//MAX((bb.max.z-bb.min.z),vectorLength);
}

// ---------------------------------------------------------------------

void GeometryViewer::DrawIndex() {

	char tmp[256];

	// Draw index number
	if(showIndex || showVertex) {

		// Get selected vertex
		Geometry *geom=work->GetGeometry();
		int nbVertex = geom->GetNbVertex();
		int nbFacet = geom->GetNbFacet();
		if( nbVertex<=0 ) return;

		int *vIdx = (int *)malloc(nbVertex*sizeof(int));
		memset(vIdx,0xFF,nbVertex*sizeof(int));
		for(int i=0;i<nbFacet;i++) {
			Facet *f = geom->GetFacet(i);
			if( f->selected ) {
				int nb = f->sh.nbIndex;
				for(int i=0;i<nb;i++) {
					vIdx[f->indices[i]] = i;
				}
			}
		}

		// Draw dot
		glPointSize(4.0f);
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
		glDisable(GL_BLEND);
		glDisable(GL_CULL_FACE);

		glColor3f(1.0f,0.2f,0.2f);

		glBegin(GL_POINTS);
		for(int i=0;i<nbVertex;i++) {
			if( vIdx[i]>=0 ) {
				VERTEX3D *v = geom->GetVertex(i);
				glVertex3d(v->x,v->y,v->z);
			}
		}
		glEnd();

		// Save contexct
		GLToolkit::DrawStringInit();
		GLToolkit::GetDialogFont()->SetTextColor(0.5f,0.9f,0.9f);

		// Draw Labels
		for(int i=0;i<nbVertex;i++) {
			int idx = vIdx[i];
			if( idx>= 0 ) {
				if(showIndex && showVertex) {
					sprintf(tmp,"%d,%d ",idx+1,i+1);
				} else if (showIndex && !showVertex) {
					sprintf(tmp,"%d ",idx+1);
				} else {
					sprintf(tmp,"%d ",i+1);
				}
				VERTEX3D *v = geom->GetVertex(i);
				GLToolkit::DrawString((float)v->x,(float)v->y,(float)v->z,tmp,GLToolkit::GetDialogFont(),2,2);
			}
		}

		//Restore
		GLToolkit::DrawStringRestore();
		free(vIdx);

	}

}

// ---------------------------------------------------------------------

void GeometryViewer::DrawRule() {

	if( showRule ) {
		// Restore large clipping plane for drawing rules
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		if( view.projMode==PERSPECTIVE_PROJ ) {
			double aspect = (double)width/(double)(height-DOWN_MARGIN);
			GLToolkit::PerspectiveLH(FOV_ANGLE,aspect,0.05,1E6);
		} else {
			glOrtho(view.vLeft,view.vRight,view.vBottom,view.vTop,-1E6,1E6);
		}
		glDisable(GL_DEPTH_TEST);
		GLToolkit::SetMaterial(&greenMaterial);
		GLToolkit::DrawRule(vectorLength,FALSE,FALSE,FALSE,arrowLength);
		GLToolkit::GetDialogFontBold()->SetTextColor(0.4f,0.8f,0.8f);
		GLToolkit::DrawStringInit();
		GLToolkit::DrawString((float)vectorLength,0.0f,0.0f,"x",GLToolkit::GetDialogFontBold());
		GLToolkit::DrawString(0.0f,(float)vectorLength,0.0f,"y",GLToolkit::GetDialogFontBold());
		GLToolkit::DrawString(0.0f,0.0f,(float)vectorLength,"z",GLToolkit::GetDialogFontBold());
		GLToolkit::DrawStringRestore();
	}

}

// ---------------------------------------------------------------------

void GeometryViewer::DrawNormal() {

	if(showNormal) {
		Geometry *geom=work->GetGeometry();
		for(int i=0;i<geom->GetNbFacet();i++) {
			Facet *f = geom->GetFacet(i);
			if( f->selected ) {
				VERTEX3D v1 = geom->GetFacetCenter(i);
				VERTEX3D v2 = f->sh.N;
				GLToolkit::SetMaterial(&blueMaterial);
				glLineWidth(2.0f);
				GLToolkit::DrawVector(v1.x,v1.y,v1.z,v1.x+v2.x*vectorLength,v1.y+v2.y*vectorLength,v1.z+v2.z*vectorLength,arrowLength);
				glPointSize(3.0f);
				glColor3f(1.0f,0.0f,0.0f);
				glBegin(GL_POINTS);
				glVertex3d(v1.x,v1.y,v1.z);
				glEnd();
			}
		}
	}

}

// ---------------------------------------------------------------------

void GeometryViewer::DrawUV() {

	if(showUV) {
		Geometry *geom=work->GetGeometry();
		for(int i=0;i<geom->GetNbFacet();i++) {
			Facet *f = geom->GetFacet(i);
			if( f->selected ) {
				VERTEX3D O = f->sh.O;
				VERTEX3D U = f->sh.U;
				VERTEX3D V = f->sh.V;
				GLToolkit::SetMaterial(&blueMaterial);
				if (antiAliasing) {
					glEnable(GL_LINE_SMOOTH);
					//glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
					//glColor4f(0.0f,0.0f,1.0f,0.5f);
					glEnable(GL_DEPTH_TEST);
				}
				glLineWidth(1.0f);
				GLToolkit::DrawVector(O.x,O.y,O.z,O.x+U.x,O.y+U.y,O.z+U.z,arrowLength);
				GLToolkit::DrawVector(O.x,O.y,O.z,O.x+V.x,O.y+V.y,O.z+V.z,arrowLength);
				if (antiAliasing) glDisable(GL_LINE_SMOOTH);
				glPointSize(3.0f);
				glColor3f(0.5f,1.0f,1.0f);
				glBegin(GL_POINTS);
				glVertex3d(O.x,O.y,O.z);
				glEnd();
				//glEnable(GL_BLEND);
				GLToolkit::GetDialogFont()->SetTextColor(0.5f,0.6f,1.0f);
				GLToolkit::DrawStringInit();
				GLToolkit::DrawString((float)(O.x+U.x),(float)(O.y+U.y),(float)(O.z+U.z),"\201",GLToolkit::GetDialogFont());
				GLToolkit::DrawString((float)(O.x+V.x),(float)(O.y+V.y),(float)(O.z+V.z),"\202",GLToolkit::GetDialogFont());
				GLToolkit::DrawStringRestore();
				//glDisable(GL_BLEND);
			}
		}
	}

}

// ---------------------------------------------------------------------

void GeometryViewer::DrawLeak() {

	// Draw leak
	if(showLeak) {

		// Retreive leak data
		LEAK pLeak[NBHLEAK];
		int nbLeak;
		work->GetLeak(pLeak,&nbLeak);
		if( nbLeak ) {

			glPointSize(4.0f);
			glDisable(GL_TEXTURE_2D);
			glDisable(GL_LIGHTING);
			glDisable(GL_BLEND);
			glDisable(GL_CULL_FACE);
			glEnable(GL_LINE_SMOOTH);
			for(int i=0;i<dispNumLeaks;i++) {

				VERTEX3D p = pLeak[i].pos;
				VERTEX3D d = pLeak[i].dir;

				glColor3f(0.9f,0.2f,0.5f);
				glBegin(GL_POINTS);
				glVertex3d(p.x , p.y , p.z);
				glEnd();

				GLToolkit::DrawVector(p.x , p.y  , p.z ,
					p.x + d.x*vectorLength , p.y + d.y*vectorLength , p.z + d.z*vectorLength , arrowLength);

			}
			glDisable(GL_LINE_SMOOTH);
		}

	}

}

// ---------------------------------------------------------------------

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
				glDisable(GL_BLEND);
				glDisable(GL_CULL_FACE);

				//Determine MAX dF value
				int count=0;
				
				int generation_mode=(work->regions.size()>0)?work->regions[0].generation_mode:SYNGEN_MODE_FLUXWISE;

				double opacityMax=0.0;
				while(count<dispNumHits && pHits[count].type!=0) {
					if (((generation_mode==SYNGEN_MODE_FLUXWISE)?pHits[count].dP:pHits[count].dF)>opacityMax) opacityMax=(generation_mode==SYNGEN_MODE_FLUXWISE)?pHits[count].dP:pHits[count].dF;
					count++;
				}

				count=0;
				

				while(count<dispNumHits && pHits[count].type!=0) {
					//if (count>0&&pHits[count].type==HIT_DES&&pHits[count-1].type!=HIT_ABS) __debugbreak(); //desorbed without being absorbed first
					float lineOpacity;
					
					if (whiteBg) { //whitebg
						lineOpacity=1.0f;
						glColor4f(0.2f,0.7f,0.2f,lineOpacity);
						
					}
					else {
						lineOpacity=(float)((generation_mode==SYNGEN_MODE_FLUXWISE)?pHits[count].dP:pHits[count].dF)/(float)opacityMax;
						glColor4f(0.5f,1.0f,0.5f,lineOpacity);
					}
					if (antiAliasing) {
						glEnable(GL_BLEND);//,glEnable(GL_LINE_SMOOTH);
						//glLineWidth(	2.0f);
					}
					
					glBegin(GL_LINE_STRIP);
					while(count<dispNumHits && pHits[count].type!=HIT_ABS) {
						//teleport routine
						if (pHits[count].type==HIT_TELEPORT) {
							glVertex3d(pHits[count].pos.x , pHits[count].pos.y , pHits[count].pos.z);
							glEnd();
							if (!whiteBg) {
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

							if (whiteBg) { //whitebg
								glColor3f(0.3f,0.9f,0.3f);
							} else {
								glColor3f(0.6f,1.0f,0.6f);
							}
							glBegin(GL_LINE_STRIP);
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

				}
				if (antiAliasing) {
					//glDisable(GL_LINE_SMOOTH);
					glDisable(GL_BLEND);//glLineWidth(	1.0f);
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
				if (whiteBg) { //whitebg
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
				pointSize=(bigDots)?4.0f:2.0f;
				glPointSize(pointSize);
				glColor3f(0.5f,1.0f,1.0f);
				glBegin(GL_POINTS);
				for(int i=0;i<dispNumHits;i++)
					if(pHits[i].type==HIT_TRANS)
						glVertex3d(pHits[i].pos.x , pHits[i].pos.y , pHits[i].pos.z);
				glEnd();

				// Teleport
				pointSize=(bigDots)?4.0f:2.0f;
				glPointSize(pointSize);
				if (!whiteBg) {
					glColor3f(1.0f,0.7f,0.2f);
				} else {
					glColor3f(1.0f,0.0f,1.0f);
				}
				glBegin(GL_POINTS);
				for(int i=0;i<dispNumHits;i++)
					if(pHits[i].type==HIT_TELEPORT)
						glVertex3d(pHits[i].pos.x , pHits[i].pos.y , pHits[i].pos.z);
				glEnd();

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

// ---------------------------------------------------------------------

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

// ---------------------------------------------------------------------

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

// ---------------------------------------------------------------------

void GeometryViewer::ComputeBB(BOOL getAll) {

	Geometry *geom = work->GetGeometry();

	GLMatrix mv; 
	float rx,ry,rz,rw;
	double dx,dy,dz;
	xMin =  1e100;
	yMin =  1e100;
	zNear = 1e100;
	xMax = -1e100;
	yMax = -1e100;
	zFar = -1e100;
	mv.LoadGL(matView);

	//Transform the AABB (fast method, but less accurate 
	//than full vertex transform)
	//AABB bbO = geom->GetBB();
	//TRANSFORMBB(min.x,min.y,min.z);
	//TRANSFORMBB(max.x,min.y,min.z);
	//TRANSFORMBB(max.x,min.y,max.z);
	//TRANSFORMBB(min.x,min.y,max.z);
	//TRANSFORMBB(min.x,max.y,min.z);
	//TRANSFORMBB(max.x,max.y,min.z);
	//TRANSFORMBB(max.x,max.y,max.z);
	//TRANSFORMBB(min.x,max.y,max.z);

	int nbV = geom->GetNbVertex();

	if( geom->viewStruct<0 || getAll ) {

		for(int i=0;i<nbV;i++) {
			VERTEX3D *p = geom->GetVertex(i);
			TRANSFORMVERTEX(p->x,p->y,p->z);
		}
		
		//regions included
		AABB bb=geom->GetBB();
		TRANSFORMVERTEX(bb.min.x,bb.min.y,bb.min.z);
		TRANSFORMVERTEX(bb.max.x,bb.min.y,bb.min.z);
		TRANSFORMVERTEX(bb.min.x,bb.max.y,bb.min.z);
		TRANSFORMVERTEX(bb.min.x,bb.min.y,bb.max.z);
		TRANSFORMVERTEX(bb.min.x,bb.max.y,bb.max.z);
		TRANSFORMVERTEX(bb.max.x,bb.min.y,bb.max.z);
		TRANSFORMVERTEX(bb.max.x,bb.max.y,bb.min.z);
		TRANSFORMVERTEX(bb.max.x,bb.max.y,bb.max.z);

	} else {

		int *refIdx = (int *)malloc( nbV * sizeof(int) );
		memset(refIdx,0,nbV * sizeof(int));

		// Get facet of the selected structure
		int nbF = geom->GetNbFacet();
		for(int i=0;i<nbF;i++) {
			Facet *f = geom->GetFacet(i);
			if( f->sh.superIdx == geom->viewStruct ) {
				for(int j=0;j<f->sh.nbIndex;j++) refIdx[f->indices[j]]=1;
			}
		}

		// Transform vertex
		for(int i=0;i<nbV;i++) {
			if( refIdx[i] ) {
				VERTEX3D *p = geom->GetVertex(i);
				TRANSFORMVERTEX(p->x,p->y,p->z);
			}
		}

		free(refIdx);

	}


}


// ---------------------------------------------------------------------

void GeometryViewer::AutoScale(BOOL reUpdateMouseCursor) {

	if(!work) return;
	Geometry *geom = work->GetGeometry();
	if(!geom) return;

	double aspect = (double)width/(double)(height-DOWN_MARGIN);
	VERTEX3D org = geom->GetCenter();

	// Reset offset, zoom
	view.camOffset.x = 0.0;
	view.camOffset.y = 0.0;
	view.camOffset.z = 0.0;
	view.camDist = 1.0f;
	UpdateMatrix();

	// Get geometry transformed BB
	ComputeBB(FALSE);

	VERTEX3D v;
	v.x = xMax - org.x;
	v.y = yMax - org.y;
	v.z = zFar - org.z;
	camDistInc = Norme(&v) / 100.0;
	view.camOffset.x = 0.0;
	view.camOffset.y = 0.0;
	view.camOffset.z = 0.0;
	if( view.projMode==PERSPECTIVE_PROJ ) {

		// Autoscale dist, at least try to ;) .stub.
		double a = 0.5 / tan((FOV_ANGLE/360.0)*PI);
		view.camDist = MAX((xMax - xMin)/aspect,
			(yMax - yMin)) * a
			+ (zFar - zNear)/1.9;

	} else {

		// Scale
		view.camDist = 1.0;

		double mDist = MAX( (xMax-xMin) , (yMax-yMin)*aspect );
		mDist *= 1.1; // 10% margin
		double dx = (xMax+xMin) / 2.0;
		double dy = (yMax+yMin) / 2.0;
		view.vLeft   =  dx - mDist/2.0;
		view.vRight  =  dx + mDist/2.0;
		view.vTop    =  dy - mDist/(2.0*aspect);
		view.vBottom =  dy + mDist/(2.0*aspect);

	}

	UpdateMatrix();
	zoomBtn->SetEnabled(view.performXY!=XYZ_NONE);
	//if(reUpdateMouseCursor) UpdateMouseCursor(MODE_SELECT);
	if(reUpdateMouseCursor) UpdateMouseCursor(mode);


}

// ---------------------------------------------------------------------

void GeometryViewer::Zoom() {

	if(abs(selX1-selX2)>=2 && abs(selY1-selY2)>=2 && work) {

		int wS = abs(selX1-selX2);
		int hS = abs(selY1-selY2);
		double sAspect = (double)wS/(double)hS;
		double aspect = (double)width / (double)(height-DOWN_MARGIN);
		double x0,y0,w0,h0;
		double dx=0.0,dy=0.0,dz=0.0;
		VERTEX3D org = work->GetGeometry()->GetCenter();
		double z;

		if( hS > wS ) {
			w0 = (double)hS*aspect;
			h0 = (double)hS;
			x0 = ((double)wS - w0)/2.0f + (double)MIN(selX1,selX2) + w0/2.0;
			y0 = (double)MIN(selY1,selY2) + h0/2.0;
			z  = (double)(height-DOWN_MARGIN)/h0;
		} else {
			w0 = (double)wS;
			h0 = (double)wS/aspect;
			x0 = (double)MIN(selX1,selX2) + w0/2.0;
			y0 = ((double)hS - h0)/2.0f + (double)MIN(selY1,selY2) + h0/2.0;
			z  = (double)width/w0;
		}

		x0 -= (double)posX;
		y0 -= (double)posY;

		switch(view.performXY) {
		case XYZ_TOP: // TopView
			dx = (0.5 - x0/(double)width)  * (view.vRight-view.vLeft);
			dz = ( 0.5 - y0/(double)(height-DOWN_MARGIN)) * (view.vBottom-view.vTop);
			break;
		case XYZ_SIDE: // Side View
			dz = ( 0.5 - x0/(double)width) * (view.vRight-view.vLeft);
			dy = ( 0.5 - y0/(double)(height-DOWN_MARGIN)) * (view.vBottom-view.vTop);
			break;
		case XYZ_FRONT: // Front View
			dx = (-0.5 + x0/(double)width)  * (view.vRight-view.vLeft) ;
			dy = ( 0.5 - y0/(double)(height-DOWN_MARGIN)) * (view.vBottom-view.vTop) ;
			break;
		}

		view.camOffset.x -= dx/view.camDist;
		view.camOffset.y += dy/view.camDist;
		view.camOffset.z += dz/view.camDist;
		view.camDist *= z;

		UpdateMatrix();

	}

}

// ---------------------------------------------------------------------

void GeometryViewer::Paint() {

	char tmp[256];

	if(!parent) return;
	GLComponent::Paint();
	//Background gradient
	int x,y,width,height;
	((GLComponent*)this)->GetBounds(&x,&y,&width,&height);
	if (!whiteBg) { glBegin(GL_QUADS);

	glColor3f(0.7f,0.4f,0.3f); //glColor3f(0.3f,0.5f,0.7f); //blue top
	glVertex2i(x,y);
	glVertex2i(x+width,y);
	glColor3f(0.05f,0.05f,0.05f); //grey bottom
	glVertex2i(x+width,y+height);
	glVertex2i(x,y+height);

	glEnd();}

	if(!work) return;
	Geometry *geom=work->GetGeometry();
	if( !geom->IsLoaded() ) {
		PaintCompAndBorder();
		return;
	}

	sprintf(tmp,"");
	topBtn->SetCheck(FALSE);
	frontBtn->SetCheck(FALSE);
	sideBtn->SetCheck(FALSE);
	if( view.performXY ) {
		// Draw coordinates on screen when aligned
		VERTEX3D org = geom->GetCenter();
		double x,y,z;
		switch(view.performXY) {
		case XYZ_TOP: // TopView
			x = -view.vLeft - (1.0-(double)mXOrg/(double)width) * (view.vRight-view.vLeft) + (org.x+view.camOffset.x)*view.camDist;
			z = -view.vTop - ((double)mYOrg/(double)(height-DOWN_MARGIN)) * (view.vBottom-view.vTop) + (org.z+view.camOffset.z)*view.camDist;
			sprintf(tmp,"X=%g, Z=%g",x/view.camDist,z/view.camDist);
			topBtn->SetCheck(TRUE);
			break;
		case XYZ_SIDE: // Side View
			z = -view.vLeft - ((double)mXOrg/(double)width) * (view.vRight-view.vLeft) + (org.z+view.camOffset.z)*view.camDist;
			y = -view.vTop - ((double)mYOrg/(double)(height-DOWN_MARGIN)) * (view.vBottom-view.vTop) + (org.y+view.camOffset.y)*view.camDist;
			sprintf(tmp,"Z=%g, Y=%g",z/view.camDist,y/view.camDist);
			sideBtn->SetCheck(TRUE);
			break;
		case XYZ_FRONT: // Front View
			x = -view.vLeft - (1.0-(double)mXOrg/(double)width) * (view.vRight-view.vLeft) + (org.x+view.camOffset.x)*view.camDist;
			y = -view.vTop - ((double)mYOrg/(double)(height-DOWN_MARGIN)) * (view.vBottom-view.vTop) + (org.y+view.camOffset.y)*view.camDist;
			sprintf(tmp,"X=%g, Y=%g",x/view.camDist,y/view.camDist);
			frontBtn->SetCheck(TRUE);
			break;
		}
	}
	coordLab->SetText(tmp);


	// Clipping and projection matrix
	GetWindow()->Clip(this,0,0,0,DOWN_MARGIN);
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(matProj);
	UpdateLight();
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(matView);
	glDisable(GL_BLEND);

	// Draw geometry
	if( showVolume || showTexture ) glEnable(GL_DEPTH_TEST);
	glClearDepth(1.0);
	glClear(GL_DEPTH_BUFFER_BIT);
	glDepthFunc(GL_LEQUAL);

	/*// Draw geometry
	glEnable(GL_DEPTH_TEST);
	glClearDepth(1.0);
	glClear(GL_DEPTH_BUFFER_BIT);
	glDepthFunc(GL_LEQUAL);*/

	int bgCol=(whiteBg)?255:0;
	SetBackgroundColor(bgCol,bgCol,bgCol);

	DrawLinesAndHits();
	geom->Render((GLfloat *)matView,showVolume,showTexture,showBack,showFilter,showHidden,showMesh,showDir);
	for (int i=0;i<(int)work->regions.size();i++)
		work->regions[i].Render(dispNumTraj,&blueMaterial,vectorLength);

	DrawIndex();
	DrawNormal();
	DrawUV();
	DrawLeak();
	DrawRule();
	geom->PaintSelectedVertices(showHiddenVertex);
	//DrawBB();

	// Restore old transformation/viewport
	GetWindow()->ClipToWindow();
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	GLWindowManager::SetDefault();

	// Draw selection rectangle
	if( (draggMode==DRAGG_SELECT || draggMode==DRAGG_SELECTVERTEX || draggMode==DRAGG_ZOOM ) && (mode==MODE_SELECT || mode==MODE_SELECTVERTEX || mode==MODE_ZOOM) ) {

		GLushort dashPattern = 0xCCCC;

		glDisable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
		glDisable(GL_BLEND);
		glDisable(GL_CULL_FACE);
		glColor3f(1.0f,0.8f,0.9f);
		glEnable(GL_LINE_STIPPLE);
		glLineStipple(1,dashPattern);
		glBegin(GL_LINE_LOOP);
		_glVertex2i(selX1,selY1);
		_glVertex2i(selX1,selY2);
		_glVertex2i(selX2,selY2);
		_glVertex2i(selX2,selY1);
		glEnd();
		glDisable(GL_LINE_STIPPLE);

	}

	PaintCompAndBorder();

	capsLockLabel->SetVisible(GetWindow()->IsCapsLockOn());
}



//-----------------------------------------------------------------------------

void GeometryViewer::PaintCompAndBorder() {

	// Components
	PaintComponents();
	GLToolkit::CheckGLErrors("Geometryviewer::PaintCompandBorder()");

	// Border
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);
	glDisable(GL_CULL_FACE);

	if( selected ) {

		glColor3f(0.5f,0.5f,1.0f);
		glLineWidth(3.0f);
		glBegin(GL_LINES);
		_glVertex2i(posX,posY);
		_glVertex2i(posX+width,posY);
		_glVertex2i(posX,posY);
		_glVertex2i(posX,posY+height);
		_glVertex2i(posX+width,posY+height);
		_glVertex2i(posX,posY+height);
		_glVertex2i(posX+width,posY+height);
		_glVertex2i(posX+width,posY-1);
		glEnd();
		glLineWidth(1.0f);

	} else {

		glColor3f(0.5f,0.5f,0.5f);
		glBegin(GL_LINES);
		_glVertex2i(posX,posY);
		_glVertex2i(posX+width,posY);
		_glVertex2i(posX,posY);
		_glVertex2i(posX,posY+height);
		glEnd();
		glColor3f(1.0f,1.0f,1.0f);
		glBegin(GL_LINES);
		_glVertex2i(posX+width,posY+height);
		_glVertex2i(posX,posY+height);
		_glVertex2i(posX+width,posY+height);
		_glVertex2i(posX+width,posY-1);
		glEnd();

	}

}

//-----------------------------------------------------------------------------

void GeometryViewer::TranslateScale(double diff) {

	// Translate or Scale view according to the projection mode
	if( view.projMode == PERSPECTIVE_PROJ ) {
		view.camDist += diff * camDistInc;
		if(view.camDist<0.01) view.camDist=0.01;
	} else {
		double r = 1.0 - (diff/100.0);
		if(r>0.01) view.camDist *= r;
	}
	UpdateMatrix();

}

//-----------------------------------------------------------------------------

void GeometryViewer::ManageEvent(SDL_Event *evt)
{

	if(!work) return;
	Geometry *geom = work->GetGeometry();
	// Key pressed
	if( evt->type == SDL_KEYDOWN ) {
		int unicode = (evt->key.keysym.unicode & 0x7F);
		if( !unicode ) unicode = evt->key.keysym.sym;

		if( unicode==SDLK_UP ) {
			if( GetWindow()->IsShiftDown() ) {
				view.camAngleOx += angleStep;
				view.performXY=XYZ_NONE;
			} else if( GetWindow()->IsCtrlDown() ) {
				// Up
				view.camOffset.x += transStep * camUp.x;
				view.camOffset.y += transStep * camUp.y;
				view.camOffset.z += transStep * camUp.z;
			} else {
				if( view.projMode == PERSPECTIVE_PROJ ) {
					// Forward
					view.camOffset.x -= transStep * camDir.x;
					view.camOffset.y -= transStep * camDir.y;
					view.camOffset.z -= transStep * camDir.z;
				} else {
					// Up
					view.camOffset.x += transStep * camUp.x;
					view.camOffset.y += transStep * camUp.y;
					view.camOffset.z += transStep * camUp.z;
				}
			}
			UpdateMatrix();
			autoScaleOn=FALSE;
			autoBtn->SetCheck(FALSE);
		}

		if( unicode==SDLK_DOWN ) {
			if( GetWindow()->IsShiftDown() ) {
				view.camAngleOx -= angleStep;
				view.performXY=XYZ_NONE;
			} else if( GetWindow()->IsCtrlDown() ) {
				// Down
				view.camOffset.x -= transStep * camUp.x;
				view.camOffset.y -= transStep * camUp.y;
				view.camOffset.z -= transStep * camUp.z;
			} else {
				if( view.projMode == PERSPECTIVE_PROJ ) {
					// Backward
					view.camOffset.x += transStep * camDir.x;
					view.camOffset.y += transStep * camDir.y;
					view.camOffset.z += transStep * camDir.z;
				} else {
					// Down
					view.camOffset.x -= transStep * camUp.x;
					view.camOffset.y -= transStep * camUp.y;
					view.camOffset.z -= transStep * camUp.z;
				}
			}
			UpdateMatrix();
			autoScaleOn=FALSE;
			autoBtn->SetCheck(FALSE);
		}

		if( unicode==SDLK_LEFT ) {
			if( GetWindow()->IsShiftDown() ) {
				view.camAngleOy += angleStep;
				view.performXY=XYZ_NONE;
			} else {
				// Strafe left
				view.camOffset.x += transStep * camLeft.x;
				view.camOffset.y += transStep * camLeft.y;
				view.camOffset.z += transStep * camLeft.z;
			}
			UpdateMatrix();
		}

		if( unicode==SDLK_RIGHT ) {
			if( GetWindow()->IsShiftDown() ) {
				view.camAngleOy -= angleStep;
				view.performXY=XYZ_NONE;
			} else {
				// Strafe right
				view.camOffset.x -= transStep * camLeft.x;
				view.camOffset.y -= transStep * camLeft.y;
				view.camOffset.z -= transStep * camLeft.z;
			}
			UpdateMatrix();
			autoScaleOn=FALSE;
			autoBtn->SetCheck(FALSE);
		}

		if( unicode == SDLK_LCTRL ||  unicode == SDLK_RCTRL ) {
			//UpdateMouseCursor(MODE_SELECT);
			UpdateMouseCursor(mode);
		}

		if( unicode == SDLK_LSHIFT ||  unicode == SDLK_RSHIFT ) {
			//UpdateMouseCursor(MODE_SELECT);
			UpdateMouseCursor(mode);
		}
		if( unicode == SDLK_LALT ||  unicode == SDLK_RALT ) {
			//UpdateMouseCursor(MODE_SELECT);
			UpdateMouseCursor(mode);
		}

		return;
	}

	if( evt->type == SDL_KEYUP) {

		int unicode = (evt->key.keysym.unicode & 0x7F);
		if( !unicode ) unicode = evt->key.keysym.sym;

		if( unicode == SDLK_LCTRL ||  unicode == SDLK_RCTRL ) {
			//UpdateMouseCursor(MODE_SELECT);
			UpdateMouseCursor(mode);
		}
		if( unicode == SDLK_LSHIFT ||  unicode == SDLK_RSHIFT ) {
			//UpdateMouseCursor(MODE_SELECT);
			UpdateMouseCursor(mode);
		}
		if( unicode == SDLK_LALT ||  unicode == SDLK_RALT ) {
			//UpdateMouseCursor(MODE_SELECT);
			UpdateMouseCursor(mode);
		}

		return;

	}

	if( !draggMode ) {
		GLContainer::ManageEvent(evt);
		GLContainer::RelayEvent(evt);
		if( evtProcessed ) {
			// Latch active cursor
			SetCursor(GLToolkit::GetCursor());
			return;
		} else {
			UpdateMouseCursor(mode);
		}
	}

	// (mX,mY) in window coorinates
	int mX = GetWindow()->GetX(this,evt) + posX;
	int mY = GetWindow()->GetY(this,evt) + posY;

	// Handle mouse events
	if( evt->type == SDL_MOUSEBUTTONDOWN || evt->type == SDL_MOUSEBUTTONDBLCLICK ) {
		mXOrg = mX;
		mYOrg = mY;
		if( evt->button.button == SDL_BUTTON_LEFT ) {
			// Selection dragging
			selX1 = selX2 = mX;
			selY1 = selY2 = mY;
			if (mode==MODE_SELECT) draggMode=DRAGG_SELECT;
			else if (mode==MODE_SELECTVERTEX) draggMode=DRAGG_SELECTVERTEX;
			else if (mode==MODE_SELECTTRAJ) draggMode=DRAGG_SELECTTRAJ;
			else if (mode==MODE_MOVE) draggMode=DRAGG_MOVE;
			else if (mode==MODE_ZOOM) draggMode=DRAGG_ZOOM;
		}
		if( evt->button.button == SDL_BUTTON_MIDDLE ) {
			// Camera translational dragging
			draggMode=DRAGG_MOVE;
			//UpdateMouseCursor(MODE_MOVE);
		}
		if( evt->button.button == SDL_BUTTON_RIGHT ) {
			// Camera rotating
			draggMode=DRAGG_ROTATE;
			//UpdateMouseCursor(MODE_MOVE);
		}
		if( evt->button.button == SDL_BUTTON_WHEELUP ) {
			if (GetWindow()->IsShiftDown()) {
				TranslateScale(-2.0); //Zoom slower when SHIFT is pressed
			}
			else if (GetWindow()->IsCtrlDown()) {
				TranslateScale(-75.0); //Zoom faster when CTRL is pressed
			} else {
				TranslateScale(-20.0);
			}
			autoScaleOn=FALSE;
			autoBtn->SetCheck(FALSE);
		}
		if( evt->button.button == SDL_BUTTON_WHEELDOWN ) {
			if (GetWindow()->IsShiftDown()) {
				TranslateScale(2.0); //Zoom slower when SHIFT is pressed
			}
			else if (GetWindow()->IsCtrlDown()) {
				TranslateScale(75.0); //Zoom slower when SHIFT is pressed
			} else {
				TranslateScale(20.0);
			}
			autoScaleOn=FALSE;
			autoBtn->SetCheck(FALSE);
		}
		UpdateMouseCursor(mode);
	}

	if( evt->type == SDL_MOUSEBUTTONUP ) {
		switch( draggMode ) {

		case DRAGG_SELECT:
		case DRAGG_SELECTVERTEX:
		case DRAGG_ZOOM:
		case DRAGG_SELECTTRAJ:
			switch(mode) {
			case MODE_ZOOM:
				Zoom();
				autoScaleOn=FALSE;
				autoBtn->SetCheck(FALSE);
				break;
			case MODE_SELECT:
				GetWindow()->Clip(this,0,0,0,DOWN_MARGIN);
				glMatrixMode(GL_PROJECTION);
				glLoadMatrixf(matProj);
				glMatrixMode(GL_MODELVIEW);
				glLoadMatrixf(matView);
				selectionChange = TRUE;
				if(abs(selX1-selX2)<=1 && abs(selY1-selY2)<=1) {
					// Simple click, select/unselect facet
					//SetCursor(CURSOR_BUSY);
					GLToolkit::SetCursor(CURSOR_BUSY);
					geom->Select(mX-posX,mY-posY,!GetWindow()->IsShiftDown(),GetWindow()->IsCtrlDown(),GetWindow()->IsCapsLockOn(),this->width,this->height);
					//UpdateMouseCursor(mode);
				} else {
					// Select region
					geom->SelectArea(selX1-posX,selY1-posY,selX2-posX,selY2-posY,
						!GetWindow()->IsShiftDown(),GetWindow()->IsCtrlDown(),GetWindow()->IsCapsLockOn());
				}
				break;
			case MODE_SELECTVERTEX:
				GetWindow()->Clip(this,0,0,0,DOWN_MARGIN);
				glMatrixMode(GL_PROJECTION);
				glLoadMatrixf(matProj);
				glMatrixMode(GL_MODELVIEW);
				glLoadMatrixf(matView);
				//selectionChange = TRUE;
				if(abs(selX1-selX2)<=1 && abs(selY1-selY2)<=1) {
					// Simple click, select/unselect vertex
					geom->SelectVertex(mX-posX,mY-posY, GetWindow()->IsShiftDown(),GetWindow()->IsCtrlDown());
					//select closest vertex
				} else {
					// Select region
					geom->SelectVertex(selX1-posX,selY1-posY,selX2-posX,selY2-posY,
						GetWindow()->IsShiftDown(),GetWindow()->IsCtrlDown());
				}
				break;
			case MODE_SELECTTRAJ:
				GetWindow()->Clip(this,0,0,0,DOWN_MARGIN);
				glMatrixMode(GL_PROJECTION);
				glLoadMatrixf(matProj);
				glMatrixMode(GL_MODELVIEW);
				glLoadMatrixf(matView);
				for (int i=0;i<(int)work->regions.size();i++)
					work->regions[i].SelectTrajPoint(mX-posX,mY-posY);
				break;
			}
			break;
		}

		draggMode=DRAGG_NONE;
		UpdateMouseCursor(mode); //Sets cursor
	}

	if( evt->type == SDL_MOUSEMOTION ) {

		double diffX = (mX - mXOrg);
		double diffY = (mY - mYOrg);
		mXOrg = mX;
		mYOrg = mY;

		UpdateMouseCursor(mode);

		switch( draggMode ) {

		case DRAGG_NONE:
			// performXY
			mXOrg = GetWindow()->GetX(this,evt);
			mYOrg = GetWindow()->GetY(this,evt);
			break;

		case DRAGG_SELECTVERTEX:
		case DRAGG_SELECT:
		case DRAGG_ZOOM:
			// Selection rectangle
			if( GetWindow()->IsAltDown() ) {
				draggMode=DRAGG_MOVE;
			}
			else{
				selX2 = mX;
				selY2 = mY;
			}
			break;

		case DRAGG_MOVE:
			
			if( view.projMode == PERSPECTIVE_PROJ ) {
				double factor=GetWindow()->IsShiftDown()?0.05:1.0;
				double tv = factor*diffX/(double)width * view.camDist * 0.75;
				double tu = factor*diffY/(double)(height-DOWN_MARGIN) * view.camDist * 0.75;
				view.camOffset.x += tu * camUp.x + tv * camLeft.x;
				view.camOffset.y += tu * camUp.y + tv * camLeft.y;
				view.camOffset.z += tu * camUp.z + tv * camLeft.z;
			} else {
				double factor=GetWindow()->IsShiftDown()?0.05:1.0;
				double tv = factor*(diffX/(double)width)  * (view.vRight-view.vLeft)/view.camDist;
				double tu = factor*(diffY/(double)(height-DOWN_MARGIN)) * (view.vBottom-view.vTop)/view.camDist;
				view.camOffset.x += tu * camUp.x + tv * camLeft.x;
				view.camOffset.y += tu * camUp.y + tv * camLeft.y;
				view.camOffset.z += tu * camUp.z + tv * camLeft.z;
			}

			UpdateMatrix();
			autoScaleOn=FALSE;
			autoBtn->SetCheck(FALSE);
			break;

		case DRAGG_ROTATE:

			if( fabs(diffX)>1.9 || fabs(diffY)>1.9 ) {
				if( GetWindow()->IsCtrlDown() ) {
					TranslateScale(diffY);
				}
				else {
					// Rotate view
					double factor=GetWindow()->IsShiftDown()?0.05:1.0;
					if (view.projMode == PERSPECTIVE_PROJ) {
						view.camAngleOx += diffY * angleStep*factor;
						view.camAngleOy += diffX * angleStep*factor;
					} else {
						view.camAngleOx -= diffY * angleStep*factor;
						view.camAngleOy -= diffX * angleStep*factor;
					}
					view.performXY=XYZ_NONE;
					zoomBtn->SetEnabled(FALSE);
					UpdateMatrix();
					if (autoScaleOn) (AutoScale(FALSE));
					//UpdateMouseCursor(MODE_SELECT); //was already MODE_SELECT
					//draggMode=DRAGG_NONE;
				}
			}
			//UpdateMatrix();
			break;
		}
	}

}

// -------------------------------------------------------------
void GeometryViewer::SelectCoplanar(double tolerance) {
	if(!work) return;
	Geometry *geom = work->GetGeometry();
	GetWindow()->Clip(this,0,0,0,DOWN_MARGIN);
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(matProj);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(matView);
	selectionChange = TRUE;
	geom->SelectCoplanar(this->width,this->height,tolerance);
}

// -------------------------------------------------------------

void GeometryViewer::ProcessMessage(GLComponent *src,int message) {

	switch(message) {
	case MSG_BUTTON:
		if( src==topBtn ) {
			ToTopView();
		} else if (src==sideBtn ) {
			ToSideView();
		} else if (src==frontBtn ) {
			ToFrontView();
		} else if (src==zoomBtn) {
			UpdateMouseCursor(MODE_ZOOM);
		} else if (src==sysBtn) {
			GetParent()->ProcessMessage(this,MSG_GEOMVIEWER_MAXIMISE);
		} else if (src==selBtn) {
			UpdateMouseCursor(MODE_SELECT);
		} else if (src==selVxBtn) {
			UpdateMouseCursor(MODE_SELECTVERTEX);
		} else if (src==selTrajBtn) {
			UpdateMouseCursor(MODE_SELECTTRAJ);
		} else if (src==autoBtn) {
			autoScaleOn=!autoScaleOn;
			//autoBtn->SetCheck(autoScaleOn);
			if (autoScaleOn) AutoScale(FALSE);
		} else if (src==handBtn) {
			UpdateMouseCursor(MODE_MOVE);
		}
		break;
	case MSG_COMBO:
		if( src==projCombo ) {
			view.projMode = projCombo->GetSelectedIndex();
			ToFrontView();
		}
		break;
	}

}
