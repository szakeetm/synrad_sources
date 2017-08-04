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

#include "Worker.h"
#include "SynradGeometry.h"
#include "Facet.h"
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "GLApp/GLMatrix.h"
#include "SynRad.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"

extern SynRad *mApp;

void SynradGeometry::BuildFacetTextures(BYTE *hits, bool renderRegularTexture, bool renderDirectionTexture) {

	SHGHITS *shGHit = (SHGHITS *)hits;
	Worker *worker = &(mApp->worker);

	GLProgress *prg = new GLProgress("Building texture", "Frame update");
	prg->SetBounds(5, 28, 300, 90);
	int startTime = SDL_GetTicks();

	if (renderRegularTexture) {
		texCMin_MC = shGHit->minHit_MC;// * dCoef;
		texCMax_MC = shGHit->maxHit_MC;// * dCoef;
		texCMin_flux = shGHit->minHit_flux / worker->no_scans;// * dCoef;
		texCMax_flux = shGHit->maxHit_flux / worker->no_scans;// * dCoef;
		texCMin_power = shGHit->minHit_power / worker->no_scans;// * dCoef;
		texCMax_power = shGHit->maxHit_power / worker->no_scans;// * dCoef;
	}

	for (int i = 0; i < sh.nbFacet; i++) {
		int time = SDL_GetTicks();
		if (!prg->IsVisible() && ((time - startTime) > 500)) {
			prg->SetVisible(true);
		}
		prg->SetProgress((double)i / (double)sh.nbFacet);
		Facet *f = facets[i];

		int profSize = (f->sh.isProfile) ? (PROFILE_SIZE*(sizeof(llong) + 2 * sizeof(double))) : 0;
		//int spectrumSize = (f->sh.hasSpectrum)?(SPECTRUM_SIZE*2*sizeof(double)):0;
		size_t nbElem = f->sh.texWidth*f->sh.texHeight;
		size_t tSize = nbElem * sizeof(double);

		if (renderRegularTexture && f->sh.isTextured) {

			GLint max_t;
			glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max_t);
			if (f->sh.texHeight > max_t || f->sh.texWidth > max_t) {
				if (!f->textureError) {
					char tmp[1024];
					sprintf(tmp, "Facet #%d has a texture of %zdx%zd cells.\n"
						"Your video card only supports texture dimensions (width or height) up to %d cells.\n"
						"Texture rendering has been disabled on this facet, but you can still read texture values\n"
						"using the Texture Plotter window. Consider using a smaller mesh resolution, or split the facet\n"
						"into smaller parts. (Use Facet/Explode... command)", i + 1, f->sh.texHeight, f->sh.texWidth, max_t);
					GLMessageBox::Display(tmp, "OpenGL Error", GLDLG_OK, GLDLG_ICONWARNING);
				}
				f->textureError = true;
				return;
			}
			else {
				f->textureError = false;
			}

		   // Retrieve texture from shared memory (every seconds)
			llong *hits_MC = (llong *)((BYTE *)shGHit + (f->sh.hitOffset + sizeof(SHHITS) + profSize));
			double *hits_flux = (double *)((BYTE *)shGHit + (f->sh.hitOffset + sizeof(SHHITS) + profSize + nbElem * sizeof(llong)));
			double *hits_power = (double *)((BYTE *)shGHit + (f->sh.hitOffset + sizeof(SHHITS) + profSize + nbElem*(sizeof(llong) + sizeof(double))));

			if (texAutoScale) {
				if (textureMode == TEXTURE_MODE_MCHITS) f->BuildTexture(hits_MC, texCMin_MC, texCMax_MC, texColormap, texLogScale);
				else if (textureMode == TEXTURE_MODE_FLUX) f->BuildTexture(hits_flux, texCMin_flux, texCMax_flux, worker->no_scans, texColormap, texLogScale);
				else if (textureMode == TEXTURE_MODE_POWER) f->BuildTexture(hits_power, texCMin_power, texCMax_power, worker->no_scans, texColormap, texLogScale);
			}
			else {
				if (textureMode == TEXTURE_MODE_MCHITS) f->BuildTexture(hits_MC, texMin_MC, texMax_MC, texColormap, texLogScale);
				else if (textureMode == TEXTURE_MODE_FLUX) f->BuildTexture(hits_flux, texMin_flux, texMax_flux, worker->no_scans, texColormap, texLogScale);
				else if (textureMode == TEXTURE_MODE_POWER) f->BuildTexture(hits_power, texMin_power, texMax_power, worker->no_scans, texColormap, texLogScale);
			}
		}
		if (renderDirectionTexture && f->sh.countDirection && f->dirCache) {
			
			size_t dSize = nbElem * sizeof(VHIT);
			
			double iDesorbed = 0.0;
			if (shGHit->total.nbDesorbed)
				iDesorbed = 1.0 / (double)shGHit->total.nbDesorbed;

			VHIT *dirs = (VHIT *)((BYTE *)shGHit + (f->sh.hitOffset + sizeof(SHHITS) + profSize + tSize));
			for (int j = 0; j < nbElem; j++) {
				f->dirCache[j].dir = dirs[j].dir * iDesorbed;
				f->dirCache[j].count = dirs[j].count;
			}
		}

	}
	prg->SetVisible(false);
	SAFE_DELETE(prg);
}