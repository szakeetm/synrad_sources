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
#include "Worker.h"
#include "SynradGeometry.h"
#include "Facet_shared.h"
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "GLApp/GLMatrix.h"
#include "SynRad.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"

extern SynRad *mApp;

void SynradGeometry::BuildFacetTextures(BYTE *hits, bool renderRegularTexture, bool renderDirectionTexture) {

	GlobalHitBuffer *shGHit = (GlobalHitBuffer *)hits;
	Worker *worker = &(mApp->worker);

	GLProgress *prg = new GLProgress("Building texture", "Frame update");
	prg->SetBounds(5, 28, 300, 90);
	int startTime = SDL_GetTicks();

	if (renderRegularTexture) {
		textureMin_auto = shGHit->hitMin;// * dCoef;
		textureMax_auto = shGHit->hitMax;// * dCoef;
		textureMin_auto.flux /= worker->no_scans;
		textureMax_auto.flux /= worker->no_scans;
		textureMin_auto.power /= worker->no_scans;
		textureMax_auto.power /= worker->no_scans;
	}

	for (int i = 0; i < sh.nbFacet; i++) {
		int time = SDL_GetTicks();
		if (!prg->IsVisible() && ((time - startTime) > 500)) {
			prg->SetVisible(true);
		}
		prg->SetProgress((double)i / (double)sh.nbFacet);
		Facet *f = facets[i];

		int profSize = (f->sh.isProfile) ? PROFILE_SIZE*sizeof(ProfileSlice) : 0;
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
			TextureCell *texture = (TextureCell *)((BYTE *)shGHit + (f->sh.hitOffset + sizeof(FacetHitBuffer) + profSize));
			f->BuildTexture(texture, textureMode, texAutoScale ? textureMin_auto : textureMin_manual, texAutoScale ? textureMax_auto : textureMax_manual, worker->no_scans, texColormap, texLogScale);
		}
		if (renderDirectionTexture && f->sh.countDirection && f->dirCache) {
			
			size_t dSize = nbElem * sizeof(DirectionCell);
			
			double iDesorbed = 0.0;
			if (shGHit->globalHits.hit.nbDesorbed)
				iDesorbed = 1.0 / (double)shGHit->globalHits.hit.nbDesorbed;

			DirectionCell *dirs = (DirectionCell *)((BYTE *)shGHit + (f->sh.hitOffset + sizeof(FacetHitBuffer) + profSize + tSize));
			for (int j = 0; j < nbElem; j++) {
				f->dirCache[j].dir = dirs[j].dir * iDesorbed;
				f->dirCache[j].count = dirs[j].count;
			}
		}

	}
	prg->SetVisible(false);
	SAFE_DELETE(prg);
}