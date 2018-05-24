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
#include "Facet_shared.h"
#include "SynradTypes.h"
#include "GlApp/MathTools.h" //IS_ZERO
#include "SynradDistributions.h" //Material

// Colormap stuff, defined in GLGradient.cpp
extern std::vector<int> colorMap;

void Facet::LoadGEO(FileReader *file, int version, size_t nbVertex) {

	file->ReadKeyword("indices"); file->ReadKeyword(":");
	for (int i = 0; i < sh.nbIndex; i++) {
		indices[i] = file->ReadInt() - 1;
		if (indices[i] >= nbVertex)
			throw Error(file->MakeError("Facet index out of bounds"));
	}

	file->ReadKeyword("sticking"); file->ReadKeyword(":");
	sh.sticking = file->ReadDouble();
	file->ReadKeyword("opacity"); file->ReadKeyword(":");
	sh.opacity = file->ReadDouble();
	file->ReadKeyword("desorbType"); file->ReadKeyword(":");
	file->ReadInt();
	if (version >= 9) {
		file->ReadKeyword("desorbTypeN"); file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("reflectType"); file->ReadKeyword(":");
	sh.reflectType = file->ReadInt();
	file->ReadKeyword("profileType"); file->ReadKeyword(":");
	file->ReadInt();
	sh.profileType = PROFILE_NONE;
	file->ReadKeyword("superDest"); file->ReadKeyword(":");
	sh.superDest = file->ReadInt();
	file->ReadKeyword("superIdx"); file->ReadKeyword(":");
	sh.superIdx = file->ReadInt();
	file->ReadKeyword("is2sided"); file->ReadKeyword(":");
	sh.is2sided = file->ReadInt();
	if (version < 8) {
		file->ReadKeyword("area"); file->ReadKeyword(":");
		sh.area = file->ReadDouble();
	}
	file->ReadKeyword("mesh"); file->ReadKeyword(":");
	hasMesh = file->ReadInt();
	if (version >= 7) {
		file->ReadKeyword("outgassing"); file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("texDimX"); file->ReadKeyword(":");
	file->ReadDouble();
	sh.texWidthD = 0;
	file->ReadKeyword("texDimY"); file->ReadKeyword(":");
	sh.texHeightD = 0;
	file->ReadDouble();
	file->ReadKeyword("countDes"); file->ReadKeyword(":");
	file->ReadInt();
	file->ReadKeyword("countAbs"); file->ReadKeyword(":");
	sh.countAbs = 0;
	file->ReadInt();
	file->ReadKeyword("countRefl"); file->ReadKeyword(":");
	sh.countRefl = 0;
	file->ReadInt();
	file->ReadKeyword("countTrans"); file->ReadKeyword(":");
	sh.countTrans = 0;
	file->ReadInt();
	file->ReadKeyword("acMode"); file->ReadKeyword(":");
	file->ReadInt();
	file->ReadKeyword("nbAbs"); file->ReadKeyword(":");
	counterCache.nbAbsEquiv = 0.0;
	file->ReadLLong();
	file->ReadKeyword("nbDes"); file->ReadKeyword(":");
	counterCache.nbDesorbed = 0;
	file->ReadLLong();
	file->ReadKeyword("nbHit"); file->ReadKeyword(":");
	counterCache.nbMCHit = 0;
	counterCache.nbHitEquiv = 0.0;
	file->ReadLLong();
	if (version >= 2) {
		// Added in GEO version 2
		file->ReadKeyword("temperature"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("countDirection"); file->ReadKeyword(":");
		sh.countDirection = 0;
		file->ReadInt();
	}
	if (version >= 4) {
		// Added in GEO version 4
		file->ReadKeyword("textureVisible"); file->ReadKeyword(":");
		textureVisible = file->ReadInt();
		file->ReadKeyword("volumeVisible"); file->ReadKeyword(":");
		volumeVisible = file->ReadInt();
	}

	if (version >= 5) {
		// Added in GEO version 5
		file->ReadKeyword("teleportDest"); file->ReadKeyword(":");
		sh.teleportDest = file->ReadInt();
	}

	if (version >= 13) {
		// Added in GEO version 13
		file->ReadKeyword("accomodationFactor"); file->ReadKeyword(":");
		/*sh.accomodationFactor =*/ file->ReadDouble();
	}

	UpdateFlags();
}

void Facet::LoadTXT(FileReader *file) {

	// Opacity parameters descripton (TXT format)
	// -4    => Pressure profile (1 sided)
	// -3    => Desorption distribution
	// -2    => Angular profile
	// -1    => Pressure profile (2 sided)
	// [0,1] => Partial opacity (1 sided)
	// [1,2] => Partial opacity (2 sided)

	// Read facet parameters from TXT format
	sh.sticking = file->ReadDouble();
	double o = file->ReadDouble();
	/*sh.area =*/ file->ReadDouble();
	counterCache.nbDesorbed = 0; file->ReadDouble();
	counterCache.nbMCHit = 0; counterCache.nbHitEquiv = 0.0; file->ReadDouble();
	counterCache.nbAbsEquiv = 0.0; file->ReadDouble();
	file->ReadDouble(); //desorbtype, unused in Synrad
	counterCache.fluxAbs = 0.0; counterCache.powerAbs = 0.0;

	// Convert opacity
	sh.profileType = PROFILE_NONE;
	if (o < 0.0) {

		sh.opacity = 0.0;
		if (IsZero(o + 1.0)) {
			//sh.profileType = REC_PRESSUREU;
			sh.is2sided = true;
		}
		if (IsZero(o + 2.0))
			//sh.profileType = REC_ANGULAR;
			if (IsZero(o + 4.0)) {
				//sh.profileType = REC_PRESSUREU;
				sh.is2sided = false;
			}

	}
	else {
		if (o >= 1.0000001) {
			sh.opacity = o - 1.0;
			sh.is2sided = true;
		}
		else
			sh.opacity = o;
	}

	// Convert desorbType
	/*switch(sh.desorbType) {
	case 0:
	sh.desorbType = DES_COSINE;
	break;
	case 1:
	sh.desorbType = DES_UNIFORM;
	break;
	case 2:
	case 3:
	case 4:
	sh.desorbType = sh.desorbType+1; // cos^n
	break;
	}
	ConvertOldDesorbType();*/
	sh.reflectType = (int)(file->ReadDouble() + 0.5);

	// Convert reflectType
	switch (sh.reflectType) {
	case 0:
		sh.reflectType = REFLECTION_DIFFUSE;
		break;
	case 1:
		sh.reflectType = REFLECTION_SPECULAR;
		break;
	default:
		sh.reflectType = REFLECTION_DIFFUSE;
		break;
	}

	file->ReadDouble(); // Unused

	if (IsTXTLinkFacet())
		sh.superDest = (int)(sh.sticking + 0.5);

	UpdateFlags();

}

void Facet::LoadXML(xml_node f, size_t nbVertex, bool isSynradFile, int vertexOffset) {
	int idx = 0;
	for (xml_node indice : f.child("Indices").children("Indice")) {
		indices[idx] = indice.attribute("vertex").as_int() + vertexOffset;
		if (indices[idx] >= nbVertex) {
			char err[128];
			sprintf(err, "Facet %d refers to vertex %d which doesn't exist", f.attribute("id").as_int() + 1, idx + 1);
			throw Error(err);
		}
		idx++;
	}
	sh.opacity = f.child("Opacity").attribute("constValue").as_double();
	sh.is2sided = f.child("Opacity").attribute("is2sided").as_int();
	sh.superIdx = f.child("Structure").attribute("inStructure").as_int();
	sh.superDest = f.child("Structure").attribute("linksTo").as_int();
	sh.teleportDest = f.child("Teleport").attribute("target").as_int();

	if (isSynradFile) {
		sh.sticking = f.child("Sticking").attribute("constValue").as_double();
		xml_node recNode = f.child("Recordings");
		sh.profileType = recNode.child("Profile").attribute("type").as_int();
		sh.recordSpectrum = recNode.child("Spectrum").attribute("record").as_bool();
		xml_node texNode = recNode.child("Texture");
		hasMesh = texNode.attribute("hasMesh").as_bool();
		sh.texWidthD = texNode.attribute("texDimX").as_double();
		sh.texHeightD = texNode.attribute("texDimY").as_double();
		sh.countAbs = texNode.attribute("countAbs").as_int();
		sh.countRefl = texNode.attribute("countRefl").as_int();
		sh.countTrans = texNode.attribute("countTrans").as_int();
		sh.countDirection = texNode.attribute("countDir").as_int();
	} //else use default values at Facet() constructor

	textureVisible = f.child("ViewSettings").attribute("textureVisible").as_int();
	volumeVisible = f.child("ViewSettings").attribute("volumeVisible").as_int();

	UpdateFlags();
}

void Facet::SaveTXT(FileWriter *file) {

	if (!sh.superDest)
		file->Write(sh.sticking, "\n");
	else {
		file->Write((double)sh.superDest, "\n");
		sh.opacity = 0.0;
	}

	if (sh.is2sided)
		file->Write(sh.opacity + 1.0, "\n");
	else
		file->Write(sh.opacity, "\n");

	file->Write(sh.area, "\n");

	file->Write(0.0, "\n"); //no desorption
	file->Write(0.0, "\n"); //nbHit
	file->Write(0.0, "\n"); //nbAbsorbed

	file->Write(0.0, "\n"); //no desorption

	switch (sh.reflectType) {
	case REFLECTION_DIFFUSE:
		file->Write(0.0, "\n");
		break;
	case REFLECTION_SPECULAR:
		file->Write(1.0, "\n");
		break;
	default:
		file->Write((double)(sh.reflectType), "\n");
		break;
	}

	file->Write(0.0, "\n"); // Unused
}

/*
void Facet::SaveGEO(FileWriter *file,int idx) {

char tmp[256];

sprintf(tmp,"facet %d {\n",idx+1);
file->Write(tmp);
file->Write("  nbIndex:");file->Write(sh.nbIndex,"\n");
file->Write("  indices:\n");
for(int i=0;i<sh.nbIndex;i++) {
file->Write("    ");
file->Write(indices[i]+1,"\n");
}
//file->Write("\n");
file->Write("  sticking:");file->Write(sh.sticking,"\n");
file->Write("  opacity:");file->Write(sh.opacity,"\n");
file->Write("  desorbType:");file->Write(0,"\n");
file->Write("  desorbTypeN:");file->Write(0.0,"\n");
file->Write("  reflectType:");file->Write(sh.reflectType,"\n");
file->Write("  profileType:");file->Write(REC_NONE,"\n");
file->Write("  superDest:");file->Write(sh.superDest,"\n");
file->Write("  superIdx:");file->Write(sh.superIdx,"\n");
file->Write("  is2sided:");file->Write(sh.is2sided,"\n");
//file->Write("  area:");file->Write(sh.area,"\n");
file->Write("  mesh:");file->Write( false ,"\n");
file->Write("  outgassing:");file->Write(0.0 ,"\n");
file->Write("  texDimX:");file->Write(0,"\n");
file->Write("  texDimY:");file->Write(0,"\n");
file->Write("  countDes:");file->Write(0,"\n");
file->Write("  countAbs:");file->Write(sh.countAbs,"\n");
file->Write("  countRefl:");file->Write(sh.countRefl,"\n");
file->Write("  countTrans:");file->Write(sh.countTrans,"\n");
file->Write("  acMode:");file->Write(0,"\n");
file->Write("  nbAbs:");file->Write(0,"\n");
file->Write("  nbDes:");file->Write(0,"\n");
file->Write("  nbHit:");file->Write(0,"\n");

// Version 2
file->Write("  temperature:");file->Write(293.15,"\n");
file->Write("  countDirection:");file->Write(sh.countDirection,"\n");

// Version 4
file->Write("  textureVisible:");file->Write(textureVisible,"\n");
file->Write("  volumeVisible:");file->Write(volumeVisible,"\n");

// Version 5
file->Write("  teleportDest:");file->Write(sh.teleportDest,"\n");

file->Write("}\n");
}
*/

size_t Facet::GetGeometrySize() {

	size_t s = sizeof(FacetProperties)
		+ (sh.nbIndex * sizeof(size_t)) //indices
		+ (sh.nbIndex * sizeof(Vector2d));

	// Size of the 'element area' array passed to the geometry buffer
	if (sh.isTextured) s += sizeof(double)*sh.texWidth*sh.texHeight;
	return s;

}

size_t Facet::GetHitsSize() {

	return   sizeof(FacetHitBuffer)
		+ (sh.isProfile ? PROFILE_SIZE*sizeof(ProfileSlice) : 0)
		+ sh.texWidth*sh.texHeight*sizeof(TextureCell)
		+ (sh.countDirection ? sh.texWidth*sh.texHeight*sizeof(DirectionCell) : 0)
		+ (sh.recordSpectrum ? SPECTRUM_SIZE * sizeof(ProfileSlice) : 0);

}

size_t Facet::GetTexRamSize() {
	//Values
	size_t sizePerCell = sizeof(TextureCell); //hits + flux + power
	if (sh.countDirection) sizePerCell += sizeof(DirectionCell); //DirectionCell: Vector3d + long
	//Mesh
	sizePerCell += sizeof(int); //CellPropertiesIds
	size_t sizePerMeshElement = sizeof(CellProperties);
	sizePerMeshElement += 4 * sizeof(Vector2d); //Estimate: most mesh elements have 4 points
	return sh.texWidth*sh.texHeight*sizePerCell + meshvectorsize*sizePerMeshElement;
}

size_t Facet::GetTexRamSizeForRatio(double ratio, bool useMesh, bool countDir) {

	double nU = sh.U.Norme();
	double nV = sh.V.Norme();
	double width = nU*ratio;
	double height = nV*ratio;

	bool dimOK = (width*height > 0.0000001);

	if (dimOK) {
		int iwidth = (int)ceil(width);
		int iheight = (int)ceil(height);
		//Values
		size_t sizePerCell = sizeof(TextureCell); //hits + flux + power
		if (sh.countDirection) sizePerCell += sizeof(DirectionCell); //DirectionCell: Vector3d + long
		//Mesh
		sizePerCell += sizeof(int); //CellPropertiesIds
		size_t sizePerMeshElement = sizeof(CellProperties);
		sizePerMeshElement += 4 * sizeof(Vector2d); //Estimate: most mesh elements have 4 points
		return iwidth*iheight*(sizePerCell + sizePerMeshElement); //Conservative: assuming all cells are non-full

	}
	else {

		return 0;

	}

}

#define LOG10(x) log10f((float)x)

void  Facet::BuildTexture(TextureCell *texture, const size_t& textureMode, const TextureCell& minVal, const TextureCell& maxVal, const double& no_scans, const bool& useColorMap, bool doLog, const bool& normalize) {
	double min, max;
	if (textureMode == TEXTURE_MODE_MCHITS) {
		min = (double)minVal.count;
		max = (double)maxVal.count;
	}
	else if (textureMode == TEXTURE_MODE_FLUX) {
		min = minVal.flux * no_scans;
		max = maxVal.flux * no_scans;
	}
	else if (textureMode == TEXTURE_MODE_POWER) {
		min = minVal.power * no_scans;
		max = maxVal.power * no_scans;
	}
	
	size_t size = sh.texWidth*sh.texHeight;
	size_t tSize = texDimW*texDimH;
	if (size == 0 || tSize == 0) return;

	double scaleFactor = 1.0;
	int val;

	glBindTexture(GL_TEXTURE_2D, glTex);
		
		// Scale
		if (min < max) {
			if (doLog) {
				if (min < 1e-20) min = 1e-20;
				scaleFactor = (useColorMap ? 65534.0 : 255.0) / (log10(max) - log10(min)); // -1 for saturation color
			}
			else {
				scaleFactor = (useColorMap ? 65534.0 : 255.0) / (max - min); // -1 for saturation color
			}
		}
		else {
			doLog = false;
			min = 0.0;
		}

		unsigned char *buff8;
		int *buff32;

		if (useColorMap) {
			buff32 = (int *)calloc(tSize,sizeof(int)); //Color
			if (!buff32) throw Error("Cannot allocate memory for texture buffer");
		}
		else {
			buff8 = (unsigned char *)calloc(tSize, sizeof(unsigned char)); //Greyscale
			if (!buff8) throw Error("Cannot allocate memory for texture buffer");
		}		

		for (size_t j = 0; j < sh.texHeight; j++) {
			for (size_t i = 0; i < sh.texWidth; i++) {
				size_t idx = i + j*sh.texWidth;
				double cellValue;
				if (textureMode == TEXTURE_MODE_MCHITS)
					cellValue = (double)texture[idx].count;
				else if (textureMode == TEXTURE_MODE_FLUX)
					cellValue = texture[idx].flux;
				else if (textureMode == TEXTURE_MODE_POWER)
					cellValue = texture[idx].power;
				if (doLog) {
					val = (int)((log10(cellValue) - log10(min))*scaleFactor + 0.5);
				}
				else {
					val = (int)((cellValue - min)*scaleFactor + 0.5);
				}
				Saturate(val, 0, useColorMap?65535:255);
				if (useColorMap) {
					buff32[(i + 1) + (j + 1)*texDimW] = colorMap[val];
					if (texture[idx].count == 0) buff32[(i + 1) + (j + 1)*texDimW] = (COLORREF)(65535 + 256 + 1); //show unset value as white
				}
				else {
					buff8[(i + 1) + (j + 1)*texDimW] = val;
				}
			}
		}
		// Perform edge smoothing (only with mesh)
		if (cellPropertiesIds) {

			for (int j = -1; j <= sh.texHeight; j++) {
				for (int i = -1; i <= sh.texWidth; i++) {
					bool doSmooth = (i < 0) || (i >= sh.texWidth) ||
						(j < 0) || (j >= sh.texHeight) ||
						GetMeshArea(i + j*sh.texWidth) == 0.0f;
					if (doSmooth) {
						if (doLog) {
							val = (int)((log10(GetSmooth(i, j, texture, textureMode, 1.0)) - log10(min))*scaleFactor + 0.5);
						}
						else {
							val = (int)((GetSmooth(i, j, texture, textureMode, 1.0f) - min)*scaleFactor + 0.5);
						}
						Saturate(val, 0, useColorMap ? 65535 : 255);
						if (useColorMap)
							buff32[(i + 1) + (j + 1)*texDimW] = colorMap[val];
						else
							buff8[(i + 1) + (j + 1)*texDimW] = val;
					}
				}
			}
		}

		GLint width, height, format;
		glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &width);
		glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &height);
		glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, &format);
		if (format == (useColorMap?GL_RGBA:GL_LUMINANCE) && width == texDimW && height == texDimH) {
			//Update texture
			glTexSubImage2D(
				GL_TEXTURE_2D,       // Type
				0,                   // No Mipmap
				0,					// X offset
				0,					// Y offset
				(int)texDimW,             // Width
				(int)texDimH,             // Height
				(useColorMap ? GL_RGBA : GL_LUMINANCE),             // Format RGBA
				GL_UNSIGNED_BYTE,    // 8 Bit/pixel
				(useColorMap ? (void*)buff32 : (void*)buff8)              // Data
			);
		}
		else {
			//Rebuild texture
			glTexImage2D(
				GL_TEXTURE_2D,       // Type
				0,                   // No Mipmap
				(useColorMap ? GL_RGBA : GL_LUMINANCE),             // Format RGBA or LUMINANCE
				(int)texDimW,             // Width
				(int)texDimH,             // Height
				0,                   // Border
				(useColorMap ? GL_RGBA : GL_LUMINANCE),             // Format RGBA or LUMINANCE
				GL_UNSIGNED_BYTE,    // 8 Bit/pixel
				(useColorMap ? (void*)buff32 : (void*)buff8)              // Data
			);
		}
		free(useColorMap ? (void*)buff32 : (void*)buff8);
		GLToolkit::CheckGLErrors("Facet::BuildTexture()");
}

inline void Facet::Weigh_Neighbor(const size_t & i, const size_t & j, const double & weight, TextureCell * texture, const size_t & textureMode, const float & scaleF, double & weighedSum, double & totalWeigh) {
	if (i >= 0 && i < sh.texWidth && j >= 0 && j < sh.texHeight) {
		size_t index = i + j * sh.texWidth;
		if (GetMeshArea(index) > 0.0) {
			double cellValue;
			if (textureMode == TEXTURE_MODE_MCHITS) {
				cellValue = (double)texture[index].count;
			}
			else if (textureMode == TEXTURE_MODE_FLUX) {
				cellValue = (double)texture[index].flux;
			}
			else if (textureMode == TEXTURE_MODE_POWER) {
				cellValue = (double)texture[index].power;
			}
			weighedSum += weight * (cellValue / GetMeshArea(index)*scaleF);
			totalWeigh += weight;
		}
	}
}

inline double Facet::GetSmooth(const int & i, const int & j, TextureCell * texture, const size_t & textureMode, const float & scaleF) {
	double totalWeigh = 0.0;
	double weighedSum = 0.0;
	if (sh.texWidth > 0 && sh.texHeight > 0) {
		Weigh_Neighbor(i - 1, j - 1, 1.0, texture, textureMode, scaleF, weighedSum, totalWeigh);
		Weigh_Neighbor(i - 1, j + 1, 1.0, texture, textureMode, scaleF, weighedSum, totalWeigh);
		Weigh_Neighbor(i + 1, j - 1, 1.0, texture, textureMode, scaleF, weighedSum, totalWeigh);
		Weigh_Neighbor(i + 1, j + 1, 1.0, texture, textureMode, scaleF, weighedSum, totalWeigh);
		Weigh_Neighbor(i, j - 1, 2.0, texture, textureMode, scaleF, weighedSum, totalWeigh);
		Weigh_Neighbor(i, j + 1, 2.0, texture, textureMode, scaleF, weighedSum, totalWeigh);
		Weigh_Neighbor(i - 1, j, 2.0, texture, textureMode, scaleF, weighedSum, totalWeigh);
		Weigh_Neighbor(i + 1, j, 2.0, texture, textureMode, scaleF, weighedSum, totalWeigh);
	}
	if (totalWeigh == 0.0)
		return 0.0;
	else
		return weighedSum / totalWeigh;
}

void Facet::SaveSYN(FileWriter *file, const std::vector<Material> &materials, int idx, bool crashSave) {

	char tmp[256];

	sprintf(tmp, "facet %d {\n", idx + 1);
	file->Write(tmp);
	file->Write("  nbIndex:"); file->Write(sh.nbIndex, "\n");
	file->Write("  indices:\n");
	for (int i = 0; i < sh.nbIndex; i++) {
		file->Write("    ");
		file->Write(indices[i] + 1, "\n");
	}
	file->Write("  reflectType:"); file->Write(sh.reflectType, "\n");
	file->Write("  sticking:"); file->Write(sh.sticking, "\n");

	if (sh.reflectType >= 10) { //material reflection
		file->Write("  materialName:"); file->Write(materials[sh.reflectType - 10].name + "\n");
	}

	file->Write("  doScattering:"); file->Write(sh.doScattering, "\n");
	file->Write("  rmsRoughness:"); file->Write(sh.rmsRoughness, "\n");
	file->Write("  autoCorrLength:"); file->Write(sh.autoCorrLength, "\n");

	file->Write("  opacity:"); file->Write(sh.opacity, "\n");
	file->Write("  profileType:"); file->Write(sh.profileType, "\n");
	file->Write("  hasSpectrum:"); file->Write(sh.recordSpectrum, "\n");
	file->Write("  superDest:"); file->Write(sh.superDest, "\n");
	file->Write("  superIdx:"); file->Write(sh.superIdx, "\n");
	file->Write("  is2sided:"); file->Write(sh.is2sided, "\n");
	file->Write("  mesh:"); file->Write((cellPropertiesIds != NULL), "\n");
	file->Write("  texDimX:"); file->Write(sh.texWidthD, "\n");
	file->Write("  texDimY:"); file->Write(sh.texHeightD, "\n");
	file->Write("  countAbs:"); file->Write(sh.countAbs, "\n");
	file->Write("  countRefl:"); file->Write(sh.countRefl, "\n");
	file->Write("  countTrans:"); file->Write(sh.countTrans, "\n");
	file->Write("  nbAbsEquiv:"); file->Write(counterCache.nbAbsEquiv, "\n");
	file->Write("  nbHit:"); file->Write(counterCache.nbMCHit, "\n");
	file->Write("  nbHitEquiv:"); file->Write(counterCache.nbHitEquiv, "\n");
	file->Write("  fluxAbs:"); file->Write(counterCache.fluxAbs, "\n");
	file->Write("  powerAbs:"); file->Write(counterCache.powerAbs, "\n");
	file->Write("  countDirection:"); file->Write(sh.countDirection, "\n");
	file->Write("  textureVisible:"); file->Write(textureVisible, "\n");
	file->Write("  volumeVisible:"); file->Write(volumeVisible, "\n");
	file->Write("  teleportDest:"); file->Write(sh.teleportDest, "\n");

	file->Write("}\n");
}
void Facet::LoadSYN(FileReader *file, const std::vector<Material> &materials, int version, size_t nbVertex) {

	file->ReadKeyword("indices"); file->ReadKeyword(":");
	for (int i = 0; i < sh.nbIndex; i++) {
		indices[i] = file->ReadInt() - 1;
		if (indices[i] >= nbVertex)
			throw Error(file->MakeError("Facet index out of bounds"));
	}

	if (version >= 9) { //new reflection model
		file->ReadKeyword("reflectType"); file->ReadKeyword(":");
		sh.reflectType = file->ReadInt();
		file->ReadKeyword("sticking"); file->ReadKeyword(":");
		sh.sticking = file->ReadDouble();

		if (sh.reflectType >= 10) { //Material reflection: update index from the material's name
			file->ReadKeyword("materialName"); file->ReadKeyword(":"); std::string matName = file->ReadWord();
			sh.reflectType = 9; //Invalid material, unless we find it below
			for (size_t i = 0; i < materials.size(); i++) {
				if (materials[i].name == matName) {
					sh.reflectType = 10 + (int)i;
					break;
				}
			}
		}
		file->ReadKeyword("doScattering"); file->ReadKeyword(":");
		sh.doScattering = file->ReadInt();
		file->ReadKeyword("rmsRoughness"); file->ReadKeyword(":");
		sh.rmsRoughness = file->ReadDouble();
		file->ReadKeyword("autoCorrLength"); file->ReadKeyword(":");
		sh.autoCorrLength = file->ReadDouble();
		file->ReadKeyword("opacity"); file->ReadKeyword(":");
		sh.opacity = file->ReadDouble();
	}
	else { //legacy reflection model
		file->ReadKeyword("sticking"); file->ReadKeyword(":");
		sh.sticking = file->ReadDouble();
		if (version >= 4) {
			file->ReadKeyword("roughness"); file->ReadKeyword(":");
			sh.autoCorrLength = 10000.0E-9; //We need a number here. Will take 10 microns so specular reflection probability is low
			sh.rmsRoughness = file->ReadDouble()*sh.autoCorrLength; // 1/Tau approximately corresponds to roughness ratio in Synrad 1.3-
		}
		file->ReadKeyword("opacity"); file->ReadKeyword(":");
		sh.opacity = file->ReadDouble();
		file->ReadKeyword("reflectType"); file->ReadKeyword(":");
		sh.reflectType = file->ReadInt();
		if (sh.reflectType >= 2) { //material reflection
			if (sh.reflectType - 2 < materials.size()) { //material exists
				sh.reflectType += 8; //Material reflections now run from 10 upwards
			}
			else { //invalid material
				sh.reflectType = 9; //invalid material
			}
			sh.doScattering = true;
		}
		else { //diffuse or mirror reflection
			sh.doScattering = false;
		}
	}
	file->ReadKeyword("profileType"); file->ReadKeyword(":");

	sh.profileType = file->ReadInt();
	file->ReadKeyword("hasSpectrum"); file->ReadKeyword(":");
	sh.recordSpectrum = file->ReadInt();
	file->ReadKeyword("superDest"); file->ReadKeyword(":");
	sh.superDest = file->ReadInt();
	file->ReadKeyword("superIdx"); file->ReadKeyword(":");
	sh.superIdx = file->ReadInt();
	file->ReadKeyword("is2sided"); file->ReadKeyword(":");
	sh.is2sided = file->ReadInt();
	file->ReadKeyword("mesh"); file->ReadKeyword(":");
	hasMesh = file->ReadInt();
	file->ReadKeyword("texDimX"); file->ReadKeyword(":");
	sh.texWidthD = file->ReadDouble();
	file->ReadKeyword("texDimY"); file->ReadKeyword(":");
	sh.texHeightD = file->ReadDouble();
	if (version < 3) {
		file->ReadKeyword("countDes"); file->ReadKeyword(":");
		file->ReadInt();
	}
	file->ReadKeyword("countAbs"); file->ReadKeyword(":");
	sh.countAbs = file->ReadInt();
	file->ReadKeyword("countRefl"); file->ReadKeyword(":");
	sh.countRefl = file->ReadInt();
	file->ReadKeyword("countTrans"); file->ReadKeyword(":");
	sh.countTrans = file->ReadInt();
	if (version>=10) file->ReadKeyword("nbAbsEquiv");
	else file->ReadKeyword("nbAbs"); file->ReadKeyword(":");
	counterCache.nbAbsEquiv = file->ReadDouble();
	if (version < 3) {
		file->ReadKeyword("nbDes"); file->ReadKeyword(":");
		counterCache.nbDesorbed = file->ReadLLong();
	}
	file->ReadKeyword("nbHit"); file->ReadKeyword(":");
	counterCache.nbMCHit = file->ReadLLong();
	if (version >= 10) {
		file->ReadKeyword("nbHitEquiv"); file->ReadKeyword(":");
		counterCache.nbHitEquiv = file->ReadDouble();
	}
	else counterCache.nbHitEquiv = static_cast<double>(counterCache.nbMCHit);

	if (version >= 3) {
		file->ReadKeyword("fluxAbs"); file->ReadKeyword(":");
		counterCache.fluxAbs = file->ReadDouble();
		file->ReadKeyword("powerAbs"); file->ReadKeyword(":");
		counterCache.powerAbs = file->ReadDouble();
	}
	file->ReadKeyword("countDirection"); file->ReadKeyword(":");
	sh.countDirection = file->ReadInt();
	file->ReadKeyword("textureVisible"); file->ReadKeyword(":");
	textureVisible = file->ReadInt();
	file->ReadKeyword("volumeVisible"); file->ReadKeyword(":");
	volumeVisible = file->ReadInt();
	file->ReadKeyword("teleportDest"); file->ReadKeyword(":");
	sh.teleportDest = file->ReadInt();

	UpdateFlags();

}

