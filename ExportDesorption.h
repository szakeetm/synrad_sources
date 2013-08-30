/*
  File:        ExportDesorption.h
  Description: Move facet by offset dialog
*/

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"

#include "Geometry.h"
#include "Worker.h"

#ifndef _EXPORTDESH_
#define _EXPORTDESH_

class ExportDesorption : public GLWindow {

public:

  // Construction
  ExportDesorption(Geometry *geom,Worker *work);

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  int LoadConvFile(char *fileName);

private:

  Geometry     *geom;
  Worker	   *work;

  GLButton    *exportButton;
  GLButton    *browseButton;
  GLButton    *cancelButton;
  GLLabel     *fileNameLabel;
  GLLabel     *fileInfoLabel;
  GLTextField *eta0Field;
  GLTextField *alphaField;
  GLToggle *toggle1,*toggle2,*toggle3,*selectedToggle;
  BOOL fileLoaded;

  Distribution2D *conversionDistr;
  double eta0,alpha;
  char* fileName;

  int nbFacetS,mode;
};

#endif /* _ExportDesorptionH_ */
