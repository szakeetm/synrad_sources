/*
  File:        GLChart.h
  Description: 2D 'scientific oriented' chart component
               C++ port of fr.esrf.tangoatk.widget.util.chart
  Author:      J-L PONS (2007)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "GLWindow.h"
#include "GLComponent.h"
#include "GLMenu.h"
#include "GLButton.h"

#ifndef _GLCHARTH_
#define _GLCHARTH_

class GLChart;
class GLAxis;
class SearchInfo;

#include "GLChart/GLChartConst.h"
#include "GLChart/DataList.h"
#include "GLChart/GLDataView.h"
#include "GLChart/GLAxis.h"
#include "GLChart/SearchInfo.h"
#include "GLChart/GLChartOptions.h"
#include "GLChart/GLDataViewOptions.h"

typedef struct {

  GLCRectangle  rect;
  GLDataView   *view;

} LabelRect;

class GLChart : public GLComponent {

public:

  // Construction
  GLChart(int compId);
  ~GLChart();

  // Components method
  GLAxis *GetXAxis();
  GLAxis *GetY1Axis();
  GLAxis *GetY2Axis();
  void SetMargin(GLCDimension d);
  GLCDimension GetMargin();
  void SetBackground(GLCColor c);
  GLCColor GetBackground();
  void SetChartBackground(GLCColor c);
  GLCColor GetChartBackground();
  void SetPaintAxisFirst(BOOL b);
  BOOL IsPaintAxisFirst();
  void SetHeaderVisible(BOOL b);
  void SetHeader(char *s);
  char *GetHeader();
  void SetDisplayDuration(double v);
  double GetDisplayDuration();
  void SetHeaderColor(GLCColor c);
  void SetLabelVisible(BOOL b);
  BOOL IsLabelVisible();
  void SetLabelPlacement(int p);
  int GetLabelPlacement();
  void ShowOptionDialog();
  void ShowDataOptionDialog(GLDataView *v);
  void SaveFile();
  void LoadFile();
  BOOL IsZoomed();
  void EnterZoom();
  void ExitZoom();
  void UnselectAll();
  void AddData(GLDataView *v, double x, double y);
  void RemoveDataView(GLDataView *view);
  void CopyAllToClipboard();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void Paint();
  void ManageEvent(SDL_Event *evt);
  void InvalidateDeviceObjects();
  void RestoreDeviceObjects();

private:

  void paintLabel(GLDataView *v,GLAxis *axis,int x,int y,int w);
  void paintLabelAndHeader();
  void measureGraphItems(int w,int h,GLDataView **views,int nbView);
  void paintZoomButton(int x,int y);
  void paintZoomSelection(int x,int y);
  GLCRectangle buildRect(int x1, int y1, int x2, int y2);
  void mouseDragged(SDL_Event *evt);
  void mouseReleased(SDL_Event *evt);
  void mousePressed(SDL_Event *evt);
  void redrawPanel();
  char **buildPanelString(SearchInfo *si);
  void showPanel(SearchInfo *si);
  int  garbageData(GLDataView *v);
  void setRect(GLCRectangle *r,int x,int y,int w,int h);
  void paintChartBackground(int xOrg,int yOrg);
  BOOL insideRect(LabelRect *r,int x,int y);
  void showChartMenu(int x,int y);

  char header[256];
  BOOL headerVisible;
  GLCColor headerColor;

  BOOL ipanelVisible;
  BOOL paintAxisFirst;
  GLCColor chartBackground;

  double displayDuration;

  GLMenu *chartMenu;
  GLMenu *dvMenuY1;
  GLMenu *dvMenuY2;

  BOOL zoomDrag;
  BOOL zoomDragAllowed;
  int zoomX;
  int zoomY;
  GLButton *zoomButton;

  int lastX;
  int lastY;
  SearchInfo lastSearch;

  // Measurements stuff
  GLCRectangle headerR;
  GLCRectangle labelR;
  GLCRectangle viewR;
  GLCDimension margin;

  // Label placement
  LabelRect labelRect[MAX_VIEWS];
  BOOL labelVisible;
  int  labelMode;
  int labelHeight;
  int labelWidth;
  int labelPerLine;
  int nbLabel;
  int headerWidth;
  int axisHeight;
  int axisWidth;
  int y1AxisThickness;
  int y2AxisThickness;
  int xAxisThickness;
  int xAxisUpMargin;

  // Axis
  GLAxis *xAxis;
  GLAxis *y1Axis;
  GLAxis *y2Axis;

  // Config dialog
  GLChartOptions *chartOptions;
  GLDataViewOptions *dvOptions;

  // File name
  char lastDir[512];
  char lastFile[64];

};

#endif /* _GLCHARTH_ */
