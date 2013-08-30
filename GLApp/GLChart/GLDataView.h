/*
  File:        GLDataView.h
  Description: 2D 'scientific oriented' chart component (Data views managenent)
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

#ifndef _GLCHARTDATAVIEWH_
#define _GLCHARTDATAVIEWH_

typedef struct {

  double x;
  double y;

} APoint2D;

typedef struct {

  APoint2D *pts;
  int length;

} Point2D;

class GLDataView {

public:

  // Construction
  GLDataView();
  ~GLDataView();

  // Dataview methods
  void SetViewType(int s);
  int GetViewType();
  void SetFillStyle(int b);
  int GetFillStyle();
  void SetFillMethod(int m);
  int GetFillMethod();
  void SetFillColor(GLCColor c);
  GLCColor GetFillColor();
  void SetColor(GLCColor c);
  GLCColor GetColor();
  BOOL IsFill();
  void SetFill(BOOL b);
  void SetClickable(BOOL b);
  BOOL IsClickable();
  void SetLabelVisible(BOOL b);
  BOOL IsLabelVisible();
  GLCColor GetLabelColor();
  void SetLabelColor (GLCColor labelColor);
  void SetInterpolationMethod(int method);
  int GetInterpolationMethod();
  void SetInterpolationStep(int step);
  int GetInterpolationStep();
  void SetHermiteTension(double tension);
  double GetHermiteTension();
  void SetHermiteBias(double bias);
  double GetHermiteBias();
  void SetSmoothingMethod(int method);
  int GetSmoothingMethod();
  void SetSmoothingNeighbors(int n);
  int GetSmoothingNeighbors();
  void SetSmoothingGaussSigma(double sigma);
  double GetSmoothingGaussSigma();
  void SetSmoothingExtrapolation(int extMode);
  int GetSmoothingExtrapolation();
  void SetMathFunction(int function);
  int GetMathFunction();
  void SetBarWidth(int w);
  int GetBarWidth();
  void SetMarkerColor(GLCColor c);
  GLCColor GetMarkerColor();
  void SetStyle(int c);
  int GetMarkerSize();
  void SetMarkerSize(int c);
  int GetStyle();
  int GetLineWidth();
  void SetLineWidth(int c);
  void SetName(char *s);
  char *GetName();
  void SetUnit(char *s);
  char *GetUnit();
  char *GetExtendedName();
  int GetMarker();
  void SetMarker(int m);
  double GetA0();
  double GetA1();
  double GetA2();
  void SetA0(double d);
  void SetA1(double d);
  void SetA2(double d);
  BOOL HasTransform();
  BOOL HasFilter();
  void SetAxis(GLAxis *a);
  GLAxis *GetAxis();
  double GetMinimum();
  double GetMaximum();
  double GetMinTime();
  double GetMaxTime();
  double GetMinXValue();
  double GetPositiveMinXValue();
  double GetPositiveMinTime();
  double GetMaxXValue();
  int GetDataLength();
  DataList *GetData();
  void CommitChange();
  void Add(double x, double y);
  void Add(double x, double y,BOOL updateFilter);
  void SetData(double *x,double *y,int nbData);
  int GarbagePointTime(double garbageLimit);
  void GarbagePointLimit(int garbageLimit);
  void ComputeTransformedMinMax(double *tmin,double *tmax);
  double ComputePositiveMin();
  double GetTransformedValue(double x);
  DataList *GetLastValue();
  void Reset();
  double GetYValueByIndex(int idx);
  double GetXValueByIndex(int idx);
  void SetUserFormat(char *format);
  char *GetUserFormat();
  char *FormatValue(double v);
  BOOL IsXDataSorted();
  void SetXDataSorted (BOOL dataSorted);
  APoint2D LinearInterpolate(APoint2D p1,APoint2D p2,double t);
  APoint2D CosineInterpolate(APoint2D p1,APoint2D p2,double t);
  APoint2D CubicInterpolate(APoint2D p0,APoint2D p1,APoint2D p2,APoint2D p3,double t);
  APoint2D HermiteInterpolate(APoint2D p0, APoint2D p1, APoint2D p2, APoint2D p3,double mu, double tension, double bias);

  // Free handle
  int userData1,userData2;

private:

  void FreeList(DataList *l);
  void computeMin();
  void computeMax();
  void computeMinXValue();
  void computeMaxXValue();
  void computeDataBounds();
  Point2D getSource(DataList *src,int nbExtra,BOOL interpNan);
  void mathop();
  int  reverse(int idx,int p);
  void doFFT(Point2D in,int mode);
  void updateFilters();
  void updateSmoothCoefs();
  void convolution();
  void addInt(APoint2D p);
  void addInt(double x,double y);
  Point2D getSegments(DataList *l);
  void interpolate();
  void addPts(Point2D *p,double x,double y);
  void addPts(Point2D *p,APoint2D pt);

  GLAxis *parentAxis;
  GLCColor lineColor;
  GLCColor fillColor;
  GLCColor markerColor;
  int lineStyle;
  int lineWidth;
  int markerType;
  int markerSize;
  int barWidth;
  int fillStyle;
  int fillMethod;
  int type;
  double A0;
  double A1;
  double A2;
  DataList *theData;
  DataList *theFilteredData;
  int dataLength;
  int filteredDataLength;
  DataList *theDataEnd;
  DataList *theFilteredDataEnd;
  double max;
  double min;
  double maxXValue;
  double minXValue;
  char name[64];
  char unit[32];
  BOOL clickable;
  BOOL labelVisible;
  char *userFormat;
  GLCColor labelColor;
  int interpMethod;
  int interpStep;
  double interpTension;
  double interpBias;
  int smoothMethod;
  double *smoothCoefs;
  int smoothNeighbor;
  double smoothSigma;
  int smoothExtrapolation;
  int mathFunction;

};

#endif /* _GLCHARTDATAVIEWH_ */
