#ifndef _TOOLS_
#define _TOOLS_

#include "Distributions.h"
typedef  unsigned __int64 LONGINT;

#define PI 3.14159265358979323846
#define proton_charge 1.602E-19 //in Coulombs
#define proton_mass 1.67E-27 //in kg

class Vector {
public:
	double x;
	double y;
	double z;
	double Norme();
	Vector();
	Vector(const Vector &src);
	Vector& operator=(const Vector &src);
	Vector(double x,double y,double z);
	~Vector();
	Vector Rotate(Vector Axis_dir,double alpha);
	Vector Normalize();
};

class Trajectory_Point {
public:
	Vector position;
	Vector direction;
	Vector rho;
	Trajectory_Point();
	~Trajectory_Point();
	double Critical_Energy(double gamma);
	double dAlpha(double dL);
};

class GenPhoton {
public:
	double dX,dY,divX,divY,dF,dP,energy;
};

class Histogram {
private:
	int number_of_bins;
	double delta,min,max;
	double *counts;
public:
	Histogram(double min,double max,int number_of_bins,bool logscale);
	~Histogram();
	double GetCount(int index);
	double GetFrequency(int index);
	double GetNormalized(int index);
	double GetX(int index);
	void Add(double x,double dY);
	bool logarithmic;
	double max_count,total_count;
	void ResetCounts();
};

class PARfileList {
public:
	int nbFiles;
	char** fileNames;
	PARfileList(int N);
	~PARfileList();
	PARfileList& operator=(const PARfileList &src);
	PARfileList(const PARfileList &src);
};

double DotProduct(Vector a,Vector b);
Vector Crossproduct(Vector v1,Vector v2);
Vector ScalarMult(Vector v,double r);
Vector Add(Vector v1,Vector v2);
double Min(double a,double b);
double Max(double a,double b);
int Min(int a,int b);
int Max(int a,int b);
double Sqr(double x);
double Gaussian(double sigma);

#endif
