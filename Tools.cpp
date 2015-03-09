#include "Tools.h"
#include <math.h>
#include <malloc.h>
#include <string>
#include "Random.h"

#define SAFE_FREE(x) if(x) { free(x);x=NULL; }
#define SAFE_DELETE(x) if(x) { delete x;x=NULL; }

double DotProduct(const Vector &a,const Vector &b) {
	double result=0.0;
	result+=a.x*b.x;
	result+=a.y*b.y;
	result+=a.z*b.z;
	return result;
}

Vector Crossproduct(const Vector &v1,const Vector &v2) {
	Vector result;
	result.x = (v1.y)*(v2.z) - (v1.z)*(v2.y);
	result.y = (v1.z)*(v2.x) - (v1.x)*(v2.z);
	result.z = (v1.x)*(v2.y) - (v1.y)*(v2.x);
	return result;
}

Vector ScalarMult(const Vector &v,const double &r) {
  Vector result;
  result.x =v.x* r;
  result.y =v.y* r;
  result.z =v.z* r;
  return result;
}

Vector RandomPerpendicularVector(const Vector &v,const double &length){
	Vector randomVector=Vector(rnd(),rnd(),rnd());
	Vector perpendicularVector=Crossproduct(randomVector,v);
	perpendicularVector=perpendicularVector.Normalize();
	return ScalarMult(perpendicularVector,length);
}

Vector Add(const Vector &v1,const Vector &v2) {
  Vector result;
  result.x = (v1.x) + (v2.x);
  result.y = (v1.y) + (v2.y);
  result.z = (v1.z) + (v2.z);
  return result;
}

Vector Subtract(const Vector &v1,const Vector &v2) {
  Vector result;
  result.x = (v1.x) - (v2.x);
  result.y = (v1.y) - (v2.y);
  result.z = (v1.z) - (v2.z);
  return result;
}

Vector Weigh(const Vector &v1,const Vector &v2,const double &weight){ //interpolate between vectors A and B. weight is [0..1]
	return Vector(
		WEIGH(v1.x,v2.x,weight),
		WEIGH(v1.y,v2.y,weight),
		WEIGH(v1.z,v2.z,weight)
		);
}

double Min(const double &a,const double &b){
	return (a>b)?b:a;
}

double Max(const double &a,const double &b){
	return (a>b)?a:b;
}

int Min(const int &a,const int &b){
	return (a>b)?b:a;
}

int Max(const int &a,const int &b){
	return (a>b)?a:b;
}

double Sqr(const double &x) {
	return x*x;
}

double Vector::Norme() {
	return sqrt(DotProduct(*this,*this));
}

Vector::Vector(const double &x,const double &y,const double &z) {
	this->x=x;
	this->y=y;
	this->z=z;
}

Vector::Vector() {
}

/*Vector::~Vector() {
}*/

Vector::Vector(const Vector &src){
	this->x=src.x;
	this->y=src.y;
	this->z=src.z;
}

Vector& Vector::operator=(const Vector &src){
	this->x=src.x;
	this->y=src.y;
	this->z=src.z;
	return *this;
}

double Trajectory_Point::Critical_Energy(const double &gamma) {
	double crit_en=2.959E-5*pow(gamma,3)/rho.Norme();
	//if (!(crit_en==crit_en)) __debugbreak();
	return crit_en;
}

double Trajectory_Point::dAlpha(const double &dL) {
	return dL/rho.Norme();
}

Vector Vector::Normalize() {
	if (Norme()==0.0) return *this;
	return ScalarMult(*this,1/Norme());
}

Vector Vector::Rotate(Vector Axis,double theta) {
	Axis=Axis.Normalize();
	Vector result;
	double u=Axis.x;
	double v=Axis.y;
	double w=Axis.z;
	double cos_theta=cos(theta);
	double sin_theta=sin(theta);
	double product=(u*x+v*y+w*z)*(1-cos_theta);
	result.x=u*product+x*cos_theta+(-w*y+v*z)*sin_theta;
	result.y=v*product+y*cos_theta+(w*x-u*z)*sin_theta;
	result.z=w*product+z*cos_theta+(-v*x+u*y)*sin_theta;
	return result;
}

Histogram::Histogram(double min_V,double max_V,int N,bool logscale){
	number_of_bins=N;
	counts = (double*)malloc(number_of_bins*sizeof(double));
	if (!counts) throw Error("Can't reserve memory for histogram");
	memset(counts,0,number_of_bins*sizeof(double));
	logarithmic=logscale;	
	this->min=min_V;
	this->max=max_V;
	
	if (!logarithmic) { //linearly distributed bins
		delta=(max-min)/number_of_bins;
	} else {
		delta=(log(max)-log(min))/number_of_bins;
	}
	max_count=0.0;
	total_count=0.0;
}

Histogram::~Histogram() {
	free(counts);
}

void Histogram::Add(const double &x,const double &dY,const double &bandwidth) {
	if (x<max && x>=min) {
		int binIndex;
		if (!logarithmic) {
			binIndex = (int)(.5 + (x - min) / delta);
		} else {
			binIndex = (int)(.5 + (log(x) - log(min)) / delta);
		}
		double binX = GetX(binIndex);
		/*if (binIndex<number_of_bins && bandwidth == -1 || (abs(x - binX) / binX) < (bandwidth / 2)) {
			counts[binIndex] += dY;
			total_count += dY;
		}*/
		double factor = 1.0;
		if (bandwidth > -1) factor = bandwidth / delta;
		if (binIndex < number_of_bins) {
			counts[binIndex] += dY*factor;
			total_count += dY*factor;
		}
		if (counts[binIndex]>max_count) max_count = counts[binIndex];
	}
}

double Histogram::GetCount(int index){
	return counts[index];
}

double Histogram::GetFrequency(int index){
	return (counts[index]/total_count);
}

double Histogram::GetX(int index){
	double X;
	if (!logarithmic) {
		X=min+index*delta;
	} else {
		X=exp(log(min)+index*delta);
	}
	return X;
}

double Histogram::GetNormalized(int index){
	return (counts[index]/max_count);
}

void Histogram::ResetCounts(){
	max_count=total_count=0.0;
	memset(counts,0,number_of_bins*sizeof(double));
}

PARfileList::PARfileList(int N){
	fileNames=new char*[N];
	for (int i=0;i<N;i++) {
		fileNames[i]=new char[512];
		*(fileNames[i])=NULL;
	}
	nbFiles=N;
}

PARfileList::~PARfileList(){
	for (int i=0;i<nbFiles;i++)
		SAFE_DELETE(fileNames[i]);
	SAFE_DELETE(fileNames);
}

PARfileList& PARfileList::operator=(const PARfileList &src){
	this->fileNames=new char*[src.nbFiles];
	this->nbFiles=src.nbFiles;
	for (int i=0;i<nbFiles;i++)
		this->fileNames[i]=_strdup(src.fileNames[i]);
	return *this;
};

PARfileList::PARfileList(const PARfileList &src) {
	this->fileNames=new char*[src.nbFiles];
	this->nbFiles=src.nbFiles;
	for (int i=0;i<nbFiles;i++)
		this->fileNames[i]=_strdup(src.fileNames[i]);
};