#include <math.h>
#include "SynradTypes.h"
#include "GLApp/MathTools.h"
#include "File.h" //Error

double Trajectory_Point::Critical_Energy(const double &gamma) {
	double crit_en=2.959E-5*pow(gamma,3)/rho.Norme(); //rho in cm...
	return crit_en;
}

double Trajectory_Point::dAlpha(const double &dL) {
	return dL/rho.Norme();
}

Histogram::Histogram(double min_V,double max_V,int N,bool logscale){
	number_of_bins=N;
	counts = (ProfileSlice*)calloc(number_of_bins,sizeof(ProfileSlice));
	if (!counts) throw Error("Can't reserve memory for histogram");
	logarithmic=logscale;	
	this->min=min_V;
	this->max=max_V;
	
	if (!logarithmic) { //linearly distributed bins
		delta=(max-min)/number_of_bins;
	} else {
		delta=(log10(max)-log10(min))/number_of_bins;
	}
}

Histogram::~Histogram() {
	free(counts);
}

void Histogram::Add(const double &x,const ProfileSlice &increment) {
	if (x<max && x>=min) {
		int binIndex;
		if (!logarithmic) {
			binIndex = (int)((x - min) / delta);
		} else {
			binIndex = (int)((log10(x) - log10(min)) / delta);
		}
		double binX = GetX(binIndex);
		counts[binIndex] += increment;
	}
}

ProfileSlice Histogram::GetCounts(size_t index){
	return counts[index];
}



double Histogram::GetX(size_t index){
	double X;
	if (!logarithmic) {
		X=min+index*delta;
	} else {
		X=Pow10(log10(min)+index*delta);
	}
	return X;
}

void Histogram::ResetCounts(){
	memset(counts,0,number_of_bins*sizeof(ProfileSlice));
}

ProfileSlice & ProfileSlice::operator+=(const ProfileSlice & rhs)
{
	this->count_absorbed += rhs.count_absorbed;
	this->count_incident += rhs.count_incident;
	this->flux_absorbed += rhs.flux_absorbed;
	this->flux_incident += rhs.flux_incident;
	this->power_absorbed += rhs.power_absorbed;
	this->power_incident += rhs.power_incident;
	return *this;
}

TextureCell & TextureCell::operator+=(const TextureCell & rhs)
{
	this->count += rhs.count;
	this->flux += rhs.flux;
	this->power += rhs.power;
	return *this;
}
