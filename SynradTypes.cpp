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
	counts = (double*)malloc(number_of_bins*sizeof(double));
	if (!counts) throw Error("Can't reserve memory for histogram");
	memset(counts,0,number_of_bins*sizeof(double));
	logarithmic=logscale;	
	this->min=min_V;
	this->max=max_V;
	
	if (!logarithmic) { //linearly distributed bins
		delta=(max-min)/number_of_bins;
	} else {
		delta=(log10(max)-log10(min))/number_of_bins;
	}
	//max_count=0.0;
	total_count=0.0;
}

Histogram::~Histogram() {
	free(counts);
}

void Histogram::Add(const double &x,const double &dY) {
	if (x<max && x>=min) {
		int binIndex;
		if (!logarithmic) {
			binIndex = (int)((x - min) / delta);
		} else {
			binIndex = (int)((log10(x) - log10(min)) / delta);
		}
		double binX = GetX(binIndex);
		/*if (binIndex<number_of_bins && bandwidth == -1 || (abs(x - binX) / binX) < (bandwidth / 2)) {
			counts[binIndex] += dY;
			total_count += dY;
		}*/
		
		//if (binIndex < number_of_bins) { //x<max
			counts[binIndex] += dY;
			total_count += dY;
		//}
		//if (counts[binIndex]>max_count) max_count = counts[binIndex];
	}
}

double Histogram::GetCount(size_t index){
	return counts[index];
}

double Histogram::GetFrequency(size_t index){
	return (counts[index]/total_count);
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
	/*max_count=*/total_count=0.0;
	memset(counts,0,number_of_bins*sizeof(double));
}