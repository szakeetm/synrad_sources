#include "Region_mathonly.h"
#include "Random.h"
#include "GLApp\GLTypes.h"
#include "GLApp\MathTools.h"
#include <ctime>
#include <vector>
#include <string>

Vector3d Region_mathonly::B(size_t pointId, const Vector3d &offset) {
	//Calculates the magnetic field at a given point
	
	Vector3d result; //return value with the magnetic field vector

	//Creating references to X,Y,Z components -> we can avoid repeating code 3 times
	double* result_components[3]={&result.x,&result.y,&result.z};
	double* Bconst_components[3]={&params.B_const.x,&params.B_const.y,&params.B_const.z};
	Vector3d* Bdir_components[3]={&params.Bx_dir,&params.By_dir,&params.Bz_dir};
	int* Bmode_components[3]={&params.Bx_mode,&params.By_mode,&params.Bz_mode};
	double* Bperiod_components[3]={&params.Bx_period,&params.By_period,&params.Bz_period};
	double* Bphase_components[3]={&params.Bx_phase,&params.By_phase,&params.Bz_phase};
	Distribution2D* distr_components[3]={&Bx_distr,&By_distr,&Bz_distr};
	double Ls_,ratio,K_,K_x,K_y;
	bool Bset=false; //if true, then when we have calculated the X component, we already set Y and Z as well. Used for rotating dipole and analytic expressions

	Vector3d position_along_beam;
	/*if (Points.size()<=1) {
		position_along_beam=Points[0].position;
	} else {
		//Interpolate point
		Saturate(pointId,0,(double)Points.size()-1.00000001); //make sure we stay within limits
		Trajectory_Point previousPoint=Points[(int)pointId];
		Trajectory_Point nextPoint=Points[(int)pointId+1];
		double overshoot=pointId-(int)pointId;
		position_along_beam=Vector3d(
			WEIGH(previousPoint.position.x,nextPoint.position.x,overshoot),
			WEIGH(previousPoint.position.y,nextPoint.position.y,overshoot),
			WEIGH(previousPoint.position.z,nextPoint.position.z,overshoot)
			);
	}*/
	position_along_beam = Points[pointId].position;
	Vector3d position_with_offset = position_along_beam + offset;
	
	for (size_t componentIndex=0;componentIndex<3&&!Bset;componentIndex++) { //X,Y,Z components
		
		switch(*Bmode_components[componentIndex]) {
		case B_MODE_CONSTANT:
			*result_components[componentIndex]=*Bconst_components[componentIndex];
			break;
		case B_MODE_DIRECTION:
			Ls_ = Dot(Bdir_components[componentIndex]->Normalized(), position_with_offset - params.startPoint);//distance towards Bx_dir direction (specified in .MAG file)
			Ls_-=(int)(Ls_/(*Bperiod_components)[componentIndex])*(*Bperiod_components)[componentIndex]; //substract filled periods
			*result_components[componentIndex]=distr_components[componentIndex]->InterpolateY(Ls_,false);
			break;
		case B_MODE_ALONGBEAM:
			Ls_=pointId*this->params.dL_cm; //distance along the beam path
			Ls_-=(int)(Ls_/(*Bperiod_components)[componentIndex])*(*Bperiod_components)[componentIndex]; //substract filled periods
			*result_components[componentIndex]=distr_components[componentIndex]->InterpolateY(Ls_,false);
			break;
		case B_MODE_SINCOS:
			Ls_ = Dot(Bdir_components[componentIndex]->Normalized(), position_with_offset - params.startPoint);//distance towards Bx_dir direction (specified in .MAG file)
			//Ls_ -= (int)(Ls_ / (*Bperiod_components)[componentIndex])*(*Bperiod_components)[componentIndex]; //substract filled periods

			*result_components[componentIndex]=0.0;
			ratio=Ls_*2*PI/(*Bperiod_components[componentIndex]);
			for (int j=0;j<distr_components[componentIndex]->GetSize();j++) { //loop through orders
				*result_components[componentIndex]+=distr_components[componentIndex]->GetX(j)*pow(sin((double)(j+1)*ratio),j+1);
				*result_components[componentIndex]+=distr_components[componentIndex]->GetY(j)*pow(cos((double)(j+1)*ratio),j+1);
			}
			break;
		case B_MODE_QUADRUPOLE: //mode 4
			Bset=true;
			result= params.quad_params.B(position_with_offset);
			break;
		case B_MODE_ANALYTIC: //mode 5
			K_=2*PI/ params.Bx_period;
			K_x=Bx_distr.GetX(0);
			K_y=sqrt(Sqr(K_)-Sqr(K_x));
			result.x=K_x/K_y*Bx_distr.GetY(0)*sinh(K_x*position_with_offset.x)*sinh(K_y*position_with_offset.y)*cos(K_*position_with_offset.z);
			result.y=Bx_distr.GetY(0)*cosh(K_x*position_with_offset.x)*cosh(K_y*position_with_offset.y)*cos(K_*position_with_offset.z);
			result.z=-K_/K_y*Bx_distr.GetY(0)*cosh(K_x*position_with_offset.x)*sinh(K_y*position_with_offset.y)*sin(K_*position_with_offset.z);
			Bset=true; //all is set, don't run again for the remaining two components
			break;
		case B_MODE_HELICOIDAL:  //helicoidal field, as per PE 31/1/2001
			Ls_ = Dot(Bdir_components[componentIndex]->Normalized(), position_with_offset - params.startPoint);//distance towards Bx_dir direction (specified in .MAG file)
			//Ls_ -= (int)(Ls_ / (*Bperiod_components)[componentIndex])*(*Bperiod_components)[componentIndex]; //substract filled periods

			*result_components[componentIndex]=0.0;
			ratio=Ls_*2*PI/(*Bperiod_components[componentIndex]);
			for (int j=0;j<distr_components[componentIndex]->GetSize();j++) { //loop through orders (sinX+sin^2X+...)
				*result_components[componentIndex]+=distr_components[componentIndex]->GetX(j)*sin((double)(j+1)*ratio)*cos(PI*(*Bphase_components[componentIndex])/(*Bperiod_components[componentIndex]));
				*result_components[componentIndex]+=distr_components[componentIndex]->GetY(j)*cos((double)(j+1)*ratio)*sin(PI*(*Bphase_components[componentIndex])/(*Bperiod_components[componentIndex]));
			}
			break;
		case B_MODE_ROTATING_DIPOLE: //rotating dipole field ( see S. Duncan's presentation on generation of 20MeV circ.pol. photons 
			Ls_ = Dot(Bdir_components[componentIndex]->Normalized(), position_with_offset - params.startPoint);//distance towards Bx_dir direction (specified in .MAG file)
			//Ls_ -= (int)(Ls_ / (*Bperiod_components)[componentIndex])*(*Bperiod_components)[componentIndex]; //substract filled periods

			*result_components[componentIndex]=0.0;
			ratio=Ls_*2*PI/(*Bperiod_components[componentIndex]);
			result.x=Bx_distr.GetX(0)*sin(ratio);
			result.y=Bx_distr.GetX(0)*cos(ratio);
			result.z=0.0;
			Bset=true;
			break;
		}
	}
	return result;
}

Region_mathonly::Region_mathonly(){

	//default values
	//generation_mode = SYNGEN_MODE_POWERWISE;
	params.nbDistr_BXY = 0;
	params.dL_cm=0.01;
	params.limits=Vector3d(1000,1000,0.1);
	params.startPoint=Vector3d(0,0,0);
	params.startDir=Vector3d(0,0,1);
	params.particleMass_GeV=-0.0005110034; //electron
	params.E_GeV=1;params.current_mA=1;
	params.betax_const_cm= params.betay_const_cm= params.eta_x_const_cm= params.eta_x_prime_const= params.energy_spread_percent= params.emittance_cm=0.0;
	params.coupling_percent=100.0;
	params.energy_low_eV=10;
	params.energy_hi_eV=1e6;
	params.enable_ort_polarization= params.enable_par_polarization=true;
	params.polarizationCompIndex = 0; //full polarization
	params.psimaxX_rad= params.psimaxY_rad=PI;
	params.Bx_mode= params.By_mode= params.Bz_mode=B_MODE_CONSTANT;
	params.B_const=Vector3d(0,0,0);
	params.showPhotons = true;
}

/*
Region_mathonly::Region_mathonly(const Region_mathonly &src) {
	this->params = src.params;
	Bx_distr = Distribution2D(src.Bx_distr);
	By_distr = Distribution2D(src.By_distr);
	Bz_distr = Distribution2D(src.Bz_distr);
	this->betaFunctions = src.betaFunctions;
	this->Points = src.Points;
}

Region_mathonly& Region_mathonly::operator=(const Region_mathonly &src) {
	this->params = src.params;
	Bx_distr = Distribution2D(src.Bx_distr);
	By_distr = Distribution2D(src.By_distr);
	Bz_distr = Distribution2D(src.Bz_distr);
	betaFunctions = src.betaFunctions;
	this->Points = src.Points;
	return *this;
}
*/