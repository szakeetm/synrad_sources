#include "Region_mathonly.h"
#include "Random.h"
#include "Tools.h"
#include "GLApp\GLTypes.h"
#include <ctime>
#include <vector>
#include <string>


Vector Region::B(const Vector &position) {
	Vector result;
	//Creating references to components
	double* result_ptr[3]={&result.x,&result.y,&result.z};
	double* Bconst_ptr[3]={&B_const.x,&B_const.y,&B_const.z};
	Vector* Bdir_ptr[3]={&Bx_dir,&By_dir,&Bz_dir};
	int* Bmode_ptr[3]={&Bx_mode,&By_mode,&Bz_mode};
	double* Bperiod_ptr[3]={&Bx_period,&By_period,&Bz_period};
	double* Bphase_ptr[3]={&Bx_phase,&By_phase,&Bz_phase};
	Distribution2D* distr_ptr[3]={Bx_distr,By_distr,Bz_distr};
	double Ls_,ratio,K_,K_x,K_y;
	Vector dL;
	bool Bset=false; //means that i=0 case sets the remaining two components as well

	for (int i=0;i<3&&!Bset;i++) { //3 components
		dL=Add(position,ScalarMult(startPoint,-1.0));
		Ls_=DotProduct(*Bdir_ptr[i],dL);//distance towards Bx_dir direction (specified in .MAG file)
		Ls_-=(int)(Ls_/(*Bperiod_ptr)[i])*(*Bperiod_ptr)[i]; //substract filled periods
		switch(*Bmode_ptr[i]) {
		case B_MODE_CONSTANT:
			*result_ptr[i]=*Bconst_ptr[i];
			break;
		case B_MODE_COORDVALUES:
			*result_ptr[i]=distr_ptr[i]->InterpolateY(Ls_);
			break;
		case B_MODE_SINCOS:
			*result_ptr[i]=0.0;
			ratio=Ls_*2*PI/(*Bperiod_ptr[i]);
			for (int j=0;j<distr_ptr[i]->size;j++) {
				*result_ptr[i]+=distr_ptr[i]->valuesX[j]*pow(sin(((double)j+1.0)*ratio),j+1);
				*result_ptr[i]+=distr_ptr[i]->valuesY[j]*pow(cos(((double)j+1.0)*ratio),j+1);
			}
			break;
		case B_MODE_QUADRUPOLE: //mode 4
			Bset=true;
			result=quad.B(position);
			break;
		case B_MODE_ANALYTIC: //mode 5
			K_=2*PI/Bx_period;
			K_x=Bx_distr->valuesX[0];
			K_y=sqrt(Sqr(K_)-Sqr(K_x));
			result.x=K_x/K_y*Bx_distr->valuesY[0]*sinh(K_x*position.x)*sinh(K_y*position.y)*cos(K_*position.z);
			result.y=Bx_distr->valuesY[0]*cosh(K_x*position.x)*cosh(K_y*position.y)*cos(K_*position.z);
			result.z=-K_/K_y*Bx_distr->valuesY[0]*cosh(K_x*position.x)*sinh(K_y*position.y)*sin(K_*position.z);
			Bset=true; //all is set, don't run again for the remaining two components
			break;
		case B_MODE_HELICOIDAL:  //mode 6
			*result_ptr[i]=0.0;
			ratio=Ls_*2*PI/(*Bperiod_ptr[i]);
			for (int j=0;j<distr_ptr[i]->size;j++) {
				*result_ptr[i]+=distr_ptr[i]->valuesX[j]*sin((double)j*ratio)*cos(PI*(*Bphase_ptr[i])/(*Bperiod_ptr[i]));
				*result_ptr[i]+=distr_ptr[i]->valuesY[j]*cos((double)j*ratio)*sin(PI*(*Bphase_ptr[i])/(*Bperiod_ptr[i]));
			}
			break;
		case B_MODE_ROTATING_DIPOLE:
			*result_ptr[i]=0.0;
			ratio=Ls_*2*PI/(*Bperiod_ptr[i]);
			result.x=Bx_distr->valuesX[0]*sin(ratio);
			result.y=Bx_distr->valuesX[0]*cos(ratio);
			result.z=0.0;
			Bset=true;
			break;
		}
	}
	return result;
}

Region::Region(){
	
	generation_mode=SYNGEN_MODE_POWERWISE;
	selectedPoint=-1;
	isLoaded=false;
	//object placeholders until MAG files are loaded
	
	Bx_distr = new Distribution2D(1);
	By_distr = new Distribution2D(1);
	Bz_distr = new Distribution2D(1);
	beta_x_distr = new Distribution2D(1);
	beta_y_distr = new Distribution2D(1);
	eta_distr = new Distribution2D(1);
	etaprime_distr = new Distribution2D(1);
	e_spread_distr = new Distribution2D(1);
}

Region::~Region(){
	Distribution2D* distr_ptr[9]={Bx_distr,By_distr,Bz_distr,beta_x_distr,beta_y_distr,eta_distr,
		etaprime_distr,e_spread_distr};
	for (int i=0;i<9;i++)
		SAFE_DELETE(distr_ptr[i]);
	//if ((int)Points.size()>0) Points.clear();
	Points=std::vector<Trajectory_Point>();
}

Region::Region(const Region &src) {
	//alfa0=src.alfa0;
	betax=src.betax;
	betay=src.betay;
	
	Bx_dir=Vector(src.Bx_dir.x,src.Bx_dir.y,src.Bx_dir.z);
	Bx_distr=new Distribution2D((int)src.nbDistr_MAG.x);
	//*Bx_distr=*(src.Bx_distr);
	this->Bx_mode=src.Bx_mode;
	this->Bx_period=src.Bx_period;
	this->Bx_phase=src.Bx_phase;

	By_dir=Vector(src.By_dir.x,src.By_dir.y,src.By_dir.z);
	By_distr=new Distribution2D((int)src.nbDistr_MAG.y);
	//*By_distr=*(src.By_distr);
	this->By_mode=src.By_mode;
	this->By_period=src.By_period;
	this->By_phase=src.By_phase;

	Bz_dir=Vector(src.Bz_dir.x,src.Bz_dir.y,src.Bz_dir.z);
	Bz_distr=new Distribution2D((int)src.nbDistr_MAG.z);
	//*Bz_distr=*(src.Bz_distr);
	this->Bz_mode=src.Bz_mode;
	this->Bz_period=src.Bz_period;
	this->Bz_phase=src.Bz_phase;

	this->nbDistr_BXY=src.nbDistr_BXY;
	beta_x_distr=new Distribution2D(src.nbDistr_BXY);
	//*beta_x_distr=*(src.beta_x_distr);
	beta_y_distr=new Distribution2D(src.nbDistr_BXY);
	//*beta_y_distr=*(src.beta_y_distr);
	eta_distr=new Distribution2D(src.nbDistr_BXY);
	//*eta_distr=*(src.eta_distr);
	etaprime_distr=new Distribution2D(src.nbDistr_BXY);
	//*etaprime_distr=*(src.etaprime_distr);
	e_spread_distr=new Distribution2D(src.nbDistr_BXY);

	this->B_const=Vector(src.B_const.x,src.B_const.y,src.B_const.z);
	this->dL=src.dL;
	this->E=src.E;
	this->current=src.current;
	this->emittance=src.emittance;
	this->coupling=src.coupling;
	this->energy_spread=src.energy_spread;
	this->enable_ort_polarization=src.enable_ort_polarization;
	this->enable_par_polarization=src.enable_par_polarization;
	this->energy_hi=src.energy_hi;
	this->energy_low=src.energy_low;
	this->eta=src.eta;
	this->etaprime=src.etaprime;
	this->gamma=src.gamma;
	this->generation_mode=src.generation_mode;
	this->isLoaded=src.isLoaded;
	//this->i_struct1=src.i_struct1;
	this->limits=Vector(src.limits.x,src.limits.y,src.limits.z);
	//this->no_scans=src.no_scans;
	//this->nregions=src.nregions;
	this->particleMass=src.particleMass;
	//this->Points=src.Points; //does it work like that?
	this->psimaxX=src.psimaxX;
	this->psimaxY=src.psimaxY;
	this->quad=src.quad; //write assignment operator...
	this->selectedPoint=src.selectedPoint;
	this->startDir=Vector(src.startDir.x,src.startDir.y,src.startDir.z);
	this->startPoint=Vector(src.startPoint.x,src.startPoint.y,src.startPoint.z);
	//this->teta0=src.teta0;
	//this->x0=src.x0;
	//this->y0=src.y0;
	//this->z0=src.z0;
	this->nbPoints=(int)src.nbPoints;
	this->nbDistr_MAG=Vector(src.nbDistr_MAG.x,src.nbDistr_MAG.y,src.nbDistr_MAG.z);
	this->beta_kind=src.beta_kind;
	this->AABBmin=Vector(src.AABBmin.x,src.AABBmin.y,src.AABBmin.z);
	this->AABBmax=Vector(src.AABBmax.x,src.AABBmax.y,src.AABBmax.z);
	//this->fileName.assign(src.fileName);
}

Region& Region::operator=(const Region &src) {
	//alfa0=src.alfa0;
	betax=src.betax;
	betay=src.betay;
	
	Bx_dir=Vector(src.Bx_dir.x,src.Bx_dir.y,src.Bx_dir.z);
	Bx_distr=new Distribution2D((int)src.nbDistr_MAG.x);
	//*Bx_distr=*(src.Bx_distr);
	this->Bx_mode=src.Bx_mode;
	this->Bx_period=src.Bx_period;
	this->Bx_phase=src.Bx_phase;

	By_dir=Vector(src.By_dir.x,src.By_dir.y,src.By_dir.z);
	By_distr=new Distribution2D((int)src.nbDistr_MAG.y);
	//*By_distr=*(src.By_distr);
	this->By_mode=src.By_mode;
	this->By_period=src.By_period;
	this->By_phase=src.By_phase;

	Bz_dir=Vector(src.Bz_dir.x,src.Bz_dir.y,src.Bz_dir.z);
	Bz_distr=new Distribution2D((int)src.nbDistr_MAG.z);
	//*Bz_distr=*(src.Bz_distr);
	this->Bz_mode=src.Bz_mode;
	this->Bz_period=src.Bz_period;
	this->Bz_phase=src.Bz_phase;

	this->nbDistr_BXY=src.nbDistr_BXY;
	beta_x_distr=new Distribution2D(src.nbDistr_BXY);
	//*beta_x_distr=*(src.beta_x_distr);
	beta_y_distr=new Distribution2D(src.nbDistr_BXY);
	//*beta_y_distr=*(src.beta_y_distr);
	eta_distr=new Distribution2D(src.nbDistr_BXY);
	//*eta_distr=*(src.eta_distr);
	etaprime_distr=new Distribution2D(src.nbDistr_BXY);
	//*etaprime_distr=*(src.etaprime_distr);
	e_spread_distr=new Distribution2D(src.nbDistr_BXY);


	this->B_const=Vector(src.B_const.x,src.B_const.y,src.B_const.z);
	this->dL=src.dL;
	this->E=src.E;
	this->current=src.current;
	this->emittance=src.emittance;
	this->coupling=src.coupling;
	this->energy_spread=src.energy_spread;
	this->enable_ort_polarization=src.enable_ort_polarization;
	this->enable_par_polarization=src.enable_par_polarization;
	this->energy_hi=src.energy_hi;
	this->energy_low=src.energy_low;
	this->eta=src.eta;
	this->etaprime=src.etaprime;
	this->gamma=src.gamma;
	this->generation_mode=src.generation_mode;
	this->isLoaded=src.isLoaded;
	//this->i_struct1=src.i_struct1;
	this->limits=Vector(src.limits.x,src.limits.y,src.limits.z);
	//this->no_scans=src.no_scans;
	//this->nregions=src.nregions;
	this->particleMass=src.particleMass;
	this->Points=src.Points; //does it work like that?
	this->psimaxX=src.psimaxX;
	this->psimaxY=src.psimaxY;
	this->quad=src.quad; //write assignment operator...
	this->selectedPoint=src.selectedPoint;
	this->startDir=Vector(src.startDir.x,src.startDir.y,src.startDir.z);
	this->startPoint=Vector(src.startPoint.x,src.startPoint.y,src.startPoint.z);
	//this->teta0=src.teta0;
	//this->x0=src.x0;
	//this->y0=src.y0;
	//this->z0=src.z0;
	this->nbPoints=nbPoints;
	this->nbDistr_MAG=Vector(src.nbDistr_MAG.x,src.nbDistr_MAG.y,src.nbDistr_MAG.z);
	this->beta_kind=src.beta_kind;
	this->AABBmin=Vector(src.AABBmin.x,src.AABBmin.y,src.AABBmin.z);
	this->AABBmax=Vector(src.AABBmax.x,src.AABBmax.y,src.AABBmax.z);
	//this->fileName.assign(src.fileName);
	return *this;
}