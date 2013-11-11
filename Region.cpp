#include "Region.h"
#include "Random.h"
#include "Tools.h"
#include <ctime>
#include <vector>
#include "GLApp\GLToolkit.h"
#include "GLApp\GLProgress.h"
#include "SynRad.h"
#include "GeometryViewer.h"
#include <string>
#include <direct.h> //for CWD
using namespace std;
//extern int antiAliasing;

vector<Trajectory_Point> Region::CalculateTrajectory(double dL,double E,Vector limits,Vector pos0,Vector dir0,int max_steps){
	//global variables. All distances in cm!
	Trajectory_Point current_point;
	extern Distribution2D polarization_distribution,integral_N_photons,integral_SR_power,g1h2_distribution;

	//initial position and speed
	current_point.position=pos0;
	current_point.direction=dir0;
	current_point.rho=Vector(0.0,0.0,1e30); //big enough so no rotation
	AABBmin=AABBmax=pos0;

	//Calculate trajectory
	std::vector <Trajectory_Point>Points;
	Points.push_back(current_point); //beam starting position

	for (int i=0;i<max_steps&&!isOutsideBoundaries(current_point.position,i==0);i++) {
		current_point=OneStep(current_point,dL,E); //step forward
		Points.push_back(current_point); //store point
		//extend trajectory limits if needed
		if (current_point.position.x<AABBmin.x) AABBmin.x=current_point.position.x;
		if (current_point.position.y<AABBmin.y) AABBmin.y=current_point.position.y;
		if (current_point.position.z<AABBmin.z) AABBmin.z=current_point.position.z;
		if (current_point.position.x>AABBmax.x) AABBmax.x=current_point.position.x;
		if (current_point.position.y>AABBmax.y) AABBmax.y=current_point.position.y;
		if (current_point.position.z>AABBmax.z) AABBmax.z=current_point.position.z;
	}
	return Points;
}

//constant B field
Vector Region::B(Vector position) {
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

Trajectory_Point Region::OneStep(Trajectory_Point p0,double dL,double E) {
	Vector rotation_axis=Crossproduct(p0.direction,p0.rho);
	//Calculate new point
	Trajectory_Point p;
	p.direction=p0.direction.Rotate(rotation_axis,dL/p0.rho.Norme());
	p.position=Add(p0.position,ScalarMult(p.direction,dL/p.direction.Norme()));

	Vector B_local=B(p.position);
	if (B_local.Norme()<VERY_SMALL) {//straight line
		p.rho=Vector(0.0,0.0,1.0e30); //no rotation
		return p;
	}
	double B_parallel=DotProduct(p.direction,B_local);
	double B_orthogonal=sqrt(Sqr(B_local.Norme())-Sqr(B_parallel));
	Vector lorentz_force_direction=Crossproduct(p.direction,B_local).Normalize();
	if (particleMass<0.0) lorentz_force_direction=ScalarMult(lorentz_force_direction,-1.0);
	double radius;
	if (B_orthogonal>VERY_SMALL)
		radius=E/0.00299792458/B_orthogonal;
	else radius=1.0E30;
	p.rho=ScalarMult(lorentz_force_direction,radius);
	return p;
}

bool Region::isOutsideBoundaries(Vector a,BOOL recalcDirs){
	static double xDir=(limits.x>startPoint.x)?1.0:-1.0;
	static double yDir=(limits.y>startPoint.y)?1.0:-1.0;
	static double zDir=(limits.z>startPoint.z)?1.0:-1.0;
	if (recalcDirs) {
		xDir=(limits.x>startPoint.x)?1.0:-1.0;
		yDir=(limits.y>startPoint.y)?1.0:-1.0;
		zDir=(limits.z>startPoint.z)?1.0:-1.0;
	}
	return ((a.x-limits.x)*xDir>0.0)||((a.y-limits.y)*yDir>0.0)||((a.z-limits.z)*zDir>0.0);
}

void Region::LoadPAR(FileReader *file,GLProgress *prg){

	//Creating references to X,Y,Z components (to avoid writing code 3 times)
	//double* result_ptr[3]={&result.x,&result.y,&result.z};
	double* Bconst_ptr[3]={&B_const.x,&B_const.y,&B_const.z};
	Vector* Bdir_ptr[3]={&Bx_dir,&By_dir,&Bz_dir};
	int* Bmode_ptr[3]={&Bx_mode,&By_mode,&Bz_mode};
	double* Bperiod_ptr[3]={&Bx_period,&By_period,&Bz_period};
	double* Bphase_ptr[3]={&Bx_phase,&By_phase,&Bz_phase};
	Distribution2D* distr_ptr[3]={Bx_distr,By_distr,Bz_distr};
	std::string* MagFileNamePtr[3]={&MAGXfileName,&MAGYfileName,&MAGZfileName};

	prg->SetMessage("Reading parameter file...");
	int no_scans=file->ReadInt(); //unused
	int nregions=file->ReadInt(); //unused
	generation_mode=file->ReadInt(); //default value
	file->JumpComment();

	particleMass=file->ReadDouble();
	E=file->ReadDouble();
	energy_low=file->ReadDouble();
	energy_hi=file->ReadDouble();
	file->JumpComment();

	gamma=abs(E/particleMass);

	dL=file->ReadDouble();
	double x0=file->ReadDouble();
	double y0=file->ReadDouble();
	double z0=file->ReadDouble();
	startPoint=Vector(x0,y0,z0);

	double teta0=file->ReadDouble();
	double alfa0=file->ReadDouble();
	startDir=Vector(-cos(alfa0)*sin(teta0),-sin(alfa0),cos(alfa0)*cos(teta0));
	file->JumpComment();

	double xmax=file->ReadDouble();
	double ymax=file->ReadDouble();
	double zmax=file->ReadDouble();
	limits=Vector(xmax,ymax,zmax);
	file->JumpComment();

	int i_struct1=file->ReadInt(); //thrown away
	file->JumpComment();

	enable_par_polarization=(file->ReadInt()==1);
	enable_ort_polarization=(file->ReadInt()==1);
	file->JumpComment();

	for (int i=0;i<3;i++) {
		*Bmode_ptr[i]=file->ReadInt();
		if (*Bmode_ptr[i]==B_MODE_CONSTANT) {
			*Bconst_ptr[i]=file->ReadDouble();
			file->JumpComment();
		} else { //parameter file to read
			int fileInt=file->ReadInt();file->JumpComment();
			char fileName[256];
			strcpy(fileName,file->GetName());
			char *filebegin= strrchr(fileName,'\\');
			if (filebegin) filebegin++;
			else filebegin=fileName;
			sprintf(filebegin,"%d.mag",fileInt); //append .MAG extension

			
				if (!FileUtils::Exist(fileName)) { //no .MAG file
					char tmperr[256];
					sprintf(tmperr,"Referenced MAG file doesn't exist:\n%s\nReferred to: ",fileName);
					Error err=file->MakeError(tmperr);
					throw err;
				}
				FileReader MAGfile(fileName);
				MagFileNamePtr[i]->assign(fileName);
				*(distr_ptr[i])=LoadMAGFile(&MAGfile,Bdir_ptr[i],Bperiod_ptr[i],Bphase_ptr[i],*Bmode_ptr[i]);
			
		}
	}

	emittance=eta=etaprime=coupling=energy_spread=betax=betay=0.0; //default values
	nbDistr_BXY=0;

	emittance=file->ReadDouble();
	betax=file->ReadDouble();
	if (!(betax<0.0)) {
		file->wasLineEnd=false;
		betay=file->ReadDouble();
		double value=file->ReadDouble();
		if (!file->wasLineEnd) {
			eta=value;
			etaprime=file->ReadDouble();
			value=file->ReadDouble();
			if (!file->wasLineEnd) {
				coupling=value;
				energy_spread=file->ReadDouble();
				value=file->ReadDouble();
			}
		}
		psimaxX=value*0.001;
		file->wasLineEnd=false;
	} else psimaxX=file->ReadDouble()*0.001;
	psimaxY=file->ReadDouble()*0.001;
	file->JumpComment();

	beta_kind=0;
	if (betax<0.0) { //read BXY file
		double fileDouble=-1.0*betax;
		char fileName[256];
		strcpy(fileName,file->GetName());
		char *filebegin= strrchr(fileName,'\\');
		if (filebegin) filebegin++;
		else filebegin=fileName;
		sprintf(filebegin,"%d.bxy",(int)(fileDouble+0.5)); //append .BXY extension
		
			if (!FileUtils::Exist(fileName)) { //no .BXY file
				char tmperr[256];
				sprintf(tmperr,"Referenced BXY file doesn't exist:\n%s\nReferred to: ",fileName);
				Error err=file->MakeError(tmperr);
				throw err;
			}
			FileReader BXYfile(fileName);
			BXYfileName.assign(fileName);
			nbDistr_BXY=LoadBXY(&BXYfile,beta_x_distr,beta_y_distr,eta_distr,etaprime_distr,coupling_distr,e_spread_distr);
		
	}

	prg->SetMessage("Calculating trajectory...");
	Points=CalculateTrajectory(dL,E,limits,startPoint,startDir,1000000); //max 1 million points
	isLoaded=true;
	extern SynRad *theApp;
	SynRad *mApp = (SynRad *)theApp;
	if (mApp->regionInfo) mApp->regionInfo->Update();
}

Distribution2D Region::LoadMAGFile(FileReader *file,Vector *dir,double *period,double *phase,int mode){
	Distribution2D result(1);
	if (mode==B_MODE_QUADRUPOLE) {
		double quadrX=file->ReadDouble();
		double quadrY=file->ReadDouble();
		double quadrZ=file->ReadDouble();
		file->JumpComment();
		quad.center=Vector(quadrX,quadrY,quadrZ);

		quad.alfa_q=file->ReadDouble();
		quad.beta_q=file->ReadDouble();
		quad.rot_q=file->ReadDouble();
		file->JumpComment();

		quad.K_q=file->ReadDouble()/100.0; //m->cm conversion
		quad.L_q=file->ReadDouble();

		quad.cosalfa_q=cos(quad.alfa_q);
		quad.sinalfa_q=sin(quad.alfa_q);
		quad.cosbeta_q=cos(quad.beta_q);
		quad.sinbeta_q=sin(quad.beta_q);
		quad.cosrot_q=cos(quad.rot_q);
		quad.sinrot_q=sin(quad.rot_q);

		quad.direction=Vector(-quad.cosalfa_q*quad.sinbeta_q,
			-quad.sinalfa_q,
			quad.cosalfa_q*quad.cosbeta_q);

	} else {
		*period=file->ReadDouble();
		if (mode==B_MODE_HELICOIDAL) *phase=file->ReadDouble();
		file->JumpComment();

		double dirx=file->ReadDouble();
		double diry=file->ReadDouble();
		double dirz=file->ReadDouble();
		*dir=Vector(dirx,diry,dirz);
		file->JumpComment();

		int N=file->ReadInt();
		result=Distribution2D(N);
		for (int i=0;i<N;i++) {
			result.valuesX[i]=file->ReadDouble();
			result.valuesY[i]=file->ReadDouble();
			file->JumpComment();
		}
	}
	return result;
}

Region::Region(){

	generation_mode=SYNGEN_MODE_FLUXWISE;
	emittance=eta=etaprime=coupling=energy_spread=betax=betay=0.0;
	selectedPoint=-1;
	isLoaded=FALSE;
	//object placeholders until MAG files are loaded

	Bx_distr = new Distribution2D(1);
	By_distr = new Distribution2D(1);
	Bz_distr = new Distribution2D(1);
	beta_x_distr = new Distribution2D(1);
	beta_y_distr = new Distribution2D(1);
	eta_distr = new Distribution2D(1);
	etaprime_distr = new Distribution2D(1);
	coupling_distr = new Distribution2D(1);
	e_spread_distr = new Distribution2D(1);

	MAGXfileName="";
	MAGYfileName="";
	MAGZfileName="";
	BXYfileName="";

	//default values for empty region
	dL=0.01;
	limits=Vector(1000,1000,0.1);
	startPoint=Vector(0,0,0);
	startDir=Vector(0,0,1);
	particleMass=-0.0005110034; //electron
	E=1;
	betax=betay=eta=etaprime=coupling=energy_spread=emittance=0,0;
	energy_low=10;
	energy_hi=1e6;
	enable_ort_polarization=enable_par_polarization=TRUE;
	psimaxX=psimaxY=PI;
	Bx_mode=By_mode=Bz_mode=B_MODE_CONSTANT;
	B_const=Vector(0,0,0);
	this->nbPoints=0;

}

Region::~Region(){
	Distribution2D* distr_ptr[9]={Bx_distr,By_distr,Bz_distr,beta_x_distr,beta_y_distr,eta_distr,
		etaprime_distr,coupling_distr,e_spread_distr};
	for (int i=0;i<9;i++)
		SAFE_DELETE(distr_ptr[i]);
	Points=std::vector<Trajectory_Point>();
}

void Region::Render(int dispNumTraj,GLMATERIAL *B_material,double vectorLength){
	if (!whiteBg) glPointSize(1.0f);
	else glPointSize(2.0f);
	//if (antiAliasing) glEnable(GL_POINT_SMOOTH);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);
	glDisable(GL_CULL_FACE);
	if (!whiteBg) glColor3f(1.0f,0.9f,0.2f);
	else glColor3f(1.0f,0.5f,0.2f);
	glBegin(GL_POINTS);
	//Region
	int N=(int)Points.size();
	int increment=Max(1,(int)(N/dispNumTraj));
	for (int i=0;i<N;i+=increment) {
		// Draw dot

		//if (hiddenVertex) glDisable(GL_DEPTH_TEST);

		Vector *v = &(Points[i].position);
		glVertex3d(v->x,v->y,v->z);
	}
	glEnd();

	//Selected trajectory point
	if (selectedPoint!=-1) {
		// Draw dot
		if (!whiteBg) glPointSize(6.0f);
		else glPointSize(7.0f);
		glEnable(GL_POINT_SMOOTH);
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
		glDisable(GL_BLEND);
		glDisable(GL_CULL_FACE);

		glColor3f(1.0f,0.0f,0.0f);

		glBegin(GL_POINTS);
		glVertex3d(Points[selectedPoint].position.x,Points[selectedPoint].position.y,Points[selectedPoint].position.z);
		glEnd();

		GLToolkit::SetMaterial(B_material);
		if (antiAliasing) {
			glEnable(GL_LINE_SMOOTH);
			//glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			//glColor4f(0.0f,0.0f,1.0f,0.5f);
			glEnable(GL_DEPTH_TEST);
		}
		glLineWidth(0.5f);
		Vector O=Points[selectedPoint].position;
		Vector B_local=B(O);
		Vector B_local_norm=B_local.Normalize();
		Vector Rho_local=Points[selectedPoint].rho;
		Vector Rho_local_norm=Rho_local.Normalize();
		Vector dir=Points[selectedPoint].direction;
		double factor=vectorLength;
		if (B_local.Norme()>0.0) 

			GLToolkit::DrawVector(O.x,O.y,O.z,O.x+factor*B_local_norm.x,O.y+factor*B_local_norm.y,O.z+factor*B_local_norm.z);
		if (Rho_local.Norme()>0.0) 
			GLToolkit::DrawVector(O.x,O.y,O.z,O.x+factor*Rho_local_norm.x,O.y+factor*Rho_local_norm.y,O.z+factor*Rho_local_norm.z);
		if (dir.Norme()>0.0) 
			GLToolkit::DrawVector(O.x,O.y,O.z,O.x+factor*dir.x,O.y+factor*dir.y,O.z+factor*dir.z);
		if (antiAliasing) glDisable(GL_LINE_SMOOTH);
		glPointSize(3.0f);
		glColor3f(0.5f,1.0f,1.0f);
		glBegin(GL_POINTS);
		glVertex3d(O.x,O.y,O.z);
		glEnd();
		//glEnable(GL_BLEND);
		char point_label[256];
		char B_label[128];
		char Rho_label[128];
		char dir_label[128];
		sprintf(point_label,"point #%d (%g,%g,%g) L=%g",selectedPoint+1,O.x,O.y,O.z,selectedPoint*dL);
		sprintf(B_label,"B (%g Tesla)",B_local.Norme());
		sprintf(Rho_label,"Rho (%g cm)",Rho_local.Norme());
		sprintf(dir_label,"Direction");
		GLToolkit::GetDialogFont()->SetTextColor(0.5f,0.6f,1.0f);
		GLToolkit::DrawStringInit();
		GLToolkit::DrawString((float)(O.x),(float)(O.y),(float)(O.z),point_label,GLToolkit::GetDialogFont());

		GLToolkit::DrawString((float)(O.x+factor*B_local_norm.x),(float)(O.y+factor*B_local_norm.y),(float)(O.z+factor*B_local_norm.z),B_label,GLToolkit::GetDialogFont());
		GLToolkit::DrawString((float)(O.x+factor*Rho_local_norm.x),(float)(O.y+factor*Rho_local_norm.y),(float)(O.z+factor*Rho_local_norm.z),Rho_label,GLToolkit::GetDialogFont());
		GLToolkit::DrawString((float)(O.x+factor*dir.x),(float)(O.y+factor*dir.y),(float)(O.z+factor*dir.z),dir_label,GLToolkit::GetDialogFont());
		GLToolkit::DrawStringRestore();
	}
}

void Region::SelectTrajPoint(int x,int y) {

	int i;
	if(!isLoaded) return;

	// Select a vertex on a mouse click in 3D perspectivce view 
	// (x,y) are in screen coordinates
	// TODO: Handle clipped polygon

	// Check intersection of the facet and a "perspective ray"
	int *allXe = (int *)malloc(Points.size() * sizeof(int));
	int *allYe = (int *)malloc(Points.size() * sizeof(int));

	// Transform points to screen coordinates
	BOOL *ok = (BOOL *)malloc(Points.size()*sizeof(BOOL));
	for(i=0;i<(int)Points.size();i++)
		ok[i] = GLToolkit::Get2DScreenCoord((float)Points[i].position.x,(float)Points[i].position.y,(float)Points[i].position.z,allXe+i,allYe+i);

	//Get Closest Point to click
	double minDist=9999;
	double distance;
	int minId=-1;
	for(i=0;i<(int)Points.size();i++) {
		if (ok[i] && !(allXe[i]<0) && !(allYe[i]<0)) { //calculate only for points on screen
			distance=pow((double)(allXe[i]-x),2)+pow((double)(allYe[i]-y),2);
			if (distance<minDist) {
				minDist=distance;
				minId=i;
			}
		}
	}


	selectedPoint=-1;

	if (minDist<250.0) {
		selectedPoint=minId;
	}

	free(allXe);
	free(allYe);
	free(ok);
	//UpdateSelection();
	extern SynRad *theApp;
	SynRad *mApp = (SynRad *)theApp;
	if (mApp->regionInfo) mApp->regionInfo->Update();
}

Region::Region(const Region &src) {
	//alfa0=src.alfa0;
	betax=src.betax;
	betay=src.betay;

	Bx_dir=Vector(src.Bx_dir.x,src.Bx_dir.y,src.Bx_dir.z);
	Bx_distr=new Distribution2D((int)src.nbDistr_MAG.x);
	*Bx_distr=*(src.Bx_distr);
	this->Bx_mode=src.Bx_mode;
	this->Bx_period=src.Bx_period;
	this->Bx_phase=src.Bx_phase;

	By_dir=Vector(src.By_dir.x,src.By_dir.y,src.By_dir.z);
	By_distr=new Distribution2D((int)src.nbDistr_MAG.y);
	*By_distr=*(src.By_distr);
	this->By_mode=src.By_mode;
	this->By_period=src.By_period;
	this->By_phase=src.By_phase;

	Bz_dir=Vector(src.Bz_dir.x,src.Bz_dir.y,src.Bz_dir.z);
	Bz_distr=new Distribution2D((int)src.nbDistr_MAG.z);
	*Bz_distr=*(src.Bz_distr);
	this->Bz_mode=src.Bz_mode;
	this->Bz_period=src.Bz_period;
	this->Bz_phase=src.Bz_phase;

	this->nbDistr_BXY=src.nbDistr_BXY;
	beta_x_distr=new Distribution2D(src.nbDistr_BXY);
	*beta_x_distr=*(src.beta_x_distr);
	beta_y_distr=new Distribution2D(src.nbDistr_BXY);
	*beta_y_distr=*(src.beta_y_distr);
	eta_distr=new Distribution2D(src.nbDistr_BXY);
	*eta_distr=*(src.eta_distr);
	etaprime_distr=new Distribution2D(src.nbDistr_BXY);
	*etaprime_distr=*(src.etaprime_distr);
	coupling_distr=new Distribution2D(src.nbDistr_BXY);
	*coupling_distr=*(src.coupling_distr);
	e_spread_distr=new Distribution2D(src.nbDistr_BXY);
	*e_spread_distr=*(src.e_spread_distr);

	this->B_const=Vector(src.B_const.x,src.B_const.y,src.B_const.z);
	this->dL=src.dL;
	this->E=src.E;
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
	this->nbPoints=(int)src.nbPoints;
	this->nbDistr_MAG=Vector(src.nbDistr_MAG.x,src.nbDistr_MAG.y,src.nbDistr_MAG.z);
	this->beta_kind=src.beta_kind;
	this->AABBmin=Vector(src.AABBmin.x,src.AABBmin.y,src.AABBmin.z);
	this->AABBmax=Vector(src.AABBmax.x,src.AABBmax.y,src.AABBmax.z);
	this->fileName.assign(src.fileName);
	this->MAGXfileName.assign(src.MAGXfileName);
	this->MAGYfileName.assign(src.MAGYfileName);
	this->MAGZfileName.assign(src.MAGZfileName);
	this->BXYfileName.assign(src.BXYfileName);
}

Region& Region::operator=(const Region &src) {
	//alfa0=src.alfa0;
	betax=src.betax;
	betay=src.betay;

	Bx_dir=Vector(src.Bx_dir.x,src.Bx_dir.y,src.Bx_dir.z);
	Bx_distr=new Distribution2D((int)src.nbDistr_MAG.x);
	*Bx_distr=*(src.Bx_distr);
	this->Bx_mode=src.Bx_mode;
	this->Bx_period=src.Bx_period;
	this->Bx_phase=src.Bx_phase;

	By_dir=Vector(src.By_dir.x,src.By_dir.y,src.By_dir.z);
	By_distr=new Distribution2D((int)src.nbDistr_MAG.y);
	*By_distr=*(src.By_distr);
	this->By_mode=src.By_mode;
	this->By_period=src.By_period;
	this->By_phase=src.By_phase;

	Bz_dir=Vector(src.Bz_dir.x,src.Bz_dir.y,src.Bz_dir.z);
	Bz_distr=new Distribution2D((int)src.nbDistr_MAG.z);
	*Bz_distr=*(src.Bz_distr);
	this->Bz_mode=src.Bz_mode;
	this->Bz_period=src.Bz_period;
	this->Bz_phase=src.Bz_phase;

	this->nbDistr_BXY=src.nbDistr_BXY;
	beta_x_distr=new Distribution2D(src.nbDistr_BXY);
	*beta_x_distr=*(src.beta_x_distr);
	beta_y_distr=new Distribution2D(src.nbDistr_BXY);
	*beta_y_distr=*(src.beta_y_distr);
	eta_distr=new Distribution2D(src.nbDistr_BXY);
	*eta_distr=*(src.eta_distr);
	etaprime_distr=new Distribution2D(src.nbDistr_BXY);
	*etaprime_distr=*(src.etaprime_distr);
	coupling_distr=new Distribution2D(src.nbDistr_BXY);
	*coupling_distr=*(src.coupling_distr);
	e_spread_distr=new Distribution2D(src.nbDistr_BXY);
	*e_spread_distr=*(src.e_spread_distr);

	this->B_const=Vector(src.B_const.x,src.B_const.y,src.B_const.z);
	this->dL=src.dL;
	this->E=src.E;
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
	//this->Points=src.Points; //trajectory points will be copied in the next step
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
	this->fileName.assign(src.fileName);
	this->MAGXfileName.assign(src.MAGXfileName);
	this->MAGYfileName.assign(src.MAGYfileName);
	this->MAGZfileName.assign(src.MAGZfileName);
	this->BXYfileName.assign(src.BXYfileName);
	return *this;
}

int Region::LoadBXY(FileReader *file,Distribution2D *beta_x_distr,Distribution2D *beta_y_distr,
					Distribution2D *eta_distr,Distribution2D *etaprime_distr,
					Distribution2D *coupling_distr,	Distribution2D *e_spread_distr)
{
	//file->wasLineEnd=false;
	int nbDistr_BXY=file->ReadInt();
	beta_kind=file->ReadInt();

	*beta_x_distr=   Distribution2D(nbDistr_BXY);
	*beta_y_distr=   Distribution2D(nbDistr_BXY);
	*eta_distr=      Distribution2D(nbDistr_BXY);
	*etaprime_distr= Distribution2D(nbDistr_BXY);
	*coupling_distr= Distribution2D(nbDistr_BXY);
	*e_spread_distr= Distribution2D(nbDistr_BXY);

	/*char* tmp=file->ReadWord();
	if (!file->wasLineEnd) {
		sscanf(tmp,"%d",&beta_kind);
		beta_x_distr->valuesX[0]=beta_y_distr->valuesX[0]=eta_distr->valuesX[0]=etaprime_distr->valuesX[0]=file->ReadDouble();
	}
	else {
		sscanf(tmp,"%g",&(etaprime_distr->valuesX[0]));
		beta_x_distr->valuesX[0]=beta_y_distr->valuesX[0]=eta_distr->valuesX[0]=etaprime_distr->valuesX[0];
	}

	double value;
	for (int i=0;i<nbDistr_BXY;i++)
	{
		beta_x_distr->valuesY[i]=file->ReadDouble();
		file->wasLineEnd=false;
		beta_y_distr->valuesY[i]=file->ReadDouble();
		if (!file->IsEof()) {
			value=file->ReadDouble();
			if (!file->wasLineEnd)
			{
				eta_distr->valuesY[i]=value;
				etaprime_distr->valuesY[i]=file->ReadDouble();
				if (i<(nbDistr_BXY-1)) beta_x_distr->valuesX[i+1]=beta_y_distr->valuesX[i+1]=eta_distr->valuesX[i+1]=etaprime_distr->valuesX[i+1]=file->ReadDouble();
			} else
				if (i<(nbDistr_BXY-1)) beta_x_distr->valuesX[i+1]=beta_y_distr->valuesX[i+1]=eta_distr->valuesX[i+1]=etaprime_distr->valuesX[i+1]=value;
		}
	}*/
	for (int i=0;i<nbDistr_BXY;i++)
	{
		beta_x_distr->valuesX[i]=beta_y_distr->valuesX[i]=eta_distr->valuesX[i]=etaprime_distr->valuesX[i]=
			coupling_distr->valuesX[i]=e_spread_distr->valuesX[i]=file->ReadDouble();
		beta_x_distr->valuesY[i]=file->ReadDouble();
		beta_y_distr->valuesY[i]=file->ReadDouble();
		eta_distr->valuesY[i]=file->ReadDouble();
		etaprime_distr->valuesY[i]=file->ReadDouble();
		coupling_distr->valuesY[i]=file->ReadDouble();
		e_spread_distr->valuesY[i]=file->ReadDouble();
	}
	return nbDistr_BXY;
}

void Region::ExportPoints(FileWriter *file,GLProgress *prg,int frequency,BOOL doFullScan){

	char delimiter[3];
	char *ext;
	ext = strrchr(file->GetName(),'.');
	BOOL isTXT = _stricmp(ext,"txt")==0;
	sprintf(delimiter,(isTXT)?"\t":",");

	static char* varNames[] = {
		"#",
		"Orbit_X (cm)",
		"Orbit_Y (cm)",
		"Orbit_Z (cm)",
		"Traj. alpha (rad)",
		"Traj. theta (rad)",
		"Natural_divx (rad)",
		"Natural_divy (rad)",
		"Offset_X (cm)",
		"Offset_Y (cm)",
		"Offset_divx (rad)",
		"Offset_divy (rad)",
		"L (cm)",
		"Radius (cm)",
		"B (Tesla)",
		"B_X (Tesla)",
		"B_Y (Tesla)",
		"B_Z (Tesla)",
		"Critical Energy (eV)",
		"E / E_critical",
		"G1/H2",
		"B_factor",
		"B_factor_power",
		"SR_flux",
		"SR_power",
		"Beta_X",
		"Beta_Y",
		"Eta",
		"Eta_Prime",
		"Coupling",
		"Energy_Spread",
		"Emittance",
		"Emittance_X",
		"Emittance_Y",
		"Sigma_X",
		"Sigma_Y",
		"Sigma_X_Prime",
		"Sigma_Y_Prime",
	}; //38 variables

	static char* varNamesFullScan[] = {
		"SR_flux.Integral",
		"SR_power.Integral",
		"B_X.Integral",
		"B_Y.Integral",
		"B_Z.Integral"
	}; //5 variables

	//write column header
	for (int i=0;i<38;i++) {
		file->Write(varNames[i]);
		file->Write(delimiter);
	}

	if (doFullScan) {
		for (int i=0;i<5;i++) {
			file->Write(varNamesFullScan[i]);
			file->Write(delimiter);
		}
	}

	//creating distributions
	Distribution2D K_1_3_distribution=Generate_K_Distribution(1.0/3.0);
	Distribution2D K_2_3_distribution=Generate_K_Distribution(2.0/3.0);
	Distribution2D integral_N_photons=Generate_Integral(LOWER_LIMIT,UPPER_LIMIT,INTEGRAL_MODE_N_PHOTONS);
	Distribution2D integral_SR_power=Generate_Integral(LOWER_LIMIT,UPPER_LIMIT,INTEGRAL_MODE_SR_POWER);
	Distribution2D polarization_distribution=Generate_Polarization_Distribution(true,true);
	Distribution2D g1h2_distribution=Generate_G1_H2_Distribution();

	double integ_BX=0.0;
	double integ_BY=0.0;
	double integ_BZ=0.0;
	double integ_flux=0.0;
	double integ_power=0.0;

	//initialize random number generator
	extern SynRad *theApp;
	DWORD seed = (DWORD)((int)(theApp->GetTick()*1000.0)*_getpid());
	rseed (seed );

	for (int i=0;i<(int)Points.size();i++) {
		if (doFullScan || ((i%frequency)==0)) {

			//progressbar
			int percent=(int)(100.0*(double)i/(double)Points.size());
			prg->SetProgress((double)percent/100.0);

			Trajectory_Point *P=&Points[i];

			//Generate a photon (Radiate routine)
			//---------------------------------------------------
			int sourceId=i;
			Region *current_region=this;

			Trajectory_Point *source=&(current_region->Points[sourceId]);
			static double last_critical_energy,last_Bfactor,last_Bfactor_power; //to speed up calculation if critical energy didn't change
			static double last_average_ans;
			//actual variables: variables at actual point
			double average_;
			double B_factor,B_factor_power;
			double x_offset,y_offset,divx_offset,divy_offset,alpha,theta;
			double emittance_x_,emittance_y_;
			double betax_,betay_,eta_,etaprime_,coupling_,e_spread_;
			double sigmax,sigmay,sigmaxprime,sigmayprime;
			double critical_energy,generated_energy;
			double SR_flux,SR_power,g1h2,natural_divx,natural_divy,radius;
			Vector B_local=Vector(0.0,0.0,0.0);

			//init
			last_critical_energy=last_Bfactor=last_Bfactor_power=last_average_ans=average_=B_factor=
				B_factor_power=x_offset=y_offset=divx_offset=divy_offset=betax_=betay_=
				eta_=etaprime_=sigmax=sigmay=sigmaxprime=sigmayprime=critical_energy=generated_energy=
				SR_flux=SR_power=g1h2=natural_divx=natural_divy=radius=coupling_=e_spread_=
				emittance_x_=emittance_y_=0.0;

			//trajectory direction angles
			if (source->direction.z==0) {
				if (source->direction.x>=0) theta=-PI/2;
				else theta=PI/2;
			} else theta=-atan(source->direction.x/source->direction.z);
			alpha=-asin(source->direction.y);

			if (current_region->emittance>0.0) { //if beam is not ideal
				//calculate non-ideal beam's offset (the four sigmas)

				if (current_region->betax<0.0) { //negative betax value: load BXY file
					double coordinate; //interpolation X value (first column of BXY file)
					if (current_region->beta_kind==0) coordinate=sourceId*current_region->dL;
					else if (current_region->beta_kind==1) coordinate=source->position.x;
					else if (current_region->beta_kind==2) coordinate=source->position.y;
					else if (current_region->beta_kind==3) coordinate=source->position.z;

					betax_=current_region->beta_x_distr->InterpolateY(coordinate);
					betay_=current_region->beta_y_distr->InterpolateY(coordinate);
					eta_=current_region->eta_distr->InterpolateY(coordinate);
					etaprime_=current_region->etaprime_distr->InterpolateY(coordinate);
					coupling_=current_region->coupling_distr->InterpolateY(coordinate);
					e_spread_=current_region->e_spread_distr->InterpolateY(coordinate);
					//the six above distributions are the ones that are read from the BXY file
					//interpolateY finds the Y value corresponding to X input

				} else { //if no BXY file, use average values
					betax_=current_region->betax;
					betay_=current_region->betay;
					eta_=current_region->eta;
					etaprime_=current_region->etaprime;
					coupling_=current_region->coupling;
				}

				emittance_x_=current_region->emittance/(1.0+coupling_);
				emittance_y_=emittance_x_*coupling_;

				sigmaxprime=sqrt(emittance_x_/betax_+Sqr(etaprime_*e_spread_));
				//{ hor lattice-dependent divergence, radians }
				sigmayprime=sqrt(emittance_y_/betay_);
				//{ same for vertical divergence, radians }
				sigmax=sqrt(emittance_x_*betax_+Sqr(eta_*e_spread_));
				//{ hor lattice-dependent beam dimension, cm }
				sigmay=sqrt(emittance_y_*betay_);
				//{ same for vertical, cm }

				x_offset=Gaussian(sigmax);
				y_offset=Gaussian(sigmay);
				divx_offset=Gaussian(sigmaxprime);
				divy_offset=Gaussian(sigmayprime);

				//Calculate local base vectors
				Vector Z_local=source->direction.Normalize(); //Z' base vector
	/*if (source->rho.Norme()<1E29) {
		X_local=source->rho.Normalize(); //X' base vector
	} else {
		X_local=RandomPerpendicularVector(Z_local,1.0);
		}*/
				Vector Y_local=Vector(0.0,1.0,0.0); //same as absolute Y - assuming that machine's orbit is in the XZ plane
				Vector X_local=Crossproduct(Z_local,Y_local);

				Vector offset_pos=source->position; //choose ideal beam as origin
				offset_pos=Add(offset_pos,ScalarMult(X_local,x_offset)); //apply dX offset
				offset_pos=Add(offset_pos,ScalarMult(Y_local,y_offset)); //apply dY offset

				B_local=current_region->B(offset_pos); //recalculate B at offset position
				double B_parallel=DotProduct(source->direction,B_local);
				double B_orthogonal=sqrt(Sqr(B_local.Norme())-Sqr(B_parallel));

				if (B_orthogonal>VERY_SMALL)
					radius=current_region->E/0.00299792458/B_orthogonal;
				else radius=1.0E30;

				critical_energy=2.959E-5*pow(current_region->gamma,3)/radius; //becomes ~1E-30 if radius is 1E30
			} else {
				//0 emittance, ideal beam
				critical_energy=source->Critical_Energy(current_region->gamma);
				radius=source->rho.Norme();
				B_local=current_region->B(P->position);
			}

			generated_energy=SYNGEN1(current_region->energy_low/critical_energy,current_region->energy_hi/critical_energy,
				current_region->generation_mode);

			if (critical_energy==last_critical_energy)
				B_factor=last_Bfactor;
			else B_factor=(exp(integral_N_photons.InterpolateY(log(current_region->energy_hi/critical_energy)))
				-exp(integral_N_photons.InterpolateY(log(current_region->energy_low/critical_energy))))
				/exp(integral_N_photons.valuesY[integral_N_photons.size-1]);

			SR_flux=current_region->dL/radius/(2*PI)*current_region->gamma*4.1289E14*B_factor;

			//if (generation_mode==SYNGEN_MODE_POWERWISE) {
			if (critical_energy==last_critical_energy)
				B_factor_power=last_Bfactor_power;
			else B_factor_power=(exp(integral_SR_power.InterpolateY(log(current_region->energy_hi/critical_energy)))
				-exp(integral_SR_power.InterpolateY(log(current_region->energy_low/critical_energy))))
				/exp(integral_SR_power.valuesY[integral_SR_power.size-1]);

			if (critical_energy==last_critical_energy)
				average_=last_average_ans;
			else average_=integral_N_photons.Interval_Mean(current_region->energy_low/critical_energy,current_region->energy_hi/critical_energy);

			double average=integral_N_photons.average;
			double average1=integral_N_photons.average1;

			if (current_region->generation_mode==SYNGEN_MODE_POWERWISE)
				SR_flux=SR_flux/generated_energy*average_;

			double f=polarization_distribution.InterpolateY(generated_energy);
			g1h2=exp(g1h2_distribution.InterpolateY(generated_energy));
			double f_times_g1h2=f*g1h2;
			natural_divy=find_psi(generated_energy,Sqr(current_region->gamma),f_times_g1h2,
				current_region->enable_par_polarization,current_region->enable_ort_polarization)/current_region->gamma;
			natural_divx=find_chi(natural_divy,current_region->gamma,f_times_g1h2,
				current_region->enable_par_polarization,current_region->enable_ort_polarization); //divided by sHandle->gamma inside the function


			if (rnd()<0.5) natural_divx=-natural_divx;
			if (rnd()<0.5) natural_divy=-natural_divy;

			if (B_factor>0.0 && average_>VERY_SMALL) {
				SR_power=SR_flux*generated_energy*critical_energy*1.602189E-19*average/average_*B_factor_power/B_factor;
			} else SR_power=0.0;

			last_critical_energy=critical_energy;
			last_Bfactor=B_factor;
			last_Bfactor_power=B_factor_power;
			last_average_ans=average_;

			//---------------------------------------------------
			//Photon generated

			integ_flux+=SR_flux;
			integ_power+=SR_power;
			integ_BX+=B_local.x;
			integ_BY+=B_local.y;
			integ_BZ+=B_local.z;

			if ((i%frequency)==0) {
				file->Write("\n");

				file->WriteInt(i+1);file->Write(delimiter);
				file->WriteDouble(P->position.x);file->Write(delimiter);
				file->WriteDouble(P->position.y);file->Write(delimiter);
				file->WriteDouble(P->position.z);file->Write(delimiter);
				file->WriteDouble(alpha);file->Write(delimiter);
				file->WriteDouble(theta);file->Write(delimiter);
				file->WriteDouble(natural_divx);file->Write(delimiter);
				file->WriteDouble(natural_divy);file->Write(delimiter);
				file->WriteDouble(x_offset);file->Write(delimiter);
				file->WriteDouble(y_offset);file->Write(delimiter);
				file->WriteDouble(divx_offset);file->Write(delimiter);
				file->WriteDouble(divy_offset);file->Write(delimiter);
				file->WriteDouble(i*this->dL);file->Write(delimiter);
				file->WriteDouble(radius);file->Write(delimiter);
				file->WriteDouble(B_local.Norme());file->Write(delimiter);
				file->WriteDouble(B_local.x);file->Write(delimiter);
				file->WriteDouble(B_local.y);file->Write(delimiter);
				file->WriteDouble(B_local.z);file->Write(delimiter);
				file->WriteDouble(critical_energy);file->Write(delimiter);
				file->WriteDouble(generated_energy);file->Write(delimiter);
				file->WriteDouble(g1h2);file->Write(delimiter);
				file->WriteDouble(B_factor);file->Write(delimiter);
				file->WriteDouble(B_factor_power);file->Write(delimiter);
				file->WriteDouble(SR_flux);file->Write(delimiter);
				file->WriteDouble(SR_power);file->Write(delimiter);
				file->WriteDouble(betax_);file->Write(delimiter);
				file->WriteDouble(betay_);file->Write(delimiter);
				file->WriteDouble(eta_);file->Write(delimiter);
				file->WriteDouble(etaprime_);file->Write(delimiter);
				file->WriteDouble(coupling_);file->Write(delimiter);
				file->WriteDouble(e_spread_);file->Write(delimiter);
				file->WriteDouble(emittance);file->Write(delimiter);
				file->WriteDouble(emittance_x_);file->Write(delimiter);
				file->WriteDouble(emittance_y_);file->Write(delimiter);
				file->WriteDouble(sigmax);file->Write(delimiter);
				file->WriteDouble(sigmay);file->Write(delimiter);
				file->WriteDouble(sigmaxprime);file->Write(delimiter);
				file->WriteDouble(sigmayprime);file->Write(delimiter);
				
				if (doFullScan) {
					file->WriteDouble(integ_flux);file->Write(delimiter);
					file->WriteDouble(integ_power);file->Write(delimiter);
					file->WriteDouble(integ_BX);file->Write(delimiter);
					file->WriteDouble(integ_BY);file->Write(delimiter);
					file->WriteDouble(integ_BZ);file->Write(delimiter);
				}
			}
		}
	}
	prg->SetMessage("Writing file to disk...");
}


void Region::SaveParam(FileWriter *file) {
	file->Write("param_file_version:");file->WriteInt(PARAMVERSION,"\n");
	file->Write("startPos_cm:");file->WriteDouble(startPoint.x);file->WriteDouble(startPoint.y);file->WriteDouble(startPoint.z,"\n");
	file->Write("startDir_cm:");file->WriteDouble(startDir.x);file->WriteDouble(startDir.y);file->WriteDouble(startDir.z,"\n");
	file->Write("dL_cm:");file->WriteDouble(dL,"\n");
	file->Write("boundaries_cm:");file->WriteDouble(limits.x);file->WriteDouble(limits.y);file->WriteDouble(limits.z,"\n");
	file->Write("particleMass_GeV:");file->WriteDouble(particleMass,"\n");
	file->Write("beamEnergy_GeV:");file->WriteDouble(E,"\n");
	file->Write("emittance:");file->WriteDouble(emittance,"\n");
	if (emittance!=0.0) { //non-ideal beam
		file->Write("beta_x:");file->WriteDouble(betax,"\n");
		if (betax<0.0) {//use BXY file
			file->Write("BXYfileName:\"");file->Write((char*)SplitFilename(BXYfileName).c_str());file->Write("\"\n"); //truncate path
		} else {//constants
			file->Write("beta_y:");file->WriteDouble(betay,"\n");
			file->Write("eta:");file->WriteDouble(eta,"\n");
			file->Write("eta_prime:");file->WriteDouble(etaprime,"\n");
			file->Write("coupling:");file->WriteDouble(coupling,"\n");
			file->Write("energy_spread:");file->WriteDouble(energy_spread,"\n");
		}
	}
	file->Write("E_min_eV:");file->WriteDouble(energy_low,"\n");
	file->Write("E_max_eV:");file->WriteDouble(energy_hi,"\n");
	file->Write("enable_par_polarization:");file->WriteInt(enable_par_polarization,"\n");
	file->Write("enable_ort_polarization:");file->WriteInt(enable_ort_polarization,"\n");
	file->Write("psiMax_rad:");file->WriteDouble(psimaxX);file->WriteDouble(psimaxY,"\n");
	
	file->Write("Bx_mode:");file->WriteInt(Bx_mode,"\n");
	if (Bx_mode==B_MODE_CONSTANT) {
		file->Write("Bx_const_Tesla:");file->WriteDouble(B_const.x,"\n");
	} else {
		file->Write("Bx_fileName:\"");file->Write((char*)SplitFilename(MAGXfileName).c_str());file->Write("\"\n"); //truncate path
	}
	file->Write("By_mode:");file->WriteInt(By_mode,"\n");
	if (By_mode==B_MODE_CONSTANT) {
		file->Write("By_const_Tesla:");file->WriteDouble(B_const.y,"\n");
	} else {
		file->Write("By_fileName:\"");file->Write((char*)SplitFilename(MAGYfileName).c_str());file->Write("\"\n"); //truncate path
	}
	file->Write("Bz_mode:");file->WriteInt(Bz_mode,"\n");
	if (Bz_mode==B_MODE_CONSTANT) {
		file->Write("Bz_const_Tesla:");file->WriteDouble(B_const.z,"\n");
	} else {
		file->Write("Bz_fileName:\"");file->Write((char*)SplitFilename(MAGZfileName).c_str());file->Write("\"\n"); //truncate path
	}
}

void Region::LoadParam(FileReader *file,GLProgress *prg){

	//Creating references to X,Y,Z components (to avoid writing code 3 times)
	double* Bconst_ptr[3]={&B_const.x,&B_const.y,&B_const.z};
	Vector* Bdir_ptr[3]={&Bx_dir,&By_dir,&Bz_dir};
	int* Bmode_ptr[3]={&Bx_mode,&By_mode,&Bz_mode};
	double* Bperiod_ptr[3]={&Bx_period,&By_period,&Bz_period};
	double* Bphase_ptr[3]={&Bx_phase,&By_phase,&Bz_phase};
	Distribution2D* distr_ptr[3]={Bx_distr,By_distr,Bz_distr};
	std::string* MagFileNamePtr[3]={&MAGXfileName,&MAGYfileName,&MAGZfileName};
	std::string Bmode_str[3]={"Bx_mode","By_mode","Bz_mode"};
	std::string Bconst_str[3]={"Bx_const_Tesla","By_const_Tesla","Bz_const_Tesla"};
	std::string MagFileNameStr[3]={"Bx_fileName","By_fileName","Bz_fileName"};
	int paramVersion;

	prg->SetMessage("Reading parameter file...");

	file->ReadKeyword("param_file_version");file->ReadKeyword(":");paramVersion=file->ReadInt();
	if (paramVersion>PARAMVERSION) {
		Error err=Error("Param file version not supported");
		throw err;
		return;
	}
	
	file->ReadKeyword("startPos_cm");file->ReadKeyword(":");
	double x=file->ReadDouble();double y=file->ReadDouble();double z=file->ReadDouble(); startPoint=Vector(x,y,z);
	file->ReadKeyword("startDir_cm");file->ReadKeyword(":");
	x=file->ReadDouble(); y=file->ReadDouble(); z=file->ReadDouble(); startDir=Vector(x,y,z);
	file->ReadKeyword("dL_cm");file->ReadKeyword(":");dL=file->ReadDouble();
	file->ReadKeyword("boundaries_cm");file->ReadKeyword(":");
	x=file->ReadDouble(); y=file->ReadDouble(); z=file->ReadDouble(); limits=Vector(x,y,z);
	file->ReadKeyword("particleMass_GeV");file->ReadKeyword(":");particleMass=file->ReadDouble();
	file->ReadKeyword("beamEnergy_GeV");file->ReadKeyword(":");E=file->ReadDouble();

	emittance=eta=etaprime=coupling=energy_spread=betax=betay=0.0; //default values
	nbDistr_BXY=0;

	file->ReadKeyword("emittance");file->ReadKeyword(":");emittance=file->ReadDouble();
	if (emittance!=0.0) { //non-ideal beam
		file->ReadKeyword("beta_x");file->ReadKeyword(":");betax=file->ReadDouble();
		beta_kind=0;
		if (betax<0.0) {//use BXY file
			file->ReadKeyword("BXYfileName");file->ReadKeyword(":");
			
				char tmp[512];
				char tmp2[512];
				char CWD [MAX_PATH];
				_getcwd( CWD, MAX_PATH );
				strcpy(tmp2,file->ReadString()); //get BXY file name
				sprintf(tmp,"%s\\tmp\\%s",CWD,tmp2); //look for file in tmp directory (extracted from syn7z)
				if (FileUtils::Exist(tmp)) //file found in tmp directory
					BXYfileName.assign(tmp);
				else {//not found, look for it in current directory (syn files)				
					std::string path=SplitPath(file->GetName());
					sprintf(tmp,"%s\\%s",path.c_str(),tmp2);
					if (FileUtils::Exist(tmp)) //file found in current directory
						BXYfileName.assign(tmp);
					else {//not in tmp, nor in current, error.
						sprintf(tmp2,"Referenced BXY file not found:\n%s",tmp);
						throw Error(tmp2);
					}
				}
				FileReader BXYfile((char*)BXYfileName.c_str());
				nbDistr_BXY=LoadBXY(&BXYfile,beta_x_distr,beta_y_distr,eta_distr,etaprime_distr,coupling_distr,e_spread_distr);
			 
		} else { //constants
			file->ReadKeyword("beta_y");file->ReadKeyword(":");betay=file->ReadDouble();
			file->ReadKeyword("eta");file->ReadKeyword(":");eta=file->ReadDouble();
			file->ReadKeyword("eta_prime");file->ReadKeyword(":");etaprime=file->ReadDouble();
			file->ReadKeyword("coupling");file->ReadKeyword(":");coupling=file->ReadDouble();
			file->ReadKeyword("energy_spread");file->ReadKeyword(":");energy_spread=file->ReadDouble();
		}
	}
	file->ReadKeyword("E_min_eV");file->ReadKeyword(":");energy_low=file->ReadDouble();
	file->ReadKeyword("E_max_eV");file->ReadKeyword(":");energy_hi=file->ReadDouble();
	file->ReadKeyword("enable_par_polarization");file->ReadKeyword(":");enable_par_polarization=file->ReadInt();
	file->ReadKeyword("enable_ort_polarization");file->ReadKeyword(":");enable_ort_polarization=file->ReadInt();
	file->ReadKeyword("psiMax_rad");file->ReadKeyword(":");psimaxX=file->ReadDouble();psimaxY=file->ReadDouble();

	for (int i=0;i<3;i++) {
		file->ReadKeyword((char*)Bmode_str[i].c_str());file->ReadKeyword(":");*Bmode_ptr[i]=file->ReadInt();
		if (*Bmode_ptr[i]==B_MODE_CONSTANT) {
			file->ReadKeyword((char*)Bconst_str[i].c_str());file->ReadKeyword(":");
			*Bconst_ptr[i]=file->ReadDouble();
		} else {
			file->ReadKeyword((char*)MagFileNameStr[i].c_str());file->ReadKeyword(":");
			
				char tmp[512];
				char tmp2[512];
				char CWD [MAX_PATH];
				_getcwd( CWD, MAX_PATH );
				strcpy(tmp2,file->ReadString()); //get BXY file name
				sprintf(tmp,"%s\\tmp\\%s",CWD,tmp2); //look for file in tmp directory (extracted from syn7z)
				if (FileUtils::Exist(tmp)) //file found in tmp directory
					MagFileNamePtr[i]->assign(tmp);
				else {//not found, look for it in current directory (syn files)
					std::string path=SplitPath(file->GetName());
					sprintf(tmp,"%s\\%s",path.c_str(),tmp2);
					if (FileUtils::Exist(tmp)) //file found in current directory
						MagFileNamePtr[i]->assign(tmp);
					else {//not in tmp, nor in current, error.
						sprintf(tmp2,"Referenced MAG file not found: %s",tmp);
						throw Error(tmp2);
					}
				}
				FileReader MAGfile((char*)MagFileNamePtr[i]->c_str());
				*(distr_ptr[i])=LoadMAGFile(&MAGfile,Bdir_ptr[i],Bperiod_ptr[i],Bphase_ptr[i],*Bmode_ptr[i]);
			
		}
	}
	gamma=abs(E/particleMass);

	prg->SetMessage("Calculating trajectory...");
	Points=CalculateTrajectory(dL,E,limits,startPoint,startDir,1000000); //max 1 million points
	isLoaded=true;
	extern SynRad *theApp;
	SynRad *mApp = (SynRad *)theApp;
	if (mApp->regionInfo) mApp->regionInfo->Update();
}
