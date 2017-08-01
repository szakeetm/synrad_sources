#include "Region_full.h"
#include "SynradGeometry.h"
#include "SynradTypes.h"
#include "Random.h"
#include "GLApp/MathTools.h"
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

extern SynRad *mApp;

void Region_full::CalculateTrajectory(int max_steps){
	//global variables. All distances in cm!
	Trajectory_Point current_point;
	//extern Distribution2D polarization_distribution,integral_N_photons,integral_SR_power,g1h2_distribution;

	//initial position and speed
	current_point.position=this->params.startPoint;
	current_point.direction=this->params.startDir;
	//current_point.rho=Vector3d(0.0,0.0,1e30); //big enough so no rotation
	AABBmin=AABBmax=this->params.startPoint;

	//Calculate trajectory
	Points.push_back(current_point); //beam starting position
	CalcPointProperties(0); //calculate magnetic field, etc. of first point

	for (int i=0;i<max_steps&&!isOutsideBoundaries(current_point.position,i==0);i++) {
		current_point=OneStep(i); //step forward
		Points.push_back(current_point); //store point on i+1st place
		CalcPointProperties(i + 1);
		//extend trajectory limits if needed
		if (current_point.position.x<AABBmin.x) AABBmin.x = current_point.position.x;
		if (current_point.position.y<AABBmin.y) AABBmin.y = current_point.position.y;
		if (current_point.position.z<AABBmin.z) AABBmin.z = current_point.position.z;
		if (current_point.position.x>AABBmax.x) AABBmax.x = current_point.position.x;
		if (current_point.position.y>AABBmax.y) AABBmax.y = current_point.position.y;
		if (current_point.position.z>AABBmax.z) AABBmax.z = current_point.position.z;
	}
}

Trajectory_Point Region_full::OneStep(int pointId) {
	Trajectory_Point* p0 = &(Points[pointId]);

	Vector3d rotation_axis = CrossProduct(p0->direction, p0->rho); //Rotate function will normalize it
	//Calculate new point
	Trajectory_Point p;
	//p.position = Rotate(p0->position, p0->position + p0->rho, rotation_axis, params.dL_cm / p0->rho.Norme());//
	p.position = p0->position + p0->direction*params.dL_cm;
	p.direction = Rotate(p0->direction, Vector3d(0,0,0), rotation_axis, params.dL_cm / p0->rho.Norme()).Normalized(); //Renormalize to prevent accumulating rounding errors
	/*
	p.direction = Rotate(p0->direction, Vector3d(0, 0, 0), rotation_axis, params.dL_cm / p0->rho.Norme());
	p.position = p0->position + p.direction * (params.dL_cm / p.direction.Norme());
	*/
	return p;
}

void Region_full::CalcPointProperties(int pointId) {
	Trajectory_Point* p = &(Points[pointId]);

	p->B = B(pointId, Vector3d(0, 0, 0)); //local magnetic field
	if (p->B.Norme() < VERY_SMALL) {//no magnetic field, beam goes in straight line
		p->rho = Vector3d(0.0, 0.0, 1.0e30); //no rotation
	}
	else {
		double B_parallel = Dot(p->direction, p->B);
		double B_orthogonal = sqrt(Sqr(p->B.Norme()) - Sqr(B_parallel));
		Vector3d lorentz_force_direction = CrossProduct(p->direction, p->B).Normalized();
		if (params.particleMass_GeV<0.0) lorentz_force_direction = -1.0 * lorentz_force_direction;
		double radius;
		if (B_orthogonal>VERY_SMALL)
			radius = params.E_GeV / 0.00299792458 / B_orthogonal; //in centimeters
		else radius = 1.0E30;
		p->rho = lorentz_force_direction * radius;
	}

	p->critical_energy = p->Critical_Energy(params.gamma);

	//Calculate local base vectors
	p->Z_local = p->direction.Normalized(); //Z' base vector
	p->Y_local = Vector3d(0.0, 1.0, 0.0); //same as absolute Y - assuming that machine's orbit is in the XZ plane
	p->X_local = CrossProduct(p->Y_local, p->Z_local); //Left-hand coord system

	if (params.emittance_cm > 0.0) { //if beam is not ideal
		//calculate non-ideal beam's offset (the four sigmas)
		if (params.betax_const_cm < 0.0) { //negative betax value: load BXY file
			double coordinate; //interpolation X value (first column of BXY file)
			if (params.beta_kind == 0) coordinate = pointId*params.dL_cm;
			else if (params.beta_kind == 1) coordinate = p->position.x;
			else if (params.beta_kind == 2) coordinate = p->position.y;
			else if (params.beta_kind == 3) coordinate = p->position.z;

			std::vector<double> betaValues = betaFunctions.GetYValues(coordinate); //interpolation

			p->beta_X = betaValues[0];    // [cm]
			p->beta_Y = betaValues[1];    // [cm]
			p->eta = betaValues[2];       // [cm]
			p->eta_prime = betaValues[3]; // [derivative -> dimensionless]
			p->alpha_X = betaValues[4];   // [derivative -> dimensionless]
			p->alpha_Y = betaValues[5];   // [derivative -> dimensionless]
			//the six above distributions are the ones that are read from the BXY file
			//GetValuesAt finds the Y values corresponding to X input

		}
		else { //if no BXY file, use average values
			p->beta_X = params.betax_const_cm;
			p->beta_Y = params.betay_const_cm;
			p->eta = params.eta_x_const_cm;
			p->eta_prime = params.eta_x_prime_const;
			//If beta is constant, alpha, which is its derivative, must be 0
			p->alpha_X = 0.0;
			p->alpha_Y = 0.0;
		}

		p->gamma_X = (1.0 + Sqr(p->alpha_X)) / p->beta_X; // [1/cm, hundred times smaller than if it was in 1/m]
		p->gamma_Y = (1.0 + Sqr(p->alpha_Y)) / p->beta_Y; // [1/cm, hundred times smaller than if it was in 1/m]

		p->emittance_X = params.emittance_cm / (1.0 + params.coupling_percent*0.01); // [cm]
		p->emittance_Y = p->emittance_X*params.coupling_percent*0.01;             // [cm]
		
		p->theta_X = 0.5 * atan(2 * p->alpha_X / (p->gamma_X - p->beta_X)); // xx' phase ellipse rotation [rad], about 100 times smaller than if X would be in m
		p->theta_Y = 0.5 * atan(2 * p->alpha_Y / (p->gamma_Y - p->beta_Y)); // yy' phase ellipse rotation [rad], about 100 times smaller than if X would be in m

		//Correcting parameters by energy spread:
		//The phase ellipse equation Beta*x'^2 + Gamma*x^2 + 2*Alpha*xx' = Emittance
		//Can be rewritten as 
		// <x^2>_rms / Emittance_rms * x'^2 + <x'^2>_rms / Emittance_rms * x^2 - 2 * <xx'>_rms / Emittance_rms * xx' = Emittance_rms
		// Where <x^2>_rms  = (sigma_x)^2  =   Beta  * Emittance_rms + (Eta_x  * dE/E)^2            := Beta_corr  * Emittance_rms
		// And   <x'^2>_rms = (sigma_x')^2 =   Gamma * Emittance_rms + (Eta_x' * dE/E)^2            := Gamma_corr * Emittance_rms
		// And   <xx'>_rms  =  sigma_xx'   = - Alpha * Emittance_rms +  Eta_x *  Eta_x' * (dE/E)^2  := Alpha_corr * Emittance_rms

		double alpha_corr = p->alpha_X  - p->eta*p->eta_prime*Sqr(params.energy_spread_percent*0.01)/p->emittance_X;
		double beta_corr =  p->beta_X   + Sqr(p->eta*params.energy_spread_percent*0.01)/p->emittance_X;
		double gamma_corr = p->gamma_X  + Sqr(p->eta_prime*params.energy_spread_percent*0.01)/p->emittance_X;

		p->a_x = sqrt(p->emittance_X / (gamma_corr*Sqr(cos(p->theta_X)) + 2 * alpha_corr*cos(p->theta_X)*sin(p->theta_X) + beta_corr*Sqr(sin(p->theta_X)))); //major axis of xx' phase ellipse (equal to sigma_x[cm] when alpha=theta=0)
		p->a_y = sqrt(p->emittance_Y / (p->gamma_Y*Sqr(cos(p->theta_Y)) + 2 * p->alpha_Y*cos(p->theta_Y)*sin(p->theta_Y) + p->beta_Y*Sqr(sin(p->theta_Y)))); //major axis of yy' phase ellipse (equal to sigma_y[rad] when alpha=theta=0)
		p->b_x = sqrt(p->emittance_X / (gamma_corr*Sqr(sin(p->theta_X)) - 2 * alpha_corr*cos(p->theta_X)*sin(p->theta_X) + beta_corr*Sqr(cos(p->theta_X)))); //minor axis of xx' phase ellipse (equal to sigma_x'[cm] when alpha=theta=0)
		p->b_y = sqrt(p->emittance_Y / (p->gamma_Y*Sqr(sin(p->theta_Y)) - 2 * p->alpha_Y*cos(p->theta_Y)*sin(p->theta_Y) + p->beta_Y*Sqr(cos(p->theta_Y)))); //minor axis of yy' phase ellipse (equal to sigma_y'[rad] when alpha=theta=0)

		p->sigma_x = sqrt(p->emittance_X*beta_corr);						//{ hor lattice-dependent beam dimension, cm }
		p->sigma_y = sqrt(p->emittance_Y*p->beta_Y);                                        						//{ same for vertical, cm }
		p->sigma_x_prime = sqrt(p->emittance_X * gamma_corr);		//{ hor lattice-dependent divergence, radians, same as if s and x would be in [m] }
		p->sigma_y_prime = sqrt(p->emittance_Y * p->gamma_Y); 														//{ same for vertical divergence, radians, same as if s and x would be in [m] }

	}
	else {
		//0 emittance, ideal beam
		p->emittance_X =
			p->emittance_Y =
			p->beta_X =
			p->beta_Y =
			p->eta =
			p->eta_prime =
			p->alpha_X =
			p->alpha_Y =
			p->gamma_X =
			p->gamma_Y =
			p->sigma_x =
			p->sigma_y =
			p->sigma_x_prime =
			p->sigma_y_prime =
			p->theta_X =
			p->theta_Y =
			p->a_x =
			p->b_x =
			p->a_y =
			p->b_y = 0.0;
	}
}

bool Region_full::isOutsideBoundaries(Vector3d a,bool recalcDirs){
	static double xDir=(params.limits.x>params.startPoint.x)?1.0:-1.0;
	static double yDir=(params.limits.y>params.startPoint.y)?1.0:-1.0;
	static double zDir=(params.limits.z>params.startPoint.z)?1.0:-1.0;
	if (recalcDirs) {
		xDir=(params.limits.x>params.startPoint.x)?1.0:-1.0;
		yDir=(params.limits.y>params.startPoint.y)?1.0:-1.0;
		zDir=(params.limits.z>params.startPoint.z)?1.0:-1.0;
	}
	return ((a.x- params.limits.x)*xDir>0.0)||((a.y- params.limits.y)*yDir>0.0)||((a.z- params.limits.z)*zDir>0.0);
}

void Region_full::LoadPAR(FileReader *file){

	//Creating references to X,Y,Z components (to avoid writing code 3 times)
	//double* result_ptr[3]={&result.x,&result.y,&result.z};
	double* Bconst_ptr[3]={&params.B_const.x,&params.B_const.y,&params.B_const.z};
	Vector3d* Bdir_ptr[3]={&params.Bx_dir,&params.By_dir,&params.Bz_dir};
	int* Bmode_ptr[3]={&params.Bx_mode,&params.By_mode,&params.Bz_mode};
	double* Bperiod_ptr[3]={&params.Bx_period,&params.By_period,&params.Bz_period};
	double* Bphase_ptr[3]={&params.Bx_phase,&params.By_phase,&params.Bz_phase};
	Distribution2D* distr_ptr[3]={&Bx_distr,&By_distr,&Bz_distr};
	std::string* MagFileNamePtr[3]={&MAGXfileName,&MAGYfileName,&MAGZfileName};

	//prg->SetMessage("Reading parameter file...");
	file->ReadInt(); //unused (no_scans)
	file->ReadInt(); //unused (nregions)
	file->ReadInt(); //unused (generation mode)
	file->JumpComment();

	params.particleMass_GeV=file->ReadDouble();
	params.E_GeV=file->ReadDouble();params.current_mA=1;
	params.energy_low_eV=file->ReadDouble();
	params.energy_hi_eV=file->ReadDouble();
	file->JumpComment();

	params.gamma=abs(params.E_GeV/ params.particleMass_GeV);

	params.dL_cm=file->ReadDouble();
	double x0=file->ReadDouble();
	double y0=file->ReadDouble();
	double z0=file->ReadDouble();
	params.startPoint=Vector3d(x0,y0,z0);

	double teta0=file->ReadDouble();
	double alfa0=file->ReadDouble();
	params.startDir=Vector3d(-cos(alfa0)*sin(teta0),-sin(alfa0),cos(alfa0)*cos(teta0));
	file->JumpComment();

	double xmax=file->ReadDouble();
	double ymax=file->ReadDouble();
	double zmax=file->ReadDouble();
	params.limits=Vector3d(xmax,ymax,zmax);
	file->JumpComment();

	int i_struct1=file->ReadInt(); //thrown away
	file->JumpComment();

	params.enable_par_polarization=(file->ReadInt()==1);
	params.enable_ort_polarization=(file->ReadInt()==1);

	//Polarization component selection
	if (params.enable_par_polarization && params.enable_ort_polarization)
		params.polarizationCompIndex = 0; //Full polarization
	else if (params.enable_par_polarization)
		params.polarizationCompIndex = 1; //Parallel polarization
	else
		params.polarizationCompIndex = 2; //Orthogonal polarization

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

	params.emittance_cm= params.eta_x_const_cm= params.eta_x_prime_const= params.energy_spread_percent= params.betax_const_cm= params.betay_const_cm=0.0; //default values
	params.coupling_percent=100.0;
	params.nbDistr_BXY=0;

	params.emittance_cm=file->ReadDouble();
	params.betax_const_cm=file->ReadDouble();
	if (!(params.betax_const_cm<0.0)) {
		file->wasLineEnd=false;
		params.betay_const_cm=file->ReadDouble();
		double value=file->ReadDouble();
		if (!file->wasLineEnd) {
			params.eta_x_const_cm=value;
			params.eta_x_prime_const=file->ReadDouble();
			value=file->ReadDouble();
			if (!file->wasLineEnd) {
				params.coupling_percent=value;
				params.energy_spread_percent=file->ReadDouble();
				value=file->ReadDouble();
			}
		}
		params.psimaxX_rad=value*0.001;
		file->wasLineEnd=false;
	} else params.psimaxX_rad=file->ReadDouble()*0.001;
	params.psimaxY_rad=file->ReadDouble()*0.001;
	file->JumpComment();

	params.beta_kind=0;
	if (params.betax_const_cm<0.0) { //read BXY file
		double fileDouble=-1.0*params.betax_const_cm;
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
			BXYfileName.assign(fileName);
			params.nbDistr_BXY=LoadBXY(fileName);
	}

	//prg->SetMessage("Calculating trajectory...");
	CalculateTrajectory(1000000); //max 1 million points
	isLoaded=true;
	if (mApp->regionInfo) mApp->regionInfo->Update();
}

Distribution2D Region_full::LoadMAGFile(FileReader *file,Vector3d *dir,double *period,double *phase,int mode){
	Distribution2D result(1);
	if (mode==B_MODE_QUADRUPOLE) {
		double quadrX=file->ReadDouble();
		double quadrY=file->ReadDouble();
		double quadrZ=file->ReadDouble();
		file->JumpComment();
		params.quad_params.center=Vector3d(quadrX,quadrY,quadrZ);

		params.quad_params.alfa_q=file->ReadDouble();
		params.quad_params.beta_q=file->ReadDouble();
		params.quad_params.rot_q=file->ReadDouble();
		file->JumpComment();

		params.quad_params.K_q=file->ReadDouble()/100.0; //T/m->T/cm conversion
		params.quad_params.L_q=file->ReadDouble();

		params.quad_params.cosalfa_q=cos(params.quad_params.alfa_q);
		params.quad_params.sinalfa_q=sin(params.quad_params.alfa_q);
		params.quad_params.cosbeta_q=cos(params.quad_params.beta_q);
		params.quad_params.sinbeta_q=sin(params.quad_params.beta_q);
		params.quad_params.cosrot_q=cos(params.quad_params.rot_q);
		params.quad_params.sinrot_q=sin(params.quad_params.rot_q);

		params.quad_params.direction=Vector3d(-params.quad_params.cosalfa_q*params.quad_params.sinbeta_q,
			-params.quad_params.sinalfa_q,
			params.quad_params.cosalfa_q*params.quad_params.cosbeta_q);

	} else {
		*period=file->ReadDouble();
		if (mode==B_MODE_HELICOIDAL) *phase=file->ReadDouble();
		file->JumpComment();

		double dirx=file->ReadDouble();
		double diry=file->ReadDouble();
		double dirz=file->ReadDouble();
		*dir=Vector3d(dirx,diry,dirz);
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

Region_full::Region_full():Region_mathonly(){

	//generation_mode=SYNGEN_MODE_POWERWISE;
	//emittance=eta=etaprime=energy_spread=betax=betay=0.0;
	//coupling=100.0;
	selectedPoint=-1;
	isLoaded=false;
	//object placeholders until MAG files are loaded

	/*Bx_distr = new Distribution2D(1);
	By_distr = new Distribution2D(1);
	Bz_distr = new Distribution2D(1);
	beta_x_distr = new Distribution2D(1);
	beta_y_distr = new Distribution2D(1);
	eta_distr = new Distribution2D(1);
	etaprime_distr = new Distribution2D(1);
	e_spread_distr = new Distribution2D(1);*/

	MAGXfileName="";
	MAGYfileName="";
	MAGZfileName="";
	BXYfileName="";

	//default values for empty region
	/*dL=0.01;
	limits=Vector3d(1000,1000,0.1);
	startPoint=Vector3d(0,0,0);
	startDir=Vector3d(0,0,1);
	particleMass=-0.0005110034; //electron
	E=1;current=1;
	betax=betay=eta=etaprime=energy_spread=emittance=0.0;
	coupling=100.0;
	energy_low=10;
	energy_hi=1e6;
	enable_ort_polarization=enable_par_polarization=true;
	psimaxX=psimaxY=PI;
	Bx_mode=By_mode=Bz_mode=B_MODE_CONSTANT;
	B_const=Vector3d(0,0,0);
	this->nbPoints=0;*/
}

Region_full::~Region_full(){
	/*Distribution2D* distr_ptr[8]={Bx_distr,By_distr,Bz_distr,beta_x_distr,beta_y_distr,eta_distr,
		etaprime_distr,e_spread_distr};
	for (int i=0;i<8;i++)
		SAFE_DELETE(distr_ptr[i]);
	Points=std::vector<Trajectory_Point>();*/
	//Already called by region_mathonly's destructor
}

void Region_full::Render(const int& regionId, const int& dispNumTraj, GLMATERIAL *B_material, const double& vectorLength) {
	if (!mApp->whiteBg) glPointSize(1.0f);
	else glPointSize(2.0f);
	//if (antiAliasing) glEnable(GL_POINT_SMOOTH);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);
	glDisable(GL_CULL_FACE);
	if (!mApp->whiteBg) glColor3f(1.0f, 0.9f, 0.2f);
	else glColor3f(1.0f,0.5f,0.2f);
	glBegin(GL_POINTS);
	//Region
	int N=(int)Points.size();
	int increment=Max(1,(int)(N/dispNumTraj));
	for (int i=0;i<N;i+=increment) {
		// Draw dot

		//if (hiddenVertex) glDisable(GL_DEPTH_TEST);

		Vector3d *v = &(Points[i].position);
		glVertex3d(v->x,v->y,v->z);
	}
	glEnd();

	//Selected trajectory point
	if (selectedPoint!=-1) {
		// Draw dot
		if (!mApp->whiteBg) glPointSize(6.0f);
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
		if (mApp->antiAliasing) {
			glEnable(GL_LINE_SMOOTH);
			//glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			//glColor4f(0.0f,0.0f,1.0f,0.5f);
			glEnable(GL_DEPTH_TEST);
		}
		glLineWidth(0.5f);

		Vector3d O=Points[selectedPoint].position;
		/*
		Vector3d B_local=Points[selectedPoint].B;
		Vector3d B_local_norm=B_local.Normalized();
		Vector3d Rho_local=Points[selectedPoint].rho;
		Vector3d Rho_local_norm=Rho_local.Normalized();
		Vector3d dir=Points[selectedPoint].direction;
		double factor=vectorLength;
		if (B_local.Norme()>0.0) 

			GLToolkit::DrawVector(O.x,O.y,O.z,O.x+factor*B_local_norm.x,O.y+factor*B_local_norm.y,O.z+factor*B_local_norm.z);
		if (Rho_local.Norme()>0.0) 
			GLToolkit::DrawVector(O.x,O.y,O.z,O.x+factor*Rho_local_norm.x,O.y+factor*Rho_local_norm.y,O.z+factor*Rho_local_norm.z);
		if (dir.Norme()>0.0) 
			GLToolkit::DrawVector(O.x,O.y,O.z,O.x+factor*dir.x,O.y+factor*dir.y,O.z+factor*dir.z);
			*/
		if (mApp->antiAliasing) glDisable(GL_LINE_SMOOTH);
		
		glPointSize(3.0f);
		glColor3f(0.5f,1.0f,1.0f);
		glBegin(GL_POINTS);
		glVertex3d(O.x,O.y,O.z);
		glEnd();

		//glEnable(GL_BLEND);
		char point_label[256];
		/*char B_label[128];
		char Rho_label[128];
		char dir_label[128];*/
		sprintf(point_label,"Region %d point #%d   [%g , %g , %g] L= %g",
			regionId + 1 ,	selectedPoint+1, O.x,O.y,O.z,selectedPoint*params.dL_cm);
		/*sprintf(B_label,"B (%g Tesla)",B_local.Norme());
		sprintf(Rho_label,"Rho (%g cm)",Rho_local.Norme());
		sprintf(dir_label,"Direction");*/
		GLToolkit::GetDialogFont()->SetTextColor(0.5f,0.6f,1.0f);
		GLToolkit::DrawStringInit();
		GLToolkit::DrawString((float)(O.x),(float)(O.y),(float)(O.z),point_label,GLToolkit::GetDialogFont());
		/*
		GLToolkit::DrawString((float)(O.x+factor*B_local_norm.x),(float)(O.y+factor*B_local_norm.y),(float)(O.z+factor*B_local_norm.z),B_label,GLToolkit::GetDialogFont());
		GLToolkit::DrawString((float)(O.x+factor*Rho_local_norm.x),(float)(O.y+factor*Rho_local_norm.y),(float)(O.z+factor*Rho_local_norm.z),Rho_label,GLToolkit::GetDialogFont());
		GLToolkit::DrawString((float)(O.x+factor*dir.x),(float)(O.y+factor*dir.y),(float)(O.z+factor*dir.z),dir_label,GLToolkit::GetDialogFont());*/
		GLToolkit::DrawStringRestore();
	}
}

void Region_full::SelectTrajPoint(int x,int y, size_t regionId) {

	int i;
	if(!isLoaded) return;

	// Select a vertex on a mouse click in 3D perspectivce view 
	// (x,y) are in screen coordinates
	// TODO: Handle clipped polygon

	// Check intersection of the facet and a "perspective ray"
	int *allXe = (int *)malloc(Points.size() * sizeof(int));
	int *allYe = (int *)malloc(Points.size() * sizeof(int));

	// Transform points to screen coordinates
	bool *ok = (bool *)malloc(Points.size()*sizeof(bool));
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
		if (mApp->trajectoryDetails && mApp->trajectoryDetails->GetRegionId()==regionId) mApp->trajectoryDetails->SelectPoint(selectedPoint);
	}

	free(allXe);
	free(allYe);
	free(ok);
	//UpdateSelection();
	
	if (mApp->regionInfo) mApp->regionInfo->Update();
}
/*
Region_full::Region_full(const Region_full &src) {
	//alfa0=src.alfa0;
	betax=src.params.betax_const_cm;
	betay=src.params.betay_const_cm;

	Bx_dir=Vector3d(src.params.Bx_dir.x,src.params.Bx_dir.y,src.params.Bx_dir.z);
	Bx_distr=new Distribution2D((int)src.params.nbDistr_MAG.x);
	*Bx_distr=*(src.Bx_distr);
	this->Bx_mode=src.params.Bx_mode;
	this->Bx_period=src.params.Bx_period;
	this->Bx_phase=src.params.Bx_phase;

	By_dir=Vector3d(src.params.By_dir.x,src.params.By_dir.y,src.params.By_dir.z);
	By_distr=new Distribution2D((int)src.params.nbDistr_MAG.y);
	*By_distr=*(src.By_distr);
	this->By_mode=src.params.By_mode;
	this->By_period=src.params.By_period;
	this->By_phase=src.params.By_phase;

	Bz_dir=Vector3d(src.params.Bz_dir.x,src.params.Bz_dir.y,src.params.Bz_dir.z);
	Bz_distr=new Distribution2D((int)src.params.nbDistr_MAG.z);
	*Bz_distr=*(src.Bz_distr);
	this->Bz_mode=src.params.Bz_mode;
	this->Bz_period=src.params.Bz_period;
	this->Bz_phase=src.params.Bz_phase;

	this->nbDistr_BXY=src.params.nbDistr_BXY;
	beta_x_distr=new Distribution2D(src.params.nbDistr_BXY);
	*beta_x_distr=*(src.beta_x_distr);
	beta_y_distr=new Distribution2D(src.params.nbDistr_BXY);
	*beta_y_distr=*(src.beta_y_distr);
	eta_distr=new Distribution2D(src.params.nbDistr_BXY);
	*eta_distr=*(src.params.eta_x_const_cm_distr);
	etaprime_distr=new Distribution2D(src.params.nbDistr_BXY);
	*etaprime_distr=*(src.params.eta_x_prime_const_distr);
	e_spread_distr=new Distribution2D(src.params.nbDistr_BXY);
	*e_spread_distr=*(src.e_spread_distr);

	this->B_const=Vector3d(src.params.B_const.x,src.params.B_const.y,src.params.B_const.z);
	this->dL=src.dL;
	this->E=src.E;
	this->current=src.params.current_mA;
	this->emittance=src.params.emittance_cm;
	this->coupling=src.params.coupling_percent;
	this->energy_spread=src.params.energy_spread_percent;
	this->enable_ort_polarization=src.params.enable_ort_polarization;
	this->enable_par_polarization=src.params.enable_par_polarization;
	this->energy_hi=src.params.energy_hi_eV;
	this->energy_low=src.params.energy_low_eV;
	this->eta=src.params.eta_x_const_cm;
	this->etaprime=src.params.eta_x_prime_const;
	this->gamma=src.params.gamma;
	//this->generation_mode=src.generation_mode;
	this->isLoaded=src.isLoaded;
	//this->i_struct1=src.i_struct1;
	this->limits=Vector3d(src.params.limits.x,src.params.limits.y,src.params.limits.z);
	//this->no_scans=src.no_scans;
	//this->nregions=src.nregions;
	this->particleMass=src.params.particleMass_GeV;
	this->Points=src.Points; //does it work like that?
	this->psimaxX=src.params.psimaxX_rad;
	this->psimaxY=src.params.psimaxY_rad;
	this->quad=src.params.quad; //write assignment operator...
	this->selectedPoint=src.selectedPoint;
	this->startDir=Vector3d(src.params.startDir.x,src.params.startDir.y,src.params.startDir.z);
	this->startPoint=Vector3d(src.params.startPoint.x,src.params.startPoint.y,src.params.startPoint.z);
	//this->teta0=src.teta0;
	//this->x0=src.x0;
	//this->y0=src.y0;
	//this->z0=src.z0;
	this->nbPointsToCopy=(int)src.params.nbPointsToCopy;
	this->nbDistr_MAG=Vector3d(src.params.nbDistr_MAG.x,src.params.nbDistr_MAG.y,src.params.nbDistr_MAG.z);
	this->beta_kind=src.params.beta_kind;
	this->AABBmin=Vector3d(src.AABBmin.x,src.AABBmin.y,src.AABBmin.z);
	this->AABBmax=Vector3d(src.AABBmax.x,src.AABBmax.y,src.AABBmax.z);
	this->fileName.assign(src.fileName);
	this->MAGXfileName.assign(src.MAGXfileName);
	this->MAGYfileName.assign(src.MAGYfileName);
	this->MAGZfileName.assign(src.MAGZfileName);
	this->BXYfileName.assign(src.BXYfileName);
}*/
/*
Region_full& Region_full::operator=(const Region_full &src) {
	//alfa0=src.alfa0;
	betax=src.params.betax_const_cm;
	betay=src.params.betay_const_cm;

	Bx_dir=Vector3d(src.params.Bx_dir.x,src.params.Bx_dir.y,src.params.Bx_dir.z);
	Bx_distr=new Distribution2D((int)src.params.nbDistr_MAG.x);
	*Bx_distr=*(src.Bx_distr);
	this->Bx_mode=src.params.Bx_mode;
	this->Bx_period=src.params.Bx_period;
	this->Bx_phase=src.params.Bx_phase;

	By_dir=Vector3d(src.params.By_dir.x,src.params.By_dir.y,src.params.By_dir.z);
	By_distr=new Distribution2D((int)src.params.nbDistr_MAG.y);
	*By_distr=*(src.By_distr);
	this->By_mode=src.params.By_mode;
	this->By_period=src.params.By_period;
	this->By_phase=src.params.By_phase;

	Bz_dir=Vector3d(src.params.Bz_dir.x,src.params.Bz_dir.y,src.params.Bz_dir.z);
	Bz_distr=new Distribution2D((int)src.params.nbDistr_MAG.z);
	*Bz_distr=*(src.Bz_distr);
	this->Bz_mode=src.params.Bz_mode;
	this->Bz_period=src.params.Bz_period;
	this->Bz_phase=src.params.Bz_phase;

	this->nbDistr_BXY=src.params.nbDistr_BXY;
	beta_x_distr=new Distribution2D(src.params.nbDistr_BXY);
	*beta_x_distr=*(src.beta_x_distr);
	beta_y_distr=new Distribution2D(src.params.nbDistr_BXY);
	*beta_y_distr=*(src.beta_y_distr);
	eta_distr=new Distribution2D(src.params.nbDistr_BXY);
	*eta_distr=*(src.params.eta_x_const_cm_distr);
	etaprime_distr=new Distribution2D(src.params.nbDistr_BXY);
	*etaprime_distr=*(src.params.eta_x_prime_const_distr);
	e_spread_distr=new Distribution2D(src.params.nbDistr_BXY);
	*e_spread_distr=*(src.e_spread_distr);

	this->B_const=Vector3d(src.params.B_const.x,src.params.B_const.y,src.params.B_const.z);
	this->dL=src.dL;
	this->E=src.E;
	this->current=src.params.current_mA;
	this->emittance=src.params.emittance_cm;
	this->coupling=src.params.coupling_percent;
	this->energy_spread=src.params.energy_spread_percent;
	this->enable_ort_polarization=src.params.enable_ort_polarization;
	this->enable_par_polarization=src.params.enable_par_polarization;
	this->energy_hi=src.params.energy_hi_eV;
	this->energy_low=src.params.energy_low_eV;
	this->eta=src.params.eta_x_const_cm;
	this->etaprime=src.params.eta_x_prime_const;
	this->gamma=src.params.gamma;
	//this->generation_mode=src.generation_mode;
	this->isLoaded=src.isLoaded;
	//this->i_struct1=src.i_struct1;
	this->limits=Vector3d(src.params.limits.x,src.params.limits.y,src.params.limits.z);
	//this->no_scans=src.no_scans;
	//this->nregions=src.nregions;
	this->particleMass=src.params.particleMass_GeV;
	//this->Points=src.Points; //trajectory points will be copied in the next step
	this->psimaxX=src.params.psimaxX_rad;
	this->psimaxY=src.params.psimaxY_rad;
	this->quad=src.params.quad; //write assignment operator...
	this->selectedPoint=src.selectedPoint;
	this->startDir=Vector3d(src.params.startDir.x,src.params.startDir.y,src.params.startDir.z);
	this->startPoint=Vector3d(src.params.startPoint.x,src.params.startPoint.y,src.params.startPoint.z);
	//this->teta0=src.teta0;
	//this->x0=src.x0;
	//this->y0=src.y0;
	//this->z0=src.z0;
	this->nbPointsToCopy=(int)src.params.nbPointsToCopy;
	this->nbDistr_MAG=Vector3d(src.params.nbDistr_MAG.x,src.params.nbDistr_MAG.y,src.params.nbDistr_MAG.z);
	this->beta_kind=src.params.beta_kind;
	this->AABBmin=Vector3d(src.AABBmin.x,src.AABBmin.y,src.AABBmin.z);
	this->AABBmax=Vector3d(src.AABBmax.x,src.AABBmax.y,src.AABBmax.z);
	this->fileName.assign(src.fileName);
	this->MAGXfileName.assign(src.MAGXfileName);
	this->MAGYfileName.assign(src.MAGYfileName);
	this->MAGZfileName.assign(src.MAGZfileName);
	this->BXYfileName.assign(src.BXYfileName);
	return *this;
}*/

int Region_full::LoadBXY(const std::string& fileName)
{
	std::fstream file(fileName);
	if (file.fail()) {
		std::ostringstream errMsg;
		errMsg << "Can't open the BXY file " << fileName;
		throw Error(errMsg.str().c_str());
	}
	int nbLines = 0;
	std::vector<std::vector<double>> distrValues;
	std::string line;
	while (std::getline(file, line)) {
		nbLines++;
		vector<string> tokens = SplitString(line);
		if (nbLines == 1) { //First line, format: "Coord type [N]"
			if ((tokens.size()==2 && Contains({ "0","1","2","3","l","x","y","z","L","X","Y","Z" },tokens[1]))
				|| (tokens.size() == 1 && Contains({ "0","1","2","3","l","x","y","z","L","X","Y","Z" }, tokens[0]))) { //Old format: "N coord_type"
				if (Contains({ "0","l","L" }, tokens.back())) {
					params.beta_kind = 0;
				} else if (Contains({ "1","x","X" }, tokens.back())) {
					params.beta_kind = 1;
				} else if (Contains({ "2","y","Y" }, tokens.back())) {
					params.beta_kind = 2;
				}
				else {
					params.beta_kind = 3;
				}
			}
			else  {
				throw Error("Couldn't parse first line of BXY file.\n"
					"Format should be: an integer between 0 and 3 (coord.type)\n"
					"Or X, Y, Z or L, lower- or uppercase\n"
					"Optionally preceded by the number of entries (old format, now ignored)");
			}
		}
		else { //Entries. Format: "Coord BetaX BetaY [EtaX EtaX' [(Energy-spread)||(AlphaX AlphaY)]]"
			std::vector<double> lineValues{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			if (Contains({ 3,5,6,7 },tokens.size())) {

				for (size_t i = 0; i < Min(tokens.size(), (size_t)7); i++) {
					std::istringstream tmp(tokens[i]);
					if (!(tmp >> lineValues[i])) {
						std::ostringstream errMsg;
						errMsg << "Can't parse token " << i << " in line " << nbLines << "of BXY file.\n";
						errMsg << "Token value: " << tokens[i];
						throw Error(errMsg.str().c_str());
					}
				}
				if (tokens.size() == 6) lineValues[5] = 0.0; //Old BXY file in [Coord BetaX BetaY EtaX EtaX' energy_spread] format -> drop energy_spread!

				//Sanitize user input
				if (lineValues[1] <= 0.0) {
					std::ostringstream errMsg;
					errMsg << "BetaX value in line " << nbLines << " of BXY file is not positive.\n";
					errMsg << "Parsed value: " << lineValues[1];
					throw Error(errMsg.str().c_str());
				}
				if (lineValues[2] <= 0.0) {
					std::ostringstream errMsg;
					errMsg << "BetaY value in line " << nbLines << " of BXY file is not positive.\n";
					errMsg << "Parsed value: " << lineValues[2];
					throw Error(errMsg.str().c_str());
				}
				/*
				if (lineValues[3] < 0.0) {
					std::ostringstream errMsg;
					errMsg << "EtaX value in line " << nbLines << " of BXY file is negative.\n";
					errMsg << "Parsed value: " << lineValues[3];
					throw Error(errMsg.str().c_str());
				}
				if (lineValues[4] < 0.0) {
					std::ostringstream errMsg;
					errMsg << "EtaX_prime value in line " << nbLines << " of BXY file is negative.\n";
					errMsg << "Parsed value: " << lineValues[4];
					throw Error(errMsg.str().c_str());
				}
				*/
				distrValues.push_back(lineValues);
			}
			else {
				if (line == "\n") {
					continue; //Empty line, probably at file end, ignore it
				}
				else {
					std::stringstream errMsg;
					errMsg << "Couldn't parse line " << nbLines << " of BXY file.\nExpecting 3, 5, 6 or 7 double values.";
					throw Error(errMsg.str().c_str());
				}
			}
		}
	}
	int nbDistr_BXY = (int)distrValues.size();

	betaFunctions.Clear();

	for (int i = 0; i < nbDistr_BXY; i++) {
		std::vector<double> newVec(distrValues[i].begin() + 1, distrValues[i].end()); //Subvector without the X coord
		betaFunctions.AddPair(distrValues[i][0], newVec);
	}
	return nbDistr_BXY;
}

void Region_full::SaveParam(FileWriter *file) {
	file->Write("param_file_version:");file->Write(PARAMVERSION,"\n");
	file->Write("startPos_cm:");file->Write(params.startPoint.x);file->Write(params.startPoint.y);file->Write(params.startPoint.z,"\n");
	file->Write("startDir_cm:");file->Write(params.startDir.x);file->Write(params.startDir.y);file->Write(params.startDir.z,"\n");
	file->Write("dL_cm:");file->Write(params.dL_cm,"\n");
	file->Write("boundaries_cm:");file->Write(params.limits.x);file->Write(params.limits.y);file->Write(params.limits.z,"\n");
	file->Write("particleMass_GeV:");file->Write(params.particleMass_GeV,"\n");
	file->Write("beamEnergy_GeV:");file->Write(params.E_GeV,"\n");
	file->Write("beamCurrent_mA:");file->Write(params.current_mA,"\n");
	file->Write("emittance_cm:");file->Write(params.emittance_cm,"\n");
	file->Write("energy_spread_percent:"); file->Write(params.energy_spread_percent, "\n");
	file->Write("coupling_percent:");file->Write(params.coupling_percent,"\n");
	if (params.emittance_cm!=0.0) { //non-ideal beam
		file->Write("const_beta_x_cm:");file->Write(params.betax_const_cm,"\n");
		if (params.betax_const_cm<0.0) {//use BXY file
			file->Write("BXYfileName:\"");file->Write(FileUtils::GetFilename(BXYfileName).c_str());file->Write("\"\n"); //truncate path
		} else {//constants
			file->Write("const_beta_y_cm:");file->Write(params.betay_const_cm,"\n");
			file->Write("const_eta_x_cm:");file->Write(params.eta_x_const_cm,"\n");
			file->Write("const_eta_x_prime:");file->Write(params.eta_x_prime_const,"\n");
		}
	}
	file->Write("E_min_eV:");file->Write(params.energy_low_eV,"\n");
	file->Write("E_max_eV:");file->Write(params.energy_hi_eV,"\n");
	file->Write("enable_par_polarization:");file->Write(params.enable_par_polarization,"\n");
	file->Write("enable_ort_polarization:");file->Write(params.enable_ort_polarization,"\n");
	file->Write("psiMax_X_Y_rad:");file->Write(params.psimaxX_rad);file->Write(params.psimaxY_rad,"\n");
	
	file->Write("Bx_mode:");file->Write(params.Bx_mode,"\n");
	if (params.Bx_mode==B_MODE_CONSTANT) {
		file->Write("Bx_const_Tesla:");file->Write(params.B_const.x,"\n");
	} else {
		file->Write("Bx_fileName:\"");file->Write(FileUtils::GetFilename(MAGXfileName).c_str());file->Write("\"\n"); //truncate path
	}
	file->Write("By_mode:");file->Write(params.By_mode,"\n");
	if (params.By_mode==B_MODE_CONSTANT) {
		file->Write("By_const_Tesla:");file->Write(params.B_const.y,"\n");
	} else {
		file->Write("By_fileName:\""); file->Write(FileUtils::GetFilename(MAGYfileName).c_str()); file->Write("\"\n"); //truncate path
	}
	file->Write("Bz_mode:");file->Write(params.Bz_mode,"\n");
	if (params.Bz_mode==B_MODE_CONSTANT) {
		file->Write("Bz_const_Tesla:");file->Write(params.B_const.z,"\n");
	} else {
		file->Write("Bz_fileName:\""); file->Write(FileUtils::GetFilename(MAGZfileName).c_str()); file->Write("\"\n"); //truncate path
	}
}

void Region_full::LoadParam(FileReader *file){

	//Creating references to X,Y,Z components (to avoid writing code 3 times)
	double* Bconst_ptr[3]={&params.B_const.x,&params.B_const.y,&params.B_const.z};
	Vector3d* Bdir_ptr[3]={&params.Bx_dir,&params.By_dir,&params.Bz_dir};
	int* Bmode_ptr[3]={&params.Bx_mode,&params.By_mode,&params.Bz_mode};
	double* Bperiod_ptr[3]={&params.Bx_period,&params.By_period,&params.Bz_period};
	double* Bphase_ptr[3]={&params.Bx_phase,&params.By_phase,&params.Bz_phase};
	Distribution2D* distr_ptr[3]={&Bx_distr,&By_distr,&Bz_distr};
	std::string* MagFileNamePtr[3]={&MAGXfileName,&MAGYfileName,&MAGZfileName};
	std::string Bmode_str[3]={"Bx_mode","By_mode","Bz_mode"};
	std::string Bconst_str[3]={"Bx_const_Tesla","By_const_Tesla","Bz_const_Tesla"};
	std::string MagFileNameStr[3]={"Bx_fileName","By_fileName","Bz_fileName"};
	int paramVersion;

	//prg->SetMessage("Reading parameter file...");

	file->ReadKeyword("param_file_version");file->ReadKeyword(":");paramVersion=file->ReadInt();
	if (paramVersion>PARAMVERSION) {
		throw Error("Param file version not supported");
	}
	
	file->ReadKeyword("startPos_cm");file->ReadKeyword(":");
	double x=file->ReadDouble();double y=file->ReadDouble();double z=file->ReadDouble(); params.startPoint=Vector3d(x,y,z);
	file->ReadKeyword("startDir_cm");file->ReadKeyword(":");
	x=file->ReadDouble(); y=file->ReadDouble(); z=file->ReadDouble(); params.startDir=Vector3d(x,y,z);
	file->ReadKeyword("dL_cm");file->ReadKeyword(":");params.dL_cm = file->ReadDouble();
	file->ReadKeyword("boundaries_cm");file->ReadKeyword(":");
	x=file->ReadDouble(); y=file->ReadDouble(); z=file->ReadDouble(); params.limits=Vector3d(x,y,z);
	file->ReadKeyword("particleMass_GeV");file->ReadKeyword(":");params.particleMass_GeV=file->ReadDouble();
	file->ReadKeyword("beamEnergy_GeV");file->ReadKeyword(":");params.E_GeV=file->ReadDouble();
	if (paramVersion>=2) {file->ReadKeyword("beamCurrent_mA");file->ReadKeyword(":");params.current_mA=file->ReadDouble();}
	else params.current_mA=1000;

	params.emittance_cm= params.eta_x_const_cm = params.eta_x_prime_const = params.energy_spread_percent =
		params.betax_const_cm = params.betay_const_cm = 0.0; //default values
	params.coupling_percent=100.0;
	params.nbDistr_BXY=0;

	(paramVersion>=3)?file->ReadKeyword("emittance_cm"):file->ReadKeyword("emittance");file->ReadKeyword(":");params.emittance_cm=file->ReadDouble();
	if (paramVersion >= 3) {
		file->ReadKeyword("energy_spread_percent");file->ReadKeyword(":");params.energy_spread_percent=file->ReadDouble();
	}
	if (paramVersion>=2) {
		(paramVersion>=3)?file->ReadKeyword("coupling_percent"):file->ReadKeyword("coupling");file->ReadKeyword(":");params.coupling_percent=file->ReadDouble();
	}
	if (params.emittance_cm!=0.0) { //non-ideal beam
		(paramVersion>=3)?file->ReadKeyword("const_beta_x_cm"):file->ReadKeyword("beta_x");file->ReadKeyword(":");params.betax_const_cm=file->ReadDouble();
		params.beta_kind=0;
		if (params.betax_const_cm<0.0) {//use BXY file
			file->ReadKeyword("BXYfileName");file->ReadKeyword(":");
				char tmp[2048];
				char tmp2[2048];
				char CWD [MAX_PATH];
				_getcwd( CWD, MAX_PATH );
				strcpy(tmp2,file->ReadString()); //get BXY file name
				sprintf(tmp,"%s\\tmp\\%s",CWD,tmp2); //look for file in tmp directory (extracted from syn7z)
				if (FileUtils::Exist(tmp)) //file found in tmp directory
					BXYfileName.assign(tmp);
				else {//not found, look for it in current directory (syn files)				
					std::string path=FileUtils::GetPath(file->GetName());
					sprintf(tmp,"%s%s",path.c_str(),tmp2);
					if (FileUtils::Exist(tmp)) //file found in current directory
						BXYfileName.assign(tmp);
					else {//not in tmp, nor in current, error.
						sprintf(tmp2,"Referenced BXY file not found:\n%s",tmp);
						throw Error(tmp2);
					}
				}
				try {
					params.nbDistr_BXY = LoadBXY(BXYfileName);
				} catch(Error &e) {
					//geom->Clear();
					sprintf(tmp, "Error loading BXY file (%s):\n%s", BXYfileName.c_str(), e.GetMsg());
					throw Error(tmp);
				}
			 
		} else { //constants
			(paramVersion>=3)?file->ReadKeyword("const_beta_y_cm"):file->ReadKeyword("beta_y");file->ReadKeyword(":");params.betay_const_cm=file->ReadDouble();
			(paramVersion>=3)?file->ReadKeyword("const_eta_x_cm"):file->ReadKeyword("eta");file->ReadKeyword(":");params.eta_x_const_cm=file->ReadDouble();
			(paramVersion>=3)?file->ReadKeyword("const_eta_x_prime"):file->ReadKeyword("eta_prime");file->ReadKeyword(":");params.eta_x_prime_const=file->ReadDouble();
			if (paramVersion==1) {file->ReadKeyword("coupling");file->ReadKeyword(":");params.coupling_percent=file->ReadDouble();}
			if (paramVersion < 3) {
				file->ReadKeyword("energy_spread");
				file->ReadKeyword(":");
				params.energy_spread_percent = file->ReadDouble();
			}
		}
	}
	file->ReadKeyword("E_min_eV");file->ReadKeyword(":");params.energy_low_eV=file->ReadDouble();
	file->ReadKeyword("E_max_eV");file->ReadKeyword(":");params.energy_hi_eV=file->ReadDouble();
	file->ReadKeyword("enable_par_polarization");file->ReadKeyword(":");params.enable_par_polarization=file->ReadInt();
	file->ReadKeyword("enable_ort_polarization");file->ReadKeyword(":");params.enable_ort_polarization=file->ReadInt();

	//Polarization component selection
	if (params.enable_par_polarization && params.enable_ort_polarization)
		params.polarizationCompIndex = 0; //Full polarization
	else if (params.enable_par_polarization)
		params.polarizationCompIndex = 1; //Parallel polarization
	else
		params.polarizationCompIndex = 2; //Orthogonal polarization

	(paramVersion>=3)?file->ReadKeyword("psiMax_X_Y_rad"):file->ReadKeyword("psiMax_rad");file->ReadKeyword(":");params.psimaxX_rad=file->ReadDouble();params.psimaxY_rad=file->ReadDouble();

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
					std::string path=FileUtils::GetPath(file->GetName());
					sprintf(tmp,"%s%s",path.c_str(),tmp2);
					if (FileUtils::Exist(tmp)) //file found in current directory
						MagFileNamePtr[i]->assign(tmp);
					else {//not in tmp, nor in current, error.
						sprintf(tmp2,"Referenced MAG file not found: %s",tmp);
						SAFE_DELETE(file);
						throw Error(tmp2);
					}
				}
				FileReader MAGfile((char*)MagFileNamePtr[i]->c_str());
				*(distr_ptr[i])=LoadMAGFile(&MAGfile,Bdir_ptr[i],Bperiod_ptr[i],Bphase_ptr[i],*Bmode_ptr[i]);
		}
	}
	params.gamma=abs(params.E_GeV/ params.particleMass_GeV);

	//prg->SetMessage("Calculating trajectory...");
	CalculateTrajectory(1000000); //max 1 million points
	isLoaded=true;
	
	if (mApp->regionInfo) mApp->regionInfo->Update();
}
