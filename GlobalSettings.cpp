/*
File:        GlobalSettings.cpp
Description: Global settings dialog (antiAliasing,units)
Program:     SynRad

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*/

#include "GLApp/GLLabel.h"
#include "GlobalSettings.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLInputBox.h"
#include "Utils.h"
#include "Synrad.h"

extern SynRad *mApp;

static const int   plWidth[] = {15,40,70,70,50,330};
static const char *plName[] = {"#","PID","Mem Usage","Mem Peak","CPU","Status"};
static const int   plAligns[] = { ALIGN_LEFT,ALIGN_CENTER,ALIGN_CENTER,ALIGN_CENTER,ALIGN_CENTER,ALIGN_LEFT };


//HANDLE synradHandle;

// --------------------------------------------------------------------

GlobalSettings::GlobalSettings():GLWindow() {

	int wD = 610;
	int hD = 525;

	SetTitle("Global Settings");
	SetIconfiable(TRUE);

	GLTitledPanel *panel = new GLTitledPanel("View settings");
	panel->SetBounds(5,2,300,78);
	Add(panel);

	chkAntiAliasing = new GLToggle(0,"Anti-Aliasing");
	chkAntiAliasing->SetBounds(10,20,80,19);
	panel->Add(chkAntiAliasing);

	chkWhiteBg = new GLToggle(0,"White Background");
	chkWhiteBg->SetBounds(10,50,80,19);
	panel->Add(chkWhiteBg);

	GLTitledPanel *panel2 = new GLTitledPanel("Low flux mode");
	panel2->SetBounds(310,2,300,78);
	Add(panel2);

	lowFluxInfo = new GLButton(0, "Info");
	lowFluxInfo->SetBounds(480, 20, 40, 19);
	panel2->Add(lowFluxInfo);

	lowFluxToggle = new GLToggle(0,"Enable low flux mode");
	lowFluxToggle->SetBounds(315, 20, 150, 19);
	panel2->Add(lowFluxToggle);

	GLLabel *cutoffLabel = new GLLabel("Cutoff ratio:");
	cutoffLabel->SetBounds(330, 50, 80, 19);
	panel2->Add(cutoffLabel);

	cutoffText = new GLTextField(0, "");
	cutoffText->SetBounds(410, 50, 70, 19);
	cutoffText->SetEditable(FALSE);
	panel2->Add(cutoffText);

	GLTitledPanel *panel4 = new GLTitledPanel("Program settings");
	panel4->SetBounds(5,85,600,90);
	Add(panel4);

	GLLabel *asLabel = new GLLabel("Autosave frequency (minutes):");
	asLabel->SetBounds(15,100,160,19);
	panel4->Add(asLabel);

	autoSaveText = new GLTextField(0,"");
	autoSaveText->SetBounds(170,100,30,19);
	panel4->Add(autoSaveText);

	chkSimuOnly = new GLToggle(0,"Autosave only when simulation is running");
	chkSimuOnly->SetBounds(10,125,160,19);
	panel4->Add(chkSimuOnly);

	chkCheckForUpdates = new GLToggle(0,"Check for updates at startup");
	chkCheckForUpdates->SetBounds(315,100,100,19);
	//chkCheckForUpdates->SetEnabled(FALSE);
	Add(chkCheckForUpdates);

	chkAutoUpdateFormulas = new GLToggle(0, "Auto update formulas");
	chkAutoUpdateFormulas->SetBounds(315, 125, 160, 19);
	Add(chkAutoUpdateFormulas);

	chkCompressSavedFiles = new GLToggle(0,"Compress saved files (use .SYN7Z format)");
	chkCompressSavedFiles->SetBounds(10,150,100,19);
	Add(chkCompressSavedFiles);

	/*chkNonIsothermal = new GLToggle(0,"Non-isothermal system (textures only, experimental)");
	chkNonIsothermal->SetBounds(315,125,100,19);
	Add(chkNonIsothermal);*/

	GLTitledPanel *panel3 = new GLTitledPanel("Subprocess control");
	panel3->SetBounds(5,180,wD-10,hD-225);
	Add(panel3);


	processList = new GLList(0);
	processList->SetHScrollVisible(TRUE);
	processList->SetSize(6,MAX_PROCESS);
	processList->SetColumnWidths((int*)plWidth);
	processList->SetColumnLabels((char **)plName);
	processList->SetColumnAligns((int *)plAligns);
	processList->SetColumnLabelVisible(TRUE);
	processList->SetBounds(10,195,wD-20,hD-305);
	panel3->Add(processList);

	char tmp[128];
	sprintf(tmp,"Number of CPU cores:     %d",mApp->numCPU);
	GLLabel *coreLabel = new GLLabel(tmp);
	coreLabel->SetBounds(10,hD-99,120,19);
	panel3->Add(coreLabel);

	GLLabel *l1 = new GLLabel("Number of subprocesses:");
	l1->SetBounds(10,hD-74,120,19);
	panel3->Add(l1);

	nbProcText = new GLTextField(0,"");
	nbProcText->SetEditable(TRUE);
	nbProcText->SetBounds(135,hD-76,30,19);
	panel3->Add(nbProcText);

	restartButton = new GLButton(0,"Restart");
	restartButton->SetBounds(170,hD-76,90,19);
	panel3->Add(restartButton);

	maxButton = new GLButton(0,"Change MAX generated photons");
	maxButton->SetBounds(wD-195,hD-76,180,19);
	panel3->Add(maxButton);

	// ---------------------------------------------------

	applyButton = new GLButton(0,"Apply");
	applyButton->SetBounds(wD-175,hD-43,80,19);
	Add(applyButton);

	cancelButton = new GLButton(0,"Dismiss");
	cancelButton->SetBounds(wD-90,hD-43,80,19);
	Add(cancelButton);

	// Center dialog
	int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;
	SetBounds(xD,yD,wD,hD);

	RestoreDeviceObjects();

	lastUpdate = 0.0f;
	for(int i=0;i<MAX_PROCESS;i++) lastCPUTime[i]=-1.0f;
	memset(lastCPULoad,0,MAX_PROCESS*sizeof(float));
}


// --------------------------------------------------------------------

void GlobalSettings::Display(Worker *w) {
	worker = w;
	char tmp[256];
	chkAntiAliasing->SetState(mApp->antiAliasing);
	chkWhiteBg->SetState(mApp->whiteBg);

	cutoffText->SetText(worker->lowFluxCutoff);
	cutoffText->SetEditable(worker->lowFluxMode);
	lowFluxToggle->SetState(worker->lowFluxMode);

	sprintf(tmp,"%g",mApp->autoSaveFrequency);
	autoSaveText->SetText(tmp);
	chkSimuOnly->SetState(mApp->autoSaveSimuOnly);
	chkCheckForUpdates->SetState(mApp->checkForUpdates);
	chkCompressSavedFiles->SetState(mApp->compressSavedFiles);
	
	int nb = worker->GetProcNumber();
	sprintf(tmp,"%d",nb);
	nbProcText->SetText(tmp);
	
	
	SetVisible(TRUE);

}

// ----------------------------------------------------------------

void GlobalSettings::SMPUpdate(float appTime) {

	if(!IsVisible() || IsIconic()) return;  
	int nb = worker->GetProcNumber();

	if( appTime-lastUpdate>1.0 && nb>0 ) {

		char tmp[512];
		PROCESS_INFO pInfo;
		int  states[MAX_PROCESS];
		char statusStr[MAX_PROCESS][64];

		memset(states,0,MAX_PROCESS*sizeof(int));
		worker->GetProcStatus(states,(char **)statusStr);

		processList->ResetValues();
		for(int i=0;i<nb;i++) {
			DWORD pid = worker->GetPID(i);
			sprintf(tmp,"%d",i+1);
			processList->SetValueAt(0,i,tmp);
			sprintf(tmp,"%d",pid);
			processList->SetValueAt(1,i,tmp);
			if( !GetProcInfo(pid,&pInfo) ) {
				processList->SetValueAt(2,i,"0 KB");
				processList->SetValueAt(3,i,"0 KB");
				processList->SetValueAt(4,i,"0 %");
				processList->SetValueAt(5,i,"Dead");
			} else {
				processList->SetValueAt(2,i,FormatMemory(pInfo.mem_use));
				processList->SetValueAt(3,i,FormatMemory(pInfo.mem_peak));

				// CPU usage
				if( lastCPUTime[i]!=-1.0f ) {
					float dTime = appTime-lastUpdate;
					float dCPUTime = (float)pInfo.cpu_time-lastCPUTime[i];
					float cpuLoad = dCPUTime/dTime;
					lastCPULoad[i] = 0.85f*cpuLoad + 0.15f*lastCPULoad[i];
					int percent = (int)(100.0f*lastCPULoad[i] + 0.5f);
					if(percent<0) percent=0;
					sprintf(tmp,"%d %%",percent);
					processList->SetValueAt(4,i,tmp);
				} else {
					processList->SetValueAt(4,i,"---");
				}
				lastCPUTime[i] = (float)pInfo.cpu_time;

				// State/Status
				char status[128];
				_snprintf(tmp,127,"%s: %s",prStates[states[i]],statusStr[i]);
				status[127]=0;
				processList->SetValueAt(5,i,tmp);

			}
		}

		lastUpdate = appTime;
	}

}
// ----------------------------------------------------------------

void GlobalSettings::RestartProc() {

	int nbProc;
	if( sscanf(nbProcText->GetText(),"%d",&nbProc)==0 ) {
		GLMessageBox::Display("Invalid process number","Error",GLDLG_OK,GLDLG_ICONERROR);
	} else {
		//char tmp[128];
		//sprintf(tmp,"Kill all running sub-process(es) and start %d new ones ?",nbProc);
		//int rep = GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);

		if( mApp->AskToReset() ) {
			if(nbProc<=0 || nbProc>MAX_PROCESS) {
				GLMessageBox::Display("Invalid process number [1..32]","Error",GLDLG_OK,GLDLG_ICONERROR);
			} else {
				try {
					worker->SetProcNumber(nbProc);
					worker->Reload();
					mApp->SaveConfig();
				} catch(Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
				}
			}
		}
	}

}
// --------------------------------------------------------------------

void GlobalSettings::ProcessMessage(GLComponent *src,int message) {

	switch(message) {
	case MSG_BUTTON:

		if(src==cancelButton) {

			GLWindow::ProcessMessage(NULL,MSG_CLOSE);

		} else if (src==restartButton) {
			RestartProc();
		} else if (src==maxButton) {
			if( worker->GetGeometry()->IsLoaded() ) {
				char tmp[128];
				sprintf(tmp,"%I64d",worker->maxDesorption);
				char *val = GLInputBox::GetInput(tmp,"Desorption max (0=>endless)","Edit MAX");
				if( val ) {
					llong maxDes;
					if( sscanf(val,"%I64d",&maxDes)==0 ) {
						GLMessageBox::Display("Invalid 'maximum desorption' number","Error",GLDLG_OK,GLDLG_ICONERROR);
					} else {
						worker->SetMaxDesorption(maxDes);
					}
				}
			} else {
				GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			}
		} else if (src==applyButton) {
			mApp->antiAliasing = chkAntiAliasing->GetState();
			mApp->whiteBg = chkWhiteBg->GetState();
			mApp->checkForUpdates = chkCheckForUpdates->GetState();
			mApp->compressSavedFiles = chkCompressSavedFiles->GetState();
			mApp->autoUpdateFormulas = chkAutoUpdateFormulas->GetState();
			mApp->autoSaveSimuOnly=chkSimuOnly->GetState();
			double autosavefreq;
			if (!autoSaveText->GetNumber(&autosavefreq) || !(autosavefreq>0.0)) {
				GLMessageBox::Display("Invalid autosave frequency","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
			mApp->autoSaveFrequency = autosavefreq;
			
			double cutoffnumber;
			if (!cutoffText->GetNumber(&cutoffnumber) || !(cutoffnumber>0.0 && cutoffnumber<1.0)) {
				GLMessageBox::Display("Invalid cutoff ratio, must be between 0 and 1","Error",GLDLG_OK,GLDLG_ICONWARNING);
				return;
			}

			if (abs(worker->lowFluxCutoff - cutoffnumber)>1e-10 || worker->lowFluxMode!=lowFluxToggle->GetState()) {
				if (mApp->AskToReset()) {
					worker->lowFluxCutoff = cutoffnumber;
					worker->lowFluxMode = lowFluxToggle->GetState();
					mApp->changedSinceSave = TRUE;
					// Send to sub process
					try { worker->Reload(); }
					catch (Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
						return;
					}
					mApp->UpdateFacetHits();
				}
				
			}

			

			GLWindow::ProcessMessage(NULL,MSG_CLOSE); 
			return;

		}
		else if (src == lowFluxInfo) {
			GLMessageBox::Display("Low flux mode helps to gain more statistics on low flux/power parts of the system, at the expense\n"
				"of higher flux/power parts. If a traced photon reflects from a low reflection probability surface, regardless of that probability,\n"
				"a reflected test photon representing a reduced flux will still be traced. Therefore test photons can reach low flux areas more easily, but\n"
				"at the same time tracng a test photon takes longer. The cutoff ratio defines what ratio of the originally generated flux/power\n"
				"can be neglected. If, for example, it is 0.001, then, when after subsequent reflections the test photon carries less than 0.1%\n"
				"of the original flux, it will be eliminated. A good advice is that if you'd like to see flux across N orders of magnitude, set it to 1E-N"
				, "Low flux mode", GLDLG_OK, GLDLG_ICONINFO);
			return;
		}
		break;

	case MSG_TEXT:
		ProcessMessage(applyButton,MSG_BUTTON);
		break;
	

	case MSG_TOGGLE:
		if (src == lowFluxToggle) {
			cutoffText->SetEditable(lowFluxToggle->GetState());
		}
		break;
	}
	GLWindow::ProcessMessage(src,message);
}
