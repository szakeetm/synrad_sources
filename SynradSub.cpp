/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/
#define NOMINMAX
#include <windows.h>

#include <tlhelp32.h>

//#include <iostream>

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Buffer_shared.h"

#include "Simulation.h"
#ifdef WIN
#include <Process.h> // For _getpid()
#endif

// Global process variables

#define WAITTIME    100  // Answer in STOP mode

static Dataport *dpControl=NULL;
static Dataport *dpHit=NULL;
static Dataport *dpLog = NULL;
static int       prIdx;
static size_t    prState;
static size_t    prParam;
static size_t    prParam2;
static DWORD     hostProcessId;
static char      ctrlDpName[32];
static char      loadDpName[32];
static char      hitsDpName[32];
static char		 logDpName[32];

bool end = false;
bool IsProcessRunning(DWORD pid);

void GetState() {
  prState = PROCESS_READY;
  prParam = 0;

  if( AccessDataport(dpControl) ) {
    SHCONTROL *master = (SHCONTROL *)dpControl->buff;
    prState = master->states[prIdx];
    prParam = master->cmdParam[prIdx];
    prParam2 = master->cmdParam2[prIdx];
    master->cmdParam[prIdx] = 0;
    master->cmdParam2[prIdx] = 0;

    ReleaseDataport(dpControl);

	if (!IsProcessRunning(hostProcessId)) {
		printf("Host synrad.exe (process id %d) not running. Closing.",hostProcessId);
		SetErrorSub("Host synrad.exe not running. Closing subprocess.");
		end = true;
	}

  } else {
	  printf("Subprocess couldn't connect to Synrad.\n");
	  SetErrorSub("No connection to main program. Closing subprocess.");
	  Sleep(5000);
	  end = true;
  }
}

size_t GetLocalState() {
  return prState;
}

void SetState(size_t state,const char *status, bool changeState, bool changeStatus) {

	prState = state;
	printf("\n setstate %zd \n",state);
	if( AccessDataport(dpControl) ) {
		SHCONTROL *master = (SHCONTROL *)dpControl->buff;
		if (changeState) master->states[prIdx] = state;
		if (changeStatus) {
			strncpy(master->statusStr[prIdx], status, 127);
			master->statusStr[prIdx][127] = 0;
		}
		ReleaseDataport(dpControl);
	}
}

void SetStatus(char *status) {
	if (AccessDataport(dpControl)) {
		SHCONTROL *master = (SHCONTROL *)dpControl->buff;
		strncpy(master->statusStr[prIdx], status, 127);
		master->statusStr[prIdx][127] = 0;
		ReleaseDataport(dpControl);
	}
}

void SetErrorSub(const char *message) {

  printf("Error: %s\n",message);
  SetState(PROCESS_ERROR,message);

}

char *GetSimuStatus() {

  static char ret[128];
  size_t count = sHandle->totalDesorbed;
  size_t max   = sHandle->ontheflyParams.desorptionLimit / sHandle->ontheflyParams.nbProcess;

      if( max!=0 ) {
        double percent = (double)(count)*100.0 / (double)(max);
        sprintf(ret,"(%s) MC %I64d/%I64d (%.1f%%)",sHandle->name,count,max,percent);
      } else {
        sprintf(ret,"(%s) MC %I64d",sHandle->name,count);
      }

  return ret;

}

void SetReady() {

  if(sHandle->loadOK)
    SetState(PROCESS_READY,GetSimuStatus());
  else
    SetState(PROCESS_READY,"(No geometry)");

}

void Load() {

  Dataport *loader;
  size_t hSize;

  // Load geometry
  loader = OpenDataport(loadDpName,prParam);
  if( !loader ) {
    char err[512];
    sprintf(err,"Failed to connect to 'loader' dataport %s (%zd Bytes)",loadDpName, prParam);
    SetErrorSub(err);
    return;
  }
  
  printf("Connected to %s\n",loadDpName);

  if( !LoadSimulation(loader) ) {
    CLOSEDP(loader);
    return;
  }
  CLOSEDP(loader);

  //Connect to log dataport
  if (sHandle->ontheflyParams.enableLogging) {
	  dpLog = OpenDataport(logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
	  if (!dpLog) {
		  char err[512];
		  sprintf(err, "Failed to connect to 'dpLog' dataport %s (%zd Bytes)", logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
		  SetErrorSub(err);
		  sHandle->loadOK = false;
		  return;
	  }
	  //*((size_t*)dpLog->buff) = 0; //Autofill with 0. Besides, we don't write without access!
  }


  // Connect to hit dataport
  hSize = GetHitsSize();
  dpHit = OpenDataport(hitsDpName,hSize);
  if( !dpHit ) {
    SetErrorSub("Failed to connect to 'hits' dataport");
    return;
  }

  printf("Connected to %s (%zd bytes)\n",hitsDpName,hSize);

}

bool UpdateParams() {

	Dataport *loader;

	// Load geometry
	loader = OpenDataport(loadDpName, prParam);
	if (!loader) {
		char err[512];
		sprintf(err, "Failed to connect to 'loader' dataport %s (%zd Bytes)", loadDpName, prParam);
		SetErrorSub(err);
		return false;
	}

	printf("Connected to %s\n", loadDpName);

	bool result = UpdateOntheflySimuParams(loader);
	CLOSEDP(loader);

	if (sHandle->ontheflyParams.enableLogging) {
		dpLog = OpenDataport(logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
		if (!dpLog) {
			char err[512];
			sprintf(err, "Failed to connect to 'dpLog' dataport %s (%zd Bytes)", logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
			SetErrorSub(err);
			return false;
		}
		//*((size_t*)dpLog->buff) = 0; //Autofill with 0, besides we would need access first
	}
	sHandle->tmpParticleLog.clear();
	sHandle->tmpParticleLog.shrink_to_fit();
	if (sHandle->ontheflyParams.enableLogging) sHandle->tmpParticleLog.reserve(sHandle->ontheflyParams.logLimit / sHandle->ontheflyParams.nbProcess);

	return result;
}

int main(int argc,char* argv[])
{
  bool eos = false;

  if(argc!=3) {
    printf("Usage: synradSub peerId index\n");
    return 1;
  }
  
  hostProcessId=atoi(argv[1]);
  prIdx = atoi(argv[2]);

  sprintf(ctrlDpName,"SNRDCTRL%s",argv[1]);
  sprintf(loadDpName,"SNRDLOAD%s",argv[1]);
  sprintf(hitsDpName,"SNRDHITS%s",argv[1]);
  sprintf(logDpName, "SNRDLOG%s", argv[1]);

  dpControl = OpenDataport(ctrlDpName,sizeof(SHCONTROL));
  if( !dpControl ) {
    printf("Usage: Cannot connect to SRDCTRL%s\n",argv[1]);
    return 1;
  }

  printf("Connected to %s (%zd bytes), synradSub.exe #%d\n",ctrlDpName,sizeof(SHCONTROL),prIdx);

  InitSimulation();

  // Sub process ready
  SetReady();

  // Main loop
  while( !end ) {

	GetState();
    switch(prState) {

      case COMMAND_LOAD:
        printf("COMMAND: LOAD (%zd,%I64d)\n",prParam,prParam2);
        Load();
        if( sHandle->loadOK ) {
          //sHandle->desorptionLimit = prParam2; // 0 for endless
          SetReady();
        }
        break;

	  case COMMAND_UPDATEPARAMS:
		  printf("COMMAND: UPDATEPARAMS (%zd,%I64d)\n", prParam, prParam2);
		  if (UpdateParams()) {
			  SetState(prParam, GetSimuStatus());
		  }
		  break;

	  case COMMAND_RELEASEDPLOG:
		  printf("COMMAND: UPDATEPARAMS (%zd,%I64d)\n", prParam, prParam2);
		  CLOSEDP(dpLog);
		  SetState(prParam, GetSimuStatus());
		  break;


      case COMMAND_START:
        printf("COMMAND: START (%zd,%I64d)\n",prParam,prParam2);
        if( sHandle->loadOK ) {
          if( StartSimulation() )
            SetState(PROCESS_RUN,GetSimuStatus());
          else {
            if( GetLocalState()!=PROCESS_ERROR )
              SetState(PROCESS_DONE,GetSimuStatus());
          }
        } else
          SetErrorSub("No geometry loaded");
        break;

      case COMMAND_PAUSE:
        printf("COMMAND: PAUSE (%zd,%I64d)\n",prParam,prParam2);
        if( !sHandle->lastHitUpdateOK ) {
          // Last update not successful, retry with a longer tomeout
			if (dpHit && (GetLocalState() != PROCESS_ERROR)) UpdateHits(dpHit,dpLog,prIdx,60000);
        }
        SetReady();
        break;

      case COMMAND_RESET:
        printf("COMMAND: RESET (%zd,%I64d)\n",prParam,prParam2);
        ResetSimulation();
        SetReady();
        break;

      case COMMAND_EXIT:
        printf("COMMAND: EXIT (%zd,%I64d)\n",prParam,prParam2);
        end = true;
        break;

      case COMMAND_CLOSE:
        printf("COMMAND: CLOSE (%zd,%I64d)\n",prParam,prParam2);
        ClearSimulation();
        CLOSEDP(dpHit);
		CLOSEDP(dpLog);
        SetReady();
        break;

      case PROCESS_RUN:
        SetStatus(GetSimuStatus()); //update hits only
        eos = SimulationRun();      // Run during 1 sec
        if(dpHit && (GetLocalState()!=PROCESS_ERROR)) UpdateHits(dpHit,dpLog,prIdx,20); // Update hit with 20ms timeout. If fails, probably an other subprocess is updating, so we'll keep calculating and try it later (latest when the simulation is stopped).
        if(eos) {
          if( GetLocalState()!=PROCESS_ERROR ) {
            // Max desorption reached
            SetState(PROCESS_DONE,GetSimuStatus());
            printf("COMMAND: PROCESS_DONE (Max reached)\n");
          }
        }
        break;

      default:
        Sleep(WAITTIME);
        break;
    }

  }

  // Release
  SetState(PROCESS_KILLED,"");
  //CLOSEDP(dpControl);
  //CLOSEDP(dpHit);
  //Why bother closing dataports? Windows will release handles automatically.
  return 0;

}

bool IsProcessRunning(DWORD pid)
{
	HANDLE process = OpenProcess(SYNCHRONIZE, false, pid);
	DWORD ret = WaitForSingleObject(process, 0);
	CloseHandle(process);
	return ret == WAIT_TIMEOUT;
}