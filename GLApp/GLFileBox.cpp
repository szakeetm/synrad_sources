/*
  File:        GLFileBox.cpp
  Description: File selection box class (SDL/OpenGL OpenGL application framework)
  Author:      J-L PONS (2007)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#include "GLWindow.h"
#include "GLFileBox.h"
#include "GLMessageBox.h"
#include "GLToolkit.h"

#include <math.h>
#include <string>
#include <vector>

#ifdef WIN
#include <direct.h>
#include <io.h>
#include <windows.h>
#include <shlobj.h>
//for windows open/save dialogs
#include <Windows.h>
#include <CommDlg.h>
#endif

#define MAX_HIST   256
#define MAX_DRIVE  64

#define PREV_BTN   1
#define OK_BTN     2
#define CANCEL_BTN 3
#define PATH_TEXT  4
#define FILE_TEXT  5

static char *cNames[] = { "Name" , "Size" , "Modified" };
static int  cWidths[] = { 220,89,115 };
static int  cAligns[] = { ALIGN_LEFT,ALIGN_RIGHT,ALIGN_RIGHT };

// Returned structure
static FILENAME retFile;

// Path history
char        *pathHistory[MAX_HIST];
int          nbPath = 0;

// File name history
char        *fileHistory[MAX_HIST];
int          nbFile = 0;

#ifdef WIN
char        *drivePaths[MAX_DRIVE];
char        *driveNames[MAX_DRIVE];
int          nbDrive = 0;
#endif

using namespace std;

// --------------------------------------------------------------

GLFileBox::GLFileBox(char *path,char *fileName,char *title,char **filters,int nbFilter,BOOL openMode):GLWindow() {

  char tmp[256];
  int wD = 450;
  int hD = 350;
  mode = openMode;

  SetTitle(title);

  // Path
  strcpy(curPath,path);
  nbPath = 0;

  // Filters
  filterLabel = new GLLabel("Filters");
  filterLabel->SetBounds(5,307,45,19);
  Add(filterLabel);

  this->nbFilter = nbFilter;
  curFilter = 0;
  filterCombo = new GLCombo(0);
  if( nbFilter==0 ) {
    filterCombo->SetSize(1);
    filterCombo->SetValueAt(0,"All files (*.*)");
  } else {
    filterCombo->SetSize(nbFilter);
    for(int i=0;i<nbFilter;i++) {
      filterName[i] = _strdup(filters[2*i+0]);
      filterExt[i]  = _strdup(filters[2*i+1]);
      sprintf(tmp,"%s (%s)",filterName[i],filterExt[i]);
      filterCombo->SetValueAt(i,tmp); 
    }
  }
  filterCombo->SetSelectedIndex(curFilter);
  filterCombo->SetBounds(55,307,wD-150,19);
  Add(filterCombo);

  // File name combo
  fileLabel = new GLLabel("File name");
  fileLabel->SetBounds(5,282,45,19);
  Add(fileLabel);

  fileText = new GLCombo(FILE_TEXT);
  fileText->SetEditable(TRUE);
  fileText->SetBounds(55,282,wD-150,19);
  fileText->SetSize(nbFile);
  for(int i=0;i<nbFile;i++) {
    fileText->SetValueAt(i,fileHistory[i]);
  }
  if(fileName) fileText->SetSelectedValue(fileName);
  Add(fileText);

  // Path name combo
  pathLabel = new GLLabel("Directoy");
  pathLabel->SetBounds(5,5,40,19);
  Add(pathLabel);

  pathText = new GLCombo(PATH_TEXT);
  AddToPathHist(path);
  pathText->SetSelectedIndex(0);
  pathText->SetEditable(TRUE);
  pathText->SetBounds(50,5,wD-150,19);
  Add(pathText);

  // Prev button
  prevButton = new GLButton(PREV_BTN,"Back");
  prevButton->SetBounds(wD-90,5,85,19);
  Add(prevButton);

  // File list
  fileList = new GLList(0);
  UpdateFileList(path);
  fileList->SetBounds(5,29,wD-10,245);
  Add(fileList);

  if( openMode )
    okButton = new GLButton(OK_BTN,"Open");
  else
    okButton = new GLButton(OK_BTN,"Save");
  okButton->SetBounds(wD-90,282,85,19);
  okButton->SetEnabled(!openMode);
  Add(okButton);

  cancelButton = new GLButton(CANCEL_BTN,"Cancel");
  cancelButton->SetBounds(wD-90,307,85,19);
  Add(cancelButton);

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  RestoreDeviceObjects();
  rCode = CANCEL_BTN;

}

// --------------------------------------------------------------

GLFileBox::~GLFileBox() {

  for(int i=0;i<nbFilter;i++) {
    SAFE_FREE(filterName[i]);
    SAFE_FREE(filterExt[i]);
  }

}

// --------------------------------------------------------------

void GLFileBox::ProcessMessage(GLComponent *src,int message) {

  char tmp[512];

  switch( message ) {

    // Buttons --------------------------------------------
    case MSG_BUTTON:
      if( src==cancelButton ) {
        GLWindow::ProcessMessage(NULL,MSG_CLOSE);
        return;
      }
      if( src==okButton ) {
        strcpy(curFile,fileText->GetSelectedValue());
        strcpy(curPath,pathText->GetSelectedValue());
        rCode = OK_BTN;
        GLWindow::ProcessMessage(NULL,MSG_CLOSE);
        return;
      }
      if( src==prevButton ) {
        Back();
        UpdateFileList(curPath);
      }
      break;

    // Texts ----------------------------------------------
    case MSG_TEXT: {
      if( src->GetId() == PATH_TEXT ) {
        // Return pressed in nameText
        char *nPath = pathText->GetSelectedValue();
        if( CheckDirectory(nPath) ) {
          UpdateFileList(nPath);
          AddToPathHist(curPath);
        } else {
          pathText->SetSelectedValue(curPath);
        }
      }
    } break;

    // Combos ---------------------------------------------
    case MSG_COMBO: {
      if( src==filterCombo ) {
        int s = filterCombo->GetSelectedIndex();
        if(s>=0) {
          curFilter = s;
          UpdateFileList(curPath);
        }
      } else if ( src==pathText ) {
        char *nPath = pathText->GetSelectedValue();
        UpdateFileList(nPath);
      }
    } break;

    // Lists ---------------------------------------------
    case MSG_LIST_DBL: {
        int row = fileList->GetSelectedRow();
        char *fName = fileList->GetValueAt(0,row);
        if( strncmp(fName,":B:",3)==0 ) {
          // Change directory
          if( strcmp(curPath,"My Computer")==0 ) {
            UpdateFileList(drivePaths[row]);
          } else if( strcmp(fName+3,"..")==0 ) {
            Back();
            UpdateFileList(curPath);
            AddToPathHist(curPath);
          } else {
            sprintf(tmp,"%s\\%s",curPath,fName+3);
            UpdateFileList(tmp);
            AddToPathHist(tmp);
          }
        } else {
          // Double click on File
          strcpy(curFile,fName);
          rCode = OK_BTN;
          GLWindow::ProcessMessage(NULL,MSG_CLOSE);
        }
      }break;
    case MSG_LIST: {
        int row = fileList->GetSelectedRow();
        if( row>=0 ) {
          char *fName = fileList->GetValueAt(0,row);
          BOOL isDir = (strncmp(fName,":B:",3)==0);
          if( isDir ) {
            if(mode) {
              fileText->SetSelectedValue("");
              okButton->SetEnabled(FALSE);
            }
          } else {
            strcpy(curFile,fName);
            fileText->SetSelectedValue(curFile);
            okButton->SetEnabled(TRUE);
          }
        }
      }break;
  }
  GLWindow::ProcessMessage(src,message);
}

// --------------------------------------------------------------

BOOL GLFileBox::CheckDirectory(char *dirName) {

#ifdef WIN

  if(strcmp(dirName,":B:My Computer")==0)
    return TRUE;

  struct _finddata_t seqfile;
  intptr_t f;
  char errMsg[512];

  f=_findfirst( dirName , &seqfile );
  if( f==-1L ) {
    sprintf(errMsg,"%s\nDirectory not found",dirName);
    GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
    return FALSE;
  }

  if((seqfile.attrib & _A_SUBDIR)==0) {
     sprintf(errMsg,"%s\nNot a directory",dirName);
     GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
     _findclose(f);
     return FALSE;
  }

  _findclose(f);

#endif

  return TRUE;

}

// --------------------------------------------------------------

FILENAME *GLFileBox::OpenFile(char *path,char *fileName,char *title,const char *filters,int nbFilter) {
	

	FILENAME *ret = NULL;

	//Windows File Open dialog
	OPENFILENAME ofn;       // common dialog box structure
	char fileNm[260];       // buffer for file name
				
	// Initialize OPENFILENAME
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.lpstrFile=fileNm;
	ofn.lpstrFile[0] = '\0';
	ofn.nMaxFile = sizeof(fileNm);
	ofn.lpstrFilter = filters;
	ofn.nFilterIndex = 1;
	ofn.lpstrFileTitle = NULL;
	ofn.nMaxFileTitle = 0;
	ofn.lpstrTitle = title;
	ofn.lpstrInitialDir = NULL;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

	// Display the Open dialog box. 

	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	//_chdir( CWD ); 				
	if (GetOpenFileName(&ofn)==TRUE) {
		_chdir( CWD ); //OpenFileName dialog changes the working driectory, so we must set it back
		strcpy(retFile.fullName,fileNm);
		ret = &retFile;
	}

	return ret;
	/* Original cross-platform OpenGl OpenFile dialog. Replaced with Windows file dialog
  if(!title) title = "Open File";

#ifdef WIN
  if(!path || strcmp(path,".")==0) path=_getcwd(NULL,0);
  GLFileBox::InitDrivePaths();
#endif

  FILENAME *ret = NULL;

  GLFileBox *f = new GLFileBox(path,fileName,title,filters,nbFilter,TRUE);
  f->DoModal();

  if( f->rCode == OK_BTN ) {
    strcpy(retFile.path , f->curPath);
    strcpy(retFile.file , f->curFile);
    sprintf(retFile.fullName,"%s\\%s",f->curPath,f->curFile);
    ret = &retFile;
    f->AddToFileHist(f->curFile);
  }

  delete f;
  return ret;
  */

}

std::vector<FILENAME> GLFileBox::OpenMultipleFiles(const char *filters, const char *title) {

	std::vector<FILENAME> ret;

	//Windows File Open dialog
	OPENFILENAME ofn;       // common dialog box structure
	char fileNm[65536];     // buffer for file name, increased size for multiple files

	// Initialize OPENFILENAME
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.lpstrFile = fileNm;
	ofn.lpstrFile[0] = '\0';
	ofn.nMaxFile = sizeof(fileNm);
	ofn.lpstrFilter = filters;
	ofn.lpstrTitle = title;
	ofn.nFilterIndex = 1;
	ofn.lpstrFileTitle = NULL;
	ofn.nMaxFileTitle = 0;
	ofn.lpstrInitialDir = NULL;
	ofn.Flags = OFN_ALLOWMULTISELECT | OFN_EXPLORER | OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

	// Display the Open dialog box. 
	char CWD[MAX_PATH];
	_getcwd(CWD, MAX_PATH);
	//_chdir(CWD);
	if (GetOpenFileName(&ofn) == TRUE) {
		_chdir(CWD); //OpenFileName dialog changes the working driectory, so we must set it back
		
		//Parsing from http://stackoverflow.com/questions/26317556/how-to-get-list-of-selected-files-when-using-getopenfilename-with-multiselect

		char* str = ofn.lpstrFile;
		std::string directory = str;
		str += (directory.length() + 1);
		BOOL moreThanOne = FALSE;
		while (*str) {
			moreThanOne = TRUE;
			std::string filename = str;
			str += (filename.length() + 1);
			// use the filename, e.g. add it to a vector
			FILENAME retFile;
			strcpy(retFile.path, directory.c_str());
			strcpy(retFile.file, filename.c_str());
			sprintf(retFile.fullName, "%s\\%s", directory.c_str(), filename.c_str());
			ret.push_back(retFile);
		}
		if (!moreThanOne) {
			FILENAME retFile;
			strcpy(retFile.file, directory.c_str());
			strcpy(retFile.fullName, directory.c_str());
			ret.push_back(retFile);
		}
	}

	return ret;
}

// --------------------------------------------------------------

FILENAME *GLFileBox::SaveFile(char *path,char *fileName,char *title,const char *filters,int nbFilter) {
		FILENAME *ret = NULL;

	//Windows File Open dialog
	OPENFILENAME ofn;       // common dialog box structure
	char fileNm[260];       // buffer for file name

	// Initialize OPENFILENAME
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.lpstrFile=fileNm;
	ofn.lpstrFile[0] = '\0';
	ofn.nMaxFile = sizeof(fileNm);
	ofn.lpstrFilter = filters;
	ofn.nFilterIndex = 1;
	ofn.lpstrFileTitle = NULL;
	ofn.nMaxFileTitle = 0;
	ofn.lpstrInitialDir = NULL;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

	// Display the Open dialog box. 

	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	_chdir( CWD ); 				
	if (GetSaveFileName(&ofn)==TRUE) {
		_chdir( CWD ); //GetSaveFileName dialog changes the working driectory, so we must set it back
		strcpy(retFile.fullName,fileNm);
		ret = &retFile;
	}

	return ret;
  /* Old cross-platform OpenGl SaveFile dialog, replaced by windows dialog
  if(!title) title = "Save File";

#ifdef WIN
  if(!path || strcmp(path,".")==0) path=_getcwd(NULL,0);
  GLFileBox::InitDrivePaths();
#endif

  FILENAME *ret = NULL;

  GLFileBox *f = new GLFileBox(path,fileName,title,filters,nbFilter,FALSE);
  f->DoModal();

  if( f->rCode == OK_BTN ) {
    strcpy(retFile.path , f->curPath);
    strcpy(retFile.file , f->curFile);
    sprintf(retFile.fullName,"%s\\%s",f->curPath,f->curFile);
    ret = &retFile;
    f->AddToFileHist(f->curFile);
  }

  delete f;
  return ret;
  */
}

// --------------------------------------------------------------

void GLFileBox::InitDrivePaths() {
#ifdef WIN

  TCHAR szPath[MAX_PATH];
  TCHAR sDrive[3];
  sDrive[2]=0;

  if( nbDrive==0 ) {

    // Special folders ----------------------------

    if(SUCCEEDED(SHGetFolderPath(NULL, 
                 CSIDL_DESKTOP|CSIDL_FLAG_CREATE, 
                 NULL, 
                 0, 
                 szPath))) 
    {
      drivePaths[nbDrive] = _strdup(szPath);
      driveNames[nbDrive] = "Desktop";
      nbDrive++;
    }  

    if(SUCCEEDED(SHGetFolderPath(NULL, 
                 CSIDL_PERSONAL|CSIDL_FLAG_CREATE, 
                 NULL, 
                 0, 
                 szPath))) 
    {
      drivePaths[nbDrive] = _strdup(szPath);
      driveNames[nbDrive] = "My Documents";
      nbDrive++;
    }  

    // Drives ----------------------------

    TCHAR buf[256];
    ::GetLogicalDriveStrings(256, buf);

    for (LPCTSTR drive=buf; *drive; drive += strlen(drive)+1) {
        
      sDrive[0] = drive[0];
      sDrive[1] = drive[1];
      switch( ::GetDriveType(drive) ) {
        case DRIVE_REMOVABLE:
          sprintf(szPath,"Removable Drive (%s)",sDrive);
          break;  
        case DRIVE_FIXED:
          sprintf(szPath,"Local Disk (%s)",sDrive);
          break; 
        case DRIVE_REMOTE:
          sprintf(szPath,"Remote Drive (%s)",sDrive);
          break; 
        case DRIVE_CDROM:
          sprintf(szPath,"CD Drive (%s)",sDrive);
          break;  
        case DRIVE_RAMDISK:
          sprintf(szPath,"RAM Disk (%s)",sDrive);
          break;  
        default:
          sprintf(szPath,"Unknown Drive (%s)",sDrive);
          break; 
      }

      drivePaths[nbDrive] = _strdup(sDrive);
      driveNames[nbDrive] = _strdup(szPath);
      nbDrive++;

    }

  }
#endif
}

// --------------------------------------------------------------

void GLFileBox::AddToFileHist(char *file) {

  // Remove first item if full
  if(nbFile>=MAX_HIST) {
    SAFE_FREE(fileHistory[0]);
    for(int i=1;i<nbFile;i++)
      fileHistory[i-1] = fileHistory[i];
    nbFile--;
  }
  //Search insertion pos
  BOOL found = FALSE;
  int i=0,cmp=-1;
  while(!found && i<nbFile) {
    cmp = _stricmp(fileHistory[i],file);
    found = (cmp>=0);
    if( !found ) i++;
  }

  // Already added
  if(cmp==0) return;

  if( !found ) {
    // Add at the end
    fileHistory[nbFile]=_strdup(file);
  } else {
    for(int j=nbFile;j>i;j--)
      fileHistory[j] = fileHistory[j-1];
    fileHistory[i]=_strdup(file);
  }
  nbFile++;

}


// --------------------------------------------------------------

void GLFileBox::AddToPathHist(char *path) {

  if( strcmp(path,"My Computer")==0 ) return;

  // Remove first item if full
  if(nbPath>=MAX_HIST) {
    SAFE_FREE(pathHistory[0]);
    for(int i=1;i<nbPath;i++)
      pathHistory[i-1] = pathHistory[i];
    nbPath--;
  }

  //Search insertion pos
  BOOL found = FALSE;
  int i=0,cmp=-1;
  while(!found && i<nbPath) {
    cmp = _stricmp(pathHistory[i],path);
    found = (cmp>=0);
    if( !found ) i++;
  }

  // Already added
  if(cmp==0) return;

  if( !found ) {
    // Add at the end
    pathHistory[nbPath]=_strdup(path);
  } else {
    for(int j=nbPath;j>i;j--)
      pathHistory[j] = pathHistory[j-1];
    pathHistory[i]=_strdup(path);
  }
  nbPath++;

  // Update the combo
#ifdef WIN
  pathText->SetSize(nbPath+1);
  pathText->SetValueAt(0,":B:My Computer");
  for(int i=1;i<=nbPath;i++) {
    pathText->SetValueAt(i,pathHistory[i-1]);
  }
#else
  pathText->SetSize(nbPath);
  for(int i=0;i<nbPath;i++) {
    pathText->SetValueAt(i,pathHistory[i]);
  }
#endif
}

// --------------------------------------------------------------

BOOL GLFileBox::MatchFilters(char *fileName) {

  if( nbFilter==0 ) return TRUE;
  if( strcmp(filterExt[curFilter],"*.*")==0 ) return TRUE;

  char *ext  = GetExtension(fileName);
  char *extF = GetExtension(filterExt[curFilter]);

  return _stricmp(ext,extF)==0;

}

// --------------------------------------------------------------

char *GLFileBox::GetTimeStr(time_t time) {

  static char ret[32];
  if( time>0 ) {
    struct tm *ts = localtime((time_t *)&time);
    sprintf(ret,"%02d/%02d/%04d %02d:%02d:%02d",ts->tm_mday,ts->tm_mon+1,ts->tm_year+1900,
                                                ts->tm_hour,ts->tm_min,ts->tm_sec);
  } else {
    strcpy(ret,"");
  }

  return ret;

}

// --------------------------------------------------------------

char *GLFileBox::GetSizeStr(DWORD size) {

  static char ret[256];
  const char *suffixStr[] = {"KB","MB","GB","TB"};
  double dSize = (double)size;
  int suffix = 0;

  while( dSize >= 1024.0 ) {
    dSize /= 1024.0;
    suffix++;
  }

  if( suffix==0 ) {
    sprintf(ret,"%u bytes",size);
  } else {
    if( fabs( dSize - floor(dSize) )<1e-3 )
      sprintf(ret,"%.0f%s",dSize,suffixStr[suffix-1]);
    else
      sprintf(ret,"%.2f%s",dSize,suffixStr[suffix-1]);
  }
  return ret;

}

// --------------------------------------------------------------

char *GLFileBox::GetExtension(char *fileName) {

  char *p = strrchr(fileName,'.');
  if( p ) return p+1;
  return "";

}

// --------------------------------------------------------------

char *GLFileBox::GetName(char *fileName) {
  if( strncmp(fileName,":B:",3)==0 )
    return fileName+3;
  else
    return fileName;
}

// --------------------------------------------------------------

void GLFileBox::Back() {

  if(strlen(curPath)==2)
    if(curPath[1]==':') {
      strcpy(curPath,":B:My Computer");
      return;
    }

  char  *p = strrchr(curPath,'\\');
  if(!p) p = strrchr(curPath,'/');
  if(p)  *p=0;

}

// --------------------------------------------------------------

void GLFileBox::UpdateFileList(char *path) {

  vector<string> files;
  int nbFile = 0;

#ifdef WIN

  if( strcmp(path,":B:My Computer")==0 ) {

    for(int i=0;i<nbDrive;i++) {
      files.push_back(string(":B:") + driveNames[i]);
      files.push_back("");
      files.push_back("");
      nbFile++;
    }
    strcpy(curPath,"My Computer");

  } else {

    struct   _finddata_t seqfile;
    intptr_t fId;
    string  filter = string(path) + "\\*.*";

    fId = _findfirst(filter.c_str(), &seqfile);
    if( fId != -1L ) {

      do {

        if( (strcmp(seqfile.name,".")!=0) ) {

          if( (seqfile.attrib & _A_SUBDIR) ) {
            files.push_back(string(":B:") + seqfile.name);
            files.push_back( string("") );
            files.push_back( string(GetTimeStr(seqfile.time_access)) );
            nbFile++;
          } else {
            if( MatchFilters(seqfile.name) ) {
              files.push_back(seqfile.name);
              files.push_back(string(GetSizeStr(seqfile.size)));
              files.push_back( string(GetTimeStr(seqfile.time_access)) );
              nbFile++;
            }
          }

        }

      } while (_findnext(fId, &seqfile) != -1L);

      _findclose(fId);

    }
    strcpy(curPath,path);

  }

#endif

  fileList->SetSize(3,nbFile);
  for(int i=0;i<nbFile;i++) {
    fileList->SetValueAt(0,i,(char *)files[3*i+0].c_str());
    fileList->SetValueAt(1,i,(char *)files[3*i+1].c_str());
    fileList->SetValueAt(2,i,(char *)files[3*i+2].c_str());
  }
  fileList->SetColumnAligns(cAligns);
  fileList->SetColumnWidths(cWidths);
  fileList->SetColumnLabels(cNames);
  fileList->SetColumnLabelVisible(TRUE);

  pathText->SetSelectedValue(curPath);
  pathText->ScrollTextToEnd();

}