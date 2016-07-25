/*
  File:        GLWindowManager.cpp
  Description: SDL/OpenGL OpenGL application framework (Window management)
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

#include "GLWindowManager.h"
#include "GLMessageBox.h"
#include "GLApp.h"
#include "GLToolkit.h"
#include "GLApp.h"
#include "GLWindow.h"
#include "../SynRad.h"

extern GLApplication *theApp;
extern float      m_fTime;             // Current time in seconds

// Definitition
typedef struct {
  GLComponent *src;
  int          keyCode;
  int          modifier;
  int          id;
} ACC;

// Global variables
#define MAX_WIN 64
#define MAX_ACC 512

static GLWindow *allWin[MAX_WIN];
static int nbWindow = 0;
static ACC allAcc[MAX_ACC];
static int nbAcc = 0;
static int modState = NO_MODIFIER;
static int keyFocus = 0;

extern GLApplication *theApp;

// ---------------------------------------------------------------

GLWindow* GLWindowManager::GetWindow(int idx) {
  return allWin[idx];
}

// ---------------------------------------------------------------

int GLWindowManager::GetNbWindow() {
  return nbWindow;
}

// ---------------------------------------------------------------

void GLWindowManager::RegisterWindow(GLWindow *wnd) {

  if(nbWindow<MAX_WIN) {
    allWin[nbWindow] = wnd;
    nbWindow++;
  }

}

// -------------------------------------------

void GLWindowManager::BringToFront(GLWindow *wnd) {

  BOOL found = FALSE;
  int i=0;

  while(!found && i<nbWindow) {
    found = (allWin[i] == wnd);
    if(!found) i++;
  }
  if( found && i>0 ) {
    // Shit
    GLWindow *winFoc = allWin[i];
    for(int j=i;j<nbWindow-1;j++) allWin[j] = allWin[j+1];
    allWin[nbWindow-1] = winFoc;
  }

}

// -------------------------------------------

GLWindow *GLWindowManager::GetTopLevelWindow() {
  if( nbWindow ) return allWin[0];
  else           return NULL;
}

// ---------------------------------------------------------------

void GLWindowManager::UnRegisterWindow(GLWindow *wnd) {

  BOOL found = FALSE;
  int i=0;

  while(!found && i<nbWindow) {
    found = (allWin[i] == wnd);
    if(!found) i++;
  }
  if( found ) {
    // Remove
    for(int j=i;j<nbWindow-1;j++)
      allWin[j] = allWin[j+1];
    nbWindow--;
  }

}

// ---------------------------------------------------------------

void GLWindowManager::FullRepaint() {

  // Done 2 times for the 2 buffer
  Repaint();
  Repaint();

}

// ---------------------------------------------------------------

void GLWindowManager::DrawStats() {

#ifdef _DEBUG
	
  // Statistics
  if( theApp ) {
    int x,y,w,h;
    GLWindowManager::NoClip();
    allWin[0]->GetBounds(&x,&y,&w,&h);
    GLToolkit::GetDialogFont()->SetTextColor(1.0f,1.0f,0.0f);
    glDisable(GL_CULL_FACE);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glEnable(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    char polys[256];
    sprintf(polys,"Quad:%d Line:%d Paint:%.2f ms   ",
    theApp->nbPoly,theApp->nbLine,theApp->fPaintTime*1000.0);
    GLToolkit::GetDialogFont()->DrawTextFast(7,h-99,polys);
    sprintf(polys,"Restore:%d FrameMove:%.2f ms   ",theApp->nbRestore,theApp->fMoveTime*1000.0);
    GLToolkit::GetDialogFont()->DrawTextFast(7,h-83,polys);
    GLToolkit::GetDialogFont()->DrawTextFast(7,h-65,theApp->m_strEventStats);
    GLToolkit::GetDialogFont()->DrawTextFast(7,h-48,theApp->m_strFrameStats);
	sprintf(polys,"m_fTime:%4.2f",m_fTime);
	GLToolkit::GetDialogFont()->DrawTextFast(150,h-48,polys);
  }
  
#endif

}

// ---------------------------------------------------------------

void GLWindowManager::SetDefault() {

  // Default OpenGL settings
  glDisable(GL_CULL_FACE);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_POLYGON_OFFSET_LINE);
  glDisable(GL_POLYGON_OFFSET_FILL);
  glDisable(GL_POLYGON_OFFSET_POINT);
  glDisable(GL_POLYGON_SMOOTH);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_FASTEST);
  GLToolkit::GetDialogFont()->SetTextColor(0.0f,0.0f,0.0f);

}

// ---------------------------------------------------------------

void GLWindowManager::NoClip() {

  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glOrtho(0,wS,hS,0,-1,1);
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
  GLToolkit::SetViewport(0,0,wS,hS);

}

// ---------------------------------------------------------------

void GLWindowManager::AnimateIconify(GLWindow *src) {

  int iPos = allWin[0]->GetWidth();
  for(int i=1;i<nbWindow;i++)
    if( allWin[i]->IsIconic() )
      iPos -= (allWin[i]->GetIconWidth() + 3);

  int xw,yw,ww,hw;
  src->GetBounds(&xw,&yw,&ww,&hw);
  int iconWidth = src->GetIconWidth();

  double xs,ys,ws;
  double xe,ye,we;

  xs  = (double)xw;
  ys  = (double)yw;
  ws  = (double)ww;
  xe  = (double)iPos;
  ye  = 1.0;
  we  = (double)iconWidth;

  double t0=(double)SDL_GetTicks()/1000.0;
  double t1 = t0;

  double tA = 0.15;

  while((t1-t0)<tA) {
    double t = (t1-t0) / tA;
    int x = (int)(xs + t*(xe-xs) + 0.5);
    int y = (int)(ys + t*(ye-ys) + 0.5);
    int w = (int)(ws + t*(we-ws) + 0.5);
    for(int i=0;i<nbWindow;i++) {
      if( src==allWin[i] ) {
        src->SetPosition(x,y);
        src->PaintTitle(w,20);
      } else {
        allWin[i]->Paint();
      }
    }
    DrawStats();
    SDL_GL_SwapBuffers();
    t1 = (double)SDL_GetTicks()/1000.0;
  }

}

// ---------------------------------------------------------------

void GLWindowManager::AnimateMaximize(GLWindow *src,int fsX,int fsY,int fsWidth,int fsHeight) {

  int xw,yw,ww,hw;
  src->GetBounds(&xw,&yw,&ww,&hw);

  double xs,ys,ws,hs;
  double xe,ye,we,he;

  xs  = (double)xw;
  ys  = (double)yw;
  ws  = (double)ww;
  hs  = (double)hw;
  xe  = (double)fsX;
  ye  = (double)fsY;
  we  = (double)fsWidth;
  he  = (double)fsHeight;

  double t0=(double)SDL_GetTicks()/1000.0;
  double t1 = t0;

  double tA = 0.15;

  while((t1-t0)<tA) {
    double t = (t1-t0) / tA;
    int x = (int)(xs + t*(xe-xs) + 0.5);
    int y = (int)(ys + t*(ye-ys) + 0.5);
    int w = (int)(ws + t*(we-ws) + 0.5);
    int h = (int)(hs + t*(he-hs) + 0.5);
    src->SetBounds(x,y,w,h);
    Repaint();
    t1 = (double)SDL_GetTicks()/1000.0;
  }
  src->SetBounds(fsX,fsY,fsWidth,fsHeight);
  Repaint();

}

// ---------------------------------------------------------------

void GLWindowManager::AnimateFocus(GLWindow *src) {

  int xw,yw,ww,hw;
  src->GetBounds(&xw,&yw,&ww,&hw);
  int rw,gw,bw;
  src->GetBackgroundColor(&rw,&gw,&bw);

  double rs,gs,bs,ws;
  double re,ge,be,we;

  rs  = 255.0;
  gs  = 60.0;
  bs  = 60.0;
  ws  = 4.0;

  re  = (double)rw;
  ge  = (double)gw;
  be  = (double)bw;
  we  = 1.0;

  double t0=(double)SDL_GetTicks()/1000.0;
  double t1 = t0;

  double tA = 0.2;

  while((t1-t0)<tA) {

    RepaintNoSwap();
    double t = (t1-t0) / tA;
    double r = (rs + t*(re-rs));
    double g = (gs + t*(ge-gs));
    double b = (bs + t*(be-bs));
    float  w = (float)(ws + t*(we-ws));
    SATURATE(r,0.0,255.0);
    SATURATE(g,0.0,255.0);
    SATURATE(b,0.0,255.0);
    SATURATE(w,1.0,3.0);

    // Paint a highlighted red border
    NoClip();
    glDisable(GL_CULL_FACE);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glColor3d(r/255.0,g/255.0,b/255.0);
    glLineWidth(w);
    glBegin(GL_LINE_LOOP);
    _glVertex2i(xw   ,yw);
    _glVertex2i(xw+ww,yw);
    _glVertex2i(xw+ww,yw+hw);
    _glVertex2i(xw   ,yw+hw);
    glEnd();
    DrawStats();
    SDL_GL_SwapBuffers();

    t1 = (double)SDL_GetTicks()/1000.0;
  }

  // Return to normal state
  Repaint();

}

// ---------------------------------------------------------------

void GLWindowManager::AnimateDeIconify(GLWindow *src) {

  int iPos = allWin[0]->GetWidth();
  int iIdx = 0;
  BOOL found = FALSE;
  while(!found && iIdx<nbWindow) {
    found = (src==allWin[iIdx]);
    if( allWin[iIdx]->IsIconic() )
      iPos -= (allWin[iIdx]->GetIconWidth() + 3);
    if(!found) iIdx++;
  }
  if(!found) return;

  int xw,yw,ww,hw;
  src->GetBounds(&xw,&yw,&ww,&hw);
  int iconWidth = src->GetIconWidth();

  double xs,ys,ws;
  double xe,ye,we;

  xe  = (double)xw;
  ye  = (double)yw;
  we  = (double)ww;
  xs  = (double)iPos;
  ys  = 1.0;
  ws  = (double)iconWidth;

  double t0=(double)SDL_GetTicks()/1000.0;
  double t1 = t0;

  double tA = 0.15;

  while((t1-t0)<tA) {
    double t = (t1-t0) / tA;
    int x = (int)(xs + t*(xe-xs) + 0.5);
    int y = (int)(ys + t*(ye-ys) + 0.5);
    int w = (int)(ws + t*(we-ws) + 0.5);
    allWin[0]->Paint();
    iPos = allWin[0]->GetWidth();
    for(int i=1;i<nbWindow;i++) {
      if( src==allWin[i] ) {
        iPos -= (int)((allWin[i]->GetIconWidth() + 3) * (1.0-t));
        src->SetPosition(x,y);
        if(src->IsVisible()) src->PaintTitle(w,20);
      } else {
        if( allWin[i]->IsIconic() ) {
          iPos -= (allWin[i]->GetIconWidth() + 3);
          allWin[i]->SetPosition(iPos,1);
          allWin[i]->Paint();
        } else {
          allWin[i]->Paint();
        }
      }
    }
    DrawStats();
    SDL_GL_SwapBuffers();
    t1 = (double)SDL_GetTicks()/1000.0;
  }

}

// ---------------------------------------------------------------

void GLWindowManager::Resize() {

  for(int i=1;i<nbWindow;i++) 
    allWin[i]->UpdateOnResize();
  SynRad *mApp = (SynRad *)theApp;
  try {
		  mApp->worker.Update(0.0f);
  } catch(Error &e) {
	  GLMessageBox::Display((char *)e.GetMsg(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
  }
}

// ---------------------------------------------------------------
void  GLWindowManager::Repaint() {
  RepaintNoSwap();
  DrawStats();
  SDL_GL_SwapBuffers();

}

// ---------------------------------------------------------------

void GLWindowManager::RepaintNoSwap() {

#ifdef _DEBUG
  double t0 = theApp->GetTick();
#endif

  // Search maximized window
  BOOL found = FALSE;
  int maxIdx = nbWindow-1;
  while(maxIdx>=0 && !found) {
    found = allWin[maxIdx]->IsMaximized() && !allWin[maxIdx]->IsIconic();
    if(!found) maxIdx--;
  }

  // Icon position
  int iPos = allWin[0]->GetWidth();

  SetDefault();
  if( maxIdx>0 ) {
    allWin[0]->PaintMenuBar();
  } else {
    allWin[0]->Paint();
  }
  GLToolkit::CheckGLErrors("GLWindowManager::RepaintNoSwap()");


  for(int i=1;i<nbWindow;i++) {
    if( allWin[i]->IsIconic() ) {
      iPos -= (allWin[i]->GetIconWidth() + 3);
      allWin[i]->SetPosition(iPos,1);
      allWin[i]->Paint();
    } else {
      if(i>=maxIdx) allWin[i]->Paint();
    }
  }
  for(int i=0;i<nbWindow;i++) allWin[i]->PaintMenu();

#ifdef _DEBUG
  theApp->fPaintTime = 0.9*theApp->fPaintTime + 0.1*(theApp->GetTick() - t0);
#endif

}

// ---------------------------------------------------------------

void GLWindowManager::RepaintRange(int w0,int w1) {

#ifdef _DEBUG
  double t0 = theApp->GetTick();
#endif

   if (!(w0<64 && w1<64)) throw Error("Buffer overrun: GLWindowManager::RepaintRange, array allWin");
  SetDefault();
  for(int i=w0;i<w1;i++) allWin[i]->Paint();  
  for(int i=w0;i<w1;i++) allWin[i]->PaintMenu();

#ifdef _DEBUG
  theApp->fPaintTime = 0.9*theApp->fPaintTime + 0.1*(theApp->GetTick() - t0);
#endif

  DrawStats();
  SDL_GL_SwapBuffers();

}

// ---------------------------------------------------------------

void GLWindowManager::RestoreDeviceObjects() {
  for(int i=0;i<nbWindow;i++) allWin[i]->RestoreDeviceObjects();
}

// ---------------------------------------------------------------

void GLWindowManager::InvalidateDeviceObjects() {
  for(int i=0;i<nbWindow;i++) allWin[i]->InvalidateDeviceObjects();
}

// ---------------------------------------------------------------

BOOL GLWindowManager::IsCtrlDown() {
  return (modState & CTRL_MODIFIER)!=0;
}

// ---------------------------------------------------------------

BOOL GLWindowManager::IsShiftDown() {
  return (modState & SHIFT_MODIFIER)!=0;
}

// ---------------------------------------------------------------

BOOL GLWindowManager::IsAltDown() {
  return (modState & ALT_MODIFIER)!=0;
}

BOOL GLWindowManager::IsSpaceDown() {
  return (modState & SPACE_MODIFIER)!=0;
}
// ---------------------------------------------------------------

BOOL GLWindowManager::IsCapsLockOn() {
 return (modState & CAPSLOCK_MODIFIER)!=0;
}

// ---------------------------------------------------------------

BOOL GLWindowManager::ProcessKey(SDL_Event *evt,BOOL processAcc) {

  BOOL accFound = FALSE;

  // Handle key modifier
  int unicode = (evt->key.keysym.unicode & 0x7F);
  if( !unicode ) unicode = evt->key.keysym.sym;

  if( evt->type == SDL_KEYDOWN )
  {
    if( unicode == SDLK_LCTRL ||  unicode == SDLK_RCTRL )
      modState |= CTRL_MODIFIER;

    if( unicode == SDLK_LSHIFT ||  unicode == SDLK_RSHIFT )
      modState |= SHIFT_MODIFIER;

    if( unicode == SDLK_RALT || unicode == SDLK_LALT )
      modState |= ALT_MODIFIER;
	
	if( unicode == SDLK_CAPSLOCK)
      modState |= CAPSLOCK_MODIFIER;

	if( unicode == SDLK_SPACE)
      modState |= SPACE_MODIFIER;
  }

  if( evt->type == SDL_KEYUP )
  {
    int altMask   =  ALT_MODIFIER; altMask   = ~altMask;
    int ctrlMask  = CTRL_MODIFIER; ctrlMask  = ~ctrlMask;
    int shiftMask = SHIFT_MODIFIER;shiftMask = ~shiftMask;
	int capsLockMask = CAPSLOCK_MODIFIER;capsLockMask = ~capsLockMask;
	int spaceMask = SPACE_MODIFIER;spaceMask = ~spaceMask;

    if( unicode == SDLK_LCTRL ||  unicode == SDLK_RCTRL )
      modState &= ctrlMask;

    if( unicode == SDLK_LSHIFT ||  unicode == SDLK_RSHIFT )
      modState &= shiftMask;

    if( unicode == SDLK_RALT || unicode == SDLK_LALT )
      modState &= altMask;

	if( unicode == SDLK_CAPSLOCK )
      modState &= capsLockMask;
  
	if (unicode == SDLK_SPACE)
		modState &= spaceMask;
  }

  // Process
  if( evt->type == SDL_KEYDOWN ) {
    // Search acc (Use keysym)
    int i = 0;
    while(i<nbAcc && !accFound) {
      accFound = (allAcc[i].keyCode == evt->key.keysym.sym) && ( (allAcc[i].modifier & modState)!=0 );
      if( !accFound ) i++;
    }
    if( accFound ) {
      GLContainer *p = allAcc[i].src->GetParent();
      if(p) p->SetFocus(allAcc[i].src);
      allAcc[i].src->ProcessAcc(allAcc[i].id);
    }
  }

  return accFound;

}

// ---------------------------------------------------------------

void GLWindowManager::RegisterAcc(GLComponent *src,int keyCode,int modifier,int accId) {

  // Search if already exists
  BOOL accFound = FALSE;
  int i = 0;
  while(i<nbAcc && !accFound) {
    accFound = (allAcc[i].keyCode == keyCode) && (allAcc[i].modifier == modifier);
    if( !accFound ) i++;
  }

  if( accFound ) {
    // Just replace id
    allAcc[i].id = accId;
  } else if( nbAcc<MAX_ACC ) {
    // Add new one
    allAcc[nbAcc].src = src;
    allAcc[nbAcc].keyCode = keyCode;
    allAcc[nbAcc].modifier = modifier;
    allAcc[nbAcc].id = accId;
    nbAcc++;
  }

}

// ---------------------------------------------------------------

BOOL GLWindowManager::ManageEvent(SDL_Event *evt) {

  if( evt->type == SDL_KEYDOWN || evt->type == SDL_KEYUP ) {
    if( !ProcessKey(evt,TRUE) ) {
      // Process key event
      if(nbWindow) allWin[keyFocus]->ManageEvent(evt);
      return (keyFocus==0);
    } else {
      // Give keyboard focus to main on accelerator
      keyFocus=0;
      return TRUE;
    }
  }

  // Searh for dragged window
  BOOL draggFound = FALSE;
  int i = 0;
  while(i<nbWindow && !draggFound) {
    draggFound = allWin[i]->IsDragging();
    if(!draggFound) i++;
  }

  if( draggFound ) {

    // Priority on dragged window
    allWin[i]->ManageEvent(evt);

  } else {

    // Process all non modal window
    // Create a temporary window list 
    // (the window Z order can change dynamically)
    GLWindow *wins[MAX_WIN];
    memcpy(wins,allWin,nbWindow*sizeof(GLWindow *));

    BOOL processed = FALSE;

    // Menus
    i = nbWindow-1;
    while(i>=0 && !processed) {
      wins[i]->ManageMenu(evt);
      processed = wins[i]->IsEventProcessed();
      if(!processed) i--;
    }
    if( !processed ) {
      // Components
      i = nbWindow-1;
      while(i>=0 && !processed) {
        wins[i]->ManageEvent(evt);
        processed = wins[i]->IsEventProcessed();
        if(!processed) i--;
      }
    }

    // Close menu
    if( processed ) {
      for(int j=0;j<nbWindow;j++) 
        if( j!=i && evt->type==SDL_MOUSEBUTTONDOWN ) wins[j]->CloseMenu();
    }

  }

  if( evt->type == SDL_MOUSEBUTTONDOWN ) {
    // Keboard focus
    keyFocus = i;
  }
  return (i==0); // Processed by top level

}

// -------------------------------------------

void GLWindowManager::RemoveAccFromStr(char *txt,char *acc,int *pos,int *width) {

  int i,j;
  char tmp[256];
  char c[2];

  if( !txt ) return;

  GLFont2D *fnt = GLToolkit::GetDialogFont();
  if(acc)   *acc   = 0;
  if(pos)   *pos   = 0;
  if(width) *width = 0;
  c[1]=0;

  for(i=0,j=0;i<(int)strlen(txt);i++) {
    if(txt[i]!='&') {
      tmp[j++]=txt[i];
    } else {
      if(acc) *acc=tolower(txt[i+1]);
      tmp[j]=0;
      if(pos) *pos=fnt->GetTextWidth(tmp);
      c[0]=txt[i+1];
      if(width) *width=fnt->GetTextWidth(c);
    }
  }

  tmp[j]=0;
  strcpy(txt,tmp);

}

// -------------------------------------------

char *GLWindowManager::GetAccStr(int keyCode,int keyModifier) {

    static char tmp[32];
    char A[32];

    switch( keyCode ) {

      case SDLK_BACKSPACE:
        strcpy(A,"Back");
        break;
	    case SDLK_TAB:
        strcpy(A,"Tab");
        break;
	    case SDLK_CLEAR:
        strcpy(A,"Clear");
        break;
	    case SDLK_RETURN:
        strcpy(A,"Return");
        break;
	    case SDLK_PAUSE:
        strcpy(A,"Pause");
        break;
	    case SDLK_ESCAPE:
        strcpy(A,"Esc");
        break;
	    case SDLK_DELETE:
        strcpy(A,"Del");
        break;
	    case SDLK_UP:
        strcpy(A,"\206");
        break;
	    case SDLK_DOWN:
        strcpy(A,"\205");
        break;
	    case SDLK_RIGHT:
        strcpy(A,"\203");
        break;
      case SDLK_LEFT:
        strcpy(A,"\204");
        break;
	    case SDLK_INSERT:
        strcpy(A,"Ins");
        break;
	    case SDLK_HOME:
        strcpy(A,"Home");
        break;
	    case SDLK_END:
        strcpy(A,"End");
        break;
	    case SDLK_PAGEUP:
        strcpy(A,"PageUp");
        break;
	    case SDLK_PAGEDOWN:
        strcpy(A,"PageDown");
        break;
	    case SDLK_F1:
        strcpy(A,"F1");
        break;
	    case SDLK_F2:
        strcpy(A,"F2");
        break;
	    case SDLK_F3:
        strcpy(A,"F3");
        break;
	    case SDLK_F4:
        strcpy(A,"F4");
        break;
	    case SDLK_F5:
        strcpy(A,"F5");
        break;
	    case SDLK_F6:
        strcpy(A,"F6");
        break;
	    case SDLK_F7:
        strcpy(A,"F7");
        break;
	    case SDLK_F8:
        strcpy(A,"F8");
        break;
	    case SDLK_F9:
        strcpy(A,"F9");
        break;
	    case SDLK_F10:
        strcpy(A,"F10");
        break;
	    case SDLK_F11:
        strcpy(A,"F11");
        break;
	    case SDLK_F12:
        strcpy(A,"F12");
        break;
	    case SDLK_F13:
        strcpy(A,"F13");
        break;
	    case SDLK_F14:
        strcpy(A,"F14");
        break;
	    case SDLK_F15:
        strcpy(A,"F15");
        break;
      default:
        A[0] = toupper(keyCode);
        A[1] = 0;
    }

    switch(keyModifier) {
      case ALT_MODIFIER:
        strcpy(tmp,"Alt+");
        break;
      case CTRL_MODIFIER:
        strcpy(tmp,"Ctrl+");
        break;
      case SHIFT_MODIFIER:
        strcpy(tmp,"Shift+");
        break;
      default:
        strcpy(tmp,"");
    }
    strcat(tmp,A);
    return tmp;

}
