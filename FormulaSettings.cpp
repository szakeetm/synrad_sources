/*
  File:        FormulaSettings.cpp
  Description: Formula edition dialog
  Program:     SynRad
  Author:      R. KERSEVAN / M ADY
  Copyright:   E.S.R.F / CERN

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "FormulaSettings.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"

FormulaSettings::FormulaSettings():GLWindow() {

	formulaId = -1;

  int wD = 500;
  int hD = 425;

  SetTitle("Formula Editor");

  
  exprL = new GLLabel("Expression");
  exprL->SetBounds(5,5,75,18);
  Add(exprL);
  exprT = new GLTextField(0,"");
  exprT->SetBounds(85,5,wD-65,18);
  Add(exprT);
  
  nameL = new GLLabel("Name (optional)");
  nameL->SetBounds(5,30,75,18);
  Add(nameL);
  nameT = new GLTextField(0,"");
  nameT->SetBounds(85,30,wD-65,18);
  Add(nameT);
  

  descL = new GLLabel(
    "MC Variables: An (Absorption on facet n), Dn (Desorption on facet n), Hn (Hit on facet n)\n"
	"                 Fn (Flux absorbed on facet n), Pn (Power absorbed on facet n)\n"
    "                 SUMABS (total absorbed), SUMDES (total desorbed), SUMHIT (total hit)\n"
	"                 SUMFLUX (total gen. flux), SUMPOWER (total gen. power), SCANS (no. of scans)\n\n"
    "Area variables: ARn (Area of facet n), ABSAR (total absorption area)\n\n"
    "Math functions: sin(), cos(), tan(), sinh(), cosh(), tanh()\n"
    "                   asin(), acos(), atan(), exp(), ln(), pow(x,y)\n"
    "                   log2(), log10(), inv(), sqrt(), abs()\n\n"
    "Utils functions: ci95(p,N) 95% confidence interval (p=prob,N=count)\n"
    "                  sum(prefix,i,j) sum variables ex: sum(F,1,10)=F1+F2+...+F10 \n"
	"                  sum(prefix,Si) sum of variables on selection group i  ex: sum(F,S1)=all flux on group 1 \n\n"
    "Constants: Kb (Boltzmann's constant [J.K\270\271]), R (Gas constant [J.K\270\271.mol\270\271])\n"
    "              Na (Avogadro's number [mol\270\271]), PI\n\n"
    "Expression example:\n"
    "  (A1+A45)/(D23+D12)\n"
    "  sqrt(A1^2+A45^2)*DESAR/SUMDES\n"    
    );
  descL->SetBounds(5,55,wD-10,hD-100);
  Add(descL);

  applyButton = new GLButton(0,"Create");
  applyButton->SetBounds(wD-300,hD-43,95,19);
  Add(applyButton);

  deleteButton = new GLButton(0,"Delete");
  deleteButton->SetEnabled(FALSE);
  deleteButton->SetBounds(wD-200,hD-43,95,19);
  Add(deleteButton);

  cancelButton = new GLButton(0,"Cancel");
  cancelButton->SetBounds(wD-100,hD-43,95,19);
  Add(cancelButton);

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  RestoreDeviceObjects();

}