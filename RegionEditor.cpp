/*
File:        RegionEditor.cpp
Description: Region editor window
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
#include "RegionEditor.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/MathTools.h" //PI
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLInputBox.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLCombo.h"

#include "Synrad.h"

extern GLApplication *theApp;

// --------------------------------------------------------------------

RegionEditor::RegionEditor():GLWindow() {

	int wD = 773;
	int hD = 683;
	beamPanel = new GLTitledPanel("Beam geometry");
	beamPanel->SetBounds(12, 12, 745, 160);
	Add(beamPanel);
	particlePanel = new GLTitledPanel("Particle settings");
	particlePanel->SetBounds(12, 178, 745, 81);
	Add(particlePanel);
	beamParamPanel = new GLTitledPanel("Additional beam parameters");
	beamParamPanel->SetBounds(12, 265, 745, 155);
	Add(beamParamPanel);
	photonGenPanel = new GLTitledPanel("Photon generation");
	photonGenPanel->SetBounds(13, 426, 744, 101);
	Add(photonGenPanel);
	magPanel = new GLTitledPanel("Magnetic Field");
	magPanel->SetBounds(13, 533, 744, 97);
	Add(magPanel);
	label6 = new GLLabel("cm");
	beamPanel->SetCompBounds(label6, 631, 27, 21, 13);
	beamPanel->Add(label6);

	label4 = new GLLabel("cm");
	beamPanel->SetCompBounds(label4, 449, 27, 21, 13);
	beamPanel->Add(label4);

	startPointYtext = new GLTextField(0, "");
	beamPanel->SetCompBounds(startPointYtext, 343, 24, 100, 20);
	beamPanel->Add(startPointYtext);

	startPointZtext = new GLTextField(0, "");
	beamPanel->SetCompBounds(startPointZtext, 525, 24, 100, 20);
	beamPanel->Add(startPointZtext);

	label7 = new GLLabel("Z0:");
	beamPanel->SetCompBounds(label7, 501, 27, 23, 13);
	beamPanel->Add(label7);

	label5 = new GLLabel("Y0:");
	beamPanel->SetCompBounds(label5, 319, 27, 23, 13);
	beamPanel->Add(label5);

	label3 = new GLLabel("cm");
	beamPanel->SetCompBounds(label3, 274, 27, 21, 13);
	beamPanel->Add(label3);

	startPointXtext = new GLTextField(0, "");
	beamPanel->SetCompBounds(startPointXtext, 168, 24, 100, 20);
	beamPanel->Add(startPointXtext);

	label2 = new GLLabel("X0:");
	beamPanel->SetCompBounds(label2, 142, 27, 23, 13);
	beamPanel->Add(label2);

	label1 = new GLLabel("Beam start position:");
	beamPanel->SetCompBounds(label1, 8, 27, 99, 13);
	beamPanel->Add(label1);

	startDirDefinitionCombo = new GLCombo(0);
	beamPanel->SetCompBounds(startDirDefinitionCombo, 168, 50, 100, 21);
	beamPanel->Add(startDirDefinitionCombo);
	startDirDefinitionCombo->SetSize(2);
	startDirDefinitionCombo->SetValueAt(0, "By vector:");
	startDirDefinitionCombo->SetValueAt(1, "By angle:");
	startDirDefinitionCombo->SetSelectedIndex(0); //Default

	label8 = new GLLabel("Beam start direction:");
	beamPanel->SetCompBounds(label8, 8, 53, 103, 13);
	beamPanel->Add(label8);

	label9 = new GLLabel("rad");
	beamPanel->SetCompBounds(label9, 631, 53, 22, 13);
	beamPanel->Add(label9);

	label10 = new GLLabel("rad");
	beamPanel->SetCompBounds(label10, 449, 53, 22, 13);
	beamPanel->Add(label10);

	theta0text = new GLTextField(0, "");
	beamPanel->SetCompBounds(theta0text, 343, 50, 100, 20);
	beamPanel->Add(theta0text);

	alpha0text = new GLTextField(0, "");
	beamPanel->SetCompBounds(alpha0text, 525, 50, 100, 20);
	beamPanel->Add(alpha0text);

	label11 = new GLLabel("Alpha0:");
	beamPanel->SetCompBounds(label11, 481, 53, 43, 13);
	beamPanel->Add(label11);

	label12 = new GLLabel("Theta0:");
	beamPanel->SetCompBounds(label12, 298, 53, 44, 13);
	beamPanel->Add(label12);

	startDirInfoButton = new GLButton(0, "Info");
	beamPanel->SetCompBounds(startDirInfoButton, 659, 49, 40, 21);
	beamPanel->Add(startDirInfoButton);

	startDirZtext = new GLTextField(0, "");
	beamPanel->SetCompBounds(startDirZtext, 525, 77, 100, 20);
	beamPanel->Add(startDirZtext);

	label15 = new GLLabel("dirZ:");
	beamPanel->SetCompBounds(label15, 496, 80, 28, 13);
	beamPanel->Add(label15);

	startDirYtext = new GLTextField(0, "");
	beamPanel->SetCompBounds(startDirYtext, 343, 77, 100, 20);
	beamPanel->Add(startDirYtext);

	label14 = new GLLabel("dirY:");
	beamPanel->SetCompBounds(label14, 314, 80, 28, 13);
	beamPanel->Add(label14);

	startDirXtext = new GLTextField(0, "");
	beamPanel->SetCompBounds(startDirXtext, 168, 77, 100, 20);
	beamPanel->Add(startDirXtext);

	label13 = new GLLabel("dirX:");
	beamPanel->SetCompBounds(label13, 137, 80, 28, 13);
	beamPanel->Add(label13);

	label19 = new GLLabel("cm");
	beamPanel->SetCompBounds(label19, 631, 132, 21, 13);
	beamPanel->Add(label19);

	label20 = new GLLabel("cm");
	beamPanel->SetCompBounds(label20, 449, 132, 21, 13);
	beamPanel->Add(label20);

	limitsYtext = new GLTextField(0, "");
	beamPanel->SetCompBounds(limitsYtext, 343, 129, 100, 20);
	beamPanel->Add(limitsYtext);

	limitsZtext = new GLTextField(0, "");
	beamPanel->SetCompBounds(limitsZtext, 525, 129, 100, 20);
	beamPanel->Add(limitsZtext);

	label21 = new GLLabel("Zmax:");
	beamPanel->SetCompBounds(label21, 488, 132, 36, 13);
	beamPanel->Add(label21);

	label22 = new GLLabel("Ymax:");
	beamPanel->SetCompBounds(label22, 306, 132, 36, 13);
	beamPanel->Add(label22);

	label23 = new GLLabel("cm");
	beamPanel->SetCompBounds(label23, 274, 132, 21, 13);
	beamPanel->Add(label23);

	limitsXtext = new GLTextField(0, "");
	beamPanel->SetCompBounds(limitsXtext, 168, 129, 100, 20);
	beamPanel->Add(limitsXtext);

	label24 = new GLLabel("Xmax:");
	beamPanel->SetCompBounds(label24, 129, 132, 36, 13);
	beamPanel->Add(label24);

	label25 = new GLLabel("Calculation boundary:");
	beamPanel->SetCompBounds(label25, 8, 132, 109, 13);
	beamPanel->Add(label25);

	label16 = new GLLabel("cm");
	beamPanel->SetCompBounds(label16, 274, 106, 21, 13);
	beamPanel->Add(label16);

	dLtext = new GLTextField(0, "");
	beamPanel->SetCompBounds(dLtext, 168, 103, 100, 20);
	beamPanel->Add(dLtext);

	label17 = new GLLabel("dL:");
	beamPanel->SetCompBounds(label17, 143, 106, 22, 13);
	beamPanel->Add(label17);

	label18 = new GLLabel("Trajectory step length:");
	beamPanel->SetCompBounds(label18, 8, 106, 112, 13);
	beamPanel->Add(label18);

	limitsInfoButton = new GLButton(0, "Info");
	beamPanel->SetCompBounds(limitsInfoButton, 659, 128, 40, 21);
	beamPanel->Add(limitsInfoButton);

	dlInfoButton = new GLButton(0, "Info");
	beamPanel->SetCompBounds(dlInfoButton, 297, 102, 40, 21);
	beamPanel->Add(dlInfoButton);

	label27 = new GLLabel("GeV/c2");
	particlePanel->SetCompBounds(label27, 274, 25, 45, 13);
	particlePanel->Add(label27);

	label26 = new GLLabel("Particle mass:");
	particlePanel->SetCompBounds(label26, 8, 25, 72, 13);
	particlePanel->Add(label26);

	particleMassText = new GLTextField(0, "");
	particlePanel->SetCompBounds(particleMassText, 168, 22, 100, 20);
	particlePanel->Add(particleMassText);

	label28 = new GLLabel("pmass:");
	particlePanel->SetCompBounds(label28, 125, 25, 40, 13);
	particlePanel->Add(label28);

	particleChargeCombo = new GLCombo(0);
	particlePanel->SetCompBounds(particleChargeCombo, 168, 48, 100, 21);
	particlePanel->Add(particleChargeCombo);
	particleChargeCombo->SetSize(2);
	particleChargeCombo->SetValueAt(0, "Positive");
	particleChargeCombo->SetValueAt(1, "Negative");

	label29 = new GLLabel("Particle charge:");
	particlePanel->SetCompBounds(label29, 9, 51, 81, 13);
	particlePanel->Add(label29);

	setParticleProtonButton = new GLButton(0, "Proton");
	particlePanel->SetCompBounds(setParticleProtonButton, 558, 47, 100, 21);
	particlePanel->Add(setParticleProtonButton);

	setParticlePositronButton = new GLButton(0, "Positron");
	particlePanel->SetCompBounds(setParticlePositronButton, 452, 47, 100, 21);
	particlePanel->Add(setParticlePositronButton);

	setParticleElectronButton = new GLButton(0, "Electron");
	particlePanel->SetCompBounds(setParticleElectronButton, 343, 47, 100, 21);
	particlePanel->Add(setParticleElectronButton);

	label30 = new GLLabel("Quick set to:");
	particlePanel->SetCompBounds(label30, 468, 25, 67, 13);
	particlePanel->Add(label30);

	label33 = new GLLabel("mA");
	beamParamPanel->SetCompBounds(label33, 355, 22, 22, 13);
	beamParamPanel->Add(label33);

	beamCurrentText = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(beamCurrentText, 275, 19, 76, 20);
	beamParamPanel->Add(beamCurrentText);

	label34 = new GLLabel("Current:");
	beamParamPanel->SetCompBounds(label34, 231, 22, 44, 13);
	beamParamPanel->Add(label34);

	label31 = new GLLabel("GeV");
	beamParamPanel->SetCompBounds(label31, 172, 23, 28, 13);
	beamParamPanel->Add(label31);

	beamEnergyText = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(beamEnergyText, 91, 19, 76, 20);
	beamParamPanel->Add(beamEnergyText);

	label32 = new GLLabel("Beam energy:");
	beamParamPanel->SetCompBounds(label32, 10, 22, 72, 13);
	beamParamPanel->Add(label32);

	idealBeamToggle = new GLToggle(0, "Ideal Beam");
	beamParamPanel->SetCompBounds(idealBeamToggle, 9, 55, 79, 17);
	beamParamPanel->Add(idealBeamToggle);

	emittanceCouplingToggle = new GLToggle(0, "Emittance:");
	beamParamPanel->SetCompBounds(emittanceCouplingToggle, 9, 78, 76, 17);
	beamParamPanel->Add(emittanceCouplingToggle);

	emittanceYtext = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(emittanceYtext, 275, 102, 76, 20);
	beamParamPanel->Add(emittanceYtext);

	emittanceXYtoggle = new GLToggle(0, "X emittance:");
	beamParamPanel->SetCompBounds(emittanceXYtoggle, 9, 104, 85, 17);
	beamParamPanel->Add(emittanceXYtoggle);

	couplingText = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(couplingText, 275, 76, 76, 20);
	beamParamPanel->Add(couplingText);

	emittanceXtext = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(emittanceXtext, 91, 102, 76, 20);
	beamParamPanel->Add(emittanceXtext);

	emittanceText = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(emittanceText, 91, 76, 76, 20);
	beamParamPanel->Add(emittanceText);

	label36 = new GLLabel("Y emittance:");
	beamParamPanel->SetCompBounds(label36, 205, 105, 66, 13);
	beamParamPanel->Add(label36);

	label35 = new GLLabel("Y/X coupling:");
	beamParamPanel->SetCompBounds(label35, 199, 79, 72, 13);
	beamParamPanel->Add(label35);

	beamsizeInfoButton = new GLButton(0, "Info");
	beamParamPanel->SetCompBounds(beamsizeInfoButton, 9, 126, 40, 21);
	beamParamPanel->Add(beamsizeInfoButton);

	label40 = new GLLabel("%");
	beamParamPanel->SetCompBounds(label40, 357, 79, 15, 13);
	beamParamPanel->Add(label40);

	label39 = new GLLabel("cm");
	beamParamPanel->SetCompBounds(label39, 357, 105, 21, 13);
	beamParamPanel->Add(label39);

	label38 = new GLLabel("cm");
	beamParamPanel->SetCompBounds(label38, 172, 105, 21, 13);
	beamParamPanel->Add(label38);

	label37 = new GLLabel("cm");
	beamParamPanel->SetCompBounds(label37, 172, 79, 21, 13);
	beamParamPanel->Add(label37);

	label49 = new GLLabel("%");
	beamParamPanel->SetCompBounds(label49, 550, 53, 15, 13);
	beamParamPanel->Add(label49);

	label47 = new GLLabel("cm");
	beamParamPanel->SetCompBounds(label47, 712, 105, 21, 13);
	beamParamPanel->Add(label47);

	label48 = new GLLabel("EtaX\':");
	beamParamPanel->SetCompBounds(label48, 591, 105, 35, 13);
	beamParamPanel->Add(label48);

	label50 = new GLLabel("En. spread:");
	beamParamPanel->SetCompBounds(label50, 404, 53, 61, 13);
	beamParamPanel->Add(label50);

	etaPrimeText = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(etaPrimeText, 631, 101, 75, 20);
	beamParamPanel->Add(etaPrimeText);

	energySpreadText = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(energySpreadText, 471, 50, 74, 20);
	beamParamPanel->Add(energySpreadText);

	label45 = new GLLabel("cm");
	beamParamPanel->SetCompBounds(label45, 712, 79, 21, 13);
	beamParamPanel->Add(label45);

	label46 = new GLLabel("EtaX:");
	beamParamPanel->SetCompBounds(label46, 591, 79, 33, 13);
	beamParamPanel->Add(label46);

	etaText = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(etaText, 632, 75, 74, 20);
	beamParamPanel->Add(etaText);

	label43 = new GLLabel("cm");
	beamParamPanel->SetCompBounds(label43, 550, 105, 21, 13);
	beamParamPanel->Add(label43);

	label44 = new GLLabel("BetaY:");
	beamParamPanel->SetCompBounds(label44, 425, 105, 39, 13);
	beamParamPanel->Add(label44);

	betaYtext = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(betaYtext, 471, 102, 74, 20);
	beamParamPanel->Add(betaYtext);

	label41 = new GLLabel("cm");
	beamParamPanel->SetCompBounds(label41, 550, 79, 21, 13);
	beamParamPanel->Add(label41);

	label42 = new GLLabel("BetaX:");
	beamParamPanel->SetCompBounds(label42, 425, 79, 39, 13);
	beamParamPanel->Add(label42);

	betaXtext = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(betaXtext, 471, 76, 74, 20);
	beamParamPanel->Add(betaXtext);

	gammaInfoButton = new GLButton(0, "Info");
	beamParamPanel->SetCompBounds(gammaInfoButton, 580, 128, 50, 21);
	beamParamPanel->Add(gammaInfoButton);

	BXYfileNameText = new GLTextField(0, "");
	beamParamPanel->SetCompBounds(BXYfileNameText, 275, 128, 298, 20);
	beamParamPanel->Add(BXYfileNameText);

	label51 = new GLLabel("BXY file:");
	beamParamPanel->SetCompBounds(label51, 224, 131, 47, 13);
	beamParamPanel->Add(label51);

	constantBXYtoggle = new GLToggle(0, "Constant beta function");
	beamParamPanel->SetCompBounds(constantBXYtoggle, 91, 130, 133, 17);
	beamParamPanel->Add(constantBXYtoggle);

	bxyEditButton = new GLButton(0, "Edit");
	beamParamPanel->SetCompBounds(bxyEditButton, 691, 128, 50, 21);
	beamParamPanel->Add(bxyEditButton);

	bxyBrowseButton = new GLButton(0, "...");
	beamParamPanel->SetCompBounds(bxyBrowseButton, 635, 128, 50, 21);
	beamParamPanel->Add(bxyBrowseButton);

	enableParPolarizationToggle = new GLToggle(0, "Parallel");
	photonGenPanel->SetCompBounds(enableParPolarizationToggle, 276, 50, 60, 17);
	photonGenPanel->Add(enableParPolarizationToggle);

	label59 = new GLLabel("Generated photon energy:");
	photonGenPanel->SetCompBounds(label59, 8, 51, 131, 13);
	photonGenPanel->Add(label59);

	label57 = new GLLabel("eV");
	photonGenPanel->SetCompBounds(label57, 576, 25, 20, 13);
	photonGenPanel->Add(label57);

	EmaxText = new GLTextField(0, "");
	photonGenPanel->SetCompBounds(EmaxText, 472, 22, 100, 20);
	photonGenPanel->Add(EmaxText);

	label58 = new GLLabel("Emax:");
	photonGenPanel->SetCompBounds(label58, 432, 25, 36, 13);
	photonGenPanel->Add(label58);

	label54 = new GLLabel("eV");
	photonGenPanel->SetCompBounds(label54, 379, 25, 20, 13);
	photonGenPanel->Add(label54);

	label56 = new GLLabel("Generated photon energy:");
	photonGenPanel->SetCompBounds(label56, 7, 25, 131, 13);
	photonGenPanel->Add(label56);

	EminText = new GLTextField(0, "");
	photonGenPanel->SetCompBounds(EminText, 276, 22, 100, 20);
	photonGenPanel->Add(EminText);

	label55 = new GLLabel("Emin:");
	photonGenPanel->SetCompBounds(label55, 241, 25, 33, 13);
	photonGenPanel->Add(label55);

	enableOrtPolarizationToggle = new GLToggle(0, "Orthogonal");
	photonGenPanel->SetCompBounds(enableOrtPolarizationToggle, 472, 50, 78, 17);
	photonGenPanel->Add(enableOrtPolarizationToggle);

	label61 = new GLLabel("rad");
	photonGenPanel->SetCompBounds(label61, 576, 76, 22, 13);
	photonGenPanel->Add(label61);

	label60 = new GLLabel("rad");
	photonGenPanel->SetCompBounds(label60, 379, 76, 22, 13);
	photonGenPanel->Add(label60);

	psiMaxXtext = new GLTextField(0, "");
	photonGenPanel->SetCompBounds(psiMaxXtext, 276, 73, 100, 20);
	photonGenPanel->Add(psiMaxXtext);

	psiMaxYtext = new GLTextField(0, "");
	photonGenPanel->SetCompBounds(psiMaxYtext, 472, 73, 100, 20);
	photonGenPanel->Add(psiMaxYtext);

	limitAngleToggle = new GLToggle(0, "Limit photon divergence:");
	photonGenPanel->SetCompBounds(limitAngleToggle, 10, 75, 142, 17);
	photonGenPanel->Add(limitAngleToggle);

	label63 = new GLLabel("divYmax:");
	photonGenPanel->SetCompBounds(label63, 418, 76, 50, 13);
	photonGenPanel->Add(label63);

	label62 = new GLLabel("divXmax:");
	photonGenPanel->SetCompBounds(label62, 224, 76, 50, 13);
	photonGenPanel->Add(label62);

	BxtypeCombo = new GLCombo(0);
	magPanel->SetCompBounds(BxtypeCombo, 35, 18, 135, 21);
	magPanel->Add(BxtypeCombo);
	BxtypeCombo->SetSize(8);
	BxtypeCombo->SetValueAt(0, "Constant field");
	BxtypeCombo->SetValueAt(1, "Coords along a direction");
	BxtypeCombo->SetValueAt(2, "Coords along the beam");
	BxtypeCombo->SetValueAt(3, "Sine / cosine");
	BxtypeCombo->SetValueAt(4, "Quadrupole");
	BxtypeCombo->SetValueAt(5, "Analytical expression");
	BxtypeCombo->SetValueAt(6, "Helicoidal");
	BxtypeCombo->SetValueAt(7, "Rotating dipole");

	label66 = new GLLabel("Bx:");
	magPanel->SetCompBounds(label66, 7, 22, 22, 13);
	magPanel->Add(label66);

	label65 = new GLLabel("T");
	magPanel->SetCompBounds(label65, 273, 21, 14, 13);
	magPanel->Add(label65);

	magxEditButton = new GLButton(0, "Edit");
	magPanel->SetCompBounds(magxEditButton, 690, 18, 48, 21);
	magPanel->Add(magxEditButton);

	constBXtext = new GLTextField(0, "");
	magPanel->SetCompBounds(constBXtext, 176, 19, 91, 20);
	magPanel->Add(constBXtext);

	magxBrowseButton = new GLButton(0, "...");
	magPanel->SetCompBounds(magxBrowseButton, 634, 18, 50, 21);
	magPanel->Add(magxBrowseButton);

	MAGfileXtext = new GLTextField(0, "");
	magPanel->SetCompBounds(MAGfileXtext, 349, 19, 279, 20);
	magPanel->Add(MAGfileXtext);

	label64 = new GLLabel("MAG file:");
	magPanel->SetCompBounds(label64, 297, 22, 50, 13);
	magPanel->Add(label64);

	BztypeCombo = new GLCombo(0);
	magPanel->SetCompBounds(BztypeCombo, 35, 70, 135, 21);
	magPanel->Add(BztypeCombo);
	BztypeCombo->SetSize(8);
	BztypeCombo->SetValueAt(0, "Constant field");
	BztypeCombo->SetValueAt(1, "Coords along a direction");
	BztypeCombo->SetValueAt(2, "Coords along the beam");
	BztypeCombo->SetValueAt(3, "Sine / cosine");
	BztypeCombo->SetValueAt(4, "Quadrupole");
	BztypeCombo->SetValueAt(5, "Analytical expression");
	BztypeCombo->SetValueAt(6, "Helicoidal");
	BztypeCombo->SetValueAt(7, "Rotating dipole");

	label70 = new GLLabel("Bz:");
	magPanel->SetCompBounds(label70, 7, 74, 22, 13);
	magPanel->Add(label70);

	label71 = new GLLabel("T");
	magPanel->SetCompBounds(label71, 273, 73, 14, 13);
	magPanel->Add(label71);

	magzEditButton = new GLButton(0, "Edit");
	magPanel->SetCompBounds(magzEditButton, 690, 70, 48, 21);
	magPanel->Add(magzEditButton);

	constBZtext = new GLTextField(0, "");
	magPanel->SetCompBounds(constBZtext, 176, 71, 91, 20);
	magPanel->Add(constBZtext);

	magzBrowseButton = new GLButton(0, "...");
	magPanel->SetCompBounds(magzBrowseButton, 634, 70, 50, 21);
	magPanel->Add(magzBrowseButton);

	MAGfileZtext = new GLTextField(0, "");
	magPanel->SetCompBounds(MAGfileZtext, 349, 71, 279, 20);
	magPanel->Add(MAGfileZtext);

	label72 = new GLLabel("MAG file:");
	magPanel->SetCompBounds(label72, 297, 74, 50, 13);
	magPanel->Add(label72);

	BytypeCombo = new GLCombo(0);
	magPanel->SetCompBounds(BytypeCombo, 35, 44, 135, 21);
	magPanel->Add(BytypeCombo);
	BytypeCombo->SetSize(8);
	BytypeCombo->SetValueAt(0, "Constant field");
	BytypeCombo->SetValueAt(1, "Coords along a direction");
	BytypeCombo->SetValueAt(2, "Coords along the beam");
	BytypeCombo->SetValueAt(3, "Sine / cosine");
	BytypeCombo->SetValueAt(4, "Quadrupole");
	BytypeCombo->SetValueAt(5, "Analytical expression");
	BytypeCombo->SetValueAt(6, "Helicoidal");
	BytypeCombo->SetValueAt(7, "Rotating dipole");

	label67 = new GLLabel("By:");
	magPanel->SetCompBounds(label67, 7, 48, 22, 13);
	magPanel->Add(label67);

	label68 = new GLLabel("T");
	magPanel->SetCompBounds(label68, 273, 47, 14, 13);
	magPanel->Add(label68);

	magyEditButton = new GLButton(0, "Edit");
	magPanel->SetCompBounds(magyEditButton, 690, 44, 48, 21);
	magPanel->Add(magyEditButton);

	constBYtext = new GLTextField(0, "");
	magPanel->SetCompBounds(constBYtext, 176, 45, 91, 20);
	magPanel->Add(constBYtext);

	magyBrowseButton = new GLButton(0, "...");
	magPanel->SetCompBounds(magyBrowseButton, 634, 44, 50, 21);
	magPanel->Add(magyBrowseButton);

	MAGfileYtext = new GLTextField(0, "");
	magPanel->SetCompBounds(MAGfileYtext, 349, 45, 279, 20);
	magPanel->Add(MAGfileYtext);

	label69 = new GLLabel("MAG file:");
	magPanel->SetCompBounds(label69, 297, 48, 50, 13);
	magPanel->Add(label69);

	applyButton = new GLButton(0, "Apply && Recalculate points");
	applyButton->SetBounds(289, 636, 223, 23);
	Add(applyButton);

	SetTitle("RegionEditor");
	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	RestoreDeviceObjects();
}


// --------------------------------------------------------------------

void RegionEditor::Display(Worker *w,int Id) {

	//set textfields and values here

	worker = w;
	this->regionId=Id;
	cr = &(worker->regions[regionId]);
	SetTitle("Editing region " + std::to_string(regionId + 1));

	FillValues();
	EnableDisableControls(NULL);

	SetVisible(true);
}

int RegionEditor::GetRegionId()
{
	return this->regionId;
}


void RegionEditor::ProcessMessage(GLComponent *src,int message) {
	SynRad *mApp = (SynRad *)theApp;
	switch(message) {
	case MSG_BUTTON:

		if (src==setParticleElectronButton) {
			particleMassText->SetText("0.0005110034");
			particleChargeCombo->SetSelectedIndex(1);
		} else if (src==setParticlePositronButton) {
			particleMassText->SetText("0.0005110034");
			particleChargeCombo->SetSelectedIndex(0);
		} else if (src==setParticleProtonButton) {
			particleMassText->SetText("0.938272046");
			particleChargeCombo->SetSelectedIndex(0);
		} else if (src==this->bxyEditButton) {
			char tmp[512];
			sprintf(tmp,"notepad.exe \"%s\"",BXYfileNameText->GetText());
			StartProc(tmp,STARTPROC_FOREGROUND);
		} else if (src==this->magxEditButton) {
			char tmp[512];
			sprintf(tmp,"notepad.exe \"%s\"",MAGfileXtext->GetText());
			StartProc(tmp, STARTPROC_FOREGROUND);
		} else if (src==this->magyEditButton) {
			char tmp[512];
			sprintf(tmp,"notepad.exe \"%s\"",MAGfileYtext->GetText());
			StartProc(tmp, STARTPROC_FOREGROUND);
		} else if (src==this->magzEditButton) {
			char tmp[512];
			sprintf(tmp,"notepad.exe \"%s\"",MAGfileZtext->GetText());
			StartProc(tmp, STARTPROC_FOREGROUND);
		} else if (src==this->bxyBrowseButton) {
			//load file dialog
			FILENAME *bxyFile=GLFileBox::OpenFile(NULL,NULL,"Open BXY file","BXY files\0*.bxy\0All files\0*.*\0",2);
			if (!bxyFile) return;
			if (!(bxyFile->fullName)) return;
			BXYfileNameText->SetText(bxyFile->fullName);
		} else if (src==this->magxBrowseButton) {
			//load file dialog
			FILENAME *magFile=GLFileBox::OpenFile(NULL,NULL,"Open MAG file","MAG files\0*.mag\0All files\0*.*\0",2);
			if (!magFile) return;
			if (!(magFile->fullName)) return;
			MAGfileXtext->SetText(magFile->fullName);
		} else if (src==this->magyBrowseButton) {
			//load file dialog
			FILENAME *magFile=GLFileBox::OpenFile(NULL,NULL,"Open MAG file","MAG files\0*.mag\0All files\0*.*\0",2);
			if (!magFile) return;
			if (!(magFile->fullName)) return;
			MAGfileYtext->SetText(magFile->fullName);
		} else if (src==this->magzBrowseButton) {
			//load file dialog
			FILENAME *magFile=GLFileBox::OpenFile(NULL,NULL,"Open MAG file","MAG files\0*.mag\0All files\0*.*\0",2);
			if (!magFile) return;
			if (!(magFile->fullName)) return;
			MAGfileZtext->SetText(magFile->fullName);
		} else if (src==this->startDirInfoButton) {
			char tmp[]="Theta [-PI..PI]: angle with Z axis in XZ plane (positive if X<0)\n"
				"Alpha [-PI/2..PI/2]: angle with XZ plane (positive if Y<0)\n\n"
				"Conversion:\n"
				"X=-cos(alpha)*sin(theta)\n"
				"Y=-sin(alpha)\n"
				"Z=cos(alpha)*cos(theta)\n";
			GLMessageBox::Display(tmp,"Beam start direction",GLDLG_OK,GLDLG_ICONINFO);
		} else if (src==this->beamsizeInfoButton) {
			char tmp[]="horizontal_emittance = emittance / (1 + coupling)\n"
				"vertical_emittance = horizontal_emittance * coupling\n\n"
				"sigma_x = sqrt( horizontal_emittance * beta_X + (eta * energy_spread)^2 )\n"
				"sigma_y = vertical_emittance * beta_Y\n"
				"sigma_x_prime = sqrt( horizontal_emittance / beta_X + (eta_prime * energy_spread)^2 )\n"
				"sigma_y_prime=sqrt( vertical_emittance / beta_Y )\n\n"
				"The beam size follows a bivariate Gaussian distribution with "
				"sigma_x and sigma_y standard deviations.\n"
				"The beam (angular) divergence follows a Gaussian distribution with "
				"sigma_x_prime and sigma_y_prime distributions.";
			GLMessageBox::Display(tmp,"Beam size info",GLDLG_OK,GLDLG_ICONINFO);
		} else if (src==this->dlInfoButton) {
			char tmp[]="The trajectory will consist of N calculated points situated 'dL' distance from each other.\n\n"
				"The calculation algorithm supposes that the trajectory is straight for 'dL' distance, after which\n"
				"it determines the local magnetic field and calculates the new direction\n\n"
				"For strong magnetic fields (sharp curves) you should define a small value, otherwise you can set it higher\n\n"
				"The maximum number of trajectory points is 1 million. Having too many points will result in more memory usage\n"
				"but won't affect simulation speed.\n";
			GLMessageBox::Display(tmp,"Trajectory step length",GLDLG_OK,GLDLG_ICONINFO);
		} else if (src==this->limitsInfoButton) {
			char tmp[]="The trajectory will be calculated starting from the 'Beam start position' defined above.\n"
				"Calculation stops when EITHER the X,Y or Z coordinate reaches the Xmax,Ymax,Zmax value.\n"
				"(If, for example, Xmax is smaller than X0, then calculation will stop when the X coordinate\n"
				"of the trajectory is equal or less than Xmax. If Xmax is larger then X0, calculation stops\n"
				"when the trajectory's X coordinate is equal or larger than Xmax.) If these limits are never reached,\n"
				"calculation stops when the maximum number of trajectory points (1 million) is reached.\n";
			GLMessageBox::Display(tmp,"Calculation boundaries",GLDLG_OK,GLDLG_ICONINFO);
		} else if (src==gammaInfoButton) {
			char tmp[] = "The gamma invariant (not to be confused with the gamma relativistic factor!) connects\n"
				"the accelerator Alpha and Beta functions through the relation:\n"
				"Gamma * Beta - Alpha^2 = 1     or     Gamma = (1 + Alpha)^2 / Beta\n"
				"Where the Alpha function describes the change of the Beta function: Alpha(s) = - 1/2 * dBeta/ds\n"
				"The importance of this is while the offset of generated photons depends only on the beam size,\n"
				"proportional to the Beta function, the angular divergence relative to the beam direction depends also\n"
				"on the change of the Beta, expressed by the Alpha function.\n"
				"If you have the Alpha and Beta function values at any point, calculate Gamma using the above expression,\n"
				"it will remain invariant through the trajectory. If you don't know the Alpha function, write 0 for an approximation.";
				
			GLMessageBox::Display(tmp,"Calculation boundaries",GLDLG_OK,GLDLG_ICONINFO);
		} else if (src==applyButton) {
			//read values, apply them and reload region
			//close the window
			
				ApplyChanges();
			

			return;

		}
		break;

	case MSG_TEXT: //redirect enter
		ProcessMessage(applyButton,MSG_BUTTON);
		break;
	case MSG_COMBO:
	case MSG_TOGGLE:
	{
		//Radio button behavior:
		std::vector<GLToggle*> radio1 = { idealBeamToggle, emittanceCouplingToggle, emittanceXYtoggle };
		if (Contains(radio1, (GLToggle*)src)) {
			for (auto comp : radio1) {
				comp->SetState(comp == src);
			}
		}
		//Enable/disable text fields:
		EnableDisableControls(src);
		break;
	}
	case MSG_TEXT_UPD:
		if (Contains({ emittanceText,couplingText }, src)) {
			double emittance, coupling;
			if (emittanceText->GetNumber(&emittance) && couplingText->GetNumber(&coupling)) {
				double emittanceX = emittance / (1.0 + coupling * 0.01);
				emittanceXtext->SetText(emittanceX);
				emittanceYtext->SetText(emittanceX * coupling * 0.01);
			}
		} else if (Contains({ emittanceXtext,emittanceYtext }, src)) {
			double emittanceX, emittanceY;
			if (emittanceXtext->GetNumber(&emittanceX) && emittanceYtext->GetNumber(&emittanceY)) {
				double coupling = 100.0 * emittanceY / emittanceX;
				emittanceText->SetText(emittanceX * (1.0 + coupling * 0.01));
				couplingText->SetText(coupling);
			}
		}
		break;
	}
	GLWindow::ProcessMessage(src,message);
}
void RegionEditor::EnableDisableControls(GLComponent* src) {
	bool defineByAngle=(startDirDefinitionCombo->GetSelectedIndex()==0);
	theta0text->SetEditable(defineByAngle);
	alpha0text->SetEditable(defineByAngle);
	startDirXtext->SetEditable(!defineByAngle);
	startDirYtext->SetEditable(!defineByAngle);
	startDirZtext->SetEditable(!defineByAngle);

	/*
	bool useMAGfile=(BxtypeCombo->GetSelectedIndex()!=0);
	constBXtext->SetEditable(!useMAGfile);
	MAGfileXtext->SetEditable(useMAGfile);
	magxBrowseButton->SetEnabled(useMAGfile);
	magxEditButton->SetEnabled(useMAGfile);

	useMAGfile=(BytypeCombo->GetSelectedIndex()!=0);
	constBYtext->SetEditable(!useMAGfile);
	MAGfileYtext->SetEditable(useMAGfile);
	magyBrowseButton->SetEnabled(useMAGfile);
	magyEditButton->SetEnabled(useMAGfile);

	useMAGfile=(BztypeCombo->GetSelectedIndex()!=0);
	constBZtext->SetEditable(!useMAGfile);
	MAGfileZtext->SetEditable(useMAGfile);
	magzBrowseButton->SetEnabled(useMAGfile);
	magzEditButton->SetEnabled(useMAGfile);*/

	
	//Control of magnetic field components
	std::vector<GLComponent*> combos = { BxtypeCombo,BytypeCombo,BztypeCombo };
	if (Contains(combos, src)) {
		std::vector<GLTextField*> constFields = { constBXtext , constBYtext, constBZtext };
		std::vector<GLTextField*> magFileFields = { MAGfileXtext , MAGfileYtext, MAGfileZtext };
		std::vector<GLButton*> browseButtons = { magxBrowseButton , magyBrowseButton, magzBrowseButton };
		std::vector<GLButton*> editButtons = { magxEditButton , magyEditButton, magzEditButton };
		
		std::vector<int> allCompModes = { B_MODE_QUADRUPOLE,B_MODE_ANALYTIC,B_MODE_ROTATING_DIPOLE };
		bool setAllComponents = Contains(allCompModes, ((GLCombo*)src)->GetSelectedIndex()+1);
		
		for (size_t i = 0; i < 3; i++) {
			int compMode = ((GLCombo*)combos[i])->GetSelectedIndex() + 1;
			bool useMagFile = (compMode != B_MODE_CONSTANT);
			
			((GLCombo*)combos[i])->SetEditable(!setAllComponents || src==combos[i]);
			constFields[i]->SetEditable((!setAllComponents || src == combos[i]) && !useMagFile);
			magFileFields[i]->SetEditable((!setAllComponents || src == combos[i]) && useMagFile);
			browseButtons[i]->SetEnabled((!setAllComponents || src == combos[i]) && useMagFile);
			editButtons[i]->SetEnabled((!setAllComponents || src == combos[i]) && useMagFile);
		}
	}

	bool useIdealBeam=idealBeamToggle->GetState();
	bool useBXYfile=(constantBXYtoggle->GetState()==0);

	energySpreadText->SetEditable(!useIdealBeam);

	constantBXYtoggle->SetEnabled(!useIdealBeam);
	BXYfileNameText->SetEditable(!useIdealBeam && useBXYfile);
	bxyBrowseButton->SetEnabled(!useIdealBeam && useBXYfile);
	bxyEditButton->SetEnabled(!useIdealBeam && useBXYfile);
	
	emittanceText->SetEditable(!useIdealBeam && emittanceCouplingToggle->GetState());
	couplingText->SetEditable(!useIdealBeam  && emittanceCouplingToggle->GetState());
	emittanceXtext->SetEditable(!useIdealBeam  && emittanceXYtoggle->GetState());
	emittanceYtext->SetEditable(!useIdealBeam  && emittanceXYtoggle->GetState());
	
	betaXtext->SetEditable(!useIdealBeam && !useBXYfile);
	betaYtext->SetEditable(!useIdealBeam && !useBXYfile);
	etaText->SetEditable(!useIdealBeam && !useBXYfile);
	etaPrimeText->SetEditable(!useIdealBeam && !useBXYfile);			

	bool limitAngle=limitAngleToggle->GetState();
	psiMaxXtext->SetEditable(limitAngle);
	psiMaxYtext->SetEditable(limitAngle);

	enableParPolarizationToggle->SetEnabled(false);
	enableOrtPolarizationToggle->SetEnabled(false);
}

void RegionEditor::FillValues() {
	startPointXtext->SetText(cr->params.startPoint.x);
	startPointYtext->SetText(cr->params.startPoint.y);
	startPointZtext->SetText(cr->params.startPoint.z);

	double t0;
	if (cr->params.startDir.z == 0) {
		if (cr->params.startDir.x >= 0) t0 = -PI / 2;
		else t0 = PI / 2;
	}
	else {
		t0 = -atan(cr->params.startDir.x / cr->params.startDir.z); //Good for -PI/2...PI/2
		if (cr->params.startDir.z <= 0) //atan out of period
			if (cr->params.startDir.x < 0) t0 += PI;
			else t0 -= PI;
	}
	theta0text->SetText(t0);
	alpha0text->SetText(-asin(cr->params.startDir.Normalized().y));
	startDirXtext->SetText(cr->params.startDir.x);
	startDirYtext->SetText(cr->params.startDir.y);
	startDirZtext->SetText(cr->params.startDir.z);

	dLtext->SetText(cr->params.dL_cm);

	limitsXtext->SetText(cr->params.limits.x);
	limitsYtext->SetText(cr->params.limits.y);
	limitsZtext->SetText(cr->params.limits.z);

	particleMassText->SetText(abs(cr->params.particleMass_GeV));
	particleChargeCombo->SetSelectedIndex(cr->params.particleMass_GeV<0);

	beamEnergyText->SetText(cr->params.E_GeV);
	beamCurrentText->SetText(cr->params.current_mA);
	idealBeamToggle->SetState(cr->params.emittance_cm==0.0);
	emittanceCouplingToggle->SetState(cr->params.emittance_cm != 0.0);

	if (cr->params.emittance_cm==0.0) { //ideal beam
		constantBXYtoggle->SetState(1);
		BXYfileNameText->SetText("");
		emittanceText->SetText("0");
		betaXtext->SetText("");
		betaYtext->SetText("");
		etaText->SetText("");
		etaPrimeText->SetText("");
		couplingText->SetText("");
		energySpreadText->SetText("");
		emittanceXtext->SetText("");
		emittanceYtext->SetText("");
	}
	
	emittanceText->SetText(cr->params.emittance_cm);
	couplingText->SetText(cr->params.coupling_percent);
	double emittanceX = cr->params.emittance_cm / (1.0 + cr->params.coupling_percent * 0.01);
	emittanceXtext->SetText(emittanceX);
	emittanceYtext->SetText(emittanceX * cr->params.coupling_percent * 0.01);
	energySpreadText->SetText(cr->params.energy_spread_percent);

	if (cr->params.betax_const_cm>=0.0) { //don't use BXY file
		constantBXYtoggle->SetState(1);
		BXYfileNameText->SetText("");
		betaXtext->SetText(cr->params.betax_const_cm);
		betaYtext->SetText(cr->params.betay_const_cm);
		etaText->SetText(cr->params.eta_x_const_cm);
		etaPrimeText->SetText(cr->params.eta_x_prime_const);
		energySpreadText->SetText(cr->params.energy_spread_percent);
	} else {
		constantBXYtoggle->SetState(0);
		BXYfileNameText->SetText(cr->BXYfileName);
	}

	EminText->SetText(cr->params.energy_low_eV);
	EmaxText->SetText(cr->params.energy_hi_eV);
	enableOrtPolarizationToggle->SetState(cr->params.enable_ort_polarization);
	enableParPolarizationToggle->SetState(cr->params.enable_par_polarization);
	limitAngleToggle->SetState(cr->params.psimaxX_rad<3.14159 || cr->params.psimaxY_rad<3.14159);
	psiMaxXtext->SetText(cr->params.psimaxX_rad);
	psiMaxYtext->SetText(cr->params.psimaxY_rad);

	BxtypeCombo->SetSelectedIndex(cr->params.Bx_mode-1);
	if (cr->params.Bx_mode==B_MODE_CONSTANT) {
		constBXtext->SetText(cr->params.B_const.x);
	} else {
		MAGfileXtext->SetText(cr->MAGXfileName);
	}
	BytypeCombo->SetSelectedIndex(cr->params.By_mode-1);
	if (cr->params.By_mode==B_MODE_CONSTANT) {
		constBYtext->SetText(cr->params.B_const.y);
	} else {
		MAGfileYtext->SetText(cr->MAGYfileName);
	}
	BztypeCombo->SetSelectedIndex(cr->params.Bz_mode-1);
	if (cr->params.Bz_mode==B_MODE_CONSTANT) {
		constBZtext->SetText(cr->params.B_const.z);
	} else {
		MAGfileZtext->SetText(cr->MAGZfileName);
	}
}

void RegionEditor::ApplyChanges() {
	double tmp;
	SynRad *mApp = (SynRad *)theApp;
	//First step:: error check
	if (!startPointXtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid X0","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (!startPointYtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid Y0","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (!startPointZtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid Z0","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (startDirDefinitionCombo->GetSelectedIndex()==0) {	//define by angle
		if (!theta0text->GetNumber(&tmp)) {GLMessageBox::Display("Invalid theta0","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (!alpha0text->GetNumber(&tmp)) {GLMessageBox::Display("Invalid alpha0","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	} else { //define by vector
		double x,y,z;
		if (!startDirXtext->GetNumber(&x)) {GLMessageBox::Display("Invalid dirX","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (!startDirYtext->GetNumber(&y)) {GLMessageBox::Display("Invalid dirY","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (!startDirZtext->GetNumber(&z)) {GLMessageBox::Display("Invalid dirZ","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		Vector3d dir=Vector3d(x,y,z);
		if (dir.Norme()==0.0) {GLMessageBox::Display("Start direction can't be a null-vector","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	}
	if (!dLtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid dL","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (!limitsXtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid Xmax","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (!limitsYtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid Ymax","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (!limitsZtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid Zmax","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}

	if (!particleMassText->GetNumber(&tmp)) {GLMessageBox::Display("Invalid particle mass","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (tmp<=0.0) {GLMessageBox::Display("Particle mass must be positive","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}

	if (!beamEnergyText->GetNumber(&tmp)) {GLMessageBox::Display("Invalid beam energy","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (tmp<=0.0) {GLMessageBox::Display("Beam energy must be positive","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}

	if (!beamCurrentText->GetNumber(&tmp)) {GLMessageBox::Display("Invalid beam current","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (tmp<=0.0) {GLMessageBox::Display("Beam current must be positive","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}

	if(!idealBeamToggle->GetState()) { //non-ideal beam
		if (!emittanceText->GetNumber(&tmp)) {GLMessageBox::Display("Invalid emittance","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (tmp<=0.0) {GLMessageBox::Display("Emittance must be positive","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (!couplingText->GetNumber(&tmp)) {GLMessageBox::Display("Invalid Coupling","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (tmp<0.0) {GLMessageBox::Display("Coupling must be non-negative","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (constantBXYtoggle->GetState()==1) { //constant values
			if (!betaXtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid BetaX","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (tmp<=0.0) {GLMessageBox::Display("BetaX must be positive","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (!betaYtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid BetaY","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (tmp<=0.0) {GLMessageBox::Display("BetaY must be positive","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (!etaText->GetNumber(&tmp)) {GLMessageBox::Display("Invalid Eta","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (tmp<0.0) {GLMessageBox::Display("Eta must be non-negative","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (!etaPrimeText->GetNumber(&tmp)) {GLMessageBox::Display("Invalid EtaPrime","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (tmp<0.0) {GLMessageBox::Display("EtaPrime must be non-negative","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (!energySpreadText->GetNumber(&tmp)) {GLMessageBox::Display("Invalid E_spread","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (tmp<0.0) {GLMessageBox::Display("E_Spread must be non-negative","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		} else { //use BXY file
			if (BXYfileNameText->GetTextLength()==0) {GLMessageBox::Display("BXY file name can't be empty","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (!FileUtils::Exist(BXYfileNameText->GetText())) {GLMessageBox::Display("BXY file doesn't exist","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		}
	}

	if (!EminText->GetNumber(&tmp)) {GLMessageBox::Display("Invalid Emin","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (tmp<=0.0) {GLMessageBox::Display("Emin must be positive","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	double tmp2;
	if (!EmaxText->GetNumber(&tmp2)) {GLMessageBox::Display("Invalid Emax","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (tmp2<=0.0) {GLMessageBox::Display("Emax must be positive","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	if (tmp>=tmp2) {GLMessageBox::Display("Emin must be smaller than Emax","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}

	if (limitAngleToggle->GetState()) {
		if (!psiMaxXtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid psiMaxX","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (tmp<=0.0) {GLMessageBox::Display("psiMaxX must be positive","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (!psiMaxYtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid psiMaxY","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (tmp<=0.0) {GLMessageBox::Display("psiMaxY must be positive","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	}
	if (BxtypeCombo->GetSelectedIndex()==0) {//const B field
		if (!constBXtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid Bx field value","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	} else { //use MAG file
		if (MAGfileXtext->GetTextLength()==0) {GLMessageBox::Display("X component: MAG file name can't be empty","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (!FileUtils::Exist(MAGfileXtext->GetText())) {GLMessageBox::Display("X component: the MAG file doesn't exist","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	}
	if (BytypeCombo->GetSelectedIndex()==0) {//const B field
		if (!constBYtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid By field value","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	} else { //use MAG file
		if (MAGfileYtext->GetTextLength()==0) {GLMessageBox::Display("Y component: MAG file name can't be empty","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (!FileUtils::Exist(MAGfileYtext->GetText())) {GLMessageBox::Display("Y component: the MAG file doesn't exist","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	}
	if (BztypeCombo->GetSelectedIndex()==0) {//const B field
		if (!constBZtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid Bz field value","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	} else { //use MAG file
		if (MAGfileZtext->GetTextLength()==0) {GLMessageBox::Display("Z component: MAG file name can't be empty","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
		if (!FileUtils::Exist(MAGfileZtext->GetText())) {GLMessageBox::Display("Z component: the MAG file doesn't exist","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
	}

	// Hallelujah, the user filled out the form without any errors
	// Let's apply his values
	if (!mApp->AskToReset()) return;
	//Deselect point, if has selected
	cr->selectedPoint=-1;

	startPointXtext->GetNumber(&cr->params.startPoint.x);
	startPointYtext->GetNumber(&cr->params.startPoint.y);
	startPointZtext->GetNumber(&cr->params.startPoint.z);

	if (startDirDefinitionCombo->GetSelectedIndex()==0) { //define by angle
		double a0,t0;
		alpha0text->GetNumber(&a0);
		theta0text->GetNumber(&t0);
		cr->params.startDir=Vector3d(sin(-t0)*cos(a0),sin(-a0),cos(a0)*cos(t0));  //left-handed coordinate system
	} else { // define by vector
		startDirXtext->GetNumber(&cr->params.startDir.x);
		startDirYtext->GetNumber(&cr->params.startDir.y);
		startDirZtext->GetNumber(&cr->params.startDir.z);
	}
	dLtext->GetNumber(&cr->params.dL_cm);
	limitsXtext->GetNumber(&cr->params.limits.x);
	limitsYtext->GetNumber(&cr->params.limits.y);
	limitsZtext->GetNumber(&cr->params.limits.z);
	particleMassText->GetNumber(&cr->params.particleMass_GeV);
	if (particleChargeCombo->GetSelectedIndex()==1) //negative charge
		cr->params.particleMass_GeV*=-1;
	beamEnergyText->GetNumber(&cr->params.E_GeV);
	beamCurrentText->GetNumber(&cr->params.current_mA);
	cr->params.gamma=abs(cr->params.E_GeV/cr->params.particleMass_GeV);
	if (idealBeamToggle->GetState()) cr->params.emittance_cm=0.0;
	else {
		emittanceText->GetNumber(&cr->params.emittance_cm);
		couplingText->GetNumber(&cr->params.coupling_percent);
		energySpreadText->GetNumber(&cr->params.energy_spread_percent);
		if (constantBXYtoggle->GetState()==1) { //constant values
			betaXtext->GetNumber(&cr->params.betax_const_cm);
			betaYtext->GetNumber(&cr->params.betay_const_cm);
			etaText->GetNumber(&cr->params.eta_x_const_cm);
			etaPrimeText->GetNumber(&cr->params.eta_x_prime_const);
		} else { //load BXY file
			cr->params.betax_const_cm=-1.0; //negative value: use BXY file
			cr->BXYfileName.assign(BXYfileNameText->GetText());
			//loadbxy
			try {
				cr->params.nbDistr_BXY=cr->LoadBXY(BXYfileNameText->GetText());
			} catch(Error &e) {
				char tmp[256];
				sprintf(tmp,"Couldn't load BXY file. Error message:\n%s",e.GetMsg());
				GLMessageBox::Display(tmp,"BXY file problem",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
		}
	}

	EminText->GetNumber(&cr->params.energy_low_eV);
	EmaxText->GetNumber(&cr->params.energy_hi_eV);
	cr->params.enable_ort_polarization=enableOrtPolarizationToggle->GetState();
	cr->params.enable_par_polarization=enableParPolarizationToggle->GetState();
	if (limitAngleToggle->GetState()) { //limit angle
		psiMaxXtext->GetNumber(&cr->params.psimaxX_rad);
		psiMaxYtext->GetNumber(&cr->params.psimaxY_rad);
	} else { 
		cr->params.psimaxX_rad=PI;
		cr->params.psimaxY_rad=PI;
	}

	//X comp.
	cr->params.Bx_mode=BxtypeCombo->GetSelectedIndex()+1;
	if (cr->params.Bx_mode==B_MODE_CONSTANT) { //constant B field
		constBXtext->GetNumber(&cr->params.B_const.x);
	} else { //use MAG file
		cr->MAGXfileName.assign(MAGfileXtext->GetText());
		try {
			FileReader MAGfile(MAGfileXtext->GetText());
			cr->Bx_distr=cr->LoadMAGFile(&MAGfile,&cr->params.Bx_dir,&cr->params.Bx_period,&cr->params.Bx_phase,cr->params.Bx_mode);
		} catch(Error &e) {
			char tmp[256];
			sprintf(tmp, "Couldn't load MAG file. Error message:\n%s", e.GetMsg());
			GLMessageBox::Display(tmp, "MAG file problem", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}

	//Y comp.
	cr->params.By_mode=BytypeCombo->GetSelectedIndex()+1;
	if (cr->params.By_mode==B_MODE_CONSTANT) { //constant B field
		constBYtext->GetNumber(&cr->params.B_const.y);
	} else { //use MAG file
		cr->MAGYfileName.assign(MAGfileYtext->GetText());
		try {
			FileReader MAGfile(MAGfileYtext->GetText());
			cr->By_distr=cr->LoadMAGFile(&MAGfile,&cr->params.By_dir,&cr->params.By_period,&cr->params.By_phase,cr->params.By_mode);
		} catch(Error &e) {
			char tmp[256];
			sprintf(tmp, "Couldn't load MAG file. Error message:\n%s", e.GetMsg());
			GLMessageBox::Display(tmp, "MAG file problem", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}

	//Z comp.
	cr->params.Bz_mode=BztypeCombo->GetSelectedIndex()+1;
	if (cr->params.Bz_mode==B_MODE_CONSTANT) { //constant B field
		constBZtext->GetNumber(&cr->params.B_const.z);
	} else { //use MAG file
		cr->MAGZfileName.assign(MAGfileZtext->GetText());
		try {
			FileReader MAGfile(MAGfileZtext->GetText());
			cr->Bz_distr=cr->LoadMAGFile(&MAGfile,&cr->params.Bz_dir,&cr->params.Bz_period,&cr->params.Bz_phase,cr->params.Bz_mode);
		} catch(Error &e) {
			char tmp[256];
			sprintf(tmp, "Couldn't load MAG file. Error message:\n%s", e.GetMsg());
			GLMessageBox::Display(tmp, "MAG file problem", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}

	worker->RecalcRegion(regionId);
	if (mApp->regionInfo) mApp->regionInfo->Update();
	if (mApp->trajectoryDetails && mApp->trajectoryDetails->IsVisible() && mApp->trajectoryDetails->GetRegionId() == regionId) mApp->trajectoryDetails->Update();
	if (mApp->spectrumPlotter) mApp->spectrumPlotter->SetScale();
	GLWindow::ProcessMessage(NULL,MSG_CLOSE);
}