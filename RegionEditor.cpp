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
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLInputBox.h"
#include "GLApp/GLTitledPanel.h"
#include "Utils.h"
#include "Synrad.h"

extern GLApplication *theApp;

// --------------------------------------------------------------------

RegionEditor::RegionEditor():GLWindow() {

	int wD=785;
	int hD=650;

	SetTitle("Region Editor");
	SetIconfiable(TRUE);
	SetMinimumSize(wD,hD);

	GLLabel*	label1 = new GLLabel("");
	GLLabel* label2 = new GLLabel("");
	startPointXtext = new GLTextField(0,"");
	startPointYtext = new GLTextField(0,"");
	GLLabel* label3 = new GLLabel("");
	startPointZtext = new GLTextField(0,"");
	GLLabel* label4 = new GLLabel("");
	GLLabel* label5 = new GLLabel("");
	startDirDefinitionCombo = new GLCombo(0);
	GLLabel* label6 = new GLLabel("");
	theta0text = new GLTextField(0,"");
	GLLabel* label7 = new GLLabel("");
	alpha0text = new GLTextField(0,"");
	GLLabel* label8 = new GLLabel("");
	startDirXtext = new GLTextField(0,"");
	GLLabel* label9 = new GLLabel("");
	startDirYtext = new GLTextField(0,"");
	GLLabel* label10 = new GLLabel("");
	startDirZtext = new GLTextField(0,"");
	GLLabel* label11 = new GLLabel("");
	GLLabel* label12 = new GLLabel("");
	dLtext = new GLTextField(0,"");
	GLLabel* label13 = new GLLabel("");
	GLLabel* label14 = new GLLabel("");
	GLLabel* label15 = new GLLabel("");
	GLLabel* label19 = new GLLabel("");
	GLLabel* label20 = new GLLabel("");
	GLLabel* label21 = new GLLabel("");
	GLLabel* label22 = new GLLabel("");
	GLLabel* label16 = new GLLabel("");
	limitsXtext = new GLTextField(0,"");
	GLLabel* label17 = new GLLabel("");
	GLLabel* label18 = new GLLabel("");
	limitsYtext = new GLTextField(0,"");
	GLLabel* label23 = new GLLabel("");
	GLLabel* label24 = new GLLabel("");
	limitsZtext = new GLTextField(0,"");
	GLLabel* label25 = new GLLabel("");
	GLLabel* label26 = new GLLabel("");
	particleMassText = new GLTextField(0,"");
	GLLabel* label27 = new GLLabel("");
	GLLabel* label28 = new GLLabel("");
	GLLabel* label29 = new GLLabel("");
	particleChargeCombo = new GLCombo(0);
	setParticleProtonButton = new GLButton(0,"");
	setParticleElectronButton = new GLButton(0,"");
	setParticlePositronButton = new GLButton(0,"");
	GLLabel* label30 = new GLLabel("");
	GLLabel* label301 = new GLLabel("Beam energy:");
	beamEnergyText = new GLTextField(0,"");
	GLLabel *currentLabel=new GLLabel("Current:");
	beamCurrentText = new GLTextField(0,"");
	GLLabel *mALabel=new GLLabel("mA");
	GLLabel* label302 = new GLLabel("GeV");
	idealBeamToggle = new GLToggle(0,"Ideal Beam");
	betaOrBXYCombo = new GLCombo(0);
	BXYfileNameText = new GLTextField(0,"");
	GLLabel* label31 = new GLLabel("");
	GLLabel* label32 = new GLLabel("");
	betaXtext = new GLTextField(0,"");
	betaYtext = new GLTextField(0,"");
	GLLabel* label34 = new GLLabel("");
	etaText = new GLTextField(0,"");
	emittanceText = new GLTextField(0,"");
	GLLabel* emittance_cm = new GLLabel("cm");
	GLLabel* betax_cm = new GLLabel("cm");
	GLLabel* betay_cm = new GLLabel("cm");
	GLLabel* coupling_percent = new GLLabel("%");
	GLLabel* eta_cm = new GLLabel("cm");
	GLLabel* etaprime_rad = new GLLabel("rad");
	GLLabel* E_spread_percent = new GLLabel("%");

	GLLabel* label35 = new GLLabel("");
	GLLabel* label36 = new GLLabel("");
	GLLabel* label37 = new GLLabel("");
	etaPrimeText = new GLTextField(0,"");
	GLLabel* label33 = new GLLabel("");
	EminText = new GLTextField(0,"");
	GLLabel* label38 = new GLLabel("");
	GLLabel* label39 = new GLLabel("");
	EmaxText = new GLTextField(0,"");
	GLLabel* label40 = new GLLabel("");
	GLLabel* label41 = new GLLabel("");
	GLLabel* label42 = new GLLabel("");
	enableParPolarizationToggle = new GLToggle(0,"Parallel");
	enableOrtPolarizationToggle = new GLToggle(0,"Orthogonal");
	limitAngleToggle = new GLToggle(0,"Limit photon divergence:");
	GLLabel* label43 = new GLLabel("");
	psiMaxXtext = new GLTextField(0,"");
	GLLabel* label44 = new GLLabel("");
	GLLabel* label45 = new GLLabel("");
	psiMaxYtext = new GLTextField(0,"");
	GLLabel* label46 = new GLLabel("");
	GLLabel* label47 = new GLLabel("");
	BxtypeCombo = new GLCombo(0);
	constBXtext = new GLTextField(0,"");
	GLLabel* label48 = new GLLabel("");
	GLLabel* label49 = new GLLabel("");
	MAGfileXtext = new GLTextField(0,"");
	GLLabel* label50 = new GLLabel("");
	BytypeCombo = new GLCombo(0);
	constBYtext = new GLTextField(0,"");
	GLLabel* label51 = new GLLabel("");
	GLLabel* label52 = new GLLabel("");
	MAGfileYtext = new GLTextField(0,"");
	GLLabel* label53 = new GLLabel("");
	BztypeCombo = new GLCombo(0);
	constBZtext = new GLTextField(0,"");
	GLLabel* label54 = new GLLabel("");
	GLLabel* label55 = new GLLabel("");
	MAGfileZtext = new GLTextField(0,"");
	applyButton = new GLButton(0,"");
	dismissButton = new GLButton(0,"");
	limitsInfoButton = new GLButton(0,"Info");
	dlInfoButton = new GLButton(0,"Info");
	startDirInfoButton = new GLButton(0,"Info");
	beamsizeInfoButton = new GLButton(0,"Info");
	bxyBrowseButton = new GLButton(0,"...");
	magxBrowseButton = new GLButton(0,"...");
	magyBrowseButton = new GLButton(0,"...");
	magzBrowseButton = new GLButton(0,"...");
	bxyEditButton = new GLButton(0,"Edit");
	magxEditButton = new GLButton(0,"Edit");
	magyEditButton = new GLButton(0,"Edit");
	magzEditButton = new GLButton(0,"Edit");
	GLTitledPanel* geomPanel = new GLTitledPanel("Beam geometry");
	GLTitledPanel* particlePanel = new GLTitledPanel("Particle settings");
	GLTitledPanel* beamPanel = new GLTitledPanel("Additional beam parameters");
	GLTitledPanel* generationPanel = new GLTitledPanel("Photon generation");
	GLTitledPanel* magPanel = new GLTitledPanel("Magnetic field");

	geomPanel->SetBounds(10,5,wD-20,150);
	Add(geomPanel);

	//---- label1 ----
	label1->SetText("Beam start position:");
	Add(label1);
	label1->SetBounds(30, 25, 125, 18);

	//---- label2 ----
	label2->SetText("X0:");
	Add(label2);
	label2->SetBounds(180, 25, 23,18);
	Add(startPointXtext);
	startPointXtext->SetBounds(240, 24, 80, 18);

	//---- label3 ----
	label3->SetText("Y0:");
	Add(label3);
	label3->SetBounds(365, 25, 23, 18);
	Add(startPointYtext);
	startPointYtext->SetBounds(425, 24, 80, 18);

	//---- label4 ----
	label4->SetText("Z0:");
	Add(label4);
	label4->SetBounds(565, 25, 23, 18);
	Add(startPointZtext);
	startPointZtext->SetBounds(615, 24, 80, 18);
	
	//---- label5 ----
	label5->SetText("Beam start direction:");
	Add(label5);
	label5->SetBounds(30, 50, 140, 18);
	Add(startDirDefinitionCombo);
	startDirDefinitionCombo->SetBounds(240, 48, 80, 19);
	startDirDefinitionCombo->SetSize(2);
	startDirDefinitionCombo->SetValueAt(0,"By angle:");
	startDirDefinitionCombo->SetValueAt(1,"By vector:");
	startDirDefinitionCombo->SetSelectedIndex(1); //Defining by direction vector is always well-defined

	//---- label6 ----
	label6->SetText("Theta0:");
	Add(label6);
	label6->SetBounds(365, 50, 50, 18);
	Add(theta0text);
	theta0text->SetBounds(425, 49, 80, 18);

	//---- label7 ----
	label7->SetText("Alpha0:");
	Add(label7);
	label7->SetBounds(565, 50, 53, 18);
	Add(alpha0text);
	alpha0text->SetBounds(615, 49, 80, 18);

	//---- startDirInfoButton ----
	Add(startDirInfoButton);
	startDirInfoButton->SetBounds(725, 49, 40, 18);

	//---- label8 ----
	label8->SetText("dirX:");
	Add(label8);
	label8->SetBounds(180, 75, 30, 18);

	//---- startDirXtext ----
	startDirXtext->SetEditable(FALSE);
	Add(startDirXtext);
	startDirXtext->SetBounds(240, 74, 80, 18);

	//---- label9 ----
	label9->SetText("dirY:");
	Add(label9);
	label9->SetBounds(365, 75, 33, 18);

	//---- startDirYtext ----
	startDirYtext->SetEditable(FALSE);
	Add(startDirYtext);
	startDirYtext->SetBounds(425, 74, 80, 18);

	//---- label10 ----
	label10->SetText("dirZ:");
	Add(label10);
	label10->SetBounds(565, 75, 28, 18);

	//---- startDirZtext ----
	startDirZtext->SetEditable(FALSE);
	Add(startDirZtext);
	startDirZtext->SetBounds(615, 74, 80, 18);

	//---- label11 ----
	label11->SetText("Trajectory step length:");
	Add(label11);
	label11->SetBounds(30, 100, 135, 18);

	//---- label12 ----
	label12->SetText("dL:");
	Add(label12);
	label12->SetBounds(180, 100, 18, 18);
	Add(dLtext);
	dLtext->SetBounds(240, 99, 80, 18);

	//---- dlInfoButton ----
	Add(dlInfoButton);
	dlInfoButton->SetBounds(725, 99, 40, 18);

	//---- label13 ----
	label13->SetText("cm");
	Add(label13);
	label13->SetBounds(325, 25, 17,18);

	//---- label14 ----
	label14->SetText("cm");
	Add(label14);
	label14->SetBounds(510, 25, 17, 18);

	//---- label15 ----
	label15->SetText("cm");
	Add(label15);
	label15->SetBounds(700, 25, 17, 18);

	//---- label19 ----
	label19->SetText("cm");
	Add(label19);
	label19->SetBounds(325, 100, 17, 18);

	//---- label20 ----
	label20->SetText("rad");
	Add(label20);
	label20->SetBounds(510, 50, 25, 18);

	//---- label21 ----
	label21->SetText("rad");
	Add(label21);
	label21->SetBounds(700, 50, 25, 18);

	//---- label22 ----
	label22->SetText("Calculation boundary:");
	Add(label22);
	label22->SetBounds(30, 125, 135, 18);

	//---- label16 ----
	label16->SetText("Xmax:");
	Add(label16);
	label16->SetBounds(180, 125, 40, 18);
	Add(limitsXtext);
	limitsXtext->SetBounds(240, 124, 80, 18);

	//---- label17 ----
	label17->SetText("cm");
	Add(label17);
	label17->SetBounds(325, 125, 17, 18);

	//---- label18 ----
	label18->SetText("Ymax:");
	Add(label18);
	label18->SetBounds(365, 125, 38, 18);
	Add(limitsYtext);
	limitsYtext->SetBounds(425, 124, 80, 18);

	//---- label23 ----
	label23->SetText("cm");
	Add(label23);
	label23->SetBounds(510, 125, 17, 18);

	//---- label24 ----
	label24->SetText("Zmax:");
	Add(label24);
	label24->SetBounds(565, 125, 43, 18);
	Add(limitsZtext);
	limitsZtext->SetBounds(615, 124, 80, 18);

	//---- label25 ----
	label25->SetText("cm");
	Add(label25);
	label25->SetBounds(700, 125, 17, 18);

	//---- limitsInfoButton ----
	Add(limitsInfoButton);
	limitsInfoButton->SetBounds(725, 124, 40, 18);

	particlePanel->SetBounds(10,160,wD-20,70);
	Add(particlePanel);

	//---- label26 ----
	label26->SetText("Particle mass:");
	Add(label26);
	label26->SetBounds(30, 175, 135, 18);
	Add(particleMassText);
	particleMassText->SetBounds(240, 174, 80, 18);

	//---- label27 ----
	label27->SetText("pmass:");
	Add(label27);
	label27->SetBounds(180, 175, 45, 18);

	//---- label28 ----
	label28->SetText("GeV/c\262");
	Add(label28);
	label28->SetBounds(325, 175, 17, 18);

	//---- label29 ----
	label29->SetText("Particle charge:");
	Add(label29);
	label29->SetBounds(30, 200, 135, 18);
	Add(particleChargeCombo);
	particleChargeCombo->SetBounds(240, 199, 80, 19);
	particleChargeCombo->SetSize(2);
	particleChargeCombo->SetValueAt(0,"Positive");
	particleChargeCombo->SetValueAt(1,"Negative");		
	particleChargeCombo->SetSelectedIndex(0);

	//---- label30 ----
	label30->SetText("Quick set to:");
	Add(label30);
	label30->SetBounds(425, 175, 60, 18);

	//---- setParticleElectronButton ----
	setParticleElectronButton->SetText("Electron");
	Add(setParticleElectronButton);
	setParticleElectronButton->SetBounds(425, 200, 80, 18);

	//---- setParticlePositronButton ----
	setParticlePositronButton->SetText("Positron");
	Add(setParticlePositronButton);
	setParticlePositronButton->SetBounds(515, 200, 80, 18);

	//---- setParticleProtonButton ----
	setParticleProtonButton->SetText("Proton");
	Add(setParticleProtonButton);
	setParticleProtonButton->SetBounds(605, 200, 80, 18);

	beamPanel->SetBounds(10,235,wD-20,150);
	Add(beamPanel);

	label301->SetBounds(165,250,100,18);//Beam energy:
	Add(label301);

	beamEnergyText->SetBounds(240,250,80,18);
	Add(beamEnergyText);

	label302->SetBounds(325,250,20,18);//GeV
	Add(label302);

	currentLabel->SetBounds(365,250,50,18);
	Add(currentLabel);

	beamCurrentText->SetBounds(425,250,80,18);
	Add(beamCurrentText);

	mALabel->SetBounds(510,250,30,18);
	Add(mALabel);

	//---- idealBeamToggle ----
	//idealBeamToggle->SetText("Ideal beam");
	idealBeamToggle->SetState(TRUE);
	Add(idealBeamToggle);
	idealBeamToggle->SetBounds(30, 285, 160, 18);

	//---- betaOrBXYCombo ----
	betaOrBXYCombo->SetEditable(FALSE);
	betaOrBXYCombo->SetEnabled(FALSE);
	Add(betaOrBXYCombo);
	betaOrBXYCombo->SetBounds(180, 284, 140, 19);
	betaOrBXYCombo->SetSize(2);
	betaOrBXYCombo->SetValueAt(0,"Constant values");
	betaOrBXYCombo->SetValueAt(1,"From BXY file");
	betaOrBXYCombo->SetSelectedIndex(0);

	//---- label31 ----
	label31->SetText("BXY file:");
	Add(label31);
	label31->SetBounds(365, 285, 48, 18);

	//---- BXYfileNameText ----
	BXYfileNameText->SetEditable(FALSE);
	Add(BXYfileNameText);
	BXYfileNameText->SetBounds(425, 284, 260, 18);
	
	//---- label56 ----
	Add(bxyBrowseButton);
	bxyBrowseButton->SetBounds(688, 284, 35, 18);

	bxyEditButton->SetBounds(727,284,40,18);
	Add(bxyEditButton);

	//---- label35 ----
	label35->SetText("Emittance:");
	Add(label35);
	label35->SetBounds(180, 309, 65, 18);
	
	//---- emittanceText ----
	emittanceText->SetEditable(FALSE);
	Add(emittanceText);
	emittanceText->SetBounds(240, 309, 80, 18);

	Add(emittance_cm);
	emittance_cm->SetBounds(325, 309, 17, 18);

	//---- label32 ----
	label32->SetText("BetaX:");
	Add(label32);
	label32->SetBounds(365, 310, 40, 18);

	//---- betaXtext ----
	betaXtext->SetEditable(FALSE);
	Add(betaXtext);
	betaXtext->SetBounds(425, 309, 80, 18);

	Add(betax_cm);
	betax_cm->SetBounds(510, 309, 17, 18);

	//---- label36 ----
	label36->SetText("BetaY:");
	Add(label36);
	label36->SetBounds(555, 310, 40, 18);

	//---- betaYtext ----
	betaYtext->SetEditable(FALSE);
	Add(betaYtext);
	betaYtext->SetBounds(605, 309, 80, 18);

	Add(betay_cm);
	betay_cm->SetBounds(690, 309, 17, 18);

	//Couplinglabel
	GLLabel *couplingLabel = new GLLabel("Coupling:");
	couplingLabel->SetBounds(180,335,60,18);
	Add(couplingLabel);

	//couplingText
	couplingText = new GLTextField(0,"");
	couplingText->SetEditable(FALSE);
	couplingText->SetBounds(240,334,80,18);
	Add(couplingText);

	Add(coupling_percent);
	coupling_percent->SetBounds(325,335,17,18);

	//---- label34 ----
	label34->SetText("Eta:");
	Add(label34);
	label34->SetBounds(365, 335, 23, 18);

	//---- etaText ----
	etaText->SetEditable(FALSE);
	Add(etaText);
	etaText->SetBounds(425, 334, 80, 18);

	Add(eta_cm);
	eta_cm->SetBounds(510,335,17,18);

	//---- label37 ----
	label37->SetText("EtaPrime:");
	Add(label37);
	label37->SetBounds(555, 335, 58, 18);

	//---- etaPrimeText ----
	etaPrimeText->SetEditable(FALSE);
	Add(etaPrimeText);
	etaPrimeText->SetBounds(605, 334, 80, 18);

	Add(etaprime_rad);
	etaprime_rad->SetBounds(690,335,17,18);

	//---- beamsizeInfoButton ----
	Add(beamsizeInfoButton);
	beamsizeInfoButton->SetBounds(240, 359, 80, 18);

	//energySpreadLabel
	GLLabel *energySpreadLabel = new GLLabel("E_Spread:");
	energySpreadLabel->SetBounds(365,360,60,18);
	Add(energySpreadLabel);

	//energySpreadText
	energySpreadText = new GLTextField(0,"");
	energySpreadText->SetEditable(FALSE);
	energySpreadText->SetBounds(425,359,80,18);
	Add(energySpreadText);

	Add(E_spread_percent);
	E_spread_percent->SetBounds(510,360,17,18);

	generationPanel->SetBounds(10,390,wD-20,100);
	Add(generationPanel);

	//---- label41 ----
	label41->SetText("Generated photon energy:");
	Add(label41);
	label41->SetBounds(30, 410, 165, 18);

	//---- label33 ----
	label33->SetText("Emin:");
	Add(label33);
	label33->SetBounds(180, 410, 35, 18);
	
	Add(EminText);
	EminText->SetBounds(240, 409, 80, 18);

	//---- label38 ----
	label38->SetText("eV");
	Add(label38);
	label38->SetBounds(325, 410, 17, 18);

	//---- label39 ----
	label39->SetText("Emax:");
	Add(label39);
	label39->SetBounds(365, 410, 38, 18);
	
	Add(EmaxText);
	EmaxText->SetBounds(425, 409, 80, 18);

	//---- label40 ----
	label40->SetText("eV");
	Add(label40);
	label40->SetBounds(510, 410, 17, 18);

	//---- label42 ----
	label42->SetText("Generated polarization modes:");
	Add(label42);
	label42->SetBounds(30, 435, 180, 18);

	enableParPolarizationToggle->SetState(TRUE);
	Add(enableParPolarizationToggle);
	enableParPolarizationToggle->SetBounds(240, 434, 100, 18);

	enableOrtPolarizationToggle->SetState(TRUE);
	Add(enableOrtPolarizationToggle);
	enableOrtPolarizationToggle->SetBounds(425, 434, 100, 18);

	Add(limitAngleToggle);
	limitAngleToggle->SetBounds(30, 460, 120,18);

	//---- label43 ----
	label43->SetText("psiMaxX:");
	Add(label43);
	label43->SetBounds(365, 460, 55, 18);

	//---- psiMaxXtext ----
	psiMaxXtext->SetEditable(FALSE);
	Add(psiMaxXtext);
	psiMaxXtext->SetBounds(425, 459, 80, 18);

	//---- label44 ----
	label44->SetText("Rad");
	Add(label44);
	label44->SetBounds(510, 460, 40, 18);

	//---- label45 ----
	label45->SetText("psiMaxY:");
	Add(label45);
	label45->SetBounds(560, 460, 55, 18);

	//---- psiMaxYtext ----
	psiMaxYtext->SetEditable(FALSE);
	Add(psiMaxYtext);
	psiMaxYtext->SetBounds(615, 459, 80, 18);

	//---- label46 ----
	label46->SetText("Rad");
	Add(label46);
	label46->SetBounds(700, 460, 40, 18);

	magPanel->SetBounds(10,495,wD-20,100);
	Add(magPanel);

	//---- label47 ----
	label47->SetText("Bx:");
	Add(label47);
	label47->SetBounds(30, 515, 25, 20);
	Add(BxtypeCombo);
	BxtypeCombo->SetBounds(60, 514, 150, 19);
	BxtypeCombo->SetSize(8);
	BxtypeCombo->SetValueAt(0,"Constant field");
	BxtypeCombo->SetValueAt(1,"Coords along a direction");
	BxtypeCombo->SetValueAt(2,"Coords along the beam");
	BxtypeCombo->SetValueAt(3,"Sine / cosine");
	BxtypeCombo->SetValueAt(4,"Quadrupole");
	BxtypeCombo->SetValueAt(5,"Analytical expression");
	BxtypeCombo->SetValueAt(6,"Helicoidal");
	BxtypeCombo->SetValueAt(7,"Rotating dipole");
	BxtypeCombo->SetSelectedIndex(0);
	Add(constBXtext);
	constBXtext->SetBounds(240, 514, 80, 18);

	//---- label48 ----
	label48->SetText("T");
	Add(label48);
	label48->SetBounds(325, 515, 17, 18);

	//---- label49 ----
	label49->SetText("MAG file:");
	Add(label49);
	label49->SetBounds(365, 515, 55, 18);

	//---- MAGfileXtext ----
	MAGfileXtext->SetEditable(FALSE);
	Add(MAGfileXtext);
	MAGfileXtext->SetBounds(425, 514, 260, 18);

		//---- label57 ----
	Add(magxBrowseButton);
	magxBrowseButton->SetBounds(688, 514, 35, 18);

	magxEditButton->SetBounds(727,514,40,18);
	Add(magxEditButton);

	//---- label50 ----
	label50->SetText("By:");
	Add(label50);
	label50->SetBounds(30, 540, 25, 18);
	Add(BytypeCombo);
	BytypeCombo->SetBounds(60, 539, 150, 19);
	BytypeCombo->SetSize(8);
	BytypeCombo->SetValueAt(0, "Constant field");
	BytypeCombo->SetValueAt(1, "Coords along a direction");
	BytypeCombo->SetValueAt(2, "Coords along the beam");
	BytypeCombo->SetValueAt(3, "Sine / cosine");
	BytypeCombo->SetValueAt(4, "Quadrupole");
	BytypeCombo->SetValueAt(5, "Analytical expression");
	BytypeCombo->SetValueAt(6, "Helicoidal");
	BytypeCombo->SetValueAt(7, "Rotating dipole");
	BytypeCombo->SetSelectedIndex(0);
	Add(constBYtext);
	constBYtext->SetBounds(240, 539, 80, 18);

	//---- label51 ----
	label51->SetText("T");
	Add(label51);
	label51->SetBounds(325, 540, 17, 18);

	//---- label52 ----
	label52->SetText("MAG file:");
	Add(label52);
	label52->SetBounds(365, 540, 55, 18);

	//---- MAGfileYtext ----
	MAGfileYtext->SetEditable(FALSE);
	Add(MAGfileYtext);
	MAGfileYtext->SetBounds(425, 539, 260, 18);

	//---- label58 ----
	Add(magyBrowseButton);
	magyBrowseButton->SetBounds(688, 539, 35, 18);

	magyEditButton->SetBounds(727,539,40,18);
	Add(magyEditButton);

	//---- label53 ----
	label53->SetText("Bz:");
	Add(label53);
	label53->SetBounds(30, 565, 25, 18);
	Add(BztypeCombo);
	BztypeCombo->SetBounds(60, 564, 150, 19);
	BztypeCombo->SetSize(8);
	BztypeCombo->SetValueAt(0, "Constant field");
	BztypeCombo->SetValueAt(1, "Coords along a direction");
	BztypeCombo->SetValueAt(2, "Coords along the beam");
	BztypeCombo->SetValueAt(3, "Sine / cosine");
	BztypeCombo->SetValueAt(4, "Quadrupole");
	BztypeCombo->SetValueAt(5, "Analytical expression");
	BztypeCombo->SetValueAt(6, "Helicoidal");
	BztypeCombo->SetValueAt(7, "Rotating dipole");
	BztypeCombo->SetSelectedIndex(0);
	Add(constBZtext);
	constBZtext->SetBounds(240, 564, 80, 18);

	//---- label54 ----
	label54->SetText("T");
	Add(label54);
	label54->SetBounds(325, 565, 17, 18);

	//---- label55 ----
	label55->SetText("MAG file:");
	Add(label55);
	label55->SetBounds(365, 565, 55, 18);

	//---- MAGfileZtext ----
	MAGfileZtext->SetEditable(FALSE);
	Add(MAGfileZtext);
	MAGfileZtext->SetBounds(425, 564, 260, 18);

	//---- label59 ----
	Add(magzBrowseButton);
	magzBrowseButton->SetBounds(688, 564, 35, 18);

	magzEditButton->SetBounds(727,564,40,18);
	Add(magzEditButton);

	//---- applyButton ----
	applyButton->SetText("Apply & Recalculate points");
	Add(applyButton);
	applyButton->SetBounds(425, 600, 225, 21);

	//---- dismissButton ----
	dismissButton->SetText("Dismiss");
	Add(dismissButton);
	dismissButton->SetBounds(655, 600, 120, 21);

	// Center dialog
	int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;
	SetBounds(xD,yD,wD,hD);

	RestoreDeviceObjects();
}


// --------------------------------------------------------------------

void RegionEditor::Display(Worker *w,int Id) {

	//set textfields and values here

	worker = w;
	this->regionId=Id;
	cr = &(worker->regions[regionId]);

	FillValues();
	EnableDisableControls();

	SetVisible(TRUE);
}


void RegionEditor::ProcessMessage(GLComponent *src,int message) {
	SynRad *mApp = (SynRad *)theApp;
	switch(message) {
	case MSG_BUTTON:

		if(src==dismissButton) {

			GLWindow::ProcessMessage(NULL,MSG_CLOSE);

		} else if (src==setParticleElectronButton) {
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
			StartProc_foreground(tmp);
		} else if (src==this->magxEditButton) {
			char tmp[512];
			sprintf(tmp,"notepad.exe \"%s\"",MAGfileXtext->GetText());
			StartProc_foreground(tmp);
		} else if (src==this->magyEditButton) {
			char tmp[512];
			sprintf(tmp,"notepad.exe \"%s\"",MAGfileYtext->GetText());
			StartProc_foreground(tmp);
		} else if (src==this->magzEditButton) {
			char tmp[512];
			sprintf(tmp,"notepad.exe \"%s\"",MAGfileZtext->GetText());
			StartProc_foreground(tmp);
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
				"sigma_x_prime = sqrt( horizontal_emittance / beta_X + (eta_prime / energy_spread)^2 )\n"
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
		EnableDisableControls();
		break;

	}
	GLWindow::ProcessMessage(src,message);
}
void RegionEditor::EnableDisableControls() {
	BOOL defineByAngle=(startDirDefinitionCombo->GetSelectedIndex()==0);
	theta0text->SetEditable(defineByAngle);
	alpha0text->SetEditable(defineByAngle);
	startDirXtext->SetEditable(!defineByAngle);
	startDirYtext->SetEditable(!defineByAngle);
	startDirZtext->SetEditable(!defineByAngle);
	BOOL useMAGfile=(BxtypeCombo->GetSelectedIndex()!=0);
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
	magzEditButton->SetEnabled(useMAGfile);
	BOOL useIdealBeam=idealBeamToggle->GetState();
	BOOL useBXYfile=(betaOrBXYCombo->GetSelectedIndex()==1);
	betaOrBXYCombo->SetEditable(!useIdealBeam);
	betaOrBXYCombo->SetEnabled(!useIdealBeam);
	BXYfileNameText->SetEditable(!useIdealBeam && useBXYfile);
	bxyBrowseButton->SetEnabled(!useIdealBeam && useBXYfile);
	bxyEditButton->SetEnabled(!useIdealBeam && useBXYfile);
	emittanceText->SetEditable(!useIdealBeam);
	couplingText->SetEditable(!useIdealBeam);
	betaXtext->SetEditable(!useIdealBeam && !useBXYfile);
	betaYtext->SetEditable(!useIdealBeam && !useBXYfile);
	etaText->SetEditable(!useIdealBeam && !useBXYfile);
	etaPrimeText->SetEditable(!useIdealBeam && !useBXYfile);			
	energySpreadText->SetEditable(!useIdealBeam && !useBXYfile);			
	BOOL limitAngle=limitAngleToggle->GetState();
	psiMaxXtext->SetEditable(limitAngle);
	psiMaxYtext->SetEditable(limitAngle);
}

void RegionEditor::FillValues() {
	char tmp[512];
	sprintf(tmp,"%.9g",cr->startPoint.x);startPointXtext->SetText(tmp);
	sprintf(tmp,"%.9g",cr->startPoint.y);startPointYtext->SetText(tmp);
	sprintf(tmp,"%.9g",cr->startPoint.z);startPointZtext->SetText(tmp);

	double t0;
	if (cr->startDir.z == 0) {
		if (cr->startDir.x >= 0) t0 = -PI / 2;
		else t0 = PI / 2;
	}
	else {
		t0 = -atan(cr->startDir.x / cr->startDir.z); //Good for -PI/2...PI/2
		if (cr->startDir.z <= 0) //atan out of period
			if (cr->startDir.x < 0) t0 += PI;
			else t0 -= PI;
	}
	theta0text->SetText(t0);
	alpha0text->SetText(-asin(cr->startDir.Normalize().y));
	sprintf(tmp,"%g",cr->startDir.x);startDirXtext->SetText(tmp);
	sprintf(tmp,"%g",cr->startDir.y);startDirYtext->SetText(tmp);
	sprintf(tmp,"%g",cr->startDir.z);startDirZtext->SetText(tmp);

	sprintf(tmp,"%g",cr->dL);dLtext->SetText(tmp);

	sprintf(tmp,"%g",cr->limits.x);limitsXtext->SetText(tmp);
	sprintf(tmp,"%g",cr->limits.y);limitsYtext->SetText(tmp);
	sprintf(tmp,"%g",cr->limits.z);limitsZtext->SetText(tmp);

	sprintf(tmp,"%g",abs(cr->particleMass));particleMassText->SetText(tmp);
	particleChargeCombo->SetSelectedIndex(cr->particleMass<0);

	sprintf(tmp,"%g",cr->E);beamEnergyText->SetText(tmp);
	sprintf(tmp,"%g",cr->current);beamCurrentText->SetText(tmp);
	idealBeamToggle->SetState(cr->emittance==0.0);

	if (cr->emittance==0.0) { //ideal beam
		betaOrBXYCombo->SetSelectedIndex(0);
		BXYfileNameText->SetText("");
		emittanceText->SetText("0");
		betaXtext->SetText("0");
		betaYtext->SetText("0");
		etaText->SetText("0");
		etaPrimeText->SetText("0");
		couplingText->SetText("0");
		energySpreadText->SetText("0");
	}
	
	sprintf(tmp,"%g",cr->emittance);emittanceText->SetText(tmp);
	sprintf(tmp,"%g",cr->coupling);couplingText->SetText(tmp);

	if (cr->betax>=0.0) { //don't use BXY file
		betaOrBXYCombo->SetSelectedIndex(0);
		BXYfileNameText->SetText("");
		sprintf(tmp,"%g",cr->betax);betaXtext->SetText(tmp);
		sprintf(tmp,"%g",cr->betay);betaYtext->SetText(tmp);
		sprintf(tmp,"%g",cr->eta);etaText->SetText(tmp);
		sprintf(tmp,"%g",cr->etaprime);etaPrimeText->SetText(tmp);
		sprintf(tmp,"%g",cr->energy_spread);energySpreadText->SetText(tmp);
	} else {
		betaOrBXYCombo->SetSelectedIndex(1);
		BXYfileNameText->SetText((char*)cr->BXYfileName.c_str());
	}

	sprintf(tmp,"%g",cr->energy_low);EminText->SetText(tmp);
	sprintf(tmp,"%g",cr->energy_hi);EmaxText->SetText(tmp);
	enableOrtPolarizationToggle->SetState(cr->enable_ort_polarization);
	enableParPolarizationToggle->SetState(cr->enable_par_polarization);
	limitAngleToggle->SetState(cr->psimaxX<3.14159 || cr->psimaxY<3.14159);
	sprintf(tmp,"%g",cr->psimaxX);psiMaxXtext->SetText(tmp);
	sprintf(tmp,"%g",cr->psimaxY);psiMaxYtext->SetText(tmp);

	BxtypeCombo->SetSelectedIndex(cr->Bx_mode-1);
	if (cr->Bx_mode==B_MODE_CONSTANT) {
		sprintf(tmp,"%g",cr->B_const.x);constBXtext->SetText(tmp);
	} else {
		MAGfileXtext->SetText((char*)cr->MAGXfileName.c_str());
	}
	BytypeCombo->SetSelectedIndex(cr->By_mode-1);
	if (cr->By_mode==B_MODE_CONSTANT) {
		sprintf(tmp,"%g",cr->B_const.y);constBYtext->SetText(tmp);
	} else {
		MAGfileYtext->SetText((char*)cr->MAGYfileName.c_str());
	}
	BztypeCombo->SetSelectedIndex(cr->Bz_mode-1);
	if (cr->Bz_mode==B_MODE_CONSTANT) {
		sprintf(tmp,"%g",cr->B_const.z);constBZtext->SetText(tmp);
	} else {
		MAGfileZtext->SetText((char*)cr->MAGZfileName.c_str());
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
		Vector dir=Vector(x,y,z);
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
		if (betaOrBXYCombo->GetSelectedIndex()==0) { //constant values
			if (!betaXtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid BetaX","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (tmp<0.0) {GLMessageBox::Display("BetaX must be non-negative","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (!betaYtext->GetNumber(&tmp)) {GLMessageBox::Display("Invalid BetaY","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
			if (tmp<0.0) {GLMessageBox::Display("BetaY must be non-negative","Invalid input",GLDLG_OK,GLDLG_ICONERROR);return;}
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

	startPointXtext->GetNumber(&cr->startPoint.x);
	startPointYtext->GetNumber(&cr->startPoint.y);
	startPointZtext->GetNumber(&cr->startPoint.z);

	if (startDirDefinitionCombo->GetSelectedIndex()==0) { //define by angle
		double a0,t0;
		alpha0text->GetNumber(&a0);
		theta0text->GetNumber(&t0);
		cr->startDir=Vector(sin(-t0)*cos(a0),sin(-a0),cos(a0)*cos(t0));  //left-handed coordinate system
	} else { // define by vector
		startDirXtext->GetNumber(&cr->startDir.x);
		startDirYtext->GetNumber(&cr->startDir.y);
		startDirZtext->GetNumber(&cr->startDir.z);
	}
	dLtext->GetNumber(&cr->dL);
	limitsXtext->GetNumber(&cr->limits.x);
	limitsYtext->GetNumber(&cr->limits.y);
	limitsZtext->GetNumber(&cr->limits.z);
	particleMassText->GetNumber(&cr->particleMass);
	if (particleChargeCombo->GetSelectedIndex()==1) //negative charge
		cr->particleMass*=-1;
	beamEnergyText->GetNumber(&cr->E);
	beamCurrentText->GetNumber(&cr->current);
	cr->gamma=abs(cr->E/cr->particleMass);
	if (idealBeamToggle->GetState()) cr->emittance=0.0;
	else {
		emittanceText->GetNumber(&cr->emittance);
		couplingText->GetNumber(&cr->coupling);
		if (betaOrBXYCombo->GetSelectedIndex()==0) { //constant values
			betaXtext->GetNumber(&cr->betax);
			betaYtext->GetNumber(&cr->betay);
			etaText->GetNumber(&cr->eta);
			etaPrimeText->GetNumber(&cr->etaprime);	
			energySpreadText->GetNumber(&cr->energy_spread);
		} else { //load BXY file
			cr->betax=-1.0; //negative value: use BXY file
			cr->BXYfileName.assign(BXYfileNameText->GetText());
			//loadbxy
			try {
				FileReader BXYfile(BXYfileNameText->GetText());
				cr->nbDistr_BXY=cr->LoadBXY(&BXYfile,cr->beta_x_distr,cr->beta_y_distr,
					cr->eta_distr,cr->etaprime_distr,cr->e_spread_distr);
			} catch(Error &e) {
				char tmp[256];
				sprintf(tmp,"Couldn't load BXY file. Error message:\n%s",e.GetMsg());
				GLMessageBox::Display(tmp,"BXY file problem",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
		}
	}

	EminText->GetNumber(&cr->energy_low);
	EmaxText->GetNumber(&cr->energy_hi);
	cr->enable_ort_polarization=enableOrtPolarizationToggle->GetState();
	cr->enable_par_polarization=enableParPolarizationToggle->GetState();
	if (limitAngleToggle->GetState()) { //limit angle
		psiMaxXtext->GetNumber(&cr->psimaxX);
		psiMaxYtext->GetNumber(&cr->psimaxY);
	} else { 
		cr->psimaxX=PI;
		cr->psimaxY=PI;
	}

	//X comp.
	cr->Bx_mode=BxtypeCombo->GetSelectedIndex()+1;
	if (cr->Bx_mode==B_MODE_CONSTANT) { //constant B field
		constBXtext->GetNumber(&cr->B_const.x);
	} else { //use MAG file
		cr->MAGXfileName.assign(MAGfileXtext->GetText());
		try {
			FileReader MAGfile(MAGfileXtext->GetText());
			*(cr->Bx_distr)=cr->LoadMAGFile(&MAGfile,&cr->Bx_dir,&cr->Bx_period,&cr->Bx_phase,cr->Bx_mode);
		} catch(Error &e) {
			throw e;
			return;
		}
	}

	//Y comp.
	cr->By_mode=BytypeCombo->GetSelectedIndex()+1;
	if (cr->By_mode==B_MODE_CONSTANT) { //constant B field
		constBYtext->GetNumber(&cr->B_const.y);
	} else { //use MAG file
		cr->MAGYfileName.assign(MAGfileYtext->GetText());
		try {
			FileReader MAGfile(MAGfileYtext->GetText());
			*(cr->By_distr)=cr->LoadMAGFile(&MAGfile,&cr->By_dir,&cr->By_period,&cr->By_phase,cr->By_mode);
		} catch(Error &e) {
			throw e;
			return;
		}
	}

	//Z comp.
	cr->Bz_mode=BztypeCombo->GetSelectedIndex()+1;
	if (cr->Bz_mode==B_MODE_CONSTANT) { //constant B field
		constBZtext->GetNumber(&cr->B_const.z);
	} else { //use MAG file
		cr->MAGZfileName.assign(MAGfileZtext->GetText());
		try {
			FileReader MAGfile(MAGfileZtext->GetText());
			*(cr->Bz_distr)=cr->LoadMAGFile(&MAGfile,&cr->Bz_dir,&cr->Bz_period,&cr->Bz_phase,cr->Bz_mode);
		} catch(Error &e) {
			throw e;
			return;
		}
	}

	worker->RecalcRegion(regionId);
	if (mApp->regionInfo) mApp->regionInfo->Update();
	if (mApp->spectrumPlotter) mApp->spectrumPlotter->SetScale();
	GLWindow::ProcessMessage(NULL,MSG_CLOSE);
}