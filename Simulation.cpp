#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include "GLApp/MathTools.h"

SuperStructure::SuperStructure()
{
	aabbTree = NULL;
}

SuperStructure::~SuperStructure(){
	SAFE_DELETE(aabbTree);
}

Simulation::Simulation(){

	memset(&tmpGlobalResult, 0, sizeof(GlobalHitBuffer));
	

    nbLeakSinceUpdate = 0;

    totalDesorbed = 0;

    currentParticle.lastHitFacet = NULL;
    currentParticle.structureId = 0;
    currentParticle.teleportedFrom = 0;

    sh.nbSuper = 0;

    // Geometry
    nbMaterials = 0;
    sourceArea = 0;
    nbDistrPoints_BXY = 0;
    sourceRegionId = 0;

    stepPerSec = 0.0;
    textTotalSize = 0;
    profTotalSize = 0;
    dirTotalSize = 0;
    spectrumTotalSize = 0;
    loadOK = false;
    lastHitUpdateOK = false;
    lastLogUpdateOK = false;
    hasVolatile = false;

    wp.nbRegion = 0;
    wp.nbTrajPoints = 0;
    wp.newReflectionModel = false;
}