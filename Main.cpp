#include <iostream>
#include <vector>

#include <omp.h>

#define STB_IMAGE_IMPLEMENTATION

#include "PhysicsEngine.h"
#include "GraphicsEngine.h"
#include "Definitions.h"

#include "ConstitutiveRelation.h"

int main (int argc, char ** argv)
{
	// physics engine initialization ------------------------------------------
	PhysicsEngine thePhysicsEngine;

	thePhysicsEngine.initializeWorld_AuxeticSwisscheeseCell();
//	thePhysicsEngine.initializeWorld_AuxeticPolygonCell();
//	thePhysicsEngine.initializeWorld_ContactBlock();
//	thePhysicsEngine.initializeWorld_DropWeight();
//	thePhysicsEngine.initializeWorld_Constitutive();
//	thePhysicsEngine.initializeWorld_AluminumTube();
//	thePhysicsEngine.initializeWorld_Billet();
//	thePhysicsEngine.initializeWorld_DiskRoll();
//	thePhysicsEngine.initializeWorld_SphereImpact();
//	thePhysicsEngine.initializeWorld_ConstitutiveRelation();
//	thePhysicsEngine.initializeWorld_Ring();
//	thePhysicsEngine.initializeWorld_DisplacementControlled();
//	thePhysicsEngine.initializeWorld_AuxeticCell();
//	thePhysicsEngine.initializeWorld_AuxeticMesh();
//	thePhysicsEngine.initializeWorld_Ring_T4();
//	thePhysicsEngine.initializeWorld_SteelAluminum();

	// graphics engine initialization -----------------------------------------
	GraphicsEngine theGraphicsEngine;
	theGraphicsEngine.initializeSystems();
	theGraphicsEngine.setPhysicsEngineReference(&thePhysicsEngine);

	theGraphicsEngine.Run();

	return(0);
}

