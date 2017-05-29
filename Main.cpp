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

//	thePhysicsEngine.initializeWorld_AuxeticSwisscheeseCell();
	thePhysicsEngine.initializeWorld_AuxeticPolygonCell();

	// graphics engine initialization -----------------------------------------
	GraphicsEngine theGraphicsEngine;
	theGraphicsEngine.initializeSystems();
	theGraphicsEngine.setPhysicsEngineReference(&thePhysicsEngine);

	theGraphicsEngine.Run();

	return(0);
}

