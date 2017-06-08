#include "GraphicsEngine.h"

// ----------------------------------------------------------------------------
void GraphicsEngine::runVisualization(PhysicsEngine *pPhysicsEngine)
{
	this->initializeSystems();
	this->setPhysicsEngineReference(pPhysicsEngine);

	this->gameLoop();
}
// ----------------------------------------------------------------------------
