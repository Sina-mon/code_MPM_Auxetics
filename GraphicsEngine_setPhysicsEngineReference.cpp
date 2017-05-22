#include "GraphicsEngine.h"

// ----------------------------------------------------------------------------
void GraphicsEngine::setPhysicsEngineReference(PhysicsEngine *mpmPhysicsEngine)
{
	mpm_PhysicsEngine = mpmPhysicsEngine;

	f_TimeSnapshot_Interval = 1.0e4*mpm_PhysicsEngine->getTime_Increment();

	glm::vec3 f3Bounds = mpm_PhysicsEngine->d3_Length_Grid;//->getGridDimensions();

	gl_Camera->f3_Center.x = 0.5*f3Bounds.x;
	gl_Camera->f3_Center.y = 0.5*f3Bounds.y;
	gl_Camera->f3_Center.z = 0.5*f3Bounds.z;

	gl_Camera->f3_Position.x = 0.5*f3Bounds.x;
	gl_Camera->f3_Position.y = 0.5*f3Bounds.y;
	gl_Camera->f3_Position.z = 1.3*(f3Bounds.x + f3Bounds.y + f3Bounds.z);

//	this->v_DEMParticles.resize(mpm_PhysicsEngine->getCount_MaterialPoint() + mpm_PhysicsEngine->getCount_GridPoint());
//	this->v_DEMParticles.resize(mpm_PhysicsEngine->getCount_MaterialPoint());
}
// ----------------------------------------------------------------------------

