#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
MaterialPoint_Kinetics *MaterialPoint_Factory::createMaterialPoint(glm::dvec3 d3Center, double dOffset)
{
	MaterialPoint_Kinetics *thisMaterialPoint = new MaterialPoint_Kinetics();

	thisMaterialPoint->d3_Position = d3Center;

	return(thisMaterialPoint);
}
// ----------------------------------------------------------------------------

