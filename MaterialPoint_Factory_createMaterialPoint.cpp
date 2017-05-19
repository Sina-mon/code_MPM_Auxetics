#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
MaterialPoint *MaterialPoint_Factory::createMaterialPoint(glm::dvec3 d3Center, double dOffset)
{
	MaterialPoint *thisMaterialPoint = new MaterialPoint();

	thisMaterialPoint->d3_Position = d3Center;

	return(thisMaterialPoint);
}
// ----------------------------------------------------------------------------

