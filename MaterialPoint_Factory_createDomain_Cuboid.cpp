#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint *> MaterialPoint_Factory::createDomain_Cuboid(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dOffset)
{
	std::vector<MaterialPoint *> allMaterialPoint;

	for(double dx = 0.5*dOffset; dx < 0.5*d3Dimension.x; dx += dOffset)
	{//create a quarter
		for(double dy = 0.5*dOffset; dy < 0.5*d3Dimension.y; dy += dOffset)
		{
			for(double dz = 0.5*dOffset; dz < 0.5*d3Dimension.z; dz += dOffset)
			{
				MaterialPoint *thisMaterialPoint;

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+dx, +dy, -dz), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-dx, +dy, -dz), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-dx, -dy, -dz), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+dx, -dy, -dz), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+dx, +dy, +dz), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-dx, +dy, +dz), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-dx, -dy, +dz), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+dx, -dy, +dz), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);
			}
		}
	}

	return(allMaterialPoint);
}
// ----------------------------------------------------------------------------

