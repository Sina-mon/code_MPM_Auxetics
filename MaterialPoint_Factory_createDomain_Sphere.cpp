#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint *> MaterialPoint_Factory::createDomain_Sphere(glm::dvec3 d3Center, double dRadius_Outer, double dRadius_Inner, double dOffset)
{
	std::vector<MaterialPoint *> allMaterialPoint;

	for(double dx = 0.5*dOffset; dx < dRadius_Outer; dx += dOffset)
	{//create an 8th
		for(double dy = 0.5*dOffset; dy < dRadius_Outer; dy += dOffset)
		{
			for(double dz = 0.5*dOffset; dz <= dRadius_Outer; dz += dOffset)
			{
				if(dRadius_Inner*dRadius_Inner < dx*dx + dy*dy + dz*dz && dx*dx + dy*dy + dz*dz < dRadius_Outer*dRadius_Outer)
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
	}
	return(allMaterialPoint);
}
// ----------------------------------------------------------------------------

