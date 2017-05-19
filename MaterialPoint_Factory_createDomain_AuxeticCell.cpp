#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint *> MaterialPoint_Factory::createDomain_AuxeticCell(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dOffset)
{
	std::vector<MaterialPoint *> allMaterialPoint;

	double dThickness = 0.0002;//2.0 * dOffset;
	double dMiddle_1 = 0.2;
	double dMiddle_2 = 0.3;
	std::vector<glm::dvec3> vVertex_Total;
	{
		vVertex_Total.push_back(glm::dvec3(0.0,0.0,0.0));
		vVertex_Total.push_back(glm::dvec3(0.5*d3Dimension.x, 0.0*d3Dimension.y, 0.0));
		vVertex_Total.push_back(glm::dvec3(0.5*d3Dimension.x, 0.5*d3Dimension.y, 0.0));
		vVertex_Total.push_back(glm::dvec3(0.0*d3Dimension.x, 0.5*d3Dimension.y, 0.0));
		vVertex_Total.push_back(glm::dvec3(0.0,0.0,0.0));
	}
	std::vector<glm::dvec3> vVertex_Top;
	{
		vVertex_Top.push_back(glm::dvec3(1.0*dThickness, dMiddle_1*d3Dimension.y + dThickness, 0.0));
		vVertex_Top.push_back(glm::dvec3(0.5*d3Dimension.x, dMiddle_2*d3Dimension.y + dThickness, 0.0));
		vVertex_Top.push_back(glm::dvec3(0.5*d3Dimension.x, 0.5*d3Dimension.y, 0.0));
		vVertex_Top.push_back(glm::dvec3(1.0*dThickness, 0.5*d3Dimension.y, 0.0));
		vVertex_Top.push_back(glm::dvec3(1.0*dThickness, dMiddle_1*d3Dimension.y + dThickness, 0.0));
	}
	std::vector<glm::dvec3> vVertex_Bottom;
	{
		vVertex_Total.push_back(glm::dvec3(0.0,0.0,0.0));
		vVertex_Bottom.push_back(glm::dvec3(0.5*d3Dimension.x-1.0*dThickness, 0.0*d3Dimension.y, 0.0));
		vVertex_Bottom.push_back(glm::dvec3(0.5*d3Dimension.x-1.0*dThickness, dMiddle_2*d3Dimension.y-dThickness, 0.0));
		vVertex_Bottom.push_back(glm::dvec3(0.0*d3Dimension.x, dMiddle_1*d3Dimension.y-dThickness, 0.0));
		vVertex_Total.push_back(glm::dvec3(0.0,0.0,0.0));
	}

	for(double dx = 0.5*dOffset; dx <= 0.5*d3Dimension.x-0.0*dOffset; dx += dOffset)
	{
		for(double dy = 0.5*dOffset; dy <= 0.5*d3Dimension.y-0.0*dOffset; dy += dOffset)
		{
			for(double dz = -0.5*d3Dimension.z; dz <= +0.5*d3Dimension.z; dz += dOffset)
			{
				MaterialPoint *thisMaterialPoint;

				glm::dvec3 d3Coordinate = glm::dvec3(dx, dy, dz);

				if(isInside(d3Coordinate, vVertex_Total) == true && isInside(d3Coordinate, vVertex_Top) == false && isInside(d3Coordinate, vVertex_Bottom) == false)
				{
					thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+d3Coordinate.y, +d3Coordinate.x, +d3Coordinate.z), dOffset);
					allMaterialPoint.push_back(thisMaterialPoint);

					thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-d3Coordinate.y, +d3Coordinate.x, +d3Coordinate.z), dOffset);
					allMaterialPoint.push_back(thisMaterialPoint);

					thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+d3Coordinate.y, -d3Coordinate.x, +d3Coordinate.z), dOffset);
					allMaterialPoint.push_back(thisMaterialPoint);

					thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-d3Coordinate.y, -d3Coordinate.x, +d3Coordinate.z), dOffset);
					allMaterialPoint.push_back(thisMaterialPoint);
				}
			}
		}
	}

	return(allMaterialPoint);
}
// ----------------------------------------------------------------------------

