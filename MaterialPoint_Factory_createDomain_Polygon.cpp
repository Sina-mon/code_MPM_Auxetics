#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint_Kinetics *> MaterialPoint_Factory::createDomain_Polygon(std::vector<glm::dvec3> vVertex, double dOffset)
{
	std::vector<MaterialPoint_Kinetics *> allMaterialPoint;

	glm::dvec3 d3Coordinates_Min = vVertex[0];
	glm::dvec3 d3Coordinates_Max = vVertex[0];

	// max,min values
	for(int index = 0; index < vVertex.size(); index++)
	{
		if(vVertex[index].x < d3Coordinates_Min.x)
			d3Coordinates_Min.x = vVertex[index].x;

		if(vVertex[index].y < d3Coordinates_Min.y)
			d3Coordinates_Min.y = vVertex[index].y;

		if(vVertex[index].z < d3Coordinates_Min.z)
			d3Coordinates_Min.z = vVertex[index].z;


		if(vVertex[index].x > d3Coordinates_Max.x)
			d3Coordinates_Max.x = vVertex[index].x;

		if(vVertex[index].y > d3Coordinates_Max.y)
			d3Coordinates_Max.y = vVertex[index].y;

		if(vVertex[index].z > d3Coordinates_Max.z)
			d3Coordinates_Max.z = vVertex[index].z;
	}

	for(double dx = d3Coordinates_Min.x; dx <= d3Coordinates_Max.x; dx += dOffset)
	{
		for(double dy = d3Coordinates_Min.y; dy <= d3Coordinates_Max.y; dy += dOffset)
		{
			for(double dz = d3Coordinates_Min.z; dz <= d3Coordinates_Max.z; dz += dOffset)
			{
				MaterialPoint_Kinetics *thisMaterialPoint;

//						double dz = vVertex[0].z;
				glm::dvec3 d3Coordinate = glm::dvec3(dx, dy, dz);

				if(isInside(d3Coordinate, vVertex) == true)
				{
					thisMaterialPoint = createMaterialPoint(d3Coordinate, dOffset);
					allMaterialPoint.push_back(thisMaterialPoint);
				}
			}
		}
	}

	return(allMaterialPoint);
}
// ----------------------------------------------------------------------------

