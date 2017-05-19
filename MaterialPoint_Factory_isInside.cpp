#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
bool MaterialPoint_Factory::isInside(glm::dvec3 d3Coordinate, std::vector<glm::dvec3> vVertex)
{
	bool bInside = true;

	for(int index = 0; index < vVertex.size() - 1; index++)
	{// vertices should be in the counter-clockwise direction
		glm::dvec3 d3Vector1 = vVertex[index+1] - vVertex[index];
		glm::dvec3 d3Vector2 = d3Coordinate - vVertex[index];

		if(glm::cross(d3Vector1, d3Vector2).z < 0.0)
		{
			bInside = false;
			break;
		}
	}

	return(bInside);
}
// ----------------------------------------------------------------------------
bool MaterialPoint_Factory::isInside(glm::dvec3 d3Coordinate, glm::dvec3 d3Center, glm::dvec2 d2Radii)
{
	bool bInside = true;

	glm::dvec3 d3Coordinate_Local = d3Coordinate - d3Center;

	double dValue = glm::pow(d3Coordinate_Local.x/d2Radii.x, 2.0) + glm::pow(d3Coordinate_Local.y/d2Radii.y, 2.0);

	if(dValue > 1.0)
		bInside = false;

	return(bInside);
}
// ----------------------------------------------------------------------------

