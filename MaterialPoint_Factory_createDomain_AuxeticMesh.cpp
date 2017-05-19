#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint *> MaterialPoint_Factory::createDomain_AuxeticMesh(glm::dvec3 d3Origin, glm::ivec2 i2Array, glm::dvec3 d3Dimension, double dOffset)
{
	std::vector<MaterialPoint *> allMaterialPoint;

	for(int index_x = 0; index_x < i2Array.x; index_x++)
	{
		for(int index_y = 0; index_y < i2Array.y; index_y++)
		{
			glm::dvec3 thisCell_Center = d3Origin + glm::dvec3(index_x*d3Dimension.y, index_y*d3Dimension.x, 0.0);

			std::vector<MaterialPoint *> thisCell = createDomain_AuxeticCell(thisCell_Center, d3Dimension, dOffset);

			allMaterialPoint.insert(std::end(allMaterialPoint), std::begin(thisCell), std::end(thisCell));
		}
	}

	return(allMaterialPoint);
}
// ----------------------------------------------------------------------------

