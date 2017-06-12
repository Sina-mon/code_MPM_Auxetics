#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint *> MaterialPoint_Factory::createDomain_AuxeticCell_Swisscheese(glm::dvec3 d3Center, glm::dvec3 d3Dimension, glm::dvec2 d2Spacing, glm::dvec2 d2Radii, double dOffset)
{
	std::vector<MaterialPoint *> allMaterialPoint;

	// create a quarter of the cell -----------------------
	std::vector<glm::dvec3> vVertex_Total;
	{// full box, quarter
		vVertex_Total.push_back(glm::dvec3(0.0,0.0,0.0));
		vVertex_Total.push_back(glm::dvec3(0.5,0.0,0.0)*d3Dimension);
		vVertex_Total.push_back(glm::dvec3(0.5,0.5,0.0)*d3Dimension);
		vVertex_Total.push_back(glm::dvec3(0.0,0.5,0.0)*d3Dimension);
		vVertex_Total.push_back(glm::dvec3(0.0,0.0,0.0));
	}

	std::vector<glm::dvec3> vEllipse_Center; // local coordinates
	std::vector<glm::dvec2> vEllipse_Radii; // local coordinates

	for(int iIndex_x = 0.0; iIndex_x < d3Dimension.x/d2Spacing.x; iIndex_x++)
	{
		for(int iIndex_y = 0.0; iIndex_y < d3Dimension.y/d2Spacing.y-2; iIndex_y++)
		{
			double dCenter_x = iIndex_x * d2Spacing.x;
			double dCenter_y = iIndex_y * d2Spacing.y;

			vEllipse_Center.push_back(glm::dvec3(dCenter_x, dCenter_y, 0.0));

			if((iIndex_x + iIndex_y) % 2 == 0)
				vEllipse_Radii.push_back(d2Radii);
			else
				vEllipse_Radii.push_back(glm::dvec2(d2Radii.y, d2Radii.x));
		}
	}

	// if the z-dimension is set to zero, create only one layer of material points
	// otherwise create even number of layers that are symmetric relative to the x-y plane
	bool bSingleLayer = false;
	double dz_Start = 0.5*dOffset;
	if(d3Dimension.z <= dOffset)
	{
		bSingleLayer = true;
		dz_Start = 0.0;
	}

	for(double dx = 0.5*dOffset; dx <= 0.5*d3Dimension.x; dx += dOffset)
	{
		for(double dy = 0.5*dOffset; dy <= 0.5*d3Dimension.y; dy += dOffset)
		{
			for(double dz = dz_Start; dz <= 0.5*d3Dimension.z; dz += dOffset)
			{
				MaterialPoint *thisMaterialPoint;

				glm::dvec3 d3Coordinate = glm::dvec3(dx, dy, dz);

				bool bCreate = true;
				for(int index = 0; index < vEllipse_Center.size(); index++)
				{
					if(isInside(d3Coordinate, vEllipse_Center[index], vEllipse_Radii[index]) == true)
					{
						bCreate = false;
						continue;
					}
				}

				if(bCreate == false)
					continue;

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+d3Coordinate.x, +d3Coordinate.y, +d3Coordinate.z), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-d3Coordinate.x, +d3Coordinate.y, +d3Coordinate.z), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+d3Coordinate.x, -d3Coordinate.y, +d3Coordinate.z), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-d3Coordinate.x, -d3Coordinate.y, +d3Coordinate.z), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				if(bSingleLayer != true)
				{
					thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+d3Coordinate.x, +d3Coordinate.y, -d3Coordinate.z), dOffset);
					allMaterialPoint.push_back(thisMaterialPoint);

					thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-d3Coordinate.x, +d3Coordinate.y, -d3Coordinate.z), dOffset);
					allMaterialPoint.push_back(thisMaterialPoint);

					thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+d3Coordinate.x, -d3Coordinate.y, -d3Coordinate.z), dOffset);
					allMaterialPoint.push_back(thisMaterialPoint);

					thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-d3Coordinate.x, -d3Coordinate.y, -d3Coordinate.z), dOffset);
					allMaterialPoint.push_back(thisMaterialPoint);
				}
			}
		}
	}

	return(allMaterialPoint);
}
// ----------------------------------------------------------------------------

