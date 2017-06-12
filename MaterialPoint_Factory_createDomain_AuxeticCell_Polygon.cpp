#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint *> MaterialPoint_Factory::createDomain_AuxeticCell_Polygon(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dDent, double dThickness, double dOffset)
{
	std::vector<MaterialPoint *> allMaterialPoint;

	double dLength_Incline = glm::sqrt(glm::pow(0.5*d3Dimension.x - 2.0*dDent, 2.0) + glm::pow(0.5*d3Dimension.y, 2.0));
	double dAngle_Sine = 0.5*d3Dimension.y / dLength_Incline;
//	dAngle_Sine *= dAngle_Sine;
	double dAngle_Cosine = glm::sqrt(1.0 - glm::pow(dAngle_Sine, 2.0));
	double dAngle_Tan = dAngle_Sine/dAngle_Cosine;
	//std::cout << dAngle_Sine << std::endl;

	// create a quarter of the cell -----------------------
	std::vector<glm::dvec3> vVertex_Total;
	{// full box, quarter
		vVertex_Total.push_back(glm::dvec3(0.0,0.0,0.0));
		vVertex_Total.push_back(glm::dvec3(0.5,0.0,0.0)*d3Dimension);
		vVertex_Total.push_back(glm::dvec3(0.5,0.5,0.0)*d3Dimension);
		vVertex_Total.push_back(glm::dvec3(0.0,0.5,0.0)*d3Dimension);
		vVertex_Total.push_back(glm::dvec3(0.0,0.0,0.0));
	}
	std::vector<glm::dvec3> vVertex_Right;
	{// right cutout, quarter
		vVertex_Right.push_back(glm::dvec3(0.5*d3Dimension.x,0.5*dThickness,0.0));
		vVertex_Right.push_back(glm::dvec3(0.5*d3Dimension.x,0.5*d3Dimension.y,0.0));
		vVertex_Right.push_back(glm::dvec3(0.5*d3Dimension.x-dDent+0.5*dThickness/dAngle_Sine,0.5*d3Dimension.y,0.0));
		vVertex_Right.push_back(glm::dvec3(dDent+0.5*dThickness/dAngle_Sine+0.5*dThickness/dAngle_Tan,0.5*dThickness,0.0));
		vVertex_Right.push_back(glm::dvec3(0.5*d3Dimension.x,0.5*dThickness,0.0));
	}
	std::vector<glm::dvec3> vVertex_Left;
	{// left cutout, quarter
		vVertex_Left.push_back(glm::dvec3(0.0,0.0,0.0));
		vVertex_Left.push_back(glm::dvec3(dDent-0.5*dThickness/dAngle_Sine,0.0,0.0));
		vVertex_Left.push_back(glm::dvec3(0.5*d3Dimension.x-dDent-0.5*dThickness/dAngle_Sine-0.5*dThickness/dAngle_Tan,0.5*d3Dimension.y-0.5*dThickness,0.0));
		vVertex_Left.push_back(glm::dvec3(0.0,0.5*d3Dimension.y-0.5*dThickness,0.0));
		vVertex_Left.push_back(glm::dvec3(0.0,0.0,0.0));
	}
	/* gabriele's version
	// create a quarter of the cell ----------------------- gabriele's version
	std::vector<glm::dvec3> vVertex_Total;
	{// full box, quarter
		vVertex_Total.push_back(glm::dvec3(0.0,0.0,0.0));
		vVertex_Total.push_back(glm::dvec3(0.5,0.0,0.0)*d3Dimension);
		vVertex_Total.push_back(glm::dvec3(0.5,0.5,0.0)*d3Dimension);
		vVertex_Total.push_back(glm::dvec3(0.0,0.5,0.0)*d3Dimension);
		vVertex_Total.push_back(glm::dvec3(0.0,0.0,0.0));
	}
	std::vector<glm::dvec3> vVertex_Right;
	{// right cutout, quarter
		vVertex_Right.push_back(glm::dvec3(0.5*d3Dimension.x,0.0*dThickness,0.0));
		vVertex_Right.push_back(glm::dvec3(0.5*d3Dimension.x,0.5*d3Dimension.y,0.0));
		vVertex_Right.push_back(glm::dvec3(0.5*d3Dimension.x-dDent+0.5*dThickness/dAngle_Sine,0.5*d3Dimension.y,0.0));
		vVertex_Right.push_back(glm::dvec3(dDent+0.5*dThickness/dAngle_Sine,0.0*dThickness,0.0));
		vVertex_Right.push_back(glm::dvec3(0.5*d3Dimension.x,0.0*dThickness,0.0));
	}
	std::vector<glm::dvec3> vVertex_Left;
	{// left cutout, quarter
		vVertex_Left.push_back(glm::dvec3(0.0,0.0,0.0));
		vVertex_Left.push_back(glm::dvec3(dDent-0.5*dThickness/dAngle_Sine,0.0,0.0));
		vVertex_Left.push_back(glm::dvec3(0.5*d3Dimension.x-dDent-0.5*dThickness/dAngle_Sine-2.0*dThickness/dAngle_Tan,0.5*d3Dimension.y-2.0*dThickness,0.0));
		vVertex_Left.push_back(glm::dvec3(0.0,0.5*d3Dimension.y-2.0*dThickness,0.0));
		vVertex_Left.push_back(glm::dvec3(0.0,0.0,0.0));
	}
	*/
	for(double dx = 0.5*dOffset; dx <= 0.5*d3Dimension.x; dx += dOffset)
	{
		for(double dy = 0.5*dOffset; dy <= 0.5*d3Dimension.y; dy += dOffset)
		{
			for(double dz = 0.5*dOffset; dz <= 0.5*d3Dimension.z; dz += dOffset)
			{
				MaterialPoint *thisMaterialPoint;

				glm::dvec3 d3Coordinate = glm::dvec3(dx, dy, dz);

				if(isInside(d3Coordinate, vVertex_Total) != true)
					continue;
				if(isInside(d3Coordinate, vVertex_Right) == true)
					continue;
				if(isInside(d3Coordinate, vVertex_Left) == true)
					continue;

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+d3Coordinate.x, +d3Coordinate.y, +d3Coordinate.z), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-d3Coordinate.x, +d3Coordinate.y, +d3Coordinate.z), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(+d3Coordinate.x, -d3Coordinate.y, +d3Coordinate.z), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

				thisMaterialPoint = createMaterialPoint(d3Center + glm::dvec3(-d3Coordinate.x, -d3Coordinate.y, +d3Coordinate.z), dOffset);
				allMaterialPoint.push_back(thisMaterialPoint);

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

	return(allMaterialPoint);
}
// ----------------------------------------------------------------------------

