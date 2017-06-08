#ifndef MATERIALPOINT_FACTORY_H
#define MATERIALPOINT_FACTORY_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "Definitions.h"
//#include <glm/glm.hpp>//sina, glm is a column major implementation
//#include <glm/gtx/transform.hpp>

//#include "Transformation.h"

#include "Definitions.h"
#include "MaterialPoint.h"

class MaterialPoint_Factory
{
	public:
		MaterialPoint_Factory(){;}
		virtual ~MaterialPoint_Factory(){;}

		bool isInside(glm::dvec3 d3Coordinate, std::vector<glm::dvec3> vVertex);
		bool isInside(glm::dvec3 d3Coordinate, glm::dvec3 d3Center, glm::dvec2 d2Radii); // is inside a ellipse

		MaterialPoint_Kinetics *createMaterialPoint(glm::dvec3 d3Center, double dOffset);

		std::vector<MaterialPoint_Kinetics *> createDomain_Cuboid(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dOffset);
		std::vector<MaterialPoint_Kinetics *> createDomain_Polygon(std::vector<glm::dvec3> vVertex, double dOffset);
		std::vector<MaterialPoint_Kinetics *> createDomain_AuxeticCell_Polygon(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dDent, double dThickness, double dOffset);
		std::vector<MaterialPoint_Kinetics *> createDomain_AuxeticCell_Swisscheese(glm::dvec3 d3Center, glm::dvec3 d3Dimension, glm::dvec2 d2Spacing, glm::dvec2 d2Radii, double dOffset);

		std::string getScript(MaterialPoint *thisMaterialPoint);

		void saveLatex(std::vector<MaterialPoint *> allMaterialPoint, std::string strPrefix, std::string strPostfix);

	protected:

	private:
};

#endif // MATERIALPOINT_FACTORY_H
