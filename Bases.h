#ifndef BASES_H
#define BASES_H

#include "Definitions.h"
////#include ".\include\glm\glm.hpp"
//#include "./include/glm/glm.hpp"

#include "MaterialPoint.h"

class Bases
{
	public:
		Bases() {;}
		virtual ~Bases() {;}

		double getShapeValue(void) {return(d_ShapeValue);}
		glm::dvec3 getShapeGradient(void) {return(d3_ShapeGradient);}

		double d_ShapeValue = 0.0;
		glm::dvec3 d3_ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);

		void calculateBases(glm::dvec3 d3Position_MP, glm::dvec3 d3Position_GP, glm::dvec3 d3Cell_Length)
		{
			glm::dvec3 d3Distance = d3Position_MP - d3Position_GP;

			glm::dvec3 d3ShapeValue = glm::dvec3(0.0, 0.0, 0.0);
			d3ShapeValue.x = 1.0 - fabs(d3Distance.x) / d3Cell_Length.x;
			d3ShapeValue.y = 1.0 - fabs(d3Distance.y) / d3Cell_Length.y;
			d3ShapeValue.z = 1.0 - fabs(d3Distance.z) / d3Cell_Length.z;

			if(d3ShapeValue.x < 0.0 || d3ShapeValue.x > 1.0)
				d3ShapeValue.x = 0.0;

			if(d3ShapeValue.y < 0.0 || d3ShapeValue.y > 1.0)
				d3ShapeValue.y = 0.0;

			if(d3ShapeValue.z < 0.0 || d3ShapeValue.z > 1.0)
				d3ShapeValue.z = 0.0;

			// results to pointers
			this->d_ShapeValue = d3ShapeValue.x * d3ShapeValue.y * d3ShapeValue.z;

			this->d3_ShapeGradient.x = -d3ShapeValue.y * d3ShapeValue.z * glm::sign(d3Distance.x) / d3Cell_Length.x;//sina, make sure these are formulated correctly
			this->d3_ShapeGradient.y = -d3ShapeValue.x * d3ShapeValue.z * glm::sign(d3Distance.y) / d3Cell_Length.y;
			this->d3_ShapeGradient.z = -d3ShapeValue.x * d3ShapeValue.y * glm::sign(d3Distance.z) / d3Cell_Length.z;
		}

	protected:

	private:
};

#endif // BASES_H
