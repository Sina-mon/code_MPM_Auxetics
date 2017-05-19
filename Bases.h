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

		double d_ShapeValue = 0.0;
		glm::dvec3 d3_ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);

		void calculateBases(glm::dvec3 d3Position_MP, glm::dvec3 d3Position_GP, glm::dvec3 d3Cell_Length)
		{
			glm::dvec3 d3Distance = d3Position_MP - d3Position_GP;

			glm::dvec3 d3ShapeValue = glm::dvec3(0.0, 0.0, 0.0);
			d3ShapeValue[0] = 1.0 - fabs(d3Distance[0]) / d3Cell_Length[0];
			d3ShapeValue[1] = 1.0 - fabs(d3Distance[1]) / d3Cell_Length[1];
			d3ShapeValue[2] = 1.0 - fabs(d3Distance[2]) / d3Cell_Length[2];

			if(d3ShapeValue[0] < 0.0 || d3ShapeValue[0] > 1.0)
				d3ShapeValue[0] = 0.0;

			if(d3ShapeValue[1] < 0.0 || d3ShapeValue[1] > 1.0)
				d3ShapeValue[1] = 0.0;

			if(d3ShapeValue[2] < 0.0 || d3ShapeValue[2] > 1.0)
				d3ShapeValue[2] = 0.0;

			// results to pointers
			this->d_ShapeValue = d3ShapeValue[0] * d3ShapeValue[1] * d3ShapeValue[2];

			this->d3_ShapeGradient[0] = -d3ShapeValue[1] * d3ShapeValue[2] * glm::sign(d3Distance[0]) / d3Cell_Length[0];//sina, make sure these are formulated correctly
			this->d3_ShapeGradient[1] = -d3ShapeValue[0] * d3ShapeValue[2] * glm::sign(d3Distance[1]) / d3Cell_Length[1];
			this->d3_ShapeGradient[2] = -d3ShapeValue[0] * d3ShapeValue[1] * glm::sign(d3Distance[2]) / d3Cell_Length[2];
		}

	protected:

	private:
};

#endif // BASES_H
