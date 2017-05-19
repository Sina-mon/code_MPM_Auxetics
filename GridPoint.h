#ifndef GRIDPOINT_H
#define GRIDPOINT_H

#include "Definitions.h"
////#include ".\include\glm\glm.hpp" // windows
//#include "./include/glm/glm.hpp" // linux

class GridPoint
{
	public:
		GridPoint() {;}
		virtual ~GridPoint() {;}

		bool b_Active = false;

		glm::ivec3 i3_Index = glm::vec3(0,0,0);

		bool b_Contact_Positive = false;// sina, this variable is only for graphic debugging
		bool b_Contact_Negative = false;// sina, this variable is only for graphic debugging

		glm::bvec3 b3_Fixed = glm::bvec3(false, false, false);
		glm::dvec3 d3_Position = glm::dvec3(0.0, 0.0, 0.0);

		glm::dvec3 d3_Mass			= glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_MassGradient	= glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_Velocity		= glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_Force			= glm::dvec3(0.0,0.0,0.0);

		glm::dvec3 d3_Force_Temp	= glm::dvec3(0.0,0.0,0.0);

		glm::dvec3 d3_Mass_Negative			= glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_MassGradient_Negative = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_Velocity_Negative		= glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_Force_Negative		= glm::dvec3(0.0,0.0,0.0);

		int i_NearestMP = 0;

		double		d_Kernel = 0.0;
		glm::dvec3	d3_Kernel_Gradient = glm::dvec3(0.0, 0.0, 0.0);
	protected:

	private:
};

#endif // GRIDPOINT_H
