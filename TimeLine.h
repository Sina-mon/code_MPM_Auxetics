#ifndef TIMELINE_H
#define TIMELINE_H

#include <vector>

#include "Definitions.h"
////#include ".\include\glm\glm.hpp" // windows
//#include "./include/glm/glm.hpp" // linux, sina, glm is a column major implementation


class TimeLine
{
	public:
		TimeLine()
		{
//			v_Time.push_back(0.0);
//			v_Velocity.push_back(glm::dvec3(0.0, 0.0, 0.0));
//
//			v_Time.push_back(0.01);
//			v_Velocity.push_back(glm::dvec3(0.1, 0.0, 0.0));
//
//			v_Time.push_back(1.0e6);
//			v_Velocity.push_back(glm::dvec3(0.1, 0.0, 0.0));
		}
		virtual ~TimeLine() {}

		void addTimePoint(double dTime, glm::dvec3 d3Velocity)
		{
			v_Time.push_back(dTime);
			v_Velocity.push_back(d3Velocity);
		}

		glm::dvec3 getVelocity(double dTime)
		{
			glm::dvec3 d3Result = glm::dvec3(0.0, 0.0, 0.0);
			for(int index = 0; index < v_Time.size()-1; index++)
			{
				if(v_Time[index] <= dTime && dTime < v_Time[index+1])
				{
					double dScale = (dTime - v_Time[index]) / (v_Time[index+1] - v_Time[index]);
					d3Result = (1.0 - dScale)*v_Velocity[index] + (dScale)*v_Velocity[index+1];
				}
			}

			return(d3Result);
		}

		std::vector<glm::dvec3> v_Velocity;
		std::vector<double>		v_Time;
	protected:

	private:
};

#endif // TIMELINE_H
