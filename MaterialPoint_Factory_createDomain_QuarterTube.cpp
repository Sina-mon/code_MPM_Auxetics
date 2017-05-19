#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint *> MaterialPoint_Factory::createDomain_QuarterTube(glm::dvec3 d3Center, glm::dvec3 d3Rotation, double dRadius_Outer, double dRadius_Inner, double dLength, double dOffset)
{
	std::vector<MaterialPoint *> allMaterialPoint;

	for(double dx = 0.5*dOffset; dx < dRadius_Outer; dx += dOffset)
	{//create a quarter
		for(double dy = 0.5*dOffset; dy < dRadius_Outer; dy += dOffset)
		{
			for(double dz = 0.0; dz < 0.5*dLength; dz += dOffset)
			{
				double dSurface = dOffset;
				double dRadialDistance = glm::sqrt(dx*dx + dy*dy);
				if(dRadius_Inner*dRadius_Inner < dx*dx + dy*dy && dx*dx + dy*dy < dRadius_Outer*dRadius_Outer)
				{
					MaterialPoint *thisMaterialPoint;

					thisMaterialPoint = createMaterialPoint(glm::dvec3(+dx, +dy, -dz), dOffset);
					allMaterialPoint.push_back(thisMaterialPoint);

					if(glm::abs(dRadius_Inner - dRadialDistance) < dOffset)
						thisMaterialPoint->b_Surface = true;
					if(glm::abs(dRadius_Outer - dRadialDistance) < dOffset)
						thisMaterialPoint->b_Surface = true;
//							if(glm::abs(0.5*dLength - dz) < dOffset)
//								thisMaterialPoint->b_Surface = true;

					if(dz != 0.0)
					{
						thisMaterialPoint = createMaterialPoint(glm::dvec3(+dx, +dy, +dz), dOffset);
						allMaterialPoint.push_back(thisMaterialPoint);

						if(glm::abs(dRadius_Inner - dRadialDistance) < dOffset)
							thisMaterialPoint->b_Surface = true;
						if(glm::abs(dRadius_Outer - dRadialDistance) < dOffset)
							thisMaterialPoint->b_Surface = true;
//								if(glm::abs(0.5*dLength - dz) < dOffset)
//									thisMaterialPoint->b_Surface = true;
					}
				}
			}
		}
	}

	for(unsigned int index = 0; index < allMaterialPoint.size(); index++)
	{
		MaterialPoint *thisMaterialPoint = allMaterialPoint[index];

		glm::dmat4 m4Transformation_Position = glm::translate(d3Center);
		glm::dmat4 m4Transformation_RotationX = glm::rotate(d3Rotation.x, glm::dvec3(1.0, 0.0, 0.0));
		glm::dmat4 m4Transformation_RotationY = glm::rotate(d3Rotation.y, glm::dvec3(0.0, 1.0, 0.0));
		glm::dmat4 m4Transformation_RotationZ = glm::rotate(d3Rotation.z, glm::dvec3(0.0, 0.0, 1.0));
//				glm::dmat4 m4Transformation_Scale = glm::scale(f3_Scale);

		glm::dmat4 m4Transformation_Combined;
		m4Transformation_Combined *= m4Transformation_Position;
		m4Transformation_Combined *= m4Transformation_RotationZ;
		m4Transformation_Combined *= m4Transformation_RotationY;
		m4Transformation_Combined *= m4Transformation_RotationX;
//				m4Transformation_Combined *= m4Transformation_Scale;

		glm::dvec4 d4Position = glm::dvec4(thisMaterialPoint->d3_Position, 1.0);

		thisMaterialPoint->d3_Position = glm::dvec3(m4Transformation_Combined * d4Position);
	}

	return(allMaterialPoint);
}
// ----------------------------------------------------------------------------

