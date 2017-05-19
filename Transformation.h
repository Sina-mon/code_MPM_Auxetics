#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include "Definitions.h"
//#include <glm/glm.hpp>
//#include <glm/gtx/transform.hpp>

#define GLEW_STATIC
#include <GL/glew.h>

#include "Vertex.h"

class Transformation
{
	public:
		Transformation(const glm::vec3 &f3Position, const glm::vec3 &f3Rotation, const glm::vec3 &f3Scale)
		{
			f3_Position = f3Position;
			f3_Rotation = f3Rotation;
			f3_Scale = f3Scale;
		}
		virtual ~Transformation() {}

		glm::mat4 GetModelMatrix(void) const
		{
			glm::mat4 m4Transformation_Position = glm::translate(f3_Position);
			glm::mat4 m4Transformation_RotationX = glm::rotate(f3_Rotation.x, glm::vec3(1.0, 0.0, 0.0));
			glm::mat4 m4Transformation_RotationY = glm::rotate(f3_Rotation.y, glm::vec3(0.0, 1.0, 0.0));
			glm::mat4 m4Transformation_RotationZ = glm::rotate(f3_Rotation.z, glm::vec3(0.0, 0.0, 1.0));
			glm::mat4 m4Transformation_Scale = glm::scale(f3_Scale);

			glm::mat4 m4Transformation_Combined;
			m4Transformation_Combined *= m4Transformation_Position;
			m4Transformation_Combined *= m4Transformation_RotationZ;
			m4Transformation_Combined *= m4Transformation_RotationY;
			m4Transformation_Combined *= m4Transformation_RotationX;
			m4Transformation_Combined *= m4Transformation_Scale;

			return(m4Transformation_Combined);
		}

		glm::vec3 f3_Position = glm::vec3(0.0, 0.0, 0.0);
		glm::vec3 f3_Rotation = glm::vec3(0.0, 0.0, 0.0);
		glm::vec3 f3_Scale = glm::vec3(1.0, 1.0, 1.0);
	protected:

	private:
};

#endif // TRANSFORMATION_H
