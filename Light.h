#ifndef LIGHT_H
#define LIGHT_H

#define GLEW_STATIC
#include "Definitions.h"
//#include <glm/glm.hpp>
//#include <glm/gtx/transform.hpp>

class Light
{
	public:
		Light(const glm::vec3& pos, const glm::vec3& f3Center)
		{
			m4_Projection = glm::ortho<float>(-0.025,0.025,-0.025,0.025,0.9,1.1);
			f3_Position = pos;
			f3_Center = f3Center;
			f3_Up = glm::vec3(0.0f, 1.0f, 0.0f);
		}
		virtual ~Light(){;}

		glm::vec3 getDirection(void) {return(f3_Center - f3_Position);}

		glm::mat4 getViewProjection(void) const
		{
			return(m4_Projection * glm::lookAt(f3_Position, f3_Center, f3_Up));
		}

		glm::mat4 m4_Projection;
		glm::vec3 f3_Position;
		glm::vec3 f3_Center;
		glm::vec3 f3_Up;

		glm::vec4 f4_Color = glm::vec4(1.0, 1.0, 1.0, 1.0);
	protected:
	private:
};

#endif // LIGHT_H
