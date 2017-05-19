#ifndef CAMERACC_H
#define CAMERACC_H

#define GLEW_STATIC

#include "Definitions.h"
//#include <glm/glm.hpp>
//#include <glm/gtx/transform.hpp>

class Camera
{
	public:
		Camera(const glm::vec3& pos, const glm::vec3& f3Center, float fov, float aspect, float zNear, float zFar)
		{
			m4_Projection = glm::perspective(fov, aspect, zNear, zFar);
			f3_Position = pos;
			f3_Center = f3Center;
			f3_Up = glm::vec3(0.0f, 1.0f, 0.0f);
		}
		virtual ~Camera(){;}

		void rotateAboutCenter(glm::vec3 f3Rotation)
		{
			glm::vec4 f4Radial = glm::vec4(this->f3_Position - this->f3_Center, 1.0);
			Transformation tTransformation(glm::vec3(0.0,0.0,0.0), f3Rotation, glm::vec3(1.0,1.0,1.0));
			f4Radial = tTransformation.GetModelMatrix() * f4Radial;
			this->f3_Position = this->f3_Center + glm::vec3(f4Radial.x, f4Radial.y, f4Radial.z);
		}

		void moveTowardCenter(float fDistance)
		{
			glm::vec3 f3Direction = this->f3_Center - this->f3_Position;
			glm::vec3 f3NewPosition = this->f3_Position + f3Direction*fDistance/glm::length(f3Direction);
			if(glm::length(f3NewPosition-this->f3_Center) > 0.02)
				this->f3_Position = f3NewPosition;
//				f3Direction = glm::normalize(gl_Camera->f3_Center - gl_Camera->f3_Position);
//			gl_Camera->f3_Position += f3Direction*(0.01f*theEvent.wheel.y);
		}

		glm::mat4 getViewProjection(void) const	{return(m4_Projection * glm::lookAt(f3_Position, f3_Center, f3_Up));}

		glm::mat4 m4_Projection;
		glm::vec3 f3_Position;
		glm::vec3 f3_Center;
		glm::vec3 f3_Up;
	protected:
	private:
};

#endif // CAMERACC_H
