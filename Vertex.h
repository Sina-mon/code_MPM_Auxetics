#ifndef VERTEX_H
#define VERTEX_H

#define GLEW_STATIC
#include <GL/glew.h>

#include "Definitions.h"
//#include <glm/glm.hpp>

struct Color
{
	GLubyte r;
	GLubyte g;
	GLubyte b;
	GLubyte a;
};

struct Vertex
{
	glm::vec3 f3Position = glm::vec3(0.0, 0.0, 0.0);
	glm::vec4 f4Color = glm::vec4(0.0, 0.0, 0.0, 0.0);
	glm::vec3 f3Normal = glm::vec3(0.0, 0.0, 0.0);
	glm::vec2 f2TextureCoordinate = glm::vec2(0.0, 0.0);
};

#endif // VERTEX_H
