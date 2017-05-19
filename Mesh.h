#ifndef MESH_H
#define MESH_H

#define GLEW_STATIC
#include <GL/glew.h>

#include "obj_loader.h"

#include "Vertex.h"

class Mesh
{
	public:
		Mesh(Vertex *pVertices, unsigned int nVertices, unsigned int *pIndices, unsigned int nIndices);
		Mesh(const std::string& fileName);
		virtual ~Mesh();

		void Initialize(Vertex *pVertices, unsigned int nVertices, unsigned int *pIndices, unsigned int nIndices);
		void Draw(void);

	protected:
		static const unsigned int n_Buffers = 2; //1 for vertices and 1 for indices
		GLuint gl_VertexArrayObject;
		GLuint gl_VertexArrayBuffers[n_Buffers];
		unsigned int n_Vertices = 0;
		unsigned int n_Indices = 0;

	private:
};

#endif // MESH_H
