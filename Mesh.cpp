#include "Mesh.h"
// ----------------------------------------------------------------------------
Mesh::Mesh(Vertex *pVertices, unsigned int nVertices, unsigned int *pIndices, unsigned int nIndices)
{// sina, this is not up do date
	this->Initialize(pVertices, nVertices, pIndices, nIndices);
}
// ----------------------------------------------------------------------------
Mesh::Mesh(const std::string& fileName)
{
	IndexedModel model = OBJModel(fileName).ToIndexedModel();

	unsigned int nVertices = model.positions.size();
	unsigned int nIndices = model.indices.size();

	Vertex 			*pVertices	= new Vertex[nVertices];
	unsigned int	*pIndices	= new unsigned int[nIndices];

	for(unsigned int index = 0; index < nVertices; index++)
	{
		pVertices[index].f3Position = model.positions[index];
		pVertices[index].f4Color = glm::vec4(0.8, 0.0, 0.2, 1.0);
		pVertices[index].f3Normal = model.normals[index];
		pVertices[index].f2TextureCoordinate = model.texCoords[index];
	}

	for(unsigned int index = 0; index < nIndices; index++)
		pIndices[index] = model.indices[index];


	this->Initialize(pVertices, nVertices, pIndices, nIndices);

	delete [] pVertices;
	delete [] pIndices;
}
// ----------------------------------------------------------------------------
Mesh::~Mesh()
{
	glDeleteVertexArrays(1, &gl_VertexArrayObject);
}
// ----------------------------------------------------------------------------
void Mesh::Initialize(Vertex *pVertices, unsigned int nVertices, unsigned int *pIndices, unsigned int nIndices)
{
	n_Vertices = nVertices;
	n_Indices = nIndices;

	// vertex array contains the vertex buffer, and information on how to read the memory
	// initial creation of the vertex array
	glGenVertexArrays(1, &gl_VertexArrayObject);

	// bind the VAO, so that everything we do applies to this particlal VOA
	glBindVertexArray(gl_VertexArrayObject);

	glGenBuffers(n_Buffers, gl_VertexArrayBuffers);// generate the memory to contain the vertex data

	// vertices buffer
	glBindBuffer(GL_ARRAY_BUFFER, gl_VertexArrayBuffers[0]); // 0 for 'vertices'
	glBufferData(GL_ARRAY_BUFFER, n_Vertices * sizeof(pVertices[0]), pVertices, GL_STATIC_DRAW); // move pVetices data to gpu memory, to currently bound buffer

	// tell the gpu how to interpret the data in the memory
	glEnableVertexAttribArray(0); // enable 0 (for position)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)offsetof(Vertex, f3Position)); // how to read in the memory data, 0 for position, 3 for x-y-z, GL_FALSE for no normalization, 0 skip on each interval, 0 to start from the very beginning

	glEnableVertexAttribArray(1); // enable 1 (for color)
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)offsetof(Vertex, f4Color)); // how to read in the memory data, 0 for position, 4 for rgba, GL_FALSE for no normalization, 0 skip on each interval, 0 to start from the very beginning

	glEnableVertexAttribArray(2); // enable 2 (for normal)
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)offsetof(Vertex, f3Normal)); // how to read in the memory data, 0 for position, 4 for rgba, GL_FALSE for no normalization, 0 skip on each interval, 0 to start from the very beginning

	glEnableVertexAttribArray(3); // enable 3 (for texture coordinate)
	glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)offsetof(Vertex, f2TextureCoordinate)); // how to read in the memory data, 0 for position, 2 for x-y, GL_FALSE for no normalization, 0 skip on each interval, 0 to start from the very beginning

	// indices buffer
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gl_VertexArrayBuffers[1]); // 1 for 'indices', GL_ELEMENT_ARRAY_BUFFER tell opengl that this buffer references another buffer
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_Indices * sizeof(pIndices[0]), pIndices, GL_STATIC_DRAW); // move pIndices data to gpu memory, to currently bound buffer

	glBindVertexArray(0);
}
// ----------------------------------------------------------------------------
void Mesh::Draw(void)
{
	glBindVertexArray(gl_VertexArrayObject); // bind the one we want to draw

//	glDrawArrays(GL_TRIANGLES, 0, n_Indices); // 0 to start from the start, and n to draw all
	glDrawElements(GL_TRIANGLES, n_Indices, GL_UNSIGNED_INT, 0);

	glBindVertexArray(0);
}
// ----------------------------------------------------------------------------
//		Vertex pV[3];
//		{
//			pV[0].f3Position = glm::vec3(0.0, 0.0, 0.0);
//			pV[1].f3Position = glm::vec3(0.05, 0.0, 0.0);
//			pV[2].f3Position = glm::vec3(0.0, 0.05, 0.0);
//
//			pV[0].f3Normal = glm::vec3(0.0, 0.0, 1.0);
//			pV[1].f3Normal = glm::vec3(0.0, 0.0, 1.0);
//			pV[2].f3Normal = glm::vec3(0.0, 0.0, 1.0);
//
//			pV[0].f2TextureCoordinate = glm::vec2(0.0, 0.0);
//			pV[1].f2TextureCoordinate = glm::vec2(0.1, 0.0);
//			pV[2].f2TextureCoordinate = glm::vec2(0.0, 0.1);
//
//			pV[0].f4Color = glm::vec4(1.0, 0.0, 0.0, 1.0);
//			pV[1].f4Color = glm::vec4(1.0, 0.0, 0.0, 1.0);
//			pV[2].f4Color = glm::vec4(1.0, 0.0, 0.0, 1.0);
//		}
//		unsigned int pI[3] = {0, 1, 2};
//		Mesh meshText(pV, 3, pI, 3);
//		glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
//		glm::mat4 m4TransformationMatrix_Model(1.0);
//		glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
//		glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
//		meshText.Draw();
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
