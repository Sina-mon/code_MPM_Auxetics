// https://www.opengl.org/wiki/Shader_Compilation
// to add uniform value, you only need to add it to the shader code and ask for its location using GLSLProgram::getUniformLocation
#ifndef GLSLPROGRAM_H
#define GLSLPROGRAM_H

#include <iostream>
#include <fstream>
#include <vector>

#define GLEW_STATIC
#include <GL/glew.h>

#include "Errors.h"

class GLSLProgram
{
	public:
		GLSLProgram();
		virtual ~GLSLProgram();

		void compileShaders(const std::string &strVertexShaderFile, const std::string &stdFragmentShaderFile);

		void linkShaders(void);

		void addAtribute(const std::string strAttributeName);

		GLuint getUniformLocation(const std::string strUniformName);

		void use(void);
		void unuse(void);

	protected:
		int i_numAttributes = 0;

		GLuint gl_ProgramID = 0;

		GLuint gl_VertexShaderID = 0;
		GLuint gl_FragmentShaderID = 0;

		void compileShader(const std::string &strShaderFile, GLuint glShaderID);

	private:
};

#endif // GLSLPROGRAM_H
