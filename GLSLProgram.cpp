#include "GLSLProgram.h"

// ----------------------------------------------------------------------------
GLSLProgram::GLSLProgram()
{
}
// ----------------------------------------------------------------------------
GLSLProgram::~GLSLProgram()
{
}
// ----------------------------------------------------------------------------
void GLSLProgram::compileShader(const std::string &strShaderFile, GLuint glShaderID)
{
	std::ifstream vertexFile(strShaderFile);
	if(vertexFile.fail())
	{
		perror(strShaderFile.c_str());
		fatalError("Error in GLSLProgram::compileShaders, failed to open " + strShaderFile);
	}

	std::string strFileContents;
	std::string strLine;
	while(std::getline(vertexFile, strLine))
	{
		strFileContents += strLine + "\n";
	}

	vertexFile.close();

	const char *pFileContents = strFileContents.c_str();
	glShaderSource(glShaderID, 1, &pFileContents, NULL);

	glCompileShader(glShaderID);

	GLint isCompiled = 0;
	glGetShaderiv(glShaderID, GL_COMPILE_STATUS, &isCompiled);

	if(isCompiled == GL_FALSE)
	{
		GLint maxLength = 0;
		glGetShaderiv(glShaderID, GL_INFO_LOG_LENGTH, &maxLength);

		// The maxLength includes the NULL character
		std::vector<GLchar> errorLog(maxLength);
		glGetShaderInfoLog(glShaderID, maxLength, &maxLength, &errorLog[0]);

		// Provide the infolog in whatever manor you deem best.
		// Exit with failure.
		glDeleteShader(glShaderID); // Don't leak the shader.

		std::printf("%s\n", &(errorLog[0]));
		fatalError("Error in GLSLProgram::compileShaders, shader " + strShaderFile + "failed to compile.");
	}
}
// ----------------------------------------------------------------------------
void GLSLProgram::compileShaders(const std::string &strVertexShaderFile, const std::string &strFragmentShaderFile)
{
	// https://www.opengl.org/wiki/Shader_Compilation

	//Vertex and fragment shaders are successfully compiled.
	//Now time to link them together into a program.
	//Get a program object.
	gl_ProgramID = glCreateProgram();

	gl_VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	if(gl_VertexShaderID == 0)
	{
		fatalError("Error in GLSLProgram::compileShaders, vertex shader failed to compile.");
	}

	gl_FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
	if(gl_FragmentShaderID == 0)
	{
		fatalError("Error in GLSLProgram::compileShaders, fragment shader failed to compile.");
	}

	this->compileShader(strVertexShaderFile, gl_VertexShaderID);
	this->compileShader(strFragmentShaderFile, gl_FragmentShaderID);
}
// ----------------------------------------------------------------------------
void GLSLProgram::linkShaders(void)
{
	//Attach our shaders to our program
	glAttachShader(gl_ProgramID, gl_VertexShaderID);
	glAttachShader(gl_ProgramID, gl_FragmentShaderID);

	//Link our program
	glLinkProgram(gl_ProgramID);

	//Note the different functions here: glGetProgram* instead of glGetShader*.
	GLint isLinked = 0;
	glGetProgramiv(gl_ProgramID, GL_LINK_STATUS, (int *)&isLinked);
	if(isLinked == GL_FALSE)
	{
		GLint maxLength = 0;
		glGetProgramiv(gl_ProgramID, GL_INFO_LOG_LENGTH, &maxLength);

		//The maxLength includes the NULL character
		std::vector<GLchar> infoLog(maxLength);
		glGetProgramInfoLog(gl_ProgramID, maxLength, &maxLength, &infoLog[0]);

		//We don't need the program anymore.
		glDeleteProgram(gl_ProgramID);
		//Don't leak shaders either.
		glDeleteShader(gl_VertexShaderID);
		glDeleteShader(gl_FragmentShaderID);

		std::printf("%s\n", &(infoLog[0]));
		fatalError("Error in GLSLProgram::compileShaders, shaders failed to compile.");
	}

	//Always detach shaders after a successful link.
	glDetachShader(gl_ProgramID, gl_VertexShaderID);
	glDetachShader(gl_ProgramID, gl_FragmentShaderID);
	glDeleteShader(gl_VertexShaderID);
	glDeleteShader(gl_FragmentShaderID);
}
// ----------------------------------------------------------------------------
void GLSLProgram::addAtribute(const std::string strAttributeName)
{
	glBindAttribLocation(gl_ProgramID, i_numAttributes, strAttributeName.c_str()); // defines the input of the shader locations
	i_numAttributes++;
}
// ----------------------------------------------------------------------------
GLuint GLSLProgram::getUniformLocation(const std::string strUniformName)
{
	GLuint location = glGetUniformLocation(gl_ProgramID, strUniformName.c_str());
	if(location == GL_INVALID_INDEX)
	{
		fatalError("Error in GLSLProgram::getUniformLocation, uniform " + strUniformName + " not found in shader.");
	}
	return(location);
}
// ----------------------------------------------------------------------------
void GLSLProgram::use(void)
{
	glUseProgram(this->gl_ProgramID);
	for(int i = 0; i < i_numAttributes; i++)
		glEnableVertexAttribArray(i);
}
// ----------------------------------------------------------------------------
void GLSLProgram::unuse(void)
{
	glUseProgram(0);
	for(int i = 0; i < i_numAttributes; i++)
		glDisableVertexAttribArray(i);
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
