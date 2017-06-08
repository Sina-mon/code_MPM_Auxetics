#include "GraphicsEngine.h"

// ----------------------------------------------------------------------------
void GraphicsEngine::initializeShaders(void)
{
	{// basic shader
		gl_BasicProgram.compileShaders("./Shader_Basic_Vertex.gl", "./Shader_Basic_Fragment.gl");
		gl_BasicProgram.addAtribute("vertexPosition");
		gl_BasicProgram.addAtribute("vertexColor");
		gl_BasicProgram.addAtribute("vertexNormal");
		gl_BasicProgram.addAtribute("vertexTextureCoordinate");
		gl_BasicProgram.linkShaders();
	}
	{// shadow shader
		gl_ShadowProgram.compileShaders("./Shader_Shadow_Vertex.gl", "./Shader_Shadow_Fragment.gl");
		gl_ShadowProgram.addAtribute("vertexPosition");
		gl_ShadowProgram.addAtribute("vertexColor");
		gl_ShadowProgram.addAtribute("vertexNormal");
		gl_ShadowProgram.addAtribute("vertexTextureCoordinate");
		gl_ShadowProgram.linkShaders();
	}
	{// final shader
		gl_FinalProgram.compileShaders("./Shader_Final_Vertex.gl", "./Shader_Final_Fragment.gl");
		gl_FinalProgram.addAtribute("vertexPosition");
		gl_FinalProgram.addAtribute("vertexColor");
		gl_FinalProgram.addAtribute("vertexNormal");
		gl_FinalProgram.addAtribute("vertexTextureCoordinate");
		gl_FinalProgram.linkShaders();
	}
	// Get the uniform variables location. You've probably already done that before...
	// http://stackoverflow.com/questions/25252512/how-can-i-pass-multiple-textures-to-a-single-shader
//	GLuint gl_UniformLocation_Diffuse	= gl_ColorProgram.getUniformLocation("diffuseTexture");
//	GLuint gl_UniformLocation_Normal	= gl_ColorProgram.getUniformLocation("normalTexture");

	// Then bind the uniform samplers to texture units:
//	gl_ColorProgram.use();
//	glUniform1i(gl_UniformLocation_Diffuse, 0);
//	glUniform1i(gl_UniformLocation_Normal,  1);
//	gl_ColorProgram.unuse();
}
// ----------------------------------------------------------------------------
