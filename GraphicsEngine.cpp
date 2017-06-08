#include "GraphicsEngine.h"

GraphicsEngine::GraphicsEngine()
{
}
// ----------------------------------------------------------------------------
GraphicsEngine::~GraphicsEngine()
{
//	if(gl_Texture_01 != NULL)	delete gl_Texture_01;

	if(gl_Diffuse_Texture != NULL)	delete gl_Diffuse_Texture;

	if(gl_Camera != NULL)		delete gl_Camera;

	if(gl_Light != NULL)		delete gl_Light;

	if(gl_Shadow_Texture != NULL)	delete gl_Shadow_Texture;

	for(int index = 0; index < (int)enum_Canvas::COUNT; index++)
	{
		if(v_Canvas_Texture[index] != NULL)
			delete v_Canvas_Texture[index];

		if(v_Canvas_Mesh[index] != NULL)
			delete v_Canvas_Mesh[index];
	}

	if(gl_Canvas_Mesh != NULL)		delete gl_Canvas_Mesh;
	if(gl_Canvas_Texture != NULL)	delete gl_Canvas_Texture;

	if(gl_Canvas_Mesh_Grid_Mass != NULL)		delete gl_Canvas_Mesh_Grid_Mass;
	if(gl_Canvas_Texture_Grid_Mass != NULL)	delete gl_Canvas_Texture_Grid_Mass;

	if(gl_Canvas_Mesh_Grid_Kernel != NULL)		delete gl_Canvas_Mesh_Grid_Kernel;
	if(gl_Canvas_Texture_Grid_Kernel != NULL)	delete gl_Canvas_Texture_Grid_Kernel;

	if(gl_Canvas_Mesh_Grid_KernelGradient != NULL)		delete gl_Canvas_Mesh_Grid_KernelGradient;
	if(gl_Canvas_Texture_Grid_KernelGradient != NULL)	delete gl_Canvas_Texture_Grid_KernelGradient;

	if(gl_Canvas_Mesh_MP_KernelGradient != NULL)	delete gl_Canvas_Mesh_MP_KernelGradient;
	if(gl_Canvas_Texture_MP_KernelGradient != NULL)	delete gl_Canvas_Texture_MP_KernelGradient;

	if(gl_Particle_Mesh != NULL)	delete gl_Particle_Mesh;
}
// ----------------------------------------------------------------------------
