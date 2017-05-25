#ifndef GRAPHICSENGINE_H
#define GRAPHICSENGINE_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <SDL2/SDL.h>

#include "Definitions.h"

#include "Errors.h"
#include "GLSLProgram.h"
#include "Vertex.h"
#include "Mesh.h"
#include "Texture.h"
#include "Transformation.h"
#include "Camera.h"
#include "Light.h"

#include "PhysicsEngine.h"
#include "MaterialPoint.h"
#include "GridPoint.h"

#define _RND(d)			((float)(rand()%d)/d)//random value between 0.0 and 1.0, with accuracy related to d
#define _RANDOM(a,b,d)	(a + _RND((int)d) * (b - a))

enum class GameState {PLAY, EXIT};

class GraphicsEngine
{
	public:
		GraphicsEngine();
		virtual ~GraphicsEngine();

		void initializeSystems(void);
		void setPhysicsEngineReference(PhysicsEngine *mpmPhysicsEngine);
		void drawGame(void);
		void Run(void);
		void saveScreenshot(int x, int y, int w, int h, const char * filename);
	protected:
		//snapshot save
		float f_TimeSnapshot_Interval = 10.0e-6;//100.0*5.0e-8;//0.1;
		float f_TimeSnapshot_LastSave = -1.0e12; // before creation
		int i_TimeSnapshotCycle = 0;

		SDL_Window *p_Window = NULL;

		float i_ScreenWidth = 600;
		float i_ScreenHeight = 600;

		glm::vec3 f3_World_Center = glm::vec3(0.0, 0.0, 0.0);
		glm::vec3 f3_World_Dimensions = glm::vec3(0.0, 0.0, 0.0);

		glm::vec3 f3_Camera_Position_Original = glm::vec3(0.0,0.0, -0.5);
		glm::vec3 f3_Camera_Target_Original = glm::vec3(0.0,0.0,0.0);

		GameState e_GameState = GameState::PLAY;

		// shaders
		GLSLProgram gl_BasicProgram;
		GLSLProgram gl_ShadowProgram;
		GLSLProgram gl_FinalProgram;

		Texture *gl_Texture_01;
		Texture *gl_Texture_02;

//		Transformation *gl_Transformation;

		Camera *gl_Camera;

		Light *gl_Light;

		Mesh *gl_Particle_Mesh;

		Texture	*gl_Shadow_Texture;

		enum class enum_Canvas : int {
			MAIN = 0,
			J2,
			COUNT
		};

		Texture	*v_Canvas_Texture[(int)enum_Canvas::COUNT];
		Mesh	*v_Canvas_Mesh[(int)enum_Canvas::COUNT];

		// canvas for every drawable item
		Texture	*gl_Canvas_Texture;
		Mesh	*gl_Canvas_Mesh;
		// canvas for MPM grid-mass data
		Texture	*gl_Canvas_Texture_Grid_Mass;
		Mesh	*gl_Canvas_Mesh_Grid_Mass;
		// canvas for MPM grid-kernel data
		Texture	*gl_Canvas_Texture_Grid_Kernel;
		Mesh	*gl_Canvas_Mesh_Grid_Kernel;
		// canvas for MPM grid-kernelGradient data
		Texture	*gl_Canvas_Texture_Grid_KernelGradient;
		Mesh	*gl_Canvas_Mesh_Grid_KernelGradient;
		// canvas for MPM materialPoint-kernelGradient data
		Texture	*gl_Canvas_Texture_MP_KernelGradient;
		Mesh	*gl_Canvas_Mesh_MP_KernelGradient;

		float f_Time = 0.0;

		PhysicsEngine *mpm_PhysicsEngine;
//		unsigned int n_Particles = 0;
//		std::vector<DEM_Particle *> v_DEMParticles;

		void initShaders(void);
		void gameLoop(void);
		void processInput(void);
	private:
};

#endif // GRAPHICSENGINE_H
