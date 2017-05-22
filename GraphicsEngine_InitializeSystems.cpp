#include "GraphicsEngine.h"

// ----------------------------------------------------------------------------
void GraphicsEngine::initializeSystems(void)
{
	SDL_Init(SDL_INIT_EVERYTHING);

	p_Window = SDL_CreateWindow("MPM Suite", 10, 30, i_ScreenWidth, i_ScreenHeight, SDL_WINDOW_OPENGL);
	if(p_Window == NULL)	fatalError("Error int MainGame::initSystems, SDL window could not be created!");

	SDL_GLContext glContext = SDL_GL_CreateContext(p_Window);// passes the context to the window, so we don;t need to keep it
	if(glContext == NULL)	fatalError("Error int MainGame::initSystems, SDL_GL context could not be created!");

	GLenum thisError = glewInit();
	if(thisError != GLEW_OK)	fatalError("Error int MainGame::initSystems, could not initialize Glew!");

	SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_BUFFER_SIZE, 32);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 32);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

//	glEnable(GL_TEXTURE_2D);

//	glClearColor(0.0f, 0.05f, 0.10f, 1.0f);// only needs to be set once
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);// only needs to be set once

	this->initShaders();

	gl_Shadow_Texture	= new Texture(2048, 2048);

	gl_Canvas_Mesh		= new Mesh("./res/square.obj");
	gl_Canvas_Texture	= new Texture(i_ScreenWidth, i_ScreenHeight);

	gl_Canvas_Mesh_Grid_Mass	= new Mesh("./res/square.obj");
	gl_Canvas_Texture_Grid_Mass	= new Texture(i_ScreenWidth, i_ScreenHeight);

	gl_Canvas_Mesh_Grid_Kernel	= new Mesh("./res/square.obj");
	gl_Canvas_Texture_Grid_Kernel	= new Texture(i_ScreenWidth, i_ScreenHeight);

	gl_Canvas_Mesh_Grid_KernelGradient	= new Mesh("./res/square.obj");
	gl_Canvas_Texture_Grid_KernelGradient	= new Texture(i_ScreenWidth, i_ScreenHeight);

	gl_Canvas_Mesh_MP_KernelGradient	= new Mesh("./res/square.obj");
	gl_Canvas_Texture_MP_KernelGradient	= new Texture(i_ScreenWidth, i_ScreenHeight);

//	gl_Camera = new Camera(glm::vec3(0.0, 0.0, 2.0), glm::vec3(0.0, 0.0, 0.0), _PI/180.0*30.0f, 0.25*(float)(i_ScreenWidth/i_ScreenHeight), 0.01f, 0.5f);
	gl_Camera = new Camera(glm::vec3(0.0, 0.0, 2.0), glm::vec3(0.0, 0.0, 0.0), _PI/180.0*30.0f, (float)(i_ScreenWidth/i_ScreenHeight), 0.01f, 0.5f);

	gl_Light = new Light(glm::vec3(0.0,0.0,0.0), glm::vec3(0.025,0.025,-0.025));
	gl_Light->f4_Color = glm::vec4(1.0, 1.0, 1.0, 1.0);

	gl_Texture_01 = new Texture("./res/bricks.jpg");
	gl_Texture_02 = new Texture("./res/Sand_01.jpg");
	gl_Particle_Mesh = new Mesh("./res/sphere.obj");
}
// ----------------------------------------------------------------------------
