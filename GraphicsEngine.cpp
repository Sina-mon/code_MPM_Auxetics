#include "GraphicsEngine.h"

GraphicsEngine::GraphicsEngine()
{
}
// ----------------------------------------------------------------------------
GraphicsEngine::~GraphicsEngine()
{
	if(gl_Texture_01 != NULL)	delete gl_Texture_01;

	if(gl_Texture_02 != NULL)	delete gl_Texture_02;

	if(gl_Camera != NULL)		delete gl_Camera;

	if(gl_Light != NULL)		delete gl_Light;

	if(gl_Shadow_Texture != NULL)	delete gl_Shadow_Texture;

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

//	for(unsigned int index = 0; index < v_DEMParticles.size(); index++)
//		delete v_DEMParticles[index];
}
// ----------------------------------------------------------------------------
void GraphicsEngine::Run(void)
{
//	this->initializeSystems();

	this->gameLoop();
}
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
void GraphicsEngine::initShaders(void)
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
void GraphicsEngine::gameLoop(void)
{
	int iSimulationStatus = 0;
	while(e_GameState != GameState::EXIT)
	{
		drawGame();

		double dTimeIncrement_Request = 1.0e2*mpm_PhysicsEngine->getTime_Increment();
		if(iSimulationStatus == 0)
		{
//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_SinglePass(dTimeIncrement_Request);
			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_SinglePass_MP(dTimeIncrement_Request);
//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_SinglePass_MP_Contact(dTimeIncrement_Request);
		}

		processInput();
		this->f_Time = mpm_PhysicsEngine->getTime_Current();

		// save snapshot
		if(f_Time - f_TimeSnapshot_LastSave > f_TimeSnapshot_Interval)
		{
			drawGame();
			drawGame();

			f_TimeSnapshot_LastSave = f_Time;
			std::string strFileName = _STR_SNAPFILE;//"./bmp/Snapshot_";
			strFileName += Script(i_TimeSnapshotCycle);
			strFileName += ".bmp";
			this->saveScreenshot(0, 0, i_ScreenWidth, i_ScreenHeight, strFileName.c_str());
			i_TimeSnapshotCycle += 1;
		}
	}
}
// ----------------------------------------------------------------------------
void GraphicsEngine::processInput(void)
{
	SDL_Event theEvent;

	while(SDL_PollEvent(&theEvent) == true)
	{
		if(theEvent.type == SDL_QUIT) e_GameState = GameState::EXIT;
		else if(theEvent.type == SDL_MOUSEMOTION)
		{ // mouse motion
			if(theEvent.motion.state == SDL_BUTTON_LMASK)
			{
				float fRotationIncrement_x = -0.002 * theEvent.motion.yrel;
				float fRotationIncrement_y = -0.002 * theEvent.motion.xrel;

				if(glm::abs(theEvent.motion.yrel) > glm::abs(theEvent.motion.xrel))
					gl_Camera->rotateAboutCenter(glm::vec3(fRotationIncrement_x,0.0,0.0));
				else
					gl_Camera->rotateAboutCenter(glm::vec3(0.0,fRotationIncrement_y,0.0));
			}
			if(theEvent.motion.state == SDL_BUTTON_MMASK)
			{
				glm::vec2 f2Speed = glm::vec2(1.0/i_ScreenWidth, 1.0/i_ScreenHeight) * gl_Camera->f3_Position.z;
				glm::vec3 f3Translation = glm::vec3(-f2Speed.x*theEvent.motion.xrel, f2Speed.y*theEvent.motion.yrel, 0.0);
				gl_Camera->f3_Position += f3Translation;
				gl_Camera->f3_Center += f3Translation;
			}
			if(theEvent.motion.state == SDL_BUTTON_RMASK)
			{
				float fRotationIncrement_x = -0.002 * theEvent.motion.yrel;
				float fRotationIncrement_y = -0.002 * theEvent.motion.xrel;

//				if(glm::abs(theEvent.motion.yrel) > glm::abs(theEvent.motion.xrel))
//					m_DE_Assembly.f3_Rotation_Assembly += glm::vec3(-fRotationIncrement_x,0.0,0.0);
//				else
//					m_DE_Assembly.f3_Rotation_Assembly += glm::vec3(0.0,-fRotationIncrement_y,0.0);
			}
		}
		else if(theEvent.type == SDL_MOUSEWHEEL)
		{ // mouse wheel
			float fScale = 0.1*glm::abs(gl_Camera->f3_Center.z - gl_Camera->f3_Position.z);
			gl_Camera->moveTowardCenter(fScale*theEvent.wheel.y);
		}
		else if(theEvent.type == SDL_KEYDOWN)
		{ // key press
			if(theEvent.key.keysym.sym == SDLK_UP)
			{
			}
			else if(theEvent.key.keysym.sym == SDLK_DOWN)
			{
			}
			else if(theEvent.key.keysym.sym == SDLK_LEFT)
			{
			}
			else if(theEvent.key.keysym.sym == SDLK_RIGHT)
			{
			}
		}
	}

	float fScrollSpeed = 0.2;
	const Uint8* currentKeyStates = SDL_GetKeyboardState( NULL );
	if(currentKeyStates[SDL_SCANCODE_UP])
	{
		glm::vec3 f3Translation = glm::vec3(0.0, fScrollSpeed, 0.0);
		gl_Camera->f3_Position += f3Translation;
		gl_Camera->f3_Center += f3Translation;
	}
	if(currentKeyStates[SDL_SCANCODE_DOWN])
	{
		glm::vec3 f3Translation = glm::vec3(0.0, -fScrollSpeed, 0.0);
		gl_Camera->f3_Position += f3Translation;
		gl_Camera->f3_Center += f3Translation;
	}
	if(currentKeyStates[SDL_SCANCODE_LEFT])
	{
		glm::vec3 f3Translation = glm::vec3(-fScrollSpeed, 0.0, 0.0);
		gl_Camera->f3_Position += f3Translation;
		gl_Camera->f3_Center += f3Translation;
	}
	if(currentKeyStates[SDL_SCANCODE_RIGHT])
	{
		glm::vec3 f3Translation = glm::vec3(fScrollSpeed, 0.0, 0.0);
		gl_Camera->f3_Position += f3Translation;
		gl_Camera->f3_Center += f3Translation;
	}
	if(currentKeyStates[SDL_SCANCODE_TAB])
	{
		glm::vec3 f3Bounds = mpm_PhysicsEngine->d3_Length_Grid;//->getGridDimensions();

		gl_Camera->f3_Center.x = 0.5*f3Bounds.x;
		gl_Camera->f3_Center.y = 0.5*f3Bounds.y;
		gl_Camera->f3_Center.z = 0.5*f3Bounds.z;

		gl_Camera->f3_Position.x = 0.5*f3Bounds.x;
		gl_Camera->f3_Position.y = 0.5*f3Bounds.y;
		gl_Camera->f3_Position.z = 1.1*f3Bounds.x + f3Bounds.z;
	}
}
// ----------------------------------------------------------------------------
void GraphicsEngine::drawGame(void)
{
	if(true)
	{// create shadow map
		gl_Shadow_Texture->bindRenderTarget();
//		gl_Canvas_Texture->bindRenderTarget(i_ScreenWidth, i_ScreenHeight);
		gl_ShadowProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		GLuint transformationShadowLocation = gl_ShadowProgram.getUniformLocation("transformationShadowMatrix");

		std::vector<MaterialPoint *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();

		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint *thisMP = vMaterialPoint[index_MP];
			// particle position
			float fSize = glm::pow(thisMP->d_Volume, 1.0/3.0);
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection() * glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// draw
			gl_Particle_Mesh->Draw();
		}

		gl_ShadowProgram.unuse();
	}

	if(true)
	{// draw all objects to the canvas (not the screen)
		gl_Canvas_Texture->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Texture_02->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		// material points ----------------------------------------------------
		std::vector<MaterialPoint *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 0.8*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			glm::vec4 f4objectColor = _RED;
			if(thisMP->b_Mark_Stress)
				f4objectColor = _GREEN;
			if(thisMP->b_Surface)
			{
				f4objectColor = _BLUE;
				fSize *= 1.0;
			}
			if(thisMP->b_DisplacementControl)
				f4objectColor = _WHITE;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}
		// grid points --------------------------------------------------------
		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();

		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(thisGP->b3_Fixed == glm::bvec3{false, false, false})	continue;

			// particle position
			float fSize = 0.02*mpm_PhysicsEngine->d3_Length_Cell.x;
			Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			// particle color
			glm::vec4 f4objectColor = _GRAY;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		gl_BasicProgram.unuse();
	}

	if(false)
	{// kernel
		gl_Canvas_Texture_Grid_Kernel->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Texture_02->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

//		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();
		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints_Kernel();

		float fKernel_Maximum = 0.0;
		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(thisGP->d_Kernel > fKernel_Maximum)	fKernel_Maximum = thisGP->d_Kernel;
		}

		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			// particle position
			float fSize = 0.002*thisGP->d_Kernel / fKernel_Maximum * mpm_PhysicsEngine->d3_Length_Cell.x;
			// particle color
			glm::vec4 f4objectColor = _WHITE;
//			if(thisGP->b_Contact_Positive == true && thisGP->b_Contact_Negative == true)
			if(thisGP->b_Contact_Negative == true)
			{
				f4objectColor = _RED;
				fSize *= 5.0;
			}

			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		gl_BasicProgram.unuse();
	}

	if(true)
	{// kernel gradient
		gl_Canvas_Texture_Grid_KernelGradient->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Texture_02->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();
//		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints_Kernel();

		float fKernelGradient_Maximum = 0.0;
		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(glm::length(thisGP->d3_Kernel_Gradient) > fKernelGradient_Maximum)
				fKernelGradient_Maximum = glm::length(thisGP->d3_Kernel_Gradient);
		}

//		std::cout << "Max kernel gradient: " << fKernelGradient_Maximum << std::endl;

		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			// particle position
			float fSize = 0.0002;// for problems with dimensions of 1m
			Transformation glTransformation(glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 0.0), glm::vec3(1.0));
			// color
			glm::vec4 f4objectColor = _GREEN;
			if(thisGP->b_Contact_Positive && thisGP->b_Contact_Negative)
				f4objectColor = _RED;
//			else
//				continue;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

//			gl_Particle_Mesh->Draw();

			if(glm::length(thisGP->d3_Kernel_Gradient) > 0.0)
			{
				glBegin(GL_LINES);
					glm::vec3 f3Start = thisGP->d3_Position;
					glm::vec3 f3End = glm::vec3(thisGP->d3_Position) + fSize/fKernelGradient_Maximum * glm::vec3(thisGP->d3_Kernel_Gradient);
					glVertex3f(f3Start.x, f3Start.y, f3Start.z);
					glVertex3f(f3End.x, f3End.y, f3End.z);
				glEnd();
			}
		}

		gl_BasicProgram.unuse();
	}

	if(true)
	{// kernel gradient
		gl_Canvas_Texture_Grid_Kernel->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Texture_02->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();
//		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints_Kernel();

		float fKernelGradient_Maximum = 0.0;
		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(glm::length(thisGP->d3_Kernel_Gradient) > fKernelGradient_Maximum)
				fKernelGradient_Maximum = glm::length(thisGP->d3_Kernel_Gradient);
		}

		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			// particle position
			float fSize = 0.0002;// for problems with dimensions of 1m
			Transformation glTransformation(glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 0.0), glm::vec3(1.0));
			// color
			glm::vec4 f4objectColor = _GREEN;
			if(thisGP->b_Contact_Positive && thisGP->b_Contact_Negative)
				f4objectColor = _RED;
			else
				continue;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

//			gl_Particle_Mesh->Draw();

			if(glm::length(thisGP->d3_Kernel_Gradient) > 0.0)
			{
				glBegin(GL_LINES);
					glm::vec3 f3Start = thisGP->d3_Position;
					glm::vec3 f3End = glm::vec3(thisGP->d3_Position) + fSize/fKernelGradient_Maximum * glm::vec3(thisGP->d3_Kernel_Gradient);
					glVertex3f(f3Start.x, f3Start.y, f3Start.z);
					glVertex3f(f3End.x, f3End.y, f3End.z);
				glEnd();
			}
		}

		gl_BasicProgram.unuse();
	}

	if(false)
	{// grid mass canvas
		gl_Canvas_Texture_Grid_Mass->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Texture_02->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();

		float fMass_Maximum = 0.0;
		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(thisGP->d3_Mass.x > fMass_Maximum)	fMass_Maximum = thisGP->d3_Mass.x;
		}

		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			// particle position
			float fSize = thisGP->d3_Mass.x / fMass_Maximum * mpm_PhysicsEngine->d3_Length_Cell.x;
			Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.5*fSize));
			// particle color
			glm::vec4 f4objectColor = _GRAY;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		gl_BasicProgram.unuse();
	}

	if(true)
	{// grid mass gradient
		gl_Canvas_Texture_Grid_Mass->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Texture_02->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();


		float fMassGradient_Maximum = 0.0;
		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(glm::length(thisGP->d3_MassGradient) > fMassGradient_Maximum)
				fMassGradient_Maximum = glm::length(thisGP->d3_MassGradient);
		}

		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			// particle position
			float fSize = 0.0002;// for problems with dimensions of 1m
			Transformation glTransformation(glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 0.0), glm::vec3(1.0));
			// color
			glm::vec4 f4objectColor = _GREEN;
//			if(thisGP->b_Contact_Positive && thisGP->b_Contact_Negative)
//				f4objectColor = _RED;
//			else
//				continue;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

//			gl_Particle_Mesh->Draw();

			if(glm::length(thisGP->d3_MassGradient) > 0.0)
			{
				glBegin(GL_LINES);
					glm::vec3 f3Start = thisGP->d3_Position;
					glm::vec3 f3End = glm::vec3(thisGP->d3_Position) + fSize/fMassGradient_Maximum * glm::vec3(thisGP->d3_MassGradient);
					glVertex3f(f3Start.x, f3Start.y, f3Start.z);
					glVertex3f(f3End.x, f3End.y, f3End.z);
				glEnd();
			}
		}

		gl_BasicProgram.unuse();
	}

	if(true)
	{// MP kernel gradient
		gl_Canvas_Texture_MP_KernelGradient->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Texture_02->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<MaterialPoint *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();

		float fKernelGradient_Maximum = 0.0;
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint *thisMP = vMaterialPoint[index_MP];

			if(glm::length(thisMP->d3_Kernel_Gradient) > fKernelGradient_Maximum)
				fKernelGradient_Maximum = glm::length(thisMP->d3_Kernel_Gradient);
		}

		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 0.0001;// for problems with dimensions of 1m
			Transformation glTransformation(glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 0.0), glm::vec3(1.0));
			// particle color
			glm::vec4 f4objectColor = _GREEN;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

//			if(glm::length(thisMP->d3_Kernel_Gradient) > 0.0)
			if(thisMP->b_Surface)
			{
				glBegin(GL_LINES);
					glm::vec3 f3Start = thisMP->d3_Position;
					glm::vec3 f3End = glm::vec3(thisMP->d3_Position) + fSize/fKernelGradient_Maximum * glm::vec3(thisMP->d3_Kernel_Gradient);
					glVertex3f(f3Start.x, f3Start.y, f3Start.z);
					glVertex3f(f3End.x, f3End.y, f3End.z);
				glEnd();
			}
		}

		gl_BasicProgram.unuse();
	}

	if(true)
	{// bind the screen for final output
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0); // drawing to the window

		// which .gl program to use
		gl_FinalProgram.use();

		// clear
		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if(true)
		{// top-left quadrant of screen
			// bind texture
			gl_Canvas_Texture->bindTextureUnit(0);
			// set viewport
			glViewport(0.0*i_ScreenWidth, 0.0*i_ScreenHeight, 1.0*i_ScreenWidth, 1.0*i_ScreenHeight);
			//draw
			gl_Canvas_Mesh->Draw();
		}
		if(false)
		{// top-right quadrant of screen
			// bind texture
			gl_Canvas_Texture_Grid_Kernel->bindTextureUnit(0);
			// set viewport
			glViewport(0.25*i_ScreenWidth, 0.0*i_ScreenHeight, 0.25*i_ScreenWidth, 1.0*i_ScreenHeight);
			//draw
			gl_Canvas_Mesh->Draw();
		}
		if(false)
		{// top-right quadrant of screen
			// bind texture
			gl_Canvas_Texture_Grid_Kernel->bindTextureUnit(0);
			// set viewport
			glViewport(0.5*i_ScreenWidth, 0.5*i_ScreenHeight, 0.5*i_ScreenWidth, 0.5*i_ScreenHeight);
			//draw
			gl_Canvas_Mesh->Draw();
		}
		if(false)
		{// top-right quadrant of screen
			// bind texture
			gl_Canvas_Texture_MP_KernelGradient->bindTextureUnit(0);
			// set viewport
			glViewport(0.5*i_ScreenWidth, 0.5*i_ScreenHeight, 0.5*i_ScreenWidth, 0.5*i_ScreenHeight);
			//draw
			gl_Canvas_Mesh->Draw();
		}
		if(false)
		{// bottom-left quadrant of screen
			// bind texture
			gl_Canvas_Texture_MP_KernelGradient->bindTextureUnit(0);
			// set viewport
			glViewport(0.50*i_ScreenWidth, 0.0*i_ScreenHeight, 0.25*i_ScreenWidth, 1.0*i_ScreenHeight);
			//draw
			gl_Canvas_Mesh->Draw();
		}
		if(false)
		{// bottom-right quadrant of screen
			// bind texture
			gl_Canvas_Texture_Grid_Mass->bindTextureUnit(0);
			// set viewport
			glViewport(0.50*i_ScreenWidth, 0.0*i_ScreenHeight, 0.25*i_ScreenWidth, 1.0*i_ScreenHeight);
			//draw
			gl_Canvas_Mesh->Draw();
		}
		if(false)
		{// bottom-left quadrant of screen
			// bind texture
			gl_Canvas_Texture_Grid_KernelGradient->bindTextureUnit(0);
			// set viewport
			glViewport(0.75*i_ScreenWidth, 0.0*i_ScreenHeight, 0.25*i_ScreenWidth, 1.0*i_ScreenHeight);
			//draw
			gl_Canvas_Mesh->Draw();
		}

		// undind texture
		gl_FinalProgram.unuse();
	}

	SDL_GL_SwapWindow(p_Window);
}
// ----------------------------------------------------------------------------
void GraphicsEngine::setPhysicsEngineReference(PhysicsEngine *mpmPhysicsEngine)
{
	mpm_PhysicsEngine = mpmPhysicsEngine;

	f_TimeSnapshot_Interval = 1.0e4*mpm_PhysicsEngine->getTime_Increment();

	glm::vec3 f3Bounds = mpm_PhysicsEngine->d3_Length_Grid;//->getGridDimensions();

	gl_Camera->f3_Center.x = 0.5*f3Bounds.x;
	gl_Camera->f3_Center.y = 0.5*f3Bounds.y;
	gl_Camera->f3_Center.z = 0.5*f3Bounds.z;

	gl_Camera->f3_Position.x = 0.5*f3Bounds.x;
	gl_Camera->f3_Position.y = 0.5*f3Bounds.y;
	gl_Camera->f3_Position.z = 1.3*(f3Bounds.x + f3Bounds.y + f3Bounds.z);

//	this->v_DEMParticles.resize(mpm_PhysicsEngine->getCount_MaterialPoint() + mpm_PhysicsEngine->getCount_GridPoint());
//	this->v_DEMParticles.resize(mpm_PhysicsEngine->getCount_MaterialPoint());
}
// ----------------------------------------------------------------------------
void GraphicsEngine::saveScreenshot(int x, int y, int w, int h, const char * filename)
{
    unsigned char * pixels_flipped = new unsigned char[w*h*4]; // 4 bytes for RGBA
    unsigned char * pixels = new unsigned char[w*h*4]; // 4 bytes for RGBA
    glReadPixels(x,y,w, h, GL_BGRA, GL_UNSIGNED_BYTE, pixels);

	// flip about horizontal
	for(int yCount = 0; yCount < h; yCount++)
	{
		for(int xCount = 0; xCount < w; xCount++)
		{
			int index_orig = (h-yCount-1)*w + xCount;
			int index_flip = yCount*w + xCount;
			pixels_flipped[4*index_flip + 0] = pixels[4*index_orig + 0];
			pixels_flipped[4*index_flip + 1] = pixels[4*index_orig + 1];
			pixels_flipped[4*index_flip + 2] = pixels[4*index_orig + 2];
			pixels_flipped[4*index_flip + 3] = pixels[4*index_orig + 3];
		}
	}

    SDL_Surface * surf = SDL_CreateRGBSurfaceFrom(pixels_flipped, w, h, 8*4, w*4, 0,0,0,0);
    SDL_SaveBMP(surf, filename);

    SDL_FreeSurface(surf);
    delete [] pixels;
    delete [] pixels_flipped;
}
// ----------------------------------------------------------------------------


