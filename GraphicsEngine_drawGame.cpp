#include "GraphicsEngine.h"

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

	if(true)
	{// MP parameteric value, e.g. mass, velocity, plastic strain
		v_Canvas_Texture[(int)enum_Canvas::J2]->bindRenderTarget();
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

		// material points ----------------------------------------------------
		ConstitutiveRelation CR;
		float fJ2_Maximum = 1.0e-12;
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint *thisMP = vMaterialPoint[index_MP];

			CR.calculateState_J2(thisMP->d6_Stress);
			float fJ2 = CR.d_J2;

			if(fJ2 > fJ2_Maximum)
				fJ2_Maximum = fJ2;
		}

		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 0.0004;
			// particle color
			CR.calculateState_J2(thisMP->d6_Stress);
			float fJ2 = CR.d_J2;
			glm::vec4 f4objectColor = (1.0f-fJ2/fJ2_Maximum) * _BLUE + fJ2/fJ2_Maximum * _RED;
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

	if(false)
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

	if(false)
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

	if(false)
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

	if(false)
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

		float fScreenRatio = 1.0 / (int)enum_Canvas::COUNT;

		if(true)
		{// top-left quadrant of screen
			// bind texture
			gl_Canvas_Texture->bindTextureUnit(0);
			// set viewport
			float x_Location = (float)enum_Canvas::MAIN / (int)enum_Canvas::COUNT;
			glViewport(x_Location*i_ScreenWidth, 0.0*i_ScreenHeight, fScreenRatio*i_ScreenWidth, 1.0*i_ScreenHeight);
			//draw
			gl_Canvas_Mesh->Draw();
		}
		if((int)enum_Canvas::COUNT > 1)
		{// top-left quadrant of screen
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::J2]->bindTextureUnit(0);
			// set viewport
			float x_Location = (float)enum_Canvas::J2 / (int)enum_Canvas::COUNT;
			glViewport(x_Location*i_ScreenWidth, 0.0*i_ScreenHeight, fScreenRatio*i_ScreenWidth, 1.0*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::J2]->Draw();
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
