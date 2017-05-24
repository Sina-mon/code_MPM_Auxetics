#include "GraphicsEngine.h"

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
		gl_Camera->f3_Position = f3_Camera_Position_Original;
		gl_Camera->f3_Center = f3_Camera_Target_Original;

//		glm::vec3 f3Bounds = mpm_PhysicsEngine->d3_Length_Grid;//->getGridDimensions();
//
//		gl_Camera->f3_Center.x = 0.5*f3Bounds.x;
//		gl_Camera->f3_Center.y = 0.5*f3Bounds.y;
//		gl_Camera->f3_Center.z = 0.5*f3Bounds.z;
//
//		gl_Camera->f3_Position.x = 0.5*f3Bounds.x;
//		gl_Camera->f3_Position.y = 0.5*f3Bounds.y;
//		gl_Camera->f3_Position.z = 1.1*f3Bounds.x + f3Bounds.z;
	}
}
// ----------------------------------------------------------------------------
