#ifndef TEXTURE_H
#define TEXTURE_H

#include <iostream>
#include <string>

#define GLEW_STATIC
#include <GL/glew.h>
#include "stb_image.h" // sina, this doesn;t work if you put it in the header file


class Texture
{
	public:
		Texture(const std::string &strFileName)
		{//constructor to create a texture from an existing file
//			stbi_set_flip_vertically_on_load(1);
			int iWidth = 0;
			int iHeight = 0;
			int nComponents = 0;

			unsigned char * pImageData = stbi_load(strFileName.c_str(), &i_Width, &i_Height, &nComponents, 4); // 4 for spedify the required components

			if(pImageData == NULL)
				std::cout << "Texture loading failed for file " << strFileName << std::endl;

			glGenTextures(1, &gl_Texture);
			glBindTexture(GL_TEXTURE_2D, gl_Texture);

			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // if reading from outside the texture boundaries
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // if reading from outside the texture boundaries

			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); // for scaling
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); // for scaling

			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, i_Width, i_Height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pImageData); // send texture to gpu, 0 for default mipmaping, 0 for border, first RGBA is for gpu, sencond RGBA is about the data that we are sending in

			stbi_image_free(pImageData);
		}

		Texture(int iScreenWidth, int iScreenHeight)
		{//constructor to create a texture as the render target
			i_Width = iScreenWidth;
			i_Height = iScreenHeight;
			glGenFramebuffers(1, &gl_FrameBuffer);
			glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl_FrameBuffer);
			//texture creation
			glGenTextures(1, &gl_Texture);
			glBindTexture(GL_TEXTURE_2D, gl_Texture);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, i_Width, i_Height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);//if texture is smaller than resolution
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);

			//depth buffer, sina, this has to be here for depth test to work properly
			glGenRenderbuffers(1, &gl_DepthBuffer);
			glBindRenderbuffer(GL_RENDERBUFFER, gl_DepthBuffer);
			glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, i_Width, i_Height);
			glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, gl_DepthBuffer);

			glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl_Texture, 0);//bind texture to frame buffer

			glDrawBuffer(GL_COLOR_ATTACHMENT0);
			//check for errors
			if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			{
				std::cout << "TextureCC::TextureCC, error: problem with creating frame buffer." << std::endl;
			}
			//set screen as default rendering target
			glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
			glViewport(0, 0, iScreenWidth, iScreenHeight);
		}

		virtual ~Texture()
		{
			if(gl_Texture != 0)		glDeleteTextures(1, &gl_Texture);
			if(gl_DepthBuffer != 0)	glDeleteRenderbuffers(1, &gl_DepthBuffer);
			if(gl_FrameBuffer != 0)	glDeleteFramebuffers(1, &gl_FrameBuffer);
		}

		void bindTextureUnit(unsigned int iTextureUnit) // to active the intended texture
		{
			//sina, be careful, unit has to be between 1 and 32
			glActiveTexture(GL_TEXTURE0 + iTextureUnit); // a texture can have multiple units, this is where we specify which unit to use
			glBindTexture(GL_TEXTURE_2D, gl_Texture);
		}

		void bindRenderTarget(void)
		{//sina, has to be same size as when created
			glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl_FrameBuffer);
			glViewport(0, 0, i_Width, i_Height);
		}

	protected:
		GLuint gl_Texture = 0; // handle to the texture given by opengl
		GLuint gl_FrameBuffer = 0;
		GLuint gl_DepthBuffer = 0;

		int i_Width = 0;
		int i_Height = 0;
	private:
};

#endif // TEXTURE_H
