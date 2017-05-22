#include "GraphicsEngine.h"

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


