#ifndef ERRORS_H
#define ERRORS_H

#include <iostream>

#include <SDL2/SDL.h>

extern void fatalError(std::string errorString);// extern keyword to tell the compiler that the function is defined elsewhere

#endif // ERRORS_H
