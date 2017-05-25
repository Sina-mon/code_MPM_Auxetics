#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <sstream>
#include<iomanip>

//#include ".\include\glm\glm.hpp" // windows
#include "./include/glm/glm.hpp"//sina, glm is a column major implementation
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#define _PI 3.14159265359

//-----------------------------------------------
#define _RED	glm::vec4(1.0, 0.0, 0.0, 1.0)
#define _GREEN	glm::vec4(0.0, 1.0, 0.0, 1.0)
#define _BLUE	glm::vec4(0.0, 0.0, 1.0, 1.0)
#define _WHITE	glm::vec4(1.0, 1.0, 1.0, 1.0)
#define _GRAY	glm::vec4(0.5, 0.5, 0.5, 1.0)
#define _BLACK	glm::vec4(0.0, 0.0, 0.0, 1.0)
//-----------------------------------------------
#define _STR_LOGFILE	 "./bmp/log.txt"
#define _STR_SNAPFILE	 "./bmp/Snapshot_"
//-----------------------------------------------
//-----------------------------------------------
static std::string Script(int i, int iWidth = 0)
{
	std::stringstream strScript;

	strScript.clear ();
	strScript.precision (0);
	strScript.setf (std::ios::dec);
	strScript.setf (std::ios::right);

	strScript.str("");
	if(iWidth != 0)
		strScript << std::setw(iWidth) << std::setfill('0');
	strScript << i;

	return (strScript.str ());
}
//-----------------------------------------------
static std::string Script(double r, int precision, std::ios::fmtflags flag=std::ios::scientific)
{
	std::stringstream strScript;

	strScript.clear ();
	strScript.precision (precision);
	strScript.setf (flag);
	strScript.setf (std::ios::right);
	strScript.setf (std::ios::showpos);

	strScript.str("");
	strScript << r;

	return (strScript.str());
}
//-----------------------------------------------
//-----------------------------------------------
//-----------------------------------------------


#endif // DEFINITIONS_H
