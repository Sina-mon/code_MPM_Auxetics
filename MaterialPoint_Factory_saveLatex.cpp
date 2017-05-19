#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
void MaterialPoint_Factory::saveLatex(std::vector<MaterialPoint *> allMaterialPoint, std::string strPrefix, std::string strPostfix)
{
	std::string strFileName = strPrefix + "MaterialPoint" + strPostfix + ".txt";
	std::ofstream OutputFile(strFileName.c_str(), std::ios_base::out);

	OutputFile << MaterialPoint_Factory::getScript(NULL);

	for(unsigned int index = 0; index < allMaterialPoint.size(); index+=1)
	{
		OutputFile << MaterialPoint_Factory::getScript(allMaterialPoint[index]);
	}

	OutputFile.close();
}
// ----------------------------------------------------------------------------

