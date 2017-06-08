#include "MaterialPoint_Factory.h"

// ----------------------------------------------------------------------------
std::string MaterialPoint_Factory::getScript(MaterialPoint *thisMaterialPoint)
{
	std::stringstream Stream;
	Stream.clear();

	if(thisMaterialPoint == NULL)
	{// to get the header column names
		Stream << "ID" << "\t";
		Stream << "mass" << "\t\t\t";
		Stream << "position_x" << "\t\t";
		Stream << "position_y" << "\t\t";
		Stream << "position_z" << "\t\t";
		Stream << std::endl;
	}
	else
	{
		Stream << Script(thisMaterialPoint->i_ID) << "\t";
		Stream << Script(thisMaterialPoint->d_Mass, 8) << "\t";
		Stream << Script(thisMaterialPoint->d3_Position.x, 8) << "\t";
		Stream << Script(thisMaterialPoint->d3_Position.y, 8) << "\t";
		Stream << Script(thisMaterialPoint->d3_Position.z, 8) << "\t";
		Stream << std::endl;
	}
	return (Stream.str ());
}
// ----------------------------------------------------------------------------

