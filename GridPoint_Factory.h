#ifndef GRIDPOINT_FACTORY_H
#define GRIDPOINT_FACTORY_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "Definitions.h"
#include "GridPoint.h"

class GridPoint_Factory
{
	public:
		GridPoint_Factory() {;}
		virtual ~GridPoint_Factory() {;}

		static int getIndex(int i_x, int i_y, int i_z, int n_x, int n_y, int n_z)
		{
			int index = i_x + i_y*(n_x) + i_z*(n_x*n_y);//   i_x + i_y * n_x;
			return(index);
		}

		static int getIndex(glm::ivec3 i3Node_Index, glm::ivec3 i3Node_Count)
		{
			int index = i3Node_Index.x + i3Node_Index.y*(i3Node_Count.x) + i3Node_Index.z*(i3Node_Count.x*i3Node_Count.y);

			int iNodes = i3Node_Count.x * i3Node_Count.y * i3Node_Count.z;
			if(0 <= index && index < iNodes)
				return(index);
			else
				return(-1);
		}

		std::vector<GridPoint *> createGrid(glm::dvec3 d3Length, glm::ivec3 i3Cells)
		{
			int nNodes[3] = {0, 0, 0};
			nNodes[0] = i3Cells[0] + 1;
			nNodes[1] = i3Cells[1] + 1;
			nNodes[2] = i3Cells[2] + 1;

			std::vector<GridPoint *> allGridPoint;//
			allGridPoint.resize(nNodes[0] * nNodes[1] * nNodes[2]);

			for(int i_x = 0; i_x < nNodes[0]; i_x++)
			{
				for(int i_y = 0; i_y < nNodes[1]; i_y++)
				{
					for(int i_z = 0; i_z < nNodes[2]; i_z++)
					{
						GridPoint *thisGridPoint = new GridPoint;

						thisGridPoint->i3_Index = glm::ivec3(i_x, i_y, i_z);

						double dx = i_x * d3Length[0]/i3Cells[0];
						double dy = i_y * d3Length[1]/i3Cells[1];
						double dz = i_z * d3Length[2]/i3Cells[2];

						thisGridPoint->d3_Position[0] = dx;
						thisGridPoint->d3_Position[1] = dy;
						thisGridPoint->d3_Position[2] = dz;

						int index = getIndex(i_x, i_y, i_z, nNodes[0], nNodes[1], nNodes[2]);
						allGridPoint[index] = thisGridPoint;
					}
				}
			}

			return(allGridPoint);
		}

		std::string getScript(GridPoint *thisGridPoint)
		{
			std::stringstream Stream;
			Stream.clear();

			if(thisGridPoint == NULL)
			{// to get the header column names
				Stream << "ID" << "\t";
				Stream << "mass_x" << "\t\t\t";
				Stream << "mass_y" << "\t\t\t";
				Stream << "mass_z" << "\t\t\t";
				Stream << "fixed_x" << "\t";
				Stream << "fixed_y" << "\t";
				Stream << "fixed_z" << "\t";
				Stream << "position_x" << "\t\t";
				Stream << "position_y" << "\t\t";
				Stream << "position_z" << "\t";
				Stream << std::endl;
			}
			else
			{
				Stream << Script(0) << "\t";
//				Stream << Script(thisGridPoint->dMass, 8) << "\t";
				Stream << Script(thisGridPoint->d3_Mass[0], 8) << "\t";
				Stream << Script(thisGridPoint->d3_Mass[1], 8) << "\t";
				Stream << Script(thisGridPoint->d3_Mass[2], 8) << "\t";
				Stream << Script(thisGridPoint->b3_Fixed[0]) << "\t";
				Stream << Script(thisGridPoint->b3_Fixed[1]) << "\t";
				Stream << Script(thisGridPoint->b3_Fixed[2]) << "\t";
				Stream << Script(thisGridPoint->d3_Position[0], 8) << "\t";
				Stream << Script(thisGridPoint->d3_Position[1], 8) << "\t";
				Stream << Script(thisGridPoint->d3_Position[2], 8) << "\t";
				Stream << std::endl;
			}

			return (Stream.str ());
		}

		void saveLatex(std::vector<GridPoint *> allGridPoint, std::string strPrefix, std::string strPostfix)
		{
			std::string strFileName = strPrefix + "GridPoint" + strPostfix + ".txt";
			std::ofstream OutputFile(strFileName.c_str(), std::ios_base::out);

			OutputFile << GridPoint_Factory::getScript(NULL);

			for(unsigned int index = 0; index < allGridPoint.size(); index+=1)
			{
				OutputFile << GridPoint_Factory::getScript(allGridPoint[index]);
			}

			OutputFile.close();
		}
	protected:

	private:
};

#endif // GRIDPOINT_FACTORY_H
