#ifndef GRIDPOINT_MEDIATOR_H
#define GRIDPOINT_MEDIATOR_H

#include <iostream>
#include <vector>
#include <math.h>

#include "Definitions.h"
#include "GridPoint.h"
#include "GridPoint_Factory.h"

class GridPoint_Mediator
{
	public:
		std::vector<unsigned int> v_adjacentGridPoints;

		GridPoint_Mediator() {;}
		virtual ~GridPoint_Mediator() {;}

		void findAdjacentGridPoints(glm::dvec3 d3Position, glm::ivec3 i3Node_Count, glm::dvec3 d3Length_Cell)
		{
			v_adjacentGridPoints.resize(0);

			glm::ivec3 i3BottomLeftFar = d3Position / d3Length_Cell;

			int index = 0;

			for(int ix = 0; ix < 2; ix++)
			{
				for(int iy = 0; iy < 2; iy++)
				{
					for(int iz = 0; iz < 2; iz++)
					{
						glm::ivec3 i3Node_Index = i3BottomLeftFar + glm::ivec3(ix, iy, iz);

						if(i3Node_Index.x < 0)	continue;
						if(i3Node_Index.y < 0)	continue;
						if(i3Node_Index.z < 0)	continue;

						if(i3Node_Index.x > i3Node_Count.x-1)	continue;
						if(i3Node_Index.y > i3Node_Count.y-1)	continue;
						if(i3Node_Index.z > i3Node_Count.z-1)	continue;

						index = GridPoint_Factory::getIndex(i3Node_Index, i3Node_Count);
						v_adjacentGridPoints.push_back(index);
					}
				}
			}
			// remove duplicates
//			std::sort( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() );
//			v_adjacentGridPoints.erase( std::unique( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() ), v_adjacentGridPoints.end() );
		}

		void findNeighborGridPoints(glm::dvec3 d3Position, glm::ivec3 i3Node_Count, glm::dvec3 d3Length_Cell, int iLayers)
		{// similar to getAdjacentGridPoint (iLayer=0), but can return extra layers of neighboring grid points
			v_adjacentGridPoints.resize(0);

			glm::ivec3 i3BottomLeftFar = d3Position / d3Length_Cell;

			int index = 0;

			for(int ix = -iLayers; ix < 2+iLayers; ix++)
			{
				for(int iy = -iLayers; iy < 2+iLayers; iy++)
				{
					for(int iz = -iLayers; iz < 2+iLayers; iz++)
					{
						glm::ivec3 i3Node_Index = i3BottomLeftFar + glm::ivec3(ix, iy, iz);

						if(i3Node_Index.x < 0)	continue;
						if(i3Node_Index.y < 0)	continue;
						if(i3Node_Index.z < 0)	continue;

						if(i3Node_Index.x > i3Node_Count.x-1)	continue;
						if(i3Node_Index.y > i3Node_Count.y-1)	continue;
						if(i3Node_Index.z > i3Node_Count.z-1)	continue;

						index = GridPoint_Factory::getIndex(i3Node_Index, i3Node_Count);
						v_adjacentGridPoints.push_back(index);
					}
				}
			}
			// remove duplicates
//			std::sort( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() );
//			v_adjacentGridPoints.erase( std::unique( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() ), v_adjacentGridPoints.end() );
		}

	protected:

	private:
};

#endif // GRIDPOINT_MEDIATOR_H
