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
		// for multi-threading ----------------------------
		std::vector<unsigned int> v_adjacentGridPoints;
		double d_ShapeValue = 0.0;
		glm::dvec3 d3_ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);

		double getShapeValue(void) {return(d_ShapeValue);}
		glm::dvec3 getShapeGradient(void) {return(d3_ShapeGradient);}
		// ------------------------------------------------
		glm::dvec3 d3_Length_Grid = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Length_Cell = glm::dvec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Cells = glm::ivec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Node_Count = glm::ivec3(0.0, 0.0, 0.0);

		GridPoint_Mediator() {;}
		virtual ~GridPoint_Mediator() {;}

		void findAdjacentGridPoints(glm::dvec3 d3Position)
		{
			v_adjacentGridPoints.resize(0);

			glm::ivec3 i3BottomLeftFar = d3Position / d3_Length_Cell;

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

						if(i3Node_Index.x > i3_Node_Count.x-1)	continue;
						if(i3Node_Index.y > i3_Node_Count.y-1)	continue;
						if(i3Node_Index.z > i3_Node_Count.z-1)	continue;

						index = GridPoint_Factory::getIndex(i3Node_Index, i3_Node_Count);
						v_adjacentGridPoints.push_back(index);
					}
				}
			}
			// remove duplicates
//			std::sort( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() );
//			v_adjacentGridPoints.erase( std::unique( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() ), v_adjacentGridPoints.end() );
		}
		void calculateBases(glm::dvec3 d3Position_MP, glm::dvec3 d3Position_GP)
		{
			glm::dvec3 d3Distance = d3Position_MP - d3Position_GP;

			glm::dvec3 d3ShapeValue = glm::dvec3(0.0, 0.0, 0.0);
			d3ShapeValue.x = 1.0 - fabs(d3Distance.x) / d3_Length_Cell.x;
			d3ShapeValue.y = 1.0 - fabs(d3Distance.y) / d3_Length_Cell.y;
			d3ShapeValue.z = 1.0 - fabs(d3Distance.z) / d3_Length_Cell.z;

			if(d3ShapeValue.x < 0.0 || d3ShapeValue.x > 1.0)
				d3ShapeValue.x = 0.0;

			if(d3ShapeValue.y < 0.0 || d3ShapeValue.y > 1.0)
				d3ShapeValue.y = 0.0;

			if(d3ShapeValue.z < 0.0 || d3ShapeValue.z > 1.0)
				d3ShapeValue.z = 0.0;

			// results to pointers
			this->d_ShapeValue = d3ShapeValue.x * d3ShapeValue.y * d3ShapeValue.z;

			this->d3_ShapeGradient.x = -d3ShapeValue.y * d3ShapeValue.z * glm::sign(d3Distance.x) / d3_Length_Cell.x;//sina, make sure these are formulated correctly
			this->d3_ShapeGradient.y = -d3ShapeValue.x * d3ShapeValue.z * glm::sign(d3Distance.y) / d3_Length_Cell.y;
			this->d3_ShapeGradient.z = -d3ShapeValue.x * d3ShapeValue.y * glm::sign(d3Distance.z) / d3_Length_Cell.z;
		}
//		void findAdjacentGridPoints(glm::dvec3 d3Position, glm::ivec3 i3Node_Count, glm::dvec3 d3Length_Cell)
//		{
//			v_adjacentGridPoints.resize(0);
//
//			glm::ivec3 i3BottomLeftFar = d3Position / d3Length_Cell;
//
//			int index = 0;
//
//			for(int ix = 0; ix < 2; ix++)
//			{
//				for(int iy = 0; iy < 2; iy++)
//				{
//					for(int iz = 0; iz < 2; iz++)
//					{
//						glm::ivec3 i3Node_Index = i3BottomLeftFar + glm::ivec3(ix, iy, iz);
//
//						if(i3Node_Index.x < 0)	continue;
//						if(i3Node_Index.y < 0)	continue;
//						if(i3Node_Index.z < 0)	continue;
//
//						if(i3Node_Index.x > i3Node_Count.x-1)	continue;
//						if(i3Node_Index.y > i3Node_Count.y-1)	continue;
//						if(i3Node_Index.z > i3Node_Count.z-1)	continue;
//
//						index = GridPoint_Factory::getIndex(i3Node_Index, i3Node_Count);
//						v_adjacentGridPoints.push_back(index);
//					}
//				}
//			}
//			// remove duplicates
////			std::sort( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() );
////			v_adjacentGridPoints.erase( std::unique( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() ), v_adjacentGridPoints.end() );
//		}

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
