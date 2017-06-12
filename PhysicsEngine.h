#ifndef PHYSICSENGINE_H
#define PHYSICSENGINE_H

#include <math.h>
#include <vector>
#include<algorithm>
#include <time.h>

#include <omp.h>

#include "Definitions.h"
////#include ".\include\glm\glm.hpp" // windows
//#include "./include/glm/glm.hpp"//sina, glm is a column major implementation

#include "Bases.h"
#include "GridPoint.h"
#include "GridPoint_Factory.h"
#include "GridPoint_Mediator.h"
#include "MaterialPoint.h"
#include "MaterialPoint_Factory.h"
#include "ConstitutiveRelation.h"
#include "TimeLine.h"

#define _MAX_N_THREADS 2

// adjacent grid point struct to calculate AGP data once
struct AGPstruct
{
	unsigned int index = 0;
	double dShapeValue = 0.0;
	glm::dvec3 d3ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);
};

class PhysicsEngine
{
	public:
		PhysicsEngine();
		virtual ~PhysicsEngine();

		void	initializeWorld_AuxeticSwisscheeseCell(void);
		void	initializeWorld_AuxeticPolygonCell(void);

		glm::dvec3 d3_Length_World = glm::dvec3(0.0, 0.0, 0.0);
		// MPM ----------------------------------------------------------------
		GridPoint_Mediator mpm_GP_Mediator_Thread[_MAX_N_THREADS];
		int		runSimulation_Classic_SinglePass_MP(double dTimeIncrement_Total);
		int		runSimulation_Classic_SinglePass_MP_Contact(double dTimeIncrement_Total);
		std::vector<std::array<AGPstruct, 8>> v_MP_AGP;

		double d_Offset = 0.0;
//		glm::dvec3 d3_Length_Grid = glm::dvec3(0.0, 0.0, 0.0);
//		glm::dvec3 d3_Length_Cell = glm::dvec3(0.0, 0.0, 0.0);
//		glm::ivec3 i3_Cells = glm::ivec3(0.0, 0.0, 0.0);
//		glm::ivec3 i3_Nodes = glm::ivec3(0.0, 0.0, 0.0);

		glm::dvec3 d3_Length_Grid_Kernel = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Length_Cell_Kernel = glm::dvec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Cells_Kernel = glm::ivec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Nodes_Kernel = glm::ivec3(0.0, 0.0, 0.0);

		// function to communicate with outside -------------------------------
		double getTime_Runtime(void) {return(d_Runtime_Total);}
		double getTime_Current(void) {return(d_Time);}
		double getTime_End(void) {return(d_TimeEnd);}
		double getTime_Increment(void) {return(d_TimeIncrement_Maximum);}
		// graphics interface -------------------------------------------------
		unsigned int	getCount_MaterialPoint(void) {return(allMaterialPoint.size());}
		unsigned int 	getCount_GridPoint(void) {return(allGridPoint.size());}
		std::vector<MaterialPoint *>	getMaterialPoints(void) {return(allMaterialPoint);}
		std::vector<GridPoint *>		getGridPoints(void) {return(allGridPoint);}
		std::vector<GridPoint *>		getGridPoints_Kernel(void) {return(v_GridPoint_Kernel);}
	protected:
		TimeLine m_TimeLine;
		double d_Mass_Minimum = 0.0;
		double d_DampingCoefficient = 0.0;

		double d_Runtime_Total = 0.0;
		std::array<double, 8> a_Runtime;

		std::vector<GridPoint *> allGridPoint;
		std::vector<GridPoint *> allGridPoint_Thread[_MAX_N_THREADS];
		std::vector<GridPoint *> v_GridPoint_Kernel;
		std::vector<MaterialPoint *> allMaterialPoint;
		std::vector<MaterialPoint *> v_MarkedMaterialPoints_Displacement_Monitor;
		std::vector<MaterialPoint *> v_MarkedMaterialPoints_Displacement_Control;
		std::vector<MaterialPoint *> v_MarkedMaterialPoints_Stress_Monitor;
		std::vector<MaterialPoint *> v_MarkedMaterialPoints_Force_Monitor;
		std::vector<MaterialPoint *> v_MarkedMaterialPoints_Momentum;

		std::vector<omp_lock_t *> v_GridPoint_Lock;

		double d_Time = 0.0; // current simulation time
		double d_TimeIncrement_Maximum = 1.0e-3;
		double d_TimeEnd = 10.0;//5.0e-4;
		int i_TimeCycle = 0;

		double d_TimeConsole_Interval = 500.0*d_TimeIncrement_Maximum;
		double d_TimeConsole_Last = -1.0e12; // before creation

		void reportConsole(std::string sDescription = "");
	private:
};

#endif // PHYSICSENGINE_H
