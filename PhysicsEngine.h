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
		void	initializeWorld_Constitutive(void);
		void	initializeWorld_ContactBlock(void);
		void	initializeWorld_DropWeight(void);
		void	initializeWorld_AluminumQuarterTube(void);
		void	initializeWorld_AluminumTube(void);
		void	initializeWorld_Billet(void);
		void	initializeWorld_DiskRoll(void);
		void	initializeWorld_SphereImpact(void);
		void	initializeWorld_ConstitutiveRelation(void);
		void	initializeWorld_Bar(void);
		void	initializeWorld_Ring(void);
		void	initializeWorld_Ring_T4(void);
		void	initializeWorld_TubeCompression(void);
		void	initializeWorld_Cantilever_T4(void);
		void	initializeWorld_SteelAluminum(void);
		void	initializeWorld_AuxeticCell(void);
		void	initializeWorld_AuxeticMesh(void);

		int		runSimulation_Classic_SinglePass(double dTimeIncrement_Total);
		int		runSimulation_Classic_SinglePass_MP(double dTimeIncrement_Total);
		int		runSimulation_Classic_SinglePass_MP_Contact(double dTimeIncrement_Total);

		// graphics interface -------------------------------------------------
		unsigned int 					getCount_MaterialPoint(void) {return(allMaterialPoint.size());}
		const MaterialPoint				*getMaterialPoint(unsigned int index) {return(allMaterialPoint[index]);}
		std::vector<MaterialPoint *>	getMaterialPoints(void) {return(allMaterialPoint);}
		std::vector<GridPoint *>		getGridPoints(void) {return(allGridPoint);}
		std::vector<GridPoint *>		getGridPoints_Kernel(void) {return(v_GridPoint_Kernel);}
		std::vector<std::array<AGPstruct, 8>> v_MP_AGP;

		unsigned int 	getCount_GridPoint(void) {return(allGridPoint.size());}
		const GridPoint	*getGridPoint(unsigned int index) {return(allGridPoint[index]);}
		// --------------------------------------------------------------------

		// function to communicate with outside
		double getTime_Runtime(void) {return(d_Runtime_Total);}
		double getTime_Current(void) {return(dTime);}
		double getTime_End(void) {return(dTimeEnd);}
		double getTime_Increment(void) {return(d_TimeIncrement_Maximum);}

		double d_Offset = 0.0;
		glm::dvec3 d3_Length_Grid = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Length_Cell = glm::dvec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Cells = glm::ivec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Nodes = glm::ivec3(0.0, 0.0, 0.0);

		glm::dvec3 d3_Length_Grid_Kernel = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Length_Cell_Kernel = glm::dvec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Cells_Kernel = glm::ivec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Nodes_Kernel = glm::ivec3(0.0, 0.0, 0.0);
	protected:
		TimeLine m_TimeLine;
		double d_Mass_Minimum = 0.0;
		double d_DampingCoefficient = 0.0;

		double d_Runtime_Total = 0.0;
		std::array<double, 8> a_Runtime;

		std::vector<GridPoint *> allGridPoint;
		std::vector<GridPoint *> v_GridPoint_Kernel;
		std::vector<MaterialPoint *> allMaterialPoint;
		std::vector<MaterialPoint *> v_MarkedMaterialPoints_Displacement_Monitor;
		std::vector<MaterialPoint *> v_MarkedMaterialPoints_Displacement_Control;
		std::vector<MaterialPoint *> v_MarkedMaterialPoints_Stress_Monitor;
		std::vector<MaterialPoint *> v_MarkedMaterialPoints_Force_Monitor;
		std::vector<MaterialPoint *> v_MarkedMaterialPoints_Momentum;

		std::vector<omp_lock_t *> v_GridPoint_Lock;

		double dTime = 0.0; // current simulation time
		double d_TimeIncrement_Maximum = 5.0e-8;
		double dTimeEnd = 10.0;//5.0e-4;
		int iTimeCycle = 0;

		double dTimeConsole_Interval = 500.0*d_TimeIncrement_Maximum;//0.0002*dTimeEnd;
		double dTimeConsole_Last = -1.0e12; // before creation

		double dTimeLatex_Interval = 10.5*dTimeEnd;
		double dTimeLatex_LastSave = 1.0e12; // before creation
		int iTimeLatexCycle = 0;

		void reportConsole(std::string sDescription = "");
		void saveLatex(void);
	private:
};

#endif // PHYSICSENGINE_H
