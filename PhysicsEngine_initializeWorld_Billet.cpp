#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Billet(void)
{
	MaterialPoint_Factory	MP_Factory;
	GridPoint_Factory		GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	{// initialize GP mediator
		d3_Length_Grid = glm::dvec3(0.5, 0.5, 0.05);
		i3_Cells = glm::ivec3(20, 20, 2);

//		d3_Length_Cell = d3_Length_Grid / glm::dvec3(i3_Cells);
		d3_Length_Cell.x = d3_Length_Grid.x / i3_Cells.x;
		d3_Length_Cell.y = d3_Length_Grid.y / i3_Cells.y;
		d3_Length_Cell.z = d3_Length_Grid.z / i3_Cells.z;

		i3_Nodes = i3_Cells + glm::ivec3(1, 1, 1);
	}
	allGridPoint = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);

	// contact kernel grid ---------------------------------------------------- contact grid
	{// initialize GP mediator
		d3_Length_Grid_ContactKernel = glm::dvec3(0.5, 0.5, 0.05);
		i3_Cells_ContactKernel = glm::ivec3(40, 40, 4);

		d3_Length_Cell_ContactKernel = d3_Length_Grid_ContactKernel / glm::dvec3(i3_Cells_ContactKernel);

		i3_Nodes_ContactKernel = i3_Cells_ContactKernel + glm::ivec3(1, 1, 1);
	}
	v_GridPoint_ContactKernel = GP_Factory.createGrid(d3_Length_Grid_ContactKernel, i3_Cells_ContactKernel);

	for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
	{// assign grid point boundary conditions
		GridPoint *thisGridPoint = allGridPoint[index_GP];

		double dx = thisGridPoint->d3_Position[0];
		double dy = thisGridPoint->d3_Position[1];
		double dz = thisGridPoint->d3_Position[2];
		double dTolerance = 0.1*d3_Length_Cell.x;

		//fixed grid points
		thisGridPoint->b3_Fixed = glm::bvec3(false, false, false);

		if(fabs(dx - 0.0) < dTolerance)
		{
			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
		}
		if(fabs(dx - d3_Length_Grid.x) < dTolerance)
		{
			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
		}
		if(fabs(dy - 0.0) < dTolerance)
		{
			thisGridPoint->b3_Fixed = glm::bvec3(false, true, false);
		}
		if(fabs(dy - d3_Length_Grid.y) < dTolerance)
		{
			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
		}
		if(fabs(dz - 0.0) < dTolerance)
		{
			thisGridPoint->b3_Fixed = glm::bvec3(false, false, true);
		}
		if(fabs(dz - d3_Length_Grid.z) < dTolerance)
		{
			thisGridPoint->b3_Fixed = glm::bvec3(false, false, true);
		}
	}

	if(true)
	{// platen material points ------------------------------------------------
		double dOffset = 1.0/3.0*d3_Length_Cell.x;
		double dGravity = 0.0;

		glm::dvec3 d3Dimension_Target = glm::dvec3(0.6,0.3,1.0) * d3_Length_Grid;
		d3Dimension_Target.z = dOffset;

		glm::dvec3 d3Center_Target = glm::dvec3(0.5,0.5,0.5) * d3_Length_Grid;

		std::vector<MaterialPoint *> thisMaterialDomain_Target = MP_Factory.createDomain_Cuboid_T2(d3Center_Target, d3Dimension_Target, dOffset);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Target.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint *thisMP = thisMaterialDomain_Target[index_MP];

			thisMP->i_MaterialType = _ELASTIC;
			thisMP->i_ID = 10;

			thisMP->d_Volume_Initial = dOffset * dOffset * dOffset;
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 7800.0 * thisMP->d_Volume;
			thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

			thisMP->d_ElasticModulus = 200.0e9;
			thisMP->d_Viscosity = 1.0e4;
			thisMP->d_PoissonRatio = 0.3;
			thisMP->d_YieldStress = 0.001*thisMP->d_ElasticModulus;

			thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			thisMP->d3_Momentum = thisMP->d3_Mass * thisMP->d3_Velocity;
			thisMP->d3_Force_External = thisMP->d3_Mass * glm::dvec3(0.0, -dGravity, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Target.size(); index_MP++)
		{// send to allMaterialPoint vector
			allMaterialPoint.push_back(thisMaterialDomain_Target[index_MP]);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Target.size(); index_MP++)
		{// mark material points
			double dTolerance = 100.0*dOffset;
			MaterialPoint *thisMP = thisMaterialDomain_Target[index_MP];

			// momentum log
			if(true)// all particles
			{
				v_MarkedMaterialPoints_Momentum.push_back(thisMP);
			}

			// disaplacement controlled
			if(0.7*d3_Length_Grid.y - thisMP->d3_Position.y < dTolerance)
			{
				thisMP->i_ID = 1;
				thisMP->b_DisplacementControl = true;
				thisMP->d3_Velocity = glm::dvec3(0.0,-20.0,0.0);
				v_MarkedMaterialPoints_Displacement.push_back(thisMP);
			}

			// Force application
//			if(fabs(thisMP->d3_Position.y - 0.45*d3_Length_Grid.y) < dTolerance)
//			{
//				thisMP->i_ID = 1;
//				thisMP->d3_Force_External = glm::dvec3(0.0,-200.0e6*dOffset*dOffset,0.0);
//				v_MarkedMaterialPoints_Force.push_back(thisMP);
//			}
		}
	}

	if(true)
	{// billet material points ------------------------------------------------
		double dOffset = 1.0/3.0*d3_Length_Cell.x;
		double dGravity = 0.0;

		glm::dvec3 d3Dimension_Target = glm::dvec3(0.3,0.3,1.0) * d3_Length_Grid;
		d3Dimension_Target.z = dOffset;

		glm::dvec3 d3Center_Target = glm::dvec3(0.5,0.15,0.5) * d3_Length_Grid;

		std::vector<MaterialPoint *> thisMaterialDomain_Target = MP_Factory.createDomain_Cuboid_T2(d3Center_Target, d3Dimension_Target, dOffset);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Target.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint *thisMP = thisMaterialDomain_Target[index_MP];

			thisMP->i_MaterialType = _PLASTIC;
			thisMP->i_ID = 10;

			thisMP->d_Volume_Initial = dOffset * dOffset * dOffset;
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 2700.0 * thisMP->d_Volume;
			thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

			thisMP->d_ElasticModulus = 70.0e9;
			thisMP->d_Viscosity = 1.0e4;
			thisMP->d_PoissonRatio = 0.3;
			thisMP->d_YieldStress = 0.001*thisMP->d_ElasticModulus;

			thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			thisMP->d3_Momentum = thisMP->d3_Mass * thisMP->d3_Velocity;
			thisMP->d3_Force_External = thisMP->d3_Mass * glm::dvec3(0.0, -dGravity, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Target.size(); index_MP++)
		{// send to allMaterialPoint vector
			allMaterialPoint.push_back(thisMaterialDomain_Target[index_MP]);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Target.size(); index_MP++)
		{// mark material points
			double dTolerance = dOffset;
			MaterialPoint *thisMP = thisMaterialDomain_Target[index_MP];

			// momentum log
			if(false)// all particles
			{
				v_MarkedMaterialPoints_Momentum.push_back(thisMP);
			}
		}
	}

	// timeline events --------------------------------------------------------
	if(false)
	{
		m_TimeLine.addTimePoint(0.0, glm::dvec3(0.0, 0.0, 0.0));
		m_TimeLine.addTimePoint(0.1, glm::dvec3(0.0001, 0.0, 0.0));
		m_TimeLine.addTimePoint(1.0e6, glm::dvec3(0.0001, 0.0, 0.0));
	}

	glm::dvec3 d3Mass_Domain = {0.0, 0.0, 0.0};
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
	{// calculate debug values
		d3Mass_Domain += allMaterialPoint[index_MP]->d3_Mass;
	}

	d_Mass_Minimum = 1.0e-16;
	d_DampingCoefficient = 0.0;

	dTimeEnd = 10.0;
	d_TimeIncrement_Maximum = 1.0e-8;
	dTimeConsole_Interval = 1.0e-5;

	std::cout << "Time increment: "	<< Script(d_TimeIncrement_Maximum, 6) << std::endl;

	std::cout << "Grid nodes: " << std::endl;
	std::cout << "	count: " << allGridPoint.size() << std::endl;

	std::cout << "Material Points: " << std::endl;
	std::cout << "	count: " << allMaterialPoint.size() << std::endl;
	std::cout << "	mass_x: " << d3Mass_Domain.x << std::endl;
	std::cout << "	mass_y: " << d3Mass_Domain.y << std::endl;
	std::cout << "	mass_z: " << d3Mass_Domain.z << std::endl;
}
// ----------------------------------------------------------------------------
