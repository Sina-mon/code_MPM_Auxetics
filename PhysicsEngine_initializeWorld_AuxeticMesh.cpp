#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_AuxeticMesh(void)
{
	MaterialPoint_Factory	MP_Factory;
	GridPoint_Factory		GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	{// initialize GP mediator
		d3_Length_Grid = glm::dvec3(0.05, 0.05, 0.05);
		i3_Cells = glm::ivec3(100, 100, 100);

		d3_Length_Cell.x = d3_Length_Grid.x / i3_Cells.x;
		d3_Length_Cell.y = d3_Length_Grid.y / i3_Cells.y;
		d3_Length_Cell.z = d3_Length_Grid.z / i3_Cells.z;

		i3_Nodes = i3_Cells + glm::ivec3(1, 1, 1);
	}

	allGridPoint = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);

	for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
	{// assign grid point initial values
		GridPoint *thisGridPoint = allGridPoint[index_GP];

		double dx = thisGridPoint->d3Position[0];
		double dy = thisGridPoint->d3Position[1];
		double dz = thisGridPoint->d3Position[2];
		double dTolerance = 1.0e-4;

		//fixed grid points
		thisGridPoint->b3Fixed = glm::bvec3(false, false, false);

		if(fabs(dy - 0.0) < dTolerance)
		{
			thisGridPoint->b3Fixed.x = true;
			thisGridPoint->b3Fixed.y = true;
			thisGridPoint->b3Fixed.z = true;
		}
		if(fabs(dy - d3_Length_Grid.y) < dTolerance)
		{
		}
	}

	// ------------------------------------------------------------------------
	// material points --------------------------------------------------------
	double dGravity = -0.0;
	double dOffset = d3_Length_Cell.x/4.0;

	// cell -------------------------------------------------------------------
	glm::dvec3 d3Mesh_Origin = glm::dvec3(0.01,0.003,0.5*d3_Length_Grid.z);//0.5*d3_Length_Grid + 1.01*dOffset*glm::dvec3(1.0, 1.0, 1.0);

	glm::dvec3 d3Dimensions_Cell = glm::dvec3(0.005, 0.01, 0.004);//0.2*d3_Length_Grid;

	std::vector<MaterialPoint *> thisMaterialDomain_Cell = MP_Factory.createDomain_AuxeticMesh(d3Mesh_Origin, glm::ivec2(4,8), d3Dimensions_Cell, dOffset);
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Cell.size(); index_MP++)
	{// assign material point initial values
		MaterialPoint *thisMP = thisMaterialDomain_Cell[index_MP];

		thisMP->i_MaterialType = _PLASTIC;
		thisMP->i_ID = 0;

		thisMP->dVolume_Initial = dOffset * dOffset * dOffset;
		thisMP->dVolume = thisMP->dVolume_Initial;

		double dMass = 2760.0 * thisMP->dVolume;
		// global value ---------------
//		d_Mass_Minimum = 0.01 * dMass;
		// ----------------------------
		thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

		thisMP->dElasticModulus = 70.0e9;
		thisMP->dViscosity = 0.0e4;
		thisMP->dPoissonRatio = 0.3;
		thisMP->dYieldStress = 320.0e6;

		thisMP->dHardening_Isotropic_C0 = 0.0;
		thisMP->dHardening_Isotropic_C1 = 0.0;

		thisMP->d3Velocity = glm::dvec3(0.0, 0.0, 0.0);

		thisMP->d3Momentum = thisMP->d3_Mass * thisMP->d3Velocity;

		thisMP->d3Force_External = dGravity * thisMP->d3_Mass;//glm::dvec3(0.0, dGravity * thisMP->dMass, 0.0);
	}
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Cell.size(); index_MP++)
	{// send to allMaterialPoint vector
		allMaterialPoint.push_back(thisMaterialDomain_Cell[index_MP]);
	}
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Cell.size(); index_MP++)
	{// mark material points
	}

	// wall -------------------------------------------------------------------
	glm::dvec3 d3Center_Wall = 0.5*d3_Length_Grid;
	d3Center_Wall.y = 0.04 + 4.0*d3_Length_Cell.y;

	glm::dvec3 d3Dimensions_Wall = 0.8*d3_Length_Grid;
	d3Dimensions_Wall.y = 1.01*dOffset;
	d3Dimensions_Wall.z = 0.005;

	std::vector<MaterialPoint *> thisMaterialDomain_Wall = MP_Factory.createDomain_Cuboid(d3Center_Wall, d3Dimensions_Wall, dOffset);
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Wall.size(); index_MP++)
	{// assign material point initial values
		MaterialPoint *thisMP = thisMaterialDomain_Wall[index_MP];

		thisMP->i_MaterialType = _ELASTIC;
		thisMP->i_ID = 0;

		thisMP->dVolume_Initial = dOffset * dOffset * dOffset;
		thisMP->dVolume = thisMP->dVolume_Initial;

		double dMass = 2760.0 * thisMP->dVolume;
		// global value ---------------
//		d_Mass_Minimum = 0.01 * dMass;
		// ----------------------------
		thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

		thisMP->dElasticModulus = 70.0e6;
		thisMP->dViscosity = 0.0e4;
		thisMP->dPoissonRatio = 0.3;
		thisMP->dYieldStress = 320.0e9;

		thisMP->dHardening_Isotropic_C0 = 0.0;
		thisMP->dHardening_Isotropic_C1 = 0.0;

		thisMP->d3Velocity = glm::dvec3(0.0, 0.0, 0.0);

		thisMP->d3Momentum = thisMP->d3_Mass * thisMP->d3Velocity;

		thisMP->d3Force_External = dGravity * thisMP->d3_Mass;//glm::dvec3(0.0, dGravity * thisMP->dMass, 0.0);
	}
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Wall.size(); index_MP++)
	{// send to allMaterialPoint vector
		allMaterialPoint.push_back(thisMaterialDomain_Wall[index_MP]);
	}
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain_Wall.size(); index_MP++)
	{// mark material points
		MaterialPoint *thisMP = thisMaterialDomain_Wall[index_MP];

		v_MarkedMaterialPoints_Displacement.push_back(thisMP);
		thisMP->i_ID = 1;
	}


	// displacement controlled timeline events --------------------------------
	if(true)
	{
		m_TimeLine.addTimePoint(0.0, glm::dvec3(0.0, 0.0, 0.0));
		m_TimeLine.addTimePoint(1.0e-6, glm::dvec3(0.0, -100.0, 0.0));
		m_TimeLine.addTimePoint(1.0e-4, glm::dvec3(0.0, -100.0, 0.0));
		m_TimeLine.addTimePoint(1.1e-4, glm::dvec3(0.0, 100.0, 0.0));
		m_TimeLine.addTimePoint(1.0e+6, glm::dvec3(0.0, 100.0, 0.0));
	}

	glm::dvec3 d3Mass_Domain = {0.0, 0.0, 0.0};
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
	{// calculate debug values
		d3Mass_Domain += allMaterialPoint[index_MP]->d3_Mass;
	}

	std::cout << "Material Points: " << std::endl;
	std::cout << "	count: " << allMaterialPoint.size() << std::endl;
	std::cout << "	mass_x: " << d3Mass_Domain.x << std::endl;
	std::cout << "	mass_y: " << d3Mass_Domain.y << std::endl;
	std::cout << "	mass_z: " << d3Mass_Domain.z << std::endl;
}
// ----------------------------------------------------------------------------
