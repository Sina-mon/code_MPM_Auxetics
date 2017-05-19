#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_ConstitutiveRelation(void)
{
	MaterialPoint_Factory	MP_Factory;
	GridPoint_Factory		GP_Factory;
//	GridPoint_Mediator	GP_Mediator;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	{// initialize GP mediator
		d3_Length_Grid = glm::dvec3(0.02, 0.004, 0.004);
		i3_Cells = glm::ivec3(20, 4, 4);

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

		if(fabs(dx - 0.0) < dTolerance)
		{
//			thisGridPoint->b3Fixed = glm::bvec3(true, true, true);
			thisGridPoint->b3Fixed.x = true;
//			thisGridPoint->b3Fixed.z = true;
		}
		if(fabs(dx - d3_Length_Grid.x) < dTolerance)
		{
//			thisGridPoint->b3Fixed = glm::bvec3(true, true, true);
			thisGridPoint->b3Fixed.x = true;
		}
		if(fabs(dy - 0.0) < dTolerance)
		{
			thisGridPoint->b3Fixed = glm::bvec3(true, true, true);
		}
		if(fabs(dy - d3_Length_Grid.y) < dTolerance)
		{
			thisGridPoint->b3Fixed = glm::bvec3(true, true, true);
		}
	}

	// ------------------------------------------------------------------------
	// material points --------------------------------------------------------
	double dGravity = -0.0;
	double dOffset = 0.0002;
	double dLength = 0.01;
	glm::dvec3 d3Center = 0.5*d3_Length_Grid;
	d3Center.x = 0.5*dLength;

	std::vector<MaterialPoint *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center, glm::dvec3(dLength, 0.0008, 0.0008), dOffset);
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
	{// assign material point initial values
		MaterialPoint *thisMP = thisMaterialDomain[index_MP];

		thisMP->i_MaterialType = _VONMISESHARDENING;
		thisMP->i_ID = 0;

		thisMP->dVolume_Initial = dOffset * dOffset * dOffset;
		thisMP->dVolume = thisMP->dVolume_Initial;

		double dMass = 2760.0 * thisMP->dVolume;
		thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

		thisMP->dElasticModulus = 70.0e9;
		thisMP->dViscosity = 0.0e4;
		thisMP->dPoissonRatio = 0.3;
		thisMP->dYieldStress = 320.0e6;

		thisMP->dHardening_Isotropic_C0 = 1.0;
		thisMP->dHardening_Isotropic_C1 = 30.0e6;

		thisMP->d3Velocity = glm::dvec3(0.0, 0.0, 0.0);

		thisMP->d3Momentum = thisMP->d3_Mass * thisMP->d3Velocity;

		thisMP->d3Force_External = dGravity * thisMP->d3_Mass;
	}
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
	{// send to allMaterialPoint vector
		allMaterialPoint.push_back(thisMaterialDomain[index_MP]);
	}
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
	{// mark material points
		double dTolerance = 2.0e-4;
		MaterialPoint *thisMP = thisMaterialDomain[index_MP];
		if(fabs(thisMP->d3Position.x - 0.5*dLength) < dTolerance)
		{
			v_MarkedMaterialPoints_Force.push_back(thisMP);
			thisMP->i_ID = 2;
		}
		if(fabs(thisMP->d3Position.x - dLength) < dTolerance)
		{
			v_MarkedMaterialPoints_Displacement.push_back(thisMP);
			thisMP->i_ID = 1;
		}
	}
	// displacement controlled timeline events --------------------------------
	if(true)
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

	std::cout << "Material Points: " << std::endl;
	std::cout << "	count: " << allMaterialPoint.size() << std::endl;
	std::cout << "	mass_x: " << d3Mass_Domain.x << std::endl;
	std::cout << "	mass_y: " << d3Mass_Domain.y << std::endl;
	std::cout << "	mass_z: " << d3Mass_Domain.z << std::endl;
}
// ----------------------------------------------------------------------------
