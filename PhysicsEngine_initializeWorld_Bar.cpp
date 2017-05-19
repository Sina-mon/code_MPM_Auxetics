#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Bar(void)
{
	MaterialPoint_Factory	MP_Factory;
	GridPoint_Factory		GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	{// initialize GP mediator
		d3_Length_Grid[0] = 0.40;
		d3_Length_Grid[1] = 0.40;
		d3_Length_Grid[2] = 0.40;

		i3_Cells[0] = 20;
		i3_Cells[1] = 20;
		i3_Cells[2] = 20;

		d3_Length_Cell[0] = d3_Length_Grid[0] / i3_Cells[0];
		d3_Length_Cell[1] = d3_Length_Grid[1] / i3_Cells[1];
		d3_Length_Cell[2] = d3_Length_Grid[2] / i3_Cells[2];

		i3_Nodes[0] = i3_Cells[0] + 1;
		i3_Nodes[1] = i3_Cells[1] + 1;
		i3_Nodes[2] = i3_Cells[2] + 1;
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
		thisGridPoint->b3Fixed[0] = false;
		thisGridPoint->b3Fixed[1] = false;
		thisGridPoint->b3Fixed[2] = false;

		if(fabs(dx - 0.0) < dTolerance)
		{
			thisGridPoint->b3Fixed.x = true;
		}
		if(fabs(dx - 10.0) < dTolerance)
		{
			thisGridPoint->b3Fixed.x = true;
		}
		if(fabs(dy - 0.0) < dTolerance)
		{
			thisGridPoint->b3Fixed.y = true;
		}
		if(fabs(dy - 10.0) < dTolerance)
		{
			thisGridPoint->b3Fixed.y = true;
		}
		if(fabs(dz - 0.0) < dTolerance)
		{
			thisGridPoint->b3Fixed.z = true;
		}
		if(fabs(dz - 10.0) < dTolerance)
		{
			thisGridPoint->b3Fixed.z = true;
		}
	}

	// ------------------------------------------------------------------------
	// material points --------------------------------------------------------
	double dGravity = -0.0;
	double dOffset = 0.001;

/*	std::vector<MaterialPoint *> thisMaterialDomain0 = MaterialPoint_Factory::createDomain_Cuboid(glm::dvec3(5.0, 6.0, 5.0), glm::vec3(3.0,0.2,3.0), dOffset);
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain0.size(); index_MP++)
	{// assign material point initial values
		MaterialPoint *thisMP = thisMaterialDomain0[index_MP];

		thisMP->i_MaterialType = _PLASTIC;
		thisMP->i_ID = 0;

		thisMP->dVolume_Initial = dOffset * dOffset * dOffset;
		thisMP->dVolume = thisMP->dVolume_Initial;
		thisMP->dMass = 105.0 * thisMP->dVolume;

		thisMP->dElasticModulus = 1.0e5;
		thisMP->dViscosity = 1.0e4;
		thisMP->dPoissonRatio = 0.3;
		thisMP->dYieldStress = 1.0e24;

		thisMP->d3Velocity += glm::dvec3(0.0, 0.0, 0.0);

		thisMP->d3Momentum = thisMP->dMass * thisMP->d3Velocity;

		thisMP->d3Force_External = glm::dvec3(0.0, dGravity * thisMP->dMass, 0.0);
	}
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain0.size(); index_MP++)
	{// send to allMaterialPoint vector
		allMaterialPoint.push_back(thisMaterialDomain0[index_MP]);
	}
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain0.size(); index_MP++)
	{// mark material points
		MaterialPoint *thisMP = thisMaterialDomain0[index_MP];

		v_MarkedMaterialPoints_Displacement.push_back(thisMP);
		thisMP->i_ID = 0;
		thisMP->dMass = 100.0;
	}
*/

//	std::vector<MaterialPoint *> thisMaterialDomain1 = MaterialPoint_Factory::createDomain_Cuboid(glm::dvec3(5.0, 2.5, 5.0), glm::vec3(1.0,5.0,0.2), dOffset);
	std::vector<MaterialPoint *> thisMaterialDomain1 = MP_Factory.createDomain_Tube(glm::dvec3(0.2, 0.2, 0.2), 0.5*0.02855, 0.5*0.02855-0.0009, 0.002, dOffset);
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain1.size(); index_MP++)
	{// assign material point initial values
		MaterialPoint *thisMP = thisMaterialDomain1[index_MP];

		thisMP->i_MaterialType = _PLASTIC;
		thisMP->i_ID = 1;

		thisMP->dVolume_Initial = dOffset * dOffset * dOffset;
		thisMP->dVolume = thisMP->dVolume_Initial;

		double dMass = 2760.0 * thisMP->dVolume;
		thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

		thisMP->dElasticModulus = 70.0e9;
		thisMP->dViscosity = 1.0e4;
		thisMP->dPoissonRatio = 0.3;
		thisMP->dYieldStress = 300.0e6;

		thisMP->d3Velocity += glm::dvec3(-50.0, 0.0, 0.0);

		thisMP->d3Momentum = thisMP->d3_Mass * thisMP->d3Velocity;

		thisMP->d3Force_External = dGravity * thisMP->d3_Mass;//glm::dvec3(0.0, dGravity * thisMP->dMass, 0.0);
	}
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain1.size(); index_MP++)
	{// send to allMaterialPoint vector
		allMaterialPoint.push_back(thisMaterialDomain1[index_MP]);
	}
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain1.size(); index_MP++)
	{// mark material points
		double dTolerance = 1.0e-1;
		MaterialPoint *thisMP = thisMaterialDomain1[index_MP];
		if(fabs(thisMP->d3Position.y - 4.0) < dTolerance)
		{
			v_MarkedMaterialPoints_Force.push_back(thisMP);
			thisMP->i_ID = 0;
		}
		if(fabs(thisMP->d3Position.y - 8.0) < dTolerance)
		{
			v_MarkedMaterialPoints_Displacement.push_back(thisMP);
			thisMP->i_ID = 0;
			thisMP->d3_Mass.y *= 100.0;
		}
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
