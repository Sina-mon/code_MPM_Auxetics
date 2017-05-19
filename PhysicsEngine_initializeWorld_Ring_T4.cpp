#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Ring_T4(void)
{
	MaterialPoint_Factory	MP_Factory;
	GridPoint_Factory		GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	{// initialize GP mediator
//		d3_Length_Grid = glm::dvec3(0.05, 0.05, 0.015);
//		i3_Cells = glm::ivec3(100, 100, 0.015*100/0.05);
		d3_Length_Grid = glm::dvec3(0.05, 0.05, 0.05/100.0);
		i3_Cells = glm::ivec3(100, 100, 1);

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
			thisGridPoint->b3Fixed.x = true;
			thisGridPoint->b3Fixed.y = true;
		}
		if(fabs(dz - 0.0) < dTolerance)
		{
			thisGridPoint->b3Fixed.z = true;
		}
		if(fabs(dz - d3_Length_Grid.z) < dTolerance)
		{
			thisGridPoint->b3Fixed.z = true;
		}
	}

	// ------------------------------------------------------------------------
	// material points --------------------------------------------------------
	double dGravity = -0.0;
	double dOffset = d3_Length_Cell.x/4.0;

	// ring -------------------------------------------------------------------
//	double dRadius_Outer = 0.002;
//	double dRadius_Inner = 0.0;
//	double dRadius_Outer = 0.5*0.02850;
//	double dRadius_Inner = 0.5*0.02850 - 0.0009;
	double dRadius_Outer = 0.5*0.02540;
	double dRadius_Inner = 0.5*0.02540 - 0.0009;
//	double dRadius_Outer = 0.5*0.02540;
//	double dRadius_Inner = 0.5*0.02540 - 0.0013;

	glm::dvec3 d3Center = 0.5*d3_Length_Grid;
	d3Center.x = 1.5*d3_Length_Cell.x + dRadius_Outer;

	std::vector<MaterialPoint *> thisMaterialDomain1 = MP_Factory.createDomain_Tube_T4(d3Center, dRadius_Outer, dRadius_Inner, 0.00, dOffset);
	for(unsigned int index_MP = 0; index_MP < thisMaterialDomain1.size(); index_MP++)
	{// assign material point initial values
		MaterialPoint *thisMP = thisMaterialDomain1[index_MP];

		thisMP->i_MaterialType = _PLASTIC;
		thisMP->i_ID = 0;

		// sina, be careful, 1.0/5.0 might not be the volume of the tetrahedra
		thisMP->dVolume_Initial = MP_Factory.getVolume(thisMP);
//		thisMP->dVolume_Initial = 1.0/6.0 * dOffset * dOffset * dOffset;
		thisMP->dVolume = thisMP->dVolume_Initial;

		double dMass = 1.0*2760.0 * thisMP->dVolume;
		// global value ---------------
		d_Mass_Minimum = 0.0001 * dMass;
		// ----------------------------
		thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

		thisMP->dElasticModulus = 70.0e9;
		thisMP->dViscosity = 0.0e4;
		thisMP->dPoissonRatio = 0.3;
		thisMP->dYieldStress = 300.0e6;

		thisMP->dHardening_Isotropic_C0 = 0.0;
		thisMP->dHardening_Isotropic_C1 = 0.0e6;

		thisMP->d3Velocity = glm::dvec3(-72.7, 0.0, 0.0);

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
			thisMP->i_ID = 2;
		}
		if(fabs(thisMP->d3Position.y - 8.0) < dTolerance)
		{
			v_MarkedMaterialPoints_Displacement.push_back(thisMP);
			thisMP->i_ID = 2;
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
