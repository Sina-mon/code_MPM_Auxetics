#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Cheese(void)
{
	MaterialPoint_Factory	MP_Factory;
	GridPoint_Factory		GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	{// initialize GP mediator
		d3_Length_Grid = glm::dvec3(0.05, 0.05, 0.05);
		i3_Cells = glm::ivec3(50, 50, 50);

		d3_Length_Cell.x = d3_Length_Grid.x / i3_Cells.x;
		d3_Length_Cell.y = d3_Length_Grid.y / i3_Cells.y;
		d3_Length_Cell.z = d3_Length_Grid.z / i3_Cells.z;

		i3_Nodes = i3_Cells + glm::ivec3(1, 1, 1);
	}

	allGridPoint = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);

	// contact kernel grid ---------------------------------------------------- contact grid
	{// initialize GP mediator
		d3_Length_Grid_Kernel = d3_Length_Grid;
//		i3_Cells_Kernel = glm::ivec3(100, 200, 100);
		i3_Cells_Kernel = 1*i3_Cells;

		d3_Length_Cell_Kernel = d3_Length_Grid_Kernel / glm::dvec3(i3_Cells_Kernel);

		i3_Nodes_Kernel = i3_Cells_Kernel + glm::ivec3(1, 1, 1);
	}
	v_GridPoint_Kernel = GP_Factory.createGrid(d3_Length_Grid_Kernel, i3_Cells_Kernel);

	// locks on grid points for atomic operations
	v_GridPoint_Lock.resize(allGridPoint.size());
	for(int index = 0; index < v_GridPoint_Lock.size(); index++)
	{
		v_GridPoint_Lock[index] = new omp_lock_t;
		omp_init_lock(v_GridPoint_Lock[index]);
	}

	for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
	{// assign grid point initial values
		GridPoint *thisGridPoint = allGridPoint[index_GP];

		double dx = thisGridPoint->d3_Position[0];
		double dy = thisGridPoint->d3_Position[1];
		double dz = thisGridPoint->d3_Position[2];
		double dTolerance = 1.0e-4;

		//fixed grid points
		thisGridPoint->b3_Fixed = glm::bvec3(false, false, false);

		if(fabs(dy - 0.0) < dTolerance)
		{
			thisGridPoint->b3_Fixed.x = true;
			thisGridPoint->b3_Fixed.y = true;
			thisGridPoint->b3_Fixed.z = true;
		}
		if(fabs(dy - d3_Length_Grid.y) < dTolerance)
		{
		}
	}

	// ------------------------------------------------------------------------
	// material points --------------------------------------------------------
	double dGravity = -0.0;
	double dOffset = d3_Length_Cell.x/4.0;

	for(float fx = 0.0; fx < d3_Length_Grid.x; fx += dOffset)
	{
		for(float fy = 0.0; fy < d3_Length_Grid.y; fy += dOffset)
		{
			for(float fz = 0.0; fz < dOffset; fz += dOffset)
			{
				glm::dvec3 d3Center;
				double dRadius_x;
				double dRadius_y;
				double dCriterion;
				glm::dvec3 d3Radial;
				// ----------------------------
				d3Center = glm::dvec3(0.25,0.25,0.0)*d3_Length_Grid;
				dRadius_x = 0.01;
				dRadius_y = 0.005;

				d3Radial = glm::dvec3(fx,fy,fz) - d3Center;
				dCriterion = glm::pow(d3Radial.x,2.0)/glm::pow(dRadius_x,2.0) + glm::pow(d3Radial.y,2.0)/glm::pow(dRadius_y,2.0);

				if(dCriterion < 1.0)
					continue;
				// ----------------------------
				d3Center = glm::dvec3(0.75,0.25,0.0)*d3_Length_Grid;
				dRadius_x = 0.005;
				dRadius_y = 0.01;

				d3Radial = glm::dvec3(fx,fy,fz) - d3Center;
				dCriterion = glm::pow(d3Radial.x,2.0)/glm::pow(dRadius_x,2.0) + glm::pow(d3Radial.y,2.0)/glm::pow(dRadius_y,2.0);

				if(dCriterion < 1.0)
					continue;
				// ----------------------------
				d3Center = glm::dvec3(0.25,0.75,0.0)*d3_Length_Grid;
				dRadius_x = 0.005;
				dRadius_y = 0.01;

				d3Radial = glm::dvec3(fx,fy,fz) - d3Center;
				dCriterion = glm::pow(d3Radial.x,2.0)/glm::pow(dRadius_x,2.0) + glm::pow(d3Radial.y,2.0)/glm::pow(dRadius_y,2.0);

				if(dCriterion < 1.0)
					continue;
				// ----------------------------
				d3Center = glm::dvec3(0.75,0.75,0.0)*d3_Length_Grid;
				dRadius_x = 0.01;
				dRadius_y = 0.005;

				d3Radial = glm::dvec3(fx,fy,fz) - d3Center;
				dCriterion = glm::pow(d3Radial.x,2.0)/glm::pow(dRadius_x,2.0) + glm::pow(d3Radial.y,2.0)/glm::pow(dRadius_y,2.0);

				if(dCriterion < 1.0)
					continue;

				MaterialPoint *thisMP = MP_Factory.createMaterialPoint(glm::dvec3(fx,fy,fz), dOffset);
				allMaterialPoint.push_back(thisMP);

				thisMP->i_MaterialType = _PLASTIC;
				thisMP->i_ID = 0;

				thisMP->d_Volume_Initial = dOffset * dOffset * dOffset;
				thisMP->d_Volume = thisMP->d_Volume_Initial;

				double dMass = 2760.0 * thisMP->d_Volume;
				// global value ---------------
		//		d_Mass_Minimum = 0.01 * dMass;
				// ----------------------------
				thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

				thisMP->d_ElasticModulus = 70.0e9;
				thisMP->d_Viscosity = 0.0e4;
				thisMP->d_PoissonRatio = 0.3;
				thisMP->d_YieldStress = 320.0e6;

				thisMP->d_Hardening_Isotropic_C0 = 0.0;
				thisMP->d_Hardening_Isotropic_C1 = 0.0;

				thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);

		//		thisMP->d3Momentum = thisMP->d3_Mass * thisMP->d3Velocity;

				thisMP->d3_Force_External = dGravity * thisMP->d3_Mass;//glm::dvec3(0.0, dGravity * thisMP->dMass, 0.0);
			}
		}
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
