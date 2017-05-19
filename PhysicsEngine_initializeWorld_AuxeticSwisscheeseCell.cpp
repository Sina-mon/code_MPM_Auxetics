#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_AuxeticSwisscheeseCell(void)
{
	MaterialPoint_Factory	MP_Factory;
	GridPoint_Factory		GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	{// initialize GP mediator
		d3_Length_Grid = glm::dvec3(0.06, 0.06, 0.001);
		i3_Cells = 1*glm::ivec3(60, 60, 1);

		d3_Length_Cell = d3_Length_Grid / glm::dvec3(i3_Cells);

		i3_Nodes = i3_Cells + glm::ivec3(1, 1, 1);
	}
	allGridPoint = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);

	// contact kernel grid ---------------------------------------------------- contact grid
	{// initialize GP mediator
		d3_Length_Grid_Kernel = d3_Length_Grid;
		i3_Cells_Kernel = 1*i3_Cells;

		d3_Length_Cell_Kernel = d3_Length_Grid_Kernel / glm::dvec3(i3_Cells_Kernel);

		i3_Nodes_Kernel = i3_Cells_Kernel + glm::ivec3(1, 1, 1);
	}
	v_GridPoint_Kernel = GP_Factory.createGrid(d3_Length_Grid_Kernel, i3_Cells_Kernel);

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
			thisGridPoint->b3_Fixed.x = true;
		}
		if(fabs(dx - d3_Length_Grid.x) < dTolerance)
		{
			thisGridPoint->b3_Fixed.x = true;
		}
		if(fabs(dy - 0.0) < dTolerance)
		{
			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
//			thisGridPoint->b3_Fixed.y = true;
		}
		if(fabs(dy - d3_Length_Grid.y) < dTolerance)
		{
//			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
		}
		if(fabs(dz - 0.0) < dTolerance)
		{
//			thisGridPoint->b3_Fixed.z = true;
		}
		if(fabs(dz - d3_Length_Grid.z) < dTolerance)
		{
//			thisGridPoint->b3_Fixed.z = true;
		}
	}

	// locks on grid points for atomic operations
	v_GridPoint_Lock.resize(allGridPoint.size());
	for(int index = 0; index < v_GridPoint_Lock.size(); index++)
	{
		v_GridPoint_Lock[index] = new omp_lock_t;
		omp_init_lock(v_GridPoint_Lock[index]);
	}

	d_Offset = 1.0/4.0*d3_Length_Cell.x;
	glm::dvec3 d3Dimensions_Cell = glm::dvec3(0.05,0.05,1.1*d_Offset);
	glm::dvec3 d3Center_Cell = glm::dvec3(0.5,0.5,0.5)*d3_Length_Grid;
	d3Center_Cell.y = 0.5*d3Dimensions_Cell.y + 1.5*d3_Length_Cell.y;
	glm::dvec2 d2Spacing = glm::dvec2(0.014, 0.014);
	glm::dvec2 d2Radii = glm::dvec2(0.0075, 0.0055);
	if(true)
	{// cell material points -------------------------------------------------- tube MP
		double dOffset = d_Offset;
		double dGravity = 0.0;

		glm::dvec3 d3Center =  d3Center_Cell;
		glm::dvec3 d3Dimensions = d3Dimensions_Cell;

		std::vector<MaterialPoint *> thisMaterialDomain = MP_Factory.createDomain_AuxeticCell_Swisscheese(d3Center, d3Dimensions, d2Spacing, d2Radii, dOffset);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint *thisMP = thisMaterialDomain[index_MP];

//			thisMP->i_MaterialType = _VONMISESHARDENING;
			thisMP->i_MaterialType = _PLASTIC;
//			thisMP->i_MaterialType = _ELASTIC;
			thisMP->i_ID = 1;

			thisMP->d_Volume_Initial = dOffset * dOffset * dOffset;
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 2700.0 * thisMP->d_Volume;
			d_Mass_Minimum = 0.001 * dMass;
			thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

			thisMP->d_ElasticModulus = 70.0e9;
			thisMP->d_Viscosity = 0.0;
			thisMP->d_PoissonRatio = 0.33;
			thisMP->d_YieldStress = 280.0e6;

			thisMP->d_Hardening_Isotropic_C0 = 4.0e+1;
			thisMP->d_Hardening_Isotropic_C1 = 50.0e6;

			thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
//			thisMP->d3_Momentum = thisMP->d3_Mass * thisMP->d3_Velocity;
			thisMP->d3_Force_External = thisMP->d3_Mass * glm::dvec3(0.0, 0.0, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// send to MP vectors
			MaterialPoint *thisMP = thisMaterialDomain[index_MP];
			// all MPs
			allMaterialPoint.push_back(thisMP);
			// moment log
			v_MarkedMaterialPoints_Momentum.push_back(thisMP);
			// mark for stress monitor
		}
	}

	if(true)
	{// top platen material points -------------------------------------------- platen MP
		double dOffset = d_Offset;
		glm::dvec3 d3Center = glm::dvec3(0.5,0.5,0.5) * d3_Length_Grid;
		d3Center.y = d3Center_Cell.y + 0.5*d3Dimensions_Cell.y + 1.0*d3_Length_Cell.y;
		glm::dvec3 d3Dimension = 1.1*d3Dimensions_Cell;
		d3Dimension.y = 1.1*dOffset;

		std::vector<MaterialPoint *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center, d3Dimension, dOffset);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint *thisMP = thisMaterialDomain[index_MP];

			thisMP->i_MaterialType = _ELASTIC;
			thisMP->i_ID = 1;

			thisMP->d_Volume_Initial = dOffset * dOffset * dOffset;
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 7800.0 * thisMP->d_Volume;
			thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

			thisMP->d_ElasticModulus = 210.0e9;
			thisMP->d_Viscosity = 0.0;
			thisMP->d_PoissonRatio = 0.33;
			thisMP->d_YieldStress = 200.0e6;

			thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
//			thisMP->d3_Momentum = thisMP->d3_Mass * thisMP->d3_Velocity;
			thisMP->d3_Force_External = thisMP->d3_Mass * glm::dvec3(0.0, 0.0, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// send to MP vectors
			MaterialPoint *thisMP = thisMaterialDomain[index_MP];
			// all MPs
			allMaterialPoint.push_back(thisMP);
			// displacement control
			if(true)
			{
				thisMP->b_DisplacementControl = true;
				thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
				v_MarkedMaterialPoints_Displacement_Control.push_back(thisMP);
				v_MarkedMaterialPoints_Displacement_Monitor.push_back(thisMP);
			}
		}
	}

	if(true)
	{// bottom platen material points ----------------------------------------- platen MP
		double dOffset = d_Offset;
		glm::dvec3 d3Center = glm::dvec3(0.5,0.5,0.5) * d3_Length_Grid;
		d3Center.y = d3Center_Cell.y - 0.5*d3Dimensions_Cell.y - 1.0*d3_Length_Cell.y;
		glm::dvec3 d3Dimension = 1.1*d3Dimensions_Cell;
		d3Dimension.y = 1.1*dOffset;

		std::vector<MaterialPoint *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center, d3Dimension, dOffset);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint *thisMP = thisMaterialDomain[index_MP];

			thisMP->i_MaterialType = _ELASTIC;
			thisMP->i_ID = 1;

			thisMP->d_Volume_Initial = dOffset * dOffset * dOffset;
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 7800.0 * thisMP->d_Volume;
			thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

			thisMP->d_ElasticModulus = 210.0e9;
			thisMP->d_Viscosity = 0.0;
			thisMP->d_PoissonRatio = 0.33;
			thisMP->d_YieldStress = 200.0e6;

			thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
//			thisMP->d3_Momentum = thisMP->d3_Mass * thisMP->d3_Velocity;
			thisMP->d3_Force_External = thisMP->d3_Mass * glm::dvec3(0.0, 0.0, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// send to MP vectors
			MaterialPoint *thisMP = thisMaterialDomain[index_MP];
			// all MPs
			allMaterialPoint.push_back(thisMP);
			// stress monitor
//			v_MarkedMaterialPoints_Stress_Monitor.push_back(thisMP);
		}
	}

	// sina, be careful, this requires the number of adjacent grid points to be exactly 8
	v_MP_AGP.resize(allMaterialPoint.size());

	if(true)
	{// timeline events -------------------------------------------------------
		m_TimeLine.addTimePoint(0.0, glm::dvec3(0.0, 0.0, 0.0));
		m_TimeLine.addTimePoint(1.0e-5, glm::dvec3(0.0, -10.0, 0.0));
//		m_TimeLine.addTimePoint(2.0e-5, glm::dvec3(0.0, -2.0, 0.0));
//		m_TimeLine.addTimePoint(3.0e-5, glm::dvec3(0.0, -4.0, 0.0));
//		m_TimeLine.addTimePoint(4.0e-5, glm::dvec3(0.0, -8.0, 0.0));
//		m_TimeLine.addTimePoint(5.0e-5, glm::dvec3(0.0, -16.0, 0.0));
//		m_TimeLine.addTimePoint(6.0e-5, glm::dvec3(0.0, -32.0, 0.0));
//		m_TimeLine.addTimePoint(7.0e-5, glm::dvec3(0.0, -64.0, 0.0));
//		m_TimeLine.addTimePoint(8.0e-5, glm::dvec3(0.0, -128.0, 0.0));
//		m_TimeLine.addTimePoint(9.0e-5, glm::dvec3(0.0, -256.0, 0.0));
//		m_TimeLine.addTimePoint(10.0e-5, glm::dvec3(0.0, -512.0, 0.0));
		m_TimeLine.addTimePoint(1.0e6, glm::dvec3(0.0, -10.0, 0.0));
	}

	glm::dvec3 d3Mass_Domain = {0.0, 0.0, 0.0};
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
	{// calculate debug values
		d3Mass_Domain += allMaterialPoint[index_MP]->d3_Mass;
	}

	a_Runtime.fill(0.0);
//	d_Mass_Minimum = 1.0e-9;
	d_DampingCoefficient = 0.0;

	dTimeEnd = 10.0e-4;
	d_TimeIncrement_Maximum = 5.0e-9;
	dTimeConsole_Interval = 1.0e-6;

	std::string sDescription = "";
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time (&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer,80,"%d-%m-%Y %H:%M:%S",timeinfo);
		std::string strTime(buffer);

		sDescription += "-------------------------------------------------------------\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Process started on: " + strTime + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Time increment: " + Script(d_TimeIncrement_Maximum, 6) + "\n";
		sDescription += "Material Point count: " + Script(allMaterialPoint.size()) + "\n";
		sDescription += "Mass: " + Script(d3Mass_Domain.x,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Grid Resolution: (" + Script(i3_Cells.x) + "," + Script(i3_Cells.y) + "," + Script(i3_Cells.z) + ")\n";
		sDescription += "Kernel Resolution: (" + Script(i3_Cells_Kernel.x) + "," + Script(i3_Cells_Kernel.y) + "," + Script(i3_Cells_Kernel.z) + ")\n";
		sDescription += "Material Point offset: " + Script(d_Offset,3) + "\n";
//        if(v_MarkedMaterialPoints_Displacement_Monitor.size() > 0)
//            sDescription += "Platen Speed: " + Script(v_MarkedMaterialPoints_Displacement_Monitor[0]->d3_Velocity.y, 3) + " m/s" + "\n";
		sDescription += "Timeline Speed: " + Script(m_TimeLine.getVelocity(1.0).y, 3) + " m/s" + "\n";
		sDescription += "Yield: " + Script(allMaterialPoint[0]->d_YieldStress, 3) + " N/m^2" + "\n";
		sDescription += "Modulus: " + Script(allMaterialPoint[0]->d_ElasticModulus, 3) + " N/m^2" + "\n";

		sDescription += "Global Damping: " + Script(d_DampingCoefficient, 3) + "\n";
		sDescription += "Non-slip contact\n";
//		sDescription += "Elastic-Perfectly plastic\n";
	}

	dTimeConsole_Last = 0.0;
	this->reportConsole(sDescription);
}
// ----------------------------------------------------------------------------
