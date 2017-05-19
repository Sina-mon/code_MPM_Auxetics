#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_DropWeight(void)
{
	MaterialPoint_Factory	MP_Factory;
	GridPoint_Factory		GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	{// initialize GP mediator
		d3_Length_Grid = glm::dvec3(0.2, 0.5, 0.2);
		i3_Cells = glm::ivec3(80, 200, 80);

		d3_Length_Cell = d3_Length_Grid / glm::dvec3(i3_Cells);

		i3_Nodes = i3_Cells + glm::ivec3(1, 1, 1);
	}
	allGridPoint = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);

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
			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
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

	// locks on grid points for atomic operations
	v_GridPoint_Lock.resize(allGridPoint.size());
	for(int index = 0; index < v_GridPoint_Lock.size(); index++)
	{
		v_GridPoint_Lock[index] = new omp_lock_t;
		omp_init_lock(v_GridPoint_Lock[index]);
	}

	double dDiameter = 0.080;
	double dLength = 0.450;
	double dThickness = 0.002;
	if(true)
	{// tube material points -------------------------------------------------- tube MP
//		double dDiameter = 0.030;
//		double dThickness = 0.001;
//		double dLength = 0.040;

		double dRadius_Outer = 0.5*dDiameter;
		double dRadius_Inner = dRadius_Outer - dThickness;
		double dOffset = 1.0/2.0*dThickness;
		d_Offset = dOffset;

		double dGravity = 0.0;

		glm::dvec3 d3Center = glm::dvec3(0.5,0.5,0.5) * d3_Length_Grid;
		d3Center.y = 0.5*dLength + 0.5*d3_Length_Cell.y;

		std::vector<MaterialPoint *> thisMaterialDomain = MP_Factory.createDomain_Tube_T2(d3Center, glm::dvec3(_PI/2.0,0.0,0.0), dRadius_Outer, dRadius_Inner, dLength, dOffset);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint *thisMP = thisMaterialDomain[index_MP];

			thisMP->i_MaterialType = _PLASTIC;
			thisMP->i_ID = 1;

			thisMP->d_Volume_Initial = dOffset * dOffset * dOffset;
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 7930.0 * thisMP->d_Volume;
			thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

			thisMP->d_ElasticModulus = 210.0e9;
			thisMP->d_Viscosity = 0.0;
			thisMP->d_PoissonRatio = 0.3;
			thisMP->d_YieldStress = 215.0e6;

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
//			if(glm::abs(thisMP->d3_Position.y - 0.999*dLength) < 0.5*d3_Length_Cell.y)
//			if(glm::abs(thisMP->d3_Position.y - 0.1*dLength) < 0.02*dLength)
			if(glm::abs(thisMP->d3_Position.y - 0.0*dLength) < 0.1*dLength)
			{
				thisMP->b_Mark_Stress = true;
				v_MarkedMaterialPoints_Stress.push_back(thisMP);
			}
		}
	}

	if(true)
	{// platen material points ------------------------------------------------ platen MP
		//double dDiamater = 0.0508;
		double dRadius_Outer = 0.5*0.110;
		double dRadius_Inner = 0.0;
		double dOffset = d_Offset;
		double dLength_Real = 3.0;
		double dLength_Plate = 10.0*dOffset;

		glm::dvec3 d3Center = glm::dvec3(0.5,0.5,0.5) * d3_Length_Grid;
		d3Center.y = dLength + 0.5*dLength_Plate + 3.0*d3_Length_Cell.y;

		std::vector<MaterialPoint *> thisMaterialDomain = MP_Factory.createDomain_Tube_T2(d3Center, glm::dvec3(_PI/2.0,0.0,0.0), dRadius_Outer, dRadius_Inner, dLength_Plate, dOffset);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint *thisMP = thisMaterialDomain[index_MP];

			thisMP->i_MaterialType = _ELASTIC;
			thisMP->i_ID = 1;

			thisMP->d_Volume_Initial = dOffset * dOffset * dOffset;
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 7800.0 * thisMP->d_Volume * dLength_Real/dLength_Plate;
			thisMP->d3_Mass = glm::dvec3(dMass, dMass, dMass);

			thisMP->d_ElasticModulus = 210.0e9;
			thisMP->d_Viscosity = 0.0;
			thisMP->d_PoissonRatio = 0.3;
			thisMP->d_YieldStress = 200.0e6;

			thisMP->d3_Velocity = glm::dvec3(0.0, -20.0, 0.0);
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
//				thisMP->d3_Velocity = glm::dvec3(0.0,-5.0,0.0);
				v_MarkedMaterialPoints_Displacement_Monitor.push_back(thisMP);
			}
		}
	}

	// sina, be careful, this requires the number of adjacent grid points to be exactly 8
	v_MP_AGP.resize(allMaterialPoint.size());

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

	d_Mass_Minimum = 1.0e-9;
	d_DampingCoefficient = 0.0;

	dTimeEnd = 10.0;
	d_TimeIncrement_Maximum = 0.2e-8;
	dTimeConsole_Interval = 1.0e-5;

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
		sDescription += "Grid Point count: " + Script(allGridPoint.size()) + "\n";
		sDescription += "Material Point count: " + Script(allMaterialPoint.size()) + "\n";
		sDescription += "Mass: " + Script(d3Mass_Domain.x,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Grid Resolution: (" + Script(i3_Cells.x) + "," + Script(i3_Cells.y) + "," + Script(i3_Cells.z) + ")\n";
		sDescription += "Material Point offset: " + Script(d_Offset,3) + "\n";
		sDescription += "Platen Speed: " + Script(v_MarkedMaterialPoints_Displacement_Monitor[0]->d3_Velocity.y, 3) + " m/s" + "\n";
		sDescription += "Target yield: " + Script(v_MarkedMaterialPoints_Stress[0]->d_YieldStress, 3) + " N/m^2" + "\n";
		sDescription += "Target Modulus: " + Script(v_MarkedMaterialPoints_Stress[0]->d_ElasticModulus, 3) + " N/m^2" + "\n";
		sDescription += "Global Damping: " + Script(d_DampingCoefficient, 3) + "\n";

		sDescription += "Non-slip contact\n";
		sDescription += "Elastic-Perfectly plastic\n";
	}

	this->reportConsole(sDescription);
}
// ----------------------------------------------------------------------------
