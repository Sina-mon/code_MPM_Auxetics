#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_Classic_SinglePass(double dTimeIncrement_Total)
{
	const int max_num_Threads = 1; // sina, only use 1 for now as you haven;t implemented atomics
	omp_set_num_threads(max_num_Threads);

	clock_t clockCurrent_Total;
	clock_t clockCurrent;

	// locks on grid points for atomic operations
//	std::vector<omp_lock_t *> vGridPoint_Lock;
//	vGridPoint_Lock.resize(allGridPoint.size());
//	for(int index = 0; index < vGridPoint_Lock.size(); index++)
//	{
//		vGridPoint_Lock[index] = new omp_lock_t;
//		omp_init_lock(vGridPoint_Lock[index]);
//	}

	double dDebug_ContactCutoff = -1000.0;

	double dTimeIncrement_Accumulated = 0.0;
	while(dTimeIncrement_Accumulated < dTimeIncrement_Total)
	{
		// default time-increment value
		double dTimeIncrement = d_TimeIncrement_Maximum;
		// unless reaching the end, where the remaining increment is less than the maximum
		if(d_TimeIncrement_Maximum > (dTimeIncrement_Total-dTimeIncrement_Accumulated))
			dTimeIncrement = dTimeIncrement_Total - dTimeIncrement_Accumulated;
		// calculate increments accumulated
		dTimeIncrement_Accumulated += dTimeIncrement;

		clockCurrent_Total = clock();
		#pragma omp parallel
		{
			GridPoint_Mediator	GP_Mediator_Thread;
			Bases BasesFunction_Thread;

			int nThreads = omp_get_num_threads();
			int	iThread = omp_get_thread_num();

			#pragma omp barrier
			// reset grid points ---------------------------------------------- reset grid points
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active == true)
				{
					thisGP->b_Active = false;
					thisGP->d3_Mass = {0.0, 0.0, 0.0};
					thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
					thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
				}
			}

			// reset grid kernel points --------------------------------------- reset grid kernel points

			#pragma omp barrier
			// find adjacent grid points for every material point ------------- find adjacent grid points
			// and the corresponding bases values -----------------------------
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3_Position, i3_Nodes, d3_Length_Cell);

				for(int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[GP_Mediator_Thread.v_adjacentGridPoints[index_AGP]];
					thisAGP->b_Active = true;

					// adjacent grid points
					v_MP_AGP[index_MP][index_AGP].index = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];

					// shape value and shape gradient value
					BasesFunction_Thread.calculateBases(thisMP->d3_Position, thisAGP->d3_Position, d3_Length_Cell);

					v_MP_AGP[index_MP][index_AGP].dShapeValue = BasesFunction_Thread.d_ShapeValue;
					v_MP_AGP[index_MP][index_AGP].d3ShapeGradient = BasesFunction_Thread.d3_ShapeGradient;
				}
			}

			// grid contact kernel -------------------------------------------- GP contact kernel
			// grid contact kernel gradient ----------------------------------- GP contact kernel gradient
			// GP contact kernel (from grid values) --------------------------- GP contact kernel and gradient
			// material point contact kernel (from grid values) --------------- MP contact kernel and gradient
			// detect contact ------------------------------------------------- detect contacts

			#pragma omp barrier
			// material point to grid, only mass ------------------------------ MP to GP (mass)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < v_MP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = v_MP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];

					double dShapeValue = v_MP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = v_MP_AGP[index_MP][index_AGP].d3ShapeGradient;

					// mass
					thisAGP->d3_Mass += dShapeValue * thisMP->d3_Mass;
				}
			}

			#pragma omp barrier
			// material point to grid, velocity and force---------------------- MP to GP (velocity and force)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < v_MP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = v_MP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];

					double dShapeValue = v_MP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = v_MP_AGP[index_MP][index_AGP].d3ShapeGradient;

					// velocity
					if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
						thisAGP->d3_Velocity += dShapeValue * (thisMP->d3_Mass * thisMP->d3_Velocity) / thisAGP->d3_Mass;

					// internal forces
					double dVolume = thisMP->d_Volume;
					thisAGP->d3_Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6_Stress[0] + d3ShapeGradient.y*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[5]);
					thisAGP->d3_Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6_Stress[1] + d3ShapeGradient.x*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[4]);
					thisAGP->d3_Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6_Stress[2] + d3ShapeGradient.x*thisMP->d6_Stress[5] + d3ShapeGradient.y*thisMP->d6_Stress[4]);

					// external forces
					thisAGP->d3_Force += dShapeValue*thisMP->d3_Force_External;
				}
			}

			#pragma omp barrier
			// update grid momentum and apply boundary conditions ------------- update GP momentum
			// and damping ----------------------------------------------------
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active)
				{
					if(glm::length(thisGP->d3_Velocity) > 1.0e-16)
						thisGP->d3_Force -= d_DampingCoefficient * glm::length(thisGP->d3_Force) * glm::normalize(thisGP->d3_Velocity);

					if(glm::length(thisGP->d3_Mass) > d_Mass_Minimum)
						thisGP->d3_Velocity += thisGP->d3_Force / thisGP->d3_Mass * dTimeIncrement;

					if(thisGP->b3_Fixed.x == true)
					{
						thisGP->d3_Velocity.x = 0.0;
						thisGP->d3_Force.x = 0.0;
					}
					if(thisGP->b3_Fixed.y == true)
					{
						thisGP->d3_Velocity.y = 0.0;
						thisGP->d3_Force.y = 0.0;
					}
					if(thisGP->b3_Fixed.z == true)
					{
						thisGP->d3_Velocity.z = 0.0;
						thisGP->d3_Force.z = 0.0;
					}
				}
			}

			#pragma omp barrier
			// sina, not thread safe
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement_Control.size(); index_MP++)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement_Control[index_MP];

				GP_Mediator_Thread.findNeighborGridPoints(thisMP->d3_Position, i3_Nodes, d3_Length_Cell, 0);

				for(int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					unsigned int index_GP = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];
					GridPoint *thisAGP = allGridPoint[index_GP];

					thisAGP->d3_Velocity = thisMP->d3_Velocity;
					thisAGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
				}

//				thisMP->d3_Position += thisMP->d3_Velocity * dTimeIncrement;
			}

			#pragma omp barrier
			// grid to material ----------------------------------------------- GP to MP
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				if(thisMP->b_DisplacementControl == true)
					continue;

				glm::dmat3 d33VelocityGradient = glm::dmat3(0.0);

				for(unsigned int index_AGP = 0; index_AGP < v_MP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[v_MP_AGP[index_MP][index_AGP].index];

					double dShapeValue = v_MP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = v_MP_AGP[index_MP][index_AGP].d3ShapeGradient;

					// velocity
					if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
						thisMP->d3_Velocity += dShapeValue * (thisAGP->d3_Force/thisAGP->d3_Mass) * dTimeIncrement;

					// position
					thisMP->d3_Position += dShapeValue * (thisAGP->d3_Velocity) * dTimeIncrement;

					// velocity gradient, to be used to calculate strains
					d33VelocityGradient += glm::outerProduct(thisAGP->d3_Velocity, d3ShapeGradient);// this glm function does the pre-transposition that we want
				}

				thisMP->d33_DeformationGradient += (d33VelocityGradient * thisMP->d33_DeformationGradient) * dTimeIncrement;

				double dDet = glm::determinant(thisMP->d33_DeformationGradient);
				thisMP->d_Volume = dDet * thisMP->d_Volume_Initial;

				glm::dmat3 d33DeformationGradientIncrement = glm::dmat3(1.0) + d33VelocityGradient * dTimeIncrement;

				double d6StrainIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				d6StrainIncrement[0] = d33DeformationGradientIncrement[0][0] - 1.0;
				d6StrainIncrement[1] = d33DeformationGradientIncrement[1][1] - 1.0;
				d6StrainIncrement[2] = d33DeformationGradientIncrement[2][2] - 1.0;
				d6StrainIncrement[3] = d33DeformationGradientIncrement[0][1] + d33DeformationGradientIncrement[1][0];
				d6StrainIncrement[4] = d33DeformationGradientIncrement[1][2] + d33DeformationGradientIncrement[2][1];
				d6StrainIncrement[5] = d33DeformationGradientIncrement[2][0] + d33DeformationGradientIncrement[0][2];

				double d6StrainRate[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				for(int index = 0; index < 6; index++)
					d6StrainRate[index] = d6StrainIncrement[index] / dTimeIncrement;

				for(int index = 0; index < 6; index++)
					thisMP->d6_Strain[index] += d6StrainIncrement[index];

				// elastic
				double dE = thisMP->d_ElasticModulus;
				// plastic
				double dNu = thisMP->d_PoissonRatio;
				double dYield = thisMP->d_YieldStress;

				double d6StressIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				double d6PlasticStrainIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

				ConstitutiveRelation vonMises_Thread;

				if(thisMP->i_MaterialType == _ELASTIC)
					vonMises_Thread.calculateIncrement_Elastic(dE, dNu, d6StrainIncrement);
				else if(thisMP->i_MaterialType == _PLASTIC)
					vonMises_Thread.calculateIncrement_PerfectlyPlastic_6D(dE, dNu, dYield, thisMP->d6_Stress, d6StrainIncrement);
				else
					vonMises_Thread.calculateIncrement_PerfectlyPlastic_6D(dE, dNu, dYield, thisMP->d6_Stress, d6StrainIncrement);

				for(int index = 0; index < 6; index++)
					d6StressIncrement[index] = vonMises_Thread.d6StressIncrement[index];

				for(int index = 0; index < 6; index++)
					d6PlasticStrainIncrement[index] = vonMises_Thread.d6PlasticStrainIncrement[index];

				for(int index = 0; index < 6; index++)
					thisMP->d6_Stress[index] += d6StressIncrement[index];

				for(int index = 0; index < 6; index++)
					thisMP->d6_Strain_Plastic[index] += d6PlasticStrainIncrement[index];
			}

			#pragma omp barrier
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement_Control.size(); index_MP++)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement_Control[index_MP];

				thisMP->d3_Position += thisMP->d3_Velocity * dTimeIncrement;
			}
		}

		d_Runtime_Total += double(double(clock()-clockCurrent_Total)/CLOCKS_PER_SEC);
		//report to console ---------------------------------------------------
		if(dTime - dTimeConsole_Last > dTimeConsole_Interval)
		{
			dTimeConsole_Last = dTime;
			this->reportConsole();
		}
		//save to latex dataset files ---------------------------------------------
		if(dTime - dTimeLatex_LastSave > dTimeLatex_Interval)
		{
			dTimeLatex_LastSave = dTime;
			this->saveLatex();
		}

		dTime += dTimeIncrement;
		iTimeCycle++;
	}

	// destroy GridPoint locks
//	for(int index = 0; index < vGridPoint_Lock.size(); index++)
//		omp_destroy_lock(vGridPoint_Lock[index]);

	if(dTime < dTimeEnd)
		return(0);
	else
	{
		std::cout << "Analysis complete ..." << std::endl;
		return(1);
	}
}
// ----------------------------------------------------------------------------
