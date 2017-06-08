#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_Classic_SinglePass_MP(double dTimeIncrement_Total)
{
//	const int max_num_Threads = 1; // sina, only use 1 for now as you haven;t implemented atomics
//	omp_set_num_threads(max_num_Threads);

	omp_set_num_threads(_MAX_N_THREADS);

	clock_t clockCurrent_Total;
	clock_t clockCurrent;

	double dRuntime_MP = 0.0;
	double dRuntime_Block = 0.0;
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
		dRuntime_MP = omp_get_wtime();
		#pragma omp parallel
		{
			GridPoint_Mediator	GP_Mediator_Thread;
			Bases BasesFunction_Thread;

			int nThreads = omp_get_num_threads();
			int	iThread = omp_get_thread_num();

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// reset grid points ---------------------------------------------- reset grid points
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active == true)
				{
					thisGP->b_Active = false;
					thisGP->d_Mass = 0.0;
					thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
					thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
					thisGP->d3_Force_Temp = glm::dvec3(0.0, 0.0, 0.0);

					for(int index_Thread = 0; index_Thread < nThreads; index_Thread++)
					{
						GridPoint *thisGP_Thread = allGridPoint_Thread[index_Thread][index_GP];

						thisGP_Thread->d_Mass = 0.0;
						thisGP_Thread->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
						thisGP_Thread->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
					}
				}
			}
			a_Runtime[0] += omp_get_wtime() - dRuntime_Block;

			// reset grid kernel points --------------------------------------- reset grid kernel points

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
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
			a_Runtime[1] += omp_get_wtime() - dRuntime_Block;

			// grid contact kernel -------------------------------------------- GP contact kernel
			// grid contact kernel gradient ----------------------------------- GP contact kernel gradient
			// GP contact kernel (from grid values) --------------------------- GP contact kernel and gradient
			// material point contact kernel (from grid values) --------------- MP contact kernel and gradient
			// detect contact ------------------------------------------------- detect contacts

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// material point to grid, only mass ------------------------------ MP to GP (mass)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				if(thisMP->b_DisplacementControl == true)
					continue;

				for(unsigned int index_AGP = 0; index_AGP < v_MP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = v_MP_AGP[index_MP][index_AGP].index;
//					GridPoint *thisAGP = allGridPoint[index_GP];
					GridPoint *thisAGP_Thread = allGridPoint_Thread[iThread][index_GP];

					double dShapeValue = v_MP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = v_MP_AGP[index_MP][index_AGP].d3ShapeGradient;

//					if(nThreads > 1)	omp_set_lock(v_GridPoint_Lock[index_GP]);
					{
						// mass
//						thisAGP->d3_Mass += dShapeValue * thisMP->d3_Mass;
						thisAGP_Thread->d_Mass += dShapeValue * thisMP->d_Mass;
					}
//					if(nThreads > 1)	omp_unset_lock(v_GridPoint_Lock[index_GP]);
				}
			}
			#pragma omp barrier
			// accumulate GP-layered values ----------------------------------- GP-layered
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];
				if(thisGP->b_Active == true)
				{
					for(int index_Thread = 0; index_Thread < nThreads; index_Thread++)
					{
						GridPoint *thisGP_Thread = allGridPoint_Thread[index_Thread][index_GP];

						thisGP->d_Mass += thisGP_Thread->d_Mass;
					}
				}
			}
			a_Runtime[2] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// material point to grid, velocity and force---------------------- MP to GP (velocity and force)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				if(thisMP->b_DisplacementControl == true)
					continue;

				for(unsigned int index_AGP = 0; index_AGP < v_MP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = v_MP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];
					GridPoint *thisAGP_Thread = allGridPoint_Thread[iThread][index_GP];

					double dShapeValue = v_MP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = v_MP_AGP[index_MP][index_AGP].d3ShapeGradient;

//					if(nThreads > 1)	omp_set_lock(v_GridPoint_Lock[index_GP]);
					{
						// velocity
						if(thisAGP->d_Mass > d_Mass_Minimum)
							thisAGP_Thread->d3_Velocity += dShapeValue * (thisMP->d_Mass * thisMP->d3_Velocity) / thisAGP->d_Mass;

						// internal forces
						double dVolume = thisMP->d_Volume;
						thisAGP_Thread->d3_Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6_Stress[0] + d3ShapeGradient.y*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[5]);
						thisAGP_Thread->d3_Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6_Stress[1] + d3ShapeGradient.x*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[4]);
						thisAGP_Thread->d3_Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6_Stress[2] + d3ShapeGradient.x*thisMP->d6_Stress[5] + d3ShapeGradient.y*thisMP->d6_Stress[4]);

						// external forces
						thisAGP_Thread->d3_Force += dShapeValue*thisMP->d3_Force_External;
					}
//					if(nThreads > 1)	omp_unset_lock(v_GridPoint_Lock[index_GP]);
				}
			}
			#pragma omp barrier
			// accumulate GP-layered values ----------------------------------- GP-layered
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];
				if(thisGP->b_Active == true)
				{
					for(int index_Thread = 0; index_Thread < nThreads; index_Thread++)
					{
						GridPoint *thisGP_Thread = allGridPoint_Thread[index_Thread][index_GP];

						thisGP->d3_Velocity	+= thisGP_Thread->d3_Velocity;
						thisGP->d3_Force	+= thisGP_Thread->d3_Force;
					}
				}
			}
			a_Runtime[3] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// update grid momentum and apply boundary conditions ------------- update GP momentum
			// and damping ----------------------------------------------------
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active)
				{
					if(glm::length(thisGP->d3_Velocity) > 1.0e-9)
						thisGP->d3_Force -= d_DampingCoefficient * glm::length(thisGP->d3_Force) * glm::normalize(thisGP->d3_Velocity);

					if(thisGP->d_Mass > d_Mass_Minimum)
						thisGP->d3_Velocity += thisGP->d3_Force / thisGP->d_Mass * dTimeIncrement;

					if(thisGP->b3_Fixed.x == true)
					{
						thisGP->d3_Velocity.x = 0.0;
						thisGP->d3_Force.x = 0.0;
					}
					if(thisGP->b3_Fixed.y == true)
					{
						thisGP->d3_Velocity.y = 0.0;
						thisGP->d3_Force_Temp.y += thisGP->d3_Force.y;
						thisGP->d3_Force.y = 0.0;
					}
					if(thisGP->b3_Fixed.z == true)
					{
						thisGP->d3_Velocity.z = 0.0;
						thisGP->d3_Force.z = 0.0;
					}
				}
			}
			a_Runtime[4] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
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

//					if(nThreads > 1)	omp_set_lock(v_GridPoint_Lock[index_GP]);
					{
						thisAGP->d3_Velocity = thisMP->d3_Velocity;
						//thisAGP->d3_Force_Temp += thisAGP->d3_Force;
						thisAGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
					}
//					if(nThreads > 1)	omp_unset_lock(v_GridPoint_Lock[index_GP]);
				}
			}
			a_Runtime[5] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
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
					if(glm::length(thisAGP->d_Mass) > d_Mass_Minimum)
						thisMP->d3_Velocity += dShapeValue * (thisAGP->d3_Force/thisAGP->d_Mass) * dTimeIncrement;

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
				else if(thisMP->i_MaterialType == _VONMISESHARDENING)
					vonMises_Thread.calculateIncrement_VonMisesHardening_6D(dE, dNu, dYield, thisMP->d_BackStress_Isotropic, thisMP->d_Hardening_Isotropic_C0, thisMP->d_Hardening_Isotropic_C1, thisMP->d6_Stress, d6StrainIncrement);
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

				thisMP->d_BackStress_Isotropic += vonMises_Thread.dBackstress_IsotropicIncrement;
			}
			a_Runtime[6] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement_Control.size(); index_MP++)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement_Control[index_MP];

				if(m_TimeLine.v_Time.size() != 0)
					thisMP->d3_Velocity = m_TimeLine.getVelocity(d_Time);
				thisMP->d3_Position += thisMP->d3_Velocity * dTimeIncrement;
			}
			a_Runtime[7] += omp_get_wtime() - dRuntime_Block;

		}

//		d_Runtime_Total += double(double(clock()-clockCurrent_Total)/CLOCKS_PER_SEC);
		d_Runtime_Total += omp_get_wtime() - dRuntime_MP;
		//report to console ---------------------------------------------------
		if(d_Time - d_TimeConsole_Last > d_TimeConsole_Interval)
		{
			d_TimeConsole_Last = d_Time;
			this->reportConsole();
		}

		d_Time += dTimeIncrement;
		i_TimeCycle++;
	}

	if(d_Time < d_TimeEnd)
		return(0);
	else
	{
		d_TimeConsole_Last = d_Time;
		this->reportConsole();
		std::cout << "Analysis complete ..." << std::endl;
		return(1);
	}
}
// ----------------------------------------------------------------------------
