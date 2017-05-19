
#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_Classic_DoublePass_New(double dTimeIncrement_Total)
{
	const int max_num_Threads = 4;
	omp_set_num_threads(max_num_Threads);

	// sina, important --------------------------------------------------------
	std::vector<GridPoint *> allGridPoint_Thread[max_num_Threads];
	for(int iThread = 0; iThread < max_num_Threads; iThread++)
	{
		allGridPoint_Thread[iThread].resize(allGridPoint.size());

		for(int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = new GridPoint;
			allGridPoint_Thread[iThread][index_GP] = thisGP;
		}
	}
	// ------------------------------------------------------------------------

	clock_t clockCurrent_Total;
	clock_t clockCurrent;

	struct AGPstruct
	{
		unsigned int index = 0;
		double dShapeValue = 0.0;
		glm::dvec3 d3ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);
	};

	std::vector<std::array<AGPstruct, 8>> vMP_AGP;
	vMP_AGP.resize(allMaterialPoint.size());

	double dTimeIncrement_Accumulated = 0.0;
	while(dTimeIncrement_Accumulated < dTimeIncrement_Total)
	{
		clockCurrent_Total = clock();

		// default time-increment value
		double dTimeIncrement = d_TimeIncrement_Maximum;
		// unless reaching the end, where the remaining increment is less than the maximum
		if(d_TimeIncrement_Maximum > (dTimeIncrement_Total-dTimeIncrement_Accumulated))
			dTimeIncrement = dTimeIncrement_Total - dTimeIncrement_Accumulated;
		// calculate increments accumulated
		dTimeIncrement_Accumulated += dTimeIncrement;

		#pragma omp parallel
		{
			GridPoint_Mediator	GP_Mediator_Thread;
			Bases BasesFunction_Thread;

			int nThreads = omp_get_num_threads();
			int	iThread = omp_get_thread_num();

			#pragma omp barrier
			// reset grid points ----------------------------------------------
//			for (unsigned int index_GP = iThread; index_GP < allGridPoint.size() && true; index_GP += nThreads)
			#pragma omp for
			for (unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP += 1)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				thisGP->b_Active = false;

				thisGP->d3_Mass = {0.0, 0.0, 0.0};
				thisGP->d3Velocity = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3Momentum = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3Force = glm::dvec3(0.0, 0.0, 0.0);
			}

			#pragma omp barrier
			// find adjacent grid points for every material point -------------
			// and the corresponding bases values -----------------------------
//			for(unsigned int index_MP = iThread; index_MP < allMaterialPoint.size() && true; index_MP += nThreads)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3Position, i3_Nodes, d3_Length_Cell);

				for(int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[GP_Mediator_Thread.v_adjacentGridPoints[index_AGP]];

					thisAGP->b_Active = true;

					// adjacent grid points
					vMP_AGP[index_MP][index_AGP].index = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];

					// shape value and shape gradient value
					BasesFunction_Thread.calculateBases(thisMP->d3Position, thisAGP->d3Position, d3_Length_Cell);

					vMP_AGP[index_MP][index_AGP].dShapeValue = BasesFunction_Thread.d_ShapeValue;
					vMP_AGP[index_MP][index_AGP].d3ShapeGradient = BasesFunction_Thread.d3_ShapeGradient;
				}
			}

			#pragma omp barrier
			// material point to grid
//			for(unsigned int index_MP = iThread; index_MP < allMaterialPoint.size() && true; index_MP += nThreads)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP_Thread = allGridPoint_Thread[iThread][vMP_AGP[index_MP][index_AGP].index];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

					// mass
					thisAGP_Thread->d3_Mass += dShapeValue * thisMP->d3_Mass;

					// momentum
					thisAGP_Thread->d3Momentum += dShapeValue * (thisMP->d3_Mass * thisMP->d3Velocity);

					// internal forces
					double dVolume = thisMP->dVolume;
					thisAGP_Thread->d3Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6Stress[0] + d3ShapeGradient.y*thisMP->d6Stress[3] + d3ShapeGradient.z*thisMP->d6Stress[5]);
					thisAGP_Thread->d3Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6Stress[1] + d3ShapeGradient.x*thisMP->d6Stress[3] + d3ShapeGradient.z*thisMP->d6Stress[4]);
					thisAGP_Thread->d3Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6Stress[2] + d3ShapeGradient.x*thisMP->d6Stress[5] + d3ShapeGradient.y*thisMP->d6Stress[4]);

					// external forces
					thisAGP_Thread->d3Force += dShapeValue*thisMP->d3Force_External;
				}
			}

			// assemble GridPoint from GridPoint_Threads ----------------------
			#pragma omp barrier
//			for (unsigned int index_GP = iThread; index_GP < allGridPoint.size() && true; index_GP += nThreads)
			#pragma omp for
			for (unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP += 1)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active == true)
				{
					for(unsigned int index_Thread = 0; index_Thread < nThreads; index_Thread++)
					{
						GridPoint *thisGP_Thread = allGridPoint_Thread[index_Thread][index_GP];

						thisGP->d3_Mass		+= thisGP_Thread->d3_Mass;
						//thisGP->d3Velocity	+= thisGP_Thread->d3Velocity;
						thisGP->d3Momentum	+= thisGP_Thread->d3Momentum;
						thisGP->d3Force		+= thisGP_Thread->d3Force;

						thisGP_Thread->d3_Mass = {0.0, 0.0, 0.0};
						//thisGP_Thread->d3Velocity = glm::dvec3(0.0, 0.0, 0.0);
						thisGP_Thread->d3Momentum = glm::dvec3(0.0, 0.0, 0.0);
						thisGP_Thread->d3Force = glm::dvec3(0.0, 0.0, 0.0);
					}
				}
			}

			#pragma omp barrier
			// displacement controlled material points ------------------------
//			for(unsigned int index_MP = iThread; index_MP < v_MarkedMaterialPoints_Displacement.size() && true; index_MP += nThreads)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement[index_MP];

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3Position, i3_Nodes, d3_Length_Cell);

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[vMP_AGP[index_MP][index_AGP].index];

					thisAGP->d3Velocity = thisMP->d3Velocity;
					thisAGP->d3Momentum = thisAGP->d3_Mass * thisMP->d3Velocity;
					thisAGP->d3Force = glm::dvec3(0.0, 0.0, 0.0);
				}
			}

			#pragma omp barrier
			// update grid momentum and apply boundary conditions -------------
			// and damping ----------------------------------------------------
//			for(unsigned int index_GP = iThread; index_GP < allGridPoint.size() && true; index_GP += nThreads)
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP += 1)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(glm::length(thisGP->d3Velocity) > 1.0e-16)
					thisGP->d3Force -= d_DampingCoefficient * glm::length(thisGP->d3Force) * glm::normalize(thisGP->d3Velocity);

				thisGP->d3Momentum += thisGP->d3Force * dTimeIncrement;

				if(thisGP->b3Fixed.x == true)
				{
					thisGP->d3Velocity.x = 0.0;
					thisGP->d3Momentum.x = 0.0;
					//thisGP->d3Force.x = 0.0;
				}
				if(thisGP->b3Fixed.y == true)
				{
					thisGP->d3Velocity.y = 0.0;
					thisGP->d3Momentum.y = 0.0;
					//thisGP->d3Force.y = 0.0;
				}
				if(thisGP->b3Fixed.z == true)
				{
					thisGP->d3Velocity.z = 0.0;
					thisGP->d3Momentum.z = 0.0;
					//thisGP->d3Force.z = 0.0;
				}
			}

			#pragma omp barrier
			// grid to material pass 1 --------------------------------------------
//			for(unsigned int index_MP = iThread; index_MP < allMaterialPoint.size() && true; index_MP += nThreads)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				if(thisMP->b_DisplacementControl == true) // this is displacement controlled material point and its velocity should not be updated
					continue;

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[vMP_AGP[index_MP][index_AGP].index];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

					if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
						thisMP->d3Velocity += dShapeValue * (thisAGP->d3Force/thisAGP->d3_Mass) * dTimeIncrement;
				}
			}

			#pragma omp barrier
			// map particle velocity/momenta back to grid -------------------------
			// mass in NOT mapped here --------------------------------------------
//			for(unsigned int index_MP = iThread; index_MP < allMaterialPoint.size() && true; index_MP += nThreads)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[vMP_AGP[index_MP][index_AGP].index];
					GridPoint *thisAGP_Thread = allGridPoint_Thread[iThread][vMP_AGP[index_MP][index_AGP].index];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

					if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
					{
						thisAGP_Thread->d3Velocity += dShapeValue * (thisMP->d3_Mass * thisMP->d3Velocity) / thisAGP->d3_Mass;

//						thisAGP_Thread->d3Velocity.x += dShapeValue * (thisMP->d3_Mass.x * thisMP->d3Velocity.x) / thisAGP->d3_Mass.x;
//						thisAGP_Thread->d3Velocity.y += dShapeValue * (thisMP->d3_Mass.y * thisMP->d3Velocity.y) / thisAGP->d3_Mass.y;
//						thisAGP_Thread->d3Velocity.z += dShapeValue * (thisMP->d3_Mass.z * thisMP->d3Velocity.z) / thisAGP->d3_Mass.z;
					}
				}
			}

			// assemble GridPoint from GridPoint_Threads ----------------------
			#pragma omp barrier
//			for (unsigned int index_GP = iThread; index_GP < allGridPoint.size() && true; index_GP += nThreads)
			#pragma omp for
			for (unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP += 1)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active == true)
				{
					for(unsigned int index_Thread = 0; index_Thread < nThreads; index_Thread++)
					{
						GridPoint *thisGP_Thread = allGridPoint_Thread[index_Thread][index_GP];

						//thisGP->d3_Mass		+= thisGP_Thread->d3_Mass;
						thisGP->d3Velocity	+= thisGP_Thread->d3Velocity;
						//thisGP->d3Momentum	+= thisGP_Thread->d3Momentum;
						//thisGP->d3Force		+= thisGP_Thread->d3Force;

						//thisGP_Thread->d3_Mass = {0.0, 0.0, 0.0};
						thisGP_Thread->d3Velocity = glm::dvec3(0.0, 0.0, 0.0);
						//thisGP_Thread->d3Momentum = glm::dvec3(0.0, 0.0, 0.0);
						//thisGP_Thread->d3Force = glm::dvec3(0.0, 0.0, 0.0);
					}
				}
			}

			#pragma omp barrier
			// displacement controlled material points ------------------------
//			for(unsigned int index_MP = iThread; index_MP < v_MarkedMaterialPoints_Displacement.size() && true; index_MP += nThreads)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[vMP_AGP[index_MP][index_AGP].index];

					thisAGP->d3Velocity = thisMP->d3Velocity;
					thisAGP->d3Momentum = thisAGP->d3_Mass * thisMP->d3Velocity;
					thisAGP->d3Force = glm::dvec3(0.0, 0.0, 0.0);
				}
			}

			#pragma omp barrier
			// apply boundary conditions --------------------------------------
//			for(unsigned int index_GP = iThread; index_GP < allGridPoint.size() && true; index_GP += nThreads)
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP += 1)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b3Fixed.x == true && thisGP->d3_Mass.x > 0.0*allMaterialPoint[0]->d3_Mass.x)
				{
					thisGP->d3Velocity.x = 0.0;
					thisGP->d3Momentum.x = 0.0;
					thisGP->d3Force.x = 0.0;
				}
				if(thisGP->b3Fixed.y == true && thisGP->d3_Mass.x > 0.0*allMaterialPoint[0]->d3_Mass.x)
				{
					thisGP->d3Velocity.y = 0.0;
					thisGP->d3Momentum.y = 0.0;
					thisGP->d3Force.y = 0.0;
				}
				if(thisGP->b3Fixed.z == true)
				{
					thisGP->d3Velocity.z = 0.0;
					thisGP->d3Momentum.z = 0.0;
					thisGP->d3Force.z = 0.0;
				}
			}

			#pragma omp barrier
			// --------------------------------------------------------------------
			// grid to material pass 2 --------------------------------------------
//			for(unsigned int index_MP = iThread; index_MP < allMaterialPoint.size() && true; index_MP += nThreads)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				if(thisMP->b_DisplacementControl == true) // displacement controlled material points
					continue;

				glm::dvec3 d3Displacement = glm::dvec3(0.0, 0.0, 0.0);
				glm::dmat3 d33DeformationGradientIncrement = glm::dmat3(1.0);

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[vMP_AGP[index_MP][index_AGP].index];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

					if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
						d3Displacement += dShapeValue * (thisAGP->d3Momentum/thisAGP->d3_Mass) * dTimeIncrement;

					glm::dvec3 d3Velocity = thisAGP->d3Velocity;
					d33DeformationGradientIncrement += glm::outerProduct(d3Velocity, d3ShapeGradient) * dTimeIncrement;// sina, this glm function does the transposition that we want
				}

				thisMP->d3Position += d3Displacement;

				thisMP->d33DeformationGradient = d33DeformationGradientIncrement * thisMP->d33DeformationGradient;

				double dDet = glm::determinant(thisMP->d33DeformationGradient);
				thisMP->dVolume = dDet * thisMP->dVolume_Initial;

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
					thisMP->d6Strain[index] += d6StrainIncrement[index];

				// elastic
				double dE = thisMP->dElasticModulus;
				// plastic
				double dNu = thisMP->dPoissonRatio;
				double dYield = thisMP->dYieldStress;
				// hardening
				double db = thisMP->dHardening_Isotropic_C0;
				double dQ = thisMP->dHardening_Isotropic_C1;

				double d6StressIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				double d6PlasticStrainIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

				double dTemperature = 293.0;
				ConstitutiveRelation vonMises_Thread;

				if(thisMP->i_MaterialType == _ELASTIC)
					vonMises_Thread.calculateIncrement_Elastic(dE, dNu, d6StrainIncrement);
				else if(thisMP->i_MaterialType == _VISCOELASTIC)
					vonMises_Thread.calculateIncrement_ViscoElastic_6D(dE, dNu, thisMP->dViscosity, thisMP->d6Stress, thisMP->d6Strain, d6StrainRate);
				else if(thisMP->i_MaterialType == _PLASTIC)
					vonMises_Thread.calculateIncrement_PerfectlyPlastic_6D(dE, dNu, dYield, thisMP->d6Stress, d6StrainIncrement);
				else if(thisMP->i_MaterialType == _VONMISESHARDENING)
					vonMises_Thread.calculateIncrement_VonMisesHardening_6D(dE, dNu, dYield, thisMP->dBackStress_Isotropic, thisMP->dHardening_Isotropic_C0, thisMP->dHardening_Isotropic_C1, thisMP->d6Stress, d6StrainIncrement);
	//			else if(thisMP->i_MaterialType == _GASS)
	//				vonMises_Thread.calculateIncrement_IdealGass(thisMP->dHeatCapacityRatio, thisMP->dSpecificHeat, 0.2, 1.0e-3, thisMP->d6Stress, (thisMP->dMass/thisMP->dVolume), dTemperature, d6StrainRate);
				else
					vonMises_Thread.calculateIncrement_PerfectlyPlastic_6D(dE, dNu, dYield, thisMP->d6Stress, d6StrainIncrement);

				for(int index = 0; index < 6; index++)
					d6StressIncrement[index] = vonMises_Thread.d6StressIncrement[index];

				for(int index = 0; index < 6; index++)
					d6PlasticStrainIncrement[index] = vonMises_Thread.d6PlasticStrainIncrement[index];

				for(int index = 0; index < 6; index++)
					thisMP->d6Stress[index] += d6StressIncrement[index];

				for(int index = 0; index < 6; index++)
					thisMP->d6Strain_Plastic[index] += d6PlasticStrainIncrement[index];

				thisMP->dBackStress_Isotropic += vonMises_Thread.dBackstress_IsotropicIncrement;

				thisMP->d3Momentum = thisMP->d3Velocity * thisMP->d3_Mass;
			}

			#pragma omp barrier
			// displacement controlled material points ------------------------
//			for(unsigned int index_MP = iThread; index_MP < v_MarkedMaterialPoints_Displacement.size() && true; index_MP += nThreads)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement[index_MP];

				thisMP->d3Velocity = glm::dvec3(0.0, 0.0, 0.0);

				if(m_TimeLine.v_Time.size() > 2)
					thisMP->d3Velocity = m_TimeLine.getVelocity(dTime);

				thisMP->d3Momentum = thisMP->d3_Mass * thisMP->d3Velocity;
				thisMP->d3Position += thisMP->d3Velocity * dTimeIncrement;
			}
		}
		dTime += dTimeIncrement;
		iTimeCycle++;

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
	}

	for(int iThread = 0; iThread < max_num_Threads; iThread++)
	{
		for(int index_GP = 0; index_GP < allGridPoint_Thread[iThread].size(); index_GP++)
		{
			GridPoint *thisGP = allGridPoint_Thread[iThread][index_GP];
			delete thisGP;
		}
	}

	if(dTime < dTimeEnd)
		return(0);
	else
	{
		std::cout << "Analysis complete ..." << std::endl;
		return(1);
	}
}
// ----------------------------------------------------------------------------
