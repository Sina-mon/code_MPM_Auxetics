#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_Classic_DoublePass_Segment(double dTimeIncrement_Total, std::vector<MaterialPoint *> *activeMaterialPoint)
{
	const int max_num_Threads = 1;
	omp_set_num_threads(max_num_Threads);

	clock_t clockCurrent_Total;
	clock_t clockCurrent;

	double dTimeIncrement_Accumulated = 0.0;
//	for(unsigned int iIncrement_Count = 0; iIncrement_Count < nIncrements; iIncrement_Count++)
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
			// reset grid points---------------------------------------------------
			for (unsigned int index_GP = iThread; index_GP < allGridPoint.size() && true; index_GP += nThreads)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				thisGP->d3_Mass = {0.0, 0.0, 0.0};
				thisGP->d3Velocity = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3Momentum = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3Force = glm::dvec3(0.0, 0.0, 0.0);
			}

			#pragma omp barrier
			// material point to grid
			for(unsigned int index_MP = iThread; index_MP < activeMaterialPoint->size() && true; index_MP += nThreads)
			{
				MaterialPoint *thisMP = (*activeMaterialPoint)[index_MP];

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3Position, i3_Nodes, d3_Length_Cell);
				for(unsigned int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[GP_Mediator_Thread.v_adjacentGridPoints[index_AGP]];

					double dShapeValue = 0.0;
					glm::dvec3 d3ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);
					BasesFunction_Thread.calculateBases(thisMP->d3Position, thisAGP->d3Position, d3_Length_Cell);
					dShapeValue = BasesFunction_Thread.d_ShapeValue;
					d3ShapeGradient = BasesFunction_Thread.d3_ShapeGradient;

					// mass
					#pragma omp atomic
						thisAGP->d3_Mass.x += dShapeValue * thisMP->d3_Mass.x;
					#pragma omp atomic
						thisAGP->d3_Mass.y += dShapeValue * thisMP->d3_Mass.y;
					#pragma omp atomic
						thisAGP->d3_Mass.z += dShapeValue * thisMP->d3_Mass.z;

					// momentum
					#pragma omp atomic
						thisAGP->d3Momentum.x += dShapeValue * thisMP->d3_Mass.x * thisMP->d3Velocity.x;
					#pragma omp atomic
						thisAGP->d3Momentum.y += dShapeValue * thisMP->d3_Mass.y * thisMP->d3Velocity.y;
					#pragma omp atomic
						thisAGP->d3Momentum.z += dShapeValue * thisMP->d3_Mass.z * thisMP->d3Velocity.z;

					// internal forces
					double dVolume = thisMP->dVolume;
					#pragma omp atomic
						thisAGP->d3Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6Stress[0] + d3ShapeGradient.y*thisMP->d6Stress[3] + d3ShapeGradient.z*thisMP->d6Stress[5]);
					#pragma omp atomic
						thisAGP->d3Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6Stress[1] + d3ShapeGradient.x*thisMP->d6Stress[3] + d3ShapeGradient.z*thisMP->d6Stress[4]);
					#pragma omp atomic
						thisAGP->d3Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6Stress[2] + d3ShapeGradient.x*thisMP->d6Stress[5] + d3ShapeGradient.y*thisMP->d6Stress[4]);

					// external forces
					#pragma omp atomic
						thisAGP->d3Force.x += dShapeValue*thisMP->d3Force_External.x;
					#pragma omp atomic
						thisAGP->d3Force.y += dShapeValue*thisMP->d3Force_External.y;
					#pragma omp atomic
						thisAGP->d3Force.z += dShapeValue*thisMP->d3Force_External.z;
				}
			}

			#pragma omp barrier
			// displacement controlled material points ------------------------
			for(unsigned int index_MP = iThread; index_MP < v_MarkedMaterialPoints_Displacement.size() && true; index_MP += nThreads)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement[index_MP];

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3Position, i3_Nodes, d3_Length_Cell);

				for(unsigned int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[GP_Mediator_Thread.v_adjacentGridPoints[index_AGP]];

					thisAGP->d3Velocity = thisMP->d3Velocity;
					thisAGP->d3Momentum = thisAGP->d3_Mass * thisMP->d3Velocity;
					thisAGP->d3Force = glm::dvec3(0.0, 0.0, 0.0);
				}
			}

			#pragma omp barrier
			// update grid momentum and apply boundary conditions -------------
			// and damping ----------------------------------------------------
			for(unsigned int index_GP = iThread; index_GP < allGridPoint.size() && true; index_GP += nThreads)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(glm::length(thisGP->d3Velocity) > 1.0e-16)
					thisGP->d3Force -= d_DampingCoefficient * glm::length(thisGP->d3Force) * glm::normalize(thisGP->d3Velocity);

				thisGP->d3Momentum += thisGP->d3Force * dTimeIncrement;

				thisGP->d3Force_Temp = thisGP->d3Force;

				if(thisGP->b3Fixed.x == true)
				{
					thisGP->d3Velocity.x = 0.0;
					thisGP->d3Momentum.x = 0.0;
					thisGP->d3Force.x = 0.0;
				}
				if(thisGP->b3Fixed.y == true)
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
			// grid to material pass 1 --------------------------------------------
			for(unsigned int index_MP = iThread; index_MP < activeMaterialPoint->size() && true; index_MP += nThreads)
			{
				MaterialPoint *thisMP = (*activeMaterialPoint)[index_MP];

				if(thisMP->b_DisplacementControl == true) // this is displacement controlled material point and its velocity should not be updated
					continue;

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3Position, i3_Nodes, d3_Length_Cell);
				for(unsigned int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[GP_Mediator_Thread.v_adjacentGridPoints[index_AGP]];
					double dShapeValue = 0.0;
					glm::dvec3 d3ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);
					BasesFunction_Thread.calculateBases(thisMP->d3Position, thisAGP->d3Position, d3_Length_Cell);
					dShapeValue = BasesFunction_Thread.d_ShapeValue;
					d3ShapeGradient = BasesFunction_Thread.d3_ShapeGradient;

					if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
						thisMP->d3Velocity += dShapeValue * (thisAGP->d3Force/thisAGP->d3_Mass) * dTimeIncrement;
				}
			}

			#pragma omp barrier
			// map particle velocity/momenta back to grid -------------------------
			// mass in NOT mapped here --------------------------------------------
			for(unsigned int index_MP = iThread; index_MP < activeMaterialPoint->size() && true; index_MP += nThreads)
			{
				MaterialPoint *thisMP = (*activeMaterialPoint)[index_MP];

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3Position, i3_Nodes, d3_Length_Cell);
				for(unsigned int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[GP_Mediator_Thread.v_adjacentGridPoints[index_AGP]];
					double dShapeValue = 0.0;
					glm::dvec3 d3ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);
					BasesFunction_Thread.calculateBases(thisMP->d3Position, thisAGP->d3Position, d3_Length_Cell);
					dShapeValue = BasesFunction_Thread.d_ShapeValue;
					d3ShapeGradient = BasesFunction_Thread.d3_ShapeGradient;

					if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
					{
						#pragma omp atomic
							thisAGP->d3Velocity.x += dShapeValue * (thisMP->d3_Mass.x * thisMP->d3Velocity.x) / thisAGP->d3_Mass.x;
						#pragma omp atomic
							thisAGP->d3Velocity.y += dShapeValue * (thisMP->d3_Mass.y * thisMP->d3Velocity.y) / thisAGP->d3_Mass.y;
						#pragma omp atomic
							thisAGP->d3Velocity.z += dShapeValue * (thisMP->d3_Mass.z * thisMP->d3Velocity.z) / thisAGP->d3_Mass.z;
					}
				}
			}

			#pragma omp barrier
			// displacement controlled material points ------------------------
			for(unsigned int index_MP = iThread; index_MP < v_MarkedMaterialPoints_Displacement.size() && true; index_MP += nThreads)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement[index_MP];

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3Position, i3_Nodes, d3_Length_Cell);

				for(unsigned int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[GP_Mediator_Thread.v_adjacentGridPoints[index_AGP]];

					thisAGP->d3Velocity = thisMP->d3Velocity;
					thisAGP->d3Momentum = thisAGP->d3_Mass * thisMP->d3Velocity;
					thisAGP->d3Force = glm::dvec3(0.0, 0.0, 0.0);
				}
			}

			#pragma omp barrier
			// apply boundary conditions --------------------------------------
			for(unsigned int index_GP = iThread; index_GP < allGridPoint.size() && true; index_GP += nThreads)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b3Fixed.x == true && thisGP->d3_Mass.x > 0.0*(*activeMaterialPoint)[0]->d3_Mass.x)
				{
					thisGP->d3Velocity.x = 0.0;
					thisGP->d3Momentum.x = 0.0;
					thisGP->d3Force.x = 0.0;
				}
				if(thisGP->b3Fixed.y == true && thisGP->d3_Mass.x > 0.0*(*activeMaterialPoint)[0]->d3_Mass.x)
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
			for(unsigned int index_MP = iThread; index_MP < activeMaterialPoint->size() && true; index_MP += nThreads)
			{
				MaterialPoint *thisMP = (*activeMaterialPoint)[index_MP];

				if(thisMP->b_DisplacementControl == true) // displacement controlled material points
					continue;

				glm::dvec3 d3Displacement = glm::dvec3(0.0, 0.0, 0.0);
				glm::dmat3 d33DeformationGradientIncrement = glm::dmat3(1.0);

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3Position, i3_Nodes, d3_Length_Cell);

				for(unsigned int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[GP_Mediator_Thread.v_adjacentGridPoints[index_AGP]];

					double dShapeValue = 0.0;
					glm::dvec3 d3ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);
					BasesFunction_Thread.calculateBases(thisMP->d3Position, thisAGP->d3Position, d3_Length_Cell);
					dShapeValue = BasesFunction_Thread.d_ShapeValue;
					d3ShapeGradient = BasesFunction_Thread.d3_ShapeGradient;

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
			for(unsigned int index_MP = iThread; index_MP < v_MarkedMaterialPoints_Displacement.size() && true; index_MP += nThreads)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement[index_MP];

				thisMP->d3Velocity = glm::dvec3(0.0, 0.0, 0.0);

				if(m_TimeLine.v_Time.size() > 2)
					thisMP->d3Velocity = m_TimeLine.getVelocity(dTime);

				thisMP->d3Momentum = thisMP->d3_Mass * thisMP->d3Velocity;
				thisMP->d3Position += thisMP->d3Velocity * dTimeIncrement;
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

	if(dTime < dTimeEnd)
		return(0);
	else
	{
		std::cout << "Analysis complete ..." << std::endl;
		return(1);
	}
}
// ----------------------------------------------------------------------------
