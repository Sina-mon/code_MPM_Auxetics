#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_Classic_DoublePass_Contact(double dTimeIncrement_Total)
{
	const int max_num_Threads = 1; // sina, only use 1 for now as you haven;t implemented atomics
	omp_set_num_threads(max_num_Threads);

	clock_t clockCurrent_Total;
	clock_t clockCurrent;

	// adjacent grid point struct to calculate AGP data once
	struct AGPstruct
	{
		unsigned int index = 0;
		double dShapeValue = 0.0;
		glm::dvec3 d3ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);
	};

	// sina, be careful, this requires the number of adjacent grid points to be exactly 8
	std::vector<std::array<AGPstruct, 8>> vMP_AGP;
	vMP_AGP.resize(allMaterialPoint.size());

	// locks on grid points for atomic operations
	std::vector<omp_lock_t *> vGridPoint_Lock;
	vGridPoint_Lock.resize(allGridPoint.size());
	for(int index = 0; index < vGridPoint_Lock.size(); index++)
	{
		vGridPoint_Lock[index] = new omp_lock_t;
		omp_init_lock(vGridPoint_Lock[index]);
	}

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
			// reset grid points ----------------------------------------------
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP += 1)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

//				thisGP->b_Active = false;
				thisGP->b_Contact_Positive = false;
				thisGP->b_Contact_Negative = false;

				thisGP->d3_Mass = {0.0, 0.0, 0.0};
				thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Momentum = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);

				thisGP->d_Kernel_Contact = 0.0;
				thisGP->d3_Kernel_ContactGradient = glm::dvec3(0.0, 0.0, 0.0);
			}

			//#pragma omp barrier
			// find adjacent grid points for every material point -------------
			// and the corresponding bases values -----------------------------
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3_Position, i3_Nodes, d3_Length_Cell);

				for(int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[GP_Mediator_Thread.v_adjacentGridPoints[index_AGP]];

					// adjacent grid points
					vMP_AGP[index_MP][index_AGP].index = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];

					// shape value and shape gradient value
					BasesFunction_Thread.calculateBases(thisMP->d3_Position, thisAGP->d3_Position, d3_Length_Cell);

					vMP_AGP[index_MP][index_AGP].dShapeValue = BasesFunction_Thread.d_ShapeValue;
					vMP_AGP[index_MP][index_AGP].d3ShapeGradient = BasesFunction_Thread.d3_ShapeGradient;
				}
			}

			#pragma omp barrier
			// grid contact kernel --------------------------------------------
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				double dRadius_Effective = glm::length(d3_Length_Cell);

				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				GP_Mediator_Thread.findNeighborGridPoints(thisMP->d3_Position, i3_Nodes, d3_Length_Cell, 2);

				for(int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					unsigned int index_GP = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];
					GridPoint *thisAGP = allGridPoint[index_GP];

					// distance of MP from this AGP
					double dDistance = glm::length(thisMP->d3_Position - thisAGP->d3_Position);
					// find kernel value
					double dDistance_Normalized = dDistance / dRadius_Effective;
					// find kernel gradient

					if(dDistance_Normalized <= 1.0)
					{
						omp_set_lock(vGridPoint_Lock[index_GP]);
						{// atomic grid operations
//							double dShapeValue = 1.0 - 3.0*glm::pow(dDistance_Normalized,2.0) + 2.0*glm::pow(dDistance_Normalized, 3.0);
//							thisAGP->d_Kernel_Contact += dShapeValue;
//
//							thisAGP->d3_Kernel_ContactGradient += d3ShapeGradient * (-6.0*dDistance_Normalized + 6.0*glm::pow(dDistance_Normalized, 2.0)) / dDistance_Normalized;
//							glm::dvec3 d3ShapeGradient = glm::vec3(0.0,0.0,0.0);
//							if(glm::length(thisMP->d3_Position - thisAGP->d3Position) > 1.0e-2 * glm::length(d3_Length_Cell))
//								d3ShapeGradient = (thisMP->d3_Position - thisAGP->d3Position);
//							thisAGP->d3_Kernel_ContactGradient += d3ShapeGradient * (-6.0 + 6.0*dDistance_Normalized);

							double dShapeValue = 1.0 - dDistance_Normalized;
							thisAGP->d_Kernel_Contact += dShapeValue;

							glm::dvec3 d3ShapeGradient = glm::vec3(0.0,0.0,0.0);
							if(glm::length(thisMP->d3_Position - thisAGP->d3_Position) > 1.0e-3 * glm::length(d3_Length_Cell))
								d3ShapeGradient = glm::normalize(thisMP->d3_Position - thisAGP->d3_Position);
							thisAGP->d3_Kernel_ContactGradient += d3ShapeGradient * -1.0;
						}
						omp_unset_lock(vGridPoint_Lock[index_GP]);
					}
				}
			}


			#pragma omp barrier
			// material point contact kernel (from grid values) ---------------
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				thisMP->d_Kernel_Contact = 0.0;
				thisMP->d3_Kernel_ContactGradient = glm::dvec3(0.0, 0.0, 0.0);

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[vMP_AGP[index_MP][index_AGP].index];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

					thisMP->d_Kernel_Contact += dShapeValue * thisAGP->d_Kernel_Contact;
					if(glm::length(thisAGP->d3_Kernel_ContactGradient) > 10.0)
						thisMP->d3_Kernel_ContactGradient += dShapeValue * thisAGP->d3_Kernel_ContactGradient;
				}

//				if(glm::length(thisMP->d3_Kernel_ContactGradient) < 10.0)
//					thisMP->d3_Kernel_ContactGradient = glm::vec3(0.0, 0.0, 0.0);
			}

			#pragma omp barrier
			// material point to grid
			#pragma omp for// schedule(dynamic,200)
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

/*				if(nThreads == 1)
				{ // avoid atomic operations if there is only 1 thread running
					for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
					{
						unsigned int index_GP = vMP_AGP[index_MP][index_AGP].index;

						GridPoint *thisAGP = allGridPoint[index_GP];
						double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
						glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

						// mass
						thisAGP->d3_Mass += dShapeValue * thisMP->d3_Mass;

						// momentum
						thisAGP->d3Momentum += dShapeValue * (thisMP->d3_Mass * thisMP->d3Velocity);

						// internal forces
						double dVolume = thisMP->dVolume;
						thisAGP->d3Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6Stress[0] + d3ShapeGradient.y*thisMP->d6Stress[3] + d3ShapeGradient.z*thisMP->d6Stress[5]);
						thisAGP->d3Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6Stress[1] + d3ShapeGradient.x*thisMP->d6Stress[3] + d3ShapeGradient.z*thisMP->d6Stress[4]);
						thisAGP->d3Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6Stress[2] + d3ShapeGradient.x*thisMP->d6Stress[5] + d3ShapeGradient.y*thisMP->d6Stress[4]);

						// external forces
						thisAGP->d3Force += dShapeValue*thisMP->d3Force_External;
					}
				}
				else
				{ // use atomics
					unsigned int iAGP_Remaining = vMP_AGP[index_MP].size();

					// sina, be careful, this algorithm requires the adjacent grid point array to be exactly of size 8
					// this may not conform to your newer implementations of getAdjacentGridPoints
					bool bAGP_Remaining[8] = {true, true, true, true, true, true, true, true};

					while(iAGP_Remaining > 0)
					{ // finish the job when all AGPs are processed
						for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
						{ // find an AGP that is not processed and is not locked, process it, and mark it as processed
							if(bAGP_Remaining[index_AGP] == true)
							{
								unsigned int index_GP = vMP_AGP[index_MP][index_AGP].index;

								if(omp_test_lock(vGridPoint_Lock[index_GP]) == true)
								{
									GridPoint *thisAGP = allGridPoint[index_GP];
									double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
									glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

									thisAGP->d3_Mass += dShapeValue * thisMP->d3_Mass;

									// momentum
									thisAGP->d3Momentum += dShapeValue * thisMP->d3_Mass * thisMP->d3Velocity;

									// internal forces
									double dVolume = thisMP->dVolume;

									thisAGP->d3Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6Stress[0] + d3ShapeGradient.y*thisMP->d6Stress[3] + d3ShapeGradient.z*thisMP->d6Stress[5]) + dShapeValue*thisMP->d3Force_External.x;
									thisAGP->d3Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6Stress[1] + d3ShapeGradient.x*thisMP->d6Stress[3] + d3ShapeGradient.z*thisMP->d6Stress[4]) + dShapeValue*thisMP->d3Force_External.y;
									thisAGP->d3Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6Stress[2] + d3ShapeGradient.x*thisMP->d6Stress[5] + d3ShapeGradient.y*thisMP->d6Stress[4]) + dShapeValue*thisMP->d3Force_External.z;

									iAGP_Remaining--;
									bAGP_Remaining[index_AGP] = false;

									omp_unset_lock(vGridPoint_Lock[index_GP]);
								}
							}
						}
					}

				}
*/

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = vMP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

					if(nThreads == 1)
					{ // avoid atomic operations if there is only 1 thread running
						// check for contacting particles
						if(glm::length(thisAGP->d3_Kernel_ContactGradient) > 10.0 && glm::length(thisMP->d3_Kernel_ContactGradient) > 10.0)
						{
							if(glm::dot(glm::normalize(thisAGP->d3_Kernel_ContactGradient), glm::normalize(thisMP->d3_Kernel_ContactGradient)) > 0.1)
								thisAGP->b_Contact_Positive = true;

							if(glm::dot(glm::normalize(thisAGP->d3_Kernel_ContactGradient), glm::normalize(thisMP->d3_Kernel_ContactGradient)) < -0.1)
								thisAGP->b_Contact_Negative = true;
						}

						// mass
						thisAGP->d3_Mass += dShapeValue * thisMP->d3_Mass;

						// momentum
						thisAGP->d3_Momentum += dShapeValue * (thisMP->d3_Mass * thisMP->d3_Velocity);

						// internal forces
						double dVolume = thisMP->d_Volume;
						thisAGP->d3_Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6_Stress[0] + d3ShapeGradient.y*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[5]);
						thisAGP->d3_Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6_Stress[1] + d3ShapeGradient.x*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[4]);
						thisAGP->d3_Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6_Stress[2] + d3ShapeGradient.x*thisMP->d6_Stress[5] + d3ShapeGradient.y*thisMP->d6_Stress[4]);

						// external forces
						thisAGP->d3_Force += dShapeValue*thisMP->d3_Force_External;
					}
					else
					{ // use atomics
						// mass
						omp_set_lock(vGridPoint_Lock[index_GP]);
						{// atomic grid operations
							thisAGP->d3_Mass += dShapeValue * thisMP->d3_Mass;

							// momentum
							thisAGP->d3_Momentum += dShapeValue * thisMP->d3_Mass * thisMP->d3_Velocity;

							// internal forces
							double dVolume = thisMP->d_Volume;

							thisAGP->d3_Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6_Stress[0] + d3ShapeGradient.y*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[5]) + dShapeValue*thisMP->d3_Force_External.x;
							thisAGP->d3_Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6_Stress[1] + d3ShapeGradient.x*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[4]) + dShapeValue*thisMP->d3_Force_External.y;
							thisAGP->d3_Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6_Stress[2] + d3ShapeGradient.x*thisMP->d6_Stress[5] + d3ShapeGradient.y*thisMP->d6_Stress[4]) + dShapeValue*thisMP->d3_Force_External.z;
						}
						omp_unset_lock(vGridPoint_Lock[index_GP]);
					}
				}

			}

			#pragma omp barrier
			// displacement controlled material points ------------------------
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement[index_MP];

				GP_Mediator_Thread.findAdjacentGridPoints(thisMP->d3_Position, i3_Nodes, d3_Length_Cell);

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{ // sina, might need mutexes over here
					GridPoint *thisAGP = allGridPoint[vMP_AGP[index_MP][index_AGP].index];

					thisAGP->d3_Velocity = thisMP->d3_Velocity;
					thisAGP->d3_Momentum = thisAGP->d3_Mass * thisMP->d3_Velocity;
					thisAGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
				}
			}

			#pragma omp barrier
			// update grid momentum and apply boundary conditions -------------
			// and damping ----------------------------------------------------
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP += 1)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(glm::length(thisGP->d3_Velocity) > 1.0e-16)
					thisGP->d3_Force -= d_DampingCoefficient * glm::length(thisGP->d3_Force) * glm::normalize(thisGP->d3_Velocity);

				thisGP->d3_Momentum += thisGP->d3_Force * dTimeIncrement;

				if(thisGP->b3_Fixed.x == true)
				{
					thisGP->d3_Velocity.x = 0.0;
					thisGP->d3_Momentum.x = 0.0;
					// sina, not sure if forces should be set to zero or not
					thisGP->d3_Force.x = 0.0;
				}
				if(thisGP->b3_Fixed.y == true)
				{
					thisGP->d3_Velocity.y = 0.0;
					thisGP->d3_Momentum.y = 0.0;
					thisGP->d3_Force.y = 0.0;
				}
				if(thisGP->b3_Fixed.z == true)
				{
					thisGP->d3_Velocity.z = 0.0;
					thisGP->d3_Momentum.z = 0.0;
					thisGP->d3_Force.z = 0.0;
				}
			}

			#pragma omp barrier
			// grid to material pass 1 --------------------------------------------
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
						thisMP->d3_Velocity += dShapeValue * (thisAGP->d3_Force/thisAGP->d3_Mass) * dTimeIncrement;
				}
			}

			#pragma omp barrier
			// map particle velocity/momenta back to grid -------------------------
			// mass in NOT mapped here --------------------------------------------
			#pragma omp for //schedule(dynamic,100)
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = vMP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

					if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
					{
						if(nThreads == 1)
						{ // avoid atomic operations if there is only 1 thread running
							thisAGP->d3_Velocity.x += dShapeValue * (thisMP->d3_Mass.x * thisMP->d3_Velocity.x) / thisAGP->d3_Mass.x;
							thisAGP->d3_Velocity.y += dShapeValue * (thisMP->d3_Mass.y * thisMP->d3_Velocity.y) / thisAGP->d3_Mass.y;
							thisAGP->d3_Velocity.z += dShapeValue * (thisMP->d3_Mass.z * thisMP->d3_Velocity.z) / thisAGP->d3_Mass.z;
						}
						else
						{
//							omp_set_lock(vGridPoint_Lock[index_GP]);
//							{// atomic grid operations
//								thisAGP->d3Velocity += dShapeValue * (thisMP->d3_Mass * thisMP->d3Velocity) / thisAGP->d3_Mass;
//							}
//							omp_unset_lock(vGridPoint_Lock[index_GP]);
							#pragma omp atomic
								thisAGP->d3_Velocity.x += dShapeValue * (thisMP->d3_Mass.x * thisMP->d3_Velocity.x) / thisAGP->d3_Mass.x;
							#pragma omp atomic
								thisAGP->d3_Velocity.y += dShapeValue * (thisMP->d3_Mass.y * thisMP->d3_Velocity.y) / thisAGP->d3_Mass.y;
							#pragma omp atomic
								thisAGP->d3_Velocity.z += dShapeValue * (thisMP->d3_Mass.z * thisMP->d3_Velocity.z) / thisAGP->d3_Mass.z;
						}
					}
				}
			}

			#pragma omp barrier
			// displacement controlled material points ------------------------
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[vMP_AGP[index_MP][index_AGP].index];

					thisAGP->d3_Velocity = thisMP->d3_Velocity;
					thisAGP->d3_Momentum = thisAGP->d3_Mass * thisMP->d3_Velocity;
					thisAGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
				}
			}

			#pragma omp barrier
			// apply boundary conditions --------------------------------------
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP += 1)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b3_Fixed.x == true && thisGP->d3_Mass.x > 0.0*allMaterialPoint[0]->d3_Mass.x)
				{
					thisGP->d3_Velocity.x = 0.0;
					thisGP->d3_Momentum.x = 0.0;
					thisGP->d3_Force.x = 0.0;
				}
				if(thisGP->b3_Fixed.y == true && thisGP->d3_Mass.x > 0.0*allMaterialPoint[0]->d3_Mass.x)
				{
					thisGP->d3_Velocity.y = 0.0;
					thisGP->d3_Momentum.y = 0.0;
					thisGP->d3_Force.y = 0.0;
				}
				if(thisGP->b3_Fixed.z == true)
				{
					thisGP->d3_Velocity.z = 0.0;
					thisGP->d3_Momentum.z = 0.0;
					thisGP->d3_Force.z = 0.0;
				}
			}

			#pragma omp barrier
			// --------------------------------------------------------------------
			// grid to material pass 2 --------------------------------------------
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
						d3Displacement += dShapeValue * (thisAGP->d3_Momentum/thisAGP->d3_Mass) * dTimeIncrement;

					glm::dvec3 d3Velocity = thisAGP->d3_Velocity;
					d33DeformationGradientIncrement += glm::outerProduct(d3Velocity, d3ShapeGradient) * dTimeIncrement;// sina, this glm function does the transposition that we want
				}

				thisMP->d3_Position += d3Displacement;

				thisMP->d33_DeformationGradient = d33DeformationGradientIncrement * thisMP->d33_DeformationGradient;

				double dDet = glm::determinant(thisMP->d33_DeformationGradient);
				thisMP->d_Volume = dDet * thisMP->d_Volume_Initial;

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
				// hardening
				double db = thisMP->d_Hardening_Isotropic_C0;
				double dQ = thisMP->d_Hardening_Isotropic_C1;

				double d6StressIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				double d6PlasticStrainIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

				double dTemperature = 293.0;
				ConstitutiveRelation vonMises_Thread;

				if(thisMP->i_MaterialType == _ELASTIC)
					vonMises_Thread.calculateIncrement_Elastic(dE, dNu, d6StrainIncrement);
				else if(thisMP->i_MaterialType == _VISCOELASTIC)
					vonMises_Thread.calculateIncrement_ViscoElastic_6D(dE, dNu, thisMP->d_Viscosity, thisMP->d6_Stress, thisMP->d6_Strain, d6StrainRate);
				else if(thisMP->i_MaterialType == _PLASTIC)
					vonMises_Thread.calculateIncrement_PerfectlyPlastic_6D(dE, dNu, dYield, thisMP->d6_Stress, d6StrainIncrement);
				else if(thisMP->i_MaterialType == _VONMISESHARDENING)
					vonMises_Thread.calculateIncrement_VonMisesHardening_6D(dE, dNu, dYield, thisMP->d_BackStress_Isotropic, thisMP->d_Hardening_Isotropic_C0, thisMP->d_Hardening_Isotropic_C1, thisMP->d6_Stress, d6StrainIncrement);
	//			else if(thisMP->i_MaterialType == _GASS)
	//				vonMises_Thread.calculateIncrement_IdealGass(thisMP->dHeatCapacityRatio, thisMP->dSpecificHeat, 0.2, 1.0e-3, thisMP->d6Stress, (thisMP->dMass/thisMP->dVolume), dTemperature, d6StrainRate);
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

				thisMP->d3_Momentum = thisMP->d3_Velocity * thisMP->d3_Mass;
			}

			#pragma omp barrier
			// displacement controlled material points ------------------------
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement.size(); index_MP += 1)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement[index_MP];

				thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);

				if(m_TimeLine.v_Time.size() > 2)
					thisMP->d3_Velocity = m_TimeLine.getVelocity(dTime);

				thisMP->d3_Momentum = thisMP->d3_Mass * thisMP->d3_Velocity;
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
	for(int index = 0; index < vGridPoint_Lock.size(); index++)
		omp_destroy_lock(vGridPoint_Lock[index]);

	if(dTime < dTimeEnd)
		return(0);
	else
	{
		std::cout << "Analysis complete ..." << std::endl;
		return(1);
	}
}
// ----------------------------------------------------------------------------
