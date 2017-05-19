#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_Classic_SinglePass_MP_Contact(double dTimeIncrement_Total)
{
	const int max_num_Threads = 2;
	omp_set_num_threads(max_num_Threads);

	clock_t clockCurrent_Total;
	clock_t clockCurrent;

	double dRuntime_MP = 0.0;
	double dRuntime_Block = 0.0;

//	double dDebug_ContactCutoff = 1.0;

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
					thisGP->d3_Mass = {0.0, 0.0, 0.0};
					thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
					thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
					thisGP->d3_Force_Temp = glm::dvec3(0.0, 0.0, 0.0);

					thisGP->b_Contact_Positive = false;// sina, this variable is only for graphic debugging
					thisGP->b_Contact_Negative = false;// sina, this variable is only for graphic debugging

					thisGP->d3_MassGradient	= glm::dvec3(0.0,0.0,0.0);

					thisGP->d3_Mass_Negative			= glm::dvec3(0.0,0.0,0.0);
					thisGP->d3_MassGradient_Negative = glm::dvec3(0.0,0.0,0.0);
					thisGP->d3_Velocity_Negative		= glm::dvec3(0.0,0.0,0.0);
					thisGP->d3_Force_Negative		= glm::dvec3(0.0,0.0,0.0);

					thisGP->i_NearestMP = 0;

					thisGP->d_Kernel = 0.0;
					thisGP->d3_Kernel_Gradient = glm::dvec3(0.0, 0.0, 0.0);
				}
			}
			a_Runtime[0] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// reset grid kernel points --------------------------------------- reset grid kernel points
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < v_GridPoint_Kernel.size(); index_GP++)
			{
				GridPoint *thisGP = v_GridPoint_Kernel[index_GP];

				if(thisGP->b_Active == true)
				{
					thisGP->b_Active = false;

					thisGP->d_Kernel = 0.0;
					thisGP->d3_Kernel_Gradient = glm::dvec3(0.0, 0.0, 0.0);
				}
			}
			a_Runtime[1] += omp_get_wtime() - dRuntime_Block;

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
			a_Runtime[2] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// kernel --------------------------------------------------------- kernel
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				double dRadius_Effective = glm::length(1.0*d3_Length_Cell.x);

				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				GP_Mediator_Thread.findNeighborGridPoints(thisMP->d3_Position, i3_Nodes_Kernel, d3_Length_Cell_Kernel, 1);

				for(int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					unsigned int index_GP = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];
					GridPoint *thisAGP = v_GridPoint_Kernel[index_GP];

					// distance of MP from this AGP
					double dDistance = glm::length(thisMP->d3_Position - thisAGP->d3_Position);

					if(dDistance < dRadius_Effective)
					{
						thisAGP->b_Active = true;

						// find kernel value
						double dDistance_Normalized = dDistance / dRadius_Effective;

						double dShapeValue = 1.0 - 3.0*dDistance_Normalized*dDistance_Normalized + 2.0*dDistance_Normalized*dDistance_Normalized*dDistance_Normalized;
//						double dShapeValue = 1.0 - dDistance_Normalized;

						if(nThreads > 1)	omp_set_lock(v_GridPoint_Lock[index_GP]);
						thisAGP->d_Kernel += dShapeValue;
						if(nThreads > 1)	omp_unset_lock(v_GridPoint_Lock[index_GP]);
					}
				}
			}
			a_Runtime[3] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// kernel gradient ------------------------------------------------ kernel gradient
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < v_GridPoint_Kernel.size(); index_GP++)
			{
				GridPoint *thisGP = v_GridPoint_Kernel[index_GP];
				if(thisGP->d_Kernel < 1.0)// skip GPs that don't have a kernel value
					continue;

				if(!(thisGP->b_Active))
					continue;

				glm::ivec3 i3Index = thisGP->i3_Index;

				unsigned int iIndex_x_Befor = GridPoint_Factory::getIndex(i3Index + glm::ivec3(-1,0,0), i3_Nodes_Kernel);
				unsigned int iIndex_x_After = GridPoint_Factory::getIndex(i3Index + glm::ivec3(+1,0,0), i3_Nodes_Kernel);
				unsigned int iIndex_y_Befor = GridPoint_Factory::getIndex(i3Index + glm::ivec3(0,-1,0), i3_Nodes_Kernel);
				unsigned int iIndex_y_After = GridPoint_Factory::getIndex(i3Index + glm::ivec3(0,+1,0), i3_Nodes_Kernel);
				unsigned int iIndex_z_Befor = GridPoint_Factory::getIndex(i3Index + glm::ivec3(0,0,-1), i3_Nodes_Kernel);
				unsigned int iIndex_z_After = GridPoint_Factory::getIndex(i3Index + glm::ivec3(0,0,+1), i3_Nodes_Kernel);

				glm::dvec3 d3KernelGradient = glm::dvec3(0.0,0.0,0.0);

				if(iIndex_x_Befor != -1) d3KernelGradient.x += +0.5/d3_Length_Cell_Kernel.x * (thisGP->d_Kernel - v_GridPoint_Kernel[iIndex_x_Befor]->d_Kernel);
				if(iIndex_x_After != -1) d3KernelGradient.x += -0.5/d3_Length_Cell_Kernel.x * (thisGP->d_Kernel - v_GridPoint_Kernel[iIndex_x_After]->d_Kernel);
				if(iIndex_y_Befor != -1) d3KernelGradient.y += +0.5/d3_Length_Cell_Kernel.y * (thisGP->d_Kernel - v_GridPoint_Kernel[iIndex_y_Befor]->d_Kernel);
				if(iIndex_y_After != -1) d3KernelGradient.y += -0.5/d3_Length_Cell_Kernel.y * (thisGP->d_Kernel - v_GridPoint_Kernel[iIndex_y_After]->d_Kernel);
				if(iIndex_z_Befor != -1) d3KernelGradient.z += +0.5/d3_Length_Cell_Kernel.z * (thisGP->d_Kernel - v_GridPoint_Kernel[iIndex_z_Befor]->d_Kernel);
				if(iIndex_z_After != -1) d3KernelGradient.z += -0.5/d3_Length_Cell_Kernel.z * (thisGP->d_Kernel - v_GridPoint_Kernel[iIndex_z_After]->d_Kernel);

				thisGP->d3_Kernel_Gradient = d3KernelGradient;

				if(glm::length(d3KernelGradient) < 1.0/d3_Length_Cell_Kernel.x)
					thisGP->d3_Kernel_Gradient = glm::dvec3(0.0,0.0,0.0);
			}
			a_Runtime[4] += omp_get_wtime() - dRuntime_Block;

			// GP kernel (from grid values) --------------------------- GP contact kernel and gradient

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// MP kernel gradient --------------------------------------------- MP kernel gradient
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				thisMP->d_Kernel = 0.0;
				thisMP->d3_Kernel_Gradient = glm::dvec3(0.0, 0.0, 0.0);

				GP_Mediator_Thread.findNeighborGridPoints(thisMP->d3_Position, i3_Nodes_Kernel, d3_Length_Cell_Kernel, 0);

				for(int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					unsigned int index_GP = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];
					GridPoint *thisAGP = v_GridPoint_Kernel[index_GP];

					// shape value and shape gradient value
					BasesFunction_Thread.calculateBases(thisMP->d3_Position, thisAGP->d3_Position, d3_Length_Cell_Kernel);
					double dShapeValue = BasesFunction_Thread.d_ShapeValue;

					thisMP->d_Kernel += dShapeValue * thisAGP->d_Kernel;
					thisMP->d3_Kernel_Gradient += dShapeValue * thisAGP->d3_Kernel_Gradient;
				}

				// find the MP nearest to each GP
				for(unsigned int index_AGP = 0; index_AGP < v_MP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = v_MP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];

					if(nThreads > 1)	omp_set_lock(v_GridPoint_Lock[index_GP]);
					if(glm::length(thisMP->d3_Position - thisAGP->d3_Position) < glm::length(allMaterialPoint[thisAGP->i_NearestMP]->d3_Position - thisAGP->d3_Position))
					{
						thisAGP->i_NearestMP = index_MP;
					}
					if(nThreads > 1)	omp_unset_lock(v_GridPoint_Lock[index_GP]);
				}
				// use material point kernel gradients to find grid point kernel gradients
				for(unsigned int index_AGP = 0; index_AGP < v_MP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = v_MP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];

					double dShapeValue = v_MP_AGP[index_MP][index_AGP].dShapeValue;

//					if(thisMP->b_Surface && glm::length(dShapeValue*thisMP->d3_Kernel_Gradient) > glm::length(thisAGP->d3_Kernel_Gradient))
//						thisAGP->d3_Kernel_Gradient = dShapeValue*thisMP->d3_Kernel_Gradient;
					if(nThreads > 1)	omp_set_lock(v_GridPoint_Lock[index_GP]);
					if(thisMP->b_Surface && glm::length(1.0*thisMP->d3_Kernel_Gradient) > glm::length(thisAGP->d3_Kernel_Gradient))
					{
						thisAGP->d3_Kernel_Gradient = 1.0*thisMP->d3_Kernel_Gradient;
					}
					if(nThreads > 1)	omp_unset_lock(v_GridPoint_Lock[index_GP]);
				}
			}
			a_Runtime[5] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// detect contact ------------------------------------------------- detect contacts
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				if(!(thisMP->b_Surface))
					continue;

				for(unsigned int index_AGP = 0; index_AGP < v_MP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[v_MP_AGP[index_MP][index_AGP].index];

					if(!(allMaterialPoint[thisAGP->i_NearestMP]->b_Surface))
						continue;

					if(glm::dot(thisAGP->d3_Kernel_Gradient, thisMP->d3_Kernel_Gradient) >= 0.0)
					{
						thisAGP->b_Contact_Positive = true;
					}
					else
					{
						thisAGP->b_Contact_Negative = true;
					}
				}
			}
			a_Runtime[6] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
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

					if(nThreads > 1)	omp_set_lock(v_GridPoint_Lock[index_GP]);
					{
						// mass
//						thisAGP->d3_Mass += dShapeValue * thisMP->d3_Mass;
						if(glm::dot(thisAGP->d3_Kernel_Gradient, thisMP->d3_Kernel_Gradient) >= 0.0)
						{
							thisAGP->d3_Mass += dShapeValue * thisMP->d3_Mass;
							thisAGP->d3_MassGradient += d3ShapeGradient * thisMP->d3_Mass;
						}
						else
						{
							thisAGP->d3_Mass_Negative += dShapeValue * thisMP->d3_Mass;
							thisAGP->d3_MassGradient_Negative += d3ShapeGradient * thisMP->d3_Mass;
						}
					}
					if(nThreads > 1)	omp_unset_lock(v_GridPoint_Lock[index_GP]);
				}
			}
			a_Runtime[7] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
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

					if(nThreads > 1)	omp_set_lock(v_GridPoint_Lock[index_GP]);
					{
//						// velocity
//						if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
//							thisAGP->d3_Velocity += dShapeValue * (thisMP->d3_Mass * thisMP->d3_Velocity) / thisAGP->d3_Mass;
//
//						// internal forces
//						double dVolume = thisMP->d_Volume;
//						thisAGP->d3_Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6_Stress[0] + d3ShapeGradient.y*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[5]);
//						thisAGP->d3_Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6_Stress[1] + d3ShapeGradient.x*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[4]);
//						thisAGP->d3_Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6_Stress[2] + d3ShapeGradient.x*thisMP->d6_Stress[5] + d3ShapeGradient.y*thisMP->d6_Stress[4]);
//
//						// external forces
//						thisAGP->d3_Force += dShapeValue*thisMP->d3_Force_External;

						if(glm::dot(thisAGP->d3_Kernel_Gradient, thisMP->d3_Kernel_Gradient) >= 0.0)
						{
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
						else
						{
							// velocity
							if(glm::length(thisAGP->d3_Mass_Negative) > d_Mass_Minimum)
								thisAGP->d3_Velocity_Negative += dShapeValue * (thisMP->d3_Mass * thisMP->d3_Velocity) / thisAGP->d3_Mass_Negative;

							// internal forces
							double dVolume = thisMP->d_Volume;
							thisAGP->d3_Force_Negative.x += -dVolume * (d3ShapeGradient.x*thisMP->d6_Stress[0] + d3ShapeGradient.y*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[5]);
							thisAGP->d3_Force_Negative.y += -dVolume * (d3ShapeGradient.y*thisMP->d6_Stress[1] + d3ShapeGradient.x*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[4]);
							thisAGP->d3_Force_Negative.z += -dVolume * (d3ShapeGradient.z*thisMP->d6_Stress[2] + d3ShapeGradient.x*thisMP->d6_Stress[5] + d3ShapeGradient.y*thisMP->d6_Stress[4]);

							// external forces
							thisAGP->d3_Force_Negative += dShapeValue*thisMP->d3_Force_External;
						}
					}
					if(nThreads > 1)	omp_unset_lock(v_GridPoint_Lock[index_GP]);
				}
			}
			a_Runtime[8] += omp_get_wtime() - dRuntime_Block;

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
					// positive contacts
					if(glm::length(thisGP->d3_Velocity) > 1.0e-9)
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

					// negative contacts
					if(glm::length(thisGP->d3_Velocity_Negative) > 1.0e-16)
						thisGP->d3_Force_Negative -= d_DampingCoefficient * glm::length(thisGP->d3_Force_Negative) * glm::normalize(thisGP->d3_Velocity_Negative);

					if(glm::length(thisGP->d3_Mass_Negative) > d_Mass_Minimum)
						thisGP->d3_Velocity_Negative += thisGP->d3_Force_Negative / thisGP->d3_Mass_Negative * dTimeIncrement;

					if(thisGP->b3_Fixed.x == true)
					{
						thisGP->d3_Velocity_Negative.x = 0.0;
						thisGP->d3_Force_Negative.x = 0.0;
					}
					if(thisGP->b3_Fixed.y == true)
					{
						thisGP->d3_Velocity_Negative.y = 0.0;
						thisGP->d3_Force_Negative.y = 0.0;
					}
					if(thisGP->b3_Fixed.z == true)
					{
						thisGP->d3_Velocity_Negative.z = 0.0;
						thisGP->d3_Force_Negative.z = 0.0;
					}

					if(thisGP->b_Contact_Negative == true)
					{// contact stuff
						glm::dvec3 d3Velocity_CM = glm::dvec3(0.0,0.0,0.0);

						d3Velocity_CM += thisGP->d3_Mass * thisGP->d3_Velocity;
						d3Velocity_CM += thisGP->d3_Mass_Negative * thisGP->d3_Velocity_Negative;
						d3Velocity_CM /= (thisGP->d3_Mass + thisGP->d3_Mass_Negative);

						glm::dvec3 d3Force_Contact = glm::dvec3(0.0,0.0,0.0);
						glm::dvec3 d3Force_Contact_Negative = glm::dvec3(0.0,0.0,0.0);

						glm::dvec3 d3Normal = glm::dvec3(0.0,0.0,0.0);
	//					if(glm::length(thisGP->d3_Kernel_ContactGradient) > 1.0)
//							d3Normal = glm::normalize(thisGP->d3_Kernel_Gradient);
						d3Normal = -0.5*(thisGP->d3_MassGradient - thisGP->d3_MassGradient_Negative);

						d3Normal = glm::normalize(d3Normal);

						if(glm::dot(thisGP->d3_Velocity,d3Normal) + glm::dot(thisGP->d3_Velocity_Negative,-d3Normal) < 0.0)
						{// there is penetration
							//std::cout << "penetration" << std::endl;
							d3Force_Contact = thisGP->d3_Mass/dTimeIncrement * glm::dot((d3Velocity_CM - thisGP->d3_Velocity), d3Normal) * d3Normal;
							d3Force_Contact_Negative = thisGP->d3_Mass_Negative/dTimeIncrement * glm::dot((d3Velocity_CM - thisGP->d3_Velocity_Negative), -d3Normal) * -d3Normal;
						}

						if(glm::length(thisGP->d3_Mass) > d_Mass_Minimum)
						{
							thisGP->d3_Force += d3Force_Contact;
							thisGP->d3_Velocity += d3Force_Contact / thisGP->d3_Mass * dTimeIncrement;
						}
						if(glm::length(thisGP->d3_Mass_Negative) > d_Mass_Minimum)
						{
							thisGP->d3_Force_Negative += d3Force_Contact_Negative;
							thisGP->d3_Velocity_Negative += d3Force_Contact_Negative / thisGP->d3_Mass_Negative * dTimeIncrement;
						}
					}

					// total, sina, this is just for debugging, remove for actual contact algorithm
//					thisGP->d3_Mass += thisGP->d3_Mass_Negative;
//					thisGP->d3_Velocity = (thisGP->d3_Mass * thisGP->d3_Velocity + thisGP->d3_Mass_Negative * thisGP->d3_Velocity_Negative) / (thisGP->d3_Mass + thisGP->d3_Mass_Negative);
//					thisGP->d3_Force += thisGP->d3_Force_Negative;
				}
			}
			a_Runtime[9] += omp_get_wtime() - dRuntime_Block;

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

					if(nThreads > 1)	omp_set_lock(v_GridPoint_Lock[index_GP]);
					{
						thisAGP->d3_Velocity = thisMP->d3_Velocity;
						thisAGP->d3_Velocity_Negative = thisMP->d3_Velocity;

						thisAGP->d3_Force_Temp += thisAGP->d3_Force;
						thisAGP->d3_Force_Temp += thisAGP->d3_Force_Negative;

						thisAGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
						thisAGP->d3_Force_Negative = glm::dvec3(0.0, 0.0, 0.0);
					}
					if(nThreads > 1)	omp_unset_lock(v_GridPoint_Lock[index_GP]);
				}
			}
			a_Runtime[10] += omp_get_wtime() - dRuntime_Block;

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

//					// velocity
//					if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
//						thisMP->d3_Velocity += dShapeValue * (thisAGP->d3_Force/thisAGP->d3_Mass) * dTimeIncrement;
//
//					// position
//					thisMP->d3_Position += dShapeValue * (thisAGP->d3_Velocity) * dTimeIncrement;
//
//					// velocity gradient, to be used to calculate strains
//					d33VelocityGradient += glm::outerProduct(thisAGP->d3_Velocity, d3ShapeGradient);// this glm function does the pre-transposition that we want

					if(glm::dot(thisAGP->d3_Kernel_Gradient, thisMP->d3_Kernel_Gradient) >= 0.0)
					{
						// velocity
						if(glm::length(thisAGP->d3_Mass) > d_Mass_Minimum)
							thisMP->d3_Velocity += dShapeValue * (thisAGP->d3_Force/thisAGP->d3_Mass) * dTimeIncrement;

						// position
						thisMP->d3_Position += dShapeValue * (thisAGP->d3_Velocity) * dTimeIncrement;

						// velocity gradient, to be used to calculate strains
						d33VelocityGradient += glm::outerProduct(thisAGP->d3_Velocity, d3ShapeGradient);// this glm function does the pre-transposition that we want
					}
					else
					{
						// velocity
						if(glm::length(thisAGP->d3_Mass_Negative) > d_Mass_Minimum)
							thisMP->d3_Velocity += dShapeValue * (thisAGP->d3_Force_Negative/thisAGP->d3_Mass_Negative) * dTimeIncrement;

						// position
						thisMP->d3_Position += dShapeValue * (thisAGP->d3_Velocity_Negative) * dTimeIncrement;

						// velocity gradient, to be used to calculate strains
						d33VelocityGradient += glm::outerProduct(thisAGP->d3_Velocity_Negative, d3ShapeGradient);// this glm function does the pre-transposition that we want
					}
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
			a_Runtime[11] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement_Control.size(); index_MP++)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement_Control[index_MP];

				thisMP->d3_Position += thisMP->d3_Velocity * dTimeIncrement;
			}
			a_Runtime[12] += omp_get_wtime() - dRuntime_Block;
		}

//		d_Runtime_Total += double(double(clock()-clockCurrent_Total)/CLOCKS_PER_SEC);
		d_Runtime_Total += omp_get_wtime() - dRuntime_MP;
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
		dTimeConsole_Last = dTime;
		this->reportConsole();
		std::cout << "Analysis complete ..." << std::endl;
		return(1);
	}
}
// ----------------------------------------------------------------------------
