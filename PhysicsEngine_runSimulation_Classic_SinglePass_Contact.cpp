#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_Classic_SinglePass_Contact(double dTimeIncrement_Total)
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

//	double dDebug_ContactCutoff = -10.0e24;// to go back to single field contact
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

				thisGP->b_Contact_Positive = false;
				thisGP->b_Contact_Negative = false;

				thisGP->d3_Mass = {0.0, 0.0, 0.0};
				thisGP->d3_MassGradient = {0.0, 0.0, 0.0};
				thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);

				thisGP->d3_Mass_Negative = {0.0, 0.0, 0.0};
				thisGP->d3_MassGradient_Negative = {0.0, 0.0, 0.0};
				thisGP->d3_Velocity_Negative = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force_Negative = glm::dvec3(0.0, 0.0, 0.0);

				thisGP->i_NearestMP = 0;

				thisGP->d_Kernel_Contact = 0.0;
				thisGP->d3_Kernel_ContactGradient = glm::dvec3(0.0, 0.0, 0.0);
			}

			#pragma omp barrier
			// reset grid kernel points --------------------------------------- reset grid kernel points
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < v_GridPoint_ContactKernel.size(); index_GP++)
			{
				GridPoint *thisGP = v_GridPoint_ContactKernel[index_GP];

				thisGP->b_Contact_Positive = false;
				thisGP->b_Contact_Negative = false;

				thisGP->d3_Mass = {0.0, 0.0, 0.0};
				thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);

				thisGP->d3_Mass_Negative = {0.0, 0.0, 0.0};
				thisGP->d3_Velocity_Negative = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force_Negative = glm::dvec3(0.0, 0.0, 0.0);

				thisGP->d_Kernel_Contact = 0.0;
				thisGP->d3_Kernel_ContactGradient = glm::dvec3(0.0, 0.0, 0.0);
			}

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

					// adjacent grid points
					vMP_AGP[index_MP][index_AGP].index = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];

					// shape value and shape gradient value
					BasesFunction_Thread.calculateBases(thisMP->d3_Position, thisAGP->d3_Position, d3_Length_Cell);

					vMP_AGP[index_MP][index_AGP].dShapeValue = BasesFunction_Thread.d_ShapeValue;
					vMP_AGP[index_MP][index_AGP].d3ShapeGradient = BasesFunction_Thread.d3_ShapeGradient;
				}
			}

			#pragma omp barrier
			// grid contact kernel -------------------------------------------- GP contact kernel
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				double dRadius_Effective = glm::length(1.0*d3_Length_Cell.x);

				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				GP_Mediator_Thread.findNeighborGridPoints(thisMP->d3_Position, i3_Nodes_ContactKernel, d3_Length_Cell_ContactKernel, 2);

				for(int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					unsigned int index_GP = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];
					GridPoint *thisAGP = v_GridPoint_ContactKernel[index_GP];

					// distance of MP from this AGP
					double dDistance = glm::length(thisMP->d3_Position - thisAGP->d3_Position);
					// find kernel value
					double dDistance_Normalized = dDistance / dRadius_Effective;

					if(dDistance_Normalized <= 1.0)
					{
//						double dShapeValue = 1.0 - 3.0*glm::pow(dDistance_Normalized,2.0) + 2.0*glm::pow(dDistance_Normalized, 3.0);
//						thisAGP->d_Kernel_Contact += dShapeValue;

						double dShapeValue = 1.0 - dDistance_Normalized;
						thisAGP->d_Kernel_Contact += dShapeValue;
					}
				}
			}

			#pragma omp barrier
			// grid contact kernel gradient ----------------------------------- GP contact kernel gradient
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < v_GridPoint_ContactKernel.size(); index_GP++)
			{
				GridPoint *thisGP = v_GridPoint_ContactKernel[index_GP];
				if(thisGP->d_Kernel_Contact < 1.0)// skip GPs that don't have a kernel value
					continue;

				glm::ivec3 i3Index = thisGP->i3_Index;

				unsigned int iIndex_x_Befor = GridPoint_Factory::getIndex(i3Index + glm::ivec3(-1,0,0), i3_Nodes_ContactKernel);
				unsigned int iIndex_x_After = GridPoint_Factory::getIndex(i3Index + glm::ivec3(+1,0,0), i3_Nodes_ContactKernel);
				unsigned int iIndex_y_Befor = GridPoint_Factory::getIndex(i3Index + glm::ivec3(0,-1,0), i3_Nodes_ContactKernel);
				unsigned int iIndex_y_After = GridPoint_Factory::getIndex(i3Index + glm::ivec3(0,+1,0), i3_Nodes_ContactKernel);
				unsigned int iIndex_z_Befor = GridPoint_Factory::getIndex(i3Index + glm::ivec3(0,0,-1), i3_Nodes_ContactKernel);
				unsigned int iIndex_z_After = GridPoint_Factory::getIndex(i3Index + glm::ivec3(0,0,+1), i3_Nodes_ContactKernel);

				glm::dvec3 d3KernelGradient = glm::dvec3(0.0,0.0,0.0);

				if(iIndex_x_Befor != -1) d3KernelGradient.x += +0.5/d3_Length_Cell_ContactKernel.x * (thisGP->d_Kernel_Contact - v_GridPoint_ContactKernel[iIndex_x_Befor]->d_Kernel_Contact);
				if(iIndex_x_After != -1) d3KernelGradient.x += -0.5/d3_Length_Cell_ContactKernel.x * (thisGP->d_Kernel_Contact - v_GridPoint_ContactKernel[iIndex_x_After]->d_Kernel_Contact);
				if(iIndex_y_Befor != -1) d3KernelGradient.y += +0.5/d3_Length_Cell_ContactKernel.y * (thisGP->d_Kernel_Contact - v_GridPoint_ContactKernel[iIndex_y_Befor]->d_Kernel_Contact);
				if(iIndex_y_After != -1) d3KernelGradient.y += -0.5/d3_Length_Cell_ContactKernel.y * (thisGP->d_Kernel_Contact - v_GridPoint_ContactKernel[iIndex_y_After]->d_Kernel_Contact);
				if(iIndex_z_Befor != -1) d3KernelGradient.z += +0.5/d3_Length_Cell_ContactKernel.z * (thisGP->d_Kernel_Contact - v_GridPoint_ContactKernel[iIndex_z_Befor]->d_Kernel_Contact);
				if(iIndex_z_After != -1) d3KernelGradient.z += -0.5/d3_Length_Cell_ContactKernel.z * (thisGP->d_Kernel_Contact - v_GridPoint_ContactKernel[iIndex_z_After]->d_Kernel_Contact);

				thisGP->d3_Kernel_ContactGradient = d3KernelGradient;

				if(glm::length(d3KernelGradient) < 1.0/d3_Length_Cell_ContactKernel.x)
					thisGP->d3_Kernel_ContactGradient = glm::dvec3(0.0,0.0,0.0);
			}

			#pragma omp barrier
			// GP contact kernel (from grid values) --------------------------- GP contact kernel and gradient
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				thisGP->d_Kernel_Contact = 0.0;
				thisGP->d3_Kernel_ContactGradient = glm::dvec3(0.0, 0.0, 0.0);

				GP_Mediator_Thread.findNeighborGridPoints(thisGP->d3_Position, i3_Nodes_ContactKernel, d3_Length_Cell_ContactKernel, 0);

				for(int index_AGP_Kernel = 0; index_AGP_Kernel < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP_Kernel++)
				{
					unsigned int index_GP_Kernel = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP_Kernel];
					GridPoint *thisAGP_Kernel = v_GridPoint_ContactKernel[index_GP_Kernel];

					// shape value and shape gradient value
					BasesFunction_Thread.calculateBases(thisGP->d3_Position, thisAGP_Kernel->d3_Position, d3_Length_Cell_ContactKernel);
					double dShapeValue = BasesFunction_Thread.d_ShapeValue;

					thisGP->d_Kernel_Contact += dShapeValue * thisAGP_Kernel->d_Kernel_Contact;
//					thisGP->d3_Kernel_ContactGradient += dShapeValue * thisAGP_Kernel->d3_Kernel_ContactGradient;
				}

//				if(glm::length(thisGP->d3_Kernel_ContactGradient) < 1.0/d3_Length_Cell.x)
//					thisGP->d3_Kernel_ContactGradient = glm::dvec3(0.0,0.0,0.0);
			}

			#pragma omp barrier
			// material point contact kernel (from grid values) --------------- MP contact kernel and gradient
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				thisMP->d_Kernel_Contact = 0.0;
				thisMP->d3_Kernel_ContactGradient = glm::dvec3(0.0, 0.0, 0.0);

				GP_Mediator_Thread.findNeighborGridPoints(thisMP->d3_Position, i3_Nodes_ContactKernel, d3_Length_Cell_ContactKernel, 0);

				for(int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					unsigned int index_GP = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];
					GridPoint *thisAGP = v_GridPoint_ContactKernel[index_GP];

					// shape value and shape gradient value
					BasesFunction_Thread.calculateBases(thisMP->d3_Position, thisAGP->d3_Position, d3_Length_Cell_ContactKernel);
					double dShapeValue = BasesFunction_Thread.d_ShapeValue;

					thisMP->d_Kernel_Contact += dShapeValue * thisAGP->d_Kernel_Contact;
					thisMP->d3_Kernel_ContactGradient += dShapeValue * thisAGP->d3_Kernel_ContactGradient;
				}

				if(5.0/d3_Length_Cell.x < glm::length(thisMP->d3_Kernel_ContactGradient))
					thisMP->i_ID = 2;

				// find the MP nearest to each GP
				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = vMP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];

					if(glm::length(thisMP->d3_Position - thisAGP->d3_Position) < glm::length(allMaterialPoint[thisAGP->i_NearestMP]->d3_Position - thisAGP->d3_Position))
					{
						thisAGP->i_NearestMP = index_MP;

						if(thisMP->i_ID == 2)
						{
//						std::cout << "here" << std::endl;
//							thisAGP->d3_Kernel_ContactGradient = thisMP->d3_Kernel_ContactGradient;
						}
					}
				}

				// use material point kernel gradients to find grid point kernel gradients
				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = vMP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;

					if(thisMP->i_ID == 2 && glm::length(dShapeValue*thisMP->d3_Kernel_ContactGradient) > glm::length(thisAGP->d3_Kernel_ContactGradient))
						thisAGP->d3_Kernel_ContactGradient = dShapeValue*thisMP->d3_Kernel_ContactGradient;
				}
			}

			#pragma omp barrier
			// detect contact ------------------------------------------------- detect contacts
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[vMP_AGP[index_MP][index_AGP].index];

					if(glm::dot(thisAGP->d3_Kernel_ContactGradient, thisMP->d3_Kernel_ContactGradient) >= dDebug_ContactCutoff)
					{
						thisAGP->b_Contact_Positive = true;
					}
//					if(glm::dot(thisAGP->d3_Kernel_ContactGradient, thisMP->d3_Kernel_ContactGradient) < -dDebug_ContactCutoff)
					else
					{
						thisAGP->b_Contact_Negative = true;
					}

//					if(thisAGP->b_Contact_Negative == true)
//						std::cout << "Negative contact detected" << std::endl;
				}
			}

			#pragma omp barrier
			// material point to grid, only mass ------------------------------ MP to GP (mass)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = vMP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

					if(nThreads == 1)
					{ // avoid atomic operations if there is only 1 thread running
						// mass
//						thisAGP->d3_Mass += dShapeValue * thisMP->d3_Mass;

						if(glm::dot(thisAGP->d3_Kernel_ContactGradient, thisMP->d3_Kernel_ContactGradient) >= dDebug_ContactCutoff)
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
				}
			}

			#pragma omp barrier
			// material point to grid, velocity and force---------------------- MP to GP (velocity and force)
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint *thisMP = allMaterialPoint[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					unsigned int index_GP = vMP_AGP[index_MP][index_AGP].index;
					GridPoint *thisAGP = allGridPoint[index_GP];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

					if(nThreads == 1)
					{ // avoid atomic operations if there is only 1 thread running

						if(glm::dot(thisAGP->d3_Kernel_ContactGradient, thisMP->d3_Kernel_ContactGradient) >= dDebug_ContactCutoff)
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
					}
				}
			}

			#pragma omp barrier
			// sina, haven't applied contact conditions on this yet
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement.size(); index_MP++)
			{
				MaterialPoint *thisMP = v_MarkedMaterialPoints_Displacement[index_MP];

				GP_Mediator_Thread.findNeighborGridPoints(thisMP->d3_Position, i3_Nodes, d3_Length_Cell, 0);

				for(int index_AGP = 0; index_AGP < GP_Mediator_Thread.v_adjacentGridPoints.size(); index_AGP++)
				{
					unsigned int index_GP = GP_Mediator_Thread.v_adjacentGridPoints[index_AGP];
					GridPoint *thisAGP = v_GridPoint_ContactKernel[index_GP];

					if(glm::dot(thisAGP->d3_Kernel_ContactGradient, thisMP->d3_Kernel_ContactGradient) >= dDebug_ContactCutoff)
					{
						thisAGP->d3_Velocity = thisMP->d3_Velocity;
						thisAGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
					}
					else
					{
						thisAGP->d3_Velocity_Negative = thisMP->d3_Velocity;
						thisAGP->d3_Force_Negative = glm::dvec3(0.0, 0.0, 0.0);
					}
				}

				thisMP->d3_Position += thisMP->d3_Velocity * dTimeIncrement;
			}

			#pragma omp barrier
			// update grid momentum and apply boundary conditions ------------- update GP momentum
			// and damping ----------------------------------------------------
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				// positive contacts
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

				// contact stuff
				if(thisGP->b_Contact_Negative == true)
				{
					glm::dvec3 d3Velocity_CM = glm::dvec3(0.0,0.0,0.0);

					d3Velocity_CM += thisGP->d3_Mass * thisGP->d3_Velocity;
					d3Velocity_CM += thisGP->d3_Mass_Negative * thisGP->d3_Velocity_Negative;
					d3Velocity_CM /= (thisGP->d3_Mass + thisGP->d3_Mass_Negative);

					glm::dvec3 d3Force_Contact = glm::dvec3(0.0,0.0,0.0);
					glm::dvec3 d3Force_Contact_Negative = glm::dvec3(0.0,0.0,0.0);

					glm::dvec3 d3Normal = glm::dvec3(0.0,0.0,0.0);
//					if(glm::length(thisGP->d3_Kernel_ContactGradient) > 1.0)
//						d3Normal = glm::normalize(thisGP->d3_Kernel_ContactGradient);
					d3Normal = -0.5*(thisGP->d3_MassGradient - thisGP->d3_MassGradient_Negative);
//					if(glm::length(d3Normal) > 1.0)
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
					// contact stuff
				}

				// total, sina, this is just for debugging, remove for actual contact algorithm
//				thisGP->d3_Mass += thisGP->d3_Mass_Negative;
//				thisGP->d3_Velocity += thisGP->d3_Velocity_Negative;
//				thisGP->d3_Force += thisGP->d3_Force_Negative;
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

				for(unsigned int index_AGP = 0; index_AGP < vMP_AGP[index_MP].size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[vMP_AGP[index_MP][index_AGP].index];

					double dShapeValue = vMP_AGP[index_MP][index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = vMP_AGP[index_MP][index_AGP].d3ShapeGradient;

					if(glm::dot(thisAGP->d3_Kernel_ContactGradient, thisMP->d3_Kernel_ContactGradient) >= dDebug_ContactCutoff)
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
				// hardening
				double db = thisMP->d_Hardening_Isotropic_C0;
				double dQ = thisMP->d_Hardening_Isotropic_C1;

				double d6StressIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				double d6PlasticStrainIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

				ConstitutiveRelation vonMises_Thread;

				if(thisMP->i_MaterialType == _ELASTIC)
					vonMises_Thread.calculateIncrement_Elastic(dE, dNu, d6StrainIncrement);
				else if(thisMP->i_MaterialType == _VISCOELASTIC)
					vonMises_Thread.calculateIncrement_ViscoElastic_6D(dE, dNu, thisMP->d_Viscosity, thisMP->d6_Stress, thisMP->d6_Strain, d6StrainRate);
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
