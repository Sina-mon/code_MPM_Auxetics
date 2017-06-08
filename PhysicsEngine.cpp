#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
PhysicsEngine::PhysicsEngine()
{
};
// ----------------------------------------------------------------------------
PhysicsEngine::~PhysicsEngine()
{
	//delete all objects created by Factory classes
//	for(unsigned int index = 0; index < allMaterialPoint.size(); index++)
//        delete allMaterialPoint[index];

	for(unsigned int index = 0; index < v_MaterialPoint_Kinetics.size(); index++)
        delete v_MaterialPoint_Kinetics[index];

	for(unsigned int index = 0; index < v_MaterialPoint_Material.size(); index++)
        delete v_MaterialPoint_Material[index];

	for(unsigned int index = 0; index < allGridPoint.size(); index++)
        delete allGridPoint[index];

	for(int iThread = 0; iThread < _MAX_N_THREADS; iThread++)
	{
		for(unsigned int index = 0; index < allGridPoint_Thread[iThread].size(); index++)
		{
			delete allGridPoint_Thread[iThread][index];
		}
	}

	// destroy GridPoint locks
	for(int index = 0; index < v_GridPoint_Lock.size(); index++)
	{
		omp_destroy_lock(v_GridPoint_Lock[index]);
		delete v_GridPoint_Lock[index];
	}
}
// ----------------------------------------------------------------------------
void PhysicsEngine::reportConsole(std::string sDescription)
{
	double dMass = 0.0;
	double dMass_Negative = 0.0;
	double dMomentum_x = 0.0;
	double dMomentum_y = 0.0;
	double dMomentum_z = 0.0;

	for(unsigned int index = 0; index < allGridPoint.size(); index++)
	{// calculate debug values
		dMass += allGridPoint[index]->d_Mass;
//		dMass_Negative += allGridPoint[index]->d3_Mass_Negative.x;
	}

	for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Momentum.size(); index_MP++)
	{// calculate debug values
		dMomentum_x += v_MarkedMaterialPoints_Momentum[index_MP]->d3_Velocity.x * v_MarkedMaterialPoints_Momentum[index_MP]->d_Mass;
		dMomentum_y += v_MarkedMaterialPoints_Momentum[index_MP]->d3_Velocity.y * v_MarkedMaterialPoints_Momentum[index_MP]->d_Mass;
	}

	glm::dvec3 d3Stress = glm::dvec3(0.0,0.0,0.0);
	glm::dvec3 d3Strain = glm::dvec3(0.0,0.0,0.0);
	for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Stress_Monitor.size(); index_MP++)
	{// calculate debug values
		MaterialPoint_Material *thisMP = v_MarkedMaterialPoints_Stress_Monitor[index_MP];

		d3Stress += glm::dvec3(thisMP->d6_Stress[0], thisMP->d6_Stress[1], thisMP->d6_Stress[2]);
		d3Strain += glm::dvec3(thisMP->d6_Strain[0], thisMP->d6_Strain[1], thisMP->d6_Strain[2]);
	}
	d3Stress /= v_MarkedMaterialPoints_Stress_Monitor.size();
	d3Strain /= v_MarkedMaterialPoints_Stress_Monitor.size();

	glm::dvec3 d3Force = glm::dvec3(0.0,0.0,0.0);
	for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
	{// calculate debug values
		GridPoint *thisGP = allGridPoint[index_GP];

		d3Force += thisGP->d3_Force_Temp;
	}

	double dKineticEnergy = 0.0;
	for(int index = 0; index < v_MaterialPoint_Kinetics.size(); index++)
	{
		MaterialPoint_Kinetics *thisMP_Kinetics = v_MaterialPoint_Kinetics[index];

		dKineticEnergy += thisMP_Kinetics->d_Mass * glm::pow(glm::length(thisMP_Kinetics->d3_Velocity),2.0);
	}

	std::string strConsole = "";
	strConsole += sDescription;
	strConsole += "\ttime: " + Script(d_Time,6);
	strConsole += "\tRuntime: " + Script(d_Runtime_Total,3);
//	strConsole += "\tMass_P: " + Script(dMass,6);
//	strConsole += "\tMass_N: " + Script(dMass_Negative,6);
//	strConsole += "\tmomentum_x: " + Script(dMomentum_x,3) + "\t momentum_y: " + Script(dMomentum_y,3) + "\t momentum_z: " + Script(dMomentum_z,3);
//	strConsole += "\tmomentum_x: " + Script(dMomentum_x,3) + "\t momentum_y: " + Script(dMomentum_y,6);
	if(v_MarkedMaterialPoints_Displacement_Monitor.size() > 0)
        strConsole += "\tPosition_y: " + Script(v_MarkedMaterialPoints_Displacement_Monitor[0]->d3_Position.y,6);
	if(v_MarkedMaterialPoints_Stress_Monitor.size() > 0)
	{
		strConsole += "\tStrain_y: " + Script(d3Strain.y,6);
		strConsole += "\tStress_y: " + Script(d3Stress.y,6);
	}
	strConsole += "\tForce_y: " + Script(d3Force.y,6);
	strConsole += "\tKinetic Energy: " + Script(dKineticEnergy,6);
	strConsole += "\n";

	if(false)
	{// runtimes
		double dRuntime_Total = 0.0;
		for(int index = 0; index < a_Runtime.size(); index++)
		{
			strConsole += "\tRuntime_Block[" + Script(index) + "]: " + Script(a_Runtime[index],3) + "(" + Script(a_Runtime[index]/d_Runtime_Total,2,std::ios::fixed) + "\%)";
			strConsole += "\n";
			dRuntime_Total += a_Runtime[index];
		}
		strConsole += "\tRuntime_Total: " + Script(dRuntime_Total,3);
		strConsole += "\n";
	}

	std::cout << strConsole;

	std::string strFileName = _STR_LOGFILE;
	std::ofstream OutputFile(strFileName.c_str(), std::ios_base::app);

	OutputFile << strConsole;

	OutputFile.close();
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
