/*
make a convenience function to set all outpur variables (stressincrement, plastic strain increment, etc to zero)
*/

#ifndef CONSTITUTIVERELATION_H
#define CONSTITUTIVERELATION_H

#include <iostream>

#include "Definitions.h"
////#include ".\include\glm\glm.hpp" // windows
//#include "./include/glm/glm.hpp" // linux

class ConstitutiveRelation
{
	public:
		ConstitutiveRelation();
		virtual ~ConstitutiveRelation();

		double d6StressIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double d6PlasticStrainIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double dBackstress_IsotropicIncrement = 0.0;
		double d_J2 = 0.0;

		void calculateIncrement_Elastic(double dE, double dNu, double d6StrainIncrement[6]);
		void calculateIncrement_NeoHookean_6D(double dE, double dNu, double d6Stress_Current[6], glm::dmat3 d33DeformationGradient);
		void calculateIncrement_ViscoElastic_6D(double dE, double dNu, double dEta, double d6Stress_Current[6], double d6Strain_New[6], double d6Strain_Rate_New[6]);
		void calculateIncrement_PerfectlyPlastic_6D(double dE, double dNu, double dYield, double d6StressCurrent[6], double d6StrainIncrement[6]);
		void calculateIncrement_VonMisesHardening_6D(double dE, double dNu, double dYield, double dBackstress_Isotropic, double dHardening_Isotropic_C0, double dHardening_Isotropic_C1, double d6StressCurrent[6], double d6StrainIncrement[6]);
		void calculateIncrement_IdealGass(double dHeatCapacityRatio, double dSpecificHeat, double dShearViscosity, double dBulkViscosity, double d6Stress_Current[6], double dDensity, double dTemperature, double d6Strain_Rate[6]);

		void calculateState_J2(double d6Stress[6]);
	protected:

	private:
};

#endif // CONSTITUTIVERELATION_H
