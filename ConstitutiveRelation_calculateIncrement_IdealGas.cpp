#include "ConstitutiveRelation.h"

// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateIncrement_IdealGass(double dHeatCapacityRatio, double dSpecificHeat, double dShearViscosity, double dBulkViscosity, double d6Stress_Current[6], double dDensity, double dTemperature, double d6Strain_Rate[6])
{
	for(int index = 0; index < 6; index++)
	{
		this->d6StressIncrement[index] = 0.0;
		this->d6PlasticStrainIncrement[index] = 0.0;
	}

	double d6Stress_New[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	double dPressure = -dBulkViscosity*(d6Strain_Rate[0] + d6Strain_Rate[1] + d6Strain_Rate[2]) - (dHeatCapacityRatio - 1.0) * dDensity * dSpecificHeat * dTemperature;
	d6Stress_New[0] = 2.0 * dShearViscosity * d6Strain_Rate[0] + dPressure;
	d6Stress_New[1] = 2.0 * dShearViscosity * d6Strain_Rate[1] + dPressure;
	d6Stress_New[2] = 2.0 * dShearViscosity * d6Strain_Rate[2] + dPressure;
	d6Stress_New[3] = 2.0 * dShearViscosity * d6Strain_Rate[3];// sina, not sure about these shear formulations
	d6Stress_New[4] = 2.0 * dShearViscosity * d6Strain_Rate[4];
	d6Stress_New[5] = 2.0 * dShearViscosity * d6Strain_Rate[5];

	for(int index = 0; index < 6; index++)
		this->d6StressIncrement[index] = d6Stress_New[index] - d6Stress_Current[index];
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
