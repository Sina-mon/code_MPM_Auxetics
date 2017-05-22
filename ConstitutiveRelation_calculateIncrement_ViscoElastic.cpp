#include "ConstitutiveRelation.h"

// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateIncrement_ViscoElastic_6D(double dE, double dNu, double dEta, double d6Stress_Current[6], double d6Strain_New[6], double d6Strain_Rate_New[6])
{
	double d6Stress_Elastic[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double d6Stress_Viscous[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double d6Stress_New[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	double dConstant = dE/(1.0 + dNu)/(1.0 - 2.0*dNu);

	d6Stress_Elastic[0] = dConstant * ((1.0-dNu)*d6Strain_New[0] + dNu*d6Strain_New[1] + dNu*d6Strain_New[2]);
	d6Stress_Elastic[1] = dConstant * ((1.0-dNu)*d6Strain_New[1] + dNu*d6Strain_New[0] + dNu*d6Strain_New[2]);
	d6Stress_Elastic[2] = dConstant * ((1.0-dNu)*d6Strain_New[2] + dNu*d6Strain_New[0] + dNu*d6Strain_New[1]);
	d6Stress_Elastic[3] = dConstant * ((0.5-dNu)*d6Strain_New[3]);
	d6Stress_Elastic[4] = dConstant * ((0.5-dNu)*d6Strain_New[4]);
	d6Stress_Elastic[5] = dConstant * ((0.5-dNu)*d6Strain_New[5]);

	for(int index = 0; index < 6; index++)
		d6Stress_Viscous[index] = dEta * d6Strain_Rate_New[index];

	for(int index = 0; index < 6; index++)
		this->d6StressIncrement[index] = (d6Stress_Elastic[index] + d6Stress_Viscous[index]) - d6Stress_Current[index];
	for(int index = 0; index < 6; index++)
		this->d6PlasticStrainIncrement[index] = 0.0;
}
// ----------------------------------------------------------------------------
