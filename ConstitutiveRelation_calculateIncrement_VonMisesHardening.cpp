#include "ConstitutiveRelation.h"

// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateIncrement_VonMisesHardening_6D(double dE, double dNu, double dYield, double dBackstress_Isotropic, double dHardening_Isotropic_C0, double dHardening_Isotropic_C1, double d6StressCurrent[6], double d6StrainIncrement[6])
{
	double d6StressIncrement_Trial[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double d6Stress_Trial[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	double dConstant = dE/(1.0 + dNu)/(1.0 - 2.0*dNu);

	d6StressIncrement_Trial[0] = dConstant * ((1.0-dNu)*d6StrainIncrement[0] + dNu*d6StrainIncrement[1] + dNu*d6StrainIncrement[2]);
	d6StressIncrement_Trial[1] = dConstant * ((1.0-dNu)*d6StrainIncrement[1] + dNu*d6StrainIncrement[0] + dNu*d6StrainIncrement[2]);
	d6StressIncrement_Trial[2] = dConstant * ((1.0-dNu)*d6StrainIncrement[2] + dNu*d6StrainIncrement[0] + dNu*d6StrainIncrement[1]);
	d6StressIncrement_Trial[3] = dConstant * ((0.5-dNu)*d6StrainIncrement[3]);
	d6StressIncrement_Trial[4] = dConstant * ((0.5-dNu)*d6StrainIncrement[4]);
	d6StressIncrement_Trial[5] = dConstant * ((0.5-dNu)*d6StrainIncrement[5]);

	for(int index = 0; index < 6; index++)
		d6Stress_Trial[index] = d6StressCurrent[index] + d6StressIncrement_Trial[index];

	double dStress_Hydrostatic_Trial = 1.0/3.0 * (d6Stress_Trial[0] + d6Stress_Trial[1] + d6Stress_Trial[2]);

	double d6Stress_Deviatoric_Trial[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	for(int index = 0; index < 6; index++)
		d6Stress_Deviatoric_Trial[index] = d6Stress_Trial[index];

	d6Stress_Deviatoric_Trial[0] -= dStress_Hydrostatic_Trial;
	d6Stress_Deviatoric_Trial[1] -= dStress_Hydrostatic_Trial;
	d6Stress_Deviatoric_Trial[2] -= dStress_Hydrostatic_Trial;

	double dJ2 = 0.0;
	dJ2 += 1.0/6.0 * (d6Stress_Trial[0] - d6Stress_Trial[1])*(d6Stress_Trial[0] - d6Stress_Trial[1]);
	dJ2 += 1.0/6.0 * (d6Stress_Trial[1] - d6Stress_Trial[2])*(d6Stress_Trial[1] - d6Stress_Trial[2]);
	dJ2 += 1.0/6.0 * (d6Stress_Trial[2] - d6Stress_Trial[0])*(d6Stress_Trial[2] - d6Stress_Trial[0]);
	dJ2 += d6Stress_Trial[3]*d6Stress_Trial[3];
	dJ2 += d6Stress_Trial[4]*d6Stress_Trial[4];
	dJ2 += d6Stress_Trial[5]*d6Stress_Trial[5];

	double dK0 = dYield / sqrt(3.0);
	double dR = dBackstress_Isotropic / sqrt(2.0);

	if(dJ2 <= (dK0+dR)*(dK0+dR))
	{
		for(int index = 0; index < 6; index++)
			d6StressIncrement[index] = d6StressIncrement_Trial[index];
		for(int index = 0; index < 6; index++)
			d6PlasticStrainIncrement[index] = 0.0;

		dBackstress_IsotropicIncrement = 0.0;
	}
	else
	{
		for(int index = 0; index < 6; index++)
			d6Stress_Deviatoric_Trial[index] *= (dK0+dR) / sqrt(dJ2);

		for(int index = 0; index < 6; index++)
			d6StressIncrement[index] = d6Stress_Deviatoric_Trial[index] - d6StressCurrent[index];

		d6StressIncrement[0] += dStress_Hydrostatic_Trial;
		d6StressIncrement[1] += dStress_Hydrostatic_Trial;
		d6StressIncrement[2] += dStress_Hydrostatic_Trial;

		double dG = 0.5*dE/(1.0+dNu);
		double d6StrainIncrement_Elastic[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};// sina, this strain increment might not be correct when assuming plain strain conditions
		d6StrainIncrement_Elastic[0] = d6StressIncrement[0]/dE - d6StressIncrement[1]*dNu/dE - d6StressIncrement[2]*dNu/dE;
		d6StrainIncrement_Elastic[1] = d6StressIncrement[1]/dE - d6StressIncrement[0]*dNu/dE - d6StressIncrement[2]*dNu/dE;
		d6StrainIncrement_Elastic[2] = d6StressIncrement[2]/dE - d6StressIncrement[0]*dNu/dE - d6StressIncrement[1]*dNu/dE;
		d6StrainIncrement_Elastic[3] = d6StressIncrement[3]/dG;
		d6StrainIncrement_Elastic[4] = d6StressIncrement[4]/dG;
		d6StrainIncrement_Elastic[5] = d6StressIncrement[5]/dG;

		for(int index = 0; index < 6; index++)
			d6PlasticStrainIncrement[index] = d6StrainIncrement[index] - d6StrainIncrement_Elastic[index];

		double dEffectivePlasticStrainIncrement = 0.0;
		for(int index = 0; index < 6; index++)
			dEffectivePlasticStrainIncrement += d6PlasticStrainIncrement[index] * d6PlasticStrainIncrement[index];
		dEffectivePlasticStrainIncrement *= 2.0/3.0;
		dEffectivePlasticStrainIncrement = glm::sqrt(dEffectivePlasticStrainIncrement);

		if((sqrt(2.0/3.0) * dHardening_Isotropic_C1 - dBackstress_Isotropic) > 0.0)
			dBackstress_IsotropicIncrement = dHardening_Isotropic_C0 * (sqrt(2.0/3.0) * dHardening_Isotropic_C1 - dBackstress_Isotropic) * dEffectivePlasticStrainIncrement;
		else
			dBackstress_IsotropicIncrement  = 0.0;
	}
}
// ----------------------------------------------------------------------------
