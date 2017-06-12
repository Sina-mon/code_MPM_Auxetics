#include "ConstitutiveRelation.h"

// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateState_J2(double d6State[6])
{

	double dStress_Hydrostatic = 1.0/3.0 * (d6State[0] + d6State[1] + d6State[2]);

	double d6Stress_Deviatoric[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	for(int index = 0; index < 6; index++)
		d6Stress_Deviatoric[index] = d6State[index];

	d6Stress_Deviatoric[0] -= dStress_Hydrostatic;
	d6Stress_Deviatoric[1] -= dStress_Hydrostatic;
	d6Stress_Deviatoric[2] -= dStress_Hydrostatic;

	double dJ2 = 0.0;
	dJ2 += 1.0/6.0 * (d6State[0] - d6State[1])*(d6State[0] - d6State[1]);
	dJ2 += 1.0/6.0 * (d6State[1] - d6State[2])*(d6State[1] - d6State[2]);
	dJ2 += 1.0/6.0 * (d6State[2] - d6State[0])*(d6State[2] - d6State[0]);
	dJ2 += d6State[3]*d6State[3];
	dJ2 += d6State[4]*d6State[4];
	dJ2 += d6State[5]*d6State[5];

	this->d_J2 = dJ2;
}
// ----------------------------------------------------------------------------
