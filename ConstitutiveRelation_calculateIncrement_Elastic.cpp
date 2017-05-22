#include "ConstitutiveRelation.h"

// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateIncrement_Elastic(double dE, double dNu, double d6StrainIncrement[6])
{
	double dConstant = dE/(1.0 + dNu)/(1.0 - 2.0*dNu);

	this->d6StressIncrement[0] = dConstant * ((1.0-dNu)*d6StrainIncrement[0] + dNu*d6StrainIncrement[1] + dNu*d6StrainIncrement[2]);
	this->d6StressIncrement[1] = dConstant * ((1.0-dNu)*d6StrainIncrement[1] + dNu*d6StrainIncrement[0] + dNu*d6StrainIncrement[2]);
	this->d6StressIncrement[2] = dConstant * ((1.0-dNu)*d6StrainIncrement[2] + dNu*d6StrainIncrement[0] + dNu*d6StrainIncrement[1]);

	this->d6StressIncrement[3] = dConstant * ((0.5-dNu)*d6StrainIncrement[3]);
	this->d6StressIncrement[4] = dConstant * ((0.5-dNu)*d6StrainIncrement[4]);
	this->d6StressIncrement[5] = dConstant * ((0.5-dNu)*d6StrainIncrement[5]);
}
// ----------------------------------------------------------------------------
