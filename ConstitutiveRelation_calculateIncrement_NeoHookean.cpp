#include "ConstitutiveRelation.h"

// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateIncrement_NeoHookean_6D(double dE, double dNu, double d6Stress_Current[6], glm::dmat3 d33DeformationGradient)
{// (2011) A convected particle domain interpolation technique to extend
	double dLame1 = 0.5*dE / (1.0 + dNu);
	double dLame2 = dNu*dE / ((1.0+dNu)*(1.0-2.0*dNu));

	double dJ = glm::determinant(d33DeformationGradient);

	glm::dmat3 d33P = dLame2*log(dJ)/dJ * glm::dmat3(1.0) + dLame1/dJ*((d33DeformationGradient)*glm::transpose(d33DeformationGradient)-glm::dmat3(1.0));

	for(int index = 0; index < 6; index++)
		this->d6StressIncrement[index] = 0.0;

	this->d6StressIncrement[0] = d33P[0][0] - d6Stress_Current[0];
	this->d6StressIncrement[1] = d33P[1][1] - d6Stress_Current[1];
	this->d6StressIncrement[2] = d33P[2][2] - d6Stress_Current[2];
	this->d6StressIncrement[3] = d33P[0][1] - d6Stress_Current[3];
	this->d6StressIncrement[4] = d33P[1][2] - d6Stress_Current[4];
	this->d6StressIncrement[5] = d33P[2][0] - d6Stress_Current[5];
}
// ----------------------------------------------------------------------------
