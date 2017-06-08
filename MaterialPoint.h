#ifndef MATERIALPOINT_H
#define MATERIALPOINT_H

#include "Definitions.h"
////#include ".\include\glm\glm.hpp" // windows
//#include "./include/glm/glm.hpp" // linux

enum MaterialType
{
	_ELASTIC,
	_NEOHOOKEAN,
	_VISCOELASTIC,
	_PLASTIC,
	_VONMISESHARDENING,
	_GASS,
};

class MaterialPoint_Kinetics
{
	public:
		MaterialPoint_Kinetics() {;}
		virtual ~MaterialPoint_Kinetics() {;}

		bool b_DisplacementControl = false;
		bool b_Mark_Force = false;
		bool b_Surface = false;

		unsigned int i_ID = 0;

		double d_Mass = 0.0;

		glm::dvec3 d3_Position = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Force_External = glm::dvec3(0.0, 0.0, 0.0);
	protected:

	private:
};

class MaterialPoint_Material
{
	public:
		MaterialPoint_Material() {;}
		virtual ~MaterialPoint_Material() {;}

		bool b_Mark_Stress = false;

		unsigned int i_ID = 0;
		unsigned int i_MaterialType = _ELASTIC;

		double d_Volume = 0.0;
		double d_Volume_Initial = 0.0;

		double d_ElasticModulus = 0.0;
		double d_PoissonRatio = 0.0;
		double d_YieldStress = 0.0;
		double d_Viscosity = 0.0;

		glm::dmat3 d33_DeformationGradient = glm::dmat3(1.0);

		double d6_Strain[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double d6_Stress[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double d6_Strain_Plastic[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

		double d_BackStress_Isotropic = 0.0;
		double d_Hardening_Isotropic_C0 = 0.0;
		double d_Hardening_Isotropic_C1 = 0.0;
	protected:

	private:
};

class MaterialPoint
{
	public:
		MaterialPoint() {;}
		virtual ~MaterialPoint() {;}

		bool b_DisplacementControl = false;
		bool b_Mark_Force = false;
		bool b_Mark_Stress = false;

		unsigned int i_ID = 0;
		unsigned int i_MaterialType = _ELASTIC;

		double d_Mass = 0.0;
		double d_Volume = 0.0;
		double d_Volume_Initial = 0.0;

		double d_ElasticModulus = 0.0;
		double d_PoissonRatio = 0.0;
		double d_YieldStress = 0.0;
		double d_Viscosity = 0.0;

		bool	b_Surface = false;
		double		d_Kernel = 0.0;
		glm::dvec3	d3_Kernel_Gradient = glm::dvec3(0.0, 0.0, 0.0);

		glm::dvec3 d3_Position = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Force_External = glm::dvec3(0.0, 0.0, 0.0);

		glm::dmat3 d33_DeformationGradient = glm::dmat3(1.0);

		double d6_Strain[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double d6_Stress[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double d6_Strain_Plastic[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

		double d_BackStress_Isotropic = 0.0;
		double d_Hardening_Isotropic_C0 = 0.0;
		double d_Hardening_Isotropic_C1 = 0.0;
	protected:

	private:
};

#endif // MATERIALPOINT_H
