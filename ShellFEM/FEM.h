#pragma once
#include "utils.h"

struct MaterialProperty {
	double E; //Young's Modulus [N/mm^2]
	double ni; //Possions Ratio
	MaterialProperty() { E = 70e9; ni = 0.3; };
	MaterialProperty(double E_, double ni_) :E(E_), ni(ni_){};
};

struct ShellProperty {
	double thickness; //[mm]
	ShellProperty() { thickness = 1.6e-3; }
	ShellProperty(double thickenss_) :thickness(thickenss_) {};
};

struct SimulationProperties{
	std::string name, output_dir;
	SimulationProperties() { name = "dummy"; output_dir = ""; };
	SimulationProperties(std::string name_, std::string output_dir_) :name(name_),output_dir(output_dir_) {};
};

struct FEMData {
	MaterialProperty matProps;
	ShellProperty shellProps;
	SimulationProperties simuProps;
};

struct FEMResults {
	MatrixXd displacements;
	VectorXd vonMisesStress;
	MatrixXd displacedVertices;
};

/**
	Performs FEM analysis with shell elements.
	@param mesh:
		V: nodes of the meshed part. each row is a node in (x,y,z) format
		F: faces of meshed part. each row is a triangular face in (v1,v2,v3)
		fixedNodes: list of fixed node indices
	@param nodalForces
	@param data material, shell and simulation data
	@return generalized displacements (rotational + axial).
*/
void Perform_FEM(Mesh const &mesh, ForcesList const &nodalForces, FEMData const &data, FEMResults &results);
