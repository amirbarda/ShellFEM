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
	Eigen::MatrixXd displacement;
	Eigen::VectorXd vonMisesStress;
};


/**
	Performs FEM analysis with shell elements.
	@param meshedV nodes of the meshed part. each row is a node in (x,y,z) format
	@param meshedF faces of meshed part. each row is a triangular face in (v1,v2,v3)
	@param nodalForces
	@param fixedNodes
	@param data material, shell and simulation data
	@return generalized displacements (rotational + axial).
*/
FEMResults Perform_FEM(Eigen::MatrixXd const &MeshV, Eigen::MatrixXi const &MeshF,
	std::vector<Force> const &nodalForces, std::vector<int> const  &fixedNode, FEMData data);
