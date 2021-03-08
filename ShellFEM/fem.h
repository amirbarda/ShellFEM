#pragma once
#include "utils.h"

#define DEFAULT_YOUNG_CNST		70e9
#define DEFAULT_POSSION_CNST	0.3
#define DEFAULT_THICKNESS_CNST	1.6e-3

struct MaterialProperty {
	double E; //Young's Modulus [N/mm^2]
	double ni; //Possions Ratio
	MaterialProperty() { E = DEFAULT_YOUNG_CNST; ni = DEFAULT_POSSION_CNST; };
	MaterialProperty(double E_, double ni_) :E(E_), ni(ni_){};
};

struct ShellProperty {
	double thickness; //[mm]
	ShellProperty() { thickness = DEFAULT_THICKNESS_CNST; }
	ShellProperty(double thickenss_) :thickness(thickenss_) {};
};

struct SimulationProperties {
	std::string name, outDir, objPath, forcesPath, fixedPath;
	bool startViewer;
	double E; //Young's Modulus
	double ni; //Possion's Ratio
	double thickness;
	SimulationProperties() : E(DEFAULT_YOUNG_CNST), ni(DEFAULT_POSSION_CNST), thickness(DEFAULT_THICKNESS_CNST), startViewer(true) {};
	SimulationProperties(std::string name_, std::string outDir_, std::string objPath_, std::string forcesPath_, std::string fixedPath_, bool startViewer_,
		double E_, double ni_, double thickness_) :
		name(name_), outDir(outDir_), objPath(objPath_), forcesPath(forcesPath_), fixedPath(fixedPath_), startViewer(startViewer_),
		E(E_), ni(ni_), thickness(thickness_) {};
};

struct FEMData {
	MaterialProperty matProps;
	ShellProperty shellProps;
	FEMData() { matProps = MaterialProperty(); shellProps = ShellProperty(); };
	FEMData(double E_, double ni_, double thickenss_) : matProps(E_, ni_), shellProps(thickenss_) {};
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
