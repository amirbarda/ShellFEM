#pragma once
#include "utils.h"
#include "s3element.h"

#define DEFAULT_YOUNG_CNST		70e9
#define DEFAULT_POSSION_CNST	0.3
#define DEFAULT_THICKNESS_CNST	1.6e-3

struct SimulationProperties {
	std::string name, outDir, objPath, forcesPath, fixedPath;
	bool startViewer;
	bool isPropsSet;
	double E; //Young's Modulus
	double ni; //Possion's Ratio
	double thickness;
	SimulationProperties() : E(DEFAULT_YOUNG_CNST), ni(DEFAULT_POSSION_CNST), thickness(DEFAULT_THICKNESS_CNST), startViewer(true), isPropsSet(false){};
	SimulationProperties(std::string name_, std::string outDir_, std::string objPath_, std::string forcesPath_, std::string fixedPath_, bool startViewer_,
		double E_, double ni_, double thickness_) :
		name(name_), outDir(outDir_), objPath(objPath_), forcesPath(forcesPath_), fixedPath(fixedPath_), startViewer(startViewer_),
		E(E_), ni(ni_), thickness(thickness_), isPropsSet(true) {};
};

struct FEMResults {
	MatrixXd displacements;
	VectorXd vonMisesStress;
	MatrixXd displacedVertices;
};

class FEMSimulation {
private:
	Mesh const mesh;
	vector3dList const nodalForces;
	SimulationProperties const simProps;
	MatrixXd DOFTranslationMap;					// [|V|x3] matrix of mapping of vertex Ux,Uy,Uz to dof #
	TriList K_triplets;							// List used to init the sparse global stiffness matrix
	ElementBuilder elementBuilder;

	void init();
	void createDOFTranslationMap();
	void calculateGlobalStiffnessMatrix();
	Element createFaceElement(IntList &verticesIndices, int faceIdx);
	void setNbrsEnvelope(Element &currElement, IntList &verticesIndices, int faceIdx, MatrixXd &TT);
	void addKeToK(Element &currElement, IntList &verticesIndices);
	bool solveSparseEquation();
	void preproccessForSolver(SparseMat &K, MatrixXd &forces);
	void calcStressFromDisplacements();
	void getDisplacedMesh();
	void printSummary();
public:
	FEMResults results;
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
	FEMSimulation(Mesh const &mesh, vector3dList const &nodalForces) : mesh(mesh), nodalForces(nodalForces) {};
	FEMSimulation(Mesh const &mesh, vector3dList const &nodalForces, SimulationProperties const &simProps): mesh(mesh), nodalForces(nodalForces), simProps(simProps){};
	void performFEM();
};