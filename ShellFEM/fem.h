#pragma once
#include "utils.h"
#include "s3element.h"

struct FEMResults {
	MatrixXd displacements;
	VectorXd vonMisesStress;
	MatrixXd displacedVertices;
};

class FEMSimulation {
private:
	void createDOFTranslationMap(Mesh &mesh, MatrixXd &DOFTranslationMap);
	void calculateGlobalStiffnessMatrix(Mesh &mesh, MatrixXd &DOFTranslationMap, TriList &K_triplets, ElementBuilder &elementBuilder);
	Element createFaceElement(Mesh &mesh, IntList &verticesIndices, int faceIdx);
	int calcNbrOppositeVrtxIndx(Mesh const &mesh, int currNbrFace, int faceIdx, int vertex);
	void setNbrsEnvelope(Mesh &mesh, Element &currElement, IntList &verticesIndices, int faceIdx, MatrixXd &TT);
	void addKeToK(Element &currElement, IntList &verticesIndices, MatrixXd &DOFTranslationMap, TriList &K_triplets);
	bool solveSparseEquation(MatrixXd &DOFTranslationMap, TriList &K_triplets, vector3dList const &nodalForces, FEMResults &results);
	void preproccessForSolver(SparseMat &K, MatrixXd &forces, MatrixXd &DOFTranslationMap, TriList &K_triplets, vector3dList const &nodalForces);
	void calcStressFromDisplacements(Mesh &mesh, MatrixXd &DOFTranslationMap, ElementBuilder &elementBuilder, FEMResults &results);
	void getDisplacedMesh(Mesh &mesh, MatrixXd &DOFTranslationMap, FEMResults &results);
	void printSummary(Mesh &mesh, FEMResults &results);
public:
	FEMSimulation() {};
	bool performFEM(Mesh &mesh, vector3dList const &nodalForces, ElementBuilder &elementBuilder, FEMResults &results);
};