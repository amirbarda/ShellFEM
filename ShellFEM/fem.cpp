#include <assert.h>
#include <igl\triangle_triangle_adjacency.h>
#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>

#include "fem.h"
#include "s3element.h"

#define ABSENT_VERTEX -1
#define FIXED_NODE	  -1

bool elementHasDOF(IntList &verticesIndices, IntList const &vertexToDOFTranslationMap) {
	for (int vertexIdx: verticesIndices) {
		if (vertexIdx != ABSENT_VERTEX && vertexToDOFTranslationMap[vertexIdx] != FIXED_NODE) return true;
	}
	return false;
}

int calcNbrOppositeVrtxIndx(Mesh const &mesh, int currNbrFace, int faceIdx, int vertex) {
	int sumOfVertices = mesh.F.row(currNbrFace).sum();
	int sumOfSharedVertices = mesh.F(faceIdx, vertex) + mesh.F(faceIdx, (vertex + 1) % 3);
	return sumOfVertices - sumOfSharedVertices;
}

void setNbrsEnvelope(Mesh const &mesh, Element &currElement, IntList &verticesIndices, int faceIdx, MatrixXd &TT) {
	for (int vertexIdx = 0; vertexIdx < 3; vertexIdx++) {        
		int currNbrFace = TT(faceIdx, vertexIdx);
		if (currNbrFace == ABSENT_VERTEX) continue; //case: no neighbour face sharing edge j of current face i .

		int currNbrOppositeVrtxIndx = calcNbrOppositeVrtxIndx(mesh, currNbrFace, faceIdx, vertexIdx);
		currElement.setNeighbour(mesh.V.row(currNbrOppositeVrtxIndx), (vertexIdx + 2) % 3);
		verticesIndices[3 + vertexIdx] = currNbrOppositeVrtxIndx;
	}
}

void addKeToK(TriList &K_triplets, Element &currElement, IntList const &vertexToDOFTranslationMap, IntList &verticesIndices) {
	for (int row = 0; row < currElement.Ke.rows(); row++) {
		int KRowVertexIdx = verticesIndices[(int)(row / 3)];
		if (KRowVertexIdx == ABSENT_VERTEX || vertexToDOFTranslationMap[KRowVertexIdx] == FIXED_NODE) continue;
		int KRowDOFIdx = vertexToDOFTranslationMap[KRowVertexIdx] * 3 + row % 3;

		for (int col = 0; col < currElement.Ke.cols(); col++) {
			int KColVertexIdx = verticesIndices[(int)(col / 3)];
			if (KColVertexIdx == ABSENT_VERTEX || vertexToDOFTranslationMap[KColVertexIdx] == FIXED_NODE) continue;
			int KColDOFIdx = vertexToDOFTranslationMap[KColVertexIdx] * 3 + col % 3;
			K_triplets.push_back(TripletXd(KRowDOFIdx, KColDOFIdx, currElement.Ke(row, col))); // TODO can sum only half of Ke since it is symmetric.
		}
	}
}

void setFixedEdges(Element &currElement, VectorXi const &face, IntList const &vertexToDOFTranslationMap) {
	for (int vertexIdx = 0; vertexIdx < 3; vertexIdx++) {
		int vertex = face(vertexIdx);
		if (vertexToDOFTranslationMap[vertex] == FIXED_NODE) {
			currElement.setFixedNode(vertexIdx); 
		}
	}
}

Element createFaceElement(Mesh const &mesh, IntList &verticesIndices, int faceIdx) {
	Vector3d vertices[3];

	for (int vertexIdx = 0; vertexIdx < 3; vertexIdx++) {
		int vertex = mesh.F(faceIdx, vertexIdx);
		vertices[vertexIdx] = mesh.V.row(vertex);
		verticesIndices[vertexIdx] = vertex;
	}
	return Element(vertices);
}

/**	1. TT  : #F by #3 adjacent matrix. Description:
	   TT(i,j) = id of the neighbour triangle, that shares the j edge of triangle i .
	For a triangle, the 1st edge is [0,1] , the 2nd edge is [1,2] , the 3rd edge [2,3]. */
void calculateGlobalStiffnessMatrix(TriList &K_triplets, ElementBuilder &elementBuilder, Mesh const &mesh, IntList const &vertexToDOFTranslationMap) {
	MatrixXd TT; 
	MatrixXd TTi; // we don't use it.

	igl::triangle_triangle_adjacency(mesh.F, TT, TTi);

	for (int faceIdx = 0; faceIdx < mesh.F.rows(); faceIdx++) {
		IntList verticesIndices(6, ABSENT_VERTEX);
		Element currElement = createFaceElement(mesh, verticesIndices, faceIdx);
		setFixedEdges(currElement, mesh.F.row(faceIdx), vertexToDOFTranslationMap);
		setNbrsEnvelope(mesh, currElement, verticesIndices, faceIdx, TT);
		if (!elementHasDOF(verticesIndices, vertexToDOFTranslationMap)) continue; // skip

		elementBuilder.calculateStiffnessMatrix(currElement);
		addKeToK(K_triplets, currElement, vertexToDOFTranslationMap, verticesIndices);
	}
}

void calculateElasticPlasticMatrix(Eigen::Matrix<double, 3, 3> &De, FEMData const &data) {
	double ni = data.matProps.ni;
	double E = data.matProps.E;
	double t = data.shellProps.thickness;
	De <<  1,   ni,  0,
		   ni,  1,   0,
		   0,   0,   (1 - ni) / 2;
	De *= E / (1 - ni*ni);
	De *= t*t*t / 12; 
}

void calculateHillPlasticStrainMatrix(Eigen::Matrix<double, 3, 3> &M, FEMData const &data) {
	double f = 0.5, g = 0.5, h = 0.5, n = 1.5; // TODO get from data
	M << g + h, -h, 0,
		-h, f + h, 0,
		0, 0, 2 * n;
}

/*
* for each vertex keep the index of DOF for Ux,Uy,Uz.
* Ux = Map[v]*3, Uy = Map[v]*3+1, Uz = Map[v]*3+2
*/
void createVertexDOFTranslationMap(Mesh const &mesh, IntList &vertexToDOFTranslationMap) {
	vertexToDOFTranslationMap = IntList(mesh.V.rows(), FIXED_NODE);

	int j = 0, offset = 0;
	for (int i = 0; i < mesh.fixedNode.size(); i++) {
		for (; j < mesh.fixedNode[i] - 1; j++) {
			vertexToDOFTranslationMap[j] = offset++; 
		}
		j++;
	}
	for (; j < mesh.V.rows(); j++) {
		vertexToDOFTranslationMap[j] = offset++;
	}
}

void preproccessForSolver(Mesh const &mesh, SparseMat &K, TriList &K_triplets, ForcesList const &nodalForces, MatrixXd &forces, IntList &vertexToDOFTranslationMap) {
	int amountOfnotFixedVertices = mesh.V.rows() - mesh.fixedNode.size();

	// prepare K sparse matrix
	K = SparseMat(3 * amountOfnotFixedVertices, 3 * amountOfnotFixedVertices);
	K.setFromTriplets(K_triplets.begin(), K_triplets.end());

	// prepare forces matrix
	forces = MatrixXd::Zero(3 * amountOfnotFixedVertices, 1);
	for (int i = 0; i < nodalForces.size(); i++) {
		int indexOfDOF = vertexToDOFTranslationMap[nodalForces[i].first-1] * 3;
		forces.block(indexOfDOF, 0, 3, 1) = nodalForces[i].second;
	}

	std::cout << "Stiffness Matrix:" << std::endl << K << std::endl;
	std::cout << "Force Vector:" << std::endl << forces << std::endl;
}

bool solveSparseEquation(Mesh const &mesh, TriList &K_triplets, ForcesList const &nodalForces, MatrixXd &displacements, IntList &vertexToDOFTranslationMap) {
	SparseMat K;								// Global Stiffness Matrix.
	MatrixXd forces;							// Forces vector for solver
	Eigen::ConjugateGradient<SparseMat> solver; // solve K*U = F
	//Eigen::LeastSquaresConjugateGradient<SparseMat> solver;
	
	preproccessForSolver(mesh, K, K_triplets, nodalForces, forces, vertexToDOFTranslationMap);

	solver.compute(K);
	displacements = solver.solve(forces);

	std::cout << "# Iterations: " << solver.iterations() << std::endl;
	std::cout << "Estimated Error: " << solver.error() << std::endl;

	if (solver.info() != Eigen::Success) {
		std::cout << "Failed to solve" << std::endl;
		return true;
	}
	return false;
}

// TODO: should keep vertices in element in a matrix of size |V|x6
void calcStressFromDisplacements(Mesh const &mesh, ElementBuilder &elementBuilder, FEMResults &results, IntList const &vertexToDOFTranslationMap) { 
	Eigen::Matrix<double, 18, 1> elementDisplacements;
	MatrixXd TT;
	MatrixXd TTi; // we don't use it.

	igl::triangle_triangle_adjacency(mesh.F, TT, TTi);
	
	results.vonMisesStress = VectorXd::Zero(mesh.F.rows());

	for (int faceIdx = 0; faceIdx < mesh.F.rows(); faceIdx++) {
		IntList verticesIndices(6, ABSENT_VERTEX);
		Element currElement = createFaceElement(mesh, verticesIndices, faceIdx);
		setNbrsEnvelope(mesh, currElement, verticesIndices, faceIdx, TT);
		if (!elementHasDOF(verticesIndices, vertexToDOFTranslationMap)) continue; // skip
		
		elementDisplacements = Eigen::MatrixXd::Zero(18, 1);
		for (int i = 0; i < 6; i++) {
			if (verticesIndices[i] == ABSENT_VERTEX) continue;
			int dof = vertexToDOFTranslationMap[verticesIndices[i]];
			if (dof != FIXED_NODE) {
				elementDisplacements.block(i * 3, 0, 3, 1) = results.displacements.block(dof * 3, 0, 3, 1);
			}
		}

		elementBuilder.calculateVonMisesStress(currElement, elementDisplacements);
		results.vonMisesStress(faceIdx) = currElement.vonMisesStress;
	}
}

void getDisplacedMesh(Mesh const &mesh, FEMResults &results, IntList const &vertexToDOFTranslationMap) {
	results.displacedVertices = mesh.V;

	for (int vertexIdx = 0; vertexIdx < mesh.V.rows(); vertexIdx++) {
		int dof = vertexToDOFTranslationMap[vertexIdx];
		if (dof != ABSENT_VERTEX) {
			results.displacedVertices.row(vertexIdx) += results.displacements.block(dof*3, 0, 3, 1).transpose();
		}
	}
}

void printSummary(Mesh const &mesh, FEMResults &results) { 
	std::cout << "Num of Vertices: " << mesh.V.rows() << std::endl;
	std::cout << "Num of Triangles: " << mesh.F.rows() << std::endl;
	std::cout << "Displacements: " << std::endl << results.displacements << std::endl;
	std::cout << "Displaced Vertices:" << std::endl << results.displacedVertices << std::endl;
	std::cout << "Von-Mises Stress:" << std::endl << results.vonMisesStress << std::endl;
}

void init(Mesh const &mesh, FEMData const &data, ElementBuilder &elementBuilder, IntList &vertexToDOFTranslationMap) {
	Eigen::Matrix<double, 3, 3> De;									// Elastic Plastic Behavior Matrix
	Eigen::Matrix<double, 3, 3> M;									// Hill's Plastic Strain Matrix

	calculateElasticPlasticMatrix(De, data);
	calculateHillPlasticStrainMatrix(M, data);
	elementBuilder = ElementBuilder(De, M, data.shellProps.thickness);
	createVertexDOFTranslationMap(mesh, vertexToDOFTranslationMap);
}

// TODO : sould run remove_duplicates before?
void Perform_FEM(Mesh const &mesh, ForcesList const &nodalForces, FEMData const &data, FEMResults &results) { // TODO should not be void

	IntList vertexToDOFTranslationMap;								// mapping from vertex Idx to DOF
	TriList K_triplets;												// List used to init the sparse global stiffness matrix
	ElementBuilder elementBuilder;
	
	init(mesh, data, elementBuilder, vertexToDOFTranslationMap);
	std::cout << "Beginnig Assembly Of Stiffness Matrix." << std::endl;
	calculateGlobalStiffnessMatrix(K_triplets, elementBuilder, mesh, vertexToDOFTranslationMap);

	std::cout << "Beginnig Solve" << std::endl;
	if (solveSparseEquation(mesh, K_triplets, nodalForces, results.displacements, vertexToDOFTranslationMap)) {
		return; // fail
	}
	
	calcStressFromDisplacements(mesh, elementBuilder, results, vertexToDOFTranslationMap); 
	getDisplacedMesh(mesh, results, vertexToDOFTranslationMap);
	printSummary(mesh, results);

	return;
}

