#include <assert.h>
#include <igl\triangle_triangle_adjacency.h>
#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>

#include "fem.h"

#define ABSENT_VERTEX -1
#define FIXED_NODE	  -1
#define NO_NIGHBOR    -1

/*
bool elementHasDOF(IntList &verticesIndices, MatrixXd const &DOFTranslationMap) {
	for (int vertexIdx: verticesIndices) {
		if (vertexIdx != ABSENT_VERTEX && DOFTranslationMap[vertexIdx] != FIXED_NODE) return true;
	}
	return false;
}
*/

int calcNbrOppositeVrtxIndx(Mesh const &mesh, int nbrFaceIdx, int faceIdx, int edgeIdx) {
	int sumOfVertices = mesh.F.row(nbrFaceIdx).sum();
	int sumOfSharedVertices = mesh.F(faceIdx, edgeIdx) + mesh.F(faceIdx, (edgeIdx + 1) % 3);
	return sumOfVertices - sumOfSharedVertices;
}

void setNbrsEnvelope(Mesh &mesh, Element &currElement, IntList &verticesIndices, int faceIdx, MatrixXd &TT) {
	for (int edgeIdx = 0; edgeIdx < 3; edgeIdx++) {
		int nbrFaceIdx = TT(faceIdx, edgeIdx);
		if (nbrFaceIdx == NO_NIGHBOR) continue; //case: no neighbour face sharing edge at index edgeIdx.

		int currNbrOppositeVrtxIndx = calcNbrOppositeVrtxIndx(mesh, nbrFaceIdx, faceIdx, edgeIdx);
		currElement.setNeighbour(mesh.V.row(currNbrOppositeVrtxIndx), edgeIdx);
		verticesIndices[3 + edgeIdx] = currNbrOppositeVrtxIndx;
	}
}

void addKeToK(Element &currElement, IntList &verticesIndices, MatrixXd &DOFTranslationMap, TriList &K_triplets) {
	for (int row = 0; row < currElement.Ke.rows(); row++) {
		int KRowVertexIdx = verticesIndices[(int)(row / 3)];
		int KRowVertexAxis = row % 3;
		int KRowDOFIdx = DOFTranslationMap(KRowVertexIdx, KRowVertexAxis);
		if (KRowVertexIdx == ABSENT_VERTEX || KRowDOFIdx == FIXED_NODE) continue;

		for (int col = 0; col < currElement.Ke.cols(); col++) { 
			int KColVertexIdx = verticesIndices[(int)(col / 3)];
			int KColVertexAxis = col % 3;
			int KColDOFIdx = DOFTranslationMap(KColVertexIdx, KColVertexAxis);
			if (KColVertexIdx == ABSENT_VERTEX || KColDOFIdx == FIXED_NODE) continue;
			K_triplets.push_back(TripletXd(KRowDOFIdx, KColDOFIdx, currElement.Ke(row, col)));
		}
	}
}

void setClampedEdges(Mesh &mesh, Element &currElement, int faceIdx) {
	for (int i = 0; i < 3; i++) {
		int edge = mesh.FE(faceIdx, i);
		currElement.isEdgeClamped[i] = mesh.isEdgeClamped[edge];
	}
}

Element createFaceElement(Mesh &mesh, IntList &verticesIndices, int faceIdx) {
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
void calculateGlobalStiffnessMatrix(Mesh &mesh, MatrixXd &DOFTranslationMap, TriList &K_triplets, ElementBuilder &elementBuilder) {
	MatrixXd TT; 
	MatrixXd TTi; // we don't use it.

	igl::triangle_triangle_adjacency(mesh.F, TT, TTi);

	for (int faceIdx = 0; faceIdx < mesh.F.rows(); faceIdx++) {
		IntList verticesIndices(6, ABSENT_VERTEX);
		Element currElement = createFaceElement(mesh, verticesIndices, faceIdx);
		setClampedEdges(mesh, currElement, faceIdx); 
		setNbrsEnvelope(mesh, currElement, verticesIndices, faceIdx, TT);
		//if (!elementHasDOF(verticesIndices, DOFTranslationMap)) continue; // skip // FIXME

		elementBuilder.calculateStiffnessMatrix(currElement);
		addKeToK(currElement, verticesIndices, DOFTranslationMap, K_triplets);
	}
}

/*
* for each vertex keep the index of DOF for Ux,Uy,Uz.
* Ux(v) = Map[v][0], Uy(v) = Map[v][1], Uz(v) = Map[v][2]
*/
void createDOFTranslationMap(Mesh &mesh, MatrixXd &DOFTranslationMap) {
	DOFTranslationMap = MatrixXd::Constant(mesh.V.rows(), 3, FIXED_NODE);

	int j = 0, offset = 0;
	for (int i = 0; i < mesh.freeDOF.size(); i++) {
		int vertexIdx = mesh.freeDOF[i].first - 1;
		for (; j < vertexIdx; j++) {
			DOFTranslationMap.row(j) = Vector3d(offset, offset+1, offset+2);
			offset += 3;
		}
		for (int k = 0; k < 3; k++) {
			if (mesh.freeDOF[i].second(k)) 
				DOFTranslationMap(j, k) = offset++;
		}
		j++;
	}
	for (; j < mesh.V.rows(); j++) {
		DOFTranslationMap.row(j) = Vector3d(offset, offset + 1, offset + 2);
		offset += 3;
	}
}

void preproccessForSolver(SparseMat &K, MatrixXd &forces, MatrixXd &DOFTranslationMap, TriList &K_triplets, vector3dList const &nodalForces) {
	int numOfDOF = 0; // TODO

	for (int i = 0; i < DOFTranslationMap.rows(); i++) { // need to save this as class parameter
		for (int j = 0; j < 3; j++) {
			if (DOFTranslationMap(i, j) > numOfDOF)
				numOfDOF = DOFTranslationMap(i, j);
		}
	}
	numOfDOF++;

	// prepare K sparse matrix
	K = SparseMat(numOfDOF, numOfDOF);
	K.setFromTriplets(K_triplets.begin(), K_triplets.end());

	// prepare forces matrix
	forces = MatrixXd::Zero(numOfDOF, 1); // TODO why not vector?
	for (int i = 0; i < nodalForces.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int indexOfDOF = DOFTranslationMap(nodalForces[i].first - 1, j);
			if (indexOfDOF != FIXED_NODE) {
				forces(indexOfDOF, 0) = nodalForces[i].second(j);
			}
		}
	}

	std::cout << "Stiffness Matrix:" << std::endl << K << std::endl;
	std::cout << "Force Vector:" << std::endl << forces << std::endl;
}

bool solveSparseEquation(MatrixXd &DOFTranslationMap, TriList &K_triplets, vector3dList const &nodalForces, FEMResults &results) {
	SparseMat K;								// Global Stiffness Matrix.
	MatrixXd forces;							// Forces vector for solver
	Eigen::ConjugateGradient<SparseMat> solver; // solve K*U = F. 
	
	preproccessForSolver(K, forces, DOFTranslationMap, K_triplets, nodalForces);

	solver.compute(K);
	results.displacements = solver.solve(forces);

	std::cout << "# Iterations: " << solver.iterations() << std::endl;
	std::cout << "Estimated Error: " << solver.error() << std::endl;

	if (solver.info() != Eigen::Success) {
		std::cout << "Failed to solve" << std::endl;
		return true;
	}
	return false;
}

// TODO: should keep vertices in element in a matrix of size |V|x6
void calcStressFromDisplacements(Mesh &mesh, MatrixXd &DOFTranslationMap, ElementBuilder &elementBuilder, FEMResults &results) { // FIXME : not the same result as before
	Eigen::Matrix<double, 18, 1> elementDisplacements;
	MatrixXd TT;
	MatrixXd TTi; // we don't use it.

	igl::triangle_triangle_adjacency(mesh.F, TT, TTi);
	
	results.vonMisesStress = VectorXd::Zero(mesh.F.rows());

	for (int faceIdx = 0; faceIdx < mesh.F.rows(); faceIdx++) {
		IntList verticesIndices(6, ABSENT_VERTEX);
		Element currElement = createFaceElement(mesh, verticesIndices, faceIdx);
		setClampedEdges(mesh, currElement, faceIdx);
		setNbrsEnvelope(mesh, currElement, verticesIndices, faceIdx, TT);
		//if (!elementHasDOF(verticesIndices, DOFTranslationMap)) continue; // skip // FIXME
		
		elementDisplacements = Eigen::MatrixXd::Zero(18, 1);
		for (int i = 0; i < 6; i++) {
			if (verticesIndices[i] == ABSENT_VERTEX) continue;
			for (int j = 0; j < 3; j++) {
				int dof = DOFTranslationMap(verticesIndices[i], j);
				if (dof != FIXED_NODE) {
					elementDisplacements(i*3+j, 0) = results.displacements(dof, 0);
				}
			}
		}

		elementBuilder.calculateVonMisesStress(currElement, elementDisplacements);
		results.vonMisesStress(faceIdx) = currElement.vonMisesStress;
	}
}

void getDisplacedMesh(Mesh &mesh, MatrixXd &DOFTranslationMap, FEMResults &results) {
	results.displacedVertices = mesh.V;

	for (int vertexIdx = 0; vertexIdx < mesh.V.rows(); vertexIdx++) {
		for (int j = 0; j < 3; j++) {
			int dof = DOFTranslationMap(vertexIdx, j);
			if (dof != ABSENT_VERTEX) {
				results.displacedVertices(vertexIdx, j) += results.displacements(dof, 0);
			}
		}
	}
}

void printSummary(Mesh &mesh, FEMResults &results) {
	//TODO print num of dof
	std::cout << "Num of Vertices: " << mesh.V.rows() << std::endl;
	std::cout << "Num of Triangles: " << mesh.F.rows() << std::endl;
	std::cout << "Displacements: " << std::endl << results.displacements << std::endl;
	std::cout << "Displaced Vertices:" << std::endl << results.displacedVertices << std::endl;
	std::cout << "Von-Mises Stress:" << std::endl << results.vonMisesStress << std::endl;
}

// TODO : sould run remove_duplicates before?
bool performFEM(Mesh &mesh, vector3dList const &nodalForces, SimulationProperties &simProps, FEMResults &results) {
	ElementBuilder elementBuilder(simProps);
	MatrixXd DOFTranslationMap;					// [|V|x3] matrix of mapping of vertex Ux,Uy,Uz to dof #
	TriList K_triplets;							// List used to init the sparse global stiffness matrix

	createDOFTranslationMap(mesh, DOFTranslationMap);
	std::cout << "Beginnig Assembly Of Stiffness Matrix." << std::endl;
	calculateGlobalStiffnessMatrix(mesh, DOFTranslationMap, K_triplets, elementBuilder);

	std::cout << "Beginnig Solve" << std::endl;
	if (solveSparseEquation(DOFTranslationMap,K_triplets, nodalForces, results)) {
		return true; // fail
	}

	calcStressFromDisplacements(mesh, DOFTranslationMap, elementBuilder, results);
	getDisplacedMesh(mesh, DOFTranslationMap, results);
	printSummary(mesh, results);

	return false;
}

