#include <assert.h>
#include <igl\triangle_triangle_adjacency.h>
#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>

#include "fem.h"

//######################################################################### Initialization ############################################################################

/*
* for each vertex keep the index of DOF for Ux,Uy,Uz.
* Ux(v) = Map[v][0], Uy(v) = Map[v][1], Uz(v) = Map[v][2]
*/
void createDOFTranslationMap(Mesh &mesh, MatrixXd &DOFTranslationMap, int &numOfDOF) {
	std::cout << "In createDOFTranslationMap" << std::endl;
	DOFTranslationMap = MatrixXd::Constant(mesh.V.rows(), 3, FIXED_NODE); // init all dof as fixed

	int j = 0, offset = 0;
	for (int i = 0; i < mesh.freeDOF.size(); i++) {
		// set dof up to next fixed node
		int vertexIdx = mesh.freeDOF[i].first;
		for (; j < vertexIdx; j++) {
			DOFTranslationMap.row(j) = Vector3d(offset, offset + 1, offset + 2);
			offset += 3;
		}
		// set the fixed node dof's according to flag for each axis 
		for (int k = 0; k < 3; k++) {
			if (mesh.freeDOF[i].second(k))
				DOFTranslationMap(j, k) = offset++;
		}
		// go to next
		j++;
	}
	// set dof for the rest of the nodes
	for (; j < mesh.V.rows(); j++) {
		DOFTranslationMap.row(j) = Vector3d(offset, offset + 1, offset + 2);
		offset += 3;
	}
	numOfDOF = offset;

	std::cout << "Number of Vertices: " << mesh.V.rows() << std::endl;
	std::cout << "Number of Faces: " << mesh.FNB.rows() << std::endl;
	std::cout << "Number of DOF: " << numOfDOF << std::endl;
	std::cout << "DOFTranslationMap" << std::endl;
	std::cout << DOFTranslationMap << std::endl;
	std::cout << DASH << std::endl;
}

//######################################################################### Stiffness Matrix ############################################################################

int getDOFIdx(int idx, VectorXi &verticesIndices, MatrixXd &DOFTranslationMap) {
	int vertexIdx = verticesIndices[(int)(idx / 3)];
	if (vertexIdx == ABSENT_VERTEX) return ABSENT_VERTEX;
	int axis = idx % 3;
	return DOFTranslationMap(vertexIdx, axis);
}

void addKeToK(Element &currElement, VectorXi &verticesIndices, MatrixXd &DOFTranslationMap, TriList &K_triplets) {
	for (int row = 0; row < currElement.Ke.rows(); row++) {
		int rowDOF = getDOFIdx(row, verticesIndices, DOFTranslationMap);
		if (rowDOF == FIXED_NODE) continue;

		for (int col = 0; col < currElement.Ke.cols(); col++) { 
			int colDOF = getDOFIdx(col, verticesIndices, DOFTranslationMap);
			if (colDOF == FIXED_NODE) continue;

			K_triplets.push_back(TripletXd(rowDOF, colDOF, currElement.Ke(row, col)));
			std::cout << "pushing to K_triplets. row: " << rowDOF << " col: " << colDOF << " value: " << currElement.Ke(row, col) << std::endl;
		}
	}
}

Element createFaceElement(Mesh &mesh, int faceIdx) {
	Vector3d vertices[3];

	// set vertices
	for (int i = 0; i < 3; i++) {
		int vertexIdx = mesh.FNB(faceIdx, i);
		vertices[i] = mesh.V.row(vertexIdx);
	}
	Element element = Element(vertices);

	// set clamped edges
	for (int i = 0; i < 3; i++) {
		int edge = mesh.FE(faceIdx, i);
		element.isEdgeClamped[i] = mesh.isEdgeClamped[edge];
		if (mesh.isEdgeClamped[edge])
			std::cout << "setting edge idx " << i << " as clamped" << std::endl;
	}

	// set neighbor vertices
	for (int i = 3; i < 6; i++) {
		int nbrVertexIdx = mesh.FNB(faceIdx, i);
		if (nbrVertexIdx != ABSENT_VERTEX)
			element.setNeighbour(mesh.V.row(nbrVertexIdx), i%3);
	}
	return element;
}

void calculateGlobalStiffnessMatrix(Mesh &mesh, MatrixXd &DOFTranslationMap, TriList &K_triplets, ElementBuilder &elementBuilder) { // TODO : skip element with 0 dof
	std::cout << "In calculateGlobalStiffnessMatrix" << std::endl;

	for (int faceIdx = 0; faceIdx < mesh.FNB.rows(); faceIdx++) {
		VectorXi verticesIndices = mesh.FNB.row(faceIdx);
		std::cout << "faceIdx: " << faceIdx << std::endl;
		std::cout << "verticesIndices" << std::endl << verticesIndices << std::endl;
		Element element = createFaceElement(mesh, faceIdx);
		elementBuilder.calculateStiffnessMatrix(element);
		addKeToK(element, verticesIndices, DOFTranslationMap, K_triplets);
		std::cout << DASH << std::endl;
	}
}

//######################################################################### Solver ############################################################################

void preproccessForSolver(SparseMat &K, MatrixXd &forces, MatrixXd &DOFTranslationMap, TriList &K_triplets, int numOfDOF, vector3dList const &nodalForces) {
	std::cout << "In preproccessForSolver" << std::endl;

	for (auto tp : K_triplets) std::cout << "row: " << tp.row() << " col: " << tp.col() << " value: " << tp.value() << std::endl;

	// prepare K sparse matrix
	K = SparseMat(numOfDOF, numOfDOF);
	K.setFromTriplets(K_triplets.begin(), K_triplets.end());

	// prepare forces matrix
	forces = MatrixXd::Zero(numOfDOF, 1); // TODO why not vector?
	for (int i = 0; i < nodalForces.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int indexOfDOF = DOFTranslationMap(nodalForces[i].first, j);
			if (indexOfDOF != FIXED_NODE) {
				forces(indexOfDOF, 0) = nodalForces[i].second(j);
			}
		}
	}

	std::cout << "Stiffness Matrix:" << std::endl << K << std::endl;
	std::cout << "Force Vector:" << std::endl << forces << std::endl;
	std::cout << DASH << std::endl;
}

void setDisplacements(Mesh &mesh, MatrixXd &DOFTranslationMap, FEMResults &results, VectorXd &disp) {
	std::cout << "In setDisplacements" << std::endl;

	results.displacements = MatrixXd::Zero(mesh.V.rows(), 3);
	for (int vertexIdx = 0; vertexIdx < mesh.V.rows(); vertexIdx++) {
		for (int j = 0; j < 3; j++) {
			int dof = DOFTranslationMap(vertexIdx, j);
			if (dof != ABSENT_VERTEX) {
				results.displacements(vertexIdx, j) = disp(dof, 0);
			}
		}
	}
	results.displacedVertices = mesh.V + results.displacements;
}

bool solveSparseEquation(Mesh &mesh, MatrixXd &DOFTranslationMap, TriList &K_triplets, vector3dList const &nodalForces, int numOfDOF, FEMResults &results) {
	SparseMat K;								// Global Stiffness Matrix.
	MatrixXd forces;							// Forces vector for solver
	Eigen::ConjugateGradient<SparseMat> solver; // solve K*U = F. 
	//Eigen::SparseLU<SparseMat> solver;
	VectorXd displacementOfDOFs;				// solution is displacments of dof

	std::cout << "In solveSparseEquation" << std::endl;
	
	preproccessForSolver(K, forces, DOFTranslationMap, K_triplets, numOfDOF, nodalForces);

	//Eigen::SparseQR <SparseMat, Eigen::COLAMDOrdering<int> > sqr;
	//sqr.compute(K);
	//std::cout << "K Rank: " << sqr.rank() << std::endl;

	solver.compute(K);
	displacementOfDOFs = solver.solve(forces); // Fixme need to set tolerance (epsilon)

	//std::cout << "# Iterations: " << solver.iterations() << std::endl;
	//std::cout << "Estimated Error: " << solver.error() << std::endl;

	if (solver.info() != Eigen::Success) {
		std::cout << "Failed to solve" << std::endl;
		if (solver.info() == Eigen::NumericalIssue) std::cout << "The provided data did not satisfy the prerequisites." << std::endl;
		if (solver.info() == Eigen::NoConvergence) std::cout << "Iterative procedure did not converge.." << std::endl;
		if (solver.info() == Eigen::InvalidInput) std::cout << "The inputs are invalid, or the algorithm has been improperly called. \
												When assertions are enabled, such errors trigger an assert.." << std::endl;
		return true;
	}
	setDisplacements(mesh, DOFTranslationMap, results, displacementOfDOFs);
	return false;
}

//#################################################################### Von Mises Stress ############################################################################

void calcStressFromDisplacements(Mesh &mesh, MatrixXd &DOFTranslationMap, ElementBuilder &elementBuilder, FEMResults &results) { 
	std::cout << "In calcStressFromDisplacements" << std::endl;
	results.vonMisesStress = VectorXd::Zero(mesh.FNB.rows());

	for (int faceIdx = 0; faceIdx < mesh.FNB.rows(); faceIdx++) {
		Eigen::Matrix<double, 18, 1> elementDisplacements = MatrixXd::Zero(18, 1);
		VectorXi verticesIndices = mesh.FNB.row(faceIdx);
		std::cout << "faceIdx: " << faceIdx << std::endl;
		std::cout << "verticesIndices" << verticesIndices << std::endl;
		std::cout << DASH << std::endl;

		Element element = createFaceElement(mesh, faceIdx);
		for (int i = 0; i < 6; i++) {
			int vertexIdx = verticesIndices[i];
			if (vertexIdx == ABSENT_VERTEX) continue;
			for (int j = 0; j < 3; j++) {
				elementDisplacements(i*3+j, 0) = results.displacements(vertexIdx, j);
			}
		}

		elementBuilder.calculateVonMisesStress(element, elementDisplacements);
		results.vonMisesStress(faceIdx) = element.vonMisesStress;
	}
}

//#################################################################### Main ############################################################################

void printSummary(Mesh &mesh, FEMResults &results) {
	std::cout << "Displacements: " << std::endl << results.displacements << std::endl;
	std::cout << DASH << std::endl;
	std::cout << "Displaced Vertices:" << std::endl << results.displacedVertices << std::endl;
	std::cout << DASH << std::endl;
	std::cout << "Von-Mises Stress:" << std::endl << results.vonMisesStress << std::endl;
	std::cout << DASH << std::endl;
}

bool performFEM(Mesh &mesh, vector3dList const &nodalForces, SimulationProperties &simProps, FEMResults &results) {
	ElementBuilder elementBuilder(simProps);
	MatrixXd DOFTranslationMap;					// [|V|x3] matrix of mapping of vertex Ux,Uy,Uz to dof #
	TriList K_triplets;							// List used to init the sparse global stiffness matrix
	int numOfDOF;

	createDOFTranslationMap(mesh, DOFTranslationMap, numOfDOF);
	std::cout << "Beginnig Assembly Of Stiffness Matrix." << std::endl;
	calculateGlobalStiffnessMatrix(mesh, DOFTranslationMap, K_triplets, elementBuilder);

	std::cout << "Beginnig Solve" << std::endl;
	if (solveSparseEquation(mesh, DOFTranslationMap, K_triplets, nodalForces, numOfDOF, results)) {
		return true; // fail
	}
	calcStressFromDisplacements(mesh, DOFTranslationMap, elementBuilder, results);
	printSummary(mesh, results);

	return false;
}

