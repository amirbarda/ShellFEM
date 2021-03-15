#include <assert.h>
#include <igl\triangle_triangle_adjacency.h>
#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>

#include "fem.h"
#include "s3element.h"

#define ABSENT_VERTEX -1
#define FIXED_NODE	  -1

/*
bool elementHasDOF(IntList &verticesIndices, MatrixXd const &DOFTranslationMap) {
	for (int vertexIdx: verticesIndices) {
		if (vertexIdx != ABSENT_VERTEX && DOFTranslationMap[vertexIdx] != FIXED_NODE) return true;
	}
	return false;
}
*/

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

void addKeToK(TriList &K_triplets, Element &currElement, MatrixXd const &DOFTranslationMap, IntList &verticesIndices) {
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

/*
void setFixedEdges(Element &currElement, VectorXi const &face, MatrixXd const &DOFTranslationMap) {
	for (int vertexIdx = 0; vertexIdx < 3; vertexIdx++) {
		int vertex = face(vertexIdx);
		if (DOFTranslationMap[vertex] == FIXED_NODE) {
			currElement.setFixedNode(vertexIdx); 
		}
	}
}
*/

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
void calculateGlobalStiffnessMatrix(TriList &K_triplets, ElementBuilder &elementBuilder, Mesh const &mesh, MatrixXd const &DOFTranslationMap) {
	MatrixXd TT; 
	MatrixXd TTi; // we don't use it.

	igl::triangle_triangle_adjacency(mesh.F, TT, TTi);

	for (int faceIdx = 0; faceIdx < mesh.F.rows(); faceIdx++) {
		IntList verticesIndices(6, ABSENT_VERTEX);
		Element currElement = createFaceElement(mesh, verticesIndices, faceIdx);
		//setFixedEdges(currElement, mesh.F.row(faceIdx), DOFTranslationMap); // TODO: need to set according to new input file
		setNbrsEnvelope(mesh, currElement, verticesIndices, faceIdx, TT);
		//if (!elementHasDOF(verticesIndices, DOFTranslationMap)) continue; // skip // FIXME

		elementBuilder.calculateStiffnessMatrix(currElement);
		addKeToK(K_triplets, currElement, DOFTranslationMap, verticesIndices);
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
* Ux = Map[v][0], Uy = Map[v][1], Uz = Map[v][2]
*/
void createDOFTranslationMap(Mesh const &mesh, MatrixXd &DOFTranslationMap) {
	DOFTranslationMap = MatrixXd::Constant(mesh.V.rows(), 3, FIXED_NODE);

	int j = 0, offset = 0;
	for (int i = 0; i < mesh.fixedDOF.size(); i++) {
		int vertexIdx = mesh.fixedDOF[i].first - 1;
		for (; j < vertexIdx; j++) {
			DOFTranslationMap.row(j) = Vector3d(offset, offset+1, offset+2);
			offset += 3;
		}
		for (int k = 0; k < 3; k++) {
			if (!mesh.fixedDOF[i].second(k)) // TODO: should be 1 or 0? what make smore sense? maybe reanme as freeDOF.
				DOFTranslationMap(j, k) = offset++;
		}
		j++;
	}
	for (; j < mesh.V.rows(); j++) {
		DOFTranslationMap.row(j) = Vector3d(offset, offset + 1, offset + 2);
		offset += 3;
	}
}

void preproccessForSolver(Mesh const &mesh, SparseMat &K, TriList &K_triplets, vector3dList const &nodalForces, MatrixXd &forces, MatrixXd &DOFTranslationMap) {
	int numOfDOF = 0; // TODO

	for (int i = 0; i < DOFTranslationMap.rows(); i++) { // need to save this as class parameter
		for (int j = 0; j < 3; j++) {
			if (DOFTranslationMap(i, j) > numOfDOF)
				numOfDOF = DOFTranslationMap(i, j);
		}
	}
	numOfDOF++;
	//std::cout << "#dof: " << numOfDOF << std::endl;

	// prepare K sparse matrix
	K = SparseMat(numOfDOF, numOfDOF);
	//std::cout << "set from triplets" << std::endl;
	//for (TripletXd t : K_triplets) std::cout << t.row() << ", " << t.col() << ", " << t.value() << std::endl;
	K.setFromTriplets(K_triplets.begin(), K_triplets.end());

	// prepare forces matrix
	forces = MatrixXd::Zero(numOfDOF, 1); // TODO why not vector?
	/*
	for (int i = 0; i < nodalForces.size(); i++) {
		int indexOfDOF = DOFTranslationMap[nodalForces[i].first-1] * 3;
		forces.block(indexOfDOF, 0, 3, 1) = nodalForces[i].second;
	}
	*/
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

bool solveSparseEquation(Mesh const &mesh, TriList &K_triplets, vector3dList const &nodalForces, MatrixXd &displacements, MatrixXd &DOFTranslationMap) {
	SparseMat K;								// Global Stiffness Matrix.
	MatrixXd forces;							// Forces vector for solver
	Eigen::ConjugateGradient<SparseMat> solver; // solve K*U = F. 
	
	preproccessForSolver(mesh, K, K_triplets, nodalForces, forces, DOFTranslationMap);

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
void calcStressFromDisplacements(Mesh const &mesh, ElementBuilder &elementBuilder, FEMResults &results, MatrixXd const &DOFTranslationMap) { 
	Eigen::Matrix<double, 18, 1> elementDisplacements;
	MatrixXd TT;
	MatrixXd TTi; // we don't use it.

	igl::triangle_triangle_adjacency(mesh.F, TT, TTi);
	
	results.vonMisesStress = VectorXd::Zero(mesh.F.rows());

	for (int faceIdx = 0; faceIdx < mesh.F.rows(); faceIdx++) {
		IntList verticesIndices(6, ABSENT_VERTEX);
		Element currElement = createFaceElement(mesh, verticesIndices, faceIdx);
		setNbrsEnvelope(mesh, currElement, verticesIndices, faceIdx, TT);
		//if (!elementHasDOF(verticesIndices, DOFTranslationMap)) continue; // skip // FIXME
		
		elementDisplacements = Eigen::MatrixXd::Zero(18, 1);
		for (int i = 0; i < 6; i++) {
			if (verticesIndices[i] == ABSENT_VERTEX) continue;
			for (int j = 0; j < 3; j++) {
				int dof = DOFTranslationMap(verticesIndices[i], j);
				if (dof != FIXED_NODE) {
					elementDisplacements(dof, 0) = results.displacements(dof, 0);
				}
			}
			/*
			int dof = DOFTranslationMap[verticesIndices[i]];
			if (dof != FIXED_NODE) {
				elementDisplacements.block(i * 3, 0, 3, 1) = results.displacements.block(dof * 3, 0, 3, 1);
			}
			*/
		}

		elementBuilder.calculateVonMisesStress(currElement, elementDisplacements);
		results.vonMisesStress(faceIdx) = currElement.vonMisesStress;
	}
}

void getDisplacedMesh(Mesh const &mesh, FEMResults &results, MatrixXd const &DOFTranslationMap) {
	results.displacedVertices = mesh.V;

	for (int vertexIdx = 0; vertexIdx < mesh.V.rows(); vertexIdx++) {
		for (int j = 0; j < 3; j++) {
			int dof = DOFTranslationMap(vertexIdx, j);
			if (dof != ABSENT_VERTEX) {
				results.displacedVertices(vertexIdx, j) += results.displacements(dof, 0);
			}
		}
		/*
		int dof = vertexToDOFTranslationMap[vertexIdx];
		if (dof != ABSENT_VERTEX) {
			results.displacedVertices.row(vertexIdx) += results.displacements.block(dof*3, 0, 3, 1).transpose();
		}
		*/
	}
}

void printSummary(Mesh const &mesh, FEMResults &results) { 
	std::cout << "Num of Vertices: " << mesh.V.rows() << std::endl;
	std::cout << "Num of Triangles: " << mesh.F.rows() << std::endl;
	std::cout << "Displacements: " << std::endl << results.displacements << std::endl;
	std::cout << "Displaced Vertices:" << std::endl << results.displacedVertices << std::endl;
	std::cout << "Von-Mises Stress:" << std::endl << results.vonMisesStress << std::endl;
}

void init(Mesh const &mesh, FEMData const &data, ElementBuilder &elementBuilder, MatrixXd &DOFTranslationMap) {
	Eigen::Matrix<double, 3, 3> De;									// Elastic Plastic Behavior Matrix
	Eigen::Matrix<double, 3, 3> M;									// Hill's Plastic Strain Matrix

	calculateElasticPlasticMatrix(De, data);
	calculateHillPlasticStrainMatrix(M, data);
	elementBuilder = ElementBuilder(De, M, data.shellProps.thickness);
	createDOFTranslationMap(mesh, DOFTranslationMap);
	std::cout << "Map:" << std::endl << DOFTranslationMap << std::endl; // FIXME  remove..
}

// TODO : sould run remove_duplicates before?
void Perform_FEM(Mesh const &mesh, vector3dList const &nodalForces, FEMData const &data, FEMResults &results) { // TODO should not be void
	MatrixXd DOFTranslationMap;					// [|V|x3] matrix of mapping of vertex Ux,Uy,Uz to dof #
	TriList K_triplets;							// List used to init the sparse global stiffness matrix
	ElementBuilder elementBuilder;
	
	init(mesh, data, elementBuilder, DOFTranslationMap);
	std::cout << "Beginnig Assembly Of Stiffness Matrix." << std::endl;
	calculateGlobalStiffnessMatrix(K_triplets, elementBuilder, mesh, DOFTranslationMap);

	std::cout << "Beginnig Solve" << std::endl;
	if (solveSparseEquation(mesh, K_triplets, nodalForces, results.displacements, DOFTranslationMap)) {
		return; // fail
	}
	
	calcStressFromDisplacements(mesh, elementBuilder, results, DOFTranslationMap);
	getDisplacedMesh(mesh, results, DOFTranslationMap);
	printSummary(mesh, results);

	return;
}

