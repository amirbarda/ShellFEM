#include <igl\triangle_triangle_adjacency.h>

#include "neighborhood.h"

int calcNbrOppositeVrtxIndx(MatrixXi &F, int nbrFaceIdx, int faceIdx, int edgeIdx) {
	int sumOfVertices = F.row(nbrFaceIdx).sum();
	int sumOfSharedVertices = F(faceIdx, edgeIdx) + F(faceIdx, (edgeIdx + 1) % 3);
	return sumOfVertices - sumOfSharedVertices;
}

/**	1. TT  : #F by #3 adjacent matrix. Description:
	   TT(i,j) = id of the neighbour triangle, that shares the j edge of triangle i .
	For a triangle, the 1st edge is [0,1] , the 2nd edge is [1,2] , the 3rd edge [2,3]. */
void calculateFaceNeighborhoodMatrix(MatrixXi &F, Eigen::MatrixXi &FNB) {
	MatrixXd TT, TTi;
	igl::triangle_triangle_adjacency(F, TT, TTi);

	FNB = MatrixXi::Constant(F.rows(), 6, ABSENT_VERTEX);
	FNB.block(0, 0, F.rows(), 3) = F;

	for (int faceIdx = 0; faceIdx < F.rows(); faceIdx++) {
		for (int edgeIdx = 0; edgeIdx < 3; edgeIdx++) {
			int nbrFaceIdx = TT(faceIdx, edgeIdx);
			if (nbrFaceIdx == ABSENT_VERTEX) continue; //case: no neighbour face sharing edge edgeIdx of current face.
			FNB(faceIdx, 3 + edgeIdx) = calcNbrOppositeVrtxIndx(F, nbrFaceIdx, faceIdx, edgeIdx);
		}
	}
}