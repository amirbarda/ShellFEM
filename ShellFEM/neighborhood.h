#pragma once

#include "utils.h"

/**
	Returns matrix of [#face X 6], which includes indices of vertices of face, and indices of neighbors.
	vertex 3 is neighbor of edge (0,1).
	vertex 4 is neighbor of edge (1,2).
	vertex 5 is neighbor of edge (2,1).
*/
void calculateFaceNeighborhoodMatrix(MatrixXi &F, MatrixXi &FNB);