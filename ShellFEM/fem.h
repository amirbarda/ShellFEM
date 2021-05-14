#pragma once
#include "utils.h"
#include "s3element.h"

struct FEMResults {
	MatrixXd displacements;
	VectorXd faceStress;		// vonMises yield criterion
	VectorXd vertexStress;		// vonMises yield criterion
	MatrixXd displacedVertices;
};

bool performFEM(Mesh &mesh, vector3dList const &nodalForces, SimulationProperties &simProps, FEMResults &results);
