#pragma once
#include "utils.h"
#include "s3element.h"

struct FEMResults {
	MatrixXd displacements;
	VectorXd vonMisesStress;
	MatrixXd displacedVertices;
};

bool performFEM(Mesh &mesh, vector3dList const &nodalForces, SimulationProperties &simProps, FEMResults &results);
