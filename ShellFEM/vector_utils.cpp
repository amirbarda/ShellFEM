#include <Eigen/Geometry>
#include <iostream>

#include "vector_utils.h"

// FIXME
#define X 0
#define Y 1
#define Z 2
#define DASH "--------------------------------------------------------------------------"

std::pair<double, Eigen::Vector3d> getAreaVector(Eigen::Vector3d const &mid, Eigen::Vector3d const &next, Eigen::Vector3d const &prev) {
	Eigen::Vector3d nextEdge = next - mid;
	Eigen::Vector3d prevEdge = prev - mid;
	Eigen::Vector3d cross = nextEdge.cross(prevEdge);
	double norm = cross.norm();
	return std::make_pair(0.5*norm, cross / norm);
}

double getCosOfAngle(Eigen::Vector3d const &mid, Eigen::Vector3d const &v1, Eigen::Vector3d const &v2) {
	Eigen::Vector3d edge1 = v1 - mid;
	Eigen::Vector3d edge2 = v2 - mid;
	double length1 = edge1.norm();
	double length2 = edge2.norm();
	double dot = edge1.dot(edge2);
	return dot / (length1*length2);
}
/*
Eigen::Matrix3d getRodrigezRotationMatrix(Eigen::Vector3d const &oldUnitVec, Eigen::Vector3d const &newUnitVec) {
	Eigen::Matrix3d rotationMat, crossProductMat;
	Eigen::Vector3d v;
	double c, s;

	v = oldUnitVec.cross(newUnitVec);
	s = v.norm();
	c = oldUnitVec.dot(newUnitVec);

	crossProductMat <<  0, -v[Z], v[Y],
						v[Z], 0, -v[X],
						-v[Y], v[X], 0;

	rotationMat << Eigen::MatrixXd::Identity(3, 3);
	rotationMat += s*crossProductMat + (1-c)*(crossProductMat * crossProductMat);

	return rotationMat;
}


void calcRotatedVertices(Eigen::Vector3d newAxes[3], Eigen::Vector3d const vertices[], Eigen::Vector3d localVertices[]) {
	Eigen::Matrix3d rotationMat;

	rotationMat = getRodrigezRotationMatrix(Eigen::Vector3d(0,0,1), newAxes[Z]);	// rotate Z axis to align to new Z
	rotationMat *= getRodrigezRotationMatrix(rotationMat.col(Y), newAxes[Y]);		// perform both rotations

	std::cout << "rotation matrix: " << std::endl << rotationMat << std::endl;
	for (int i = 0; i < 3; i++) {
		localVertices[i] = rotationMat * vertices[i];
		std::cout << "rotated: " << vertices[i] << std::endl;
	}
	std::cout << DASH << std::endl;
}

*/