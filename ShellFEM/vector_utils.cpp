#include "vector_utils.h"

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