#pragma once

#include <Eigen/Core>

std::pair<double, Eigen::Vector3d> getAreaVector(Eigen::Vector3d const &v1, Eigen::Vector3d const &v2, Eigen::Vector3d const &v3);
double getCosOfAngle(Eigen::Vector3d const &v1, Eigen::Vector3d const &v2, Eigen::Vector3d const &v3);

Eigen::Matrix3d getRotationMatrix(Eigen::Vector3d const &oldAxis, Eigen::Vector3d const &newAxis, int axisToAlign);

/**
	newBase is a Matrix3d, where each coloumn is a unit vector. (ux, uy, uz)
	Output is localVertices[i] = Rotation*vertices[i]
*/
//void calcRotatedVertices(Eigen::Matrix3d &newBase, Eigen::Vector3d const vertices[], Eigen::Vector3d localVertices[]);