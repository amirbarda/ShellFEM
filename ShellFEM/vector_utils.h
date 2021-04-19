#pragma once

#include <Eigen/Core>

std::pair<double, Eigen::Vector3d> getAreaVector(Eigen::Vector3d const &v1, Eigen::Vector3d const &v2, Eigen::Vector3d const &v3);
double getCosOfAngle(Eigen::Vector3d const &v1, Eigen::Vector3d const &v2, Eigen::Vector3d const &v3);
void calcRotatedVertices(Eigen::Matrix3d &newBase, Eigen::Vector3d const vertices[], Eigen::Vector3d localVertices[]);