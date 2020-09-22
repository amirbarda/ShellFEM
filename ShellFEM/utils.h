#pragma once
#include <iostream>
#include <Eigen\core>
#include <vector>
#include <fstream>
#include <igl\readSTL.h>

typedef std::pair<int, Eigen::Vector3d> Force;

std::vector<std::string> split_string_by_space(std::string s);
std::vector<Force> nodal_forces_from_txt(std::string path);
std::vector<int> fixed_nodes_from_txt(std::string path);
Eigen::MatrixXd displacements_from_txt(std::string path, int nodeNum);
Eigen::VectorXd vonmises_from_txt(std::string path, int nodeNum);