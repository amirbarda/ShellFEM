#pragma once
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <fstream>
#include <iostream>
#include <igl\readSTL.h>

typedef std::vector<int> IntList;
typedef std::pair<int, Eigen::Vector3d> Force;
typedef std::vector<Force> ForcesList;
typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::MatrixXi MatrixXi;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::VectorXi VectorXi;
typedef Eigen::SparseMatrix<double> SparseMat;
typedef Eigen::Triplet<double> TripletXd;
typedef std::vector<TripletXd> TriList;

struct Mesh {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	std::vector<int> fixedNode;
	Mesh(Eigen::MatrixXd const &V_, Eigen::MatrixXi const &F_, std::vector<int> const &fixedNode_) : V(V_), F(F_), fixedNode(fixedNode_){};
}; 

std::vector<std::string> split_string_by_space(std::string s);
std::vector<Force> nodal_forces_from_txt(std::string path);
std::vector<int> fixed_nodes_from_txt(std::string path);
Eigen::MatrixXd displacements_from_txt(std::string path, int nodeNum);
Eigen::VectorXd vonmises_from_txt(std::string path, int nodeNum);
void saveOBJ(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::string filepath);