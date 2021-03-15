#pragma once
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <fstream>
#include <iostream>
#include <igl\readSTL.h>

typedef std::vector<int> IntList;
typedef std::pair<int, Eigen::Vector3d> indexedVector3d;
typedef std::vector<indexedVector3d> vector3dList;
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
	Eigen::MatrixXi EV; //Ex2 matrix storing the edge description as pair of indices to
	Eigen::MatrixXi FE; //Fx3 matrix storing the Triangle - Edge relation
	vector3dList fixedDOF;
	IntList clampedEdges;
	Mesh(MatrixXd const &V_, MatrixXi const &F_, vector3dList const &fixedDOF_, IntList const &clampedEdges_) :
		V(V_), F(F_), fixedDOF(fixedDOF_), clampedEdges(clampedEdges){};
}; 

std::vector<std::string> split_string_by_space(std::string s);
vector3dList vector3d_from_txt(std::string path);
IntList clamped_from_txt(std::string path);
MatrixXd displacements_from_txt(std::string path, int nodeNum);
VectorXd vonmises_from_txt(std::string path, int nodeNum);
void saveOBJ(MatrixXd &V, MatrixXi &F, std::string filepath);