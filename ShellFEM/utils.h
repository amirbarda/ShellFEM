#pragma once
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <fstream>
#include <iostream>
#include <igl\readSTL.h>
#include <igl\edge_topology.h>

typedef std::vector<int> IntList;
typedef std::vector<bool> BoolList;
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
	MatrixXd V;
	MatrixXi F;
	MatrixXi EV; //Ex2 matrix storing the edge description as pair of indices to
	MatrixXi FE; //Fx3 matrix storing the Triangle - Edge relation
	vector3dList fixedDOF;
	BoolList isEdgeClamped; // list of size |E|, for each edge true iff edge is clamped 
	Mesh(MatrixXd &V, MatrixXi &F, MatrixXi &EV, MatrixXi &FE, vector3dList &fixedDOF, BoolList &isEdgeClamped) :
		V(V), F(F), EV(EV), FE(FE), fixedDOF(fixedDOF), isEdgeClamped(isEdgeClamped){};
}; 

struct JobProperties {
	std::string name, outDir, objPath, forcesPath, fixedPath, clampedPath;
	bool startViewer;
	JobProperties() {};
	JobProperties(std::string name, std::string outDir, std::string objPath, std::string forcesPath, std::string fixedPath, std::string clampedPath, bool startViewer) :
		name(name), outDir(outDir), objPath(objPath), forcesPath(forcesPath), fixedPath(fixedPath), clampedPath(clampedPath), startViewer(startViewer) {};
};

std::vector<std::string> split_string_by_space(std::string s);
vector3dList vector3d_from_txt(std::string path);
BoolList clamped_from_txt(std::string path, int edgeCount);
MatrixXd displacements_from_txt(std::string path, int nodeNum);
VectorXd vonmises_from_txt(std::string path, int nodeNum);
void saveOBJ(MatrixXd &V, MatrixXi &F, std::string filepath);