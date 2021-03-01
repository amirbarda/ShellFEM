#pragma once
#include "FEM.h"
#include "Viewer.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <iostream>


void saveOBJ(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::string filepath) {
	std::ofstream file;
	file.open(filepath);
	for (int i = 0; i < V.rows(); i++) {
		file << std::fixed << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
	}
	file << std::endl;
	for (int i = 0; i < F.rows(); i++) {
		file << "f " << F(i, 0) + 1 << " " << F(i, 1) + 1 << " " << F(i, 2) + 1 << std::endl;
	}
	file.close();
}

int main() {

	Eigen::MatrixXd V;
	Eigen::MatrixXd TC;
	Eigen::MatrixXd N;
	Eigen::MatrixXi F;
	Eigen::MatrixXi FTC;
	Eigen::MatrixXi FN;
	std::string objPath = "..\\..\\tests\\test2\\pyramid.obj";
	std::string outputObjPath = "..\\..\\tests\\test2\\pyramid_output.obj";
	std::string stdoutPath = "..\\..\\tests\\test2\\test2.log";
	std::string nodalForcesPath = "..\\..\\tests\\test2\\load_nodes.txt";
	std::string fixedNodesPath = "..\\..\\tests\\test2\\fixed_nodes.txt";
	FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

	igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
	auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
	auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
	FEMData data; //gets default data (for now)

	FEMResults result;
	Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
	saveOBJ(result.displacedVertices, F, outputObjPath);

	Viewer viewer;
	viewer.startView(result.displacedVertices, F, result.vonMisesStress);
	fclose(file);

}