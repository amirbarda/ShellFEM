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
	std::string objPath = "..\\..\\tests\\test4\\plate_A11.obj";
	std::string outputObjPath = "..\\..\\tests\\test4\\plate_A11_output.obj";
	std::string stdoutPath = "..\\..\\tests\\test4\\test4.log";
	std::string nodalForcesPath = "..\\..\\tests\\test4\\load_nodes_A11.txt";
	std::string fixedNodesPath = "..\\..\\tests\\test4\\fixed_nodes_A11.txt";
	FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

	igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
	//V = V;
	auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
	auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
	FEMData data; //gets default data (for now)
	data.shellProps.thickness = 1e-3;
	data.matProps.E = 200e9;

	FEMResults result;
	Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
	std::cout << "Uzc: "<<(result.displacedVertices.row(2) - V.row(2)).norm() << std::endl;
	saveOBJ(result.displacedVertices, F, outputObjPath);

	Viewer viewer;
	viewer.startView(result.displacedVertices, F, result.vonMisesStress);
	fclose(file);

}