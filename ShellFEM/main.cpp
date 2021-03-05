#include "fem.h"
#include "options.h"
#include "gui.h"

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

int main(int argc, char** argv) {
	Eigen::MatrixXd V;
	Eigen::MatrixXd TC;
	Eigen::MatrixXd N;
	Eigen::MatrixXi F;
	Eigen::MatrixXi FTC;
	Eigen::MatrixXi FN;

	FEMResults result;
	SimulationProperties simProps = parse_arguments(argc, argv); 
	std::string saveObjPath = simProps.outDir + "/" + simProps.name + ".obj";
	FEMData data(simProps.E, simProps.ni, simProps.thickness);
	Viewer viewer;

	igl::readOBJ(simProps.objPath, V, TC, N, F, FTC, FN);
	auto nodalForces = nodal_forces_from_txt(simProps.forcesPath);
	auto fixedNodes = fixed_nodes_from_txt(simProps.fixedPath);

	Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
	saveOBJ(result.displacedVertices, F, saveObjPath);

	if (simProps.startViewer) {
		viewer.startView(result.displacedVertices, F, result.vonMisesStress);
	}
}