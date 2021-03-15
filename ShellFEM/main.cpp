#include "fem.h"
#include "options.h"
#include "gui.h"
#include "utils.h"

void runInBatchMode(int argc, char** argv) {
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

	igl::readOBJ(simProps.objPath, V, TC, N, F, FTC, FN);
	//auto nodalForces = nodal_forces_from_txt(simProps.forcesPath);
	//auto fixedNodes = fixed_nodes_from_txt(simProps.fixedPath);

	//Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
	saveOBJ(result.displacedVertices, F, saveObjPath);

	if (simProps.startViewer) {
		Viewer viewer;
		viewer.startView(result.displacedVertices, F, result.vonMisesStress);
	}
}

int main(int argc, char** argv) {
	if (argc > 1) {
		runInBatchMode(argc, argv);
	}
	else {
		startProgramGUI();
	}
}