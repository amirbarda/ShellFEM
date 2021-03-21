#include "fem.h"
#include "options.h"
#include "gui.h"
#include "utils.h"

void runInBatchMode(int argc, char** argv) {
	Eigen::MatrixXd V, TC, N;
	Eigen::MatrixXi F, FTC, FN, EV, FE, EF;

	FEMResults results;
	JobProperties jobProps;
	SimulationProperties simProps;
	
	parse_arguments(argc, argv, simProps, jobProps);
	std::string saveObjPath = jobProps.outDir + "/" + jobProps.name + ".obj";

	igl::readOBJ(jobProps.objPath, V, TC, N, F, FTC, FN);
	igl::edge_topology(V, F, EV, FE, EF);
	auto nodalForces = vector3d_from_txt(jobProps.forcesPath);
	auto fixedNodes = vector3d_from_txt(jobProps.fixedPath);
	auto isEdgeClamped = clamped_from_txt(jobProps.clampedPath, EV.rows());

	Mesh mesh(V, F, EV, FE, fixedNodes, isEdgeClamped);
	performFEM(mesh, nodalForces, simProps, results);
	saveOBJ(results.displacedVertices, F, saveObjPath);

	if (jobProps.startViewer) {
		Viewer viewer;
		viewer.startView(results.displacedVertices, F, results.vonMisesStress);
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