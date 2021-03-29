#include <igl\edge_topology.h>
#include <igl\remove_duplicate_vertices.h>
#include <igl\readSTL.h>
#include <igl\readOBJ.h>

#include "job.h"
#include "gui.h"
#include "utils.h"

void runFEMJob(JobProperties &jobProps, SimulationProperties &simProps, bool createLogFile = false) {
	Eigen::MatrixXd V, TC, N, SV;
	Eigen::MatrixXi F, FTC, FN, EV, FE, EF, SVI, SVJ, SF;
	FILE *logFile = NULL;

	std::string saveObjPath = jobProps.outDir + "/" + jobProps.name + ".obj";
	std::string logPath = jobProps.outDir + "/" + jobProps.name + ".log";
	std::string saveDisplacementsPath = jobProps.outDir + "/displacements.txt";
	std::string saveStressesPath = jobProps.outDir + "/stresses.txt";

	if (createLogFile) {
		logFile = freopen(logPath.c_str(), "w", stdout); // setting stdout
	}

	if (jobProps.format == STL_T) {
		FILE *stlFile = fopen(jobProps.objPath.c_str(), "rb");
		igl::readSTL(stlFile, V, F, N);
		fclose(stlFile);
	}

	igl::remove_duplicate_vertices(V, F, 1e-7, SV, SVI, SVJ, SF);
	igl::edge_topology(SV, SF, EV, FE, EF);

	auto nodalForces = vector3d_from_txt(jobProps.forcesPath);
	auto freeDOF = vector3d_from_txt(jobProps.fixedPath);
	auto isEdgeClamped = clamped_from_txt(jobProps.clampedPath, EV.rows());

	std::sort(nodalForces.begin(), nodalForces.end(), [](auto &left, auto &right) {
		return left.first < right.first;
	});
	std::sort(freeDOF.begin(), freeDOF.end(), [](auto &left, auto &right) {
		return left.first < right.first;
	});

	FEMResults results;
	Mesh mesh(SV, SF, EV, FE, freeDOF, isEdgeClamped);
	performFEM(mesh, nodalForces, simProps, results);
	saveOBJ(results.displacedVertices, SF, saveObjPath);
	saveMatrix(results.displacements, saveDisplacementsPath);
	saveMatrix(results.vonMisesStress, saveStressesPath);

	if (jobProps.startViewer) {
		Viewer viewer;
		viewer.startView(results.displacedVertices, SF, results.vonMisesStress);
	}
	if (createLogFile) {
		fclose(logFile);
	}
}