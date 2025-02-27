#include <igl\edge_topology.h>
#include <igl\remove_duplicate_vertices.h>
#include <igl\readSTL.h>
#include <igl\readOBJ.h>

#include "job.h"
#include "gui.h"
#include "utils.h"
#include "neighborhood.h"

void runFEMJob(JobProperties &jobProps, SimulationProperties &simProps, bool createLogFile = false) {
	Eigen::MatrixXd V, TC, N, SV;
	Eigen::MatrixXi F, FTC, FN, EV, FE, EF, SVI, SVJ, SF, FNB;
	FILE *logFile = NULL;

	std::string saveObjPath = jobProps.outDir + "\\" + jobProps.name + ".obj";
	std::string logPath = jobProps.outDir + "\\" + jobProps.name + ".log";
	std::string saveDisplacementsPath = jobProps.outDir + "\\displacements.txt";
	std::string saveStressesPath = jobProps.outDir + "\\stresses.txt";

	if (createLogFile) {
		logFile = freopen(logPath.c_str(), "w", stdout); // setting stdout
	}

	if (jobProps.format == STL_T) {
		FILE *stlFile = fopen(jobProps.objPath.c_str(), "rb");
		igl::readSTL(stlFile, V, F, N);
		fclose(stlFile);
		igl::remove_duplicate_vertices(V, F, 1e-7, SV, SVI, SVJ, SF);
		std::cout << "SVI:" << std::endl << SVI << std::endl;
		std::cout << "SVJ:" << std::endl << SVJ << std::endl;
	}
	else {
		igl::readOBJ(jobProps.objPath, V, TC, N, F, FTC, FN);
		SV = V;
		SF = F;
	}

	SV *= jobProps.scale;
	std::cout << "SV: " << std::endl << SV << std::endl;
	std::cout << "SF: " << std::endl << SF << std::endl;

	calculateFaceNeighborhoodMatrix(SF, FNB);

	igl::edge_topology(SV, SF, EV, FE, EF);
	std::cout << "EV:" << std::endl;
	for (int i=0; i<EV.rows();i++) std::cout << i << ": " << EV.row(i) << std::endl;
	std::cout << DASH << std::endl;

	auto nodalForces = vector3d_from_txt(jobProps.forcesPath);
	auto freeDOF = vector3d_from_txt(jobProps.fixedPath);
	auto isEdgeClamped = clamped_from_txt(jobProps.clampedPath, EV);

	std::sort(nodalForces.begin(), nodalForces.end(), [](auto &left, auto &right) {
		return left.first < right.first;
	});
	std::sort(freeDOF.begin(), freeDOF.end(), [](auto &left, auto &right) {
		return left.first < right.first;
	});

	FEMResults results;
	Mesh mesh(SV, FNB, EV, FE, freeDOF, isEdgeClamped);
	performFEM(mesh, nodalForces, simProps, results);
	saveOBJ(results.displacedVertices, SF, saveObjPath);
	saveMatrix(results.displacements, saveDisplacementsPath);
	saveMatrix(results.vertexStress, saveStressesPath);

	if (jobProps.startViewer) {
		Viewer viewer;
		viewer.startView(results.displacedVertices, SF, results.faceStress);
	}
	if (createLogFile) {
		fclose(logFile);
	}
}