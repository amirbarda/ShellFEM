#include <igl\jet.h>

#include "gui.h"
#include "fem.h"
#include "utils.h"

void Viewer::startView(Eigen::MatrixXd const &V, Eigen::MatrixXi const &F, VectorXd &vonMisesStress) {
	MatrixXd C;

	viewer.core().lighting_factor = 0.4;
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
	viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
	igl::jet(vonMisesStress, true, C);
	viewer.data().set_colors(C);
	viewer.launch();
}

void startProgramGUI() {
	Eigen::MatrixXd V;
	Eigen::MatrixXd TC;
	Eigen::MatrixXd N;
	Eigen::MatrixXi F;
	Eigen::MatrixXi FTC;
	Eigen::MatrixXi FN;

	FEMResults result;
	SimulationProperties simProps;
	std::string saveObjPath = simProps.outDir + "/" + simProps.name + ".obj";
	FEMData data(simProps.E, simProps.ni, simProps.thickness);

	igl::readOBJ(simProps.objPath, V, TC, N, F, FTC, FN);
	//auto nodalForces = nodal_forces_from_txt(simProps.forcesPath);
	//auto fixedNodes = fixed_nodes_from_txt(simProps.fixedPath);

	//Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
	saveOBJ(result.displacedVertices, F, saveObjPath);

	// TODO : add imgui and call startView
}