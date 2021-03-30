#include <igl\jet.h>

#include "gui.h"
#include "utils.h"

void Viewer::startView(Eigen::MatrixXd const &V, Eigen::MatrixXi const &F, VectorXd &vonMisesStress) {
	MatrixXd C;

	viewer.core().lighting_factor = 0.4;
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
	viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
	igl::jet(vonMisesStress, true, C); // FIXME
	viewer.data().set_colors(C);
	viewer.launch();
}

void startProgramGUI() {
	std::cout << "Not implemented yet" << std::endl;
	// TODO : add imgui and call startView
}

