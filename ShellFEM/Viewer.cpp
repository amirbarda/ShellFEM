#include "Viewer.h"
#include <igl\jet.h>
//#include <igl\colormap.h>

void Viewer::startView(Eigen::MatrixXd const &V, Eigen::MatrixXi const &F, VectorXd &vonMisesStress) {
	MatrixXd C;

	viewer.core().lighting_factor = 0;
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
	viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
	igl::jet(vonMisesStress, true, C);
	//igl::colormap(igl::COLOR_MAP_TYPE_PLASMA, vonMisesStress, true, C);
	viewer.data().set_colors(C);
	viewer.launch();
}