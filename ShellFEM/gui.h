#pragma once
#include <igl/opengl/glfw/Viewer.h> 

#include "utils.h"

class Viewer {
private:
	igl::opengl::glfw::Viewer viewer;
public:
	Viewer() {};
	void startView(Eigen::MatrixXd const &V, Eigen::MatrixXi const &F, VectorXd &vonMisesStress);
};

void startProgramGUI();