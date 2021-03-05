#pragma once

#include <string>

typedef struct SimulationProperties {
	std::string name, outDir, objPath, forcesPath, fixedPath;
	bool startViewer;
	double E; //Young's Modulus
	double ni; //Possion's Ratio
	double thickness; 
	SimulationProperties(std::string name_, std::string outDir_, std::string objPath_, std::string forcesPath_, std::string fixedPath_, bool startViewer_,
			double E_, double ni_, double thickness_) :
		name(name_), outDir(outDir_), objPath(objPath_), forcesPath(forcesPath_), fixedPath(fixedPath_), startViewer(startViewer_),
			E(E_), ni(ni_), thickness(thickness_){};
} SimulationProperties;

SimulationProperties parse_arguments(int argc, char** argv);