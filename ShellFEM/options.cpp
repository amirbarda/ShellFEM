#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>
#include <ctime>

#include "options.h"

std::string getDateAndTime() {
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, sizeof(buffer), "%d_%m_%Y_%H_%M_%S", timeinfo);
	std::string str(buffer);

	return str;
}


SimulationProperties parse_arguments(int argc, char** argv) {
	try {
		// Init
		TCLAP::CmdLine cmd("Command description message", ' ', "0.9");

		// Add commands
		TCLAP::ValueArg<std::string> nameArg("n", "name", "Name of simulation", false, getDateAndTime(), "string");
		cmd.add(nameArg);

		TCLAP::ValueArg<std::string> outputDirArg("o", "out_dir", "Path of output logs", false, ".", "string");
		cmd.add(outputDirArg);

		TCLAP::ValueArg<std::string> objPathArg("f", "obj_path", "Path of input obj", true, " ", "string");
		cmd.add(objPathArg);

		TCLAP::ValueArg<std::string> forcesPathArg("c", "forces_path", "Path of forces file", true, " ", "string");
		cmd.add(forcesPathArg);

		TCLAP::ValueArg<std::string> fixedPathArg("x", "fixed_path", "Path of fixed nodes file", true, " ", "string");
		cmd.add(fixedPathArg);

		TCLAP::ValueArg<double> youngModulusArg("e", "young_modulus", "Value of Young's Modulus", false, 70e9, "double"); 
		cmd.add(youngModulusArg);

		TCLAP::ValueArg<double> possionRatioArg("p", "possion_ratio", "Value of Possion's ratio", false, 0.3, "double");
		cmd.add(possionRatioArg);

		TCLAP::ValueArg<double> thicknessArg("t", "thickness", "Value of sheet thickness", false, 1.6e-3, "double");
		cmd.add(thicknessArg);

		TCLAP::SwitchArg viewerSwitch("g", "gui", "start viewer", cmd, false);

		// Parse the argv array.
		cmd.parse(argc, argv);

		std::cout << "nameArg: " << nameArg.getValue() << std::endl;
		std::cout << "outputDirArg: " << outputDirArg.getValue() << std::endl;
		std::cout << "objPathArg: " << objPathArg.getValue() << std::endl;
		std::cout << "forcesPathArg: " << forcesPathArg.getValue() << std::endl;
		std::cout << "fixedPathArg: " << fixedPathArg.getValue() << std::endl;
		std::cout << "youngModulusArg: " << youngModulusArg.getValue() << std::endl;
		std::cout << "possionRatioArg: " << possionRatioArg.getValue() << std::endl;
		std::cout << "thicknessArg: " << thicknessArg.getValue() << std::endl;
		std::cout << "viewerSwitch: " << viewerSwitch.getValue() << std::endl;

		// Get the value parsed by each arg. 
		SimulationProperties simProps(nameArg.getValue(), outputDirArg.getValue(), objPathArg.getValue(), forcesPathArg.getValue(), fixedPathArg.getValue(),
			viewerSwitch.getValue(), youngModulusArg.getValue(), possionRatioArg.getValue(), thicknessArg.getValue());
		
		return simProps;
	} catch (TCLAP::ArgException &e) { // catch exceptions
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		exit(1);
	}
}