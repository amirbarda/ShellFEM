#include "utils.h"

std::vector<std::string> split_string_by_space(std::string s) { 
	std::vector<std::string> result;
	std::istringstream iss(s);
	for (std::string s; iss >> s; ) {
		result.push_back(s);
	}
	return result;
}

std::vector<Force> nodal_forces_from_txt(std::string path) {
	std::vector<Force> nodalForces;
	std::ifstream nodalForcesFile(path);
	if (nodalForcesFile.is_open()) {
		std::string line;
		while (std::getline(nodalForcesFile, line)) {
			Force f;
			auto tokens = split_string_by_space(line);
			f.first = std::stoi(tokens[0]);
			f.second = Eigen::Vector3d(std::stod(tokens[1]), std::stod(tokens[2]), std::stod(tokens[3]));
			nodalForces.push_back(f);
		}
	}
	nodalForcesFile.close();
	return nodalForces;
}

std::vector<int> fixed_nodes_from_txt(std::string path) {
	std::vector<int> fixedNodes;
	std::ifstream fixedNodesFile(path);
	if (fixedNodesFile.is_open()) {
		std::string line;
		while (std::getline(fixedNodesFile, line)) {
			fixedNodes.push_back(std::stod(line));
		}
	}
	fixedNodesFile.close();
	return fixedNodes;
}

Eigen::VectorXd vonmises_from_txt(std::string path, int nodeNum) {
	//get results
	Eigen::VectorXd stresses(nodeNum);
	std::ifstream results(path);
	std::string line;
	std::string delimiter = " ";
	while (getline(results, line)) {
		int node_id = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + delimiter.length());
		double misesStress = std::stod(line);
		stresses(node_id) = misesStress;
	}
	results.close();
	return stresses;
}

Eigen::MatrixXd displacements_from_txt(std::string path, int nodeNum) {
	std::string line;
	std::ifstream results_disp;
	std::string delimiter = " ";
	Eigen::MatrixXd displacements(nodeNum,3);
	results_disp.open(path);
	while (getline(results_disp, line)) {
		int node_id = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + delimiter.length());
		double x = std::stod(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + delimiter.length());
		double y = std::stod(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + delimiter.length());
		double z = std::stod(line.substr(0, line.find(delimiter)));
		displacements.row(node_id) << x, y, z;
	}
	results_disp.close();
	return displacements;
}