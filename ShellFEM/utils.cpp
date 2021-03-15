#include "utils.h"

std::vector<std::string> split_string_by_space(std::string s) { 
	std::vector<std::string> result;
	std::istringstream iss(s);
	for (std::string s; iss >> s; ) {
		result.push_back(s);
	}
	return result;
}

vector3dList vector3d_from_txt(std::string path) {
	vector3dList iVecList;
	std::ifstream inFile(path);
	std::string line;
	while (std::getline(inFile, line)) {
		indexedVector3d iVec;
		auto tokens = split_string_by_space(line);
		iVec.first = std::stod(tokens[0]);
		iVec.second = Eigen::Vector3d(std::stod(tokens[1]), std::stod(tokens[2]), std::stod(tokens[3]));
		iVecList.push_back(iVec);
	}
	inFile.close();
	return iVecList;
}

IntList clamped_from_txt(std::string path) {
	IntList clampedEdges;
	std::ifstream inFile(path);
	std::string line;
	while (std::getline(inFile, line)) {
		auto tokens = split_string_by_space(line);
		clampedEdges.push_back(std::stod(tokens[0]));
	}
	inFile.close();
	return clampedEdges;
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

void saveOBJ(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::string filepath) {
	std::ofstream file;
	file.open(filepath);
	for (int i = 0; i < V.rows(); i++) {
		file << std::fixed << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
	}
	file << std::endl;
	for (int i = 0; i < F.rows(); i++) {
		file << "f " << F(i, 0) + 1 << " " << F(i, 1) + 1 << " " << F(i, 2) + 1 << std::endl;
	}
	file.close();
}