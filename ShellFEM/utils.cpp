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
		if (tokens.size() > 1) iVec.second = Vector3d(std::stod(tokens[1]), std::stod(tokens[2]), std::stod(tokens[3]));
		else iVec.second = Vector3d(0, 0, 0);
		iVecList.push_back(iVec);
	}
	inFile.close();
	return iVecList;
}

void initMap(Eigen::MatrixXi &EV, std::vector<std::vector<int>> &mapEV) {
	for (int i = 0; i < EV.rows(); i++) {
		mapEV[EV(i, 0)][EV(i, 1)] = i;
		mapEV[EV(i, 1)][EV(i, 0)] = i;
	}
}

BoolList clamped_from_txt(std::string path, Eigen::MatrixXi &EV) {
	std::vector<std::vector<int>> mapEV(EV.rows(), std::vector<int>(EV.rows(), -1));
	BoolList isEdgeClamped(EV.rows(), false);
	std::ifstream inFile(path);
	std::string line;
	int edgeIdx;
	
	initMap(EV, mapEV);
	while (std::getline(inFile, line)) {
		auto tokens = split_string_by_space(line);
		if (tokens.size() > 1) {
			edgeIdx = mapEV[std::stoi(tokens[0])][std::stoi(tokens[1])];
		}
		else edgeIdx = std::stoi(tokens[0]);
		isEdgeClamped[edgeIdx] = true;
	}
	inFile.close();
	return isEdgeClamped;
}

VectorXd vonmises_from_txt(std::string path, int nodeNum) {
	//get results
	VectorXd stresses(nodeNum);
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

MatrixXd displacements_from_txt(std::string path, int nodeNum) {
	std::string line;
	std::ifstream results_disp;
	std::string delimiter = " ";
	MatrixXd displacements(nodeNum,3);
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

void saveOBJ(MatrixXd &V, MatrixXi &F, std::string filepath) {
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
