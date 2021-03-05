#pragma once
#include <vector>


class Logger {

private:
	bool logElementParameters = false;
	std::vector<int> logElementParametersIdsArr; 	// empty (default) vector logs all elements

public:
	void log_element_parameters(std::vector<int> elementIds = std::vector<int>());
};


extern Logger logger;