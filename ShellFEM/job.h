#pragma once

#include "utils.h"
#include "fem.h"
#include "s3element.h"

typedef enum { OBJ_T, STL_T } FileFormat;

struct JobProperties {
	std::string name, outDir, objPath, forcesPath, fixedPath, clampedPath;
	FileFormat format;
	bool startViewer;
	JobProperties() {};
	JobProperties(std::string name, std::string outDir, std::string objPath, std::string forcesPath, std::string fixedPath, std::string clampedPath, bool startViewer, FileFormat format) :
		name(name), outDir(outDir), objPath(objPath), forcesPath(forcesPath), fixedPath(fixedPath), clampedPath(clampedPath), startViewer(startViewer), format(format) {};
};

void runFEMJob(JobProperties &jobProps, SimulationProperties &simProps, bool createLogFile);