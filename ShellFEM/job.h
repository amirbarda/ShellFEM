#pragma once

#include "utils.h"
#include "fem.h"
#include "s3element.h"

typedef enum { OBJ_T, STL_T } FileFormat;

struct JobProperties { // TODO : add print representation..
	std::string name, outDir, objPath, forcesPath, fixedPath, clampedPath;
	FileFormat format;
	bool startViewer;
	JobProperties() {};
	JobProperties(const std::string& name, const std::string& outDir, const std::string& objPath, const std::string& forcesPath, const std::string& fixedPath,
		const std::string& clampedPath, bool startViewer, FileFormat format) :
		name(name), outDir(outDir), objPath(objPath), forcesPath(forcesPath), fixedPath(fixedPath), clampedPath(clampedPath), startViewer(startViewer), format(format) {};
};

void runFEMJob(JobProperties &jobProps, SimulationProperties &simProps, bool createLogFile);