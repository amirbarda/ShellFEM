#include "fem.h"
#include "options.h"
#include "job.h"
#include "gui.h"
#include "utils.h"

void runInBatchMode(int argc, char** argv) {
	JobProperties jobProps;
	SimulationProperties simProps;
	parse_arguments(argc, argv, simProps, jobProps);
	runFEMJob(jobProps, simProps, true); // FIXME : log file creation flag
}

int main(int argc, char** argv) {
	if (argc > 1) {
		runInBatchMode(argc, argv);
	}
	else { // FIXME
		//startProgramGUI();
		std::string name = "test";
 		std::string objPath = "test\\mesh.obj";
		std::string outDir = "test";
		std::string forcesPath = "test\\load_nodes.txt";
		std::string fixedPath = "test\\fixed_nodes.txt";
		std::string clampedPath = "test\\clamped_edges.txt";
		JobProperties jobProps(name, outDir, objPath, forcesPath, fixedPath, clampedPath, true, OBJ_T);
		SimulationProperties simProps;

		runFEMJob(jobProps, simProps, true);
	}
}
