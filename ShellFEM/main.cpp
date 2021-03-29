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
	else {
		startProgramGUI();
	}
}