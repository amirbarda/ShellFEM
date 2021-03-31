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
 		std::string objPath = "C:\\Users\\Michael\\Documents\\GoogleDrive\\GoogleDrive\\VisualStudioC\\ShellFEM\\tests\\test3\\plate.obj";
		std::string outDir = "C:\\Users\\Michael\\Documents\\GoogleDrive\\GoogleDrive\\VisualStudioC\\ShellFEM\\tests\\test3";
		std::string forcesPath = "C:\\Users\\Michael\\Documents\\GoogleDrive\\GoogleDrive\\VisualStudioC\\ShellFEM\\tests\\test3\\load_nodes.txt";
		std::string fixedPath = "C:\\Users\\Michael\\Documents\\GoogleDrive\\GoogleDrive\\VisualStudioC\\ShellFEM\\tests\\test3\\fixed_nodes.txt";
		std::string clampedPath = "C:\\Users\\Michael\\Documents\\GoogleDrive\\GoogleDrive\\VisualStudioC\\ShellFEM\\tests\\test3\\clamped_edges.txt";
		JobProperties jobProps(name, outDir, objPath, forcesPath, fixedPath, clampedPath, true, OBJ_T);
		SimulationProperties simProps;

		runFEMJob(jobProps, simProps, true);
	}
}
