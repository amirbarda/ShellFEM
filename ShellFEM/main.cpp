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
 		std::string objPath = "C:\\Users\\Michael\\Documents\\GoogleDrive\\GoogleDrive\\VisualStudioC\\ShellFEM\\x64\\Release\\test\\mesh.obj";
		std::string outDir = "C:\\Users\\Michael\\Documents\\GoogleDrive\\GoogleDrive\\VisualStudioC\\ShellFEM\\x64\\Release\\test";
		std::string forcesPath = "C:\\Users\\Michael\\Documents\\GoogleDrive\\GoogleDrive\\VisualStudioC\\ShellFEM\\x64\\Release\\test\\load_nodes.txt";
		std::string fixedPath = "C:\\Users\\Michael\\Documents\\GoogleDrive\\GoogleDrive\\VisualStudioC\\ShellFEM\\x64\\Release\\test\\fixed_nodes.txt";
		std::string clampedPath = "C:\\Users\\Michael\\Documents\\GoogleDrive\\GoogleDrive\\VisualStudioC\\ShellFEM\\x64\\Release\\test\\clamped_edges.txt";
		JobProperties jobProps(name, outDir, objPath, forcesPath, fixedPath, clampedPath, true, OBJ_T, 1);
		SimulationProperties simProps;

		runFEMJob(jobProps, simProps, true);
	}
}
