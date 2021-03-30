#include "stdafx.h"

#include <utility>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <fstream>

#include "CppUnitTest.h"

#include "utils.cpp"
#include "vector_utils.cpp"
#define NOMINMAX
#include "glad.c"
#include "s3element.cpp"
#include "fem.cpp"
#include "logger.cpp"
#include "gui.cpp"
#include "job.cpp" 

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace ShellFEMUnitTests
{
	TEST_CLASS(FEM)
	{
	public:
		/*
		TEST_METHOD(Perform_FEM_test3)
		{
			std::string objPath = "..\\..\\tests\\test3\\plate.obj";
			std::string outDir = "..\\..\\tests\\test3\\";
			std::string forcesPath = "..\\..\\tests\\test3\\load_nodes.txt";
			std::string fixedPath = "..\\..\\tests\\test3\\fixed_nodes.txt";
			std::string clampedPath = "..\\..\\tests\\test3\\clamped_edges.txt";
			std::string fixedPath = "..\\..\\tests\\test3\\fixed_nodes.txt";
			JobProperties jobProps("test", outDir, objPath, forcesPath, fixedPath, clampedPath, true, OBJ_T);
			SimulationProperties simProps;

			runFEMJob(jobProps, simProps, true);
		}

		TEST_METHOD(Perform_FEM_test12)	{
			std::string objPath = "..\\..\\tests\\test12\\shell.obj";
			std::string outDir = "..\\..\\tests\\test12\\";
			std::string forcesPath = "..\\..\\tests\\test12\\load_nodes.txt";
			std::string fixedPath = "..\\..\\tests\\test12\\fixed_nodes.txt";
			std::string clampedPath = "..\\..\\tests\\test12\\clamped_edges.txt";
			std::string fixedPath = "..\\..\\tests\\test12\\fixed_nodes.txt";
			JobProperties jobProps("test", outDir, objPath, forcesPath, fixedPath, clampedPath, true, OBJ_T);
			SimulationProperties simProps(6.825e7, DEFAULT_POSSION_CNST, 0.04);

			runFEMJob(jobProps, simProps, true);
		}

		TEST_METHOD(Perform_FEM_test13) {
			std::string objPath = "..\\..\\tests\\test13\\remeshed_surface.STL";
			std::string outDir = "..\\..\\tests\\test13\\";
			std::string forcesPath = "..\\..\\tests\\test13\\load_nodes.txt";
			std::string fixedPath = "..\\..\\tests\\test13\\fixed_nodes.txt";
			std::string clampedPath = "..\\..\\tests\\test13\\clamped_edges.txt";
			std::string fixedPath = "..\\..\\tests\\test13\\fixed_nodes.txt";
			JobProperties jobProps("test", outDir, objPath, forcesPath, fixedPath, clampedPath, true, STL_T);
			SimulationProperties simProps;

			runFEMJob(jobProps, simProps, true);
		}

		TEST_METHOD(Perform_FEM_test14) {
			std::string objPath = "..\\..\\tests\\test14\\remeshed_surface.STL";
			std::string outDir = "..\\..\\tests\\test14\\";
			std::string forcesPath = "..\\..\\tests\\test14\\load_nodes.txt";
			std::string fixedPath = "..\\..\\tests\\test14\\fixed_nodes.txt";
			std::string clampedPath = "..\\..\\tests\\test14\\clamped_edges.txt";
			std::string fixedPath = "..\\..\\tests\\test14\\fixed_nodes.txt";
			JobProperties jobProps("test", outDir, objPath, forcesPath, fixedPath, clampedPath, true, STL_T);
			SimulationProperties simProps;

			runFEMJob(jobProps, simProps, true);
		}
		*/










		// old tests:



		/*
		TEST_METHOD(Perform_FEM_test1)
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXd TC;
			Eigen::MatrixXd N;
			Eigen::MatrixXi F;
			Eigen::MatrixXi FTC;
			Eigen::MatrixXi FN;
			std::string objPath = "..\\tests\\test1\\shelve.obj";
			std::string outputObjPath = "..\\tests\\test1\\shelve_output.obj";
			std::string stdoutPath = "..\\tests\\test1\\test1.log";
			std::string nodalForcesPath = "..\\..\\tests\\test1\\load_nodes.txt";
			std::string fixedNodesPath = "..\\tests\\test1\\fixed_nodes.txt";
			FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

			igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
			auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
			auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
			FEMData data; //gets default data (for now)

			FEMResults result;
			Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
			saveOBJ(result.displacedVertices, F, outputObjPath);

			Viewer viewer;
			viewer.startView(result.displacedVertices, F, result.vonMisesStress);
			fclose(file);
		}

		TEST_METHOD(Perform_FEM_test2)
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXd TC;
			Eigen::MatrixXd N;
			Eigen::MatrixXi F;
			Eigen::MatrixXi FTC;
			Eigen::MatrixXi FN;
			std::string objPath = "..\\..\\tests\\test2\\pyramid.obj";
			std::string outputObjPath = "..\\..\\tests\\test2\\pyramid_output.obj";
			std::string stdoutPath = "..\\..\\tests\\test2\\test2.log";
			std::string nodalForcesPath = "..\\..\\tests\\test2\\load_nodes.txt";
			std::string fixedNodesPath = "..\\..\\tests\\test2\\fixed_nodes.txt";
			FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

			igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
			auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
			auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
			FEMData data; //gets default data (for now)

			FEMResults result;
			Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
			saveOBJ(result.displacedVertices, F, outputObjPath);

			Viewer viewer;
			viewer.startView(result.displacedVertices, F, result.vonMisesStress);
			fclose(file);
		}

		TEST_METHOD(Perform_FEM_test3)
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXd TC;
			Eigen::MatrixXd N;
			Eigen::MatrixXi F;
			Eigen::MatrixXi FTC;
			Eigen::MatrixXi FN;
			std::string objPath = "..\\..\\tests\\test3\\plate.obj";
			std::string outputObjPath = "..\\..\\tests\\test3\\plate_output.obj";
			std::string stdoutPath = "..\\..\\tests\\test3\\test3.log";
			std::string nodalForcesPath = "..\\..\\tests\\test3\\load_nodes.txt";
			std::string fixedNodesPath = "..\\..\\tests\\test3\\fixed_nodes.txt";
			FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

			igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
			auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
			auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
			FEMData data; //gets default data (for now)

			FEMResults result;
			Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
			saveOBJ(result.displacedVertices, F, outputObjPath);

			Viewer viewer;
			viewer.startView(result.displacedVertices, F, result.vonMisesStress);
			fclose(file);
		}

		TEST_METHOD(Perform_FEM_test5)
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXd TC;
			Eigen::MatrixXd N;
			Eigen::MatrixXi F;
			Eigen::MatrixXi FTC;
			Eigen::MatrixXi FN;
			std::string objPath = "..\\..\\tests\\test5\\rod.obj";
			std::string outputObjPath = "..\\..\\tests\\test5\\rod_output.obj";
			std::string stdoutPath = "..\\..\\tests\\test5\\test5.log";
			std::string nodalForcesPath = "..\\..\\tests\\test5\\load_nodes.txt";
			std::string fixedNodesPath = "..\\..\\tests\\test5\\fixed_nodes.txt";
			FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

			igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
			auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
			auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
			FEMData data; //gets default data (for now)

			FEMResults result;
			Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
			saveOBJ(result.displacedVertices, F, outputObjPath);

			Viewer viewer;
			viewer.startView(result.displacedVertices, F, result.vonMisesStress);
			fclose(file);
		}

		TEST_METHOD(Perform_FEM_test6)
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXd TC;
			Eigen::MatrixXd N;
			Eigen::MatrixXi F;
			Eigen::MatrixXi FTC;
			Eigen::MatrixXi FN;
			std::string objPath = "..\\..\\tests\\test6\\box.obj";
			std::string outputObjPath = "..\\..\\tests\\test6\\box_output.obj";
			std::string stdoutPath = "..\\..\\tests\\test6\\test6.log";
			std::string nodalForcesPath = "..\\..\\tests\\test6\\load_nodes.txt";
			std::string fixedNodesPath = "..\\..\\tests\\test6\\fixed_nodes.txt";
			FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

			igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
			auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
			auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
			FEMData data; //gets default data (for now)

			FEMResults result;
			Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
			saveOBJ(result.displacedVertices, F, outputObjPath);

			Viewer viewer;
			viewer.startView(result.displacedVertices, F, result.vonMisesStress);
			fclose(file);
		}

		TEST_METHOD(Perform_FEM_test7)
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXd TC;
			Eigen::MatrixXd N;
			Eigen::MatrixXi F;
			Eigen::MatrixXi FTC;
			Eigen::MatrixXi FN;
			std::string objPath = "..\\..\\tests\\test7\\sphere.obj";
			std::string outputObjPath = "..\\..\\tests\\test7\\sphere_output.obj";
			std::string stdoutPath = "..\\..\\tests\\test7\\test7.log";
			std::string nodalForcesPath = "..\\..\\tests\\test7\\load_nodes.txt";
			std::string fixedNodesPath = "..\\..\\tests\\test7\\fixed_nodes.txt";
			FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

			igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
			auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
			auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
			FEMData data; //gets default data (for now)

			FEMResults result;
			Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
			saveOBJ(result.displacedVertices, F, outputObjPath);

			Viewer viewer;
			viewer.startView(result.displacedVertices, F, result.vonMisesStress);
			fclose(file);
		}

		TEST_METHOD(Perform_FEM_test8)
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXd TC;
			Eigen::MatrixXd N;
			Eigen::MatrixXi F;
			Eigen::MatrixXi FTC;
			Eigen::MatrixXi FN;
			std::string objPath = "..\\..\\tests\\test8\\rod.obj";
			std::string outputObjPath = "..\\..\\tests\\test8\\rod_twist_output.obj";
			std::string stdoutPath = "..\\..\\tests\\test8\\test8.log";
			std::string nodalForcesPath = "..\\..\\tests\\test8\\load_nodes.txt";
			std::string fixedNodesPath = "..\\..\\tests\\test8\\fixed_nodes.txt";
			FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

			igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
			auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
			auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
			FEMData data; //gets default data (for now)

			FEMResults result;
			Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
			saveOBJ(result.displacedVertices, F, outputObjPath);

			Viewer viewer;
			viewer.startView(result.displacedVertices, F, result.vonMisesStress);
			fclose(file);
		}

		TEST_METHOD(Perform_FEM_test9)
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXd TC;
			Eigen::MatrixXd N;
			Eigen::MatrixXi F;
			Eigen::MatrixXi FTC;
			Eigen::MatrixXi FN;
			std::string objPath = "..\\..\\tests\\test9\\plate_fold.obj";
			std::string outputObjPath = "..\\..\\tests\\test9\\plate_fold_output.obj";
			std::string stdoutPath = "..\\..\\tests\\test9\\test9.log";
			std::string nodalForcesPath = "..\\..\\tests\\test9\\load_nodes.txt";
			std::string fixedNodesPath = "..\\..\\tests\\test9\\fixed_nodes.txt";
			FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

			igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
			auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
			auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
			FEMData data; //gets default data (for now)

			FEMResults result;
			Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
			saveOBJ(result.displacedVertices, F, outputObjPath);

			Viewer viewer;
			viewer.startView(result.displacedVertices, F, result.vonMisesStress);
			fclose(file);
		}

		TEST_METHOD(Perform_FEM_test10)
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXd TC;
			Eigen::MatrixXd N;
			Eigen::MatrixXi F;
			Eigen::MatrixXi FTC;
			Eigen::MatrixXi FN;
			std::string objPath = "..\\..\\tests\\test10\\Shelf_ElementSizing_30.obj";
			std::string outputObjPath = "..\\..\\tests\\test10\\shelf_30_output.obj";
			std::string stdoutPath = "..\\..\\tests\\test10\\test10.log";
			std::string nodalForcesPath = "..\\..\\tests\\test10\\load_nodes_30.txt";
			std::string fixedNodesPath = "..\\..\\tests\\test10\\fixed_nodes_30.txt";
			FILE *file = freopen(stdoutPath.c_str(), "w", stdout); // setting stdout

			igl::readOBJ(objPath, V, TC, N, F, FTC, FN);
			auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
			auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
			FEMData data; //gets default data (for now)

			FEMResults result;
			Perform_FEM(Mesh(V, F, fixedNodes), nodalForces, data, result);
			saveOBJ(result.displacedVertices, F, outputObjPath);

			Viewer viewer;
			viewer.startView(result.displacedVertices, F, result.vonMisesStress);
			fclose(file);
		}

		*/
	};
}