#include "stdafx.h"
#include <Eigen/core>
#include "utils.cpp"
#include "FEM.cpp"
#include "CppUnitTest.h"

#define EpsVal 1e-6

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace ShellFEMUnitTests
{		
	TEST_CLASS(FEM)
	{
	public:
		TEST_METHOD(Perform_FEM_test1)
		{
			Eigen::MatrixXd meshedV, meshedN;
			Eigen::MatrixXi meshedF;
			std::string meshPath =				"..\\..\\tests\\test1\\input\\remeshed_surface.stl";
			std::string nodalForcesPath =		"..\\..\\tests\\test1\\input\\load_nodes.txt";
			std::string fixedNodesPath =		"..\\..\\tests\\test1\\input\\fixed_nodes.txt";
			std::string testDisplacmentsPath =	"..\\..\\tests\\test1\\output\\displacements.txt";
			std::string testVonMisesPath =		"..\\..\\tests\\test1\\output\\stresses.txt";
			igl::readSTL(meshPath, meshedV, meshedF, meshedN);
			auto nodalForces = nodal_forces_from_txt(nodalForcesPath);
			auto fixedNodes = fixed_nodes_from_txt(fixedNodesPath);
			auto testDisplacements = displacements_from_txt(testDisplacmentsPath, meshedV.rows());
			auto testVonMisesStresses = vonmises_from_txt(testDisplacmentsPath, meshedV.rows());
			FEMData data; //gets default data (for now)
			auto FEMResults = Perform_FEM(meshedV, meshedF, nodalForces, fixedNodes, data);
			//get calculated FEM results
			auto resDisplacements = FEMResults.displacement;
			auto resVonMises = FEMResults.vonMisesStress;
			//check result dimensions
			Assert::AreEqual((size_t)testDisplacements.rows(),		(size_t)resDisplacements.rows());
			Assert::AreEqual((size_t)testVonMisesStresses.rows(),	(size_t)resVonMises.rows());
			Assert::AreEqual((size_t)testDisplacements.cols(),		(size_t)resDisplacements.cols());
			Assert::AreEqual((size_t)testVonMisesStresses.cols(),	(size_t)resVonMises.cols());
			//check displacements
			Assert::IsTrue((testDisplacements - resDisplacements).norm() < EpsVal);
			//check stresses
			Assert::IsTrue((testVonMisesStresses - resVonMises).norm() < EpsVal);
		}

	};
}