#pragma once

#include "utils.h"
#include <utility>

#define NXT(i) ((i + 1) % 3)	// next index in counter-clockwise direction with relation to normal
#define PRV(i) ((i + 2) % 3)	// previous index in counter-clockwise direction with relation to normal
#define DEFAULT_YOUNG_CNST		70e9
#define DEFAULT_POSSION_CNST	0.3
#define DEFAULT_THICKNESS_CNST	1.6e-3

struct SimulationProperties {
	double E; //Young's Modulus
	double ni; //Possion's Ratio
	double thickness;
	SimulationProperties() : E(DEFAULT_YOUNG_CNST), ni(DEFAULT_POSSION_CNST), thickness(DEFAULT_THICKNESS_CNST) {};
	SimulationProperties(double E, double ni, double thickness) : E(E), ni(ni), thickness(thickness) {};
};

/**
	Data Structure For Representing S3 Triangle Element.
	For Face with vertices [0,1,2] we order edges/neighbors according to [(0,1), (1,2), (2,0)].
	i.e. neighbor 0 is sharing the edge (0,1) with the element.
*/
struct Element {
	Vector3d vertices[3];					// List of vertices (x,y,z).
	bool neighborExists[3];					// Sometimes there is no neighbor.
	bool isEdgeClamped[3];					// Some edges need to be set as clamped.
	Vector3d neighborVertices[3];			// List of neighbor vertices (x,y,z).
	Eigen::Matrix<double, 18, 18> Ke;		// Element Stiffness Matrix.
	double vonMisesStress;					// Von Mises yield criterion calculated from displacements.

	Element(Vector3d _vertices[]) {
		for (int i = 0; i < 3; i++) {
			vertices[i] = _vertices[i];
			neighborExists[i] = false;
			isEdgeClamped[i] = false;
			neighborVertices[i] << 0, 0, 0;
		}
	}

	void setNeighbour(Vector3d v, int i) {
		neighborExists[i] = true;
		neighborVertices[i] = v;
	}

};

/**
	Data Structures For Holding Temporary Results Used to Calculate the Stiffness Matrix
*/
	struct NeighborParameters {
		double height = 0;					// Height down from vertex
		Vector3d normal = Vector3d(0,0,0);	// Normal to triangle face
		double heightArr[3] = {0};			// Index is according to adjacent element vertex index
		double cosineArr[3] = {0};			// Index is according to adjacent element vertex index
	};

	struct ElementParameters {
		double heights[3];						// Height down from vertex i
		double cosine[3];						// Cosine of angle at vertex i
		double c[3];							// TODO : find definition 
		double s[3];							// TODO : find definition
		double area;							// Element triangle area
		Vector3d axes[3];						// Local plane axes. normal to plane is Z direction.
		NeighborParameters neighborParam[3];	// List of neighbor triangles' fields
	};

/**
	This Class Receives an Element Struct and Calculates It's Stiffness Matrix
*/
class ElementBuilder { 
private:					
	Eigen::Matrix<double, 3, 3> De; // Elastic Plastic Behavior Matrix
	Eigen::Matrix<double, 3, 3> M;  // Hill's plastic strain matrix
	SimulationProperties simProps;
	
	void calculateElasticPlasticMatrix();
	void calculateHillPlasticStrainMatrix();
	void getUnitVectors(Element const &element, ElementParameters &elemParam);
	void calculateParameters(Element const &element, ElementParameters &elemParam);
	void buildRMatrix(ElementParameters const &elemParam, Eigen::Matrix<double, 3, 3> &R);		// TODO change these methods names
	void buildHMatrix(Element const &element, ElementParameters const &elemParam, Eigen::Matrix<double, 3, 6>  &H);
	void buildCMatrixElementRow(Eigen::Matrix<double, 1, 18> &row, ElementParameters const &elemParam,
									int idx1, int idx2, int idx3, int pos1, int pos2, int pos3);
	void buildCMatrixNeighborRow(Eigen::Matrix<double, 1, 18> &row, NeighborParameters const &param,
									int idx1, int idx2, int pos1, int pos2, int pos3);
	void buildCMatrix(Element const &element, ElementParameters const &elemParam, Eigen::Matrix<double, 6, 18> &C);
	void calculateBmMatrix(ElementParameters const &elemParam, Eigen::Matrix<double, 3, 9> &Bm);
	void calculateBMatrix(Element &element, ElementParameters const &elemParam, Eigen::Matrix<double, 3, 18> &B);
public:
	ElementBuilder(SimulationProperties const &_simProps);
	void calculateStiffnessMatrix(Element &element);
	void calculateVonMisesStress(Element &element, Eigen::Matrix<double, 18, 1> const globalDisplacement);
};