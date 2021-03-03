#pragma once

#include "utils.h"
#include "vector_utils.h"
#include <utility>

#define NXT(i) ((i + 1) % 3)	// next index in counter-clockwise direction with relation to normal
#define PRV(i) ((i + 2) % 3)	// previous index in counter-clockwise direction with relation to normal

/**
	Data Structure For Representing S3 Triangle Element 
*/
struct Element {
	Vector3d vertices[3];					// List of vertices (x,y,z)
	bool neighborExists[3];					// Sometimes there is no neighbor
	bool isEdgeFixed[3];					// Some vertices are fixed and cannot move.
	Vector3d neighborVertices[3];			// List of neighbor vertices (x,y,z). Neighbor at the opposite side of vertex i or vector3d
	Eigen::Matrix<double, 18, 18> Ke;		// Element Stiffness Matrix
	double vonMisesStress;					// Von Mises yield criterion calculated from dispositions

	Element(Vector3d _vertices[]) {
		for (int i = 0; i < 3; i++) {
			vertices[i] = _vertices[i];
			neighborExists[i] = false;
			isEdgeFixed[i] = false;
		}
	}

	void setNeighbour(Vector3d v, int i) {
		neighborExists[i] = true;
		neighborVertices[i] = v;
	}

	void setFixedNode(int i) {
		isEdgeFixed[NXT(i)] = false;//true;
		isEdgeFixed[PRV(i)] = false;// true;
	}
};

/**
	Data Structures For Holding Temporary Results Used to Calculate the Stiffness Matrix
*/
	struct NeighborParameters {
		double height;							// Height down from vertex
		Vector3d normal;						// Normal to triangle face
		double heightArr[3];					// Index is according to adjacent element vertex index
		double cosineArr[3];					// Index is according to adjacent element vertex index
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
	double thickness;				// Shell thickness
	
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
	ElementBuilder() {};
	ElementBuilder(Eigen::Matrix<double, 3, 3> const &_De, Eigen::Matrix<double, 3, 3> const &_M, double _thickness) : De(_De), M(_M), thickness(_thickness) {};
	void calculateStiffnessMatrix(Element &element);
	void calculateVonMisesStress(Element &element, Eigen::Matrix<double, 18, 1> const globalDisplacement);
};