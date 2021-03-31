#include "s3element.h"
#define DEBUG
#define X 0
#define Y 1
#define Z 2

//######################################################################### Auxiliary functions ############################################################################

std::ostream& operator<<(std::ostream& os, const Element& elem) {
	std::stringstream repr;

	for (int i = 0; i < 3; i++) repr << "vertices[" << i << "]: " << elem.vertices[i] << std::endl;
	repr << "centroid: " << (elem.vertices[0] + elem.vertices[1] + elem.vertices[2]) / 3 << std::endl;
	for (int i = 0; i < 3; i++) repr << "neighborExists[" << i << "]: " << elem.neighborExists[i] << std::endl;
	for (int i = 0; i < 3; i++) repr << "neighborVertices[" << i << "]: " << elem.neighborVertices[i] << std::endl;
	return os << repr.str();
}

std::ostream& operator<<(std::ostream& os, const NeighborParameters& param) {
	std::stringstream repr;

	for (int i = 0; i < 3; i++) repr << "heightArr[" << i << "]: " << param.heightArr[i] << std::endl;
	for (int i = 0; i < 3; i++) repr << "cosineArr[" << i << "]: " << param.cosineArr[i] << std::endl;
	repr << "height: " << param.height << std::endl;
	repr << "normal: " << std::endl << param.normal << std::endl;

	return os << repr.str();
}

std::ostream& operator<<(std::ostream& os, const ElementParameters& param) {
	std::stringstream repr;

	repr << DASH << std::endl << "Element: " << std::endl;
	for (int i = 0; i < 3; i++) repr << "heights["<<i<<"]: " << param.heights[i] << std::endl;
	for (int i = 0; i < 3; i++) repr << "cosine["<<i<<"]: " << param.cosine[i]<< std::endl;
	for (int i = 0; i < 3; i++) repr << "c["<<i<<"]: " << param.c[i] << std::endl;
	for (int i = 0; i < 3; i++) repr << "s["<<i<<"]: " << param.s[i] << std::endl;
	for (int i = 0; i < 3; i++) repr << "axes["<<i<<"]: " << std::endl << param.axes[i] << std::endl;
	repr << "area: " << param.area << std::endl;
	repr << DASH << std::endl << "Neighbors: " << std::endl;
	for (int i = 0; i < 3; i++) repr << "neighbor[" << i << "]" << std::endl << param.neighborParam[i] << std::endl;
	repr << DASH << std::endl;

	return os << repr.str();
}

//################################################################# Element Builder Initialization ############################################################################

void ElementBuilder::calculateElasticPlasticMatrix() {
	double ni = simProps.ni;
	double E = simProps.E;
	double t = simProps.thickness;
	De << 1, ni, 0,
		ni, 1, 0,
		0, 0, (1 - ni) / 2;
	De *= E / (1 - ni * ni);
	De *= t * t*t / 12;
}

void ElementBuilder::calculateHillPlasticStrainMatrix() {
	double f = 0.5, g = 0.5, h = 0.5, n = 1.5; // set to values which calculate von mises yield criterion.
	M << g + h, -h, 0,
		-h, f + h, 0,
		0, 0, 2 * n;
}

ElementBuilder::ElementBuilder(SimulationProperties const &_simProps) {
	simProps = _simProps;
	calculateElasticPlasticMatrix();
	calculateHillPlasticStrainMatrix();
}; 

//################################################################# Element parameters Calculation ############################################################################

void ElementBuilder::getUnitVectors(Element const &element, ElementParameters &elemParam) {
	std::pair<double, Eigen::Vector3d> areaVector;
	Eigen::Vector3d centroid;
	Eigen::Vector2d planeVertices[3];

	areaVector = getAreaVector(element.vertices[0], element.vertices[NXT(0)], element.vertices[PRV(0)]);
	elemParam.area = areaVector.first;
	elemParam.axes[Z] = areaVector.second; 
	centroid = (element.vertices[0] + element.vertices[1] + element.vertices[2]) / 3;
	elemParam.axes[Y] = (element.vertices[0] - centroid).normalized();
	elemParam.axes[X] = elemParam.axes[Y].cross(elemParam.axes[Z]);
}

void ElementBuilder::calculateParameters(Element const &element, ElementParameters &elemParam) {	
	
	std::pair<double, Eigen::Vector3d> areaVector;
	Eigen::Vector3d centroid;
	Eigen::Vector2d planeVertices[3]; //needs to be 3d?
	
	getUnitVectors(element, elemParam);
	
	for (int i = 0; i < 3; i++) {
		planeVertices[i] << Eigen::Vector2d(elemParam.axes[X].dot(element.vertices[i]), elemParam.axes[Y].dot(element.vertices[i]));
	}
	
	for (int i=0; i<3; i++) {
		// Element Fields
		double length = (element.vertices[NXT(i)] - element.vertices[PRV(i)]).norm();
		elemParam.c[i] = + (planeVertices[PRV(i)][Y] - planeVertices[NXT(i)][Y]) / length;
		elemParam.s[i] = - (planeVertices[PRV(i)][X] - planeVertices[NXT(i)][X]) / length;
		elemParam.cosine[i] = getCosOfAngle(element.vertices[i], element.vertices[NXT(i)], element.vertices[PRV(i)]);
		elemParam.heights[i] = 2 * elemParam.area / length;

		// Neighbor Fields
		if (element.neighborExists[i]) {
			areaVector = getAreaVector(element.neighborVertices[i], element.vertices[NXT(i)], element.vertices[i]);
			double neighborArea = areaVector.first;
			elemParam.neighborParam[i].normal = areaVector.second;
			elemParam.neighborParam[i].cosineArr[i] = getCosOfAngle(element.vertices[i], element.neighborVertices[i], element.vertices[NXT(i)]);
			elemParam.neighborParam[i].cosineArr[NXT(i)] = getCosOfAngle(element.vertices[NXT(i)], element.vertices[i], element.neighborVertices[i]);
			elemParam.neighborParam[i].height = 2 * neighborArea / (element.vertices[NXT(i)] - element.vertices[i]).norm();
			elemParam.neighborParam[i].heightArr[i] = 2 * neighborArea / (element.vertices[NXT(i)] - element.neighborVertices[i]).norm();
			elemParam.neighborParam[i].heightArr[NXT(i)] = 2 * neighborArea / (element.vertices[i] - element.neighborVertices[i]).norm();
		}
	}
}

//################################################################# Element Matrices Calculation ############################################################################

void ElementBuilder::buildRMatrix(ElementParameters const &elemParam, Eigen::Matrix<double, 3, 3> &R) {
	R << (1 * elemParam.c[0] * elemParam.c[0]), (1 * elemParam.c[1] * elemParam.c[1]), (1 * elemParam.c[2] * elemParam.c[2]),
		 (1 * elemParam.s[0] * elemParam.s[0]), (1 * elemParam.s[1] * elemParam.s[1]), (1 * elemParam.s[2] * elemParam.s[2]),
		 (2 * elemParam.c[0] * elemParam.s[0]), (2 * elemParam.c[1] * elemParam.s[1]), (2 * elemParam.c[2] * elemParam.s[2]);
}

void ElementBuilder::buildHMatrix(Element const &element, ElementParameters const &elemParam, Eigen::Matrix<double, 3, 6>  &H) {
	double tmp[3];			// Temps

	for (int i = 0; i < 3; i++) {
		if (element.neighborExists[NXT(i)])
			tmp[i] = -2 / (elemParam.heights[i] + elemParam.neighborParam[NXT(i)].height);
		else tmp[i] = -1 / (elemParam.heights[i]);
	}
	
	H << tmp[0], 0,		 0,		 0,		 tmp[0], 0,
		 0,		 tmp[1], 0,		 0,		 0,		 tmp[1],
		 0,		 0,		 tmp[2], tmp[2], 0,		 0;
}

void ElementBuilder::buildCMatrixElementRow(Eigen::Matrix<double, 1, 18> &row, ElementParameters const &elemParam,
												int idx1, int idx2, int idx3, int pos1, int pos2, int pos3) {
	double t1, t2; // Temps

	t1 = elemParam.cosine[idx1] / elemParam.heights[idx2];
	t2 = elemParam.cosine[idx2] / elemParam.heights[idx1];

	row.block(0, pos1 * 3, 1, 3) << - (t1 * elemParam.axes[Z][X]), \
									- (t1 * elemParam.axes[Z][Y]), \
									- (t1 * elemParam.axes[Z][Z]);
	row.block(0, pos2 * 3, 1, 3) << - (t2 * elemParam.axes[Z][X]), \
									- (t2 * elemParam.axes[Z][Y]), \
									- (t2 * elemParam.axes[Z][Z]);
	row.block(0, pos3 * 3, 1, 3) << (elemParam.axes[Z][X] / elemParam.heights[idx3]), \
									(elemParam.axes[Z][Y] / elemParam.heights[idx3]), \
									(elemParam.axes[Z][Z] / elemParam.heights[idx3]);
	row.block(0, 9, 1, 9) << 0, 0, 0, 0, 0, 0, 0, 0, 0;
}

void ElementBuilder::buildCMatrixNeighborRow(Eigen::Matrix<double, 1, 18> &row, NeighborParameters const &param,
												int idx1, int idx2, int pos1, int pos2, int pos3) {
	double t1, t2; // Temps
	row = Eigen::Matrix<double, 1, 18>::Zero();

	t1 = param.cosineArr[idx1] / param.heightArr[idx2];
	t2 = param.cosineArr[idx2] / param.heightArr[idx1];

	row.block(0, pos1 * 3, 1, 3) << - (t1 * param.normal[X]), \
									- (t1 * param.normal[Y]), \
									- (t1 * param.normal[Z]);

	row.block(0, pos2 * 3, 1, 3) << - (t2 * param.normal[X]), \
									- (t2 * param.normal[Y]), \
									- (t2 * param.normal[Z]);

	row.block(0, pos3 * 3, 1, 3) << (param.normal[X] / param.height), \
									(param.normal[Y] / param.height), \
									(param.normal[Z] / param.height);
}

void ElementBuilder::buildCMatrix(Element const &element, ElementParameters const &elemParam, Eigen::Matrix<double, 6, 18> &C) {
	Eigen::Matrix<double, 1, 18> row[6];	
	
	buildCMatrixElementRow(row[0], elemParam, 2, 1, 0, 1, 2, 0);
	buildCMatrixElementRow(row[1], elemParam, 2, 0, 1, 0, 2, 1);
	buildCMatrixElementRow(row[2], elemParam, 1, 0, 2, 0, 1, 2);

	std::cout << "clamped 0: " << element.isEdgeClamped[0] << std::endl;
	std::cout << "clamped 1: " << element.isEdgeClamped[1] << std::endl;
	std::cout << "clamped 2: " << element.isEdgeClamped[2] << std::endl;

	if (element.neighborExists[0])
		buildCMatrixNeighborRow(row[3], elemParam.neighborParam[0], 1, 0, 0, 1, 3);
	else if (element.isEdgeClamped[0]) row[3] = row[2];
	else row[3] = -row[2];

	if (element.neighborExists[1])
		buildCMatrixNeighborRow(row[4], elemParam.neighborParam[1], 2, 1, 1, 2, 4);
	else if (element.isEdgeClamped[1]) row[4] = row[0];
	else row[4] = -row[0];

	if (element.neighborExists[2])
		buildCMatrixNeighborRow(row[5], elemParam.neighborParam[2], 2, 0, 0, 2, 5);
	else if (element.isEdgeClamped[2]) row[5] = row[1];
	else row[5] = -row[1];

	for (int i = 0; i < 6; i++) {
		C.block(i, 0, 1, 18) = row[i];
	}
}

void ElementBuilder::calculateBmMatrix(ElementParameters const &elemParam, Eigen::Matrix<double, 3, 9> &Bm) {
	double t1, t2; // Temps

	for (int i = 0; i < 3; i++) {
		t1 = elemParam.c[i] / elemParam.heights[i];
		t2 = elemParam.s[i] / elemParam.heights[i];
		Eigen::Matrix<double, 3, 3> BmBlock;
		BmBlock <<  (t1*elemParam.axes[X][X]), (t1*elemParam.axes[X][Y]), (t1*elemParam.axes[X][Z]),
					(t2*elemParam.axes[Y][X]), (t2*elemParam.axes[Y][Y]), (t2*elemParam.axes[Y][Z]),
					(t1*elemParam.axes[Y][X]+t2* elemParam.axes[X][X]), \
						(t1*elemParam.axes[Y][Y] + t2 * elemParam.axes[X][Y]), \
						(t1*elemParam.axes[Y][Z] + t2 * elemParam.axes[X][Z]);
		Bm.block(0, i * 3, 3, 3) = BmBlock;
	}
	Bm *= -1;
}

void ElementBuilder::calculateBMatrix(Element &element, ElementParameters const &elemParam, Eigen::Matrix<double, 3, 18> &B) {
	Eigen::Matrix<double, 3, 3>  R;			// TODO : find description
	Eigen::Matrix<double, 3, 6>  H;			// TODO : find description
	Eigen::Matrix<double, 6, 18> C;			// TODO : find description

	buildRMatrix(elemParam, R);
	buildHMatrix(element, elemParam, H);
	buildCMatrix(element, elemParam, C);

	B = R * H * C;					// Bending Effect
#ifdef DEBUG
	std::cout << "R: " << std::endl << R << std::endl;
	std::cout << "H: " << std::endl << H << std::endl;
	std::cout << "C: " << std::endl << C << std::endl;
	std::cout << "HC: " << std::endl << H * C << std::endl;
	std::cout << "B: " << std::endl << B << std::endl;
	std::cout << DASH << std::endl;
#endif
}

//################################################################# Element Main Methods ############################################################################

void ElementBuilder::calculateStiffnessMatrix(Element &element) {
	Eigen::Matrix<double, 3, 18> B;			// Bending effect matrix
	Eigen::Matrix<double, 3, 9>  Bm;		// Membrane Effect Matrix
	Eigen::Matrix<double, 9, 9>  Kem;
	ElementParameters elemParam;

	calculateParameters(element, elemParam);
	calculateBMatrix(element, elemParam, B);
	calculateBmMatrix(elemParam, Bm); 

	Kem = simProps.thickness * elemParam.area * Bm.transpose() * De * Bm;
	element.Ke = simProps.thickness * elemParam.area * B.transpose() * De * B;
	element.Ke.block(0, 0, 9, 9) += Kem;	

#ifdef DEBUG
	std::cout << element << std::endl;
	std::cout << elemParam << std::endl;
	std::cout << "De: " << std::endl << De << std::endl;
	std::cout << "Bm: " << std::endl << Bm << std::endl;
	std::cout << "Kem: " << std::endl << Kem << std::endl;
	std::cout << "Ke: " << std::endl << element.Ke << std::endl;
	std::cout << DASH << std::endl;
#endif
}

void ElementBuilder::calculateVonMisesStress(Element &element, Eigen::Matrix<double, 18, 1> const displacements) {
	Eigen::Matrix<double, 3, 18> B;			// Bending effect matrix
	Eigen::Matrix<double, 3, 9> Bm;		    // Membrane Effect Matrix	
	Eigen::Matrix<double, 3, 1> strain;
	Eigen::Matrix<double, 3, 1> stress;
	ElementParameters elemParam;
	
	calculateParameters(element, elemParam);	
	calculateBMatrix(element, elemParam, B);
	calculateBmMatrix(elemParam, Bm);

	strain = B * displacements + Bm * displacements.block(0,0,9,1);
	stress = De * strain;
	element.vonMisesStress = stress.transpose() * M * stress;
#ifdef DEBUG
	if (displacements.sum() != 0) {
		std::cout << "calculateVonMisesStress" << std::endl;
		std::cout << "elemParam " << std::endl << elemParam << std::endl;
		std::cout << "Bm " << std::endl << Bm << std::endl;
		std::cout << "displacements " << std::endl << displacements << std::endl;
		std::cout << "Strain: " << std::endl << strain << std::endl;
		std::cout << "Stress: " << std::endl << stress << std::endl;
		std::wcout << "vonMisesStress: " << element.vonMisesStress << std::endl;
		std::cout << DASH << std::endl;
	}	
#endif
}

