#include "hybrid_composition3.h"
#include "../rng.h"
#include "../problem.h"
#include "../utils.h"

#include "rastrigin.h"
#include "schafferf6.h"
#include "griewank.h"
#include "griewank_rosenbrock.h"
#include "weierstrass.h"

#include <cmath>
#include <string>
#include <sstream>
#include <cfloat>

using namespace std;

double HC_3::sigma[NUM_OF_FUNC] = {
		1.0,	1.0,	1.0,	1.0,	1.0,
		2.0,	2.0,	2.0,	2.0,	2.0

};

double HC_3::lambda[NUM_OF_FUNC] = {
		5.0*5.0/100.0,	5.0/100.0,	5.0*1.0,	1.0,			5.0*1.0,
		1.0,			5.0*10.0,	10.0,		5.0*5.0/200.0,	5.0/200.0

};

double HC_3::biases[NUM_OF_FUNC] = {
		0.0,	100.0,	200.0,	300.0,	400.0,
		500.0,	600.0,	700.0,	800.0,	900.0
};


long double HC_3::hybrid_composition(int dim, const double* x) {

	// Get the raw weights
	double wMax = DBL_MIN;

	for (int i = 0 ; i < num_func ; i ++) {
		double sumSqr = 0.0;
		shift(dim, z2d[i], x, shift_vector2D[i]);
		for (int j = 0 ; j < dim ; j ++) {
			sumSqr += (z2d[i][j] * z2d[i][j]);
		}
		w[i] = exp(-1.0 * sumSqr / (2.0 * dim * sigma[i] * sigma[i]));
		if (wMax < w[i])
			wMax = w[i];
	}

	// Modify the weights
	double wSum = 0.0;
	double w1mMaxPow = 1.0 - pow(wMax, 10.0);
	for (int i = 0 ; i < num_func ; i ++) {
		if (w[i] != wMax) {
			w[i] *= w1mMaxPow;
		}
		wSum += w[i];
	}

	// Normalize the weights
	for (int i = 0 ; i < num_func ; i ++) {
		w[i] /= wSum;
	}

	long double sumF = 0.0;
	for (int i = 0 ; i < num_func ; i ++) {
		for (int j = 0 ; j < dim ; j ++) {
			z2d[i][j] /= lambda[i];
		}

		rotate(dim, zM[i], z2d[i], rotation_matrix3D[i]);
		sumF +=	w[i] * (C * basic_func(i, dim, zM[i]) / fmax[i] + biases[i]);
	}
	return (sumF);
}

HC_3::HC_3(Configuration* config, int variantID):Problem(config, variantID){


	C = 2000.0;
	num_func = NUM_OF_FUNC;

	zM = allocateMemory2D(dimension,NUM_OF_FUNC);
	z2d = allocateMemory2D(dimension,NUM_OF_FUNC);
	rotation_matrix3D = allocateMemory3D(dimension,NUM_OF_FUNC,dimension);

	shift_vector2D = allocateMemory2D(dimension,NUM_OF_FUNC);

	string file_data = "supportData/hybrid_func3_data.txt";

	string file_m;
	stringstream dim_name;
	dim_name << dimension;

	if(vID == ROTATED_WITH_HIGH_CONDITION_NUMBER_MATRIX)
		file_m = "supportData/hybrid_func3_HM_D" + dim_name.str() + ".txt";
	else if (vID == ROTATED || vID == NONCONTINUOUS_ROTATED)
		file_m = "supportData/hybrid_func3_M_D" + dim_name.str() + ".txt";

	Utils::loadMatrixFromFile(file_data, NUM_OF_FUNC, dimension, shift_vector2D);

	Utils::loadMatrixFromFile(file_m, NUM_OF_FUNC, dimension, dimension, rotation_matrix3D);

	w = new double[NUM_OF_FUNC];

	double testPoint[dimension];
	double testPointM[dimension];

	fmax = new double[NUM_OF_FUNC];

	// Calculate/estimate the fmax for all the functions involved
	for (int i = 0 ; i < NUM_OF_FUNC ; i ++) {
		for (int j = 0 ; j < dimension ; j ++) {
			testPoint[j] = (5.0 / lambda[i]);
		}
		rotate(dimension, testPointM, testPoint, rotation_matrix3D[i]);
		fmax[i] = fabs(basic_func(i, dimension, testPointM));
	}
}


HC_3::~HC_3(){
	deallocateMemory3D(rotation_matrix3D, NUM_OF_FUNC, dimension);
	deallocateMemory2D(zM, NUM_OF_FUNC);
	deallocateMemory2D(z2d, NUM_OF_FUNC);
	deallocateMemory2D(shift_vector2D, NUM_OF_FUNC);
	delete [] fmax;
	delete [] w;
}


long double HC_3::basic_func(int func_no, int dim, const double* x){
	long double result = 0.0;

	switch(func_no) {
	case 0:
	case 1:{
		SchafferF6 * problem = new SchafferF6(config, EXPANDED);
		result = problem->evaluate(dim, x);
		delete problem;
		break;}
	case 2:
	case 3:{
		Rastrigin * problem = new Rastrigin(config, BASIC);
		result = problem->evaluate(dim, x);
		delete problem;
		break;}
	case 4:
	case 5:{
		GriewankRosenbrock * problem = new GriewankRosenbrock(config, BASIC);
		result = problem->evaluate(dim, x);
		delete problem;
		break;}
	case 6:
	case 7:{
		Weierstrass * problem = new Weierstrass(config, BASIC);
		result = problem->evaluate(dim, x);
		delete problem;
		break;}
	case 8:
	case 9:{
		Griewank * problem = new Griewank(config, BASIC);
		result = problem->evaluate(dim, x);
		delete problem;
		break;}
	default:
		cout<<"func_no is out of range."<<endl;
		exit(-1);
	}

	return (result);
}


long double HC_3::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case ROTATED:
	case ROTATED_WITH_HIGH_CONDITION_NUMBER_MATRIX:
		result = hybridCompositionFunction3(dim, x);
		break;
	case NONCONTINUOUS_ROTATED:
		result = noncontinuousHybridCompositionFunction3(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;
}



long double HC_3::hybridCompositionFunction3(int dim, const double* x){
	long double result = 0.0;

	result = hybrid_composition(dim, x);

	return result;
}

long double HC_3::noncontinuousHybridCompositionFunction3(int dim, const double* x){
	long double result = 0.0;
	double xx[dim];
	for (int i = 0 ; i < dim ; i++) {
		xx[i] = xRound(x[i], shift_vector2D[0][i]);
	}

	result = hybrid_composition(dim, xx);

	return result;

}

