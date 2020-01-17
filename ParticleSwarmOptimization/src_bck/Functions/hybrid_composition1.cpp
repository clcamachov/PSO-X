#include "hybrid_composition1.h"
#include "../rng.h"
#include "../problem.h"
#include "../utils.h"

#include "rastrigin.h"
#include "griewank.h"
#include "ackley.h"
#include "sphere.h"
#include "weierstrass.h"

#include <string>
#include <sstream>
#include <cmath>
#include <cfloat>

using namespace std;

double HC_1::sigma[NUM_OF_FUNC] = {
		1.0,	1.0,	1.0,	1.0,	1.0,
		1.0,	1.0,	1.0,	1.0,	1.0
};

double HC_1::lambda[NUM_OF_FUNC] = {
		1.0,	1.0,	10.0,	10.0,	5.0/60.0,
		5.0/60.0,	5.0/32.0,	5.0/32.0,	5.0/100.0,	5.0/100.0
};

double HC_1::biases[NUM_OF_FUNC] = {
		0.0,	100.0,	200.0,	300.0,	400.0,
		500.0,	600.0,	700.0,	800.0,	900.0
};

void HC_1::loadIdentityMatrix(double *** M, int w, int h, int d){

	// Generate identity matrices
	for (int i = 0 ; i < h ; i ++) {
		for (int j = 0 ; j < w ; j ++) {
			for (int k = 0 ; k < d ; k ++) {
				M[i][j][k] = 0.0;
			}
		}
		for (int j = 0 ; j < w ; j ++) {
			M[i][j][j] = 1.0;
		}
	}
}

HC_1::HC_1(Configuration* config, int variantID):Problem(config, variantID){

	C = 2000.0;
	num_func = NUM_OF_FUNC;

	zM = allocateMemory2D(dimension, NUM_OF_FUNC);
	z2d = allocateMemory2D(dimension, NUM_OF_FUNC);

	rotation_matrix3D = allocateMemory3D( dimension,NUM_OF_FUNC,dimension);
	shift_vector2D = allocateMemory2D(dimension, NUM_OF_FUNC);

	string file_data = "supportData/hybrid_func1_data.txt";

	Utils::loadMatrixFromFile(file_data, NUM_OF_FUNC, dimension, shift_vector2D);

	w = new double[NUM_OF_FUNC];

	if(variantID == ROTATED || vID == NOISE_ROTATED){
		stringstream dim_name;
		dim_name << dimension;
		string file_m = "supportData/hybrid_func1_M_D" + dim_name.str() + ".txt";
		Utils::loadMatrixFromFile(file_m, NUM_OF_FUNC, dimension, dimension, rotation_matrix3D);
	}else if(vID == BASIC)
		loadIdentityMatrix(rotation_matrix3D,dimension,NUM_OF_FUNC, dimension);

	double testPointM[dimension];
	fmax = new double[NUM_OF_FUNC];
	double testP[dimension];

	// Calculate/estimate the fmax for all the functions involved
	for (int i = 0 ; i < NUM_OF_FUNC ; i ++) {
		for (int j = 0 ; j < dimension ; j ++) {
			testP[j] = (5.0 / lambda[i]);
		}
		rotate(dimension, testPointM, testP, rotation_matrix3D[i]);
		fmax[i] = fabs(basic_func(i, dimension, testPointM));
	}

}


HC_1::~HC_1(){
	deallocateMemory3D(rotation_matrix3D, NUM_OF_FUNC, dimension);
	deallocateMemory2D(zM, NUM_OF_FUNC);
	deallocateMemory2D(z2d, NUM_OF_FUNC);
	deallocateMemory2D(shift_vector2D, NUM_OF_FUNC);
	delete [] fmax;
	delete [] w;

}

long double HC_1::hybrid_composition(int dim, const double* x) {

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

//void HC_1::printShift_vector2D(){
//	//cout << "Optimum solution:" << endl;
//	cout << "[ " ;
//	for(long int i=0; i<config->getProblemDimension(); i++){
//		cout << scientific << shift_vector2D[0][i] << "  ";
//	}
//	cout << " ]" << endl;
//	//shift_vector2D
//}

long double HC_1::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
	case ROTATED:
		result = hybridCompositionFunction1(dim, x);
		break;
	case NOISE_ROTATED:
		result = noiseHybridCompositionFunction1(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;
}


long double HC_1::basic_func(int func_no, int dim, const double* x){
	long double result = 0.0;

	switch(func_no) {
	case 0:
	case 1:{
		Rastrigin * problem = new Rastrigin(config, BASIC);
		result = problem->evaluate(dim, x);
		delete problem;
		break;}
	case 2:
	case 3:{
		Weierstrass * problem = new Weierstrass(config, BASIC);
		result = problem->evaluate(dim, x);
		delete problem;
		break;}
	case 4:
	case 5:{
		Griewank * problem = new Griewank(config, BASIC);
		result = problem->evaluate(dim, x);
		delete problem;
		break;}
	case 6:
	case 7:{
		Ackley * problem = new Ackley(config, BASIC);
		result = problem->evaluate(dim, x);
		delete problem;
		break;}
	case 8:
	case 9:{
		Sphere * problem = new Sphere(config, BASIC);
		result = problem->evaluate(dim, x);
		delete problem;
		break;}
	default:{
		cout<<"func_no is out of range."<<endl;
		exit(-1);
	}
	}

	return (result);
}


long double HC_1::hybridCompositionFunction1(int dim, const double* x){
	long double result = 0.0;

	result = hybrid_composition(dim, x);

	return result;
}

long double HC_1::noiseHybridCompositionFunction1(int dim, const double * x){
	long double result = 0.0;

	result = hybrid_composition(dim, x);

	result *= (1.0 + 0.2 * fabs(RNG::randGauss(1)));

	return result;
}
