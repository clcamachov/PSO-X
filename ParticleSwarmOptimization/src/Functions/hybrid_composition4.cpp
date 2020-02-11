#include "hybrid_composition4.h"
#include "../rng.h"
#include "../problem.h"
#include "../utils.h"

#include "rastrigin.h"
#include "schafferf6.h"
#include "griewank.h"
#include "griewank_rosenbrock.h"
#include "weierstrass.h"
#include "ackley.h"
#include "elliptic.h"
#include "sphere.h"

#include <cmath>
#include <string>
#include <sstream>
#include <cfloat>

using namespace std;



double HC_4::sigma[NUM_OF_FUNC] = {
		2.0,	2.0,	2.0,	2.0,	2.0,
		2.0,	2.0,	2.0,	2.0,	2.0
};

double HC_4::lambda[NUM_OF_FUNC] = {
		10.0,		5.0/20.0,	1.0,	5.0/32.0,	1.0,
		5.0/100.0,	5.0/50.0,	1.0,	5.0/100.0,	5.0/100.0

};

double HC_4::biases[NUM_OF_FUNC] = {
		0.0,	100.0,	200.0,	300.0,	400.0,
		500.0,	600.0,	700.0,	800.0,	900.0
};


HC_4::HC_4(Configuration* config, int variantID):Problem(config,variantID){



	C = 2000.0;
	num_func = NUM_OF_FUNC;

	zM = allocateMemory2D(dimension,NUM_OF_FUNC);
	z2d = allocateMemory2D(dimension,NUM_OF_FUNC);
	rotation_matrix3D = allocateMemory3D(dimension,NUM_OF_FUNC,dimension);
	shift_vector2D = allocateMemory2D(dimension,NUM_OF_FUNC);

	string file_data = "../supportData/hybrid_func4_data.txt";

	stringstream dim_name;
	dim_name << dimension;
	string file_m; // = "../supportData/hybrid_func4_M_D" + dim_name.str() + ".txt";

	if (dimension > 2 && dimension < 10)
		file_m = "../supportData/hybrid_func4_M_D10.txt";
	else if (dimension > 10 && dimension < 30)
		file_m = "../supportData/hybrid_func4_M_D30.txt";
	else if (dimension > 30 && dimension < 50)
		file_m = "../supportData/hybrid_func4_M_D50.txt";
	else
		file_m = "../supportData/hybrid_func4_M_D" + dim_name.str() + ".txt";

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

	if( vID == ROTATED_WITHOUT_BOUNDS)
		hasBounds = false;
}


HC_4::~HC_4(){
	deallocateMemory3D(rotation_matrix3D, NUM_OF_FUNC, dimension);
	deallocateMemory2D(zM, NUM_OF_FUNC);
	deallocateMemory2D(z2d, NUM_OF_FUNC);
	deallocateMemory2D(shift_vector2D, NUM_OF_FUNC);
	delete [] fmax;
	delete [] w;
}


long double HC_4::basic_func(int func_no, int dim, const double* x){
	long double result = 0.0;

	switch(func_no) {
	case 0:
	{Weierstrass * problem = new Weierstrass(config, BASIC);
	result = problem->evaluate(dim, x);
	delete problem;
	break;}
	case 1:
	{SchafferF6 * problem = new SchafferF6(config, EXPANDED);
	result = problem->evaluate(dim, x);
	delete problem;
	break;}
	case 2:
	{	GriewankRosenbrock * problem = new GriewankRosenbrock(config, BASIC);
	result = problem->evaluate(dim, x);
	delete problem;
	break;}
	case 3:
	{	Ackley * problem = new Ackley(config, BASIC);
	result = problem->evaluate(dim,x);
	delete problem;
	break;}
	case 4:
	{	Rastrigin * problem = new Rastrigin(config, BASIC);
	result = problem->evaluate(dim, x);
	delete problem;
	break;}
	case 5:
	{	Griewank * problem = new Griewank(config, BASIC);
	result = problem->evaluate(dim, x);
	delete problem;
	break;}
	case 6:
	{	SchafferF6 * problem = new SchafferF6(config, NONCONTINUOUS_EXPANDED);
	result = problem->evaluate(dim, x);
	delete problem;
	break;}
	case 7:
	{	Rastrigin * problem = new Rastrigin(config, NONCONTINUOUS_BASIC);
	result = problem->evaluate(dim, x);
	delete problem;
	break;}
	case 8:
	{	Elliptic * problem = new Elliptic(config, BASIC);
	result = problem->evaluate(dim, x);
	delete problem;
	break;}
	case 9:
	{	Sphere * problem = new Sphere(config, NOISE_BASIC);
	result = problem->evaluate(dim, x);
	delete problem;
	break;}
	default:
		cout<<"func_no is out of range."<<endl;
		exit(-1);
	}

	return (result);
}


long double HC_4::hybrid_composition(int dim, const double* x) {

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


long double HC_4::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case ROTATED:
	case ROTATED_WITHOUT_BOUNDS:
		result = hybridCompositionFunction4(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;
}



long double HC_4::hybridCompositionFunction4(int dim, const double* x){
	long double result = 0.0;

	result = hybrid_composition(dim, x);

	return result;
}
