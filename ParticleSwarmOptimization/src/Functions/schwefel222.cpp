#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "schwefel222.h"

#include <float.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

Schwefel222::Schwefel222(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == SHIFTED){
		shift_vector = new double[dimension];
		string file_data = "../supportData/schwefel222_CEC08_func_data.txt";
		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

		z = new double[dimension];
	}

}

Schwefel222::~Schwefel222(){

	if(vID == SHIFTED){
		delete [] z;
		delete [] shift_vector;
	}

}

long double Schwefel222::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = schwefel222Function(dim, x);
		break;
	case SHIFTED:
		result = shiftedSchwefel222Function(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}


long double Schwefel222::shiftedSchwefel222Function(int dim, const double* x){

	long double F = 0.0;

	shift(dim, z, x, shift_vector);

	F = schwefel222Function(dim, z);

	return F;

}

long double Schwefel222::schwefel222Function(int dim, const double* x){

	long double sum, currentGen, prod;

	sum = 0.0;
	prod = 1.0;

	for (int i = 0; i < dim; i++)
	{
		currentGen = fabs(x[i]);
		sum += currentGen;
		prod *= currentGen;
	}

	return sum + prod;
}
