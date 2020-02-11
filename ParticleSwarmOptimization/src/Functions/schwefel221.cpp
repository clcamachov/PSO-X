#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "schwefel221.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

Schwefel221::Schwefel221(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == SHIFTED){
		shift_vector = new double[dimension];
		string file_data = "../supportData/schwefel221_CEC08_func_data.txt";
		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

		z = new double[dimension];
	}
}

Schwefel221::~Schwefel221(){  
	if(vID == SHIFTED){
		delete [] z;
		delete [] shift_vector;
	}

}

long double Schwefel221::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = schwefel221Function(dim, x);
		break;
	case SHIFTED:
		result = shiftedSchwefel221Function(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}

long double Schwefel221::shiftedSchwefel221Function(int dim, const double* x){
	long double F = 0.0;

	shift(dim, z, x, shift_vector);

	F = schwefel221Function(dim, z);

	return F;
}

long double Schwefel221::schwefel221Function(int dim, const double* x){
	int i;

	long double F = abss(x[0]);
	for(i=1;i<dim;i++){
		F = max((double)F ,abss(x[i]));
	}
	return F;
}
