#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "schwefel12.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

Schwefel12::Schwefel12(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;
	if(vID == SHIFTED || vID == NOISE_SHIFTED){ 
		shift_vector = new double[dimension];
		z = new double[dimension];
		if(config->getCompetitionID() == CEC05){
			string file_data = "../supportData/schwefel_102_CEC05_func_data.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		}else{
			string file_data = "../supportData/schwefel_102_CEC08_func_data.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		}
	}	

}

Schwefel12::~Schwefel12(){
	if(vID == SHIFTED || vID == NOISE_SHIFTED){ 
		delete [] shift_vector;
		delete [] z;
	}
}

long double Schwefel12::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = schwefel12Function(dim, x);
		break;
	case SHIFTED:
		result = shiftedSchwefel12Function(dim, x);
		break;
	case NOISE_SHIFTED:
		result = noiseShiftedSchwefel12Function(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}


long double Schwefel12::schwefel12Function(int dim, const double* x){
	long double Sum=0.0, Val=0.0;

	for (int i = 0; i < dim; i++){
		Val += x[i];
		Sum += Val * Val;
	}

	return Sum;
}


long double Schwefel12::shiftedSchwefel12Function(int dim, const double* x){
	shift(dim, z, x, shift_vector);
	long double F = schwefel12Function(dim, z);
	return F;
}

long double Schwefel12::noiseShiftedSchwefel12Function(int dim, const double* x){
	shift(dimension, z, x, shift_vector);
	long double F = schwefel12Function(dim, z);
	F *= (1.0 + 0.4 * fabs(RNG::randGauss(1)));
	return F;
}
