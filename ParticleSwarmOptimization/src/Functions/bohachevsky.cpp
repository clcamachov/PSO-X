#include "../rng.h"
#include "../problem.h"
#include "bohachevsky.h"
#include "../utils.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

Bohachevsky::Bohachevsky(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == SHIFTED){
		shift_vector = new double[dimension];
		string file_data = "supportData/bohachevsky_CEC08_func_data.txt";
		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

		z = new double[dimension];
	}

}

Bohachevsky::~Bohachevsky(){
	if(vID == SHIFTED){
		delete [] z;
		delete shift_vector;
	}
}

long double Bohachevsky::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = bohachevskyFunction(dim, x);
		break;
	case SHIFTED:
		result = shiftedBohachevskyFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}


long double Bohachevsky::bohachevskyFunction(int dim, const double* x){
	long double sum = 0.0;
	int i;
	double currentGen;
	double nextGen;

	currentGen = x[0];

	for (i = 1; i < dim; i++)
	{
		nextGen = x[i];
		sum += currentGen * currentGen + 2.0 * nextGen * nextGen;
		sum += -0.3 * cos(3.0 * PI * currentGen) -0.4 * cos(4.0 * PI * nextGen) + 0.7;
		currentGen = nextGen;
	}

	return sum;
}



long double Bohachevsky::shiftedBohachevskyFunction(int dim, const double* x){
	long double F = 0.0;

	shift(dim, z, x, shift_vector);
	F = bohachevskyFunction(dim, z);

	return F;

}
