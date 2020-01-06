#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "schaffer.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

Schaffer::Schaffer(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == SHIFTED){
		shift_vector = new double[dimension];
		z = new double[dimension];
		string file_data = "supportData/schaffer_CEC08_func_data.txt";
		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

	}

}

Schaffer::~Schaffer(){
	if(vID == SHIFTED){
		delete [] shift_vector;
		delete [] z;
	}
}


long double Schaffer::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = schafferFunction(dim, x);
		break;
	case SHIFTED:
		result = shiftedSchafferFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}


long double Schaffer::schafferFunction(int dim, const double* x){
	int i;
	long double sum, aux, aux2;
	double currentGen, nextGen;

	sum = 0.0;
	currentGen = x[0];
	currentGen = currentGen * currentGen;

	for (i = 1; i < dim; i++)
	{
		nextGen = x[i];
		nextGen = nextGen * nextGen;
		aux = currentGen + nextGen;
		currentGen = nextGen;
		aux2 = sin(50. * pow((double)aux, 0.1));
		sum += pow((double)aux, 0.25) * (aux2 * aux2 + 1.0);
	}

	return sum;
}

long double Schaffer::shiftedSchafferFunction(int dim, const double* x){
	shift(dim, z, x, shift_vector);
	return schafferFunction(dim, z);
}
