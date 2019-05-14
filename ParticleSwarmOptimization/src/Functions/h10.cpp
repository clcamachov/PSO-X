#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "schwefel222.h"
#include "bohachevsky.h"
#include "h10.h"

#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

Hybrid10::Hybrid10(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;
	problem1 = new Bohachevsky(config, BASIC);
	problem2 = new Schwefel222(config, BASIC);


	if(vID == SHIFTED){
		shift_vector = new double[dimension];
		string file_data = "supportData/hybrid10_CEC08_func_data.txt";
		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

		z = new double[dimension];
	}

}

Hybrid10::~Hybrid10(){
	delete problem1;
	delete problem2;
	if(vID == SHIFTED){
		delete [] z;
		delete [] shift_vector;
	}


}

long double Hybrid10::evaluate(int dim, const double* x) {
	long double result = 0.0;

	switch(vID){
	case BASIC:
		result = h10Function(dim, x);
		break;
	case SHIFTED:
		result = shiftedH10Function(dim,x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}

long double Hybrid10::h10Function(int dim, const double* x){
	double part1[dim], part2[dim];
	int size1, size2;
	long double f1, f2;

	divideFunctions(x, part1, part2, 0.75, &size1, &size2);

	f1 = problem1->evaluate(size1,part1);
	f2 = problem2->evaluate(size2,part2);

	return f1+f2;
}


long double Hybrid10::shiftedH10Function(int dim, const double* x){

	shift(dim, z, x, shift_vector);

	return h10Function(dim, z);
}

