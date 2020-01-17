#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "schafferf6.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;


SchafferF6::SchafferF6(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == SHIFTED_ROTATED_EXPANDED){
		if(config->getCompetitionID() == CEC05){
			rotation_matrix  = allocateMemory2D(dimension, dimension);
			shift_vector = new double[dimension];
			z = new double[dimension];
			m_r = new double[dimension];
			string file_data = "supportData/E_ScafferF6_func_data.txt";
			stringstream dim_name;
			dim_name << dimension;
			string file_m = "supportData/E_ScafferF6_M_D" + dim_name.str() + ".txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
			Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);
		}
		if(config->getCompetitionID() == MIXTURE){
			string file_data = "supportData/input_data/shift_data_16.txt";
			stringstream dim_name;
			dim_name << dimension;
			string file_m = "supportData/input_data/M_16_D" + dim_name.str() + ".txt";

			rotation_matrix  = allocateMemory2D(dimension, dimension);
			shift_vector = new double[dimension];
			z = new double[dimension];
			m_r = new double[dimension];

			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
			Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);


		}
	}
}

SchafferF6::~SchafferF6(){
	if(vID == SHIFTED_ROTATED_EXPANDED){
		deallocateMemory2D(rotation_matrix,dimension);
		delete [] shift_vector;
		delete [] z;
		delete [] m_r;
	}
}


long double SchafferF6::evaluate(int dim, const double* x) {
	long double result;
	switch(vID){
	case EXPANDED:
		result = expandedSchafferF6Function(dim, x);
		break;
	case SHIFTED_ROTATED_EXPANDED:
		result = shiftedRotatedExpandedSchafferF6Function(dim, x);
		break;
	case NONCONTINUOUS_EXPANDED:
		result = noncontinuousExpandedSchafferF6Function(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}

double SchafferF6::schafferF6Function(double x, double y){
	double temp1 = x*x + y*y;
	double temp2 = sin(sqrt(temp1));
	double temp3 = 1.0 + 0.001 * temp1;
	return (0.5 + ((temp2 * temp2 - 0.5)/(temp3 * temp3)));

}

long double SchafferF6::expandedSchafferF6Function(int dim, const double* x){
	long double sum = 0.0;

	for (int i = 1 ; i < dim ; i ++) {
		sum += schafferF6Function(x[i-1], x[i]);
	}

	sum += schafferF6Function(x[dim-1], x[0]);

	return (sum);

}

long double SchafferF6::shiftedRotatedExpandedSchafferF6Function(int dim, const double* x){
	shift(dim,z, x,shift_vector);
	rotate(dim, m_r, z, rotation_matrix);
	return expandedSchafferF6Function(dim, m_r);
}


long double SchafferF6::noncontinuousExpandedSchafferF6Function(int dim, const double* x){
	long double sum = 0.0;

	for (int i = 1 ; i < dim ; i ++) {
		sum += schafferF6Function(xRound(x[i-1]), xRound(x[i]));
	}

	sum += schafferF6Function(xRound(x[dim-1]), xRound(x[0]));

	return (sum);

}

