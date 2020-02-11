#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "schwefel26.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cfloat>

using namespace std;


Schwefel26::Schwefel26(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	double ** data_matrix  = allocateMemory2D(dimension, dimension+1);
	v_B = new double[dimension]; 
	m_A = allocateMemory2D(dimension, dimension);
	z = new double[dimension];

	shift_vector = new double[dimension]; 
	string file_data = "../supportData/schwefel_206_data.txt";
	Utils::loadMatrixFromFile(file_data, dimension+1, dimension, data_matrix);


	for (int i = 0 ; i < dimension ; i ++) {
		if ((i+1) <= ceil(dimension / 4.0))
			shift_vector[i] = -100.0;
		else if ((i+1) >= floor((3.0 * dimension) / 4.0))
			shift_vector[i] = 100.0;
		else
			shift_vector[i] = data_matrix[0][i];
	}

	for (int i = 0 ; i < dimension ; i ++) {
		for (int j = 0 ; j < dimension ; j ++) {
			m_A[i][j] = data_matrix[i+1][j];
		}
	}

	multiplyMatrix(v_B, m_A, shift_vector);

	deallocateMemory2D(data_matrix, dimension);

}

Schwefel26::~Schwefel26(){
	delete [] v_B;
	deallocateMemory2D(m_A, dimension);
	delete [] shift_vector;
	delete [] z;
}

long double Schwefel26::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC_GLOBAL_OPTIMUM_ON_BOUNDS:
		result = schwefel26GlobalOptimumOnBoundsFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}


long double Schwefel26::schwefel26GlobalOptimumOnBoundsFunction(int dim, const double* x){
	long double max = LDBL_MIN;

	multiplyMatrix(z, m_A, x);

	for (int i = 0 ; i < dim ; i ++) {
		long double temp = fabs(z[i] - v_B[i]);
		if (max < temp)
			max = temp;
	}

	return max;
}
