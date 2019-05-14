#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "schwefel213.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

Schwefel213::Schwefel213(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;


	double m_o[dimension];

	double ** data_matrix = allocateMemory2D(dimension,201);
	string file_data = "supportData/schwefel_213_data.txt";
	Utils::loadMatrixFromFile(file_data, 201, dimension, data_matrix);

	m_a = allocateMemory2D(dimension,dimension);
	m_b = allocateMemory2D(dimension,dimension);

	m_A = new double[dimension];
	m_B = new double[dimension];

	for (int i = 0 ; i < dimension ; i ++) {
		for (int j = 0 ; j < dimension ; j ++) {
			m_a[i][j] = data_matrix[i][j];
			m_b[i][j] = data_matrix[100+i][j];
		}
		m_o[i] = data_matrix[100+100][i];
	}

	for (int i = 0 ; i < dimension ; i ++) {
		m_A[i] = 0.0;
		for (int j = 0 ; j < dimension ; j ++) {
			m_A[i] += (m_a[i][j] * sin(m_o[j]) + m_b[i][j] * cos(m_o[j]));
		}
	}

	deallocateMemory2D(data_matrix,dimension);	
}

Schwefel213::~Schwefel213(){
	deallocateMemory2D(m_a,dimension);
	deallocateMemory2D(m_b,dimension);

	delete [] m_A;
	delete [] m_B;
}


long double Schwefel213::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = schwefel213Function(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}


long double Schwefel213::schwefel213Function(int dim, const double* x){
	double sum = 0.0;

	for (int i = 0 ; i < dim ; i ++) {
		m_B[i] = 0.0;
		for (int j = 0 ; j < dim ; j ++) {
			m_B[i] += (m_a[i][j] * sin(x[j]) + m_b[i][j] * cos(x[j]));
		}

		double temp = m_A[i] - m_B[i];
		sum += (temp * temp);
	}

	return sum;

}
