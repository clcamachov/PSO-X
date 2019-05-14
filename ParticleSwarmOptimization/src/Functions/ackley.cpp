/*
 * ackley.cpp
 *
 *  Created on: Jun 6, 2018
 *      Author: leonardo
 */

#include "../rng.h"
#include "../problem.h"
#include "ackley.h"
#include "../utils.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

Ackley::Ackley(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(config->getCompetitionID() == CEC05 ){
		if(vID == SHIFTED || vID == SHIFTED_ROTATED
				|| vID == SHIFTED_ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS){
			shift_vector = new double[dimension];
			z = new double[dimension];
			string file_data = "supportData/ackley_CEC05_func_data.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		}

		if(vID == SHIFTED_ROTATED || vID == SHIFTED_ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS){
			rotation_matrix  = allocateMemory2D(dimension, dimension);
			m_r = new double[dimension];
			stringstream dim_name;
			dim_name << dimension;
			string file_m = "supportData/ackley_M_D" + dim_name.str() + ".txt";
			Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);
		}
	}else if(config->getCompetitionID() == CEC14 || config->getCompetitionID() == SOFT_COMPUTING
			|| config->getCompetitionID() == MIXTURE){
		if(vID == SHIFTED){
			shift_vector = new double[dimension];
			z = new double[dimension];
			string file_data = "supportData/ackley_CEC08_func_data.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
			//27.00775710883052500000 -16.13168918554272900000 6.41055501250719570000
		}
		else if(vID == SHIFTED_ROTATED){

			shift_vector = new double[dimension];
			z = new double[dimension];
			string file_data = "supportData/input_data/shift_data_5.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

			rotation_matrix  = allocateMemory2D(dimension, dimension);
			m_r = new double[dimension];
			stringstream dim_name;
			dim_name << dimension;
			string file_m = "supportData/input_data/ackley_M_D" + dim_name.str() + ".txt";
			Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);
		}
	}

	if(vID == SHIFTED_ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS){
		for (int i = 0 ; i < dimension ; i += 2) {
			shift_vector[i] = -32.0;
		}
	}

}

Ackley::~Ackley(){
	if(vID == SHIFTED || vID == SHIFTED_ROTATED || vID == SHIFTED_ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS){
		delete [] shift_vector;
		delete [] z;
	}

	if(vID == SHIFTED_ROTATED || vID == SHIFTED_ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS){
		delete [] m_r;
		deallocateMemory2D(rotation_matrix, dimension);
	}

}

long double Ackley::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = ackleyFunction(dim, x);
		break;
	case SHIFTED:
		result = shiftedAckleyFunction(dim, x);
		break;
	case SHIFTED_ROTATED:
	case SHIFTED_ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS:
		result = shiftedRotatedAckleyFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	}

	return result;

}

long double Ackley::ackleyFunction(int dim, const double* x){
	int i;
	//long double z;
	long double Sum1 = 0;
	long double Sum2 = 0;

	for(i=0;i<dim;i++){
		Sum1 = Sum1 + pow(x[i] , 2 );
		Sum2 = Sum2 + cos(2*PI*x[i]);
	}

	return -20*exp(-0.2*sqrt(Sum1/(double)dimension)) -exp(Sum2/(double)dimension) + 20 + exp(1);
}


long double Ackley::shiftedAckleyFunction(int dim, const double * x){
	shift(dim, z, x, shift_vector);
	return ackleyFunction(dim, z);
}


long double Ackley::shiftedRotatedAckleyFunction(int dim, const double * x){
	shift(dim, z, x, shift_vector);
	rotate(dim, m_r, z, rotation_matrix);
	return ackleyFunction(dim, m_r);
}


