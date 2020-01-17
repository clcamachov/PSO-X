#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "rosenbrock.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;


Rosenbrock::Rosenbrock(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == SHIFTED){
		shift_vector = new double[dimension];
		z = new double[dimension];
		if(config->getCompetitionID() == CEC05){
			string file_data = "supportData/rosenbrock_CEC05_func_data.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

		}else{
			string file_data = "supportData/rosenbrock_CEC08_func_data.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		}

		for(int i = 0; i < dimension; i++){
			shift_vector[i] = shift_vector[i] - 1.0;
		}
	}
	else if (config->getCompetitionID() == MIXTURE || config->getCompetitionID() == CEC14){
		if(vID == SHIFTED_ROTATED){

			shift_vector = new double[dimension];
			z = new double[dimension];
			string file_data = "supportData/input_data/shift_data_4.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

			rotation_matrix  = allocateMemory2D(dimension, dimension);
			m_r = new double[dimension];
			stringstream dim_name;
			dim_name << dimension;
			string file_m = "supportData/input_data/M_4_D" + dim_name.str() + ".txt";
			Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);
		}
	}

}


Rosenbrock::~Rosenbrock(){
	if(vID == SHIFTED){
		delete [] shift_vector;
		delete [] z;
	}
}


long double Rosenbrock::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = rosenbrockFunction(dim, x);
		break;
	case SHIFTED:
		result = shiftedRosenbrockFunction(dim, x);
		break;
	case SHIFTED_ROTATED:
		result = shiftedRotatedRosenbrockFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}

long double Rosenbrock::rosenbrockFunction(int dim, const double * x){

	long double F = 0.0;

	for(int i=0;i<dim-1;i++){
		F = F + 100*( pow((pow(x[i],2)-x[i+1]) , 2) ) + pow((x[i]-1) , 2);

	}
	return F;

}

long double Rosenbrock::shiftedRosenbrockFunction(int dim, const double* x){
	long double F = 0;

	shift(dim, z, x, shift_vector);

	F = rosenbrockFunction(dim, z);

	return F;
}

long double Rosenbrock::shiftedRotatedRosenbrockFunction(int dim, const double* x){
	long double F = 0;

	shift(dim, z, x, shift_vector);
	if (config->getCompetitionID() == MIXTURE || config->getCompetitionID() == CEC14){
		double shiftRate = 2.048/100.0;
		for(int i = 0; i < dim; i++)
			z[i]*=shiftRate;
	}

	rotate(dim, m_r, z, rotation_matrix);

	F = rosenbrockFunction(dim, m_r);

	return F;
}
