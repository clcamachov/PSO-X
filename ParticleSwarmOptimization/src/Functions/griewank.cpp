#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "griewank.h"

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>

using namespace std;


Griewank::Griewank(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(config->getCompetitionID() == CEC05){
		if(vID == SHIFTED || vID == SHIFTED_ROTATED || vID == ROTATED_WITHOUT_BOUNDS){
			shift_vector = new double[dimension];
			z = new double[dimension];
			string file_data = "supportData/griewank_CEC05_func_data.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		}
		if (vID == SHIFTED_ROTATED || vID == ROTATED_WITHOUT_BOUNDS){
			rotation_matrix  = allocateMemory2D(dimension, dimension);
			m_r = new double[dimension];
			stringstream dim_name;
			dim_name<< dimension;
			string file_m; // = "supportData/griewank_M_D" + dim_name.str() + ".txt";

			if (dimension > 2 && dimension < 10)
				file_m = "supportData/griewank_M_D10.txt";
			else if (dimension > 10 && dimension < 30)
				file_m = "supportData/griewank_M_D30.txt";
			else if (dimension > 30 && dimension < 50)
				file_m = "supportData/griewank_M_D50.txt";
			else
				file_m = "supportData/griewank_M_D" + dim_name.str() + ".txt";

			Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);
			if( vID == ROTATED_WITHOUT_BOUNDS)
				hasBounds = false;
		}
	}else if (config->getCompetitionID() == MIXTURE || config->getCompetitionID() == CEC14){
		if(vID == SHIFTED_ROTATED){

			shift_vector = new double[dimension];
			z = new double[dimension];
			string file_data = "supportData/input_data/shift_data_7.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

			rotation_matrix  = allocateMemory2D(dimension, dimension);
			m_r = new double[dimension];
			stringstream dim_name;
			dim_name << dimension;
			string file_m; // = "supportData/input_data/M_7_D" + dim_name.str() + ".txt";

			if (dimension > 2 && dimension < 10)
				file_m = "supportData/input_data/M_7_D10.txt";
			else if (dimension > 10 && dimension < 20)
				file_m = "supportData/input_data/M_7_D20.txt";
			else if (dimension > 20 && dimension < 30)
				file_m = "supportData/input_data/M_7_D30.txt";
			else if (dimension > 30 && dimension < 50)
				file_m = "supportData/input_data/M_7_D50.txt";
			else if (dimension > 50 && dimension < 100)
				file_m = "supportData/input_data/M_7_D100.txt";
			else
				file_m = "supportData/input_data/M_7_D" + dim_name.str() + ".txt";

			Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);
		}
		if (vID == SHIFTED){

			shift_vector = new double[dimension];
			z = new double[dimension];
			string file_data = "supportData/griewank_CEC08_func_data.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		}
	}
	else {
		if (vID == SHIFTED){

			shift_vector = new double[dimension];
			z = new double[dimension];
			string file_data = "supportData/griewank_CEC08_func_data.txt";
			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		}
	}

}

Griewank::~Griewank(){
	if(vID == SHIFTED || vID == SHIFTED_ROTATED || vID == ROTATED_WITHOUT_BOUNDS){
		delete [] shift_vector;
		delete [] z;
	}

	if(vID == SHIFTED_ROTATED || vID == ROTATED_WITHOUT_BOUNDS){
		delete [] m_r;
		deallocateMemory2D(rotation_matrix, dimension);
	}
}

long double Griewank::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = griewankFunction(dim, x);
		break;
	case SHIFTED:
		result = shiftedGriewankFunction(dim, x);
		break;
	case SHIFTED_ROTATED:
	case ROTATED_WITHOUT_BOUNDS:
		result = shiftedRotatedGriewankFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}

long double Griewank::griewankFunction(int dim, const double* x){
	int i;
	long double F1 = 0.0;
	long double F2 = 1.0;

	for(i=0;i<dim;i++){
		F1 = F1 + ( pow(x[i],2.0) / 4000.0 );
		F2 = F2 * ( cos(x[i]/sqrt(i+1)));

	}

	return F1 - F2 + 1.0;
}


long double Griewank::shiftedGriewankFunction(int dim, const double* x){
	shift(dim, z, x, shift_vector);
	return griewankFunction(dim, z);
}

long double Griewank::shiftedRotatedGriewankFunction(int dim, const double* x){
	shift(dim, z, x, shift_vector);
	if (config->getCompetitionID() == MIXTURE || config->getCompetitionID() == CEC14){
		double shiftRate = 600.0/100.0;
		for(int i = 0; i < dim; i++)
			z[i]*=shiftRate;
	}

	rotate(dim, m_r, z, rotation_matrix);
	return griewankFunction(dim, m_r);

}
