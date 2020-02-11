#include "../problem.h"
#include "discus.h"
#include "../utils.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;


Discus::Discus(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == MIXTURE || vID == SHIFTED_ROTATED){
		string file_data = "../supportData/input_data/shift_data_2.txt";
		stringstream dim_name;
		dim_name << dimension;
		string file_m; // = "../supportData/input_data/M_2_D" + dim_name.str() + ".txt";

		if (dimension > 2 && dimension < 10)
			file_m = "../supportData/input_data/M_2_D10.txt";
		else if (dimension > 10 && dimension < 20)
			file_m = "../supportData/input_data/M_2_D20.txt";
		else if (dimension > 20 && dimension < 30)
			file_m = "../supportData/input_data/M_2_D30.txt";
		else if (dimension > 30 && dimension < 50)
			file_m = "../supportData/input_data/M_2_D50.txt";
		else if (dimension > 50 && dimension < 100)
			file_m = "../supportData/input_data/M_2_D100.txt";
		else
			file_m = "../supportData/input_data/M_2_D" + dim_name.str() + ".txt";

		rotation_matrix = allocateMemory2D(dimension, dimension);
		shift_vector = new double[dimension];
		z = new double[dimension];
		m_r = new double[dimension];

		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);
	}

}

Discus::~Discus(){
	if(vID == SHIFTED_ROTATED || vID == MIXTURE){
		deallocateMemory2D(rotation_matrix, dimension);
		delete [] shift_vector;
		delete [] z;
		delete [] m_r;
	}
}

long double Discus::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = discusFunction(dim, x);
		break;
	case SHIFTED_ROTATED:
	case MIXTURE:
		result = shiftedRotatedDiscusFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;
}


long double Discus::discusFunction(int dim, const double* x){
	long double sum = 1.0e6*x[0]*x[0];
	for (int i=1; i<dim; i++)
	{
		sum += x[i]*x[i];
	}
	return sum;
}

long double Discus::shiftedRotatedDiscusFunction(int dim, const double* x){

	shift(dim, z, x, shift_vector);
	rotate(dim, m_r, z, rotation_matrix);

	long double sum = 1.0e6*m_r[0]*m_r[0];

	for (int i=1; i<dim; i++)
	{
		sum += m_r[i]*m_r[i];
	}

	return sum;
}
