#include "../problem.h"
#include "elliptic.h"
#include "../utils.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;


Elliptic::Elliptic(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == SHIFTED_ROTATED_HIGH_CONDITIONED){
		if(config->getCompetitionID() == CEC05){
			string file_data = "supportData/high_cond_elliptic_rot_data.txt";
			stringstream dim_name;
			dim_name << dimension;
			string file_m; //= "supportData/elliptic_M_D" + dim_name.str() + ".txt";

			if (dimension > 2 && dimension < 10)
				file_m = "supportData/elliptic_M_D10.txt";
			else if (dimension > 10 && dimension < 30)
				file_m = "supportData/elliptic_M_D30.txt";
			else if (dimension > 30 && dimension < 50)
				file_m = "supportData/elliptic_M_D50.txt";
			else
				file_m = "supportData/elliptic_M_D" + dim_name.str() + ".txt";

			rotation_matrix  = allocateMemory2D(dimension, dimension);
			shift_vector = new double[dimension];
			z = new double[dimension];
			m_r = new double[dimension];

			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
			Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);
		}

		if(config->getCompetitionID() == MIXTURE || config->getCompetitionID() == CEC14){
			string file_data = "supportData/input_data/shift_data_1.txt";
			stringstream dim_name;
			dim_name << dimension;
			string file_m; //= "supportData/input_data/M_1_D" + dim_name.str() + ".txt";

			if (dimension > 2 && dimension < 10)
				file_m = "supportData/input_data/M_1_D10.txt";
			else if (dimension > 10 && dimension < 20)
				file_m = "supportData/input_data/M_1_D20.txt";
			else if (dimension > 20 && dimension < 30)
				file_m = "supportData/input_data/M_1_D30.txt";
			else if (dimension > 30 && dimension < 50)
				file_m = "supportData/input_data/M_1_D50.txt";
			else if (dimension > 50 && dimension < 100)
				file_m = "supportData/input_data/M_1_D100.txt";
			else
				file_m = "supportData/input_data/M_1_D" + dim_name.str() + ".txt";

			rotation_matrix  = allocateMemory2D(dimension, dimension);
			shift_vector = new double[dimension];
			z = new double[dimension];
			m_r = new double[dimension];

			Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
			Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);

		}

	}	

}

Elliptic::~Elliptic(){
	if(vID == SHIFTED_ROTATED_HIGH_CONDITIONED){
		deallocateMemory2D(rotation_matrix, dimension);
		delete [] shift_vector;
		delete [] z;
		delete [] m_r;
	}
}

long double Elliptic::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = ellipticFunction(dim, x);
		break;
	case SHIFTED_ROTATED_HIGH_CONDITIONED:
	case MIXTURE:
		result = shiftedRotatedHighConditionedEllipticFunction(dim, x);
		break;
	case SHIFTED_ROTATED_HIGH_CONDITIONED_CEC2013:
		result = shiftedRotatedHighConditionedEllipticFunction2(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;
}


long double Elliptic::ellipticFunction(int dim, const double* x){
	long double sum = 0.0;
	double a = 1e6;

	for (int i = 0 ; i < dim ; i ++) {
		sum += pow(a, (((double )i)/((double )(dim-1)))) * x[i] * x[i];
	}

	return sum;

}

long double Elliptic::shiftedRotatedHighConditionedEllipticFunction(int dim, const double* x){

	shift(dim, z, x, shift_vector);
	rotate(dim, m_r, z, rotation_matrix);

	long double sum = 0.0;
	const double constant = pow(1.0e6, 1.0/((double)dim-1.0));

	for (int i = 0 ; i < dim ; i++) {

		sum += pow(constant, (double)i) * m_r[i] * m_r[i];
	}

	return sum;
}

long double Elliptic::shiftedRotatedHighConditionedEllipticFunction2(int dim, const double* x){
	//cout<<"yesf"<<endl;
	shift(dim, z, x, shift_vector);
	rotate(dim, m_r, z, rotation_matrix);
	oszfunc(m_r,z,dim);

	long double sum = 0.0;
	//const double constant = pow(1.0e6, 1.0/((double)dim-1.0));

	for (int i = 0 ; i < dim ; i++) {


		sum += pow(10.0, 6.0*i/(dim-1.0)) * z[i] * z[i];
	}

	return sum;
}

