#include "../problem.h"
#include "../utils.h"
#include "hgbat.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;


Hgbat::Hgbat(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == MIXTURE || vID == SHIFTED_ROTATED){
		string file_data = "supportData/input_data/shift_data_14.txt";
		stringstream dim_name;
		dim_name << dimension;
		string file_m = "supportData/input_data/M_14_D" + dim_name.str() + ".txt";

		rotation_matrix  = allocateMemory2D(dimension, dimension);
		shift_vector = new double[dimension];
		z = new double[dimension];
		m_r = new double[dimension];

		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);

	}


}

Hgbat::~Hgbat(){
	if(vID == SHIFTED_ROTATED || vID == MIXTURE){
		deallocateMemory2D(rotation_matrix, dimension);
		delete [] shift_vector;
		delete [] z;
		delete [] m_r;
	}
}

long double Hgbat::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = hgbatFunction(dim, x);
		break;
	case SHIFTED_ROTATED:
	case MIXTURE:
		result = shiftedRotatedHgbatFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;
}


long double Hgbat::hgbatFunction(int dim, const double* x){

	int i;
	double alpha,r2,sum_z;
	alpha=1.0/4.0;

	r2 = 0.0;
	sum_z=0.0;
	for (i=0; i<dim; i++)
	{
		//x[i]=x[i]-1.0;//shift to orgin
		double val = x[i]-1.0;
		r2 += val*val;
		sum_z += val;
	}

	return pow(fabs(pow(r2,2.0)-pow(sum_z,2.0)),2*alpha) + (0.5*r2 + sum_z)/dim + 0.5;

}

long double Hgbat::shiftedRotatedHgbatFunction(int dim, const double* x){

	shift(dim, z, x, shift_vector);

	if (config->getCompetitionID() == MIXTURE || config->getCompetitionID() == CEC14){
		double shiftRate = 5.0/100.0;
		for(int i = 0; i < dim; i++)
			z[i]*=shiftRate;
	}

	rotate(dim, m_r, z, rotation_matrix);
	return 	hgbatFunction(dim, m_r);
}
