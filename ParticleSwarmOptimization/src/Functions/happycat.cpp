#include "../problem.h"
#include "../utils.h"
#include "happycat.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;


Happycat::Happycat(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == MIXTURE || vID == SHIFTED_ROTATED){
		string file_data = "../supportData/input_data/shift_data_13.txt";
		stringstream dim_name;
		dim_name << dimension;
		string file_m; // = "../supportData/input_data/M_13_D" + dim_name.str() + ".txt";

		if (dimension > 2 && dimension < 10)
			file_m = "../supportData/input_data/M_13_D10.txt";
		else if (dimension > 10 && dimension < 20)
			file_m = "../supportData/input_data/M_13_D20.txt";
		else if (dimension > 20 && dimension < 30)
			file_m = "../supportData/input_data/M_13_D30.txt";
		else if (dimension > 30 && dimension < 50)
			file_m = "../supportData/input_data/M_13_D50.txt";
		else if (dimension > 50 && dimension < 100)
			file_m = "../supportData/input_data/M_13_D100.txt";
		else
			file_m = "../supportData/input_data/M_13_D" + dim_name.str() + ".txt";

		rotation_matrix  = allocateMemory2D(dimension, dimension);
		shift_vector = new double[dimension];
		z = new double[dimension];
		m_r = new double[dimension];

		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);

	}


}

Happycat::~Happycat(){
	if(vID == SHIFTED_ROTATED || vID == MIXTURE){
		deallocateMemory2D(rotation_matrix, dimension);
		delete [] shift_vector;
		delete [] z;
		delete [] m_r;
	}
}

long double Happycat::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = happycatFunction(dim, x);
		break;
	case SHIFTED_ROTATED:
	case MIXTURE:
		result = shiftedRotatedHappycatFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;
}


long double Happycat::happycatFunction(int dim, const double* x){

	int i;
	double alpha,r2,sum_z;
	alpha=1.0/8.0;

	r2 = 0.0;
	sum_z=0.0;
	for (i=0; i<dim; i++)
	{
		double val=x[i]-1.0;//shift to orgin
		r2 += val*val;
		sum_z += val;
	}

	return pow(fabs(r2-dim),2*alpha) + (0.5*r2 + sum_z)/dim + 0.5;

}

long double Happycat::shiftedRotatedHappycatFunction(int dim, const double* x){

	shift(dim, z, x, shift_vector);

	if (config->getCompetitionID() == MIXTURE || config->getCompetitionID() == CEC14){
		double shiftRate = 5.0/100.0;
		for(int i = 0; i < dim; i++)
			z[i]*=shiftRate;
	}

	rotate(dim, m_r, z, rotation_matrix);
	return 	happycatFunction(dim, m_r);
}
