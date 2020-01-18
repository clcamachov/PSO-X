#include "../problem.h"
#include "../utils.h"
#include "katsuura.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;


Katsuura::Katsuura(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == MIXTURE || vID == SHIFTED_ROTATED){
		string file_data = "supportData/input_data/shift_data_12.txt";
		stringstream dim_name;
		dim_name << dimension;
		string file_m; // = "supportData/input_data/M_12_D" + dim_name.str() + ".txt";

		if (dimension > 2 && dimension < 10)
			file_m = "supportData/input_data/M_12_D10.txt";
		else if (dimension > 10 && dimension < 20)
			file_m = "supportData/input_data/M_12_D20.txt";
		else if (dimension > 20 && dimension < 30)
			file_m = "supportData/input_data/M_12_D30.txt";
		else if (dimension > 30 && dimension < 50)
			file_m = "supportData/input_data/M_12_D50.txt";
		else if (dimension > 50 && dimension < 100)
			file_m = "supportData/input_data/M_12_D100.txt";
		else
			file_m = "supportData/input_data/M_12_D" + dim_name.str() + ".txt";

		rotation_matrix  = allocateMemory2D(dimension, dimension);
		shift_vector = new double[dimension];
		z = new double[dimension];
		m_r = new double[dimension];

		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);

	}


}

Katsuura::~Katsuura(){
	if(vID == SHIFTED_ROTATED || vID == MIXTURE){
		deallocateMemory2D(rotation_matrix, dimension);
		delete [] shift_vector;
		delete [] z;
		delete [] m_r;
	}
}

long double Katsuura::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = katsuuraFunction(dim, x);
		break;
	case SHIFTED_ROTATED:
	case MIXTURE:
		result = shiftedRotatedKatsuuraFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;
}


long double Katsuura::katsuuraFunction(int dim, const double* x){
	long double sum = 1.0;
	int i,j;
	double temp,tmp1,tmp2,tmp3;

	tmp3=pow(1.0*dim,1.2);

	for (i=0; i<dim; i++)
	{
		temp=0.0;
		for (j=1; j<=32; j++)
		{
			tmp1=pow(2.0,j);
			tmp2=tmp1*x[i];
			temp += fabs(tmp2-floor(tmp2+0.5))/tmp1;
		}
		sum *= pow(1.0+(i+1)*temp,10.0/tmp3);
	}
	tmp1=10.0/dim/dim;
	sum=sum*tmp1-tmp1;
	return sum;
}

long double Katsuura::shiftedRotatedKatsuuraFunction(int dim, const double* x){

	shift(dim, z, x, shift_vector);

	if (config->getCompetitionID() == MIXTURE || config->getCompetitionID() == CEC14){
		double shiftRate = 5.0/100.0;
		for(int i = 0; i < dim; i++)
			z[i]*=shiftRate;
	}

	rotate(dim, m_r, z, rotation_matrix);
	return katsuuraFunction(dim, m_r);
}
