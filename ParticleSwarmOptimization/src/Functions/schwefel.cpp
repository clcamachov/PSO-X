#include "../problem.h"
#include "../utils.h"
#include "schwefel.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;


Schwefel::Schwefel(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == MIXTURE || vID == SHIFTED_ROTATED){
		string file_data = "../supportData/input_data/shift_data_11.txt";
		stringstream dim_name;
		dim_name << dimension;
		string file_m; // = "../supportData/input_data/M_11_D" + dim_name.str() + ".txt";


		if (dimension > 2 && dimension < 10)
			file_m = "../supportData/input_data/M_11_D10.txt";
		else if (dimension > 10 && dimension < 20)
			file_m = "../supportData/input_data/M_11_D20.txt";
		else if (dimension > 20 && dimension < 30)
			file_m = "../supportData/input_data/M_11_D30.txt";
		else if (dimension > 30 && dimension < 50)
			file_m = "../supportData/input_data/M_11_D50.txt";
		else if (dimension > 50 && dimension < 100)
			file_m = "../supportData/input_data/M_11_D100.txt";
		else
			file_m = "../supportData/input_data/M_11_D" + dim_name.str() + ".txt";

		rotation_matrix  = allocateMemory2D(dimension, dimension);
		shift_vector = new double[dimension];
		z = new double[dimension];
		m_r = new double[dimension];

		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);

	}

	if(vID == MIXTURE || vID == SHIFTED){
		string file_data = "../supportData/input_data/shift_data_10.txt";
		shift_vector = new double[dimension];
		z = new double[dimension];
		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

	}

}

Schwefel::~Schwefel(){
	if(vID == SHIFTED_ROTATED || vID == MIXTURE){
		deallocateMemory2D(rotation_matrix, dimension);
		delete [] shift_vector;
		delete [] z;
		delete [] m_r;
	}
	if(vID == SHIFTED || vID == MIXTURE){
		delete [] shift_vector;
		delete [] z;
	}

}

long double Schwefel::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = schwefelFunction(dim, x);
		break;
	case SHIFTED:
		result = shiftedSchwefelFunction(dim, x);
		break;
	case SHIFTED_ROTATED:
		result = shiftedRotatedSchwefelFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;
}


long double Schwefel::schwefelFunction(int dim, const double* x){
	long double sum= 0.0;
	int i;
	double tmp;

	for (i=0; i<dim; i++)
	{
		sum += 4.209687462275036e+002;
		if (x[i]>500)
		{
			sum-=(500.0-fmod(x[i],500))*sin(pow(500.0-fmod(x[i],500),0.5));
			tmp=(x[i]-500.0)/100;
			sum+= tmp*tmp/dim;
		}
		else if (x[i]<-500)
		{
			sum-=(-500.0+fmod(fabs(x[i]),500))*sin(pow(500.0-fmod(fabs(x[i]),500),0.5));
			tmp=(x[i]+500.0)/100;
			sum+= tmp*tmp/dim;
		}
		else
			sum-=x[i]*sin(pow(fabs(x[i]),0.5));
	}
	sum +=4.189828872724338e+002*dim;
	return sum;
}

long double Schwefel::shiftedRotatedSchwefelFunction(int dim, const double* x){

	shift(dim, z, x, shift_vector);
	if (config->getCompetitionID() == MIXTURE || config->getCompetitionID() == CEC14){
		double shiftRate = 1000.0/100.0;
		for(int i = 0; i < dim; i++)
			z[i]*=shiftRate;
	}

	rotate(dim, m_r, z, rotation_matrix);

	return schwefelFunction(dim,m_r);

}


long double Schwefel::shiftedSchwefelFunction(int dim, const double* x){

	shift(dim, z, x, shift_vector);
	if (config->getCompetitionID() == MIXTURE || config->getCompetitionID() == CEC14){
		double shiftRate = 1000.0/100.0;
		for(int i = 0; i < dim; i++)
			z[i]*=shiftRate;
	}

	return schwefelFunction(dim,z);

}
