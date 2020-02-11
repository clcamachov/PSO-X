#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "extendedF10.h"


#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

ExtendedF10::ExtendedF10(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == SHIFTED){
		shift_vector = new double[dimension];
		string file_data = "../supportData/extendedF10_CEC08_func_data.txt";
		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

		z = new double[dimension];
	}

}

ExtendedF10::~ExtendedF10(){
	if(vID == SHIFTED){
		delete [] z;
		delete [] shift_vector;
	}
}

long double ExtendedF10::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = extendedF10Function(dim, x);
		break;
	case SHIFTED:
		result = shiftedExtendedF10Function(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}


long double ExtendedF10::f_10(double x, double y)
{
	double p, z, t;

	p=(x*x+y*y);

	z=pow(p, 0.25);
	t=sin(50.0*pow(p, 0.1));
	t=t*t+1.0;

	return z*t;
}

long double ExtendedF10::shiftedExtendedF10Function(int dim, const double* x){
	long double F=0.0;

	shift(dim, z, x, shift_vector);

	F = extendedF10Function(dim, z);

	return F;
}

long double ExtendedF10::extendedF10Function(int dim, const double* x){
	long double suma=0.0;

	for(int i=0; i<dim-1; i++)
		suma+=f_10(x[i], x[i+1]);

	suma+=f_10(x[dim-1], x[0]);

	return suma;
}
