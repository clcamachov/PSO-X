#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "sphere.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;


Sphere::Sphere(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == SHIFTED){
		z = new double[dimension];
		shift_vector = new double[dimension];
		string file_data;
		if(config->getCompetitionID() == CEC05){
			file_data = "supportData/sphere_CEC05_func_data.txt";

		}
		else{
			file_data = "supportData/sphere_CEC08_func_data.txt";
		}

		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
	}
}

Sphere::~Sphere(){
	if(vID == SHIFTED){
		delete [] shift_vector;
		delete [] z;
	}
}

long double Sphere::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = sphereFunction(dim, x);
		break;
	case SHIFTED:
		result = shiftedSphereFunction(dim, x);
		break;
	case NOISE_BASIC:
		result = noiseSphereFunction(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}


long double Sphere::sphereFunction(int dim, const double * x){
	long double F = 0;
	for(int i=0;i<dim;i++){
		F += x[i]*x[i];
	}
	return F;
}


long double Sphere::shiftedSphereFunction(int dim, const double* x){

	shift(dim, z, x, shift_vector);

	long double F = sphereFunction(dim, z);
	return F;
}


long double Sphere::noiseSphereFunction(int dim, const double* x){
	long double F = sphereFunction(dim, x);
	F *= (1.0 + 0.1 * fabs(RNG::randGauss(config->getRNGSeed())));
	return F;
}
