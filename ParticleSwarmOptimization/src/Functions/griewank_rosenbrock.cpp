#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "griewank_rosenbrock.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

GriewankRosenbrock::GriewankRosenbrock(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;

	if(vID == SHIFTED_EXPANDED){
		shift_vector = new double[dimension];
		z = new double[dimension];
		string file_data = "../supportData/EF8F2_func_data.txt";
		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);

		for(int i = 0; i < dimension; i++){
			shift_vector[i] = shift_vector[i] - 1.0;
		}
	}
	if(vID == SHIFTED_ROTATED){
		string file_data = "../supportData/input_data/shift_data_15.txt";
		stringstream dim_name;
		dim_name << dimension;
		string file_m; // = "../supportData/input_data/M_15_D" + dim_name.str() + ".txt";

		if (dimension > 2 && dimension < 10)
			file_m = "../supportData/input_data/M_15_D10.txt";
		else if (dimension > 10 && dimension < 20)
			file_m = "../supportData/input_data/M_15_D20.txt";
		else if (dimension > 20 && dimension < 30)
			file_m = "../supportData/input_data/M_15_D30.txt";
		else if (dimension > 30 && dimension < 50)
			file_m = "../supportData/input_data/M_15_D50.txt";
		else if (dimension > 50 && dimension < 100)
			file_m = "../supportData/input_data/M_15_D100.txt";
		else
			file_m =  "../supportData/input_data/M_15_D" + dim_name.str() + ".txt";

		rotation_matrix  = allocateMemory2D(dimension, dimension);
		shift_vector = new double[dimension];
		z = new double[dimension];
		m_r = new double[dimension];

		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);
	}
}

GriewankRosenbrock::~GriewankRosenbrock(){
	if(vID == SHIFTED_EXPANDED){
		delete [] shift_vector;
		delete [] z;
	}
}

long double GriewankRosenbrock::evaluate(int dim, const double* x) {
	long double result;

	switch(vID){
	case BASIC:
		result = griewankRosenbrockFunction(dim, x);
		break;
	case SHIFTED_EXPANDED:
		result = shiftedExpandedGriewankRosenbrockFunction(dim, x);
		break;
	case SHIFTED_ROTATED:
		result =shiftedRotatedExpandedGriewankRosenbrockFunction(dim,x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}

long double GriewankRosenbrock::griewankRosenbrockFunction(int dim, const double * x){

	double sum = 0.0;

	for (int i = 1 ; i < dim ; i ++) {
		sum += griewank1DFunction(rosenbrock2DFunction(x[i-1], x[i]));
	}

	sum += griewank1DFunction(rosenbrock2DFunction(x[dim-1], x[0]));

	return sum;

}

long double GriewankRosenbrock::griewankRosenbrockFunctionCEC14(int dim, double * x){

	double sum = 0.0;

	int i;
	double temp,tmp1,tmp2;

	x[0] += 1.0;//shift to orgin
	for (i=0; i<dim-1; i++)
	{
		x[i+1] += 1.0;//shift to orgin
		tmp1 = x[i]*x[i]-x[i+1];
		tmp2 = x[i]-1.0;
		temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
		sum += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}

	tmp1 = x[dim-1]*x[dim-1]-x[0];
	tmp2 = x[dim-1]-1.0;
	temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
	sum += (temp*temp)/4000.0 - cos(temp) + 1.0 ;
	return sum;
}

long double GriewankRosenbrock::shiftedRotatedExpandedGriewankRosenbrockFunction(int dim, const double* x){

	long double F = 0;

	shift(dim, z, x, shift_vector);

	if (config->getCompetitionID() == MIXTURE){
		double shiftRate = 5.0/100.0;
		for(int i = 0; i < dim; i++)
			z[i]*=shiftRate;
	}
	rotate(dim, m_r, z, rotation_matrix);
	F = griewankRosenbrockFunctionCEC14(dim, m_r);

	return F;
}


long double GriewankRosenbrock::shiftedExpandedGriewankRosenbrockFunction(int dim, const double* x){

	long double F = 0;

	shift(dim, z, x, shift_vector);

	F = griewankRosenbrockFunction(dim, z);

	return F;
}

double GriewankRosenbrock::griewank1DFunction(const double x){
	return (((x * x) / 4000.0) - cos(x) + 1.0);
}

double GriewankRosenbrock::rosenbrock2DFunction(double x, double y){
	double temp1 = (x * x) - y;
	double temp2 = x - 1.0;
	return ((100.0 * temp1 * temp1) + (temp2 * temp2));
}
