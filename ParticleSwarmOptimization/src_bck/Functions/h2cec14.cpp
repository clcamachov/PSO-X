#include "../rng.h"
#include "../problem.h"
#include "../utils.h"
#include "rastrigin.h"
#include "cigar.h"
#include "hgbat.h"
#include "h2cec14.h"



#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

Hybrid2CEC14::Hybrid2CEC14(Configuration* config, int variantID):Problem(config, variantID){
	optimum = 0.0;
	problem1 = new Cigar(config, BASIC);
	problem2 = new Hgbat(config, BASIC);
	problem3 = new Rastrigin(config, BASIC);

	if(vID == SHIFTED_ROTATED){
		string file_data = "supportData/input_data/shift_data_18.txt";
		stringstream dim_name;
		dim_name << dimension;
		string file_m = "supportData/input_data/M_18_D" + dim_name.str() + ".txt";

		rotation_matrix  = allocateMemory2D(dimension, dimension);
		shift_vector = new double[dimension];
		z = new double[dimension];
		m_r = new double[dimension];

		Utils::loadRowVectorFromFile(file_data, dimension, shift_vector);
		Utils::loadMatrixFromFile(file_m, dimension, dimension, rotation_matrix);

	}
}

Hybrid2CEC14::~Hybrid2CEC14(){

	delete problem1;

	delete problem2;

	delete problem3;

	if(vID == SHIFTED || vID == SHIFTED_ROTATED || vID == ROTATED_WITHOUT_BOUNDS){
		delete [] shift_vector;
		delete [] z;
	}

	if(vID == SHIFTED_ROTATED || vID == ROTATED_WITHOUT_BOUNDS){
		delete [] m_r;
		deallocateMemory2D(rotation_matrix, dimension);
	}
}

long double Hybrid2CEC14::evaluate(int dim, const double* x) {
	long double result = 0.0;

	switch(vID){
	case SHIFTED_ROTATED:
		result = h2Function(dim, x);
		break;
	default: cout<<"There is no function such that !!"<< endl;
	exit(0);
	break;
	} 

	return result;

}



long double Hybrid2CEC14::h2Function(int dim, const double* x){

	//int size1, size2;
	long double f = 0.0;

	int i,tmp,cf_num=3;
	//double fit[3];
	int G[3],G_nx[3];
	double Gp[3]={0.3,0.3,0.4};

	tmp=0;
	//find the starting and end points of each parts
	for (i=0; i<cf_num-1; i++)
	{
		G_nx[i] = ceil(Gp[i]*dim);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1]=dim-tmp;
	G[0]=0;
	for (i=1; i<cf_num; i++)
	{
		G[i] = G[i-1]+G_nx[i-1];
	}

	shift(dim, z, x, shift_vector);
	rotate(dim, m_r, z, rotation_matrix);

	double part1[G_nx[0]];
	for(int j= 0; j < G_nx[0]; j++)// divide into parts
			part1[j]=m_r[G[0]+j];
	f += problem1->evaluate(G_nx[0],part1);

	double part2[G_nx[1]];
	for(int j= 0; j < G_nx[1]; j++)// divide into parts
		part2[j]=m_r[G[1]+j];
	f += problem2->evaluate(G_nx[1],part2);

	double part3[G_nx[2]];
	for(int j= 0; j < G_nx[2]; j++)// divide into parts
		part3[j]=m_r[G[2]+j];
	f += problem2->evaluate(G_nx[2],part3);

	return f;
}
