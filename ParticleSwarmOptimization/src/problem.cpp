/*
 * problem.cpp
 *
 *  Created on: May 31, 2019
 *      Author: Christian L. Camacho Villalón
 */

#include "gsl/gsl_rng.h"
#include "problem.h"
#include "config.h"
#include "rng.h"
#include <iostream>
#include <cstdlib>
#include <new>

#include <values.h>
#include <float.h>
//Timers
#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

Problem::Problem(Configuration* config, int variantID){
	this->config = config;
	evaluations = 0;
	dimension = config->getProblemDimension();
	maxEvaluations = config->getMaxFES();
	firstEvaluation = true;
	bestSolutionValue = LDBL_MAX; //MAXLONGDOUBLE;
	vID = variantID;
	bestSoFarSolution = new double[config->getProblemDimension()];
	hasBounds = true;
	shift_vector = 0;
	oldResult = 0;
	optimum = 0;
	rotation_matrix = allocateMemory2D(0, 0);
}

Problem::~Problem(){
	delete [] bestSoFarSolution;
}

void Problem::printProblem(){
	cout << "\nProblem parameters:\n"
			<< "  dimensions:         " << getProblemDimension() << "\n"
			<< "  evaluations:        " << getFunctionEvaluations() << "\n"
			<< "  maxEvaluations:     " << maxEvaluations << "\n"
			<< "  optimum:            " << scientific << getProblemOptimum() << "\n"
			<< "  oldResult:          " << scientific << oldResult << "\n"
			<< "  firstEvaluation:    " << firstEvaluation << "\n"
			<< "  bestSolutionValue:  " << getBestSolutionValue() << "\n"
			<< "  hasBounds:          " << hasBoundConstraint() << "\n"
			<< "  vID:                " << vID << "\n"
			<< "  lowerBound:         " << config->getMinInitBound() << "\n"
			<< "  upperBound:         " << config->getMaxInitBound() << "\n";
}

double* Problem::getShift_vector(){
	return (shift_vector);
}

double** Problem::getRotation_matrix(){
	return (rotation_matrix);
}

void Problem::printProblem_results(){
	cout << "Optimum:\t" << scientific << getProblemOptimum() << endl;
}

bool Problem::hasBoundConstraint(){
	return (hasBounds);
}

int Problem::getProblemDimension() {
	return (dimension);
}

double Problem::getProblemOptimum() {
	return (optimum);
}

unsigned int Problem::getFunctionEvaluations(){
	return (evaluations);
}

void Problem::setFunctionEvaluations(unsigned int evals){
	evaluations = evals;
}

long double Problem::getBestSolutionValue(){
	return (bestSolutionValue);
}

double* Problem::getBestSoFarSolution(){
	return (bestSoFarSolution);
}

void Problem::setBestSoFarSolution(const double* solution){
	for(int i = 0; i < config->getProblemDimension(); i++)
		bestSoFarSolution[i] = solution[i];
}

void Problem::setBestSolutionValue(long double value){
	bestSolutionValue = value;
}

long double Problem::getFunctionValue(const double* x) {
	long double result;

	if(evaluations < maxEvaluations){

		result = evaluate(dimension, x);
		//Update number of function evaluations
		evaluations++;

		oldResult = result;

		//Update and print best solution value
		if( bestSolutionValue > result){
			bestSolutionValue = result;
			setBestSoFarSolution(x);
		}

		return (result);
	}
	else
		return (oldResult);

}

void Problem::printProgress(){
	double etime;
	static struct timeval tp;
	gettimeofday( &tp, NULL );
	etime =(double) tp.tv_sec + (double) tp.tv_usec/1000000.0;
	cout.precision(20);
	if(bestSolutionValue > EPSILON)
		cout << evaluations << "\t" <<  scientific <<  bestSolutionValue << "\t" << etime-config->getStartTime()<< endl;
	else
		cout << evaluations << "\t" <<  0.0 << "\t" << etime-config->getStartTime()<< endl;
}

void Problem::divideFunctions(const double *s, double *part1, double *part2, double m, int *psize1, int *psize2) {
	int shared;
	int rest;
	int i;
	int total;
	double *partrest;

	if (m <= 0.5) {
		partrest = part2;
	}
	else {
		partrest = part1;
		m = 1-m;
	}

	shared = (int) floor(dimension*m);
	rest = 2*shared;

	for (i = 0; i < shared; i++) {
		part1[i] = s[i*2];
		part2[i] = s[i*2+1];
	}
	total = dimension-shared;

	for (i = 0; i < total-shared; i++) {
		partrest[i+shared] = s[i+rest];
	}

	*psize1 = shared;
	*psize2 = dimension-shared;

	if (partrest == part1) {
		int temp = *psize1;
		*psize1 = *psize2;
		*psize2 = temp;
	}
}


void Problem::shift(int dim, double * result, const double* x, double* shift_vector){
	for (int i = 0; i < dim; i++){
		result[i] = x[i] - shift_vector[i];
	}
}

void Problem::rotate(int dim, double* result, double* x, double** rotation_matrix){

	for (int i = 0 ; i < dim ; i ++) {
		result[i] = 0.0;
		for (int j = 0 ; j < dim ; j ++) {
			result[i] += (x[j] * rotation_matrix[j][i]);
		}
	}

}

double ** Problem::allocateMemory2D(unsigned int w, unsigned int h){

	double **array2D;

	array2D = new double*[h];
	for (unsigned int i = 0; i < h; ++i)
		array2D[i] = new double[w];

	return (array2D);
}

void Problem::printMatrix(unsigned int dimensions, double** rotation_matrix){
	cout << "\n\nRotationMatrix::\n[ " << endl;
	for (unsigned int i=0; i<dimensions; ++i) {
		for (unsigned int j=0; j<dimensions; ++j) {
			cout << " " << rotation_matrix[i][j];
		}
		cout << endl;
	}
	cout << " ]" << endl;
}


double *** Problem::allocateMemory3D(unsigned int w, unsigned int h, unsigned int d){
	double ***array3D;

	array3D = new double**[h];
	for (unsigned int i = 0; i < h; ++i) {
		array3D[i] = new double*[w];

		for (unsigned int j = 0; j < w; ++j)
			array3D[i][j] = new double[d];
	}
	return (array3D);

}

void Problem::deallocateMemory2D(double ** array2D, unsigned int h){

	for (unsigned int i = 0; i < h; ++i)
		delete [] array2D[i];
	delete [] array2D;
}

void Problem::deallocateMemory3D(double *** array3D, unsigned int h, unsigned int w){

	for (unsigned int i = 0; i < h; ++i) {
		for (unsigned int j = 0; j < w; ++j)
			delete [] array3D[i][j];

		delete [] array3D[i];
	}
	delete [] array3D;

}

void Problem::oszfunc (double *x, double *xosz, int dim)
{
	int i;
	int sx;
	double c1;
	double c2;
	double xx=0;
	for (i=0; i<dim; i++)
	{
		if (i==0||i==dim-1)
		{
			if (x[i]!=0)
				xx=log(fabs(x[i]));
			if (x[i]>0)
			{
				c1=10;
				c2=7.9;
			}
			else
			{
				c1=5.5;
				c2=3.1;
			}
			if (x[i]>0)
				sx=1;
			else if (x[i]==0)
				sx=0;
			else
				sx=-1;
			xosz[i]=sx*exp(xx+0.049*(sin(c1*xx)+sin(c2*xx)));
		}
		else
			xosz[i]=x[i];
	}
}

double Problem::xRound(double x, double o){
	return ((fabs(x - o) < 0.5) ? x : (roundVal(2.0 * x) / 2.0));
}

double Problem::xRound(double x){
	return ((fabs(x) < 0.5) ? x : (roundVal(2.0 * x) / 2.0));
}

double Problem::roundVal(double x){
	return (signnum(x) * floor(abs(x) + 0.5));
}

double Problem::signnum(double x){
	if (x > 0) return (1);
	if (x < 0) return (-1);
	return (0);
}

void Problem::multiplyMatrix(double* result, double* vector, double** matrix){

	for (int i = 0 ; i < dimension ; i ++) {
		result[i] = 0.0;
		for (int j = 0 ; j < dimension ; j ++) {
			result[i] += (vector[j] * matrix[j][i]);
		}
	}
}

void Problem::multiplyMatrix(double* result, double** matrix, double* vector){
	for (int i = 0 ; i < dimension ; i ++) {
		result[i] = 0.0;
		for (int j = 0 ; j < dimension ; j ++) {
			result[i] += (vector[i] * matrix[i][j]);
		}
	}

}

void Problem::multiplyMatrix(double* result, double** matrix, const double* vector){
	for (int i = 0 ; i < dimension ; i ++) {
		result[i] = 0.0;
		for (int j = 0 ; j < dimension ; j ++) {
			result[i] += (vector[i] * matrix[i][j]);
		}
	}

}
