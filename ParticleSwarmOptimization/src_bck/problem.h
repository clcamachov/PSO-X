/*
 * problem.h
 *
 *  Created on: May 31, 2018
 *      Author: leonardo
 */

#ifndef PROBLEM_H_
#define PROBLEM_H_

#include "config.h"
#include <cmath>
#include <algorithm>

#define abss(a)     (a<0 ? (-a) : a)

#define BASIC											0
#define SHIFTED 										1
#define ROTATED 										2
#define	EXPANDED										3
#define SHIFTED_EXPANDED								4
#define SHIFTED_ROTATED 								5
#define NOISE_BASIC										6
#define NOISE_SHIFTED									7
#define NOISE_ROTATED									8
#define BASIC_GLOBAL_OPTIMUM_ON_BOUNDS					9
#define SHIFTED_ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS		10
#define	SHIFTED_ROTATED_HIGH_CONDITIONED				11
#define ROTATED_WITH_NORROW_BASIN_GLOBAL_OPTIMUM		12
#define ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS				13
#define ROTATED_WITH_HIGH_CONDITION_NUMBER_MATRIX		14
#define NONCONTINUOUS_ROTATED							15
#define NONCONTINUOUS_BASIC								16
#define SHIFTED_ROTATED_EXPANDED						17
#define NONCONTINUOUS_EXPANDED							18
#define ROTATED_WITHOUT_BOUNDS							19
#define SHIFTED_ROTATED_HIGH_CONDITIONED_CEC2013		20

class Problem {

protected:
	Configuration* config;
	int dimension;
	unsigned int evaluations;
	unsigned int maxEvaluations;
	double optimum;
	long double oldResult;

	bool firstEvaluation;
	long double bestSolutionValue;
	double* bestSoFarSolution;

	double * shift_vector;
	double ** rotation_matrix;

	bool hasBounds;

	int vID;	//Variant ID: from 1 to 20

	void shift(int dim, double * result, const double* x, double* shift_vector);
	void rotate(int dim, double * result, double* x, double** rotation_matrix);

	void oszfunc (double *x, double *xosz, int dim);

	double ** allocateMemory2D(unsigned int w, unsigned int h);
	double *** allocateMemory3D(unsigned int w, unsigned int h, unsigned int d);

	void deallocateMemory2D(double ** array2D, unsigned int h);
	void deallocateMemory3D(double *** array3D, unsigned int h, unsigned int w);

	double xRound(double x, double o);
	double xRound(double x);

	void multiplyMatrix(double* result, double* vector, double** matrix);
	void multiplyMatrix(double* result, double** matrix, double* vector);
	void multiplyMatrix(double* result, double** matrix, const double* vector);

private:
	double roundVal(double x);
	double signnum(double x);

public:

	Problem(Configuration* config, int variantID);
	virtual ~Problem();

	int getProblemDimension();
	double getProblemOptimum();

	unsigned int getFunctionEvaluations();

	void setFunctionEvaluations(unsigned int evals);
	//virtual void printShift_vector2D();

	//Evaluation call
	virtual long double evaluate(int dim, const double* x) = 0;

	long double getFunctionValue(const double* x);

	long double getBestSolutionValue();
	double* getBestSoFarSolution();
	void setBestSoFarSolution(const double* solution);

	void setBestSolutionValue(long double value);

	bool hasBoundConstraint();

	void printProgress();

	void divideFunctions(const double *s,  double *part1, double *part2, double m, int *psize1, int *psize2);

	double* getShift_vector();
	double** getRotation_matrix();
	void printProblem_results();


	/*
	 * Function from Optim.h
	 */
	double getRandomX();
	double getRandomX(double lowerBound, double upperBound);
	double getRandom01();
	void printProblem();
};


#endif /* PROBLEM_H_ */
