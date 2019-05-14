#ifndef ROSENBROCK_H
#define ROSENBROCK_H

#include "../config.h"
#include "../problem.h"

class Rosenbrock: public Problem {

private:
	double * z;
	double * m_r;
	long double shiftedRosenbrockFunction(int dim, const double* x);
	long double shiftedRotatedRosenbrockFunction(int dim, const double* x);
	long double rosenbrockFunction(int dim, const double* x);

public:

	Rosenbrock(Configuration* config, int variantID);
	~Rosenbrock();

	long double evaluate(int dim, const double* x);

};

#endif
