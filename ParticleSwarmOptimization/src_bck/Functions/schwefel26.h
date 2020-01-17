#ifndef SCHWEFEL_2_6_H
#define SCHWEFEL_2_6_H

#include "../config.h"
#include "../problem.h"

class Schwefel26: public Problem {

private:
	double * v_B;
	double ** m_A;
	double * z;

	long double schwefel26GlobalOptimumOnBoundsFunction(int dim, const double* x);

public:

	Schwefel26(Configuration* config, int variantID);
	~Schwefel26();

	long double evaluate(int dim, const double* x);

};

#endif
