#ifndef SCHWEFEL_2_21_H
#define SCHWEFEL_2_21_H

#include "../config.h"
#include "../problem.h"

class Schwefel221: public Problem {

private:
	double * z;

	long double schwefel221Function(int dim, const double* x);
	long double shiftedSchwefel221Function(int dim, const double* x);

public:

	Schwefel221(Configuration* config, int variantID);
	~Schwefel221();

	long double evaluate(int dim, const double* x);
};

#endif
