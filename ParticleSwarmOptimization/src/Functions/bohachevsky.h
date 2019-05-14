#ifndef BOHACHEVSKY_H
#define BOHACHEVSKY_H

#include "../config.h"
#include "../problem.h"

class Bohachevsky: public Problem {

private:
	double * z;

	long double bohachevskyFunction(int dim, const double* x);
	long double shiftedBohachevskyFunction(int dim, const double* x);

public:

	Bohachevsky(Configuration* config, int variantID);
	~Bohachevsky();

	long double evaluate(int dim, const double* x);
};

#endif
