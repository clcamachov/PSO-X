#ifndef SCHWEFEL_2_22_H
#define SCHWEFEL_2_22_H

#include "../config.h"
#include "../problem.h"

class Schwefel222: public Problem {

private:

	double * z;
	long double schwefel222Function(int dim, const double* x);
	long double shiftedSchwefel222Function(int dim, const double* x);
public:

	Schwefel222(Configuration* config, int variantID);
	~Schwefel222();

	long double evaluate(int dim, const double* x);

};

#endif
