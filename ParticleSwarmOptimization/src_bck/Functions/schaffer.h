#ifndef SCHAFFER_H
#define SCHAFFER_H

#include "../config.h"
#include "../problem.h"

class Schaffer: public Problem {

private:
	double * z;

	long double schafferFunction(int dim, const double* x);
	long double shiftedSchafferFunction(int dim, const double* x);
public:

	Schaffer(Configuration* config, int variantID);
	~Schaffer();

	long double evaluate(int dim, const double* x);

};

#endif
