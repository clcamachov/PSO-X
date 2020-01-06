#ifndef HYBRID4_H
#define HYBRID4_H

#include "../config.h"
#include "../problem.h"

class Hybrid4: public Problem {

private:
	Bohachevsky* problem1;
	Schwefel222* problem2;
	double * z;

	long double h4Function(int dim, const double* x);
	long double shiftedH4Function(int dim, const double* x);
public:

	Hybrid4(Configuration* config, int variantID);
	~Hybrid4();

	long double evaluate(int dim, const double* x);

};

#endif
