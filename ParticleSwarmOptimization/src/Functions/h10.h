#ifndef HYBRID10_H
#define HYBRID10_H

#include "../config.h"
#include "../problem.h"

class Hybrid10: public Problem {

private:
	Bohachevsky* problem1;
	Schwefel222* problem2;
	double * z;

	long double h10Function(int dim, const double* x);
	long double shiftedH10Function(int dim, const double * x);

public:

	Hybrid10(Configuration* config, int variantID);
	~Hybrid10();

	long double evaluate(int dim, const double* x);

};

#endif
