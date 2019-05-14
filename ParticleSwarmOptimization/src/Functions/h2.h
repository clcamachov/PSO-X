#ifndef HYBRID2_H
#define HYBRID2_H

#include "../config.h"
#include "../problem.h"

class Hybrid2: public Problem {

private:
	ExtendedF10* problem1;
	Rosenbrock* problem2;

	long double h2Function(int dim, const double* x);
public:

	Hybrid2(Configuration* config, int variantID);
	~Hybrid2();

	long double evaluate(int dim, const double* x);

};

#endif
