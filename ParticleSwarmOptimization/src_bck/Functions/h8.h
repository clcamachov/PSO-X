#ifndef HYBRID8_H
#define HYBRID8_H

#include "../config.h"
#include "../problem.h"

class Hybrid8: public Problem {

private:
	ExtendedF10* problem1;
	Rosenbrock* problem2;

	long double h8Function(int dim, const double* x);
public:

	Hybrid8(Configuration* config, int variantID);
	~Hybrid8();

	long double evaluate(int dim, const double* x);

};

#endif
