#ifndef HYBRID1_H
#define HYBRID1_H

#include "../config.h"
#include "../problem.h"

class Hybrid1: public Problem {

private:
	ExtendedF10* problem1;
	Sphere* problem2;

	long double h1Function(int dim, const double* x);

public:

	Hybrid1(Configuration* config, int variantID);
	~Hybrid1();

	long double evaluate(int dim, const double* x);

};

#endif
