#ifndef HYBRID7_H
#define HYBRID7_H

#include "../config.h"
#include "../problem.h"

class Hybrid7: public Problem {

private:
	ExtendedF10* problem1;
	Sphere* problem2;

	long double h7Function(int dim, const double* x);

public:

	Hybrid7(Configuration* config, int variantID);
	~Hybrid7();

	long double evaluate(int dim, const double* x);

};

#endif
