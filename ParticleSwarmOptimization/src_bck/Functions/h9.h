#ifndef HYBRID9_H
#define HYBRID9_H

#include "../config.h"
#include "../problem.h"

class Hybrid9: public Problem {

private:
	ExtendedF10* problem1;
	Rastrigin* problem2;

	long double h9Function(int dim, const double* x);

public:

	Hybrid9(Configuration* config, int variantID);
	~Hybrid9();

	long double evaluate(int dim, const double* x);

};

#endif
