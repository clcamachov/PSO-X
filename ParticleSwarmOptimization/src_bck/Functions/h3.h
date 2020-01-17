#ifndef HYBRID3_H
#define HYBRID3_H

#include "../config.h"
#include "../problem.h"

class Hybrid3: public Problem {

private:
	ExtendedF10* problem1;
	Rastrigin* problem2;

	long double h3Function(int dim, const double* x);

public:

	Hybrid3(Configuration* config, int variantID);
	~Hybrid3();

	long double evaluate(int dim, const double* x);

};

#endif
