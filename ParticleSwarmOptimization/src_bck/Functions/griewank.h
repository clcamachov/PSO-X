#ifndef SHIFTED_GRIEWANK_H
#define SHIFTED_GRIEWANK_H

#include "../config.h"
#include "../problem.h"

class Griewank: public Problem {

private:
	double * z;
	double * m_r;
	long double griewankFunction(int dim, const double* x);
	long double shiftedGriewankFunction(int dim, const double* x);
	long double shiftedRotatedGriewankFunction(int dim, const double* x);

public:

	Griewank(Configuration* config, int variantID);
	~Griewank();

	long double evaluate(int dim, const double* x);
};

#endif
