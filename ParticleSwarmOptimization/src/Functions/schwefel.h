#ifndef SCHWEFEL_H
#define SCHWEFEL_H

#include "../config.h"
#include "../problem.h"

class Schwefel: public Problem {

private:
	double * z;
	double * m_r;

	long double shiftedRotatedSchwefelFunction(int dim, const double* x);
	long double shiftedSchwefelFunction(int dim, const double* x);
	long double schwefelFunction(int dim, const double * x);

public:

	Schwefel(Configuration* config, int variantID);
	~Schwefel();

	long double evaluate(int dim, const double* x);

};

#endif
