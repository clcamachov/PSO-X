#ifndef SCHWEFEL_1_2_H
#define SCHWEFEL_1_2_H

#include "../config.h"
#include "../problem.h"

class Schwefel12: public Problem {

private:
	double * z;
	//double * m_r;
	long double schwefel12Function(int dim, const double* x);
	long double shiftedSchwefel12Function(int dim, const double* x);
	long double noiseShiftedSchwefel12Function(int dim, const double* x);
public:

	Schwefel12(Configuration* config, int variantID);
	~Schwefel12();

	long double evaluate(int dim, const double* x);
};

#endif
