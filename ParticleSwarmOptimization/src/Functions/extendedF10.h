#ifndef EXTENDED_F10_H
#define EXTENDED_F10_H

#include "../config.h"
#include "../problem.h"

class ExtendedF10: public Problem {

private:
	double * z;

	long double extendedF10Function(int dim, const double* x);
	long double shiftedExtendedF10Function(int dim, const double* x);
	long double f_10(double x, double y);

public:

	ExtendedF10(Configuration* config, int variantID);
	~ExtendedF10();

	long double evaluate(int dim, const double* x);
};

#endif
