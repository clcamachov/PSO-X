#ifndef SCHAFFER_F6_H
#define SCHAFFER_F6_H

#include "../config.h"
#include "../problem.h"

class SchafferF6: public Problem {

private:
	double * z;
	double * m_r;

	double schafferF6Function(double x, double y);
	long double expandedSchafferF6Function(int dim, const double* x);
	long double shiftedRotatedExpandedSchafferF6Function(int dim, const double* x);
	long double noncontinuousExpandedSchafferF6Function(int dim, const double* x);

public:

	SchafferF6(Configuration* config, int variantID);
	~SchafferF6();

	long double evaluate(int dim, const double* x);
};

#endif
