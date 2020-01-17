#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include "../config.h"
#include "../problem.h"

class Elliptic: public Problem {

private:
	double * z;
	double * m_r;
	long double ellipticFunction(int dim, const double* x);
	long double shiftedRotatedHighConditionedEllipticFunction(int dim, const double* x);
	long double shiftedRotatedHighConditionedEllipticFunction2(int dim, const double* x);
public:

	Elliptic(Configuration* config, int variantID);
	~Elliptic();

	long double evaluate(int dim, const double* x);

};

#endif
