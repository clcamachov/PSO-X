#ifndef WEIERSTRASS_H
#define WEIERSTRASS_H

#include "../config.h"
#include "../problem.h"

class Weierstrass: public Problem {

private:
	double * z;
	double * m_r;

	long double weierstrassFunction(int dim, const double* x);
	long double shiftedRotatedWeierstrassFunction(int dim, const double* x);

public:

	Weierstrass(Configuration* config, int variantID);
	~Weierstrass();

	long double evaluate(int dim, const double* x);

};

#endif
