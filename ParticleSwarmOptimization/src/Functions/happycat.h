#ifndef HAPPYCAT_H
#define HAPPYCAT_H

#include "../config.h"
#include "../problem.h"

class Happycat: public Problem {

private:
	double * z;
	double * m_r;

	long double shiftedRotatedHappycatFunction(int dim, const double* x);
	long double happycatFunction(int dim, const double * x);

public:

	Happycat(Configuration* config, int variantID);
	~Happycat();

	long double evaluate(int dim, const double* x);

};

#endif
