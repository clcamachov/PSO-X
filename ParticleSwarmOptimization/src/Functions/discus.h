#ifndef DISCUS_H
#define DISCUS_H

#include "../config.h"
#include "../problem.h"

class Discus: public Problem {

private:
	double * z;
	double * m_r;

	long double shiftedRotatedDiscusFunction(int dim, const double* x);
	long double discusFunction(int dim, const double * x);

public:

	Discus(Configuration* config, int variantID);
	~Discus();

	long double evaluate(int dim, const double* x);

};

#endif
