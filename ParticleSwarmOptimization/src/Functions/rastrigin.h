#ifndef RASTRIGIN_H
#define RASTRIGIN_H

#include "../config.h"
#include "../problem.h"

class Rastrigin: public Problem {

private:
	double * z;
	double * m_r;

	long double rastriginFunction(int dim, const double* x);
	long double noncontinuousRastriginFunction(int dim, const double* x);
	long double shiftedRastriginFunction(int dim, const double* x);
	long double shiftedRotatedRastriginFunction(int dim, const double* x);

public:

	Rastrigin(Configuration* config, int variantID);
	~Rastrigin();

	long double evaluate(int dim, const double* x);

};

#endif
