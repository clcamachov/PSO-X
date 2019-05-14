#ifndef KATSUURA_H
#define KATSUURA_H

#include "../config.h"
#include "../problem.h"

class Katsuura: public Problem {

private:
	double * z;
	double * m_r;

	long double shiftedRotatedKatsuuraFunction(int dim, const double* x);
	long double katsuuraFunction(int dim, const double * x);

public:

	Katsuura(Configuration* config, int variantID);
	~Katsuura();

	long double evaluate(int dim, const double* x);

};

#endif
