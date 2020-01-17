#ifndef HGBAT_H
#define HGBAT_H

#include "../config.h"
#include "../problem.h"

class Hgbat: public Problem {

private:
	double * z;
	double * m_r;

	long double shiftedRotatedHgbatFunction(int dim, const double* x);
	long double hgbatFunction(int dim, const double * x);

public:

	Hgbat(Configuration* config, int variantID);
	~Hgbat();

	long double evaluate(int dim, const double* x);

};

#endif
