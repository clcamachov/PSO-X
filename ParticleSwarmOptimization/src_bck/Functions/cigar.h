#ifndef CIGAR_H
#define CIGAR_H

#include "../config.h"
#include "../problem.h"

class Cigar: public Problem {

private:
	double * z;
	double * m_r;

	long double shiftedRotatedCigarFunction(int dim, const double* x);
	long double cigarFunction(int dim, const double * x);

public:

	Cigar(Configuration* config, int variantID);
	~Cigar();

	long double evaluate(int dim, const double* x);

};

#endif
