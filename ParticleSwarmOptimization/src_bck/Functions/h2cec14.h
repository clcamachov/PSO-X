#ifndef HYBRID2CEC14_H
#define HYBRID2CEC14_H

#include "../config.h"
#include "../problem.h"

class Hybrid2CEC14: public Problem {

private:
	Cigar* problem1;
	Rastrigin* problem3;
	Hgbat* problem2;
	double * z;
	double * m_r;

	long double h2Function(int dim, const double* x);

public:
	Hybrid2CEC14(Configuration* config, int variantID);
	~Hybrid2CEC14();

	long double evaluate(int dim, const double* x);

};

#endif
