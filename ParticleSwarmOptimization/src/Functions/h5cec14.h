#ifndef HYBRID5CEC14_H
#define HYBRID5CEC14_H

#include "../config.h"
#include "../problem.h"

class Hybrid5CEC14: public Problem {

private:
	SchafferF6* problem1;
	Hgbat* problem2;
	Rosenbrock* problem3;
	Schwefel* problem4;
	Elliptic* problem5;
	double * z;
	double * m_r;

	long double h2Function(int dim, const double* x);

public:
	Hybrid5CEC14(Configuration* config, int variantID);
	~Hybrid5CEC14();

	long double evaluate(int dim, const double* x);

};

#endif
