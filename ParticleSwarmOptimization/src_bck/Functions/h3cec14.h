#ifndef HYBRID3CEC14_H
#define HYBRID3CEC14_H

#include "../config.h"
#include "../problem.h"

class Hybrid3CEC14: public Problem {

private:
	Griewank* problem1;
	Weierstrass* problem2;
	Rosenbrock* problem3;
	SchafferF6* problem4;
	double * z;
	double * m_r;

	long double h2Function(int dim, const double* x);

public:
	Hybrid3CEC14(Configuration* config, int variantID);
	~Hybrid3CEC14();

	long double evaluate(int dim, const double* x);

};

#endif
