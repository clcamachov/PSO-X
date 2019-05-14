#ifndef HYBRID6CEC14_H
#define HYBRID6CEC14_H

#include "../config.h"
#include "../problem.h"

class Hybrid6CEC14: public Problem {

private:
	Katsuura* problem1;
	Happycat* problem2;
	GriewankRosenbrock* problem3;
	Schwefel* problem4;
	Ackley* problem5;
	double * z;
	double * m_r;

	long double h2Function(int dim, const double* x);

public:
	Hybrid6CEC14(Configuration* config, int variantID);
	~Hybrid6CEC14();

	long double evaluate(int dim, const double* x);

};

#endif
