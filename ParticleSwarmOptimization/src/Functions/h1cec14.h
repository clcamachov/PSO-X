#ifndef HYBRID1CEC14_H
#define HYBRID1CEC14_H

#include "../config.h"
#include "../problem.h"

class Hybrid1CEC14: public Problem {

private:
	Schwefel* problem1;
	Rastrigin* problem2;
	Elliptic* problem3;

	double * z;
	double * m_r;

	long double h2Function(int dim, const double* x);

public:
	Hybrid1CEC14(Configuration* config, int variantID);
	~Hybrid1CEC14();

	long double evaluate(int dim, const double* x);

};

#endif
