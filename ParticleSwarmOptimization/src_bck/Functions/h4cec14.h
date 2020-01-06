#ifndef HYBRID4CEC14_H
#define HYBRID4CEC14_H

#include "../config.h"
#include "../problem.h"

class Hybrid4CEC14: public Problem {

private:
	Hgbat* problem1;
	Discus* problem2;
	GriewankRosenbrock* problem3;
	Rastrigin* problem4;
	double * z;
	double * m_r;

	long double h2Function(int dim, const double* x);

public:
	Hybrid4CEC14(Configuration* config, int variantID);
	~Hybrid4CEC14();

	long double evaluate(int dim, const double* x);

};

#endif
