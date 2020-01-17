#ifndef SCHWEFEL_2_13_H
#define SCHWEFEL_2_13_H

#include "../config.h"
#include "../problem.h"

class Schwefel213: public Problem {

private:
	double * m_A;
	double * m_B;
	double ** m_a;
	double ** m_b;
	//double ** data_matrix;

	long double schwefel213Function(int dim, const double* x);

public:

	Schwefel213(Configuration* config, int variantID);
	~Schwefel213();

	long double evaluate(int dim, const double* x);

};

#endif
