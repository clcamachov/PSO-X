#ifndef GRIEWANK_ROSENBROCK_H
#define GRIEWANK_ROSENBROCK_H

#include "../config.h"
#include "../problem.h"

class GriewankRosenbrock: public Problem {

private:

	double * z;
	double * m_r;
	double griewank1DFunction(const double x);
	double rosenbrock2DFunction(double x, double y);

	long double shiftedRotatedExpandedGriewankRosenbrockFunction(int dim, const double* x);
	long double shiftedExpandedGriewankRosenbrockFunction(int dim, const double* x);
	long double griewankRosenbrockFunction(int dim, const double* x);
	long double griewankRosenbrockFunctionCEC14(int dim, double* x);

public:

	GriewankRosenbrock(Configuration* config, int variantID);
	~GriewankRosenbrock();

	long double evaluate(int dim, const double* x);
};

#endif
