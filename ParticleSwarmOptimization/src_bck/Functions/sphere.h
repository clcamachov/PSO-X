#ifndef SHIFTED_SPHERE_H
#define SHIFTED_SPHERE_H

#include "../config.h"
#include "../problem.h"
#include "../rng.h"

class Sphere: public Problem {


private:
	double * z;

	long double sphereFunction(int dim, const double* x);
	long double shiftedSphereFunction(int dim, const double* x);
	long double noiseSphereFunction(int dim, const double * x);

public:

	Sphere(Configuration* config, int variantID);
	~Sphere();

	long double evaluate(int dim, const double* x);

};

#endif
