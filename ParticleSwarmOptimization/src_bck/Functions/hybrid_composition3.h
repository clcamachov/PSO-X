#ifndef HYBRID_COMPOSITION_3_H
#define HYBRID_COMPOSITION_3_H

#include "../problem.h"

#define NUM_OF_FUNC 10

class HC_3: public Problem{

public:
	HC_3(Configuration* config , int variantID);
	~HC_3();

	long double evaluate(int dim, const double * x);

protected:
	long double basic_func(int func_no, int dim, const double* x);
	long double hybridCompositionFunction3(int dim, const double* x);
	long double noncontinuousHybridCompositionFunction3(int dim, const double* x);
	long double hybrid_composition(int dim, const double* x);

	// Number of basic functions
	int num_func;

	// Predefined constant
	double C;
	// Coverage range for each basic function
	static double sigma[];
	// Biases for each basic function
	static double biases[];
	// Stretch / compress each basic function
	static double lambda[];
	// Estimated fmax
	double* fmax;
	// Shift global optimum for each basic function
	double** shift_vector2D;

	// Linear transformation matrix for each basic function
	double*** rotation_matrix3D;

	// Working areas to avoid memory allocation operations
	double* w;
	double** z2d;
	double** zM;


};


#endif
