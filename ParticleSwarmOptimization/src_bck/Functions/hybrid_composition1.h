#ifndef HYBRID_COMPOSITION_1_H
#define HYBRID_COMPOSITION_1_H

#include "../problem.h"

#define NUM_OF_FUNC 10

class HC_1: public Problem{

public:
	HC_1(Configuration* config, int variantID);
	~HC_1();

	long double evaluate(int dim, const double* x);
//    void printShift_vector2D();

protected:
	long double hybridCompositionFunction1(int dim, const double* x);
	long double noiseHybridCompositionFunction1(int dim, const double * x);

	long double basic_func(int func_no, int dim, const double* x);
	long double hybrid_composition(int dim, const double* x);
	void loadIdentityMatrix(double *** M, int w, int h, int d);

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

	double * testP;
	double * testPointM;

};

#endif
