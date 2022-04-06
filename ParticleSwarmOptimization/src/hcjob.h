/*
 * hcjob.h
 *
 *  Created on: May 31, 2019
 *      Author: Christian L. Camacho Villalón
 */

#ifndef HCJOB_H_
#define HCJOB_H_


#include "problem.h"

class HCJob: public Problem {

protected:

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
	double** z;
	double** zM;


	virtual long double basic_func(int func_no, int dim, const double* x);
	long double hybrid_composition(int dim, const double* x);
	void loadIdentityMatrix(double *** M, int w, int h, int d);

//Delegating constructor
public:
	//Delegating constructor
	HCJob(Configuration* config, int variantID): Problem(config, variantID){
			num_func = 0;
			fmax = 0;
			w = 0;
			C = 0;
			rotation_matrix3D = 0;
			shift_vector2D = 0;
			zM = 0;
			z = 0;
	}
	virtual ~HCJob(){};
};



#endif /* HCJOB_H_ */
