/*
 * hcjob.cpp
 *
 *  Created on: May 31, 2019
 *      Author: christian
 */

#include "hcjob.h"

#include <cfloat>
#include <math.h>

using namespace std;

//HCJob::~HCJob(){}

// Coverage range for each basic function
double HCJob::sigma[] = {} ;
// Biases for each basic function
double HCJob::biases[] = {} ;
// Stretch / compress each basic function
double HCJob::lambda[] = {} ;

long double HCJob::hybrid_composition(int dim, const double* x) {

	// Get the raw weights
	double wMax = DBL_MIN;

	for (int i = 0 ; i < num_func ; i ++) {
		double sumSqr = 0.0;
		shift(dim, z[i], x, shift_vector2D[i]);
		for (int j = 0 ; j < dim ; j ++) {
			sumSqr += (z[i][j] * z[i][j]);
		}
		w[i] = exp(-1.0 * sumSqr / (2.0 * dim * sigma[i] * sigma[i]));
		if (wMax < w[i])
			wMax = w[i];
	}

	// Modify the weights
	double wSum = 0.0;
	double w1mMaxPow = 1.0 - pow(wMax, 10.0);
	for (int i = 0 ; i < num_func ; i ++) {
		if (w[i] != wMax) {
			w[i] *= w1mMaxPow;
		}
		wSum += w[i];
	}

	// Normalize the weights
	for (int i = 0 ; i < num_func ; i ++) {
		w[i] /= wSum;
	}

	long double sumF = 0.0;
	for (int i = 0 ; i < num_func ; i ++) {
		for (int j = 0 ; j < dim ; j ++) {
			z[i][j] /= lambda[i];
		}

		rotate(dim, zM[i], z[i], rotation_matrix3D[i]);
		sumF +=	w[i] * (C * basic_func(i, dim, zM[i]) / fmax[i] + biases[i]);
	}
	return (sumF);
}

void HCJob::loadIdentityMatrix(double *** M, int w, int h, int d){


	// Generate identity matrices
	for (int i = 0 ; i < w ; i ++) {
		for (int j = 0 ; j < h ; j ++) {
			for (int k = 0 ; k < d ; k ++) {
				M[i][j][k] = 0.0;
			}
		}
		for (int j = 0 ; j < h ; j ++) {
			M[i][j][j] = 1.0;
		}
	}
}

