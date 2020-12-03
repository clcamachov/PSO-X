/*
 * rng.h
 *
 *  Created on: May 31, 2019
 *      Author: christian
 */

#ifndef RNG_H_
#define RNG_H_

#include "gsl/gsl_rng.h"
#include "gsl/gsl_permutation.h"

class RNG{
public:
static gsl_rng* R;
static gsl_permutation* P;
static void initializeRNG( unsigned long int rngSeed );
static void initializePermutation(unsigned int size);

static double randVal(double min,double max);
static double randGauss(double sd);
static double randGaussWithMean(double sd, double mean);
static double randCauchy(double gamma);
static double randBeta(double alpha, double beta);
static double randLevy(const double c, const double alpha);
static double randLevySkew(const double c, const double alpha, const double beta);


static unsigned int randBernoulli(double p);

static void randHypersphere( double * G);

static void shufflePermutation();
static int getPermutationElement(int index);

static void deallocateRNG();
static void deallocatePermutation();
};



#endif /* RNG_H_ */
