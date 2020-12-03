/*
 * rng.cpp
 *
 *  Created on: May 31, 2019
 *      Author: christian
 */

#include "rng.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"


using namespace std;

gsl_rng* RNG::R = gsl_rng_alloc(gsl_rng_env_setup());
gsl_permutation* RNG::P = gsl_permutation_alloc(1);

void RNG::initializeRNG(unsigned long int rngSeed){
    gsl_rng_set(RNG::R, rngSeed);
}

void RNG::initializePermutation(unsigned int size){
	//Deallocate default
	gsl_permutation_free( RNG::P );

	//Allocate correct one
	RNG::P = gsl_permutation_alloc(size);
	gsl_permutation_init (RNG::P);
}

double RNG::randVal(double min,double max){
	return (gsl_ran_flat( RNG::R, min, max));
}

double RNG::randGauss(double sd){ //mean = 0
	return (gsl_ran_gaussian( RNG::R, sd));
}

double RNG::randGaussWithMean(double sd, double mean){
	return (gsl_ran_gaussian( RNG::R, sd) + mean);
}

double RNG::randCauchy(double gamma){
	return (gsl_ran_cauchy(RNG::R, gamma));
}

double RNG::randBeta(double alpha, double beta){
	return (gsl_ran_beta(RNG::R, alpha, beta));
}

double RNG::randLevy(const double c, const double alpha){
	return (gsl_ran_levy(RNG::R, c, alpha));
	//alpha = 1 = Cauchy distribution
	//alpha = 2 = Gaussian distribution
}

double RNG::randLevySkew(const double c, const double alpha, const double beta){
	return (gsl_ran_levy_skew(RNG::R, c, alpha, beta));
	//alpha = 1 = Cauchy distribution
	//alpha = 2 = Gaussian distribution
}

unsigned int RNG::randBernoulli(double p) {
	return (gsl_ran_bernoulli (RNG::R, p));
}

void RNG::randHypersphere( double * G){
	gsl_ran_dir_nd (RNG::R, sizeof(size_t), G);
}

void RNG::shufflePermutation(){
	gsl_ran_shuffle (RNG::R, RNG::P->data, gsl_permutation_size(RNG::P), sizeof(size_t)); // @suppress("Invalid arguments")
}

int RNG::getPermutationElement(int index){
	return (gsl_permutation_get (RNG::P,index));
}

void RNG::deallocateRNG(){
	 gsl_rng_free( RNG::R );
}

void RNG::deallocatePermutation(){
	 gsl_permutation_free( RNG::P );
}
