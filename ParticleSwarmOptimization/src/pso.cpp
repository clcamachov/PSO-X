//============================================================================
// Name        : Particle.cpp
// Author      : Christian Leonardo Camacho Villalón
// Version     :
// Copyright   :
// Description :
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <float.h>
#include <iostream>
#include <vector>

#include "config.h"
#include "problem.h"
#include "particle.h"
#include "swarm.h"
#include "hcjob.h"

//Timers
#include <sys/time.h>
#include <sys/resource.h>

//Random numbers
#include "rng.h"

//Functions
#include "Functions/sphere.h"
#include "Functions/schwefel213.h"
#include "Functions/schwefel221.h"
#include "Functions/schwefel222.h"
#include "Functions/schwefel12.h"
#include "Functions/schwefel26.h"
#include "Functions/rosenbrock.h"
#include "Functions/rastrigin.h"
#include "Functions/griewank.h"
#include "Functions/ackley.h"
#include "Functions/extendedF10.h"
#include "Functions/bohachevsky.h"
#include "Functions/schaffer.h"
#include "Functions/schafferf6.h"
#include "Functions/elliptic.h"
#include "Functions/weierstrass.h"
#include "Functions/griewank_rosenbrock.h"
#include "Functions/schwefel.h"
#include "Functions/hgbat.h"
#include "Functions/happycat.h"
#include "Functions/cigar.h"
#include "Functions/discus.h"
#include "Functions/katsuura.h"
#include "Functions/h1.h"
#include "Functions/h2.h"
#include "Functions/h3.h"
#include "Functions/h4.h"
#include "Functions/h7.h"
#include "Functions/h8.h"
#include "Functions/h9.h"
#include "Functions/h10.h"
#include "Functions/h1cec14.h"
#include "Functions/h2cec14.h"
#include "Functions/h3cec14.h"
#include "Functions/h4cec14.h"
#include "Functions/h5cec14.h"
#include "Functions/h6cec14.h"
#include "Functions/hybrid_composition1.h"
#include "Functions/hybrid_composition2.h"
#include "Functions/hybrid_composition3.h"
#include "Functions/hybrid_composition4.h"

using namespace std;

/* This software has three classes containing:
 * 1) a bunch of options -- Configuration: parameters given by the user
 * 2) a bunch of continuous functions -- Problem: continuous optimization problem
 * 3) a bunch of PSO variants -- Swarm: swarm of particles
 *
 * Extending the components of this framework:
 * 	- More neighborhoods structures (topologies) can be added in class Swarm
 * 	- More velocity update strategies can be added in class Particle
 * 	- To use and/or configure components extending the framework is needed to included
 * 	  parameters, arguments, etc., in class Configure and other referenced classes.
 * */

Configuration* config;
Problem* problem;
Swarm* swarm;

long int iterations=0;			/* counter of iterations */
long int evaluations=0;			/* counter of evaluations */

// Benchmark functions
Problem* initializeProblem() {

	//Available problems
	//ABC-X (CEC05, CEC14 and SOCO)
	if(config->getCompetitionID() == MIXTURE){
		//Problem initialization
		switch(config->getProblemID()){
		/////////UNIMODAL FUNCTIONS //////////////////////////////////////////////////
		case SHIFTED_SPHERE:{
			problem = new Sphere(config, SHIFTED);
		}break;
		case SHIFTED_ROTATED_HIGH_CONDITIONED_ELLIPTIC:{
			problem = new Elliptic(config, SHIFTED_ROTATED_HIGH_CONDITIONED);
		}break;
		case ROTATED_BENT_CIGER:{
			problem = new Cigar(config, SHIFTED_ROTATED);
		}break;
		case ROTATED_DISCUS:{
			problem = new Discus(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_SCHWEFEL221:{
			problem = new Schwefel221(config, SHIFTED);
		}break;
		case SHIFTED_SCHWEFEL12:{
			problem = new Schwefel12(config, SHIFTED);
			config->setMinInitRange(-65.536);	//lower bound of the function
			config->setMaxInitRange(65.536);	//upper bound of the function
		}break;
		case SHIFTED_SCHWEFELS12_NOISE_IN_FITNESS:{
			problem = new Schwefel12(config, NOISE_SHIFTED);
		}break;
		case SHIFTED_SCHWEFEL222:{
			problem = new Schwefel222(config, SHIFTED);
			config->setMinInitRange(-10);	//lower bound of the function
			config->setMaxInitRange(10);	//upper bound of the function
		}break;
		case SHIFTED_EXTENDED_F10:{
			problem = new ExtendedF10(config, SHIFTED);
		}break;
		case SHIFTED_BOHACHEVSKY:{
			problem = new Bohachevsky(config, SHIFTED);
		}break;
		case SHIFTED_SCHAFFER:{
			problem = new Schaffer(config, SHIFTED);
		}break;
		case SCHWEFEL26_GLOBAL_OPTIMUM_ON_BOUNDS:{
			problem = new Schwefel26(config, BASIC_GLOBAL_OPTIMUM_ON_BOUNDS);
		}break;
		///////////////////////////////////////////MULTIMODAL FUNCTIONS
		case SHIFTED_ACKLEY:{
			problem = new Ackley(config, SHIFTED);
			config->setMinInitRange(-32);	//lower bound of the function
			config->setMaxInitRange(32);	//upper bound of the function
		}break;
		case SHIFTED_ROTATED_ACKLEY:{
			problem = new Ackley(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROSENBROCK:{
			problem = new Rosenbrock(config, SHIFTED);
		}break;
		case SHIFTED_ROTATED_ROSENBROCK:{
			problem = new Rosenbrock(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_GRIEWANK:{
			problem = new Griewank(config, SHIFTED);
			config->setMinInitRange(-600);	//lower bound of the function
			config->setMaxInitRange(600);	//upper bound of the function
		}break;
		case SHIFTED_ROTATED_GRIEWANK:{
			problem = new Griewank(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_RASTRIGIN:{
			problem = new Rastrigin(config, SHIFTED);
		}break;
		case SHIFTED_ROTATED_RASTRIGIN:{
			problem = new Rastrigin(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_SCHWEFEL:{
			problem = new Schwefel(config, SHIFTED);
		}break;
		case SHIFTED_ROTATED_SCHWEFEL:{
			problem = new Schwefel(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_WEIERSTRASS:{
			problem = new Weierstrass(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_KATSUURA:{
			problem = new Katsuura(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_HAPPYCAT:{
			problem = new Happycat(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_HGBAT:{
			problem = new Hgbat(config, SHIFTED_ROTATED);
		}break;
		/////////////////////////////////////HYBRID FUNCTIONS
		case H1_SOCO:{
			problem = new Hybrid1(config, BASIC);
		}break;
		case H2_SOCO:{
			problem = new Hybrid2(config, BASIC);
		}break;
		case H3_SOCO:{
			problem = new Hybrid3(config, BASIC);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case H4_SOCO:{
			problem = new Hybrid4(config, SHIFTED);
			config->setMinInitRange(-10);	//lower bound of the function
			config->setMaxInitRange(10);		//upper bound of the function
		}break;
		case H7_SOCO:{
			problem = new Hybrid7(config, BASIC);
		}break;
		case H8_SOCO:{
			problem = new Hybrid8(config, BASIC);
		}break;
		case H9_SOCO:{
			problem = new Hybrid9(config, BASIC);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case H10_SOCO:{
			problem = new Hybrid10(config, SHIFTED);
			config->setMinInitRange(-10);	//lower bound of the function
			config->setMaxInitRange(10);		//upper bound of the function
		}break;
		case H1_CEC14:{
			problem = new Hybrid1CEC14(config, SHIFTED_ROTATED);
		}break;
		case H2_CEC14:{
			problem = new Hybrid2CEC14(config, SHIFTED_ROTATED);
		}break;
		case H3_CEC14:{
			problem = new Hybrid3CEC14(config, SHIFTED_ROTATED);
		}break;
		case H4_CEC14:{
			problem = new Hybrid4CEC14(config, SHIFTED_ROTATED);
		}break;
		case H5_CEC14:{
			problem = new Hybrid5CEC14(config, SHIFTED_ROTATED);
		}break;
		case H6_CEC14:{
			problem = new Hybrid6CEC14(config, SHIFTED_ROTATED);
		}break;
		/////////////////////////// COMPOSITION FUNCTIONS
		case COMPOSITION_F1:{
			problem = new HC_1(config, BASIC);
		}break;
		case COMPOSITION_F2:{
			problem = new HC_1(config, ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case COMPOSITION_F3:{
			problem = new HC_1(config, NOISE_ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case COMPOSITION_F4:{
			problem = new HC_2(config, ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case COMPOSITION_F5:{
			problem = new HC_2(config, ROTATED_WITH_NORROW_BASIN_GLOBAL_OPTIMUM);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case COMPOSITION_F6:{
			problem = new HC_2(config, ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case COMPOSITION_F7:{
			problem = new HC_3(config, ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case COMPOSITION_F8:{
			problem = new HC_3(config, ROTATED_WITH_HIGH_CONDITION_NUMBER_MATRIX);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case COMPOSITION_F9:{
			problem = new HC_3(config, NONCONTINUOUS_ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case COMPOSITION_F10:{
			problem = new HC_4(config, ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		}
	}
	if(config->getCompetitionID() == SOFT_COMPUTING){
		//Problem initialization
		switch(config->getProblemID()){
		case SHIFTED_SPHERE_SOCO:{
			problem = new Sphere(config, SHIFTED);
		}break;
		case SHIFTED_SCHWEFEL221_SOCO:{
			problem = new Schwefel221(config, SHIFTED);
		}break;
		case SHIFTED_ROSENBROCK_SOCO:{
			problem = new Rosenbrock(config, SHIFTED);
		}break;
		case SHIFTED_RASTRIGIN_SOCO:{
			problem = new Rastrigin(config, SHIFTED);
		}break;
		case SHIFTED_GRIEWANK_SOCO:{
			problem = new Griewank(config, SHIFTED);
			config->setMinInitRange(-600);	//lower bound of the function
			config->setMaxInitRange(600);	//upper bound of the function
		}break;
		case SHIFTED_ACKLEY_SOCO:{
			problem = new Ackley(config, SHIFTED);
			config->setMinInitRange(-32);	//lower bound of the function
			config->setMaxInitRange(32);	//upper bound of the function
		}break;
		case SHIFTED_SCHWEFEL222_SOCO:{
			problem = new Schwefel222(config, SHIFTED);
			config->setMinInitRange(-10);	//lower bound of the function
			config->setMaxInitRange(10);	//upper bound of the function
		}break;
		case SHIFTED_SCHWEFEL12_SOCO:{
			problem = new Schwefel12(config, SHIFTED);
			config->setMinInitRange(-65.536);	//lower bound of the function
			config->setMaxInitRange(65.536);	//upper bound of the function
		}break;
		case SHIFTED_EXTENDED_F10_SOCO:{
			problem = new ExtendedF10(config, SHIFTED);
		}break;
		case SHIFTED_BOHACHEVSKY_SOCO:{
			problem = new Bohachevsky(config, SHIFTED);
		}break;
		case SHIFTED_SCHAFFER_SOCO:{
			problem = new Schaffer(config, SHIFTED);
		}break;
		case H1:{
			problem = new Hybrid1(config, BASIC);
		}break;
		case H2:{
			problem = new Hybrid2(config, BASIC);
		}break;
		case H3:{
			problem = new Hybrid3(config, BASIC);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case H4:{
			problem = new Hybrid4(config, SHIFTED);
			config->setMinInitRange(-10);	//lower bound of the function
			config->setMaxInitRange(10);	//upper bound of the function
		}break;
		case H7:{
			problem = new Hybrid7(config, BASIC);
		}break;
		case H8:{
			problem = new Hybrid8(config, BASIC);
		}break;
		case H9:{
			problem = new Hybrid9(config, BASIC);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case H10:{
			problem = new Hybrid10(config, SHIFTED);
			config->setMinInitRange(-10);	//lower bound of the function
			config->setMaxInitRange(10);	//upper bound of the function
		}break;
		}
	}
	if(config->getCompetitionID() == CEC05){
		//Problem initialization
		switch(config->getProblemID()){
		case SHIFTED_SPHERE_CEC05:{
			problem = new Sphere(config, SHIFTED);
		}break;
		case SHIFTED_SCHWEFEL12_CEC05:{
			problem = new Schwefel12(config, SHIFTED);
		}break;
		case SHIFTED_ROTATED_HIGH_CONDITIONED_ELLIPTIC_CEC05:{
			problem = new Elliptic(config, SHIFTED_ROTATED_HIGH_CONDITIONED);
		}break;
		case NOISE_SHIFTED_SCHWEFEL12_CEC05:{
			problem = new Schwefel12(config, NOISE_SHIFTED);
		}break;
		case SCHWEFEL26_GLOBAL_OPTIMUM_ON_BOUNDS_CEC05:{
			problem = new Schwefel26(config, BASIC_GLOBAL_OPTIMUM_ON_BOUNDS);
		}break;
		case SHIFTED_ROSENBROCK_CEC05:{
			problem = new Rosenbrock(config, SHIFTED);
		}break;
		case SHIFTED_ROTATED_GRIEWANK_CEC05:{
			problem = new Griewank(config, ROTATED_WITHOUT_BOUNDS);
			config->setMinInitRange(LDBL_MAX*-1.0);	//lower bound of the function
			config->setMaxInitRange(LDBL_MAX);	//upper bound of the function
		}break;
		case SHIFTED_ROTATED_ACKLEY_GOOB_CEC05:{
			problem = new Ackley(config, SHIFTED_ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS);
		}break;
		case SHIFTED_RASTRIGIN_CEC05:{
			problem = new Rastrigin(config, SHIFTED);
		}break;
		case SHIFTED_ROTATED_RASTRIGIN_CEC05:{
			problem = new Rastrigin(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_WEIERSTRASS_CEC05:{
			problem = new Weierstrass(config, SHIFTED_ROTATED);
		}break;
		case BASIC_SCHWEFEL213_CEC05:{
			problem = new Schwefel213(config, BASIC);
		}break;
		case SHIFTED_EXPANDED_GRIEWANKROSENBROCK_CEC05:{
			problem = new GriewankRosenbrock(config, SHIFTED_EXPANDED);
		}break;
		case SHIFTED_ROTATED_EXPANDED_SCHAFFERF6_CEC05:{
			problem = new SchafferF6(config, SHIFTED_ROTATED_EXPANDED);
		}break;
		case BASIC_HYBRIDCOMPOSITION1:{
			problem = new HC_1(config, BASIC);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case ROTATED_HYBRIDCOMPOSITION1:{
			problem = new HC_1(config, ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case NOISE_ROTATED_HYBRIDCOMPOSITION1:{
			problem = new HC_1(config, NOISE_ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case ROTATED_HYBRIDCOMPOSITION2:{
			problem = new HC_2(config, ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case ROTATED_HYBRIDCOMPOSITION2_NBGO:{
			problem = new HC_2(config, ROTATED_WITH_NORROW_BASIN_GLOBAL_OPTIMUM);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case ROTATED_HYBRIDCOMPOSITION2_GOOB:{
			problem = new HC_2(config, ROTATED_GLOBAL_OPTIMUM_ON_BOUNDS);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case ROTATED_HYBRIDCOMPOSITION3:{
			problem = new HC_3(config, ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case ROTATED_HYBRIDCOMPOSITION3_HCNM:{
			problem = new HC_3(config, ROTATED_WITH_HIGH_CONDITION_NUMBER_MATRIX);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case NONCONTINUOUS_ROTATED_HYBRIDCOMPOSITION3:{
			problem = new HC_3(config, NONCONTINUOUS_ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case ROTATED_HYBRIDCOMPOSITION4:{
			problem = new HC_4(config, ROTATED);
			config->setMinInitRange(-5);	//lower bound of the function
			config->setMaxInitRange(5);		//upper bound of the function
		}break;
		case ROTATED_HYBRIDCOMPOSITION4_NO_BOUNDS:{
			problem = new HC_4(config, ROTATED_WITHOUT_BOUNDS);
			config->setMinInitRange(LDBL_MAX*-1.0);	//lower bound of the function
			config->setMaxInitRange(LDBL_MAX);	//upper bound of the function
		}break;
		}
	}
	if(config->getCompetitionID() == CEC14){
		//Problem initialization
		switch(config->getProblemID()){
		case SHIFTED_ROTATED_HIGH_CONDITIONED_ELLIPTIC_CEC14:{
			problem = new Elliptic(config, SHIFTED_ROTATED_HIGH_CONDITIONED);
		}break;
		case ROTATED_BENT_CIGER_CEC14:{
			problem = new Cigar(config, SHIFTED_ROTATED);
		}break;
		case ROTATED_DISCUS_CEC14:{
			problem = new Discus(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_ACKLEY_CEC14:{
			problem = new Ackley(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_ROSENBROCK_CEC14:{
			problem = new Rosenbrock(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_GRIEWANK_CEC14:{
			problem = new Griewank(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_RASTRIGIN_CEC14:{
			problem = new Rastrigin(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_SCHWEFEL_CEC14:{
			problem = new Schwefel(config, SHIFTED);
		}break;
		case SHIFTED_ROTATED_SCHWEFEL_CEC14:{
			problem = new Schwefel(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_WEIERSTRASS_CEC14:{
			problem = new Weierstrass(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_KATSUURA_CEC14:{
			problem = new Katsuura(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_HAPPYCAT_CEC14:{
			problem = new Happycat(config, SHIFTED_ROTATED);
		}break;
		case SHIFTED_ROTATED_HGBAT_CEC14:{
			problem = new Hgbat(config, SHIFTED_ROTATED);
		}break;
		case H1_CEC14_CEC14:{
			problem = new Hybrid1CEC14(config, SHIFTED_ROTATED);
		}break;
		case H2_CEC14_CEC14:{
			problem = new Hybrid2CEC14(config, SHIFTED_ROTATED);
		}break;
		case H3_CEC14_CEC14:{
			problem = new Hybrid3CEC14(config, SHIFTED_ROTATED);
		}break;
		case H4_CEC14_CEC14:{
			problem = new Hybrid4CEC14(config, SHIFTED_ROTATED);
		}break;
		case H5_CEC14_CEC14:{
			problem = new Hybrid5CEC14(config, SHIFTED_ROTATED);
		}break;
		case H6_CEC14_CEC14:{
			problem = new Hybrid6CEC14(config, SHIFTED_ROTATED);
		}break;
		}
	}

	return problem;
}

/* Termination condition */
bool terminationCondition() {
	//Budget based on iterations
	if (config->getMaxIterations() != 0 && iterations > config->getMaxIterations()-1) {
		return true;
	}
	//Budget based on evaluations of the objective function
	if (config->getMaxFES() != 0 && evaluations > config->getMaxFES()) {
		return true;
	}
	return false;
}

/*Free memory used*/
void freeMemory(){
	//Memory release
	RNG::deallocateRNG();
	RNG::deallocatePermutation();
	delete problem;
	delete config;
	//~Swarm() is called automatically when we use hierarchical topology
	if (!swarm->isHierarchical())
		delete swarm;
}

int main(int argc, char *argv[] ){

	double stime;
	static struct timeval tp;

	config = new Configuration();
	if(!config->getConfig(argc, argv)){
		exit(-1);
	}

	//Random number generator
	RNG::initializeRNG(config->getRNGSeed());

	initializeProblem();
	//problem->printProblem();

	swarm = new Swarm(problem, config);
	config->printParameters();

	//Start timer
	gettimeofday( &tp, NULL );
	stime =(double) tp.tv_sec + (double) tp.tv_usec/1000000.0;
	config->setStartTime(stime);

	//Iterations loop
	//cout << "\nMoving swarm...\n"<< endl;
	cout << "\n";
	while(!terminationCondition()){
		iterations++;

		if (config->getTopology() == TOP_TIMEVARYING) {
			swarm->updateTimeVaryingTopology(config, iterations) ;
		}

		swarm->moveSwarm(config, iterations, config->getMinInitBound(),config->getMaxInitBound());
		evaluations=evaluations + config->getSwarmSize();
	}

	cout << "\nBest value found:\t" << scientific << swarm->getGlobalBest().eval << endl;
	swarm->printGbest(config->getProblemDimension());
	problem->printProblem_results();
	freeMemory();   // Free memory.
}
