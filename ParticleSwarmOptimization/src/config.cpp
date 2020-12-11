/*
 * config.cpp
 *
 *  Created on: May 31, 2019
 *      Author: christian
 */

#include "config.h"
#include "gsl/gsl_rng.h"
#include <float.h>
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <time.h>       /* time */

/****************************************************************************
 * The Configuration class contains the methods for importing the parameters
 * required for the program and the methods for exporting them.
 * **************************************************************************/
//TCLAP makes it easy to import parameters from a certain format

using namespace std;

Configuration::~Configuration(){}

/* Read arguments from command line */
Configuration::Configuration(){
	Configuration::setDefaultParameters();
}

bool Configuration::getConfig(int argc, char *argv[]){
	for(int i=1; i<argc; i++){
		if(strcmp(argv[i], "--competition") == 0){
			competitionID = atoi(argv[i+1]);
			i++;
			//cout << "\n competition has been received \n";
		} else if(strcmp(argv[i], "--problem") == 0){
			problemID = atoi(argv[i+1]);
			i++;
			//cout << "\n problem has been received \n";
		} else if(strcmp(argv[i], "--dimensions") == 0){
			problemDimension = atoi(argv[i+1]);
			i++;
			//cout << "\n number of dimensions has been received \n";
		} else if (strcmp(argv[i], "--seed") == 0) {
			rngSeed = atol(argv[i+1]);
			i++;
			//cout << "\n random seed has been received \n";
		} else if (strcmp(argv[i], "--evaluations") == 0){
			maxFES = atoi(argv[i+1]);
			i++;
			//cout << "\n number of evaluations has been received \n";
		} else if (strcmp(argv[i], "--iterations") == 0){
			max_iterations = atol(argv[i+1]);
			i++;
			//cout << "\n number of iterations has been received \n";
		} else if (strcmp(argv[i], "--particles") == 0){
			particles = atol(argv[i+1]);
			i++;
			//cout << "\n number of particles has been received \n";
		} else if (strcmp(argv[i], "--clamped") == 0){
			useVelClamping = true;
			//cout << "\n velocity clamping has been set to true \n";
		} else if (strcmp(argv[i], "--reinitialized") == 0){
			reinitializePosition = true;
			//cout << "\n reinitializePosition has been set to true \n";
		} else if (strcmp(argv[i], "--perturbed-lb") == 0){
			perturbedlBest = true;
			//cout << "\n use indStrategies has been set to true \n";
		} else if (strcmp(argv[i], "--unstuck") == 0){
			detectStagnation = true;
			//cout << "\n use detectStagnation has been set to true \n";
		} else if (strcmp(argv[i], "--populationCS") == 0) {
			populationCS = atoi(argv[i+1]);
			i++;
			//cout << "\n populationCS has been received \n";
		} else if (strcmp(argv[i], "--initialPopSize") == 0){
			initialPopSize = atol(argv[i+1]);
			i++;
			//cout << "\n initialPopSize has been received \n";
		} else if (strcmp(argv[i], "--finalPopSize") == 0){
			finalPopSize = atol(argv[i+1]);
			i++;
			//cout << "\n finalPopSize has been received \n";
		} else if (strcmp(argv[i], "--particlesToAdd") == 0){
			particlesToAdd = atoi(argv[i+1]);
			i++;
			//cout << "\n finalPopSize has been received \n";
		} else if (strcmp(argv[i], "--popTViterations") == 0){
			popTViterations = atoi(argv[i+1]);
			i++;
			//cout << "\n popTViterations has been received \n";
		} else if (strcmp(argv[i], "--pIntitType") == 0){
			p_intitType = atoi(argv[i+1]);
			i++;
			//cout << "\n finalPopSize has been received \n";
		} else if (strcmp(argv[i], "--phi1") == 0) {
			phi_1 = atof(argv[i+1]);
			i++;
			//cout << "\n phi1 has been received \n";
		} else if (strcmp(argv[i], "--phi2") == 0) {
			phi_2 = atof(argv[i+1]);
			i++;
			//cout << "\n phi2 has been received \n";
		} else if (strcmp(argv[i], "--accelCoeffCS") == 0) {
			accelCoeffCS = atoi(argv[i+1]);
			i++;
			//cout << "\n accelCoeffCS has been received \n";
		} else if (strcmp(argv[i], "--initialPhi1") == 0) {
			initialPhi1 = atof(argv[i+1]);
			i++;
			//cout << "\n initialPhi1 has been received \n";
		} else if (strcmp(argv[i], "--finalPhi1") == 0) {
			finalPhi1 = atof(argv[i+1]);
			i++;
			//cout << "\n finalPhi1 has been received \n";
		} else if (strcmp(argv[i], "--initialPhi2") == 0) {
			initialPhi2 = atof(argv[i+1]);
			i++;
			//cout << "\n initPhi2 has been received \n";
		} else if (strcmp(argv[i], "--finalPhi2") == 0) {
			finalPhi2 = atof(argv[i+1]);
			i++;
			//cout << "\n finalPhi2 has been received \n";
		} else if (strcmp(argv[i], "--topology") == 0){
			topology = atoi(argv[i+1]);
			i++;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--fullyConnected") == 0){
			topology = TOP_FULLYCONNECTED;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--star") == 0){
			topology = TOP_FULLYCONNECTED;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--ring") == 0){
			topology = TOP_RING;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--wheel") == 0){
			topology = TOP_WHEEL;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--random") == 0){
			topology = TOP_RANDOM;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--vonneumann") == 0){
			topology = TOP_VONNEUMANN;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--timevarying") == 0){
			topology = TOP_TIMEVARYING;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--hierarchical") == 0){
			topology = TOP_HIERARCHICAL;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--modInfluence") == 0) {
			modelOfInfluence = atoi(argv[i+1]);
			i++;
			//cout << "\n model of influence has been received \n";
		} else if (strcmp(argv[i], "--tSchedule") == 0) {
			tSchedule = atol(argv[i+1]);
			i++;
			//cout << "\n topology schedule has been received \n";
		} else if (strcmp(argv[i], "--branching") == 0) {
			branching = atoi(argv[i+1]);
			i++;
			//cout << "\n branching degree has been received \n";
		} else if (strcmp(argv[i], "--omega1CS") == 0) {
			omega1CS = atoi(argv[i+1]);
			i++;
			//cout << "\n inertia control strategy has been received \n";
		} else if (strcmp(argv[i], "--inertia") == 0){
			inertia = atof(argv[i+1]);
			i++;
			//cout << "\n inertia has been received \n";
		} else if (strcmp(argv[i], "--initialIW") == 0) {
			initialIW = atof(argv[i+1]);
			i++;
			//cout << "\n initial inertia value has been received \n";
		} else if (strcmp(argv[i], "--finalIW") == 0) {
			finalIW = atof(argv[i+1]);
			i++;
			//cout << "\n final inertia value has been received \n";
		} else if (strcmp(argv[i], "--iwSchedule") == 0) {
			iwSchedule = atol(argv[i+1]);
			i++;
			//cout << "\n inertia weight schedule has been received \n";
		} else if (strcmp(argv[i], "--omega2CS") == 0) {
			omega2CS = atoi(argv[i+1]);
			i++;
			//cout << "\n omega2 strategy has been received \n";
		} else if (strcmp(argv[i], "--omega2") == 0) {
			omega2 = atof(argv[i+1]);
			i++;
			//cout << "\n final inertia value has been received \n";
		} else if (strcmp(argv[i], "--omega3CS") == 0) {
			omega3CS = atoi(argv[i+1]);
			i++;
			//cout << "\n omega3 strategy has been received \n";
		} else if (strcmp(argv[i], "--omega3") == 0) {
			omega3 = atof(argv[i+1]);
			i++;
			//cout << "\n final inertia value has been received \n";
		} else if (strcmp(argv[i], "--perturbation1") == 0) {
			perturbation1CS = atoi(argv[i+1]);
			i++;
			//cout << "\n perturbation value has been received \n";
		} else if (strcmp(argv[i], "--perturbation2") == 0) {
			perturbation2CS = atoi(argv[i+1]);
			i++;
			//cout << "\n perturbation value has been received \n";
		} else if (strcmp(argv[i], "--magnitude1CS") == 0) {
			magnitude1CS = atoi(argv[i+1]);
			i++;
			//cout << "\n magnitude1 strategy has been received \n";
		} else if (strcmp(argv[i], "--magnitude2CS") == 0) {
			magnitude2CS = atoi(argv[i+1]);
			i++;
			//cout << "\n magnitude2 strategy has been received \n";
		} else if (strcmp(argv[i], "--magnitude1") == 0) {
			magnitude1 = atof(argv[i+1]);
			i++;
			//cout << "\n magnitude1 value has been received \n";
		} else if (strcmp(argv[i], "--mag1_parm_l_CS") == 0) {
			mag1_parm_l_CS = atoi(argv[i+1]);
			i++;
			//cout << "\n magnitude1 parameter has been received \n";
		} else if (strcmp(argv[i], "--mag1_par_l") == 0) {
			mag1_parm_l = atof(argv[i+1]);
			i++;
			//cout << "\n magnitude1 parameter has been received \n";
		} else if (strcmp(argv[i], "--mag1_par_m") == 0) {
			mag1_parm_m = atof(argv[i+1]);
			i++;
			//cout << "\n magnitude1 parameter has been received \n";
		} else if (strcmp(argv[i], "--mag1_parm_success") == 0) {
			mag1_parm_success = atoi(argv[i+1]);
			i++;
			//cout << "\n magnitude1 parameter has been received \n";
		} else if (strcmp(argv[i], "--mag1_parm_failure") == 0) {
			mag1_parm_failure = atoi(argv[i+1]);
			i++;
			//cout << "\n magnitude1 parameter has been received \n";
		} else if (strcmp(argv[i], "--magnitude2") == 0) {
			magnitude2 = atof(argv[i+1]);
			i++;
			//cout << "\n magnitude2 value has been received \n";
		} else if (strcmp(argv[i], "--mag2_parm_l_CS") == 0) {
			mag2_parm_l_CS = atoi(argv[i+1]);
			i++;
			//cout << "\n magnitude2 parameter has been received \n";
		} else if (strcmp(argv[i], "--mag2_par_l") == 0) {
			mag2_parm_l = atof(argv[i+1]);
			i++;
			//cout << "\n magnitude2 parameter has been received \n";
		} else if (strcmp(argv[i], "--mag2_par_m") == 0) {
			mag2_parm_m = atof(argv[i+1]);
			i++;
			//cout << "\n magnitude2 parameter has been received \n";
		} else if (strcmp(argv[i], "--mag2_parm_success") == 0) {
			mag2_parm_success = atoi(argv[i+1]);
			i++;
			//cout << "\n magnitude2 parameter has been received \n";
		} else if (strcmp(argv[i], "--mag2_parm_failure") == 0) {
			mag2_parm_failure = atoi(argv[i+1]);
			i++;
			//cout << "\n magnitude2 parameter has been received \n";
		} else if (strcmp(argv[i], "--randomMatrix") == 0) {
			randomMatrix = atoi(argv[i+1]);
			i++;
			//cout << "\n randomMatrix type has been received \n";
		} else if (strcmp(argv[i], "--angleCS") == 0) {
			angleCS = atoi(argv[i+1]);
			i++;
			//cout << "\n operator_q has been received \n";
		} else if (strcmp(argv[i], "--angleSD") == 0) {
			angleSD = atof(argv[i+1]);
			i++;
			//cout << "\n final inertia value has been received \n";
		} else if (strcmp(argv[i], "--angle_par_alpha") == 0) {
			angle_par_alpha = atof(argv[i+1]);
			i++;
			//cout << "\n final inertia value has been received \n";
		} else if (strcmp(argv[i], "--angle_par_beta") == 0) {
			angle_par_beta = atof(argv[i+1]);
			i++;
			//cout << "\n final inertia value has been received \n";
		} else if (strcmp(argv[i], "--rotationAngle") == 0) {
			rotation_angle = atof(argv[i+1]);
			i++;
			//cout << "\n final inertia value has been received \n";
		} else if (strcmp(argv[i], "--DNPP") == 0) {
			distributionNPP = atoi(argv[i+1]);
			i++;
			//cout << "\n operator_q has been received \n";
		} else if (strcmp(argv[i], "--operator_q") == 0) {
			operator_q = atoi(argv[i+1]);
			i++;
			//cout << "\n operator_q has been received \n";
		} else if (strcmp(argv[i], "--operatorCG_parm_r") == 0) {
			operatorCG_parm_r = atof(argv[i+1]);
			i++;
			//cout << "\n operatorCG_parm_r has been received \n";
		} else if (strcmp(argv[i], "--randNeighbor") == 0) {
			if (atoi(argv[i+1]) == 0)
				randNeighbor = false;
			if (atoi(argv[i+1]) == 1)
				randNeighbor = true;
			i++;
			//cout << "\n operator_q will be used with randNeighbor\n";
		} else if (strcmp(argv[i], "--noLogs") == 0) {
			useLogs = false;
			//cout << "\n no logs - suppress logs \n";
		} else if (strcmp(argv[i], "--quiet") == 0) {
			verbose = false;
			//cout << "\n no logs - suppress logs \n";
		} else if (strcmp(argv[i], "--verbose") == 0) {
			verbose = true;
			//cout << "\n no logs - suppress logs \n";
		} else if (strcmp(argv[i], "--output-path") == 0) {
			outputPath = argv[i+1];
			i++;
			//cout << "\n operator_q has been received \n";
		} else if (strcmp(argv[i], "--iw_par_eta") == 0) {
			iw_par_eta = atof(argv[i+1]);
			i++;
			//cout << "\n parameter iw_par_eta for IW_SELF_REGULATING - 11 has been received \n";
		} else if (strcmp(argv[i], "--iw_par_deltaOmega") == 0) {
			iw_par_lambda = atof(argv[i+1]);
			i++;
			//cout << "\n parameter iw_par_deltaOmega for IW_VELOCITY_BASED - 12 has been received \n";
		} else if (strcmp(argv[i], "--iw_par_alpha_2") == 0) {
			iw_par_alpha_2 = atof(argv[i+1]);
			i++;
			//cout << "\n parameter iw_par_alpha_2 for IW_CONVERGE_BASED - 16 has been received \n";
		} else if (strcmp(argv[i], "--iw_par_beta_2") == 0) {
			iw_par_beta_2 = atof(argv[i+1]);
			i++;
			//cout << "\n parameter iw_par_beta_2 for IW_CONVERGE_BASED - 16 has been received \n";
		} else if (strcmp(argv[i], "--help") == 0) {
			Configuration::printUsage();
			return(false);
		}
		else {
			cerr << "\nError: Parameter " << argv[i] << "not recognized.\n";
			return(false);
		}
	}

	//Check the random seed
	if(rngSeed == 0){
		cerr << "\nError: Parameter rngSeed has to be and positive integer, the value received was (" << rngSeed << ")\n";
		return(false);
	}

	//Check problems and competitions range
	if (competitionID == CEC05 && problemID > 24){
		cerr << "\nERROR: the range of problems for CEC05 is from 0 to 24.\n";
		return(false);
	}
	if (competitionID == CEC14 && problemID > 18){
		cerr << "\nERROR: the range of problems for CEC14 is from 0 to 18.\n";
		return(false);
	}
	if (competitionID == SOFT_COMPUTING && problemID > 18){
		cerr << "\nERROR: the range of problems for SOFT COMPUTING is from 0 to 18.\n";
		return(false);
	}
	if (competitionID == MIXTURE && problemID > 49){
		cerr << "\nERROR: the range of problems for MIXTURE is from 0 to 49.\n";
		return(false);
	}

	if (particles == -1)
		particles = problemDimension;

	if (initialPopSize > finalPopSize)
		cerr << "\nError: Wrong initial (or) final population size..." << endl;

	//Check that initialPopSize is at least of the size of the branching degree
	if (topology == TOP_HIERARCHICAL && populationCS != POP_CONSTANT){
		if (branching > (int)initialPopSize)
			initialPopSize = branching;
	}
	if (populationCS == POP_INCREMENTAL)
		particles = initialPopSize;
	if (populationCS == POP_TIME_VARYING){
		if (particles < initialPopSize) particles =  initialPopSize;
		if (particles > finalPopSize) particles =  finalPopSize;
	}

	//Check the bounds of the acceleration coefficients
	//initialPhi1 finalPhi1 initialPhi2 finalPhi2
	// > <
	if (accelCoeffCS == AC_TIME_VARYING){
		//Check order
		if (initialPhi1 < finalPhi1){
			double tmp = initialPhi1;
			initialPhi1 = finalPhi1;
			finalPhi1 = tmp;
			//swap(initialPhi1, finalPhi1);
		}
		if (initialPhi2 > finalPhi2){
			double tmp = initialPhi2;
			initialPhi2 = finalPhi2;
			finalPhi2 = tmp;
			//			swap(initialPhi2,finalPhi2);
		}
		//Check they are different, otherwise use default values
		if (initialPhi1 == finalPhi1){
			initialPhi1 = 2.5;
			finalPhi1 = 0.5;
		}
		if (initialPhi2 == finalPhi2){
			initialPhi2 = 0.5;
			finalPhi2 = 2.5;
		}
	}
	if (accelCoeffCS == AC_RANDOM){
		//In this case initialPhi1 is the lower bound and finalPhi1 the upper bound
		if (finalPhi1 < initialPhi1){
			double tmp = initialPhi1;
			initialPhi1 = finalPhi1;
			finalPhi1 = tmp;
			//			swap(initialPhi1, finalPhi1);
		}
		if (initialPhi1 == finalPhi1){
			initialPhi1 = 0.5;
			finalPhi1 = 2.5;
		}
		initialPhi2 = initialPhi1;
		finalPhi2 = finalPhi1;
	}

	//The topology schedule should be maximum six times the swarm size, i.e. it goes from n to 6n
	if (tSchedule > 10)
		tSchedule = 10;
	if (tSchedule < 1)
		tSchedule = 1;
	if (populationCS == POP_CONSTANT)
		tSchedule = particles*tSchedule;
	else {
		tSchedule = finalPopSize*tSchedule;
		esteps = finalPopSize-3;
	}

	//Check constraints of the Hierarchical Topology
	if (topology == TOP_HIERARCHICAL){
		if( modelOfInfluence == MOI_RANKED_FI )
			modelOfInfluence = MOI_FI;
		if (branching < 2)
			branching = 2;
		if (branching > particles){
			//branching = floor(particles/2);
			cerr << "\nError: The branching degree of a hierarchical topology (" << branching <<  ") should be <= that the population size (" << particles << ").\n";
			return(false);
		}
	}

	//The inertia weight schedule
	if (iwSchedule > 10)
		iwSchedule = 10;
	if (iwSchedule <= 0)
		iwSchedule = 0;
	populationCS == POP_CONSTANT ? iwSchedule = iwSchedule*pow(particles,2) : iwSchedule = iwSchedule*pow(finalPopSize,2);

	if (initialIW > finalIW){
		double tmp = initialIW;
		initialIW = finalIW;
		finalIW = tmp;
	}

	//Check DNPPs
	if (distributionNPP == DIST_ADD_STOCH){
		randomMatrix = MATRIX_NONE;
	}

	//Check valid combinations of Perturbation1 and Magnitude1
	if (perturbation1CS == PERT1_NONE)
		magnitude1CS = MAGNITUDE_NONE;
	else{
		if (magnitude1CS == MAGNITUDE_NONE){
			magnitude1CS = MAGNITUDE_CONSTANT;
			magnitude1 = 0.001;
		}
		if (magnitude1CS == MAGNITUDE_EUC_DISTANCE)
			magnitude1 = 1;
	}
	//Check valid combinations of Perturbation2 and Magnitude2
	if (perturbation2CS == PERT2_NONE){
		omega3CS = O3_ZERO;
		magnitude2CS = MAGNITUDE_NONE;
	}
	else{
		if (magnitude2CS == MAGNITUDE_NONE){
			magnitude2CS = MAGNITUDE_CONSTANT;
			magnitude2 = 0.001;
		}
		if (magnitude2CS == MAGNITUDE_NONE)
			magnitude2 = 1;
	}

	//Check Perturbation magnitude
	if (magnitude1CS == MAGNITUDE_CONSTANT && magnitude1 < 0)
		magnitude1 = 0.001;
	if (magnitude2CS == MAGNITUDE_CONSTANT && magnitude2 < 0)
		magnitude2 = 0.001;


	//Check problem dimensions
	if (problemDimension < 2 || problemDimension > 1000) {
		cerr << "\nError: Dimension should be between 2 and 100.\n";
		return(false);
	}

	if (populationCS == POP_TIME_VARYING) {
		if (popTViterations < 2)
			popTViterations = 2;
		if (popTViterations > 100)
			popTViterations = 100;
	}

	//Termination criteria
	if (maxFES == -1)	//Use a small budget for the number of FEs
		maxFES = 1000 * problemDimension;
	if (maxFES == -2)	//Use the budget of CEC competitions for the number of FEs
		maxFES = 10000 * problemDimension;
	if (maxFES == -3)	//Use half the budget of CEC competitions for the number of FEs
		maxFES = 5000 * problemDimension;

	if (max_iterations < 0)
		populationCS != POP_CONSTANT ? max_iterations = maxFES/initialPopSize : max_iterations = maxFES/particles;

	//Check parameters value
	if(mag1_parm_l_CS == MAG_PARAM_L_INDEPENDENT){
		//scaling factor for MAGNITUDE_EUC_DISTANCE a1=0.911, a2=0.21, a3=0.51, a4=0.58
		mag1_parm_l = 0.911*0.51/(pow((particles/problemDimension),0.21)*pow(problemDimension,0.58));
		if (mag1_parm_l > 0.03) mag1_parm_l = 0.03;
		if (mag1_parm_l < 1E-15) mag1_parm_l = 1E-15;
	}
	if(mag2_parm_l_CS == MAG_PARAM_L_INDEPENDENT){
		//scaling factor for MAGNITUDE_EUC_DISTANCE a1=0.911, a2=0.21, a3=0.51, a4=0.58
		mag2_parm_l = 0.911*0.51/(pow((particles/problemDimension),0.21)*pow(problemDimension,0.58));
		if (mag2_parm_l > 0.03) mag2_parm_l = 0.03;
		if (mag2_parm_l < 1E-15) mag2_parm_l = 1E-15;
	}
	return(true);
}

//TODO: Create complete man page
void Configuration::printUsage(){
	cout << "" << endl;
	cout << "PSO-2020: A flexible and configurable particle swarm optimization framework" << endl;
	cout << "" << endl;
	cout << "General parameters:" << endl;
	cout << "\t--competition <competitionID>" << endl;
	cout << "\t\t We have three benchmark options:" << endl;
	cout << "\t\t CEC05 or 0, CEC14 or 1, SOFT_COMPUTING or 2, MIXTURE or 3." << endl;
	cout << "\t--problem <problemID>" << endl;
	cout << "\t\t For CEC05 can be from 0 to 24, " << endl;
	cout << "\t\t for CEC14 can be from 0 to 18, " << endl;
	cout << "\t\t for SOFT_COMPUTING can be from 0 to 18 and" << endl;
	cout << "\t\t for MIXTURE can be from 0 to 49." << endl;
	cout << "\t--dimensions <problemDimension>" << endl;
	cout << "\t\t To run tests: 10, 30 and 50." << endl;
	cout << "\t\t To tune the algorithm from 10 to 59." << endl;
	cout << "\t--seed <rngSeed>" << endl;
	cout << "\t\t Seed for random numbers." << endl;
	cout << "\t--evaluations <maxFES>" << endl;
	cout << "\t\t Maximum number of function evaluation, for running tests is recommended (10,000 * dimensions)." << endl;
	cout << "\t--iterations <max_iterations>" << endl;
	cout << "\t\t Maximum number of iterations." << endl;
	cout << "\t--particles" << endl;
	cout << "\t\t This is the swarm size, for running test the number of particles recommended are 20, 40 or 60." << endl;
	cout << "\t--inertia" << endl;
	cout << "\t\t Value of the inertia weight use when updating the velocity" << endl;
	cout << "\t--phi1" << endl;
	cout << "\t\t ValuinitialIWe of the personal coefficient" << endl;
	cout << "\t--phi2" << endl;
	cout << "\t\t Value of the social coefficient" << endl;
	cout << "" << endl;
	cout << "Available neighborhood structures:" << endl;
	cout << "\t--fullyConnected" << endl;
	cout << "\t\t TOP_FULLYCONNECTED or 0." << endl;
	cout << "\t--ring" << endl;
	cout << "\t\t TOP_RING or 1." << endl;
	cout << "\t--star" << endl;
	cout << "\t\t TOP_STAR or 2." << endl;
	cout << "\t--random" << endl;
	cout << "\t\t TOP_RANDOM or 3." << endl;
	cout << "\t--vonneumann" << endl;
	cout << "\t\t TOP_VONNEUMANN or 4." << endl;
	cout << "\t--timevarying" << endl;
	cout << "\t\t TOP_TIMEVARYING or 5." << endl;
	cout << "\t--tSchedule" << endl;
	cout << "\t\t When using --timevarying topology is necessary to set the topology update schedule. The value is computed as tSchedule * n, where n=particles." << endl;
	cout << "\t\t The chosen value will be multiplied by n=particles" << endl;
	cout << "\t--scalefree" << endl;
	cout << "\t\t TOP_SCALEFREE or 6." << endl;
	cout << "" << endl;
	cout << "\t--clamped" << endl;
	cout << "\t\t Use velocity (step size) and position clamping." << endl;
	cout << "\t--iwSchedule" << endl;
	cout << "\t\t When using --inertiaCS is necessary to set the the inertia weight schedule. The value used is computed as initialIW * n^2, where n=particles." << endl;
	cout << "\t--initialIW" << endl;
	cout << "\t\t To be described." << endl;
	cout << "\t--finalIW" << endl;
	cout << "\t\t To be described." << endl;
	cout << "\t--vRule" << endl;
	cout << "\t\t To be described." << endl;
}

/*Default parameters (KISS)*/
void Configuration::setDefaultParameters(){
	/** General parameters **/
	/** Problem parameters **/
	competitionID = CEC05;					//CEC05, CEC14, SOFT_COMPUTING, MIXTURE
	problemID = 0;							//25, 19, 19 and 50, respectively
	problemDimension = 10; 					//dimensions
	minInitRange = -100;					//lower bound of the function
	maxInitRange = 100;						//upper bound of the function
	//srand (time(NULL));
	maxFES = 5000*problemDimension;			//max function evaluation
	max_iterations = maxFES;				//max iterations
	useVelClamping = false;					//clamp velocity (step size)
	perturbedlBest = false;
	detectStagnation = false;
	reinitializePosition = false;			//reinitialize particles' position with precision of 10^-5

	/** Population **/
	particles = 10;							//particles (swarm size)
	populationCS = POP_CONSTANT;			//population control strategy
	initialPopSize = 2;						//initial population
	finalPopSize= 1000;						//final or maximum number of individuals allowed
	particlesToAdd = 1;						//number of particles added in a non-constant PopCS
	p_intitType = PARTICLE_INIT_RANDOM;      //type of initialization of particles in a non-constant PopCS
	popTViterations = 2;

	/** Acceleration coefficients **/
	accelCoeffCS = AC_CONSTANT;				//acceleration coefficients control strategy
	phi_1 = 1.496180;						//personal coefficient
	phi_2 = 1.496180;						//social coefficient
	initialPhi1 = 2.5;						//initial personal coefficient value
	initialPhi2 = 0.5;						//initial social coefficient value
	finalPhi1 = 0.5;						//final personal coefficient value
	finalPhi2 = 2.5;						//final social coefficient value

	/** Topology parameters **/
	topology = TOP_TIMEVARYING;				//topology
	tSchedule = 2*particles;				//topology update schedule
	branching = 4;							//branching degree for the hierarchical topology
	esteps = 0;

	/** Model of influence **/
	modelOfInfluence = MOI_BEST_OF_N;		//self-explanatory

	/** Inertia control parameters (omega1 in the GVU) **/
	omega1CS = IW_RANDOM;					//inertia control strategy
	inertia = 1.0;							//inertia weight
	initialIW =  0.9;						//initial inertia value
	finalIW = 0.4;							//final inertia value
	iwSchedule = 0;							//inertia weight schedule
	omega2CS = O2_EQUAL_TO_O1;				//omega2 control strategy (see GVU formula)
	omega2 = 1.0;							//omega2 value when used constant
	omega3CS = O3_ZERO;						//omega3 control strategy (see GVU formula)
	omega3 = 1.0;							//omega2 value when used constant
	iw_par_eta = 1;							//from 0.1 to 1 in IW_SELF_REGULATING - 11
	iw_par_lambda = 0.1;					//from 0.1 to 1 small positive constant in IW_VELOCITY_BASED - 12
	iw_par_alpha_2 = 0.5;					//from 0 to  1 in IW_CONVERGE_BASED - 16
	iw_par_beta_2 = 0.5;					//from 0 to  1 in IW_CONVERGE_BASED - 16

	/** Perturbation and Magnitude**/
	perturbation1CS = PERT1_NONE;		//informed perturbation
	perturbation2CS = PERT2_NONE;		//random (additive) perturbation
	magnitude1CS = MAGNITUDE_CONSTANT;	//strategy to compute the magnitude of the "informed" perturbation
	magnitude2CS = MAGNITUDE_CONSTANT;	//strategy to compute the magnitude of the "informed" perturbation

	//Magnitude 1
	magnitude1 = 0.01;
	mag1_parm_l_CS = MAG_PARAM_L_INDEPENDENT;
	mag1_parm_l = 0.085;
	mag1_parm_m = 0.085;
	mag1_parm_success = 15;
	mag1_parm_failure = 5;

	//Magnitude 2
	magnitude2 = 0.01;
	mag2_parm_l_CS = MAG_PARAM_L_INDEPENDENT;
	mag2_parm_l = 0.085;
	mag2_parm_m = 0.085;
	mag2_parm_success = 15;
	mag2_parm_failure = 5;

	//Magnitude1 and Magnitude2 variables
	mag1_sc = 0; mag1_fc = 0;	//success and failure counters
	mag2_sc = 0; mag2_fc = 0;	//success and failure counters

	/** Matrix **/
	randomMatrix = MATRIX_NONE;			//random matrix
	angleCS = ANGLE_NORMAL;
	rotation_angle = 5;						//rotation angle of RRMs
	angleSD = 20;							//standard deviation of the angle
	angle_par_alpha = 30;
	angle_par_beta = 0.01;

	/** DNPP **/
	distributionNPP = DIST_RECTANGULAR;		//distribution of next possible positions
	operator_q = Q_STANDARD;				//q_operator in simple dynamics PSO
	randNeighbor = true;					//chose a random neighbor as p2 in operator_q
	operatorCG_parm_r = 0.5;				//probability for the Cauchy distribution

	/** Logs **/
	useLogs = true;							//create a folder an log the execution of the algorithm
	verbose = false;
	outputPath = "../";
}

/*Print parameters */
void Configuration::printParameters(){

	cout << "\nPSO-X parameters:\n";
	switch (getCompetitionID()){
	case CEC05: 			cout	<< "  competition:       CEC05\n"; break;
	case CEC14: 			cout	<< "  competition:       CEC14\n"; break;
	case SOFT_COMPUTING:	cout	<< "  competition:       SOFT_COMPUTING\n";	break;
	case MIXTURE:   		cout	<< "  competition:       MIXTURE\n"; break;
	}
	cout    << "  problem:           " << getProblemID() << "\n"
			<< "  dimensions:        " << getProblemDimension() << "\n"
			<< "  minInitRange:      " << getMinInitBound() << "\n"
			<< "  maxInitRange:      " << getMaxInitBound() << "\n"
			<< "  seed:              " << getRNGSeed() << "\n"
			<< "  evaluations:       " << getMaxFES() << "\n"
			<< "  iterations:        " << getMaxIterations() << "\n"
			<< "  particles:         " << getSwarmSize() << "\n";
	switch(useReinitialization()){
	case true:		cout	<< "  reinitilize_pos:   YES\n"; break;
	case false:		cout	<< "  reinitilize_pos:   NO\n"; break;
	}
	switch(detectParticleStagnated()){
	case true:		cout	<< "  detect_stagnation: YES\n"; break;
	case false:		cout	<< "  detect_stagnation: NO\n"; break;
	}
	switch(usePerturbedlBest()){
	case true:		cout	<< "  perturbedlBest:    YES\n"; break;
	case false:		cout	<< "  perturbedlBest:    NO\n"; break;
	}
	switch (getPopulationCS()){
	case POP_CONSTANT: 		cout	<< "  populationCS:      POP_CONSTANT\n"; break;
	case POP_TIME_VARYING: 		cout	<< "  populationCS:      POP_TIME_VARYING\n"
			<< "  minimumPopSize:    " << initialPopSize << "\n"
			<< "  maximumPopSize:    " << finalPopSize << "\n"
			<< "  popTViterations:   " << getPopTViterations() << "\n"; break;
	case POP_INCREMENTAL:	cout	<< "  populationCS:      POP_INCREMENTAL\n"
			<< "  initialPopSize:    " << initialPopSize << "\n"
			<< "  finalPopSize:      " << finalPopSize << "\n"
			<< "  particlesToAdd:    " << getParticlesToAdd() << "\n"
			<< "  p_intitType:       " << getParticleInitType() << "\n"; break;
	}
	switch (getModelOfInfluence()){
	case MOI_BEST_OF_N: 	cout	<< "  modelOfInfluence:  BEST_OF_NEIGHBORHOOD\n"; break;
	case MOI_FI: 			cout	<< "  modelOfInfluence:  FI\n"; break;
	case MOI_RANKED_FI:		cout	<< "  modelOfInfluence:  RANKED_FI\n"; break;
	}
	switch (getTopology()){
	case TOP_FULLYCONNECTED:	cout	<< "  topology:          TOP_FULLYCONNECTED\n";	break;
	case TOP_RING: 				cout	<< "  topology:          TOP_RING\n"; break;
	case TOP_WHEEL:				cout	<< "  topology:          TOP_WHEEL\n"; break;
	case TOP_RANDOM:			cout	<< "  topology:          TOP_RANDOM\n"; break;
	case TOP_VONNEUMANN:		cout	<< "  topology:          TOP_VONNEUMANN\n"; break;
	case TOP_TIMEVARYING:		cout	<< "  topology:          TOP_TIMEVARYING\n";
	if (getTopologySchedule() <= 3)
		cout << "  tSchedule          " << getIWSchedule() << " (fast)\n";
	else if (getTopologySchedule() > 3 && getTopologySchedule() <= 6)
		cout << "  tSchedule          " << getIWSchedule() << " (medium speed)\n";
	else
		cout << "  tSchedule          " << getTopologySchedule() << " (slow)\n";
	break;
	case TOP_HIERARCHICAL:		cout	<< "  topology:          TOP_HIERARCHICAL\n"
			<< "  branching          " << getBranchingDegree() << "\n"; break;
	}
	switch (getRandNeighbor()){
	case true: cout <<  "  randNeighbor:      YES\n"; break;
	case false: cout << "  randNeighbor:      NO\n"; break;
	}
	switch(useVelocityClamping()){
	case true:		cout	<< "  velocity clamped:  YES\n"; break;
	case false:		cout	<< "  velocity clamped:  NO\n"; break;
	}
	switch(getOmega1CS()){
	case IW_CONSTANT:
		cout	<< "  omega1CS:          CONSTANT\n"
		<< "  inertia:           " << getOmega1() << "\n"; break;
	case IW_L_INC:
		cout	<< "  omega1CS:          LINEAR_INCREASING\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n";
		if (getIWSchedule() != 0){
			if (getIWSchedule() <= 3)
				cout << "  iwSchedule         " << getIWSchedule() << " (fast)\n";
			else if (getIWSchedule() > 3 && getIWSchedule() <= 6)
				cout << "  iwSchedule         " << getIWSchedule() << " (medium speed)\n";
			else
				cout << "  iwSchedule         " << getIWSchedule() << " (slow)\n";
		}
		break;
	case IW_L_DEC:
		cout	<< "  omega1CS:          LINEAR_DECREASING\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n";
		if (getIWSchedule() != 0){
			if (getIWSchedule() <= 3)
				cout << "  iwSchedule         " << getIWSchedule() << " (fast)\n";
			else if (getIWSchedule() > 3 && getIWSchedule() <= 6)
				cout << "  iwSchedule         " << getIWSchedule() << " (medium speed)\n";
			else
				cout << "  iwSchedule         " << getIWSchedule() << " (slow)\n";
		}
		break;
	case IW_RANDOM:
		cout	<< "  omega1CS:          RANDOM\n"; break;
	case IW_NONL_DEC:
		cout	<< "  omega1CS:          NON_LINEAR_DECREASING\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_NONL_DEC_IMP:
		cout	<< "  omega1CS:          NON_LINEAR_DECREASING_IMPROVED\n"; break;
	case IW_NONL_DEC_TIME:
		cout	<< "  omega1CS:          NON_LINEAR_DECREASING_TIME\n"; break;
	case IW_CHAOTIC_DEC:
		cout	<< "  omega1CS:          CHAOTIC_DECREASING\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_EXP_DEC:
		cout	<< "  omega1CS:          EXPONENTIAL_DECREASING\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_OSCILLATING:
		cout	<< "  omega1CS:          OSCILLATING\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_LOG_DEC:
		cout	<< "  omega1CS:          LOGARITHMIC_DECREASING\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_SELF_REGULATING:
		cout	<< "  omega1CS:          (adaptive) SELF_REGULATING\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  iw_par_eta         " << get_iw_par_eta() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_VELOCITY_BASED:
		cout	<< "  omega1CS:          (adaptive) VELOCITY_BASED\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  iw_par_deltaOmega  " << get_iw_par_deltaOmega() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_DOUBLE_EXP:
		cout	<< "  omega1CS:         (adaptive) DOUBLE_EXPONENTIAL\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_RANKS_BASED:
		cout	<< "  omega1CS:          (adaptive) RANKS_BASED\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_SUCCESS_BASED:
		cout	<< "  omega1CS:          (adaptive) SUCCESS_BASED\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n";  break;
	case IW_CONVERGE_BASED:
		cout	<< "  omega1CS:          (adaptive) CONVERGE_BASED\n"
		<< "  iw_par_alpha_2     " << get_iw_par_alpha_2() << "\n"
		<< "  iw_par_beta_2      " << get_iw_par_beta_2() << "\n"; break;
	}
	switch (getOmega2CS()){
	case O2_EQUAL_TO_O1: 	cout	<< "  omega2CS:          EQUAL_TO_omega1\n"; break;
	case O2_ZERO: 			cout	<< "  omega2CS:          ZERO (no DNNP)\n"; break;
	case O2_ONE: 			cout	<< "  omega2CS:          ONE\n"; break;
	case O2_RANDOM: 		cout	<< "  omega2CS:          RANDOM\n"; break;
	case O2_CONSTANT: 		cout	<< "  omega2CS:          CONSTANT\n"
			<< "  omega2:            " << getOmega2() << "\n"; break;
	}
	switch (getOmega3CS()){
	case O3_EQUAL_TO_O1: 	cout	<< "  omega3CS:          EQUAL_TO_omega1\n"; break;
	case O3_ZERO:			cout	<< "  omega3CS:          ZERO (no additive perturbation)\n"; break;
	case O3_ONE:			cout	<< "  omega3CS:          ONE\n"; break;
	case O3_RANDOM:			cout	<< "  omega3CS:          RANDOM\n"; break;
	case O3_CONSTANT: 		cout	<< "  omega3CS:          CONSTANT\n"
			<< "  omega3:            " << getOmega3() << "\n"; break;
	}
	switch (getAccelCoeffCS()){
	case AC_CONSTANT:		cout << "  accelCoeffCS:      CONSTANT\n"
			<< "  phi_1:             " << getPhi1() << "\n"
			<< "  phi_2:             " << getPhi2() << "\n"; break;
	case AC_TIME_VARYING: 	cout << "  accelCoeffCS:      TIME_VARYING\n"
			<< "  initialPhi1:       " << getInitialPhi1() << "\n"
			<< "  finalPhi1:         " << getFinalPhi1() << "\n"
			<< "  initialPhi2:       " << getInitialPhi2() << "\n"
			<< "  finalPhi2:         " << getFinalPhi2() << "\n"; break;
	case AC_EXTRAPOLATED:	cout << "  accelCoeffCS:      EXTRAPOLATED\n"; break;
	case AC_RANDOM:			cout << "  accelCoeffCS:      RANDOM\n"
			<< "  initialPhi1:       " << getInitialPhi1() << "\n"
			<< "  finalPhi1:         " << getFinalPhi1() << "\n"
			<< "  initialPhi2:       " << getInitialPhi2() << "\n"
			<< "  finalPhi2:         " << getFinalPhi2() << "\n"; break;
	}
	switch (getPerturbation1CS()){
	case PERT1_NONE: 		cout	<< "  perturbation1CS:   NONE\n"; break;
	case PERT1_GAUSSIAN:	cout	<< "  perturbation1CS:   GAUSSIAN\n"; break;
	case PERT1_LEVY: 		cout	<< "  perturbation1CS:   CAUCHY\n"; break;
	case PERT1_UNIFORM: 	cout	<< "  perturbation1CS:   UNIFORM\n"; break;
	}
	switch (getMagnitude1CS()){
	case MAGNITUDE_CONSTANT: 		cout 	<< "  magnitude1CS:      CONSTANT\n"
			<< "  magnitude1:        " << getMagnitude1() << "\n"; break;
	case MAGNITUDE_EUC_DISTANCE:	cout 	<< "  magnitude1CS:      EUCLIDEAN_DISTANCE\n";
	switch (getMagnitude1_parm_l_CS()){
	case MAG_PARAM_L_INDEPENDENT: 	cout << "  mag1_parm_l_CS:    INDEPENDENT\n"; break;
	case MAG_PARAM_L_USER_SUPPLIED:	cout << "  mag1_parm_l_CS:    USER_SUPPLIED\n"; break;
	}
	cout << "  mag1_parm_l:       " << getMag1_parm_l() << "\n"; break;
	case MAGNITUDE_OBJ_F_DISTANCE:	cout 	<< "  magnitude1CS:      OBJECTIVE_F_DISTANCE\n"
			<< "  mag1_parm_m:       " << getMag1_parm_m() << "\n"; break;
	case MAGNITUDE_SUCCESS:			cout 	<< "  magnitude1CS:      SUCCESS\n"
			<< "  magnitude1:        " << getMagnitude1() << "\n"
			<< "  mag1_parm_success: " << getMag1_parm_success() << "\n"
			<< "  mag1_parm_failure: " << getMag1_parm_failure() << "\n"; break;
	}
	switch (getPerturbation2CS()){
	case PERT2_NONE: 	 	 cout	<< "  perturbation2CS:   NONE\n"; break;
	case PERT2_RECTANGULAR:	 cout	<< "  perturbation2CS:   RECTANGULAR\n"; break;
	case PERT2_NOISY:	 	 cout	<< "  perturbation2CS:   NOISY\n"; break;
	}
	switch (getMagnitude2CS()){
	case MAGNITUDE_CONSTANT: 		cout 	<< "  magnitude2CS:      CONSTANT\n"
			<< "  magnitude2:        " << getMagnitude2() << "\n"; break;
	case MAGNITUDE_EUC_DISTANCE:	cout 	<< "  magnitude2CS:      EUCLIDEAN_DISTANCE\n";
	switch (getMagnitude2_parm_l_CS()){
	case MAG_PARAM_L_INDEPENDENT: 	cout << "  mag2_parm_l_CS:    INDEPENDENT\n"; break;
	case MAG_PARAM_L_USER_SUPPLIED:cout << "  mag2_parm_l_CS:    USER_SUPPLIED\n"; break;
	}
	cout << "  mag2_parm_l:       " << getMag2_parm_l() << "\n"; break;
	case MAGNITUDE_OBJ_F_DISTANCE:	cout 	<< "  magnitude2CS:      OBJECTIVE_F_DISTANCE\n"
			<< "  mag2_parm_m:       " << getMag2_parm_m() << "\n"; break;
	case MAGNITUDE_SUCCESS:			cout 	<< "  magnitude2CS:      SUCCESS\n"
			<< "  magnitude2:        " << getMagnitude2() << "\n"
			<< "  mag2_parm_success: " << getMag2_parm_success() << "\n"
			<< "  mag2_parm_failure: " << getMag2_parm_failure() << "\n"; break;
	}
	switch (getRandomMatrix()){
	case MATRIX_NONE:				cout	<< "  randomMatrix:      NONE\n"; break;
	case MATRIX_DIAGONAL: 			cout	<< "  randomMatrix:      DIAGONAL\n"; break;
	case MATRIX_LINEAR:				cout	<< "  randomMatrix:      LINEAR\n"; break;
	case MATRIX_RRM_EXP_MAP:		cout	<< "  randomMatrix:      RRM_EXP_MAP\n";
	switch (getAngleCS()){
	case ANGLE_CONSTANT: cout << "  angleCS:           CONSTANT\n"
			<< "  rotation_angle:    " << getRotationAgle() << "\n"; break;
	case ANGLE_NORMAL:   cout << "  angleCS:           NORMAL\n"
			<< "  angleSD:           " << getAngleSD() << "\n"; break;
	case ANGLE_ADAPTIVE: cout << "  angleCS:           ADAPTIVE\n"
			<< "  angle_par_alpha:   " << get_angle_par_alpha() << "\n"
			<< "  angle_par_beta:    " << get_angle_par_beta() << "\n"; break;
	break;
	}
	break;
	case MATRIX_RRM_EUCLIDEAN_ONE:	cout	<< "  randomMatrix:      RRM_EUCLIDEAN_ONE\n";
	switch (getAngleCS()){
	case ANGLE_CONSTANT: cout << "  angleCS:           CONSTANT\n"
			<< "  rotation_angle:    " << getRotationAgle() << "\n"; break;
	case ANGLE_NORMAL:   cout << "  angleCS:           NORMAL\n"
			<< "  angleSD:           " << getAngleSD() << "\n"; break;
	case ANGLE_ADAPTIVE: cout << "  angleCS:           ADAPTIVE\n"
			<< "  angle_par_alpha:   " << get_angle_par_alpha() << "\n"
			<< "  angle_par_beta:    " << get_angle_par_beta() << "\n"; break;
	break;
	}
	break;
	case MATRIX_RRM_EUCLIDEAN_ALL: 	cout	<< "  randomMatrix:      RRM_EUCLIDEAN_ALL\n";
	switch (getAngleCS()){
	case ANGLE_CONSTANT: cout << "  angleCS:           CONSTANT\n"
			<< "  rotation_angle:    " << getRotationAgle() << "\n"; break;
	case ANGLE_NORMAL:   cout << "  angleCS:           NORMAL\n"
			<< "  angleSD:           " << getAngleSD() << "\n"; break;
	case ANGLE_ADAPTIVE: cout << "  angleCS:           ADAPTIVE\n"
			<< "  angle_par_alpha:   " << get_angle_par_alpha() << "\n"
			<< "  angle_par_beta:    " << get_angle_par_beta() << "\n"; break;
	break;
	}
	break;
	}
	switch (getDistributionNPP()){
	case DIST_RECTANGULAR: 		cout	<< "  DistributionNPP:   RECTANGULAR\n"; break;
	case DIST_SPHERICAL: 		cout	<< "  DistributionNPP:   SPHERICAL\n"; break;
	case DIST_ADD_STOCH:		cout	<< "  DistributionNPP:   ADD_STOCH\n";
	switch (getOperator_q()){
	case Q_STANDARD: 		cout	<< "  operator_q:        STANDARD\n"; break;
	case Q_GAUSSIAN: 		cout	<< "  operator_q:        GAUSSIAN\n"; break;
	case Q_DISCRETE_2: 		cout	<< "  operator_q:        DISCRETE_2\n"; break;
	case Q_CAUCHY_NORMAL:	cout	<< "  operator_q:        CAUCHY_NORMAL\n"
			<< "  operatorCG_parm_r: " << getOperatorCG_parm_r() << "\n"; break;
	}
	break;
	}
	cout << endl;
}
bool Configuration::logOutput(){
	return (useLogs);
}
bool Configuration::verboseMode(){
	return (verbose);
}
std::string Configuration::getOutputPath(){
	return (outputPath);
}

//Problem
int Configuration::getCompetitionID(){
	return (competitionID);
}
int Configuration::getProblemID(){
	return (problemID);
}
int Configuration::getProblemDimension(){
	return (problemDimension);
}
unsigned long Configuration::getRNGSeed(){
	return (rngSeed);
}
long int Configuration::getMaxFES(){
	return (maxFES);
}
long int Configuration::getMaxIterations(){
	return (max_iterations);
}
void Configuration::setMinInitRange(double lowerlimit) {
	minInitRange = lowerlimit;
}
double Configuration::getMinInitRange(){
	return (minInitRange);
}
double Configuration::getMinInitBound(){
	return (minInitRange);
}
void Configuration::setMaxInitRange(double upperlimit){
	maxInitRange = upperlimit;
}
double Configuration::getMaxInitRange(){
	return (maxInitRange);
}
double Configuration::getMaxInitBound(){
	return (maxInitRange);
}

//PSO
long int Configuration::getSwarmSize(){
	return (particles);
}
void Configuration::setSwarmSize(long int new_size){
	particles = new_size;
}
int Configuration::getPopulationCS(){
	return (populationCS);
}
long int Configuration::getInitialPopSize(){
	return (initialPopSize);
}
long int Configuration::getFinalPopSize(){
	return (finalPopSize);
}
void Configuration::setParticlesToAdd(int new_pool_size){
	particlesToAdd = new_pool_size;
}
int Configuration::getParticlesToAdd(){
	return (particlesToAdd);
}
int Configuration::getParticleInitType(){
	return (p_intitType);
}
int Configuration::getPopTViterations(){
	return (popTViterations);
}
short Configuration::getOmega1CS(){
	return (omega1CS);
}
short Configuration::getOmega2CS(){
	return (omega2CS);
}
short Configuration::getOmega3CS(){
	return (omega3CS);
}
double Configuration::getOmega1(){
	return (inertia);
}
void Configuration::setOmega1(double new_inertia){
	inertia = new_inertia;
}
double Configuration::getInitialIW(){
	return (initialIW);
}
double Configuration::getFinalIW(){
	return (finalIW);
}
int Configuration::getIWSchedule(){
	return (iwSchedule);
}
double Configuration::get_iw_par_eta(){
	return (iw_par_eta);
}
double Configuration::get_iw_par_deltaOmega(){
	return (iw_par_lambda);
}
double Configuration::get_iw_par_alpha_2(){
	return (iw_par_alpha_2);
}
double Configuration::get_iw_par_beta_2(){
	return (iw_par_beta_2);
}
double Configuration::getOmega2(){
	return (omega2);
}
void Configuration::setOmega2(double new_omega2){
	omega2 = new_omega2;
}
double Configuration::getOmega3(){
	return (omega3);
}
void Configuration::setOmega3(double new_omega3){
	omega3 = new_omega3;
}
bool Configuration::useVelocityClamping(){
	return (useVelClamping);
}
bool Configuration::usePerturbedlBest(){
	return (perturbedlBest);
}
//Position
bool Configuration::useReinitialization(){
	return (reinitializePosition);
}
void Configuration::setReinitialization(bool reinit){
	reinitializePosition = reinit;
}
bool Configuration::detectParticleStagnated(){
	return (detectStagnation);
}
short Configuration::getModelOfInfluence(){
	return (modelOfInfluence);
}
//Perturbation and Magnitudes Control strategies
short Configuration::getPerturbation1CS(){
	return (perturbation1CS);
}
short Configuration::getPerturbation2CS(){
	return (perturbation2CS);
}
short Configuration::getMagnitude1CS(){
	return (magnitude1CS);
}
short Configuration::getMagnitude2CS(){
	return (magnitude2CS);
}
//Parameters Magnitude1
double Configuration::getMagnitude1(){
	return (magnitude1);
}
void Configuration::setMagnitude1(double mag_1){
	magnitude1 = mag_1;
}
short Configuration::getMagnitude1_parm_l_CS(){
	return (mag1_parm_l_CS);
}
double Configuration::getMag1_parm_l(){
	return (mag1_parm_l);
}
double Configuration::getMag1_parm_m(){
	return (mag1_parm_m);
}
int Configuration::getMag1_parm_success(){
	return (mag1_parm_success);
}
void Configuration::setMag1_parm_success(int new_mag1_succ){
	mag1_parm_success = new_mag1_succ;
}
int Configuration::getMag1_parm_failure(){
	return (mag1_parm_failure);
}
void Configuration::setMag1_parm_failure(int new_mag1_fail){
	mag1_parm_failure = new_mag1_fail;
}
//Setters/Getters variables magnitude1
void Configuration::set_mag1_sc(int sc){
	mag1_sc=sc;
}
void Configuration::set_mag1_fc(int fc){
	mag1_fc=fc;
}
int Configuration::get_mag1_sc(){
	return (mag1_sc);
}
int Configuration::get_mag1_fc(){
	return (mag1_fc);
}

//Parameters Magnitude2
double Configuration::getMagnitude2(){
	return (magnitude2);
}
void Configuration::setMagnitude2(double mag_2){
	magnitude2 = mag_2;
}
short Configuration::getMagnitude2_parm_l_CS(){
	return (mag2_parm_l_CS);
}
double Configuration::getMag2_parm_l(){
	return (mag2_parm_l);
}
double Configuration::getMag2_parm_m(){
	return (mag2_parm_m);
}
int Configuration::getMag2_parm_success(){
	return (mag2_parm_success);
}
int Configuration::getMag2_parm_failure(){
	return (mag2_parm_failure);
}

//Setters/Getters variables magnitude2
void Configuration::set_mag2_sc(int sc){
	mag2_sc = sc;
}
void Configuration::set_mag2_fc(int fc){
	mag2_fc = fc;
}
int Configuration::get_mag2_sc(){
	return (mag2_sc);
}
int Configuration::get_mag2_fc(){
	return (mag2_fc);
}
//Random matrices
short Configuration::getRandomMatrix(){
	return (randomMatrix);
}
short Configuration::getAngleCS(){
	return (angleCS);
}
double Configuration::getAngleSD(){
	return (angleSD);
}
void Configuration::setAngleSD(double angle_sd){
	angleSD = angle_sd;
}
double Configuration::get_angle_par_alpha(){
	return (angle_par_alpha);
}
double Configuration::get_angle_par_beta(){
	return (angle_par_beta);
}

double Configuration::getRotationAgle(){
	return (rotation_angle);
}
void Configuration::setRotationAgle(double angle){
	rotation_angle = angle;
}
short Configuration::getDistributionNPP(){
	return (distributionNPP);
}
short Configuration::getOperator_q(){
	return (operator_q);
}
bool Configuration::getRandNeighbor(){
	return (randNeighbor);
}
double Configuration::getOperatorCG_parm_r(){
	return (operatorCG_parm_r);
}
short Configuration::getTopology(){
	return (topology);
}
unsigned int Configuration::getTopologySchedule(){
	return (tSchedule);
}
int Configuration::getBranchingDegree(){
	return (branching);
}
short Configuration::getAccelCoeffCS(){
	return (accelCoeffCS);
}
double Configuration::getPhi1(){
	return (phi_1);
}
double Configuration::getPhi2(){
	return (phi_2);
}
double Configuration::getInitialPhi1(){
	return (initialPhi1);
}
double Configuration::getFinalPhi1(){
	return (finalPhi1);
}
double Configuration::getInitialPhi2(){
	return (initialPhi2);
}
double Configuration::getFinalPhi2(){
	return (finalPhi2);
}
void Configuration::setEsteps(unsigned int num_esteps){
	esteps = num_esteps;
}
unsigned int Configuration::getEsteps(){
	return (esteps);
}
void Configuration::setTopologyUpdatePeriod(int period){
	topologyUpdatePeriod = period;
}
int Configuration::getTopologyUpdatePeriod(){
	return (topologyUpdatePeriod);
}
void Configuration::setStartTime(double stime){
	startTime = stime;
}

double Configuration::getStartTime(){
	return (startTime);
}

