/*
 * config.cpp
 *
 *  Created on: Jun 1, 2018
 *      Author: leonardo
 */

#include "config.h"
#include "gsl/gsl_rng.h"
#include <float.h>
#include <string.h>
#include <iostream>
#include <stdio.h>

/****************************************************************************
 * The Configuration class contains the methods for importing the parameters
 * required for the program and the methods for exporting them.
 * **************************************************************************/
//TCLAP makes it easy to import parameters from a certain format

using namespace std;

void printer(char *name, int value) {
	printf("name: %s\tvalue: %d\n", name, value);
}


Configuration::~Configuration(){}

/* Read arguments from command line */
Configuration::Configuration(){
	Configuration::setDefaultParameters();
}

//Program
//TODO: Update with all the parameters
//TODO: Check the integrity of all the parameters here (values and combination with other parameters)
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
		} else if (strcmp(argv[i], "--inertia") == 0){
			inertia = atof(argv[i+1]);
			i++;
			//cout << "\n inertia has been received \n";
		} else if (strcmp(argv[i], "--phi1") == 0) {
			phi_1 = atof(argv[i+1]);
			i++;
			//cout << "\n phi1 has been received \n";
		} else if (strcmp(argv[i], "--phi2") == 0) {
			phi_2 = atof(argv[i+1]);
			i++;
			//cout << "\n phi2 has been received \n";
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
		} else if (strcmp(argv[i], "--clamped") == 0){
			useVelClamping = true;
			//cout << "\n velocity clamping has been set to true \n";
		} else if (strcmp(argv[i], "--omega1CS") == 0) {
			omega1CS = atoi(argv[i+1]);
			i++;
			//cout << "\n inertia control strategy has been received \n";
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
		} else if (strcmp(argv[i], "--omega3CS") == 0) {
			omega3CS = atoi(argv[i+1]);
			i++;
			//cout << "\n omega3 strategy has been received \n";
		} else if (strcmp(argv[i], "--perturbation1") == 0) {
			perturbation1 = atoi(argv[i+1]);
			i++;
			//cout << "\n perturbation strategy has been received \n";
		} else if (strcmp(argv[i], "--perturbation2") == 0) {
			perturbation2 = atoi(argv[i+1]);
			i++;
			//cout << "\n perturbation strategy has been received \n";
		} else if (strcmp(argv[i], "--randomMatrix") == 0) {
			randomMatrix = atoi(argv[i+1]);
			i++;
			//cout << "\n randomMatrix type has been received \n";
		} else if (strcmp(argv[i], "--DNPP") == 0) {
			distributionNPP = atoi(argv[i+1]);
			i++;
			//cout << "\n operator_q has been received \n";
		} else if (strcmp(argv[i], "--operator_q") == 0) {
			operator_q = atoi(argv[i+1]);
			i++;
			//cout << "\n operator_q has been received \n";
		} else if (strcmp(argv[i], "--vRule") == 0) {
			vRule = atol(argv[i+1]);
			i++;
			//cout << "\n velocity rule has been received \n";
		} else if (strcmp(argv[i], "--help") == 0) {
			Configuration::printUsage();
			return(false);
		}
		else {
			cerr << "\nError: Parameter " << argv[i] << "not recognized.\n";
			return(false);
		}
	}


	//Some velocity update rules have special requirements
	if (vRule == 0 ){ //do not require inertia
		inertia = 1.0;
		omega1CS = 0;
	}

	//The topology schedule should be maximum six times the swarm size, i.e. it goes from n to 6n
	if (tSchedule > 6)
		tSchedule = 6;
	if (tSchedule < 1)
		tSchedule = 1;
	tSchedule = particles*tSchedule;

	if (branching < 2)
		branching = 2;
	if (branching > particles)
		branching = floor(particles/2);

	//Topology and model of influence check
	if (topology == TOP_HIERARCHICAL)
		modelOfInfluence = MOI_HIERARCHICAL;

	if (topology == TOP_WHEEL)
		modelOfInfluence = MOI_BEST_OF_N;

	if (modelOfInfluence == MOI_HIERARCHICAL && topology != TOP_HIERARCHICAL)
		modelOfInfluence = MOI_BEST_OF_N; //default

	//The inertia weight schedule
	if (iwSchedule > 4)
		iwSchedule = 4;
	if (iwSchedule < 1)
		iwSchedule = 1;
	iwSchedule = iwSchedule*pow(particles,2);

	//Check problem dimensions
	if (problemDimension < 2 || problemDimension > 50) {
		cerr << "\nError: Dimension should be between 2 and 50.\n";
		return(false);
	}
	//--competition 0 --problem 20 --dimensions 30 --seed 123456 --particles 20 --iterations 2003 --evaluations 20000 --inertiaCS 12 --iwSchedule 2 --initialIW 0.4 --finalIW 0.9 --inertia 0.75 --phi1 0.5 --phi2 0.5 --ring
	//--competition 2 --problem 14 --dimensions 60 --seed 34591 --evaluations 100000 --iterations 100000 --particles 60  --inertia 0.75 --phi1 0.5 --phi2 0.5 --timevarying --tSchedule 240 --clamped
	//Configuration::printParameters();
	//This is for the termination criterion:
	maxFES = 10000 * problemDimension;
	//max_iterations = maxFES;

	return(true);
}

void Configuration::setStartTime(double stime){
	startTime = stime;
}

double Configuration::getStartTime(){
	return startTime;
}

//TODO: Update with all the parameters
void Configuration::printUsage(){
	cout << "" << endl;
	cout << "PSO-something: A flexible and configurable particle swarm optimization framework" << endl;
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
	competitionID = 1;						//CEC05, CEC14, SOFT_COMPUTING, MIXTURE
	problemID = 18;							//25, 19, 19 and 50, respectively
	problemDimension = 2; 					//dimensions
	rngSeed = 12345;						//seed for random numbers
	particles = 10;							//particles (swarm size)
	maxFES = 1000*problemDimension;			//max function evaluation
	max_iterations = 1000*problemDimension;	//max iterations
	inertia = 0.42;							//inertia weight
	phi_1 = 1.55;							//personal coefficient
	phi_2 = 1.55;							//social coefficient
	topology = TOP_TIMEVARYING;				//topology
	minInitRange = -100;					//lower bound of the function
	maxInitRange = 100;						//upper bound of the function
	useVelClamping = true;					//clamp velocity (step size)
	omega1CS = 0;							//inertia control strategy
	initialIW =  0.9;						//initial inertia value
	finalIW = 0.4;							//final inertia value
	iwSchedule = 2*pow(particles,2);		//inertia weight schedule
	omega2CS = 0;
	omega3CS = 0;
	modelOfInfluence = 0;
	tSchedule = 2*particles;				//topology update schedule
	branching = 4;							//branching degree for the hierarchical topology
	perturbation1 = 0;
	perturbation2 = 0;
	randomMatrix = 0;
	distributionNPP = 0;
	operator_q = 0;
	vRule = 1;

	//When the maxInitRange and minInitRange are different from 100
	//the range is updated after instantiating the problem.
	//Also for velocity clamping the bound depends on the function bounds
	//the maxVelLimit and minVelLimit are updated after instantiating
	//the problem.
}

/*Print parameters */
//TODO: Update with all the parameters
void Configuration::printParameters(){

	cout << "\nPSO parameters:\n";
	//cout	<< "  competition:     " << getCompetitionID() << "\n"
	switch (getCompetitionID()){
	case CEC05: 			cout	<< "  competition:       CEC05\n"; break;
	case CEC14: 			cout	<< "  competition:       CEC14\n"; break;
	case SOFT_COMPUTING:	cout	<< "  competition:       SOFT_COMPUTING\n";	break;
	case MIXTURE:   		cout	<< "  competition:       MIXTURE\n"; break;
	}
	cout    << "  problem:           " << getProblemID() << "\n"
			<< "  dimensions:        " << getProblemDimension() << "\n";
	//		<< "  minInitRange:      " << getMinInitBound() << "\n"
	//		<< "  maxInitRange:      " << getMaxInitBound() << "\n"
	//		<< "  seed:              " << getRNGSeed() << "\n"
	cout	<< "  evaluations:       " << getMaxFES() << "\n"
			<< "  iterations:        " << getMaxIterations() << "\n"
			<< "  particles:         " << getSwarmSize() << "\n";
	//		<< "  phi_1:             " << getPhi1() << "\n"
	//		<< "  phi_2:             " << getPhi2() << "\n"
	//		<< "  topology:          " << getTopology() << "\n"
	switch (getTopology()){
	case TOP_FULLYCONNECTED:	cout	<< "  topology:          TOP_FULLYCONNECTED\n";	break;
	case TOP_RING: 				cout	<< "  topology:          TOP_RING\n"; break;
	case TOP_WHEEL:				cout	<< "  topology:          TOP_WHEEL\n"; break;
	case TOP_RANDOM:			cout	<< "  topology:          TOP_RANDOM\n"; break;
	case TOP_VONNEUMANN:		cout	<< "  topology:          TOP_VONNEUMANN\n"; break;
	case TOP_TIMEVARYING:		cout	<< "  topology:          TOP_TIMEVARYING\n"
			<< "  tSchedule          " << getTopologySchedule() << "\n"; break;
	case TOP_HIERARCHICAL:		cout	<< "  topology:          TOP_HIERARCHICAL\n"
			<< "  branching          " << getBranchingDegree() << "\n"; break;
	}
	//<< "  useVelClamping:  " << useVelocityClamping() << "\n"
	switch(useVelocityClamping()){
	case true:		cout	<< "  velocity clamped:  YES\n"; break;
	case false:		cout	<< "  velocity clamped:  NO\n"; break;
	}
	//cout	<< "  omega1CS:        " << getinertiaCS() << "\n"
	switch(getinertiaCS()){
	case IW_CONSTANT:
		cout	<< "  omega1CS:          CONSTANT\n"
		<< "  inertia:           " << getInertia() << "\n"; break;
	case IW_L_INC:
		cout	<< "  omega1CS:          L_INC\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"
		<< "  iwSchedule         " << getIWSchedule() << "\n"; break;
	case IW_L_DEC:
		cout	<< "  omega1CS:          L_DEC\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"
		<< "  iwSchedule         " << getIWSchedule() << "\n"; break;
	case IW_RANDOM:
		cout	<< "  omega1CS:          RANDOM\n"; break;
	case IW_NONL_DEC:
		cout	<< "  omega1CS:          NONL_DEC\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_NONL_DEC_IMP:
		cout	<< "  omega1CS:          NONL_DEC_IMP\n"; break;
	case IW_NONL_DEC_TIME:
		cout	<< "  omega1CS:          NONL_DEC_TIME\n"; break;
	case IW_CHAOTIC_DEC:
		cout	<< "  omega1CS:          CHAOTIC_DEC\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_EXP_DEC:
		cout	<< "  omega1CS:          EXP_DEC\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_OSCILLATING:
		cout	<< "  omega1CS:          OSCILLATING\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_LOG_DEC:
		cout	<< "  omega1CS:          LOG_DEC\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_SELF_REGULATING:
		cout	<< "  omega1CS:          (adaptive) SELF_REGULATING\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_VELOCITY_BASED:
		cout	<< "  omega1CS:          (adaptive) VELOCITY_BASED\n"
		<< "  initialIW          " << getInitialIW() << "\n"
		<< "  finalIW            " << getFinalIW() << "\n"; break;
	case IW_DOUBLE_EXP:
		cout	<< "  omega1CS:         (adaptive) DOUBLE_EXP\n"
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
		cout	<< "  omega1CS:          (adaptive) CONVERGE_BASED\n"; break;
	}
	//<< "  omega2CS          " << getomega2CS() << "\n"
	switch (getomega2CS()){
	case O2_EQUALS_IW: 	cout	<< "  omega2:            EQUALS_IW\n"; break;
	case O2_ZERO: 		cout	<< "  omega2:            ZERO\n"; break;
	case O2_ONE: 		cout	<< "  omega2:            ONE\n"; break;
	case O2_RANDOM: 	cout	<< "  omega2:            RANDOM\n"; break;
	}
	//<< "  omega3CS          " << getomega3CS() << "\n"
	switch (getomega3CS()){
	case O3_EQUALS_IW: 	cout	<< "  omega3:            EQUALS_IW\n"; break;
	case O3_ZERO:		cout	<< "  omega3:            ZERO\n"; break;
	case O3_ONE:		cout	<< "  omega3:            ONE\n"; break;
	case O3_RANDOM:		cout	<< "  omega3:            RANDOM\n"; break;
	}
	//<< "  modelOfInfluence  " << getModelOfInfluence() << "\n"
	switch (getModelOfInfluence()){
	case MOI_BEST_OF_N: 	cout	<< "  modelOfInfluence:  BEST_OF_N\n"; break;
	case MOI_FI: 			cout	<< "  modelOfInfluence:  FI\n"; break;
	case MOI_RANKED_FI:		cout	<< "  modelOfInfluence:  RANKED_FI\n"; break;
	case MOI_HIERARCHICAL:	cout	<< "  modelOfInfluence:  HIERARCHICAL\n"; break;
	}
	//<< "  perturbation1     " << getPerturbation1() << "\n"
	switch (getPerturbation1()){
	case PERT1_NONE: 			cout	<< "  perturbation1:     NONE\n"; break;
	case PERT1_ADD_RECT: 		cout	<< "  perturbation1:     ADD_RECT\n"; break;
	case PERT1_ADD_NOISY:		cout	<< "  perturbation1:     ADD_NOISY\n"; break;
	case PERT1_DIST_NORMAL:		cout	<< "  perturbation1:     DIST_NORMAL\n"; break;
	case PERT1_DIST_SUCCESS: 	cout	<< "  perturbation1:     DIST_SUCCESS\n"; break;
	}
	//<< "  perturbation2     " << getPerturbation2() << "\n"
	switch (getPerturbation2()){
	case PERT2_NONE: 	 	 cout	<< "  perturbation2:     NONE\n"; break;
	case PERT2_ADD_RECT:	 cout	<< "  perturbation2:     ADD_RECT\n"; break;
	case PERT2_ADD_NOISY:	 cout	<< "  perturbation2:     ADD_NOISY\n"; break;
	}
	//<< "  randomMatrix      " << getRandomMatrix() << "\n"
	switch (getRandomMatrix()){
	case MATRIX_NONE:				cout	<< "  randomMatrix:      NONE\n"; break;
	case MATRIX_DIAGONAL: 			cout	<< "  randomMatrix:      DIAGONAL\n"; break;
	case MATRIX_LINEAR:				cout	<< "  randomMatrix:      LINEAR\n"; break;
	case MATRIX_RRM_EXP_MAP:		cout	<< "  randomMatrix:      RRM_EXP_MAP\n"; break;
	case MATRIX_RRM_EUCLIDEAN_ONE:	cout	<< "  randomMatrix:      RRM_EUCLIDEAN_ONE\n"; break;
	case MATRIX_RRM_EUCLIDEAN_ALL: 	cout	<< "  randomMatrix:      RRM_EUCLIDEAN_ALL\n"; break;
	}
	//<< "  DistributionNPP    " << getDistributionNPP() << "\n"
	switch (getDistributionNPP()){
	case DIST_RECTANGULAR: 		cout	<< "  DistributionNPP:   RECTANGULAR\n"; break;
	case DIST_SPHERICAL: 		cout	<< "  DistributionNPP:   SPHERICAL\n"; break;
	case DIST_MULTISPHERICAL:	cout	<< "  DistributionNPP:   MULTISPHERICAL\n"; break;
	case DIST_ADD_STOCH:		cout	<< "  DistributionNPP:   ADD_STOCH\n"; break;
	}
	//<< "  operator_q        " << getOperator_q() << "\n"
	switch (getOperator_q()){
	case Q_STANDARD: cout	<< "  operator_q:        STANDARD\n"; break;
	case Q_GAUSSIAN: cout	<< "  operator_q:        GAUSSIAN\n"; break;
	case Q_DISCRETE: cout	<< "  operator_q:        DISCRETE\n"; break;
	case Q_NORMAL:	 cout	<< "  operator_q:        NORMAL\n"; break;
	}
	//<< "  vRule             " <<  << "\n"
	switch (getVelocityRule()){
	case VEL_BASIC:					cout	<< "  vRule:             BASIC\n"; break;
	case VEL_STANDARD:				cout	<< "  vRule:             STANDARD\n"; break;
	case VEL_LINEAR:				cout	<< "  vRule:             LINEAR\n"; break;
	case VEL_CONSTRICTED:			cout	<< "  vRule:             CONSTRICTED\n"; break;
	case VEL_GUARAN_CONVERG:		cout	<< "  vRule:             GUARAN_CONVERG\n"; break;
	case VEL_FULLY_INFORMED:		cout	<< "  vRule:             FULLY_INFORMED\n"; break;
	case VEL_LOC_CON_TRANS_INV:		cout	<< "  vRule:             LOC_CON_TRANS_INV\n"; break;
	case VEL_STANDARD2011:			cout	<< "  vRule:             STANDARD2011\n";break;
	case VEL_ROTATION_INV:			cout	<< "  vRule:             ROTATION_INV\n"; break;
	}
	cout << endl;
}

//Problem
unsigned int Configuration::getCompetitionID(){
	return competitionID;
}

unsigned int Configuration::getProblemID(){
	return problemID;
}

unsigned int Configuration::getProblemDimension(){
	return problemDimension;
}

unsigned long Configuration::getRNGSeed(){
	return rngSeed;
}

unsigned int Configuration::getMaxFES(){
	return maxFES;
}

unsigned int Configuration::getMaxIterations(){
	return max_iterations;
}

double Configuration::getMinInitRange(){
	if((getProblemID() == 6 || getProblemID() == 24) && (getCompetitionID() == 0))
		return LDBL_MAX*-1.0;

	return minInitRange;
}

double Configuration::getMaxInitRange(){
	if((getProblemID() == 6 || getProblemID() == 24) && (getCompetitionID() == 0))
		return LDBL_MAX;

	return maxInitRange;
}

void Configuration::setMinInitRange(double lowerlimit) {
	minInitRange = lowerlimit;
}

void Configuration::setMaxInitRange(double upperlimit){
	maxInitRange = upperlimit;
}

double Configuration::getMinInitBound(){
	if((getProblemID() == 6 || getProblemID() == 24) && (getCompetitionID() == 0))
		return LDBL_MAX*-1.0;

	return minInitRange;
}

double Configuration::getMaxInitBound(){
	if((getProblemID() == 6 || getProblemID() == 24) && (getCompetitionID() == 0))
		return LDBL_MAX;

	return maxInitRange;
}


//PSO
bool Configuration::useVelocityClamping(){
	return useVelClamping;
}
long int Configuration::getSwarmSize(){
	return particles;
}
double Configuration::getInertia(){
	return inertia;
}
void Configuration::setInertia(double new_inertia){
	inertia = new_inertia;
}
double Configuration::getInitialIW(){
	return initialIW;
}
double Configuration::getFinalIW(){
	return finalIW;
}
unsigned int Configuration::getIWSchedule(){
	return iwSchedule;
}
short Configuration::getinertiaCS(){
	return omega1CS;
}
short Configuration::getomega2CS(){
	return omega2CS;
}
short Configuration::getomega3CS(){
	return omega3CS;
}
bool Configuration::isVelocityClamped(){
	return useVelClamping;
}
void Configuration::setVelocityClamped(bool clamping){
	useVelClamping = clamping;
}
short Configuration::getModelOfInfluence(){
	return modelOfInfluence;
}
short Configuration::getPerturbation1(){
	return perturbation1;
}
short Configuration::getPerturbation2(){
	return perturbation2;
}
short Configuration::getRandomMatrix(){
	return randomMatrix;
}
short Configuration::getDistributionNPP(){
	return distributionNPP;
}
short Configuration::getOperator_q(){
	return operator_q;
}
short Configuration::getTopology(){
	return topology;
}
unsigned int Configuration::getTopologySchedule(){
	return tSchedule;
}
int Configuration::getBranchingDegree(){
	return branching;
}
double Configuration::getPhi1(){
	return phi_1;
}
double Configuration::getPhi2(){
	return phi_2;
}
void Configuration::setEsteps(unsigned int num_esteps){
	esteps = num_esteps;
}
void Configuration::setTopologyUpdatePeriod(int period){
	topologyUpdatePeriod = period;
}
unsigned int Configuration::getEsteps(){
	return esteps;
}
int Configuration::getTopologyUpdatePeriod(){
	return topologyUpdatePeriod;
}

//One can explore the possibility of changing the velocity rule
//during the execution of the algorithm
void Configuration::setVelocityRule(int rule){
	vRule = rule;
}
int Configuration::getVelocityRule(){
	return vRule;
}
