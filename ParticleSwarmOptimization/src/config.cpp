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

//Program
//TODO: Update with all the parameters
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
		} else if (strcmp(argv[i], "--fullyConnected") == 0){
			topology = TOP_FULLYCONNECTED;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--ring") == 0){
			topology = TOP_RING;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--star") == 0){
			topology = TOP_STAR;
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
		} else if (strcmp(argv[i], "--scalefree") == 0){
			topology = TOP_SCALEFREE;
			//cout << "\n topology has been received \n";
		} else if (strcmp(argv[i], "--clamped") == 0){
			useVelClamping = true;
			//cout << "\n velocity clamping has been set to true \n";
		} else if (strcmp(argv[i], "--inertiaCS") == 0) {
			inertiaCS = atoi(argv[i+1]);
			i++;
			//cout << "\n initial inertia weight has been received \n";
		} else if (strcmp(argv[i], "--initialIW") == 0) {
			initialIW = atof(argv[i+1]);
			i++;
			//cout << "\n initial inertia weight has been received \n";
		} else if (strcmp(argv[i], "--finalIW") == 0) {
			finalIW = atof(argv[i+1]);
			i++;
			//cout << "\n final inertia weight has been received \n";
		} else if (strcmp(argv[i], "--iwSchedule") == 0) {
			iwSchedule = atol(argv[i+1]);
			i++;
			//cout << "\n inertia weight schedule has been received \n";
		} else if (strcmp(argv[i], "--tSchedule") == 0) {
			tSchedule = atol(argv[i+1]);
			i++;
			//cout << "\n topology schedule has been received \n";
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
		inertiaCS = 0;
	}

	//The topology schedule should be maximum six times the swarm size, i.e. it goes from n to 6n
	if (tSchedule > 6)
		tSchedule = 6;
	if (tSchedule < 1)
		tSchedule = 1;
	tSchedule = particles*tSchedule;

	//The inertia weight schedule
	if (iwSchedule > 4)
		iwSchedule = 4;
	if (iwSchedule < 1)
		iwSchedule = 1;
	iwSchedule = iwSchedule*pow(particles,2);

	//Check problem dimensions
	if (problemDimension < 2 || problemDimension > 100) {
		cerr << "\nError: Dimension should be between 2 and 100.\n";
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
	cout << "PSO-X: A flexible and configurable particle swarm optimization framework" << endl;
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
	inertiaCS = 0;							//inertia control strategy
	initialIW =  0.9;						//initial inertia weight
	finalIW = 0.4;							//final inertia weight
	iwSchedule = 2*pow(particles,2);		//inertia weight schedule
	tSchedule = 2*particles;				//topology update schedule
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
	cout << "\nPSO parameters:\n"
			<< "  competition:     " << getCompetitionID() << "\n"
			<< "  problem:         " << getProblemID() << "\n"
			<< "  dimensions:      " << getProblemDimension() << "\n"
			<< "  seed:            " << getRNGSeed() << "\n"
			<< "  evaluations:     " << getMaxFES() << "\n"
			<< "  iterations:      " << getMaxIterations() << "\n"
			<< "  particles:       " << getSwarmSize() << "\n"
			<< "  inertia:         " << getInertia() << "\n"
			<< "  phi_1:           " << getPhi1() << "\n"
			<< "  phi_2:           " << getPhi2() << "\n"
			<< "  topology:        " << getTopology() << "\n"
			<< "  minInitRange:    " << getMinInitBound() << "\n"
			<< "  maxInitRange:    " << getMaxInitBound() << "\n"
			<< "  useVelClamping:  " << useVelocityClamping() << "\n"
			<< "  inertiaCS:       " << getinertiaCS() << "\n"
			<< "  initialIW        " << getInitialIW() << "\n"
			<< "  finalIW          " << getFinalIW() << "\n"
			<< "  iwSchedule       " << getIWSchedule() << "\n"
			<< "  tSchedule        " << getTopologySchedule() << "\n"
			<< "  vRule            " << getVelocityRule() << "\n"
			<< endl;
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
	return inertiaCS;
}

bool Configuration::isVelocityClamped(){
	return useVelClamping;
}
void Configuration::setVelocityClamped(bool clamping){
	useVelClamping = clamping;
}

short Configuration::getTopology(){
	return topology;
}
unsigned int Configuration::getTopologySchedule(){
	return tSchedule;
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
