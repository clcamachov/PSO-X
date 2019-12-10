/*
 * swarm.cpp
 *
 *  Created on: Jun 20, 2018
 *      Author: leonardo
 */
#include "iostream"
#include "float.h"
#include "utils.h"
#include <bits/stdc++.h>
#include <new>

#include "problem.h"
#include "config.h"
#include "swarm.h"

#include "rng.h"

#define ALPHA_T_PRECISION  0.00001

using namespace std;

/* Variables to compute the inertia weight using the available control strategies */
double alpha = 1/pow(M_PI,2); 				//small positive constant 		IW_NONL_DEC - 4
double omega = 0.3;							//value between [0,1] 			IW_NONL_DEC_IMP - 5
double u = 1.0002;							//value between [1.0001,1.0005] IW_NONL_DEC_IMP - 5
double omegaXu = omega * u;					//								IW_NONL_DEC_IMP - 5
double zFunction;							//								IW_CHAOTIC_DEC - 7
int	   k = 7;								//positive integer constant		IW_OSCILLATING - 9
double simNumOfCos = 2.0*M_PI*((4*k)+6);	//								IW_OSCILLATING - 9
double a = 1;								//small positive constant		IW_LOG_DEC - 10
double omega_2 = 0;							//								IW_SELF_REGULATING - 11
double eta = 1;								//								IW_SELF_REGULATING - 11
double idealVelocity;						//								IW_VELOCITY_BASED - 12
double avVel;								//								IW_VELOCITY_BASED - 12
double deltaOmega = 0.1;					//small positive constant		IW_VELOCITY_BASED - 12
struct SimplifySwarm simpSwarm;				//								IW_RANKS_BASED - 14
double alpha_2 = 0.5;						//small constant in [0,1]		IW_CONVERGE_BASED
double beta_2 = 0.5;						//small constant in [0,1]		IW_CONVERGE_BASED
/* Variable for the topology and model of influence */
int** hierarchy;							//tree structure for the hierarchical model of influence
int lastLevelComplete = 0;					//global variable for the hierarchy
int* Informants;							//variable length array that contains the ID of informants of a particle
struct SimplifySwarm rankedSwarm;			//ranked FI model of influence
bool modInfRanked = false;					//flag to indicate that the rankedSwarm structure was used and delete it
/* Variable for the perturbation strategies */
int success = 15, failure = 5;				//success and failure thresholds for the additive rectangular perturbation
int sc = 0, fc = 0;							//success and failure counters
double alpha_t = 1.0;						//side length of the rectangle for the success-rate perturbation
double delta = 1.0;							//side length of the rectangle for the uniform random perturbation
double l = 0.01;							//scaling factor for the perturbation


bool sortcol(const vector<int>& v1, const vector<int>& v2) {
	return v1[1] < v2[1];
}

//Default constructor
Swarm::Swarm(){
	problem = 0;
	size = 0;
	best_particle = 0;
	init = false;
	ranked = false;
	hierarchical = false;
}

Swarm::~Swarm(){
	if (init) {
		// Memory allocated dynamically
		for (long int i=0;i<size;i++)
			delete swarm.at(i);

		delete [] global_best.x;
		delete [] Informants;

		if (ranked){
			//free memory reserved to the use rankings
			delete [] simpSwarm.id;
			delete [] simpSwarm.eval;
		}
		if (hierarchical){
			for (int r=0; r<=lastLevelComplete+1; r++){
				delete [] hierarchy[r];
			}
			delete [] hierarchy;
		}
		if (modInfRanked){
			delete [] rankedSwarm.id;
			delete [] rankedSwarm.eval;
		}
	}
	init=false;
}

Swarm::Swarm (Problem* problem, Configuration* config){
	//cout << "Creating swarm.\n\n";

	this->problem = problem;
	size = config->getSwarmSize();

	/*Initialize global best*/
	global_best.x = new double[config->getProblemDimension()];
	for(unsigned int i=0;i<config->getProblemDimension();i++)
		global_best.x[i] = 0;
	global_best.eval=LDBL_MAX;

	for (long int i=0; i<size; i++) {
		Particle* aParticle = new Particle(problem, config, i);
		swarm.push_back(aParticle);

		if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
			updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
			best_particle = swarm.at(i);
		}
	}

	hierarchical = false;

	//Select one of the available topologies
	if (config->getTopology() == TOP_FULLYCONNECTED) {
		createFullyConnectedTopology();
	} else if (config->getTopology() == TOP_HIERARCHICAL) {
		hierarchical = true;
		createHierarchical(config->getBranchingDegree());
	} else if (config->getTopology() == TOP_RING) {
		createRingTopology();
	} else if (config->getTopology() == TOP_WHEEL) {
		createWheelTopology();
	} else if (config->getTopology() == TOP_RANDOM) {
		createRandomEdge();
	} else if (config->getTopology() == TOP_TIMEVARYING) {
		createFullyConnectedTopology();
		config->setTopologyUpdatePeriod(
				(int)floor((double)config->getTopologySchedule()/(config->getSwarmSize()-3)));
		RNG::initializePermutation(config->getSwarmSize());
	} else if (config->getTopology() == TOP_VONNEUMANN) {
		createVonNeumannTopology();
	}
	else {
		cerr << "Wrong topology" << endl;
		exit (-1);
	}

	//rankings
	if (config->getOmega1CS() == IW_RANKS_BASED || config->getOmega1CS() == IW_SUCCESS_BASED
			|| config->getOmega1CS() == IW_CONVERGE_BASED)
		ranked = true;

	init=true;
}

/* Copy constructor */
Swarm::Swarm (const Swarm &s, Configuration* config){
	problem = s.problem;
	size = s.size;

	/*Initialize global best*/
	if (! init) {
		global_best.x = new double[config->getProblemDimension()];
		for(unsigned int i=0;i<config->getProblemDimension();i++)
			global_best.x[i] = 0;
		global_best.eval = LDBL_MAX;
		for (long int i=0; i<size; i++) {
			Particle* aParticle = new Particle(problem, config, i);
			swarm.push_back(aParticle);

			if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
				updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
				best_particle = swarm.at(i);
			}
		}
		hierarchical = false;

		//Select one of the available topologies
		if (config->getTopology() == TOP_FULLYCONNECTED) {
			createFullyConnectedTopology();
		} else if (config->getTopology() == TOP_HIERARCHICAL) {
			hierarchical = true;
			createHierarchical(config->getBranchingDegree());
		} else if (config->getTopology() == TOP_RING) {
			createRingTopology();
		} else if (config->getTopology() == TOP_WHEEL) {
			createWheelTopology();
		} else if (config->getTopology() == TOP_RANDOM) {
			createRandomEdge();
		} else if (config->getTopology() == TOP_TIMEVARYING) {
			createFullyConnectedTopology();
			config->setTopologyUpdatePeriod(
					(int)floor((double)config->getTopologySchedule()/(config->getSwarmSize()-3)));
			RNG::initializePermutation(config->getSwarmSize());
		} else if (config->getTopology() == TOP_VONNEUMANN) {
			createVonNeumannTopology();
		}
		else {
			cerr << "Wrong topology" << endl;
			exit (-1);
		}

		//rankings
		if (config->getOmega1CS() == IW_RANKS_BASED || config->getOmega1CS() == IW_SUCCESS_BASED
				|| config->getOmega1CS() == IW_CONVERGE_BASED)
			ranked = true;
	}
	else {
		//Copy swarm
		for (long int i=0; i<size; i++) {
			swarm.push_back(s.swarm.at(i));
		}
		//Copy global_best
		for(unsigned int i=0;i<config->getProblemDimension();i++)
			global_best.x[i] = s.global_best.x[i];
		global_best.eval=s.global_best.eval;

		best_particle = s.best_particle;
		hierarchical = s.hierarchical;

		for (long int i=0; i<size; i++) {
			for (unsigned int j= 0; j<swarm.at(i)->neighbours.size(); j++)
				swarm.at(i)->neighbours[j]=s.swarm.at(i)->neighbours[j];
		}
		ranked = s.ranked;
	}
	init = true;
}

/*Update global best solution found */
void Swarm::updateGlobalBest(double* new_x, double eval){
	for (int i=0;i<problem->getProblemDimension();i++){
		global_best.x[i]=new_x[i];
	}
	global_best.eval=eval;
}

void Swarm::printGbest(unsigned int dimensions){
	//print best solution
	cout << "[ " ;
	//	for(unsigned int i=0; i< sizeof(getGlobalBest().x); i++){
	for(unsigned int i=0; i< dimensions; i++){
		cout << getGlobalBest().x[i] << "  ";
	}
	cout << " ]\n" << endl;
}

/*Move the swarm to new solutions */
void Swarm::moveSwarm(Configuration* config, long int iteration, const double minBound, const double maxBound) {
	//Update gBest of each particle
	for (unsigned int i=0;i<swarm.size();i++){
		//The id of the gBest particle depends on the topology and model of influence
		swarm.at(i)->getBestOfNeibourhood();
	}

	//For the cases in which the entire swarm is using omega1 with the same value
	computeOmega1(config, iteration, -1, true); // -1 is a place holder for the id of the particle

	//We call this method here because some models of influence need to initialize things
	getInformants(config,-1,iteration);	// -1 is a place holder for the id of the particle

	//Compute the acceleration coefficients of the entire swarm
	computeAccelerationCoefficients(config, iteration);

	//Move particles
	cout << "iteration: " << iteration << endl; //remove

	for (unsigned int i=0;i<swarm.size();i++){

		cout << "\tParticle [" << i << "] -- gBestID [" << swarm.at(i)->getgBestID() << "] -- "; //<< endl; //remove
		int sizeInformants = getInformants(config, i, iteration); //Get the informants of i

		//print all neighbors
		cout << "\tNeighbors ids:  [ ";
		for (unsigned int j=0;j<swarm.at(i)->neighbours.size();j++){
			cout << swarm.at(i)->neighbours[j]->getID() << " ";
		}
		cout << "]" << endl;
		//print all neighbors
		cout << "\tInformants pos: [ ";
		for (int j=0;j<sizeInformants;j++){
			cout << Informants[j] << " ";
		}
		cout << "]" << endl;

		//If using a self-adaptive strategy compute omega1, otherwise this function simply returns the value already computed
		double omega1 = computeOmega1(config, iteration, i, false);

		//When FI or RankedFI phi_2 has to be decomposed according the number of informants
		decomposePhi2(config->getModelOfInfluence(), i, sizeInformants);

		//Note that here computeOmega1 receives the number of the particle and the flag = false
		swarm.at(i)->move(config, minBound, maxBound, iteration,
				omega1,
				computeOmega2(config),
				computeOmega3(config),
				sizeInformants,
				Informants,
				lastLevelComplete,
				alpha_t, l, delta);
	}

	long double prev_Gbest_eval = global_best.eval; //best solution at iteration t-1

	//Update best_particle
	for (unsigned int i=0;i<swarm.size();i++){
		if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
			//Update best_particle (a.k.a Gbest) of the Swarm (not to be confused with the gBest of a particle)
			Swarm::updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
			best_particle = swarm[i];
		}
	}
	updatePerturbationVariables(config, prev_Gbest_eval, global_best.eval, iteration);
}

int Swarm::getInformants(Configuration* config, int particleID, long int iteration){
	if (particleID != -1){
		//Best of neighborhood
		if (config->getModelOfInfluence() == MOI_BEST_OF_N){
			if (config->getTopology() == TOP_HIERARCHICAL ){
				int Array_size = 0;		//variable to the the size of Informants
				int * TMP_Array;
				TMP_Array = new int[lastLevelComplete];
				getParticleParentsIDs(particleID, TMP_Array); //Get parents of the particle
				for (int i=0; i<=lastLevelComplete; i++){
					if (TMP_Array[i] != -2){ //-2 indicates an empty position
						Array_size++;
					}
					else
						break;
				}
				//This is the actual array with the ID of the informants
				if (iteration == 1 && particleID == 0)
					Informants = new int[Array_size];
				else {
					delete [] Informants;
					Informants = new int[Array_size];
				}

				//We need to find the position index of the IDs in TMP_Array in the particles neighbours vector
				bool exitFlag;
				for (unsigned int i=0; i<swarm.at(particleID)->neighbours.size(); i++){
					for (int j=0; j<Array_size; j++)
						if (swarm.at(particleID)->neighbours.at(i)->getID() == TMP_Array[j]){
							Informants[j] = i;
							exitFlag = true;
							break;
						}
					if (exitFlag)
						break;
				}
				cout << "Size of Informants of " << particleID << " is " << 1 << endl;
				return 1;
			}
			else {
				if (iteration == 1 && particleID == 0) //allocate memory only first time
					Informants = new int [1];
				else {
					delete [] Informants;
					Informants = new int [1];
				}
				int bestID = swarm.at(particleID)->getgBestID();
				//We need to find the position index of the IDs in TMP_Array in the particles neighbours vector
				for (unsigned int i=0; i<swarm.at(particleID)->neighbours.size(); i++){
					if (swarm.at(particleID)->neighbours.at(i)->getID() == bestID)
						Informants[0] = i;
				}
				cout << "Size of Informants: " << 1 << endl;
				return 1;

			}
		}
		//Fully informed
		else if (config->getModelOfInfluence() == MOI_FI) {
			if (config->getTopology() == TOP_HIERARCHICAL ){
				int Array_size = 0;		//variable to the the size of Informants
				int * TMP_Array;
				TMP_Array = new int[lastLevelComplete];
				getParticleParentsIDs(particleID, TMP_Array); //Get parents of the particle
				for (int i=0; i<=lastLevelComplete; i++){
					if (TMP_Array[i] != -2){ //-2 indicates an empty position
						Array_size++;
					}
					else
						break;
				}
				//This is the actual array with the ID of the informants
				if (iteration == 1 && particleID == 0)
					Informants = new int[Array_size];
				else {
					delete [] Informants;
					Informants = new int[Array_size];
				}
				//We need to find the position index of the IDs in TMP_Array in the particles neighbours vector
				for (unsigned int i=0; i<swarm.at(particleID)->neighbours.size(); i++){
					for (int j=0; j<Array_size; j++)
						if (swarm.at(particleID)->neighbours.at(i)->getID() == TMP_Array[j])
							Informants[j] = i;
				}
				//delete [] TMP_Array;
				//swarm.at(particleID)->getBestOfNeibourhood(); //update particle's gbest
				cout << "Size of Informants of " << particleID << " is " << Array_size << endl;
				return Array_size;
			}
			else {
				//Since some topologies are dynamic, the size of informants may change from iteration to iteration
				if (iteration == 1 && particleID == 0)
					Informants = new int[swarm.at(particleID)->neighbours.size()];
				else {
					delete [] Informants;
					Informants = new int[swarm.at(particleID)->neighbours.size()];
				}
				for (unsigned int i=0;i<swarm.at(particleID)->neighbours.size();i++){
					Informants[i] = i; //we use the indexes of neighbours
				}
				cout << "Size of Informants is: " << swarm.at(particleID)->neighbours.size() << endl;
				return swarm.at(particleID)->neighbours.size();
			}
		}
		//Ranked fully informed
		else if (config->getModelOfInfluence() == MOI_RANKED_FI) {
			if (iteration == 1 && particleID == 0)
				Informants = new int[swarm.at(particleID)->neighbours.size()];
			else {
				delete [] Informants;
				Informants = new int[swarm.at(particleID)->neighbours.size()];
			}

			//Container to sort the neighbors
			vector< vector<int> > TMP_vect;
			//Resize vector
			TMP_vect.resize((swarm.at(particleID)->neighbours.size()), vector<int>(2));

			//This is the same as Fully informed
			for (unsigned int i=0;i<swarm.at(particleID)->neighbours.size();i++){
				TMP_vect.at(i).at(0) = i; ///we use the indexes of neighbours
				TMP_vect.at(i).at(1) = swarm.at(particleID)->neighbours.at(i)->getRanking();
			}

			//we sort by the second column (see sortcol driver function above )
			sort(TMP_vect.begin(), TMP_vect.end(), sortcol);

			//copy informants ID to Informants sorted
			for (unsigned int i=0; i<TMP_vect.size(); i++){ //rows
				Informants[i] = TMP_vect[i][0];
				TMP_vect[i].clear();
			}
			TMP_vect.clear();

			cout << "Size of Informants is: " << swarm.at(particleID)->neighbours.size() << endl;
			return swarm.at(particleID)->neighbours.size();
		}
		else {
			cerr << "No model of influence matches the available options" << endl;
			exit (-1);
		}
	}
	else{
		//Ranked fully informed
		if (config->getModelOfInfluence() == MOI_RANKED_FI){
			//Implement ranks if particles are not using them already
			if (config->getOmega1CS() != IW_RANKS_BASED){
				if (iteration == 1) {
					if (!modInfRanked) //flag to delete the structure at the end
						modInfRanked = true;
					rankedSwarm.eval = new long double [config->getSwarmSize()];
					rankedSwarm.id = new int [config->getSwarmSize()];
					//Rank particles
					rankParticles(&rankedSwarm);
					cout << "\nRanking swarm..." << endl;
				}
				else
					rankParticles(&rankedSwarm);
			}
		}
		return 0;
	}
}

void Swarm::decomposePhi2(int modelOfInflu, int part_id, int numInformants){
	if (modelOfInflu == MOI_FI)
		swarm.at(part_id)->setPhi2(swarm.at(part_id)->getPhi2()/numInformants);
}

void Swarm::computeAccelerationCoefficients(Configuration* config, long int iteration){
	//If the strategy involves all particles using the same acceleration coefficients value,
	//it is more efficient to compute it once at the beginning of the iteration
	//Constant value
	switch (config->getAccelCoeffCS()) {
	case AC_CONSTANT:{
		for (unsigned int i=0; i<swarm.size(); i++){
			swarm.at(i)->setPhi1(config->getPhi1());
			swarm.at(i)->setPhi2(config->getPhi2());
		}
	} break;
	//Random values within bounds
	case AC_RANDOM:{
		for (unsigned int i=0; i<swarm.size(); i++){
			swarm.at(i)->setPhi1(problem->getRandomX(config->getInitialPhi1(),config->getFinalPhi1()));
			swarm.at(i)->setPhi2(problem->getRandomX(config->getFinalPhi2(),config->getInitialPhi2()));
		}
	} break;
	case AC_EXTRAPOLATED:{
		for (unsigned int i=0; i<swarm.size(); i++){
			double varPhi1 = exp((double)(-iteration/config->getMaxIterations()));
			double distanceToGbest = (swarm.at(swarm.at(i)->getgBestID())->getCurrentEvaluation()-swarm.at(i)->getCurrentEvaluation())/
					swarm.at(swarm.at(i)->getgBestID())->getCurrentEvaluation();
			swarm.at(i)->setPhi1(varPhi1);
			swarm.at(i)->setPhi2(exp(varPhi1*distanceToGbest));
		}
	} break;
	case AC_TIME_VARYING:{
		double varPhi1 = (config->getInitialPhi1()-config->getFinalPhi1()) * (double)(-iteration/config->getMaxIterations()) + config->getInitialPhi1();
		double varPhi2 = (config->getInitialPhi2()-config->getFinalPhi2()) * (double)(-iteration/config->getMaxIterations()) + config->getInitialPhi2();
		for (unsigned int i=0; i<swarm.size(); i++){
			swarm.at(i)->setPhi1( varPhi1 );
			swarm.at(i)->setPhi2( varPhi2 );
		}
	} break;
	}
}


void Swarm::updatePerturbationVariables(Configuration* config, double previousGbest_eval, double currentGbest_eval, long int iteration){
	if(		config->getPerturbation1Type() == PERT1_NORMAL_SUCCESS ||
			config->getPerturbation2Type() == PERT2_ADD_RECT ||
			config->getPerturbation1Type() == PERT1_CAUCHY_SUCCESS ){
		if (previousGbest_eval != currentGbest_eval){ //success
			sc++;
			fc=0;
			if (sc > success){
				alpha_t = alpha_t * 2.0;
				//cout << "\n DOUBLE alpha_t for success" << alpha_t << endl;
			}
		}
		else{ //failure
			fc++;
			sc=0;
			if (fc > failure){
				alpha_t = alpha_t * 0.5;
				//cout << "\n HALF alpha_t update for failure: " << alpha_t << endl;
			}
		}
		//if alpha_t becomes too small, we reinitialize it to 0.15, if we are in the first half of the
		//iterations, or to 0.001 if we are in the second half of the iterations.
		if (alpha_t < ALPHA_T_PRECISION){
			if (iteration < config->getMaxIterations()/2) //first half
				alpha_t = 0.15;
			else
				alpha_t = 0.001;
		}
	}
}

/* Topologies */
void Swarm::createFullyConnectedTopology(){       //All particles are neighbor among them
	for(unsigned int i=0;i<swarm.size();i++){
		for(unsigned int j=0;j<swarm.size();j++){
			swarm.at(i)->addNeighbour(swarm.at(j));
		}
	}
}

void Swarm::createRingTopology(){       //Every particle is neighbor of the adjacent particles
	int a,b;
	for(unsigned int i=0;i<swarm.size();i++){
		a=i-1;
		b=i+1;
		if(i==0)
			a=swarm.size()-1;
		if(i==(swarm.size()-1))
			b=0;

		swarm.at(i)->addNeighbour(swarm.at(a));
		swarm.at(i)->addNeighbour(swarm.at(b));
	}
}

void Swarm::createWheelTopology(){        //Particles are neighbors of one central particle
	for(unsigned int i=1;i<swarm.size();i++){
		swarm.at(i)->addNeighbour(swarm.at(0));
		swarm.at(0)->addNeighbour(swarm.at(i));
	}
}

void Swarm::createRandomEdge(){			//Or random edge topology
	long int randomEdge;
	for(unsigned int i=0;i<swarm.size();i++){
		randomEdge = (int)floor(RNG::randVal(0.0,(double)swarm.size()));
		if (randomEdge == i){
			randomEdge = (int)floor(RNG::randVal(0.0,(double)swarm.size()));
			swarm.at(i)->addNeighbour(swarm.at(randomEdge));
		}
		else
			swarm.at(i)->addNeighbour(swarm.at(randomEdge));
	}
}

void Swarm::createVonNeumannTopology(){
	int a,b,c;
	for(unsigned int i=0;i<swarm.size();i++){
		a=i-1;
		b=a-1;
		c=i+1;
		if(i==0){
			a=swarm.size()-1;
			b=a-1;
		}
		if (i==1){
			a=swarm.size()-1;
			b=0;
		}
		if(i==(swarm.size()-1))
			c=0;

		swarm.at(i)->addNeighbour(swarm.at(a));
		swarm.at(i)->addNeighbour(swarm.at(b));
		swarm.at(i)->addNeighbour(swarm.at(c));
	}
}

void Swarm::updateTimeVaryingTopology(Configuration* config, long int iterations){
	//Topology update
	if( (iterations > 0) && (config->getEsteps() < swarm.size()-3) && (iterations%config->getTopologyUpdatePeriod() == 0)){
		unsigned int removals = 0;
		//cout << " -- esteps " << config->getEsteps() << endl;
		//cout << " -- Update topology at iteration: " << iterations << " Target: " << swarm.size()-(2+config->getEsteps()) << endl;
		RNG::shufflePermutation();
		while(removals < swarm.size()-(2+config->getEsteps())){
			for(unsigned int i=0;i<swarm.size();i++){
				int particleIndex = RNG::getPermutationElement(i);
				if( swarm.at(particleIndex)->getNeighborhoodSize() > 3 ){ //3 because a particle is a neighbor to itself
					int neighborID = swarm.at(particleIndex)->getRandomNonAdjacentNeighborID(config);

					//cout << " -- Erasing edge " << particleIndex << " <---> " << neighborID << endl;

					swarm.at(particleIndex)->eraseNeighborbyID(neighborID);
					swarm.at(neighborID)->eraseNeighborbyID(particleIndex);

					removals++;
				}
				if(removals == swarm.size()-(2+config->getEsteps()))
					break;
			}
		}
		config->setEsteps(config->getEsteps()+1);
		//cout << "Removals " << removals << endl;
		//for(unsigned int i=0;i<swarm.size();i++)
		//	cout << "particleIndex: " << i << " -- Neighbors: " << swarm.at(i)->getNeighborhoodSize() << endl;
	}
}

void Swarm::createHierarchical(int branching){
	//cout  << "\nCreating hierarchical topology..." << endl;
	//The topology is fully-connected, but we use a hierarchy as a model of influence
	long int firstPart = 1;		//id of the last particle that can be added in a complete level without exceeding size
	int actualSpace=1;			//actual number of nodes than can be stored in a given level. This depends on the swarm size and the branching degree

	//Find last level that can be filled in the tree with max width and the id of the last particle
	while (size - pow(branching,lastLevelComplete) >= pow(branching,lastLevelComplete+1)) {
		lastLevelComplete++;
		firstPart+=pow(branching,lastLevelComplete);
	}
	//It can be the case that for a very large branching degree the firstPart variable
	//is larger than the size of the swarm. In those cases, we find the last level complete
	//without exceeding the swarm size
	while (firstPart > size){
		lastLevelComplete--;
		firstPart-=pow(branching,lastLevelComplete);
	}
	//Compute actual space
	for (int i=1; i<=lastLevelComplete; ++i)
		actualSpace+=pow(branching,i);
	if (firstPart > actualSpace)
		firstPart = actualSpace;
	int maxWidth = pow(branching,lastLevelComplete+1);

	//cout << "Last level complete: \t" << lastLevelComplete << endl;
	//cout << "Actual space: \t\t" << actualSpace << endl;
	//cout << "First part: \t\t" << firstPart << endl;
	//cout << "Max width: \t\t" << maxWidth << endl;
	//cout << "Branching degree: \t" << branching << endl << endl;

	//Create a fully connected topology first
	for(unsigned int i=0;i<swarm.size()	;i++){
		for(unsigned int j=0;j<swarm.size();j++){
			swarm.at(i)->addNeighbour(swarm.at(j));
		}
	}

	//Initialize tree structure
	hierarchy = new int*[lastLevelComplete+1];	//these are the level (depth of the tree)
	for (int i=0; i<=lastLevelComplete+1; i++)
		hierarchy[i] = new int[maxWidth];		//these are the nodes in a level (max width of the tree)

	for (int i=0; i<=lastLevelComplete+1; i++){
		for (int j=0; j<maxWidth; j++)
			hierarchy[i][j] = -2;				//-2 indicates that the position is not in use
	}

	//tree variables
	int d = 0;				//node's degree
	int nodesCounter = 0; 	//every time we reach the max branching degree, we move to the next node
	int h = 1;				//height of the tree
	int width = 0;			//width of the level

	//root node
	int parent = hierarchy[0][0] = swarm.at(0)->getID();
	swarm.at(0)->setParent(-1); //the parent of root is NULL
	//We first fill all the levels that will be complete, i.e., the levels with maximum width
	for(unsigned int i=1;i<firstPart;i++){
		if (width < pow(branching,h)) {
			hierarchy[h][width] = swarm.at(i)->getID();
			swarm.at(i)->setParent(parent);
			width++;	//tree's degree
			d++;		//node's degree
			if (d==branching) {
				nodesCounter++;
				parent = hierarchy[h-1][nodesCounter];
				d=0;
			}
			if (width == pow (branching,h)){
				h++;
				width=0;
				nodesCounter=0;
				parent = hierarchy[h-1][nodesCounter];
			}
		}
	}
	//reset variables of the tree for the rest of the particles
	width = 0;
	nodesCounter=0;
	int iterCount = 0;
	parent = hierarchy[lastLevelComplete][nodesCounter];
	//Now we add one d at a time to each node in the last level complete
	for(unsigned int i=firstPart;i<size;i++){
		hierarchy[lastLevelComplete+1][width+(nodesCounter*(branching-1))+iterCount] = swarm.at(i)->getID();
		swarm.at(i)->setParent(parent);
		width++;		//last level degree
		nodesCounter++; //we add one d to each node to have a fair distribution
		parent = hierarchy[lastLevelComplete][nodesCounter];
		if (width == pow (branching,lastLevelComplete)){
			width=0;
			nodesCounter=0;
			iterCount++;
			parent = hierarchy[lastLevelComplete][nodesCounter];
		}
	}
	//updateTree(branching);
	//printTree(branching);
	//printAllParentNodes();
}

void Swarm::getParticleParentsIDs(int particleID, int *ParentsArray1D){
	//Initialize the array to -2
	for (int i=0; i<=lastLevelComplete; i++)
		ParentsArray1D[i] = -2;

	//Find particle index in the swarm by its ID
	long int index;
	for(unsigned int i=0; i<swarm.size(); i++){
		if (particleID == swarm.at(i)->getID()){
			index = i;
			break;
		}
	}
	//cout << "\nParticle " << swarm.at(index)->getID() << ": " ;

	//Get direct parentNode ID
	int aParent = swarm.at(index)->getParent();
	if (aParent==-1){
		//cout << "Root node";
		//ParentsArray1D[0] = swarm.at(index)->getRandomNeighbor();
		ParentsArray1D[0] = index;
		swarm.at(particleID)->setgBestID(index);	//we force the use of reinitialization to model
	}
	else{
		int pos = 0;
		while (aParent != -1){
			ParentsArray1D[pos] = aParent;
			//cout << ParentsArray1D[pos] << " ";
			pos++;
			aParent = swarm.at(aParent)->getParent();
		}
	}
	//cout << endl ;
	//return ParentsArray1D;
}

void Swarm::updateTree(int branching){
	//Traverse the tree and update the position of the particles
	//reset indexes to traverse the entire tree
	//cout << endl << "Traversing the tree: " << endl;
	int h = 1;
	int width = 0;
	int iterCount = 0;
	int childCounter = 1;
	int newParent = 0;
	int newH = 0;
	int newWidth = 0;
	int parentNode = swarm.at(hierarchy[0][0])->getParent();	//start at the root top-bottom
	int parentH = 0;
	int parentWidth = -1;
	for(unsigned int i=1;i<size;i++){
		if (h <= lastLevelComplete){
			if (width < pow(branching,h)){ //max level width
				if (parentNode!=swarm.at(hierarchy[h][width])->getParent()){
					parentNode=swarm.at(hierarchy[h][width])->getParent();
					//cout << "Parent: " << parentNode << endl;	//Parent coordinates
					parentH = h-1;
					parentWidth++;
				}
				//cout << "\t\t" << hierarchy[h][width] << "->" << swarm.at(hierarchy[h][width])->getParent() << endl;
				if (swarm.at(hierarchy[h][width])->getCurrentEvaluation() < swarm.at(parentNode)->getCurrentEvaluation()){
					newParent = swarm.at(hierarchy[h][width])->getID();	//newParent coordinates (this is the child that becomes parent)
					newH = h;
					newWidth = width;
					//we change places only when newParent is different from getParent() of the last particle
					//this means that at least one of the children will become parent.
					//The id of this particle is saved in newParent
					if (width > 0 && (width+1) % branching == 0 && (newParent != swarm.at(hierarchy[h][width])->getParent())){
						swapNodes(newParent, newH, newWidth, parentNode, parentH, parentWidth, branching, h, width, iterCount);
					}
				}
				width++;
				childCounter++;
			}
			if (width == pow(branching,h)){
				h++;
				width=0;
				parentWidth=-1;
				parentNode=-1;
				//cout << endl << "Level " << h << ": " << endl;
			}
		}
		else {
			parentNode = hierarchy[h-1][iterCount];
			parentH = h-1;
			parentWidth++;
			newParent = parentNode;
			//cout << "Parent: " << parentNode << endl;

			for(int j=iterCount*branching;j<(iterCount*branching)+branching;j++){
				try{
					if (swarm.at(hierarchy[h][j])->getParent() == parentNode) {
						//cout << "\t\t" << hierarchy[h][j] << "->" << swarm.at(hierarchy[h][j])->getParent() << endl;
						childCounter++;
					}
					//Check if there is a new parent
					if (swarm.at(hierarchy[h][j])->getCurrentEvaluation() < swarm.at(parentNode)->getCurrentEvaluation()){
						newParent = swarm.at(hierarchy[h][j])->getID();	//newParent coordinates (this is the  that becomes parent)
						newH = h;
						newWidth = j;
					}
					//Switch parent and child after comparing with all children
					if (j == (iterCount*branching)+branching-1 && (newParent != swarm.at(hierarchy[h][j])->getParent())){
						swapNodes(newParent, newH, newWidth, parentNode, parentH, parentWidth, branching, h, width, iterCount);
					}
				}
				catch (const std::exception& e) {
					break;
				}
			}
			iterCount++;
		}
		if (childCounter==size)
			break;
	}
}

void Swarm::swapNodes(int newParent, int newH, int newWidth, int parentNode, int parentH, int parentWidth, int branching, int h, int width, int iterCount){
	//cout << "\t\t New parent: " << newParent << " coordinates: "  << newH << " "<< newWidth << endl;
	//cout << "\t\t Replaces: " << parentNode << " coordinates: "  << parentH << " "<< parentWidth << endl;

	//Swap child and parent IDs
	swarm.at(hierarchy[newH][newWidth])->setParent(swarm.at(hierarchy[parentH][parentWidth])->getParent());
	swarm.at(hierarchy[parentH][parentWidth])->setParent(newParent); //parent

	//Swap places
	hierarchy[parentH][parentWidth] = newParent;
	hierarchy[newH][newWidth] = parentNode;

	//cout << "\t\t s" << hierarchy[newH][newWidth] << "->" << swarm.at(hierarchy[newH][newWidth])->getParent() << endl;
	//cout << "\t\t s" << hierarchy[parentH][parentWidth] << "->" << swarm.at(hierarchy[parentH][parentWidth])->getParent() << endl;

	if (h <= lastLevelComplete) {
		//Update brothers parentID
		for (int k=width-(branching-1); k<=width; k++){
			swarm.at(hierarchy[h][k])->setParent(newParent);
			//cout << "\t\t b" << hierarchy[h][k] << "->" << swarm.at(hierarchy[h][k])->getParent() << endl;

		}
		//Update any children that the old child (now parent) may have had
		for (int k=branching*width; k<(branching*width)+branching; k++){
			try {
				swarm.at(hierarchy[h+1][k])->setParent(parentNode);
				//cout << "\t\t c" << hierarchy[h+1][k] << "->" << swarm.at(hierarchy[h+1][k])->getParent() << endl;
			}
			catch (const std::exception& e) {
				break;
			}
		}
	}
	if (h == lastLevelComplete+1) {
		//Update brothers parentID
		for(int k=iterCount*branching;k<(iterCount*branching)+branching;k++){
			try{
				swarm.at(hierarchy[h][k])->setParent(newParent);
				//cout << "\t\t b" << hierarchy[h][k] << "->" << swarm.at(hierarchy[h][k])->getParent() << endl;
			}
			catch (const std::exception& e) {
				break;
			}
		}
	}
}

void Swarm::printAllParentNodes(){
	int parentNode = 0;
	for(unsigned int i=0; i<swarm.size(); i++){
		cout << "Particle " << i << ": " ;
		parentNode = swarm.at(i)->getParent();

		if (parentNode==-1){
			cout << "Root node";
		}
		else{
			while (parentNode != -1){
				cout << parentNode << " ";
				parentNode = swarm.at(parentNode)->getParent();
			}
		}
		cout << endl ;
	}
}

void Swarm::printTree(int branching) {
	//Raw printing of the structure
	cout << endl << endl;
	for (int i=0; i<=lastLevelComplete+1; i++){
		cout << "Level " << i << ": ";
		for (int j=0; j<pow(branching,lastLevelComplete+1); j++)
			cout << hierarchy[i][j] << " ";
		cout << endl;
	}

	//Print tree -- this is the same algorithm we use to traverse the tree and update the position of the particles
	cout << "\n\nRoot node: " << endl;
	cout << "\t\t" << hierarchy[0][0] << "->" << swarm.at(hierarchy[0][0])->getParent() << endl << "Level 1: " << endl;
	//reset indexes to traverse the entire tree
	int h = 1;
	int width = 0;
	int iterCount = 0;
	int childCounter = 1;
	int parentNode = swarm.at(hierarchy[0][0])->getParent();
	for(unsigned int i=1;i<size;i++){
		if (h <= lastLevelComplete){
			if (width < pow(branching,h)){ //max level width
				if (parentNode!=swarm.at(hierarchy[h][width])->getParent()){
					parentNode=swarm.at(hierarchy[h][width])->getParent();
					cout << "Parent: " << parentNode << endl;
				}
				cout << "\t\t" << hierarchy[h][width] << "->" << swarm.at(hierarchy[h][width])->getParent() << endl;
				width++;
				childCounter++;
			}
			if (width == pow(branching,h)){
				h++;
				width=0;
				cout << "Level " << h << ": " << endl;
			}
		}
		else {
			parentNode = hierarchy[h-1][iterCount];
			cout << "Parent: " << parentNode << endl;
			for(int j=iterCount*branching;j<(iterCount*branching)+branching;j++){
				try{
					if (swarm.at(hierarchy[h][j])->getParent() == parentNode) {
						cout << "\t\t" << hierarchy[h][j] << "->" << swarm.at(hierarchy[h][j])->getParent() << endl;
						childCounter++;
					}
				}
				//This exception is thrown when hierarchy[][] = -2 and we try to access a position of swarm.at(-2), which does not exist
				catch (const std::exception& e) {
					//cout << "h[" << h << "] j[" << j << "] iter " << iterCount << endl;
					break;
				}
			}
			iterCount++;
		}
		if (childCounter==size)
			break;
	}
}

Solution Swarm::getGlobalBest(){
	return global_best;
}

bool Swarm::isHierarchical(){
	return hierarchical;
}

double Swarm::computeAvgVelocity(Configuration* config){
	double sumVelocity = 0.0;
	for (unsigned int i=0;i<swarm.size();i++){
		for (unsigned int j=0;j<config->getProblemDimension();j++){
			sumVelocity = sumVelocity + abs(swarm.at(i)->getCurrentVelocity()[j]);
		}
	}
	return (sumVelocity/(swarm.size()*config->getProblemDimension()));
}

void Swarm::rankParticles(SimplifySwarm* sS){
	//The ranking are obtained via the mergesort algorithm implemented in Utils.cpp

	//Copy particle's id and evaluation in simpSwarm
	for (unsigned int i=0;i<swarm.size();i++){
		sS->id[i] = swarm.at(i)->getID();
		sS->eval[i] = swarm.at(i)->getCurrentEvaluation();
	}
	//Sort simpSwarm by the value of particles' evaluation
	Utils::mergeSort(sS, 0, swarm.size()-1); //mergeSort(array, left(LOWER) index, right (UPPER) index);

	//Set particles' rank in the swarm
	for (unsigned int i=0;i<swarm.size();i++){
		swarm.at(sS->id[i])->setRanking(i+1);
	}
}

double Swarm::computeOmega2(Configuration* config){
	switch(config->getOmega2CS()){
	case O2_EQUALS_IW:
		return config->getOmega1();
	case O2_RANDOM:
		return 0.5 * (problem->getRandom01()/2.0);
	case O2_ZERO:
		return 0.0; //the component is not used
	default :
		return 1.0; //no strategy, set the value to one
	}
}

double Swarm::computeOmega3(Configuration* config){
	//Same as Omega1
	if (config->getOmega3CS() == O3_EQUALS_IW)
		return config->getOmega1();
	//Random value
	else if (config->getOmega3CS() == O3_RANDOM)
		return 0.5 * (problem->getRandom01()/2.0);
	//Zero -- component is not being used
	else if (config->getOmega3CS() == O3_ZERO)
		return 0.0;
	//One -- no strategy in particular
	else
		return 1.0;
}

// This function computes the inertia weight (omega1 in the GVU formula) according to the selected strategy
double Swarm::computeOmega1(Configuration* config, long int iteration, long int id, bool newIteration){
	//If all particle use the same inertia value, it is more efficient
	//to compute it once at the beginning of the iteration
	if (newIteration) {
		/* Non-adaptive strategies */
		//IW_CONSTANT - 0 - Constant value
		if (config->getOmega1CS() == IW_CONSTANT){
			if (config->getOmega1() < -1 || config->getOmega1() > 1) //check convergence bounds
				config->setOmega1(CONSTRICTION_COEFFICIENT);
		}
		//IW_L_INC - 1 - Linear increasing
		if (config->getOmega1CS() == IW_L_INC) {
			//from Frankenstein's PSO
			if(iteration <= config->getIWSchedule()){
				config->setOmega1(
						((double)(config->getIWSchedule() - iteration)/config->getIWSchedule())*
						(config->getFinalIW() - config->getInitialIW()) + config->getInitialIW()
				);
				//cout << config->getInertia() << endl;
			}
			else
				config->setOmega1(config->getFinalIW());
		}
		//IW_L_DEC - 2 - Linear decreasing
		else if (config->getOmega1CS() == IW_L_DEC) {
			//from Frankenstein's PSO
			if(iteration <= config->getIWSchedule()){
				config->setOmega1(
						((double)(config->getIWSchedule() - iteration)/config->getIWSchedule())*
						(config->getInitialIW() - config->getFinalIW()) + config->getFinalIW()
				);
				//cout << config->getInertia() << endl;
			}
			else
				config->setOmega1(config->getFinalIW());
		}
		//IW_RANDOM - 3 - Random
		else if (config->getOmega1CS() == IW_RANDOM) {
			config->setOmega1( 0.5 * (problem->getRandom01()/2.0));
			//cout << config->getInertia() << endl;
		}
		//IW_NONL_DEC - 4 - Nonlinear decreasing
		else if (config->getOmega1CS() == IW_NONL_DEC) {
			config->setOmega1(
					config->getFinalIW() - (config->getFinalIW()-config->getInitialIW())*
					pow((double)(iteration)/config->getMaxIterations(),alpha)
			);
			//cout << iteration << endl;
			//cout << config->getInertia() << endl;
		}
		//IW_NONL_DEC_IMP - 5 - Nonlinear decreasing improved
		else if (config->getOmega1CS() == IW_NONL_DEC_IMP) {
			config->setOmega1(
					pow(omegaXu, iteration)
			);
			//cout << iteration << endl;
			//cout << config->getInertia() << endl;
		}
		//IW_NONL_DEC_TIME - 6 - Nonlinear decreasing time-dependent
		else if (config->getOmega1CS() == IW_NONL_DEC_TIME) {
			config->setOmega1(
					pow((2.0/iteration), 0.3)
			);
			//cout << iteration << endl;
			//cout << config->getInertia() << endl;
		}
		//IW_CHAOTIC_DEC - 7 Chaotic decreasing
		else if (config->getOmega1CS() == IW_CHAOTIC_DEC) {
			iteration == 1 ? zFunction = problem->getRandom01() : zFunction = 4*zFunction*(1-zFunction);
			config->setOmega1(
					(zFunction*config->getInitialIW()) + (config->getFinalIW()-config->getInitialIW()) *
					((double)(config->getMaxIterations())-iteration)/config->getMaxIterations()
			);
			//cout << iteration << endl;
			//cout << config->getInertia() << endl;
		}
		//IW_EXP_DEC - 8 - Natural exponential decreasing
		else if (config->getOmega1CS() == IW_EXP_DEC) {
			config->setOmega1(
					config->getInitialIW() + (config->getFinalIW()-config->getInitialIW())*
					exp((-10.0*iteration)/config->getMaxIterations())
			);
			//cout << iteration << endl;
			//cout << config->getInertia() << endl;
		}
		//IW_OSCILLATING - 9 - Oscillating
		else if (config->getOmega1CS() == IW_OSCILLATING) {
			if (iteration < (3*config->getMaxIterations())/4)
				config->setOmega1(
						((config->getInitialIW() + config->getFinalIW()) /2.0) +
						((config->getFinalIW() + config->getInitialIW()) /2.0) *
						cos((simNumOfCos*iteration)/(config->getMaxIterations()*3.0))
				);
			else
				config->setOmega1(config->getInitialIW());
			//cout << iteration << endl;
			//cout << config->getInertia() << endl;
		}
		//IW_LOG_DEC - 10 - Logarithm decreasing
		else if (config->getOmega1CS() == IW_LOG_DEC) {
			config->setOmega1(
					config->getFinalIW() + (config->getInitialIW()-config->getFinalIW())*
					log10(((10.0*iteration)/config->getMaxIterations())+ a )
			);
			//cout << iteration << endl;
			//cout << config->getInertia() << endl;
		}

		/* ****************************************************************************************************************/
		/* ****************************************************************************************************************/

		/* Adaptive strategies */
		//IW_SELF_REGULATING - 11 - Self-regulating
		else if (config->getOmega1CS() == IW_SELF_REGULATING) {
			static double deltaOmega = (((double)config->getFinalIW()) - config->getInitialIW())/config->getMaxIterations();
			iteration == 1 ? omega_2 = config->getFinalIW() : omega_2 = omega_2-deltaOmega;
			config->setOmega1(omega_2);
			//cout << iteration << endl;
			//cout << config->getInertia() << endl;
		}
		//IW_VELOCITY_BASED - 12 - Based on velocity information
		else if (config->getOmega1CS() == IW_VELOCITY_BASED) {

			static double T_0_95 = 95*config->getMaxIterations()/100; //iteration at which 95% of search process is completed

			if (iteration == 1)
				config->setOmega1(config->getFinalIW());

			idealVelocity = swarm.at(0)->getMaxVelLimit() * ((1.0 + cos(M_PI*(iteration/T_0_95)))/2);
			avVel=computeAvgVelocity(config);	//average absolute velocity of the swarm

			if (avVel >= idealVelocity){
				(config->getOmega1()-deltaOmega) >= config->getInitialIW() ?
						config->setOmega1(config->getOmega1()-deltaOmega) : config->setOmega1(config->getInitialIW());
			}
			else{
				(config->getOmega1()+deltaOmega) >= config->getFinalIW() ?
						config->setOmega1(config->getFinalIW()) : config->setOmega1(config->getOmega1()+deltaOmega);
			}
			//cout << iteration << endl;
			//cout << avVel << " -- " << idealVelocity << endl;
			//cout << config->getInertia() << endl;
		}
		//IW_RANKS_BASED - 14 - Rank-based
		else if (config->getOmega1CS() == IW_RANKS_BASED) {
			if (iteration == 1) {
				//Memory allocation to rank particles
				simpSwarm.eval = new long double [config->getSwarmSize()];
				simpSwarm.id = new int [config->getSwarmSize()];
			}
			//Rank particles
			rankParticles(&simpSwarm);
			//The computation of the inertia weight is implemented in Particle::move
		}
		//IW_SUCCESS_BASED - 15 Success-based
		else if (config->getOmega1CS() == IW_SUCCESS_BASED) {
			if (iteration == 1) {
				//simpSwarm contains a simplified copy of the swarm at t-1
				simpSwarm.eval = new long double [config->getSwarmSize()];
				simpSwarm.id = new int [config->getSwarmSize()];
				//Copy particle's id and evaluation in simpSwarm
				for (unsigned int i=0;i<swarm.size();i++){
					simpSwarm.id[i] = swarm.at(i)->getID();
					simpSwarm.eval[i] = swarm.at(i)->getCurrentEvaluation();
				}
				config->setOmega1(config->getFinalIW()); //set inertia weight to its maximum value for the first iteration
			}
			else{
				int S_i = 0; //Number of solutions that improved after the last iteration
				for (unsigned int i=0; i<swarm.size(); i++){
					for (unsigned int j=0; j<sizeof(simpSwarm.id); j++){
						if (swarm.at(i)->getID() == simpSwarm.id[j]){
							//evaluate if the solution improved
							if (swarm.at(i)->getCurrentEvaluation() < simpSwarm.eval[j]){
								S_i++;
								break;
							}
							else
								break;
						}
					}
				}
				//set the inertia weight
				config->setOmega1( config->getInitialIW() + ((config->getFinalIW()-config->getInitialIW()) *
						((double)S_i/swarm.size()))
				);
				//Copy particle's id and evaluation in simpSwarm for the next iteration
				for (unsigned int i=0;i<swarm.size();i++){
					simpSwarm.id[i] = swarm.at(i)->getID();
					simpSwarm.eval[i] = swarm.at(i)->getCurrentEvaluation();
				}
				//cout << iteration << endl;
				//cout << "P[" << id << "]" << "rank: "<< ranking << " eval:" << current.eval << endl;
				//cout << config->getInertia() << endl;
			}

		}
		//IW_CONVERGE_BASED - 16 Convergence-based
		else if (config->getOmega1CS() == IW_CONVERGE_BASED) {
			if (iteration == 1) {
				//simpSwarm contains a simplified copy of the swarm at t-1
				simpSwarm.eval = new long double [config->getSwarmSize()];
				simpSwarm.id = new int [config->getSwarmSize()];
				//Copy particle's id and evaluation in simpSwarm
				for (unsigned int i=0;i<swarm.size();i++){
					simpSwarm.id[i] = swarm.at(i)->getID();
					simpSwarm.eval[i] = swarm.at(i)->getCurrentEvaluation();
				}
				config->setOmega1(1 - abs(alpha_2/(1+beta_2)));
			}
			else{
				//Copy particle's id and evaluation in simpSwarm for the next iteration
				for (unsigned int i=0;i<swarm.size();i++){
					simpSwarm.id[i] = swarm.at(i)->getID();
					simpSwarm.eval[i] = swarm.at(i)->getCurrentEvaluation();
				}
			}
		}
		//No strategy, inertia is constant during the execution
		else {
			config->setOmega1(config->getOmega1()); //kind of unnecessary, but ensures integrity
			//cout << inertia << endl;
			//return inertia;
		}
	}
	else {
		//These are the strategies that need to compute a independent inertia value for each particle
		//Note that the variables used here are computed/allocated when (newIteration == true).
		if (id != -1){
			//IW_SELF_REGULATING - 11 - Self-regulating
			if (config->getOmega1CS() == IW_SELF_REGULATING) {
				//The best particle has a special inertia value
				if (id == best_particle->getID()) {
					//compute the inertia weight for the global_best particle
					//cout << id << endl;
					//cout << omega_2 + eta * ((config->getFinalIW() - config->getInitialIW()) / config->getMaxIterations()) << endl;
					return omega_2 + eta * ((config->getFinalIW() - config->getInitialIW()) / config->getMaxIterations());
				}
				else {
					//cout << id << endl;
					//cout << config->getInertia() << endl;
					return config->getOmega1();
				}
			}
			//IW_DOUBLE_EXP - 13 - Double exponential self-adaptive
			else if (config->getOmega1CS() == IW_DOUBLE_EXP) {
				double R =0.0;
				swarm.at(id)->getBestOfNeibourhood();  //update particle's gbest
				if (iteration == 1){
					config->setOmega1(config->getFinalIW());
					//cout << id << endl;
					//cout << config->getInertia() << endl;
				}
				else {
					R = swarm.at(id)->computeDistPbestGbest()*((((double)config->getMaxIterations())-iteration)/config->getMaxIterations());
					config->setOmega1( exp(-1*exp((R*-1))));
					//cout << id << endl;
					//cout << config->getInertia() << endl;
				}
				return config->getOmega1();
			}
			//IW_RANKS_BASED - 14 - Rank-based
			else if (config->getOmega1CS() == IW_RANKS_BASED){
				config->setOmega1( config->getInitialIW() + ((config->getFinalIW()-config->getInitialIW()) *
						((double) swarm.at(id)->getRanking()/config->getSwarmSize()))
				);
				//cout << id << endl;
				//cout << config->getInertia() << endl;
				return config->getOmega1();
			}
			//IW_CONVERGE_BASED - 16 Convergence-based
			else if (config->getOmega1CS() == IW_CONVERGE_BASED) {
				if (iteration > 1) {
					//Find the particle in simpSwarm
					for (unsigned int j=0; j<sizeof(simpSwarm.id); j++){
						if (simpSwarm.id[j] == swarm.at(id)->getID()){
							//convergence factor
							long double c_i = abs(simpSwarm.eval[j] - swarm.at(id)->getCurrentEvaluation())/
									(simpSwarm.eval[j] + swarm.at(id)->getCurrentEvaluation());
							//diffusion factor
							long double d_i = abs(swarm.at(id)->getCurrentEvaluation() - global_best.eval) /
									(swarm.at(id)->getCurrentEvaluation() + global_best.eval);
							//set the inertia weight
							config->setOmega1(1 - abs(alpha_2*(1-c_i)) / (1+d_i)*(1+beta_2));
							break;
						}
					}
					//cout << id << endl;
					//cout << config->getInertia() << endl;
					return config->getOmega1();
				}
				else {
					return config->getOmega1();
				}
			}
			else {
				return config->getOmega1();
			}
		}
		else{
			cerr << "Unexpected error while computing the value of inertia" << endl;
			exit (-1);
		}
	}
	return config->getOmega1();
}
