/*
 * swarm.cpp
 *
 *  Created on: May 31, 2019
 *      Author: Christian L. Camacho Villalón
 */

#include "iostream"
#include <sstream>
#include "float.h"
#include "utils.h"
#include <bits/stdc++.h>
#include <new>
#include "problem.h"
#include "config.h"
#include "swarm.h"
#include "rng.h"

using namespace std;

bool sortcol(const vector<int>& v1, const vector<int>& v2) {
	return (v1[1] < v2[1]);
}

//Default constructor
Swarm::Swarm(){
	problem = NULL;
	size = 0;
	best_particle = NULL;
	init = false;
	lastLevelComplete = 0;
}

Swarm::~Swarm(){
	if (init) {
		// Memory allocated dynamically
		delete [] global_best.x;
		for (unsigned long int i=0;i<swarm.size();i++)
			delete swarm.at(i);
	}
	init=false;
}

Swarm::Swarm (Problem* problem, Configuration* config){
	this->problem = problem;
	size = config->getSwarmSize();

	/*Initialize global best*/
	global_best.x = new double[config->getProblemDimension()];
	for(int i=0;i<config->getProblemDimension();i++)
		global_best.x[i] = 0;
	global_best.eval=LDBL_MAX;

	for (long int i=0; i<size; i++) {
		Particle* aParticle = new Particle(problem, config, i, 0);
		swarm.push_back(aParticle);

		if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
			updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
			best_particle = swarm.at(i);
		}
	}
	if (best_particle == NULL){
		cerr << "The program could not determine an initial solution quality. " << endl;
		cerr << "HINT : Make sure that all input files can be accessed by the program" << endl;
		cerr << "     : This error is typically produced when the o.f. cannot be correctly evaluated" << endl;
		exit (-1);
	}

	if (config->getTopology() == TOP_FULLYCONNECTED) {
		createFullyConnectedTopology();
	} else if (config->getTopology() == TOP_HIERARCHICAL) {
		(config->getPopulationCS() == POP_CONSTANT) ? createHierarchical(config->getBranchingDegree(), swarm.size()) :
				createHierarchical(config->getBranchingDegree(), config->getFinalPopSize());
	} else if (config->getTopology() == TOP_RING) {
		createRingTopology();
	} else if (config->getTopology() == TOP_WHEEL) {
		createWheelTopology();
	} else if (config->getTopology() == TOP_RANDOM) {
		createRandomEdge();
	} else if (config->getTopology() == TOP_TIMEVARYING) {
		createFullyConnectedTopology();
		if (config->getPopulationCS() != POP_INCREMENTAL)
			config->setTopologyUpdatePeriod((int)floor((double)config->getTopologySchedule()/(swarm.size()-3)));
		else
			config->setTopologyUpdatePeriod((int)floor((double)config->getTopologySchedule()/(config->getFinalPopSize()-3)));
	} else if (config->getTopology() == TOP_VONNEUMANN) {
		createVonNeumannTopology();
	}
	else {
		cerr << "Wrong topology" << endl;
		exit (-1);
	}

	init=true;
}

int Swarm::countImprovedSolutions(Configuration* config, long int iteration){
	//Create simpSwarm if it does not exist and return 0
	if (iteration == 1) { //define simpSwarm if it does not exist
		simpSwarm.clear();
		simpSwarm.resize(swarm.size());
		for (unsigned int i=0;i<swarm.size();i++){
			simpSwarm.at(i).id = swarm.at(i)->getID();
			simpSwarm.at(i).eval = swarm.at(i)->getCurrentEvaluation();
		}
		return (0);
	}
	else {
		int S_i = 0; //Number of solutions that improved after the last iteration
		for (unsigned int i=0; i<swarm.size(); i++){
			for (unsigned int j=0; j<simpSwarm.size(); j++){
				if (swarm.at(i)->getID() == simpSwarm.at(j).id){
					//evaluate if the solution improved
					if (swarm.at(i)->getCurrentEvaluation() < simpSwarm.at(j).eval){
						S_i++;
						break;
					}
					else
						break;
				}
			}
		}
		return (S_i);
	}
}

/*Move the swarm to create new solutions */
void Swarm::moveSwarm(Configuration* config, const long int iteration) {
	static int sol_improved = 0;

	if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE)){
		cout << "\niteration: " << iteration << endl; //remove
		cout << "\tvar::swarm.size() = " << swarm.size() << endl;
	}

	//Update lBest of each particle
	for (unsigned int i=0;i<swarm.size();i++){
		if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_COMPUTATIONS))
			cout << "\tParticle [" << swarm.at(i)->getID() << "] -- lBest.at(t-1) [" << swarm.at(i)->getlBestID() << "] -- "; //<< endl; //remove
		//The id of the gBest particle depends on the topology and model of influence
		swarm.at(i)->getBestOfNeibourhood();
		if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_COMPUTATIONS))
			cout << "lBest.at(t) [" << swarm.at(i)->getlBestID() << "]" << endl;
	}
	if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_COMPUTATIONS)){
		if (config->getTopology() == TOP_HIERARCHICAL)
			printTree(config->getBranchingDegree());
		cout << endl;
	}

	//Self-explanatory
	clearResizeSimpSwarm(config, iteration);

	//The entire swarm is using the same omega1
	computeOmega1(config, iteration, -1, true); // -1 is a place holder for the id of the particle

	//Call this method here cause' some models of influence need to initialize things
	getInformants(config, -1, iteration);	// -1 is a place holder for the id of the particle

	//Compute the acceleration coefficients of the entire swarm
	computeAccelerationCoefficients(config, iteration);

	//Compute the random diagonal matrix using the increasing stochastic scaling grouping
	computeRandomMatrixGroupedIncreasing(config, iteration);

	if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE))
		cout << "\n\t------------------------------------------" << endl;

	for (unsigned int i=0; i<swarm.size(); i++){
		//Get the informants of i
		getInformants(config, i, iteration);

		if (config->verboseMode()){
			if (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE) {
				cout << "\tParticle [" << swarm.at(i)->getID() << "]" << " -- neighbors.size() = "
						<< swarm.at(i)->neighbors.size() << " -- informants.size() = " << swarm.at(i)->informants.size() << endl;
			}
			//print all neighbors
			if (config->verboseLevel() >= VERBOSE_LEVEL_COMPUTATIONS){
				cout << "\tNeighbors ids:  [ ";
				for (unsigned int j=0; j<swarm.at(i)->neighbors.size(); j++){
					cout << swarm.at(i)->neighbors.at(j)->getID() << " ";
				} cout << "]" << endl;
				//print the informants positions in swarm
				cout << "\tInformants pos: [ ";
				for (unsigned int j=0; j<swarm.at(i)->informants.size(); j++)
					cout << swarm.at(i)->informants.at(j) << " ";
				cout << "]" << endl;
				//print the informants IDs
				cout << "\tInformants ids: [ ";
				for (unsigned int j=0; j<swarm.at(i)->informants.size(); j++)
					cout << swarm.at(i)->neighbors.at(swarm.at(i)->informants.at(j))->getID() << " ";
				cout << "]"<< endl;
			}
		}
		//Here, computeOmega1 always receives the number of the particle and the flag = false
		swarm.at(i)->move(config, iteration,
				computeOmega1(config, iteration, i, false), //Compute self-adaptive omega1, otherwise this function returns the value already computed
				computeOmega2(config),
				computeOmega3(config),
				lastLevelComplete,
				sol_improved);
	}
	long double prev_Gbest_eval = global_best.eval; //best solution at iteration t-1
	//Update best_particle
	for (unsigned int i=0;i<swarm.size();i++){
		if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
			//Update best_particle of the Swarm (gBest) -- not to be confused with the lbest of a particle, which depends on the topology and MoI
			updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
			best_particle = swarm.at(i);
		}
	}

	if(config->useReinitialization())
		reinitializeParticlePosition(config, iteration); //There are two strategies available (one is commented)

	updatePerturbationVariables(config, prev_Gbest_eval, global_best.eval, iteration);

	//Compute the number of solutions that improved in the last iteration
	sol_improved = countImprovedSolutions(config, iteration);
}

void Swarm::computeRandomMatrixGroupedIncreasing(Configuration* config, long int iteration){
	double * Matrix;
	Matrix= new double[config->getProblemDimension()];

	if (config->getRandomMatrix() == MATRIX_GROUPED_INCREASING){

		int groups = ((double)(config->getProblemDimension()-1) / (double)(config->getMaxIterations()-1)) * (iteration-1) + 1;
		if (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE) {
			cout << " -- groups = " << groups << endl;
		}
		int dCounter = 0;
		double lastRandom = 0;

		//check that groups has a valid value
		if (groups<1) groups = 1;
		if (groups > config->getProblemDimension()) groups = config->getProblemDimension();

		//Create a random vector divided in groups
		for (int j=0; j<groups; j++){
			double randonNumber = RNG::randVal(0,1);
			for (int i=j*(config->getProblemDimension()/groups); i<(j*(config->getProblemDimension()/groups))+(config->getProblemDimension()/groups); i++){
				Matrix[i] =randonNumber;
				dCounter++;
				lastRandom=randonNumber;
			}
		}
		for (int k=dCounter; k<config->getProblemDimension(); k++)
			Matrix[k] =lastRandom;
	}
	else{
		for (int i=0; i<config->getProblemDimension(); i++)
			Matrix[i]=1.0;
	}
	if (config->verboseLevel() >= VERBOSE_LEVEL_COMPUTATIONS){
		cout << "\t Random diagonal matrix [";
		for (int i=0; i<config->getProblemDimension(); i++)
			cout << " " << Matrix[i];
		cout << "]" << endl;
	}
	config->setMatrixGI(Matrix);
}

void Swarm::getInformants(Configuration* config, int pPosinSwarm, long int iteration){
	if (pPosinSwarm == -1){
		//Only implement ranks if particles are not using them already
		if (config->getModelOfInfluence() == MOI_RANKED_FI && config->getOmega1CS() != IW_RANKS_BASED)
			rankParticles(simpSwarm);
	}
	else {
		//Best of neighborhood
		if (config->getModelOfInfluence() == MOI_BEST_OF_N){
			if (config->getTopology() == TOP_HIERARCHICAL){
				int Array_size = 0;		//variable to get the the size of Informants
				int TMP_Array[lastLevelComplete];
				getParticleParentsIDs(pPosinSwarm, TMP_Array); //Get parents of the particle

				for (int i=0; i<=lastLevelComplete; i++){
					if (TMP_Array[i] != -2){ //-2 indicates an empty position
						Array_size++;
					}
					else
						break;
				}

				//Clear the vector
				swarm.at(pPosinSwarm)->informants.clear();

				double bestEval = swarm.at(pPosinSwarm)->neighbors.at(0)->getPbestEvaluation();
				int bestPos = 0;
				int iter = 0;
				bool itselfAdded = false;
				//Find the position of the particle in neighbors vector that has the best PbestEvaluation
				//using the corresponding IDs in TMP_Array and neighbours.at(i)->getID()
				for (unsigned int i=0; i<swarm.at(pPosinSwarm)->neighbors.size(); i++){
					for (int j=0; j<Array_size; j++){
						if (swarm.at(pPosinSwarm)->neighbors.at(i)->getID() == TMP_Array[j]){
							if (iter == 0){
								bestPos = i;
								bestEval = swarm.at(pPosinSwarm)->neighbors.at(i)->getPbestEvaluation();
								iter++;
							}
							else{
								if (bestEval > swarm.at(pPosinSwarm)->neighbors.at(i)->getPbestEvaluation()){
									bestPos = i;
									bestEval = swarm.at(pPosinSwarm)->neighbors.at(i)->getPbestEvaluation();
									iter++;
								}
							}
						}
						//Add the particle itself to the the informants vector
						if (swarm.at(pPosinSwarm)->neighbors.at(i)->getID() == swarm.at(pPosinSwarm)->getID() && itselfAdded == false){
							swarm.at(pPosinSwarm)->informants.push_back(i);
							itselfAdded = true;
						}
					}
				}
				swarm.at(pPosinSwarm)->informants.push_back(bestPos);
			}
			else{
				//Clear the vector
				swarm.at(pPosinSwarm)->informants.clear();

				int bestID = swarm.at(pPosinSwarm)->getlBestID();

				//We need to find the position index of the IDs in TMP_Array in the particles neighbors vector
				for (unsigned int i=0; i<swarm.at(pPosinSwarm)->neighbors.size(); i++){
					if (swarm.at(pPosinSwarm)->neighbors.at(i)->getID() == bestID ||
							swarm.at(pPosinSwarm)->neighbors.at(i)->getID() == swarm.at(pPosinSwarm)->getID())
						swarm.at(pPosinSwarm)->informants.push_back(i);
				}
				if(swarm.at(pPosinSwarm)->informants.size() == 1)
					swarm.at(pPosinSwarm)->informants.push_back(swarm.at(pPosinSwarm)->informants[0]);
			}
		}
		//Fully informed
		if (config->getModelOfInfluence() == MOI_FI){
			if (config->getTopology() == TOP_HIERARCHICAL ){
				int Array_size = 0;		//variable to the the size of Informants
				int TMP_Array[lastLevelComplete];
				getParticleParentsIDs(pPosinSwarm, TMP_Array); //Get parents of the particle
				for (int i=0; i<=lastLevelComplete; i++){
					if (TMP_Array[i] != -2){ //-2 indicates an empty position
						Array_size++;
					}
					else
						break;
				}
				//Clear the vector
				swarm.at(pPosinSwarm)->informants.clear();

				//Find the index of the IDs in TMP_Array in the particles neighbors vector
				for (unsigned int i=0; i<swarm.at(pPosinSwarm)->neighbors.size(); i++){
					for (int j=0; j<Array_size; j++){
						if (swarm.at(pPosinSwarm)->neighbors.at(i)->getID() == TMP_Array[j])
							swarm.at(pPosinSwarm)->informants.push_back(i);
					}
					if (swarm.at(pPosinSwarm)->neighbors.at(i)->getID() == swarm.at(pPosinSwarm)->getID())
						swarm.at(pPosinSwarm)->informants.push_back(i);
				}
			}
			else {
				//Since some topologies are dynamic, the size of informants may change from iteration to iteration
				//Clear the vector
				swarm.at(pPosinSwarm)->informants.clear();

				for (unsigned int i=0;i<swarm.at(pPosinSwarm)->neighbors.size();i++){
					swarm.at(pPosinSwarm)->informants.push_back(i); //we use the indexes of neighbors
				}
			}
		}
		//Ranked fully informed
		if (config->getModelOfInfluence() == MOI_RANKED_FI){
			//Clear the vector
			swarm.at(pPosinSwarm)->informants.clear();

			//Container to sort the neighbors
			vector< vector<int> > TMP_vect;
			//Resize vector
			TMP_vect.resize((swarm.at(pPosinSwarm)->neighbors.size()), vector<int>(2));

			//This is the same as Fully informed
			for (unsigned int i=0;i<swarm.at(pPosinSwarm)->neighbors.size();i++){
				TMP_vect.at(i).at(0) = i; ///we use the indexes of neighbors
				TMP_vect.at(i).at(1) = swarm.at(pPosinSwarm)->neighbors.at(i)->getRanking();
			}

			//Sort by the second column (see sortcol prototype function above)
			sort(TMP_vect.begin(), TMP_vect.end(), sortcol);

			//copy informants ID to Informants sorted
			for (unsigned int i=0; i<TMP_vect.size(); i++){ //rows
				swarm.at(pPosinSwarm)->informants.push_back(TMP_vect[i][0]);
				TMP_vect[i].clear();
			}
			TMP_vect.clear();
		}
	}
}

void Swarm::resizeSwarm(Problem* problem, Configuration* config, long int iteration){
	static int iterWithoutImprvCounter=0;
	static int iterWithImprvCounter=0;
	static double prevGbest = LDBL_MAX;
	static long partRemoved = 0;
	static long partAdded = 0;
	int particleToAdd = 0;
	long int previous_size = swarm.size();

	if (config->getPopulationCS() == POP_TIME_VARYING){
		//First time: Just save the global_best value
		if (iteration == 1){
			prevGbest = global_best.eval;
		}
		//Check if it is needed to increase or decrease the swarm size
		else {
			//THERE IS AN IMPROVEMENT
			if (prevGbest > global_best.eval){
				prevGbest = global_best.eval;
				iterWithoutImprvCounter = 0;
				iterWithImprvCounter++;
				if (iterWithImprvCounter == config->getPopTViterations()){
					//Try to remove one particle from the swarm
					if ((long)swarm.size()-1 > config->getInitialPopSize()){
						removeOneParticle(config);
						updateTopologyConnections(config, -1, iteration);
						RNG::initializePermutation(swarm.size());
						partRemoved++;
					}
				}
			}
			//THERE IS NOT IMPROMEVENT
			else{
				iterWithImprvCounter=0;
				iterWithoutImprvCounter++;
				if (iterWithoutImprvCounter == config->getPopTViterations()){
					//Try to add one particle to the swarm
					if ((long)swarm.size()+1 < config->getFinalPopSize()){
						//If there is space in the swarm: add one particle
						addParticles(problem, config, 1, iteration);
						updateTopologyConnections(config, previous_size, iteration);
						RNG::initializePermutation(swarm.size());
						partAdded++;
					}
					else {
						//Reinitialize the worst particle
						int worstP_pos = 0;
						double worstP_eval = swarm.at(0)->getPbestEvaluation();
						//Find the position in the swarm of the worst particle
						for(unsigned int i=1; i<swarm.size(); i++){
							if (swarm.at(i)->getPbestEvaluation() >= worstP_eval){
								worstP_pos = i;
								worstP_eval = swarm.at(i)->getPbestEvaluation();
							}
						}
						swarm.at(worstP_pos)->initializePosition(config,100, false);
					}
				}
			}
		}
	}
	if (config->getPopulationCS() == POP_INCREMENTAL){
		particleToAdd = config->getParticlesToAdd();
		//Add k particle per iteration
		if (previous_size+particleToAdd <= config->getFinalPopSize()){
			addParticles(problem, config, particleToAdd, iteration);
			updateTopologyConnections(config, previous_size, iteration);
			RNG::initializePermutation(swarm.size());
		}
		else{
			//See if we can add the difference
			particleToAdd = config->getFinalPopSize()-previous_size;
			if ( particleToAdd > 0){
				addParticles(problem, config, particleToAdd, iteration);
				updateTopologyConnections( config, previous_size, iteration);
				RNG::initializePermutation(swarm.size());
			}
		}
	}
}

void Swarm::updateTopologyConnections(Configuration* config, long previous_size, long int iteration){

	if (config->getTopology() == TOP_FULLYCONNECTED){
		//Clear connections
		for(unsigned int i=0; i<swarm.size(); i++)
			swarm.at(i)->neighbors.clear();
		//Create new connections
		createFullyConnectedTopology();
	}
	if (config->getTopology() == TOP_RING){
		//Reconnect the whole swarm
		for(unsigned int i=0; i<swarm.size(); i++)
			swarm.at(i)->neighbors.clear();
		createRingTopology();
	}
	if (config->getTopology() == TOP_VONNEUMANN){
		//Reconnect the whole swarm
		for(unsigned int i=0; i<swarm.size(); i++)
			swarm.at(i)->neighbors.clear();
		createVonNeumannTopology();
	}

	if (previous_size > 0) { //PARTICLES WERE ADDED -- EXTEND THE TOPOLOGY
		if (config->getTopology() == TOP_HIERARCHICAL){
			for(unsigned int i=previous_size; i<swarm.size(); i++){
				for(unsigned int j=0; j<swarm.size(); j++){
					swarm.at(i)->addNeighbour(swarm.at(j));
					if (j<previous_size)
						swarm.at(j)->addNeighbour(swarm.at(i));
				}
			}
			updateHierarchical(config->getBranchingDegree(),previous_size);
		}
		if (config->getTopology() == TOP_WHEEL){
			for(unsigned int i=previous_size; i<swarm.size(); i++){
				swarm.at(i)->addNeighbour(swarm.at(i));
				swarm.at(i)->addNeighbour(swarm.at(0));
				swarm.at(0)->addNeighbour(swarm.at(i));
			}
		}
		if (config->getTopology() == TOP_RANDOM){
			long int randomEdge;

			//First: make all particles neighbors of themselves
			for(unsigned int i=previous_size; i<swarm.size(); i++)
				swarm.at(i)->addNeighbour(swarm.at(i));

			//Second: Add a random neighbor
			for(unsigned int i=previous_size; i<swarm.size(); i++){
				for(unsigned int j=0;j<swarm.size();j++){
					randomEdge = (int)floor(RNG::randVal(0.0,(double)swarm.size()));

					if (randomEdge == i){
						continue;
					}
					else{
						swarm.at(i)->addNeighbour(swarm.at(randomEdge));
						break;
					}
				}
			}
		}
		if (config->getTopology() == TOP_TIMEVARYING){
			if (iteration < config->getTopologySchedule()){
				//We connect a particle randomly with swarm ensuring that the particle is connected to the adjacent neighbors (RING)

				//Get average connections number of the swarm.
				int averageConnections = 0;
				for(unsigned int i=0; i<previous_size; i++){
					averageConnections += swarm.at(i)->neighbors.size();
				}
				averageConnections = (int)floor((double)averageConnections/previous_size);

				unsigned int a;
				unsigned int b;
				for(unsigned int i=previous_size; i<swarm.size(); i++){
					//First connect the new particles in RING
					a=i-1;
					b=i+1;
					if(i==0)
						a=swarm.size()-1;
					if(i==(swarm.size()-1)){
						b=0;
						//Reconnect first particle with the last particle
						swarm.at(0)->addNeighbour(swarm.at(i));
						//Reconnect previous last particle with the next one
						swarm.at(previous_size-1)->addNeighbour(swarm.at(previous_size));
					}
					swarm.at(i)->addNeighbour(swarm.at(i));
					swarm.at(i)->addNeighbour(swarm.at(a));
					swarm.at(i)->addNeighbour(swarm.at(b));

					//Create #averageConnections random connections
					for(int j=0; j<averageConnections-3; j++){
						unsigned int randomEdge = (unsigned int)floor(RNG::randVal(0.0,(double)swarm.size()));
						if (randomEdge != a && randomEdge != b && randomEdge != i){
							swarm.at(i)->addNeighbour(swarm.at(randomEdge));
							swarm.at(randomEdge)->addNeighbour(swarm.at(i));
						}
					}
				}
			}
			else { //RING -- After the topology update period is finished
				for(unsigned int i=0; i<swarm.size(); i++){
					swarm.at(i)->neighbors.clear();
				}
				createRingTopology();
			}
		}
	}
	//PARTICLES WERE REMOVED -- RECONNECT THE TOPOLOGY
	else {
		//Reconnect the swarm
		if (config->getTopology() == TOP_HIERARCHICAL){
			//Clear old connections
			for(unsigned int i=0; i<swarm.size(); i++){
				swarm.at(i)->neighbors.clear();
			}
			//Create new fully connected topology
			createHierarchical(config->getBranchingDegree(), config->getFinalPopSize());
			for (int j=0; j<=lastLevelComplete; j++)
				updateTree(config->getBranchingDegree());
		}
		if (config->getTopology() == TOP_WHEEL){
			//Clear old connections
			for(unsigned int i=0; i<swarm.size(); i++){
				swarm.at(i)->neighbors.clear();
			}
			//Create new wheel topology
			createWheelTopology();
		}
		if (config->getTopology() == TOP_RANDOM){
			//Clear old connections
			for(unsigned int i=0; i<swarm.size(); i++){
				swarm.at(i)->neighbors.clear();
			}
			//Create new random edge topology
			createRandomEdge();
		}
		if (config->getTopology() == TOP_TIMEVARYING){
			if (iteration < config->getTopologySchedule()){
				//Connect a particle randomly with swarm ensuring that the particle is connected to the adjacent neighbors (RING)
				//Get average connections number of the swarm.
				int averageConnections = 0;
				for(unsigned int i=0; i<swarm.size(); i++){
					averageConnections += swarm.at(i)->neighbors.size();
				}
				averageConnections = (int)floor((double)averageConnections/swarm.size()-1);

				//Clear old connections
				for(unsigned int i=0; i<swarm.size(); i++){
					swarm.at(i)->neighbors.clear();
				}

				unsigned int a;
				unsigned int b;
				for(unsigned int i=0; i<swarm.size(); i++){
					//First connect the new particles in RING
					a=i-1;
					b=i+1;
					if(i==0)
						a=swarm.size()-1;
					if(i==(swarm.size()-1))
						b=0;

					swarm.at(i)->addNeighbour(swarm.at(i));
					swarm.at(i)->addNeighbour(swarm.at(a));
					swarm.at(i)->addNeighbour(swarm.at(b));

					bool noAdd = false;
					//After connect the new particles with #averageConnections random neighbors
					for(int j=0; j<averageConnections-3; j++){
						unsigned int randomEdge = (unsigned int)floor(RNG::randVal(0.0,(double)swarm.size()));
						if (randomEdge != a && randomEdge != b && randomEdge !=i){
							for(unsigned int j=0; j<swarm.at(i)->neighbors.size(); j++){
								if (swarm.at(randomEdge)->getID() == swarm.at(i)->neighbors.at(j)->getID())
									noAdd = true;
							}
							if (noAdd == false){
								swarm.at(i)->addNeighbour(swarm.at(randomEdge));
								swarm.at(randomEdge)->addNeighbour(swarm.at(i));
							}
						}
					}
				}
			}
			else { //RING -- After the topology update period is finished
				//Clear old connections
				for(unsigned int i=0; i<swarm.size(); i++){
					swarm.at(i)->neighbors.clear();
				}
				createRingTopology();
			}
		}
	}
}

void Swarm::removeOneParticle(Configuration* config){
	//Add particles to the swarm
	int pID = swarm.at(0)->getID();
	double pQuality = swarm.at(0)->getPbestEvaluation();

	//Find worst particle in the last level if TOP_HIERARCHICAL is being used
	if (config->getTopology() == TOP_HIERARCHICAL){
		int childsInLastLevel = 0;
		//Get the number of particles in the last level complete
		for (int i=0; i<pow(config->getBranchingDegree(),lastLevelComplete+1); i++){
			if (hierarchy.at(lastLevelComplete+1).at(i) != -2) {
				//Since hierarchy contains the ID of particles, look for the particle with that ID
				for (unsigned int j=0;j<swarm.size();j++){
					if (swarm.at(j)->getID() == hierarchy.at(lastLevelComplete+1).at(i)){
						//If its the first iteration there's nothing to compare, just keep the quality in pQuality
						if (i>0){
							pQuality = swarm.at(j)->getPbestEvaluation();
							pID = hierarchy.at(lastLevelComplete+1).at(i);
						}
						else {
							if (swarm.at(j)->getPbestEvaluation() >= pQuality){
								pID = hierarchy.at(lastLevelComplete+1).at(i);
								pQuality = swarm.at(j)->getPbestEvaluation();
							}
						}
					}
				}
				childsInLastLevel++;
			}
		}
		//Remove particle ID from the hierarchy
		for (int i=0; i<pow(config->getBranchingDegree(),lastLevelComplete+1); i++){
			if (hierarchy.at(lastLevelComplete+1).at(i) == pID) {
				hierarchy.at(lastLevelComplete+1).at(i) = -2;
				break;
			}
		}

		//Update lastLevelComplete in case there are no more particles in that level
		if(childsInLastLevel == 1)
			lastLevelComplete--;
	}
	//Find worst particle
	else {
		for (unsigned int i=1;i<swarm.size();i++){
			if (swarm.at(i)->getPbestEvaluation() >=  pQuality){
				pID = swarm.at(i)->getID();
				pQuality = swarm.at(i)->getPbestEvaluation();
			}
		}
	}
	//Remove particle from the neighbors set of other particles
	for (unsigned int i=0;i<swarm.size();i++){
		swarm.at(i)->eraseNeighborbyID(pID);
	}
	//Remove particle from the swarm
	for(unsigned int i=0;i<swarm.size();i++){
		if(pID == swarm.at(i)->getID()){
			swarm.erase(swarm.begin()+i);
			break;
		}
	}
	//Update variable size
	config->setSwarmSize(swarm.size()); //long int particles


	for (unsigned int i=0;i<swarm.size();i++){
		swarm.at(i)->setID(i);
	}
}


void Swarm::addParticles(Problem* problem, Configuration* config, int numOfParticles, long int iteration){
	if ((long)swarm.size() <= config->getFinalPopSize()){
		//Add particles to the swarm
		int current_size = swarm.size();
		long int partID = 0;

		//Find the particle with the largest ID number
		for(unsigned int j=0;j<swarm.size();j++){
			if (swarm.at(j)->getID() > partID)
				partID = swarm.at(j)->getID();
		}

		for (long int i=current_size; i<current_size+numOfParticles; i++){
			partID++;
			Particle* aParticle = new Particle(problem, config, partID, iteration);
			swarm.push_back(aParticle);
			if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
				updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
				best_particle = swarm.at(i);
			}
		}
		//Update variable size
		config->setSwarmSize(swarm.size()); //long int particles
	}
}

void Swarm::reinitializeParticlePosition(Configuration* config, long int iteration){
	/*
	 * Reinitialize particles that are too close to gbest
	 */

	//First strategy: using standard deviation and overall change in the O.F.
	static long int reinitSchedule = 10*config->getProblemDimension()/swarm.size();
	static long double previousgBestEval = best_particle->getPbestEvaluation();
	double stdSwarm = 0;
	double mean = 0;
	if (reinitSchedule == 0)
		reinitSchedule=1;

	//Mean value of the O.F. value of the swarm
	for (unsigned int i=0;i<swarm.size();i++)
		mean += swarm.at(i)->getPbestEvaluation();
	mean = mean/swarm.size();

	//Standard deviation of the O.F. value of the swarm
	for (unsigned int i=0;i<swarm.size();i++)
		stdSwarm += pow(swarm.at(i)->getPbestEvaluation()-mean,2);
	stdSwarm = sqrt(stdSwarm/swarm.size());

	//Reinitialize particles if stdSwarm is lower than REINIT_PRECISION
	if (stdSwarm < REINIT_PRECISION){
		if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE))
			cout << "\n\tRestarting particles using random values in bounds\n";
		for (unsigned int i=0;i<swarm.size();i++){
			if (swarm.at(i)->getID() != best_particle->getID()) //exclude gBest
				swarm.at(i)->setRandomPositionInBoundsWithProbability(config);
		}
	}

	//Check if Overall change in the O.F: is below OBJECTIVE_FUNCTION_CHANGE_THRESHOLD
	if (iteration%reinitSchedule == 0){
		//Reinitialize particles if the best solution has not improved after reinitSchedule iterations
		if (fabs(previousgBestEval - best_particle->getPbestEvaluation()) < OBJECTIVE_FUNCTION_CHANGE_THRESHOLD){
			//cout << "\nSo far, so good"<< endl;

			if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_COMPUTATIONS)){
				cout << "\n\tRestarting particles using gbest derivative\n";
				cout << "\t\tvar::iteration " << iteration << endl;
				cout << "\t\tvar::reinitSchedule " << reinitSchedule << endl;
				cout << "\t\tvar::previousgBestEval " << previousgBestEval << endl;
			}
			for (unsigned int i=0;i<swarm.size();i++){
				if (swarm.at(i)->getID() != best_particle->getID()){ //exclude gBest
					swarm.at(i)->setRandomPositiongBestDerivative(config, best_particle->getPbestPosition());
				}
			}
		}
		previousgBestEval = best_particle->getPbestEvaluation();
	}

	// Second strategy: using similarity (It works, you can uncomment it!)
	//	for (unsigned int i=0;i<swarm.size();i++){
	//		//Check that the particle is different to itself before reinitializing its position
	//		if (swarm.at(i)->getID() != swarm.at(i)->getgBestID()){
	//			double* gBestPos = 0;
	//			double* particlePos = swarm.at(i)->getCurrentPosition();
	//			double similarity = 0; //distance between x and gBest
	//			//Find gBest of particle by id
	//			for (unsigned int r=0;r<swarm.size();r++){
	//				if (swarm.at(i)->getgBestID() == swarm.at(r)->getID()){
	//					gBestPos = swarm.at(r)->getCurrentPosition();
	//				}
	//			}
	//			//Compute similarity between particle and gBest
	//			for (int j=0;j<config->getProblemDimension();j++){
	//				similarity +=  particlePos[j] / gBestPos[j];
	//			}
	//			similarity = similarity/config->getProblemDimension();
	//			if (fabs(1-similarity) < REINIT_PRECISION){
	//				//Reinitialize particle to a random position and set v=0
	//				swarm.at(i)->initializePosition(config,100, false);
	//			}
	//		}
	//	}
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
	for(unsigned int i=0; i< dimensions; i++){
		cout << fixed << getGlobalBest().x[i] << "  ";
	}
	cout << " ]\n" << endl;
}

Solution Swarm::getGlobalBest(){
	return (global_best);
}

void Swarm::decomposePhi2(int modelOfInflu, int part_id, int numInformants){
	if (modelOfInflu == MOI_FI)
		swarm.at(part_id)->setPhi2(swarm.at(part_id)->getPhi2()/numInformants);
}

void Swarm::computeAccelerationCoefficients(Configuration* config, long int iteration){
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
			swarm.at(i)->setPhi1(RNG::randVal(config->getInitialPhi1(),config->getFinalPhi1()));
			swarm.at(i)->setPhi2(RNG::randVal(config->getInitialPhi2(),config->getFinalPhi2()));
		}
	} break;
	case AC_EXTRAPOLATED:{
		for (unsigned int i=0; i<swarm.size(); i++){
			double varPhi1 = exp((double)(-iteration/config->getMaxIterations()));
			double distanceToGbest = (swarm.at(swarm.at(i)->getlBestID())->getCurrentEvaluation()-swarm.at(i)->getCurrentEvaluation())/
					swarm.at(swarm.at(i)->getlBestID())->getCurrentEvaluation();
			swarm.at(i)->setPhi1(varPhi1);
			swarm.at(i)->setPhi2(exp(varPhi1*distanceToGbest));
		}
	} break;
	case AC_TIME_VARYING:{
		double varPhi1 = (config->getInitialPhi1()-config->getFinalPhi1()) *
				(double)(iteration/config->getMaxIterations()) + config->getInitialPhi1();
		double varPhi2 = (config->getInitialPhi2()-config->getFinalPhi2()) *
				(double)(iteration/config->getMaxIterations()) + config->getInitialPhi2();
		for (unsigned int i=0; i<swarm.size(); i++){
			swarm.at(i)->setPhi1( varPhi1 );
			swarm.at(i)->setPhi2( varPhi2 );
		}
	} break;
	}

	if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE))
		cout << "\tvarPhi1 = " << config->getPhi1() << " -- varPhi2 = " << config->getPhi2();
}

void Swarm::updatePerturbationVariables(Configuration* config, double previousGbest_eval, double currentGbest_eval, long int iteration){

	if (config->getMagnitude1CS() == MAGNITUDE_SUCCESS){
		//Success
		if (previousGbest_eval != currentGbest_eval){
			config->set_mag1_sc(config->get_mag1_sc()+1); //increase successes
			config->set_mag1_fc(0); //set failures to 0
		}
		//Failure
		else{
			config->set_mag1_fc(config->get_mag1_fc()+1); //increase failures
			config->set_mag1_sc(0); //set successes to 0
		}

		//Update magnitude
		if (config->get_mag1_sc() > config->getMag1_parm_success()){
			config->setMagnitude1(config->getMagnitude1() * 2.0);
		}
		if (config->get_mag1_fc() > config->getMag1_parm_failure()){
			config->setMagnitude1(config->getMagnitude1() * 0.5);
		}
	}

	if (config->getMagnitude2CS() == MAGNITUDE_SUCCESS){
		//Success
		if (previousGbest_eval != currentGbest_eval){
			config->set_mag2_sc(config->get_mag2_sc()+1); //increase successes
			config->set_mag2_fc(0); //set failures to 0
		}
		//Failure
		else{
			config->set_mag2_fc(config->get_mag2_fc()+1); //increase failures
			config->set_mag2_sc(0); //set successes to 0
		}

		//Update value of alpha
		//Increase magnitude
		if (config->get_mag2_sc() > config->getMag2_parm_success()){
			if (config->getPerturbation2CS() == PERT2_RECTANGULAR || config->getPerturbation2CS() == PERT2_NOISY)
				config->setMagnitude2(config->getMagnitude2()*2.0);
		}
		//Decrease magnitude
		if (config->get_mag2_fc() > config->getMag2_parm_failure()){
			if (config->getPerturbation2CS() == PERT2_RECTANGULAR || config->getPerturbation2CS() == PERT2_NOISY)
				config->setMagnitude2(config->getMagnitude2()*0.5);
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
	if (swarm.size() > 2){
		int a;
		int b;
		for(unsigned int i=0;i<swarm.size();i++){
			a=i-1;
			b=i+1;
			if(i==0)
				a=swarm.size()-1;
			if(i==(swarm.size()-1))
				b=0;

			swarm.at(i)->addNeighbour(swarm.at(i));
			swarm.at(i)->addNeighbour(swarm.at(a));
			swarm.at(i)->addNeighbour(swarm.at(b));
		}
	}
	else {
		for(unsigned int i=0;i<swarm.size();i++){
			for(unsigned int j=0;j<swarm.size();j++){
				swarm.at(i)->addNeighbour(swarm.at(j));
			}
		}
	}
}

void Swarm::createWheelTopology(){        //Particles are neighbors of one central particle
	swarm.at(0)->addNeighbour(swarm.at(0));
	for(unsigned int i=1;i<swarm.size();i++){
		swarm.at(i)->addNeighbour(swarm.at(i));
		swarm.at(i)->addNeighbour(swarm.at(0));
		swarm.at(0)->addNeighbour(swarm.at(i));
	}
}

void Swarm::createRandomEdge(){			//Random edge topology
	long int randomEdge;

	//First: make each particle a neighbor of itself
	for(unsigned int i=0; i<swarm.size(); i++)
		swarm.at(i)->addNeighbour(swarm.at(i));

	//Second: Add a random neighbor
	for(unsigned int i=0; i<swarm.size(); i++){
		for(unsigned int j=0;j<swarm.size();j++){
			randomEdge = (int)floor(RNG::randVal(0.0,(double)swarm.size()));

			if (randomEdge == i){
				continue;
			}
			else{
				swarm.at(i)->addNeighbour(swarm.at(randomEdge));
				break;
			}
		}
	}
}

void Swarm::createVonNeumannTopology(){
	if (swarm.size() > 3){
		int a;
		int b;
		int c;
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
			swarm.at(i)->addNeighbour(swarm.at(i));
			swarm.at(i)->addNeighbour(swarm.at(a));
			swarm.at(i)->addNeighbour(swarm.at(b));
			swarm.at(i)->addNeighbour(swarm.at(c));
		}
	}
	else {
		for(unsigned int i=0;i<swarm.size();i++){
			for(unsigned int j=0;j<swarm.size();j++){
				swarm.at(i)->addNeighbour(swarm.at(j));
			}
		}
	}
}

void Swarm::updateTimeVaryingTopology(Configuration* config, long int iterations){
	long int swarm_size;
	bool emergyExit = false;
	(config->getPopulationCS() == POP_CONSTANT) ? swarm_size = swarm.size() : swarm_size = config->getFinalPopSize();

	if (config->getTopologyUpdatePeriod() == 0) //Swarm too small
		return;

	//Topology update: Following the progression n-2, n-3, ..., 2. (see esteps variable)
	if ((iterations > 0) &&	(iterations%config->getTopologyUpdatePeriod() == 0) &&
			((config->getEsteps() < swarm_size-3 && (config->getPopulationCS() == POP_CONSTANT))
					|| ((config->getEsteps() > 0 && swarm.size() > 4 ) && config->getPopulationCS() != POP_CONSTANT))){

		unsigned int removals = 0;
		RNG::shufflePermutation();
		unsigned int target = swarm_size-(2+config->getEsteps());

		double averageConnections = 0;
		for(unsigned int i=0; i<swarm.size(); i++){
			averageConnections += swarm.at(i)->neighbors.size();
		}
		averageConnections = averageConnections/swarm.size();

		while(removals < target && averageConnections > 3 && emergyExit == false){
			try {// vector::at throws an out-of-range
				for(unsigned int i=0;i<swarm.size();i++){
					int particleIndex = RNG::getPermutationElement(i);
					int neighborID = -1;

					if (swarm.at(particleIndex)->neighbors.size() > 3 ){ //3 because a particle is a neighbor to itself
						neighborID = swarm.at(particleIndex)->getRandomNonAdjacentNeighborID(config);
						swarm.at(particleIndex)->eraseNeighborbyID(neighborID);
						swarm.at(neighborID)->eraseNeighborbyID(particleIndex);
						removals++;
					}
					if (i== swarm.size()-1)
						emergyExit = true;
					if (removals == target )
						break;
				}
			}
			catch (const std::out_of_range& oor) {
				break;
			}
		}

		(config->getPopulationCS() == POP_CONSTANT || config->getPopulationCS() == POP_TIME_VARYING) ? config->setEsteps(config->getEsteps()+1) : config->setEsteps(config->getEsteps()-1);
	}
}

void Swarm::createHierarchical(int branching, int finalPopSize){
	//The topology is fully-connected, but we use a hierarchy as a model of influence
	long int firstPart = 1;		//id of the last particle that can be added in a complete level without exceeding size
	int actualSpace=1;			//actual number of nodes than can be stored in a given level. This depends on the swarm size and the branching degree
	int maxWidth = 0;
	int temp_LastLevelComplete = 0;
	long int temp_firstPart = 1;
	lastLevelComplete = 0;

	//Find last level that can be filled in the tree with max width and the id of the last particle
	while (swarm.size() - pow(branching,lastLevelComplete) >= pow(branching,lastLevelComplete+1)) {
		lastLevelComplete++;
		firstPart+=pow(branching,lastLevelComplete);
	}
	//It can be the case that for a very large branching degree the firstPart variable
	//is larger than the size of the swarm. In that case, find the last level complete
	//without exceeding the swarm size
	while ((long unsigned)firstPart > swarm.size()){
		lastLevelComplete--;
		firstPart-=pow(branching,lastLevelComplete);
	}
	//Compute actual space
	for (int i=1; i<=lastLevelComplete; ++i)
		actualSpace+=pow(branching,i);
	if (firstPart > actualSpace)
		firstPart = actualSpace;
	//Constant population size
	if ((long unsigned)finalPopSize == swarm.size()){
		maxWidth = pow(branching,lastLevelComplete+1);
		temp_LastLevelComplete = lastLevelComplete;
		temp_firstPart = firstPart;
	}
	//Dynamic population size
	else {
		//recompute everything to get the maxWidth, which is used to allocate memory
		while (finalPopSize - pow(branching,temp_LastLevelComplete) >= pow(branching,temp_LastLevelComplete+1)){
			temp_LastLevelComplete++;
			temp_firstPart+=pow(branching,temp_LastLevelComplete);
		}
		while (temp_firstPart > finalPopSize){
			temp_LastLevelComplete--;
			temp_firstPart-=pow(branching,temp_LastLevelComplete);
		}
		maxWidth = pow(branching,temp_LastLevelComplete+1);
	}

	//Create a fully connected topology first
	for(unsigned int i=0;i<swarm.size()	;i++){
		for(unsigned int j=0;j<swarm.size();j++){
			swarm.at(i)->addNeighbour(swarm.at(j));
		}
	}

	//Initialize tree structure
	//Note that even when the population is dynamic, hierarchy already has the right size to allocate all
	//the particles in the population
	//temp_LastLevelComplete+1 indicates the levels (depth of the tree)
	//maxWidth indicates maximum number of nodes in a level (max width of the tree)
	hierarchy.resize(temp_LastLevelComplete+2, vector<int>(maxWidth));

	for (int i=0; i<=temp_LastLevelComplete+1; i++){
		for (int j=0; j<maxWidth; j++)
			hierarchy.at(i).at(j) = -2;				//-2 indicates that the position is not being used
	}

	//tree variables
	int d = 0;				//node's degree
	int nodesCounter = 0; 	//every time we reach the max branching degree, we move to the next node
	int h = 1;				//height of the tree
	int width = 0;			//width of the level

	//root node
	int parent = hierarchy.at(0).at(0) = swarm.at(0)->getID();
	swarm.at(0)->setParent(-1); //the root particle has parent -1
	//We first fill all the levels that will be complete, i.e., the levels with maximum width
	for(unsigned int i=1;i<firstPart;i++){
		if (width < pow(branching,h)) {
			hierarchy.at(h).at(width) = swarm.at(i)->getID();
			swarm.at(i)->setParent(parent);
			width++;	//tree's degree
			d++;		//node's degree
			if (d==branching) {
				nodesCounter++;
				parent = hierarchy.at(h-1).at(nodesCounter);
				d=0;
			}
			if (width == pow (branching,h)){
				h++;
				width=0;
				nodesCounter=0;
				parent = hierarchy.at(h-1).at(nodesCounter);
			}
		}
	}
	//reset variables of the tree for the rest of the particles
	width = 0;
	nodesCounter=0;
	int iterCount = 0;
	parent = hierarchy.at(lastLevelComplete).at(nodesCounter);
	//Now we add one d at a time to each node in the last level complete
	for(unsigned int i=firstPart;i<swarm.size();i++){
		hierarchy.at(lastLevelComplete+1).at(width+(nodesCounter*(branching-1))+iterCount) = swarm.at(i)->getID();
		swarm.at(i)->setParent(parent);
		width++;		//last level degree
		nodesCounter++; //we add one d to each node to have a fair distribution
		parent = hierarchy.at(lastLevelComplete).at(nodesCounter);
		if (width == pow (branching,lastLevelComplete)){
			width=0;
			nodesCounter=0;
			iterCount++;
			parent = hierarchy.at(lastLevelComplete).at(nodesCounter);
		}
	}
}

void Swarm::updateHierarchical(int branching, long previous_size){
	int childsInLastLevel = 0;
	int nodesToAdd = swarm.size()-previous_size;

	//Get the number of particles in the last level complete
	for (int i=0; i<pow(branching,lastLevelComplete+1); i++){
		if (hierarchy.at(lastLevelComplete+1).at(i) != -2) {
			childsInLastLevel++;
		}
	}

	//The number of particles to be added fits in lastLevelComplete
	if (pow(branching,lastLevelComplete+1)-childsInLastLevel >= nodesToAdd){
		addParticlesInLastLevel(previous_size, swarm.size(), branching);
		if (nodesToAdd+childsInLastLevel == pow(branching,lastLevelComplete+1)){
			lastLevelComplete++;
		}
	}
	//The number of particles to add is larger than the space in lastLevelComplete
	else {
		//Fill the last level
		long firstPart = (pow(branching,lastLevelComplete+1)-childsInLastLevel);
		addParticlesInLastLevel(previous_size, previous_size+firstPart, branching);
		lastLevelComplete++;

		//Compute the number of nodes that still have to be added
		nodesToAdd = nodesToAdd-firstPart;
		int newLastLevel = lastLevelComplete;

		//Find last level that can be filled in the tree with max width and the id of the last particle
		while (nodesToAdd - pow(branching,newLastLevel) >= pow(branching,newLastLevel+1))
			newLastLevel++;

		//Verify that the tree size does not exceed the number of particles
		long unsigned int treeSize =0;
		for (int c=0; c<=newLastLevel; c++)
			treeSize += pow(branching,c);
		while (treeSize > swarm.size()){
			treeSize -= pow(branching,newLastLevel);
			newLastLevel--;
		}

		//First and last are the ids of the particles that are going to be added in the level
		long first = previous_size+firstPart;
		long last = 0;
		if (nodesToAdd < pow(branching,lastLevelComplete+1))
			last = first+nodesToAdd;
		else
			last = first+pow(branching,lastLevelComplete+1);


		for (int j=lastLevelComplete; j<=newLastLevel; j++){
			addParticlesInLastLevel(first, last, branching);
			if (j+1 <=newLastLevel){
				lastLevelComplete++;
				nodesToAdd = nodesToAdd-(last-first);
				first = last;
				if (nodesToAdd < pow(branching,lastLevelComplete+1))
					last = first + nodesToAdd;
				else
					last = first+pow(branching,lastLevelComplete+1);
			}
		}
		lastLevelComplete = newLastLevel;
	}
}

void Swarm::addParticlesInLastLevel(int first, int last, int branching){
	//tree variables
	int width = 0;
	int nodesCounter=0;
	int iterCount = 0;
	int parent = hierarchy.at(lastLevelComplete).at(nodesCounter);
	bool particleAdded = false;

	//Add one node at a time to each node in the last level complete
	for (int i=first; i<last; ++i){	//particles to add
		for (int j=0; j<pow(branching,lastLevelComplete+1); j++){ //maxWidth of the level
			if (hierarchy.at(lastLevelComplete+1).at(width+(nodesCounter*(branching-1))+iterCount) == -2){ //look for free positions
				hierarchy.at(lastLevelComplete+1).at(width+(nodesCounter*(branching-1))+iterCount) = swarm.at(i)->getID();
				swarm.at(i)->setParent(parent);
				particleAdded = true;
			}
			width++;		//last level degree
			nodesCounter++; //we add one node to each parent node to have a fair distribution
			parent = hierarchy.at(lastLevelComplete).at(nodesCounter);
			if (width == pow(branching,lastLevelComplete)){
				width=0;
				nodesCounter=0;
				iterCount++;
				parent = hierarchy.at(lastLevelComplete).at(nodesCounter);
			}
			if (particleAdded){
				break;
				particleAdded = false;
			}
		}
	}
}

void Swarm::getParticleParentsIDs(int particleID, int *ParentsArray1D){
	//Initialize the array to -2
	for (int i=0; i<=lastLevelComplete; i++)
		ParentsArray1D[i] = -2;

	//Find particle's index in the swarm by its ID
	long int index=0;
	for(unsigned int i=0; i<swarm.size(); i++){
		if (particleID == swarm.at(i)->getID()){
			index = i;
			break;
		}
	}

	//Get direct parentNode ID
	int aParent = swarm.at(index)->getParent();
	//This is the case of the particle in the root node
	if (aParent==-1){
		ParentsArray1D[0] = swarm.at(index)->getID();
		swarm.at(particleID)->setlBestID(swarm.at(index)->getID());
	}
	else{
		int pos = 0;
		while (aParent != -1){
			ParentsArray1D[pos] = aParent;
			pos++;
			aParent = swarm.at(aParent)->getParent();
		}
	}
}

int Swarm::getParticleNumParents(int particleID){
	int partents_counter=0;

	//Find particle index in the swarm by its ID
	long int index=0;
	for(unsigned int i=0; i<swarm.size(); i++){
		if (particleID == swarm.at(i)->getID()){
			index = i;
			break;
		}
	}

	int aParent = swarm.at(index)->getParent();
	while (aParent != -1){
		partents_counter++;
		aParent = swarm.at(aParent)->getParent();
	}
	return (partents_counter);
}


void Swarm::updateTree(int branching){
	//Traverse the tree and update the position of the particles
	//reset indexes to traverse the entire tree
	int h = 1;
	int width = 0;
	int iterCount = 0;
	int childCounter = 1;
	int newParent = 0;
	int newH = 0;
	int newWidth = 0;
	int parentNode = swarm.at(hierarchy.at(0).at(0))->getParent();	//start at the root top-bottom
	int parentH = 0;
	int parentWidth = -1;
	for(unsigned int i=1;i<swarm.size();i++){
		if (h <= lastLevelComplete){
			if (width < pow(branching,h)){ //max level width
				if (parentNode!=swarm.at(hierarchy.at(h).at(width))->getParent()){
					parentNode=swarm.at(hierarchy.at(h).at(width))->getParent();
					parentH = h-1;
					parentWidth++;
				}
				if (swarm.at(hierarchy.at(h).at(width))->getCurrentEvaluation() < swarm.at(parentNode)->getCurrentEvaluation()){
					newParent = swarm.at(hierarchy.at(h).at(width))->getID();	//newParent coordinates (this is the child that becomes parent)
					newH = h;
					newWidth = width;
					//we change places only when newParent is different from getParent() of the last particle
					//this means that at least one of the children will become parent.
					//The id of this particle is saved in newParent
					if (width > 0 && (width+1) % branching == 0 && (newParent != swarm.at(hierarchy.at(h).at(width))->getParent())){
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
			}
		}
		else {
			parentNode = hierarchy.at(h-1).at(iterCount);
			parentH = h-1;
			parentWidth++;
			newParent = parentNode;

			for(int j=iterCount*branching;j<(iterCount*branching)+branching;j++){
				try{
					if (swarm.at(hierarchy.at(h).at(j))->getParent() == parentNode) {
						childCounter++;
					}
					//Check if there is a new parent
					if (swarm.at(hierarchy.at(h).at(j))->getCurrentEvaluation() < swarm.at(parentNode)->getCurrentEvaluation()){
						newParent = swarm.at(hierarchy.at(h).at(j))->getID();	//newParent coordinates (this is the  that becomes parent)
						newH = h;
						newWidth = j;
					}
					//Switch parent and child after comparing with all children
					if (j == (iterCount*branching)+branching-1 && (newParent != swarm.at(hierarchy.at(h).at(j))->getParent())){
						swapNodes(newParent, newH, newWidth, parentNode, parentH, parentWidth, branching, h, width, iterCount);
					}
				}
				catch (const std::exception& e) {
					break;
				}
			}
			iterCount++;
		}
		if ((unsigned int)childCounter==swarm.size())
			break;
	}
}

void Swarm::swapNodes(int newParent, int newH, int newWidth, int parentNode, int parentH, int parentWidth, int branching, int h, int width, int iterCount){
	//Swap child and parent IDs
	swarm.at(hierarchy.at(newH).at(newWidth))->setParent(swarm.at(hierarchy.at(parentH).at(parentWidth))->getParent());
	swarm.at(hierarchy.at(parentH).at(parentWidth))->setParent(newParent); //parent

	//Swap places
	hierarchy.at(parentH).at(parentWidth) = newParent;
	hierarchy.at(newH).at(newWidth) = parentNode;

	if (h <= lastLevelComplete) {
		//Update siblings' parentID
		for (int k=width-(branching-1); k<=width; k++){
			swarm.at(hierarchy.at(h).at(k))->setParent(newParent);
		}
		//Update any children that the old child (now parent) may have had
		for (int k=branching*width; k<(branching*width)+branching; k++){
			try {
				swarm.at(hierarchy.at(h+1).at(k))->setParent(parentNode);
			}
			catch (const std::exception& e) {
				break;
			}
		}
	}
	if (h == lastLevelComplete+1) {
		//Update siblings' parentID
		for(int k=iterCount*branching;k<(iterCount*branching)+branching;k++){
			try{
				swarm.at(hierarchy.at(h).at(k))->setParent(newParent);
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
	cout << "\nvar::lastLevelComplete " << lastLevelComplete << endl << endl;
	for (int i=0; i<=lastLevelComplete+1; i++){
		cout << "Level " << i << ": ";
		for (int j=0; j<pow(branching,i); j++)
			cout << hierarchy.at(i).at(j) << " ";
		cout << endl;
	}

	//Print tree -- this is the same algorithm we use to traverse the tree and update the position of the particles
	cout << "\n\nRoot node: " << endl;
	cout << "\t\t" << hierarchy.at(0).at(0) << "->" << swarm.at(hierarchy.at(0).at(0))->getParent() << endl << "Level 1: " << endl;
	//reset indexes to traverse the entire tree
	int h = 1;
	int width = 0;
	int iterCount = 0;
	int childCounter = 1;
	int parentNode = swarm.at(hierarchy.at(0).at(0))->getParent();
	for(unsigned int i=1;i<swarm.size();i++){
		if (h <= lastLevelComplete){
			if (width < pow(branching,h)){ //max level width
				if (parentNode!=swarm.at(hierarchy.at(h).at(width))->getParent()){
					parentNode=swarm.at(hierarchy.at(h).at(width))->getParent();
					cout << "Parent: " << parentNode << endl;
				}
				cout << "\t\t" << hierarchy.at(h).at(width) << "->" << swarm.at(hierarchy.at(h).at(width))->getParent() << endl;
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
			parentNode = hierarchy.at(h-1).at(iterCount);
			if (parentNode > -1){
				cout << "Parent: " << parentNode << endl;
				for(int j=iterCount*branching;j<(iterCount*branching)+branching;j++){
					try{
						if (swarm.at(hierarchy.at(h).at(j))->getParent() == parentNode) {
							cout << "\t\t" << hierarchy.at(h).at(j) << "->" << swarm.at(hierarchy.at(h).at(j))->getParent() << endl;
							childCounter++;
						}
					}
					catch (const std::exception& e) {
						break;
					}
				}
			}
			iterCount++;
		}
		if ((unsigned)childCounter==swarm.size())
			break;
	}
}

void Swarm::clearResizeSimpSwarm(Configuration* config, long int iteration){
	//Copy particle's id and evaluation in simpSwarm for the next iteration
	simpSwarm.clear();
	simpSwarm.resize(swarm.size());
	for (unsigned int i=0;i<simpSwarm.size();i++){
		simpSwarm.at(i).id = swarm.at(i)->getID();
		simpSwarm.at(i).eval = swarm.at(i)->getCurrentEvaluation();
	}
}

// This function computes the inertia weight (omega1 in the GVU formula) according to the selected strategy
double Swarm::computeOmega1(Configuration* config, long int iteration, long int posIndex, bool newIteration){
	/*TODO: Some of these static variables could become parameters of the framework and be optimized with irace*/
	static double alpha = 1/pow(M_PI,2); 			//small positive constant 		IW_NONL_DEC - 4
	static double omega = 0.3;						//value between [0,1] 			IW_NONL_DEC_IMP - 5
	static double u = 1.0002;						//value between [1.0001,1.0005] IW_NONL_DEC_IMP - 5
	static double omegaXu = omega * u;				//								IW_NONL_DEC_IMP - 5
	static double zFunction;						//								IW_CHAOTIC_DEC - 7
	static int    k = 7;							//positive integer constant		IW_OSCILLATING - 9 This should be a configurable parameter (
	static double simNumOfCos = 2.0*M_PI*((4*k)+6);	//								IW_OSCILLATING - 9
	static double a = 1;							//small positive constant		IW_LOG_DEC - 10
	static double omega_2 = 0;						//								IW_SELF_REGULATING - 11
	static double deltaOmega = (((double)config->getFinalIW()) - config->getInitialIW())/config->getMaxIterations(); //IW_SELF_REGULATING - 11 - Self-regulating
	static double T_0_95 = 95*config->getMaxIterations()/100; //iteration at which 95% of search process is completed IW_VELOCITY_BASED - 12 - Based on velocity information

	double OmegaToReturn = config->getOmega1();

	//If all particle use the same inertia value, it is more efficient
	//to compute it once at the beginning of the iteration
	if (newIteration) {
		/* Non-adaptive strategies */
		//IW_CONSTANT - 0 - Constant value
		if (config->getOmega1CS() == IW_CONSTANT){
			if (config->getOmega1() < -1 || config->getOmega1() > 1) //check convergence bounds
				config->setOmega1(CONSTRICTION_COEFFICIENT);
			OmegaToReturn =config->getOmega1();
		}
		//IW_L_INC - 1 - Linear increasing
		else if (config->getOmega1CS() == IW_L_INC) {
			if (config->getIWSchedule() > 0 ){
				//from Frankenstein's PSO
				if(iteration <= config->getIWSchedule()){
					config->setOmega1(config->getInitialIW() +
							(((double)(config->getIWSchedule() - iteration)/config->getIWSchedule())*
									(config->getFinalIW() - config->getInitialIW()))
					);
				}
				else
					config->setOmega1(config->getFinalIW());
			}
			else{
				config->setOmega1( config->getInitialIW() +
						(((double)(iteration)/config->getMaxIterations())*
								(config->getFinalIW()-config->getInitialIW()))
				);
			}
			OmegaToReturn =config->getOmega1();
		}
		//IW_L_DEC - 2 - Linear decreasing
		else if (config->getOmega1CS() == IW_L_DEC) {
			if (config->getIWSchedule() > 0 ){
				//from Frankenstein's PSO
				if(iteration <= config->getIWSchedule()){
					config->setOmega1( config->getFinalIW() +
							(((double)(config->getIWSchedule() - iteration)/config->getIWSchedule())*
									(config->getInitialIW() - config->getFinalIW()))
					);
				}
				else
					config->setOmega1(config->getFinalIW());
			}
			else {
				config->setOmega1( config->getFinalIW() +
						((config->getInitialIW()-config->getFinalIW())*iteration)/config->getMaxIterations()
				);
			}
			OmegaToReturn =config->getOmega1();
		}
		//IW_RANDOM - 3 - Random
		else if (config->getOmega1CS() == IW_RANDOM) {
			config->setOmega1( 0.5 * (RNG::randVal(0,1)/2.0));
			OmegaToReturn =config->getOmega1();
		}
		//IW_NONL_DEC - 4 - Nonlinear decreasing
		else if (config->getOmega1CS() == IW_NONL_DEC) {
			config->setOmega1(
					config->getFinalIW() + ((config->getInitialIW()-config->getFinalIW())*
							pow((double)(iteration)/config->getMaxIterations(),alpha))
			);
			OmegaToReturn =config->getOmega1();
		}
		//IW_NONL_DEC_IMP - 5 - Nonlinear decreasing improved
		else if (config->getOmega1CS() == IW_NONL_DEC_IMP) {
			config->setOmega1(
					pow(omegaXu, iteration)
			);
			OmegaToReturn =config->getOmega1();
		}
		//IW_NONL_DEC_TIME - 6 - Nonlinear decreasing time-dependent
		else if (config->getOmega1CS() == IW_NONL_DEC_TIME) {
			config->setOmega1(
					pow((2.0/iteration), 0.3)
			);
			OmegaToReturn =config->getOmega1();
		}
		//IW_CHAOTIC_DEC - 7 Chaotic decreasing
		else if (config->getOmega1CS() == IW_CHAOTIC_DEC) {
			iteration == 1 ? zFunction = RNG::randVal(0,1) : zFunction = 4*zFunction*(1-zFunction);
			config->setOmega1(
					(zFunction*config->getInitialIW()) + (config->getFinalIW()-config->getInitialIW()) *
					((double)(config->getMaxIterations())-iteration)/config->getMaxIterations()
			);
			OmegaToReturn =config->getOmega1();
		}
		//IW_EXP_DEC - 8 - Natural exponential decreasing
		else if (config->getOmega1CS() == IW_EXP_DEC) {
			config->setOmega1(
					config->getInitialIW() + (config->getFinalIW()-config->getInitialIW())*
					exp((-10.0*iteration)/config->getMaxIterations())
			);
			OmegaToReturn =config->getOmega1();
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
			OmegaToReturn =config->getOmega1();
		}
		//IW_LOG_DEC - 10 - Logarithm decreasing
		else if (config->getOmega1CS() == IW_LOG_DEC) {
			config->setOmega1(
					config->getFinalIW() + (config->getInitialIW()-config->getFinalIW())*
					log10(((10.0*iteration)/config->getMaxIterations())+ a )
			);
			OmegaToReturn =config->getOmega1();
		}

		/* ****************************************************************************************************************/
		/* ****************************************************************************************************************/

		/* Adaptive strategies */
		//IW_SELF_REGULATING - 11 - Self-regulating
		else if (config->getOmega1CS() == IW_SELF_REGULATING) {
			iteration == 1 ? omega_2 = config->getFinalIW() : omega_2 = omega_2-deltaOmega;
			config->setOmega1(omega_2);
			OmegaToReturn =config->getOmega1();
		}
		//IW_VELOCITY_BASED - 12 - Based on velocity information
		else if (config->getOmega1CS() == IW_VELOCITY_BASED) {
			double idealVelocity;
			double avVel;
			double lambda = config->get_iw_par_deltaOmega();
			if (iteration == 1)
				config->setOmega1(config->getFinalIW());

			idealVelocity = swarm.at(0)->getMaxVelLimit() * ((1.0 + cos(M_PI*(iteration/T_0_95)))/2);
			avVel=computeAvgVelocity(config);	//average absolute velocity of the swarm

			if (avVel >= idealVelocity){
				(config->getOmega1()-lambda) >= config->getInitialIW() ?
						config->setOmega1(config->getOmega1()-lambda) : config->setOmega1(config->getInitialIW());
			}
			else{
				(config->getOmega1()+lambda) >= config->getFinalIW() ?
						config->setOmega1(config->getFinalIW()) : config->setOmega1(config->getOmega1()+lambda);
			}
			OmegaToReturn =config->getOmega1();
		}
		//IW_RANKS_BASED - 14 - Rank-based
		else if (config->getOmega1CS() == IW_RANKS_BASED) {
			rankParticles(simpSwarm);
			OmegaToReturn = CONSTRICTION_COEFFICIENT;
		}
		//IW_SUCCESS_BASED - 15 Success-based
		else if (config->getOmega1CS() == IW_SUCCESS_BASED) {
			if (iteration == 1){
				config->setOmega1(config->getFinalIW()); //set inertia weight to its maximum value for the first iteration
			}
			if (iteration > 1){
				int S_i = 0; //Number of solutions that improved after the last iteration
				for (unsigned int i=0; i<swarm.size(); i++){
					for (unsigned int j=0; j<simpSwarm.size(); j++){
						if (swarm.at(i)->getID() == simpSwarm.at(j).id){
							//evaluate if the solution improved
							if (swarm.at(i)->getCurrentEvaluation() < simpSwarm.at(j).eval){
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
			}
			OmegaToReturn =config->getOmega1();
		}
		//IW_CONVERGE_BASED - 16 Convergence-based
		else if (config->getOmega1CS() == IW_CONVERGE_BASED) {
			if (iteration == 1 ){
				double alpha_2 = config->get_iw_par_alpha_2();
				double beta_2 = config->get_iw_par_beta_2();
				//				if (config->getPopulationCS() == POP_CONSTANT)
				config->setOmega1(1 - fabs(alpha_2/(1+beta_2)));
			}
			OmegaToReturn =config->getOmega1();
		}
		//Keep inertia constant during the execution
		else {
			config->setOmega1(CONSTRICTION_COEFFICIENT);
			OmegaToReturn =config->getOmega1();
		}
	}
	//These are the strategies that need to compute a independent inertia value for each particle
	if (posIndex != -1){
		//IW_L_INC - 1 - Linear increasing
		if (config->getOmega1CS() == IW_L_INC){
			if (config->getTopology() == TOP_HIERARCHICAL){
				//From Hierarchical PSO
				double k=((double)getParticleNumParents(posIndex));
				config->setOmega1( config->getFinalIW() +
						( ((config->getInitialIW()-config->getFinalIW()) * k) / (lastLevelComplete+1))
				);
				OmegaToReturn = config->getOmega1();
			}
			else
				OmegaToReturn = config->getOmega1();
		}
		//IW_L_DEC - 2 - Linear decreasing
		else if (config->getOmega1CS() ==  IW_L_DEC) {
			if (config->getTopology() == TOP_HIERARCHICAL){
				//From Hierarchical PSO
				double k=((double)getParticleNumParents(posIndex));
				config->setOmega1( config->getInitialIW() +
						( ((config->getFinalIW()-config->getInitialIW()) * k)/(lastLevelComplete+1))
				);
				OmegaToReturn = config->getOmega1();
			}
			else
				OmegaToReturn = config->getOmega1();
		}
		//IW_SELF_REGULATING - 11 - Self-regulating
		else if (config->getOmega1CS() ==  IW_SELF_REGULATING){
			double eta = config->get_iw_par_eta();
			//The best particle has a special inertia value
			if (posIndex == best_particle->getID())
				OmegaToReturn = omega_2 + eta * ((config->getFinalIW() - config->getInitialIW()) / config->getMaxIterations());
			else
				OmegaToReturn = config->getOmega1();
		}
		//IW_DOUBLE_EXP - 13 - Double exponential self-adaptive
		else if (config->getOmega1CS() ==  IW_DOUBLE_EXP){
			double R =0.0;
			if (iteration == 1){
				config->setOmega1(config->getFinalIW());
			}
			else {
				R = swarm.at(posIndex)->computeDistPbestGbest()*((((double)config->getMaxIterations())-iteration)/config->getMaxIterations());
				config->setOmega1( exp(-1*exp((R*-1))));
			}
			OmegaToReturn = config->getOmega1();
		}
		//IW_RANKS_BASED - 14 - Rank-based
		else if (config->getOmega1CS() ==  IW_RANKS_BASED){
			config->setOmega1( config->getInitialIW() + ((config->getFinalIW()-config->getInitialIW()) *
					((double) swarm.at(posIndex)->getRanking()/swarm.size()))
			);
			OmegaToReturn = config->getOmega1();
		}
		//IW_CONVERGE_BASED - 16 Convergence-based
		else if (config->getOmega1CS() ==  IW_CONVERGE_BASED){
			if (iteration > 1) {
				double alpha_2 = config->get_iw_par_alpha_2();
				double beta_2 = config->get_iw_par_beta_2();

				//Find the particle in simpSwarm
				for (unsigned int j=0; j<simpSwarm.size(); j++){
					if (simpSwarm.at(j).id == swarm.at(posIndex)->getID()){
						//convergence factor
						long double c_i = fabs(simpSwarm.at(j).eval - swarm.at(posIndex)->getCurrentEvaluation())/
								(simpSwarm.at(j).eval + swarm.at(posIndex)->getCurrentEvaluation());
						//diffusion factor
						long double d_i = fabs(swarm.at(posIndex)->getCurrentEvaluation() - global_best.eval) /
								(swarm.at(posIndex)->getCurrentEvaluation() + global_best.eval);
						//set the inertia weight
						config->setOmega1(1 - fabs(alpha_2*(1-c_i)) / (1+d_i)*(1+beta_2));
						break;
					}
				}
				OmegaToReturn = config->getOmega1();
			}
		}
		else{
			OmegaToReturn = config->getOmega1();
		}
	}
	if ((posIndex != -1) && config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE))
		cout << "\tvar::omega1 = " << fixed << OmegaToReturn;

	return (OmegaToReturn);
}


double Swarm::computeOmega2(Configuration* config){
	double Omega2ToReturn = 0;
	switch (config->getOmega2CS()) {
	case O2_EQUAL_TO_O1:
		Omega2ToReturn = config->getOmega1();
		break;
	case O2_RANDOM:
		Omega2ToReturn = 0.5 + (RNG::randVal(0,1)/2.0);
		break;
	case O2_ZERO:
		Omega2ToReturn =  0; //the component is not used
		break;
	case O2_CONSTANT:
		Omega2ToReturn = config->getOmega2();
		break;
	default:
		Omega2ToReturn =  1; //no strategy, set the value to one
		break;
	}

	if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE))
		cout << "\tvar::omega2 = " << fixed << Omega2ToReturn;

	return (Omega2ToReturn);
}

double Swarm::computeOmega3(Configuration* config){
	double Omega3ToReturn = 0;
	switch (config->getOmega3CS()) {
	case O3_EQUAL_TO_O1:
		Omega3ToReturn = config->getOmega1();
		break;
	case O3_RANDOM:
		Omega3ToReturn = 0.5 + (RNG::randVal(0,1)/2.0);
		break;
	case O3_ZERO:
		Omega3ToReturn = 0;
		break;
	case O3_CONSTANT:
		Omega3ToReturn = config->getOmega3();
		break;
	default:
		Omega3ToReturn = 1;
		break;
	}

	if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE))
		cout << "\tvar::omega3 = " << fixed << Omega3ToReturn;

	return (Omega3ToReturn);
}

double Swarm::computeAvgVelocity(Configuration* config){
	double sumVelocity = 0.0;
	for (unsigned int i=0;i<swarm.size();i++){
		for (int j=0;j<config->getProblemDimension();j++){
			sumVelocity = sumVelocity + fabs(swarm.at(i)->getCurrentVelocity()[j]);
		}
	}
	return (sumVelocity/(swarm.size()*config->getProblemDimension()));
}

void Swarm::rankParticles(vector<SimplifySwarm> &sSwarm){
	//The ranking are obtained via the mergesort algorithm implemented in Utils.cpp

	ssArray sS;
	sS.id = new int [sSwarm.size()];
	sS.eval = new long double [sSwarm.size()];

	//Copy particle's id and evaluation in simpSwarm
	for (unsigned int i=0;i<sSwarm.size();i++){
		sS.id[i] = swarm.at(i)->getID();
		sS.eval[i] = swarm.at(i)->getCurrentEvaluation();
	}
	//Sort simpSwarm by the value of particles' evaluation
	mergeSort(&sS, 0, sSwarm.size()-1); //mergeSort(array, left(LOWER) index, right (UPPER) index);


	//Set particles' rank in the swarm
	for (unsigned int i=0;i<sSwarm.size();i++){
		for (unsigned int j=0;j<swarm.size();j++)
			if (swarm.at(j)->getID() == sS.id[i]){
				swarm.at(i)->setRanking(i+1);
				sSwarm.at(i).id = sS.id[i];
				sSwarm.at(i).eval = swarm.at(j)->getCurrentEvaluation();
			}
	}
	delete [] sS.id;
	delete [] sS.eval;
}

//mergeSort(arr, 0, arr_size-1);
void Swarm::mergeSort(ssArray* arr, int l, int r){
	if (l < r) {
		// Middle point of the array
		int m = l+(r-l)/2;

		// Sort first and second halves
		mergeSort(arr, l, m);
		mergeSort(arr, m+1, r);
		merge(arr, l, m, r);
	}
}

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void Swarm::merge(ssArray* arr, int l, int m, int r){
	int i;
	int j;
	int k;
	int n1 = m - l + 1;
	int n2 =  r - m;

	/* create temp arrays */
	//struct SimplifySwarm L[n1], R[n2];
	//struct SimplifySwarm L, R;
	long double L_eval [n1];
	int L_id [n1];
	long double R_eval [n2];
	int R_id [n2];

	/* Copy data to temp arrays L[] and R[] */
	for (i = 0; i < n1; i++){
		L_eval[i] = arr->eval[l + i];
		L_id[i] = arr->id[l + i];
	}
	for (j = 0; j < n2; j++){
		R_eval[j] = arr->eval[m + 1+ j];
		R_id[j] = arr->id[m + 1+ j];
	}

	/* Merge the temp arrays back into arr[l..r]*/
	i = 0; // Initial index of first subarray
	j = 0; // Initial index of second subarray
	k = l; // Initial index of merged subarray
	while (i < n1 && j < n2) {
		if (L_eval[i] <= R_eval[j]) {
			arr->eval[k] = L_eval[i];
			arr->id[k] = L_id[i];
			i++;
		}
		else {
			arr->eval[k] = R_eval[j];
			arr->id[k] = R_id[j];
			j++;
		}
		k++;
	}

	/* Copy the remaining elements of L[], if there
       are any */
	while (i < n1) {
		arr->eval[k] = L_eval[i];
		arr->id[k] = L_id[i];
		i++;
		k++;
	}

	/* Copy the remaining elements of R[], if there
       are any */
	while (j < n2) {
		arr->eval[k] = R_eval[j];
		arr->id[k] = R_id[j];
		j++;
		k++;
	}
}
