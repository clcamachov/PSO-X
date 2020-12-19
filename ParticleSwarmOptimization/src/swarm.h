/*
 * swarm.h
 *
 *  Created on: May 31, 2019
 *      Author: christian
 */

#ifndef SWARM_H_
#define SWARM_H_

#include <vector>
#include "problem.h"
#include "config.h"
#include "particle.h"
#include "utils.h"

#include <math.h>

/* Simplified swarm structure */
struct ssArray {
	long double* eval;   /* value of the solution */
	int* id;  		     /* particle id */
};

class Swarm {

protected:

	Problem* problem;				  /* the optimization problem */

	vector<Particle*> swarm;		  /* the swarm */
	vector<SimplifySwarm> simpSwarm;  /* a vector containing only ids and eval values*/
	vector<vector< int > > hierarchy; /* the hierarchical topology is implemented as a vector*/

	struct Solution global_best;      /* best solution in the swarm*/
	Particle *best_particle;          /* best particle */

	long int size;					  /* number of particles*/
	bool init;
	int lastLevelComplete;			  /* global variable for the hierarchy */

public:

	Swarm ();
	~Swarm();
	Swarm (Problem* problem, Configuration* config); /*Create swarm structure*/

	void updateGlobalBest(double* new_x, double eval);
	void moveSwarm(Configuration* config, const long int iteration);
	void getInformants(Configuration* config, int particleID, long int iteration);	//returns the number of informants of a particle

	// Inertia and acceleration coefficients
	void computeAccelerationCoefficients(Configuration* config, long int iteration);
	void decomposePhi2(int modelOfInflu, int part_id, int numInformants);
	double computeOmega1(Configuration* config, long int iteration, long int id, bool newIteration);
	double computeOmega2(Configuration* config);
	double computeOmega3(Configuration* config);
	double computeAvgVelocity(Configuration* config);

	// GlobalBest particle
	Solution getGlobalBest();
	void printGbest(unsigned int dimensions);

	/* Available topologies*/
	void createFullyConnectedTopology();
	void createRingTopology();
	void createWheelTopology();
	void createRandomEdge();
	void createVonNeumannTopology();
	void updateTimeVaryingTopology(Configuration* config, long int iterations);

	//Hierarchical topology
	void createHierarchical(int branching, int finalPopSize);
	void printTree(int branching);
	void swapNodes(int newParent, int newH, int newWidth, int parentNode, int parentH, int parentWidth, int branching, int h, int width, int iterCount);
	void updateTree(int branching);
	bool isHierarchical();
	void getParticleParentsIDs(int particleID, int *ParentsArray1D);
	void printAllParentNodes();
	int getParticleNumParents(int particleID);
	void updateIDinTree(int branching, long swarm_size, int OldAndNewIDs[][1]);

	//Perturbation
	void updatePerturbationVariables(Configuration* config, double previousGbest_eval, double currentGbest_eval, long int iteration);
	int countImprovedSolutions(Configuration* config, long int iteration);

	//Population size
	void resizeSwarm(Problem* problem, Configuration* config, long int iteration);
	void addParticles(Problem* problem, Configuration* config, int numOfParticles, long int iteration);
	void updateTopologyConnections(Configuration* config, long previous_size, long int iteration);
	void addParticlesInLastLevel(int first, int last, int branching);
	void updateHierarchical(int branching, long previous_size);
	void clearResizeSimpSwarm(Configuration* config, long int iteration);
	void removeOneParticle(Configuration* config);

	//Utils
	void rankParticles(vector<SimplifySwarm> &sSwarm);
	void mergeSort(ssArray* arr, int l, int r);
	void merge(ssArray* arr, int l, int m, int r);
	void reinitializeParticlePosition(Configuration* config, long int iteration);
};

#endif /* SWARM_H_ */
