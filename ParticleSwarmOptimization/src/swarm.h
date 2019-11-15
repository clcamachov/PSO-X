/*
 * swarm.h
 *
 *  Created on: Jun 20, 2018
 *      Author: leonardo
 */

#ifndef SWARM_H_
#define SWARM_H_

#include <vector>
#include "problem.h"
#include "config.h"
#include "particle.h"
#include "utils.h"

#include <math.h>

class Swarm {

protected:

	Problem* problem;				/* the optimization problem */
	vector<Particle*> swarm;		/* the swarm */
	long int size;					/* number of particles*/

	struct Solution global_best;    /* best solution in the swarm*/
	Particle *best_particle;        /* best particle */
	bool init;
	bool ranked;					/* implement particle rankings */
	bool hierarchical;				/* hierarchical topology flag */

public:
	Swarm ();
	~Swarm();
	Swarm (Problem* problem, Configuration* config); /*Create swarm structure*/
	Swarm (const Swarm &s, Configuration* config);

	Swarm& operator= (const Swarm& s);

	void updateGlobalBest(double* new_x, double eval);
	void moveSwarm(Configuration* config, long int iteration, const double minBound, const double maxBound);
	void getInformants(Configuration* config, int particleID, long int iteration);

	// Inertia
	double computeOmega1(Configuration* config, long int iteration, long int id, bool newIteration);
	double computeOmega2(Configuration* config);
	double computeOmega3(Configuration* config);
	double computeAvgVelocity(Configuration* config);
	void rankParticles(SimplifySwarm* particles);

	// GlobalBest value to print
	Solution getGlobalBest();
	void printGbest(unsigned int dimensions);

	/* Available topologies*/
	void createFullyConnectedTopology();
	void createRingTopology();
	void createWheelTopology();
	void createRandomTopology();
	void createVonNeumannTopology();
	void updateTimeVaryingTopology(Configuration* config, long int iterations);

	//Hierarchical topology
	void createHierarchical(int branching);
	void printTree(int branching);
	void swapNodes(int newParent, int newH, int newWidth, int parentNode, int parentH, int parentWidth, int branching, int h, int width, int iterCount);
	void updateTree(int branching);
	bool isHierarchical();
	void getParticleParentsIDs(int particleID, int *ParentsArray1D);
	void printAllParentNodes();
};

#endif /* SWARM_H_ */
