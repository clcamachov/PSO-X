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
	bool ranked;					/* implement particle rankings */

public:
	Swarm ();
	~Swarm();
	Swarm (Problem* problem, Configuration* config); /*Create swarm structure*/

	void updateGlobalBest(double* new_x, double eval);
	void moveSwarm(Configuration* config, long int iteration, const double minBound, const double maxBound);
	//void (*setNeighborhood)();		// pointer to neighborhood function

	// Inertia
	void computeInertia(Configuration* config, long int iteration);
	double computeAvgVelocity(Configuration* config);
	void rankParticles(SimplifySwarm* particles);

	// GlobalBest value to print
	Solution getGlobalBest();
	void printGbest(unsigned int dimensions);

	/* Available topologies*/
	void createFullyConnectedTopology();
	void createRingTopology();
	void createStarTopology();
	void createRandomTopology();
	void createVonNeumannTopology();
	void updateTimeVaryingTopology(Configuration* config, long int iterations);
};

#endif /* SWARM_H_ */
