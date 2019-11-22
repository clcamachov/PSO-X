/*
 * particle.h
 *
 *  Created on: May 31, 2018
 *      Author: leonardo
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <vector>
#include "problem.h"
#include "config.h"
#include "utils.h"

/* Solution structure */
struct Solution {
	double *x;    /* coefficients of each dimension */
	long double eval;  /* value of the solution */
};

using namespace std;

class Particle {

protected:

	Problem* problem;		/* the optimization problem */
	long int size; 			/* problem size: number of dimensions of the functions */

	bool init;
	int  id;
	bool hasVelocitybounds;
	int ranking;
	int parent;
	int stereotype;
	double perturbationVal;

	/*Solution variables*/
	//Each particle has to remember this three vectors at each iteration
	struct Solution current;  /* current solution */
	struct Solution pbest;    /* personal best solution */
	struct Solution gbest;    /* global best solution (According to the topology) */

	/*Velocity variables*/
	double* velocity;               /* velocity */
	double phi_1;
	double phi_2;
	double inertia; 				/* parameters */

	double minVelLimit;
	double maxVelLimit;

public:

	vector < Particle* > neighbours;  /* vector of neighbors particles */
	Particle ();  														/* empty constructor */
	~Particle();  														/* destructor */
	Particle (Problem* problem, Configuration* config, int identifier);	/* constructor */
	Particle (const Particle &p);  										/* copy constructor */
	Particle& operator= (const Particle& p);  							/* overriding of '=' */

	void move(Configuration* config, double minBound, double maxBound, long int iteration,
			double omega1, double omega2, double omega3, int numInformants, int *theInformants, int lastLevelComplete, double alpha_t, double l, double delta);
	double computeNewVelocity(Configuration* config, double vel, double rand1, double rand2,double perInf, double socInf, double pos, double additionalVal);

	double* getCurrentPosition();
	long double getCurrentEvaluation();
	double* getPbestPosition();
	long double getPbestEvaluation();
	double* getCurrentVelocity();
	void computeEvaluation();

	void addNeighbour(Particle* p);
	void checkNeibourhood();
	int getBestOfNeibourhood();

	void updateGlobalBest(double* x, double eval);
	unsigned int getNeighborhoodSize();

	double computeDistPbestGbest();
	double computeDistance(double * x, double * p);

	//Velocity
	void setVelocityLimits(Configuration* config);
	double getMinVelLimit();
	double getMaxVelLimit();

	void initializeUniform();
	void printPosition();
	void printNeighborByID(int id);

	void getHypersphericalVector(double* H, double* V1);
	int getRandomNeighbor();

	//Perturbation
	double computePerturbation(Configuration* config, double * pos_x, double * pbest_x, double alpha_t,
			double l, double delta, bool newIteration);

	//Frankenstein's members
	int getID();
	void eraseNeighborbyID(int nid);
	int getRandomNonAdjacentNeighborID(Configuration* config);

	//inertia
	int getRanking();
	void setRanking(int rank);

	//hierarchical topology
	int getParent();
	void setParent(int papa);

};

#endif /* PARTICLE_H_ */
