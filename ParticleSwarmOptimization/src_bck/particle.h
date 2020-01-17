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
	double perturbMagnitud1;
	double perturbMagnitud2;
	int gBestID;

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
	int getBestOfNeibourhood();
	int getgBestID();
	void setgBestID(int gB_ID);

	void updateGlobalBest(double* x, double eval);
	unsigned int getNeighborhoodSize();

	double computeDistPbestGbest();
	double computeDistance(double * x, double * p);

	bool ispBestIntheInformants(int numInformants, int *theInformants);

	//Velocity
	void setVelocityLimits(Configuration* config);
	double getMinVelLimit();
	double getMaxVelLimit();

	void initializeUniform();
	void printPosition();
	void printNeighborByID(int id);

	double ** computeSubtractionPerturbationRotation(int DNPP, double ** vect_PbestMinusPosition, int &numInformants, int *theInformants, double *** rndMatrix, int RmatrixType,
				bool pBestIntheInformants, int pertubType, double alpha_t, double l_value);
	void getRectangularDNPP(double * vect_distribution, int numInformants, int *theInformants, bool pBestIntheInformants, double ** vect_PbestMinusPosition, int modelOfInflu);
	void getSphericalDNPP(double * vect_distribution, int numInformants, int *theInformants, bool pBestIntheInformants, double ** vect_PbestMinusPosition, int modelOfInflu);
	void getAdditiveStochasticDNPP(double * vect_distribution, int numInformants, int *theInformants, bool pBestIntheInformants, double ** vect_PbestMinusPosition, int modelOfInflu, bool randNeighbor, int operatorQ);
	int getRandomNeighbor();
	int getRandomInformantPosition(int numInformants, int *theInformants, bool pBestIntheInformants);
	int getPositionOfpBest(int numInformants, int *theInformants, bool pBestIntheInformants);


	//Perturbation
	void setPerturbationMagnitude(int pertubType, double * pos_x, double * pbest_x, double alpha_t, double l);
	void setPerturbationMagnitude(int pertubType, double alpha_t);
	double getPerturbationMagnitude(int pertubType, double alpha_t, double delta);
	double applyPerturbation(int pertubType, double pos_xi);

	//Random Matrix
	void computeRndMatrix(double *** rndMatrix, int RmatrixType);
	double * multiplyVectorByRndMatrix(double * aVector, double *** rndMatrix, int RmatrixType);

	//Frankenstein's members
	int getID();
	void eraseNeighborbyID(int nid);
	int getRandomNonAdjacentNeighborID(Configuration* config);

	//inertia
	int getRanking();
	void setRanking(int rank);

	//hierarchical topology
	int getParent();
	void setParent(int papa); //some Spanish :P

	//acceleration coefficients
	void setPhi1(double new_phi_1);
	void setPhi2(double new_phi_2);
	double getPhi1();
	double getPhi2();

};

#endif /* PARTICLE_H_ */