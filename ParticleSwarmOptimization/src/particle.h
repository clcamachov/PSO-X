/*
 * particle.h
 *
 *  Created on: May 31, 2019
 *      Author: christian
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <vector>
#include "problem.h"
#include "config.h"
#include "utils.h"

/* Solution structure */
struct Solution {
	double* x;    /* coefficients of each dimension */
	long double eval;  /* value of the solution */
};

using namespace std;

class Particle {

protected:

	Problem* problem;	/* the optimization problem */

	int size; 			/* problem size: number of dimensions of the functions */
	int id;
	int ranking;
	int parent;
	int stereotype;
	int lbestID;

	/*Solution variables*/
	//Each particle has to remember this three vectors at each iteration
	struct Solution current;  /* current solution */
	struct Solution pbest;    /* personal best solution */
	struct Solution lbest;    /* global best solution (According to the topology) */

	/*Velocity and acceleration coefficients variables*/
	double* velocity;
	double minVelLimit;
	double maxVelLimit;
	double inertia;
	double phi_1;
	double phi_2;

	bool init;

public:

	vector < Particle* > neighbors;  /* vector of neighbors particles */
	vector < int > informants;
	Particle ();  														/* empty constructor */
	~Particle();  														/* destructor */
	Particle (Problem* problem, Configuration* config, int identifier, long int iteration);	/* constructor */
	Particle (const Particle &p);  										/* copy constructor */
	Particle& operator= (const Particle& p);  							/* overriding of '=' */

	void move(Configuration* config, long int iteration, double omega1, double omega2, double omega3,
			int lastLevelComplete, int solImproved);

	double* getCurrentPosition();
	long double getCurrentEvaluation();
	double* getPbestPosition();
	long double getPbestEvaluation();
	double* getCurrentVelocity();
	void setRandomPositionInBoundsWithProbability(Configuration* config);
	void setRandomPositiongBestDerivative(Configuration* config, double* global_bestPos);
	void evaluateSolution();

	void addNeighbour(Particle* p);
	int getBestOfNeibourhood();
	int getlBestID();
	void setlBestID(int gB_ID);

	void updatelBestParticle(double* x, double eval);
	unsigned int getNeighborhoodSize();
	double computeDistPbestGbest();

	bool ispBestIntheInformants(int numInformants);

	//Velocity
	void setVelocityLimits(Configuration* config);
	double getMinVelLimit();
	double getMaxVelLimit();

	//Initial position
	void initUniform(Configuration* config);
	void initToModel();
	void initializePosition(Configuration* config, long int initialPopSize, bool updatePbest);
	void initializeVelocity(Configuration* config);
	void printPosition();
	void printNeighborByID(int id);

	//Velocity and position computation
	void computeSubtractionPerturbationRotation(Configuration* config, vector< vector<double> > &vect_PbestMinusPosition,
			long int iteration, int solImprov);
	void getRectangularDNPP(Configuration* config, double vect_distribution[], vector< vector< double> > &vect_PbestMinusPosition);
	void getSphericalDNPP(Configuration* config, double vect_distribution[], vector< vector< double> > &vect_PbestMinusPosition);
	void getAdditiveStochasticDNPP(Configuration* config, double vect_distribution[], vector< vector< double> > &vect_PbestMinusPosition);

	void computeAC(Configuration* config, double &c1, double &c2);
	int getRandomInformantPosition();
	int getPositionOfpBest();
	void detectStagnation(Configuration* config);

	//Perturbation
	void setPerturbation1Magnitude(Configuration* config, double pertMagnitude[], double * pos_x, double * pbest_x);
	double applyInformedPerturbation(Configuration* config, double pertMagnitude[], double pos_xi, int index, long int iteration);
	void getRandomAdditivePerturbation(Configuration* config, double vect_perturbation[]);

	//Random Matrix
	void computeRndMatrix(Configuration* config, double ** rndMatrix, int RmatrixType, double angle);
	void multiplyVectorByRndMatrix(Configuration* config, vector<vector< double> > &vect_PbestMinusPosition, int informant, double ** rndMatrix,
			int RmatrixType, int solImprov, long int iteration);
	double getAnAngle(Configuration* config, int solImprov, long int iteration);

	//Frankenstein's members
	int getID();
	void setID(int newID);
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
