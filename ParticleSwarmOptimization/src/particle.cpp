/*
 * particle.cpp
 *
 *  Created on: May 31, 2018
 *      Author: leonardo
 */

#include "float.h"
#include "particle.h"
#include "iostream"
#include "rng.h"
#include <math.h>
#include "utils.h"

using namespace std;

double additionalVal= 0.0; //variable in case some additional value has to be use to update the formula

/* Default constructor*/
Particle::Particle(){
	id = 0;
	phi_1 = 0;
	phi_2 = 0;
	inertia = 0;
	//velocity = {0};
	velocity = 0;
	size = 0;
	problem = 0;
	init=false;
	hasVelocitybounds = false;
	minVelLimit = -100;					//minimum velocity -v_max = (x_max-x_min)/2
	maxVelLimit = 100;					//maximum velocity  v_max
	ranking = 0;
	parent =0;
	stereotype = 0;
	perturbMagnitud1 = 0;
	perturbMagnitud2 = 0;
	gBestID = -1;
}

/* Constructor*/
Particle::Particle (Problem* problem, Configuration* config, int identifier){

	this->problem = problem;
	size= problem->getProblemDimension();
	id = identifier;
	ranking = 0;
	parent = 0;
	stereotype = 0;
	perturbMagnitud1 = 0.01;
	perturbMagnitud2 = 0.01;
	gBestID = -1;

	hasVelocitybounds = config->useVelocityClamping();
	setVelocityLimits(config); //Set velocity clamping limits

	current.x = new double[size];
	pbest.x = new double[size];
	gbest.x = new double[size];
	current.eval = LDBL_MAX;
	pbest.eval = LDBL_MAX;
	gbest.eval = LDBL_MAX;

	velocity = new double[size];

	phi_1 = config->getPhi1();
	phi_2 = config->getPhi1();
	inertia = config->getInertia();

	//cout << "Particle " << id << " : [ ";
	initializeUniform();

	//printPosition();
	init = true;
}

/* Copy constructor */
Particle::Particle (const Particle &p){
	problem = p.problem;
	size = p.size;
	id = p.id;
	ranking = p.ranking;
	parent = p.parent;
	stereotype = p.stereotype;
	perturbMagnitud1 = p.perturbMagnitud1;
	perturbMagnitud2 = p.perturbMagnitud2;
	gBestID = p.gBestID;

	hasVelocitybounds = p.hasVelocitybounds;
	minVelLimit = p.minVelLimit;
	maxVelLimit = p.maxVelLimit;

	if(!init){
		current.x = new double[size];
		pbest.x = new double[size];
		gbest.x = new double[size];
		velocity = new double[size];
	}

	for(int i=0; i<size;i++){
		current.x[i] = p.current.x[i];
		pbest.x[i] = p.pbest.x[i];
		gbest.x[i] = p.gbest.x[i];
		velocity[i] = p.velocity[i];
	}

	current.eval = p.current.eval;
	pbest.eval = p.pbest.eval;
	gbest.eval = p.gbest.eval;

	inertia=p.inertia;
	phi_1 = p.phi_1;
	phi_2 = p.phi_2;
	init = true;

}

/* destructor */
Particle::~Particle(){
	if(init){
		delete[] current.x;
		delete[] pbest.x;
		delete[] gbest.x;
		delete[] velocity;
	}
	init=false;
}

/* overriding of '=' operator for particles
   (now 'p1 = p2;' does what one would expect) */
Particle& Particle::operator= (const Particle& p){
	if (this != &p) {
		problem = p.problem;
		size = p.size;
		id = p.id;
		ranking = p.ranking;
		parent = p.parent;
		stereotype = p.stereotype;
		perturbMagnitud1 = p.perturbMagnitud1;
		perturbMagnitud2 = p.perturbMagnitud2;
		gBestID = p.gBestID;

		hasVelocitybounds = p.hasVelocitybounds;
		minVelLimit = p.minVelLimit;
		maxVelLimit = p.maxVelLimit;

		if(!init){
			current.x = new double[size];
			pbest.x = new double[size];
			gbest.x = new double[size];
			velocity = new double[size];
		}

		for(int i=0; i<size;i++){
			current.x[i] = p.current.x[i];
			pbest.x[i] = p.pbest.x[i];
			gbest.x[i] = p.gbest.x[i];
			velocity[i] = p.velocity[i];
		}

		current.eval = p.current.eval;
		pbest.eval = p.pbest.eval;
		gbest.eval = p.gbest.eval;

		inertia=p.inertia;
		phi_1 = p.phi_1;
		phi_2 = p.phi_2;
		init = true;
	}
	return *this;
}

/* Initialize particle */
void Particle::initializeUniform(){
	for(int i=0; i<size; i++){
		current.x[i] = problem->getRandomX(); //random values within the bounds of the function
		pbest.x[i]=current.x[i];
		velocity[i]=0;
	}
	computeEvaluation();
}

void Particle::setVelocityLimits(Configuration* config){
	if((config->getProblemID() == 6 || config->getProblemID() == 24) && (config->getCompetitionID() == 0)){
		minVelLimit= (LDBL_MAX)/-2;
		maxVelLimit= (LDBL_MAX)/2;
	}
	else {
		minVelLimit=((config->getMaxInitRange()-config->getMinInitRange())/2.0)*-1;
		maxVelLimit=((config->getMaxInitRange()-config->getMinInitRange())/2.0);
	}
}

double Particle::getMinVelLimit(){
	return minVelLimit;
}

double Particle::getMaxVelLimit(){
	return maxVelLimit;
}

int Particle::getRanking(){
	return ranking;
}

void Particle::setRanking(int rank){
	ranking = rank;
}

int Particle::getParent(){
	return parent;
}

void Particle::setParent(int node){
	parent = node;
}

bool Particle::isPbestInformant(int numInformants, int *theInformants){
	bool result = false;
	for (int j=0; j<numInformants; j++)
		if (this->id == neighbours.at(theInformants[j])->getID() )
			result = true;

	return result;
}

/* Generate a new solution by updating the particle's position */
void Particle::move(Configuration* config, double minBound, double maxBound, long int iteration,
		double omega1, double omega2, double omega3, int numInformants, int *theInformants, int lastLevelComplete,
		double alpha_t, double l, double delta){

	double V2[size], V1[size]; 				 //used for DIST_SPHERICAL
	//double * vect_PosToPbest[numInformants]; //used for DIST_RECTANGULAR and DIST_MULTISPHERICAL

	/**
	 * ROTATION MATRICES	--->	has to be computed per Informant
	 * **/
	//number of rotation matrices to rotate in all possible planes
	int numMatrices =(int)floor(((size*(size-1))/2));
	double **rndMatrix[numMatrices];
	//Compute the random matrix if its going to be used
	for(int i=0; i<numMatrices; i++){
		rndMatrix[i] = new double*[size];
		for (unsigned int j=0; j<size; j++)
			rndMatrix[i][j] = new double[size];
	}
	computeRndMatrix(rndMatrix, config->getRandomMatrix()); //we need to compute this for each informant

	/**
	 * PERTURBATION 1	--->	has to be computed per Informant and/or per Dimension
	 * **/
	//For the distribution-based perturbation strategies we need to compute the std. deviation in advance,
	//If the last flag is true, we are computing the magnitude of the perturbation for a new informant
	perturbMagnitud1 = computePerturbation(config->getPerturbation1(), current.x, neighbours[theInformants[0]]->pbest.x, alpha_t, l, delta, true);

	/**
	 * PERTURBATION 2	--->	has to be computed per dimension
	 * **/
	perturbMagnitud2 = computePerturbation(config->getPerturbation2(), alpha_t, delta); //Just additive

	switch (config->getDistributionNPP()) {
	case DIST_RECTANGULAR:
		for (int h=0; h<numInformants; h++){

		}
		break;
	case DIST_SPHERICAL:
		getHypersphericalVector(V2, V1, numInformants, theInformants, rndMatrix, config->getRandomMatrix());
		//additionalVal = V2[i] + V1[i];
		break;
	case DIST_MULTISPHERICAL: // has to be computed per Informant
		break;
	case DIST_ADD_STOCH:
		break;
	}

	//Compute new position
	for (int i=0;i<size;i++) {

		//computePerturbation1 can be directly multiply or added (depending on the strategy) to the vector dimension
		perturbMagnitud1 = computePerturbation(config->getPerturbation1(), current.x, neighbours[0]->current.x, alpha_t, l, delta, false);
		//cout << "\n Normal_perturbed: " << RNG::randGaussWithMean(pow(perturbMagnitud1,2), current.x[i]) << endl;

		double PersonalInfluence = pbest.x[i]-current.x[i];
		double SocialInfluence = gbest.x[i]-current.x[i];

		velocity[i] = computeNewVelocity(config, velocity[i],
				1.0, 1.0,
				PersonalInfluence,
				SocialInfluence,
				current.x[i],
				additionalVal);

		//Check velocity clamping bounds
		if (hasVelocitybounds) {
			if (velocity[i] > maxVelLimit)
				current.x[i] = current.x[i] + maxVelLimit;
			else if (velocity[i] < minVelLimit)
				current.x[i] = current.x[i] + minVelLimit;
			else
				current.x[i] = current.x[i] + velocity[i];
		}
		else
			current.x[i] = current.x[i] + velocity[i];

		//values assigned to the variables cannot be outside the function bounds
		if(current.x[i] < minBound)
			current.x[i] = minBound;
		if(current.x[i] > maxBound)
			current.x[i]= maxBound;
	}
	computeEvaluation(); //evaluate the objective function and update pbest if a new pbest has been found

	for(int i=0; i<numMatrices; i++){
		for (unsigned int j=0; j<size; j++)
			delete [] rndMatrix[i][j];
		delete [] rndMatrix[i];
	}
}

double Particle::computeNewVelocity(Configuration* config, double vel, double u1, double u2, double perInf,
		double socInf, double pos, double additionalVal){
	double new_vel = 0.0;

	//Configurable PSO velocity equation
	//Constriction coefficient PSO
	//if (config->getVelocityRule() == VEL_CONSTRICTED){
	new_vel = CONSTRICTION_COEFFICIENT * (
			(inertia * vel) +
			(phi_1 * problem->getRandom01() * perInf) +
			(phi_2 * problem->getRandom01() * socInf)
	);
	//}
	return new_vel;
}

// The computation of the radius and the random point in the HyperSphere
// was taken from the publicly available code of Maurice Clerc - Standard PSO 2011
// https://www.particleswarm.info/Programs.html
void Particle::getHypersphericalVector(double V2[], double V1[], int numInformants, int *theInformants, double *** rndMatrix, int RmatrixType){
	double G[size];	//center of the sphere
	double l[size]; //particle's Gbest
	double radius = 0.0;	//radius G-X
	double pw=1./(double)size;
	double * v_PosToInfPbest[numInformants];
	bool pBestInformant = isPbestInformant(numInformants, theInformants); //flag to know if pBest is included in theInformants

	//When the hierarchical topology is used, this verification might be irrelevant most of the time
	//since the gBest particle can be located in a lower level in the hierarchy.
	//1.- Check if the particle is gBest
	if (this->id == this->gBestID){
		//use a random neighbor as Gbest
		int randNeighbor = getRandomNeighbor();
		if (randNeighbor != this->id){
			cout << "\t\t ... using a random neighbor's pbest--" << endl;
			for (int i=0; i<size; i++)
				l[i] = neighbours.at(randNeighbor)->pbest.x[i]; //Use random neighbor's pBest
		}
		//initialization rule of Incremental PSO -- reinitialize position to model
		else{
			cout << "\t\t ... using reinitialization to model --" << endl;
			for (int i=0; i<size; i++)
				l[i] = problem->getRandomX() + problem->getRandom01()*(gbest.x[i]-current.x[i]);
		}
	}
	else
		for (int i=0; i<size; i++)
			l[i] = gbest.x[i];

	//2.- Rotate if some type of rotation is chosen
	//If gBest is not included in theInformant array, it will be added at the end
	if (pBestInformant == false)
		numInformants = numInformants+1;
	//Create a temporary structure
	for (int i=0; i<numInformants; i++)
		v_PosToInfPbest[i] = new double [size];

	//Get the p^k-x^i of all Informants and then add pBest as the end of the Array
	if (pBestInformant == false){
		//Copy the informants in a temporary structure
		for (int j=0; j<numInformants-1; j++){
			if (this->id == this->gBestID){ //use the right gBest
				for (int i=0; i<size; i++)
					v_PosToInfPbest[j][i] = (l[i]-current.x[i]); //vector to be rotated
			}
			else
				for (int i=0; i<size; i++)
					v_PosToInfPbest[j][i] = (neighbours.at(theInformants[j])->pbest.x[i]-current.x[i]); //vector to be rotated
		}
		for (int i=0; i<size; i++)
			v_PosToInfPbest[numInformants-1][i] = (pbest.x[i]-current.x[i]); //vector to be rotated
	}
	else {
		//Copy the informants in a temporary structure
		for (int j=0; j<numInformants; j++){
			if (this->id == this->gBestID){ //use the right gBest
				for (int i=0; i<size; i++)
					v_PosToInfPbest[j][i] = (l[i]-current.x[i]);  //vector to be rotated
			}
			else
				for (int i=0; i<size; i++)
					v_PosToInfPbest[j][i] = (neighbours.at(theInformants[j])->pbest.x[i]-current.x[i]);  //vector to be rotated
		}
	}
	//Rotate each vector
	for (int j=0; j<numInformants; j++){ //the rest of informants
		multiplyVectorByRndMatrix(v_PosToInfPbest[j], rndMatrix, RmatrixType);
		if (j<numInformants) //avoid computing a new random matrix in the last iteration
			computeRndMatrix(rndMatrix, RmatrixType); //compute a random matrix for the next informant
	}

	//Compute G (center of the sphere) and V1 (radius of each dimension)
	for (int i=0; i<size; i++){
		double R = 0.0;
		for (int j=0; j<numInformants; j++){
			//cout << "\t\t  informant at [" << j << "] is " <<  theInformants[j] << endl; //remove
			R += current.x[i] + (phi_1 * v_PosToInfPbest[j][i]); //v_PosToInfPbest is the vector already rotated
		}
		G[i] = (current.x[i] + R )/(numInformants + 1);
		radius += pow(abs(current.x[i] - G[i]), 2);
		V1[i] = G[i] - current.x[i];
	}
	radius = sqrt(radius); //this is the actual radius of the hyper-sphere

	// Get a random vector in the hyper-sphere H(G||G-X||)
	// ----------------------------------- Step 1.  Direction
	double length=0.0;
	for (int i=0;i<size;i++) {
		V2[i] = RNG::randGauss(1.0);
		length += pow(V2[i],2);
	}
	length=sqrt(length);
	//----------------------------------- Step 2. Random radius
	// Random uniform distribution on the sphere
	double r = pow(RNG::randVal(0.0,1.0),pw);

	//----------------------------------- Step 3. Random vector
	for (int i=0;i<size;i++) {
		V2[i] = radius*r*V2[i]/length;
	}

	//if (pBestInformant == false && RmatrixType != MATRIX_NONE)
	if (pBestInformant == false)
		numInformants-=1;
}

//Perturbation 1
double Particle::computePerturbation(int pertubType, double * pos_x, double * pbest_x, double alpha_t, double l, double delta, bool newInformant){
	switch(pertubType){
	case PERT1_NONE: //Do not apply perturbation
		return 1.0;
	case PERT1_ADD_RECT || PERT1_DIST_SUCCESS: //Additional rectangular
	return alpha_t*(1-(2*problem->getRandomX(0,1)));
	case PERT1_ADD_NOISY : //Additional noisy
		return problem->getRandomX(-delta/2,delta/2);
	case PERT1_DIST_NORMAL: //Normally distributed (here, we only compute the std. deviation
		if (newInformant){
			double distance = computeDistance(pos_x, pbest_x);
			if (distance == 0)
				return (this->perturbMagnitud1);
			else
				return (l*distance);
		}
		else
			return (this->perturbMagnitud1);
	default:
		return 0.001;
	}
}

//Perturbation 2, a.k.a. additional perturbation
double Particle::computePerturbation(int pertubType, double alpha_t, double delta){
	switch(pertubType){
	case PERT2_NONE: //Do not apply perturbation
		return 1.0;
	case PERT2_ADD_RECT : //Additional rectangular
		return alpha_t*(1-(2*problem->getRandomX(0,1)));
	case PERT2_ADD_NOISY: //Additional noisy
		return problem->getRandomX(-delta/2,delta/2);
	default:
		return 0.001;
	}
}

//Multiply by the random matrix
double * Particle::multiplyVectorByRndMatrix(double * aVector, double *** rndMatrix, int RmatrixType){
	double resultvxM[size]; //working space variable

	switch(RmatrixType){
	case MATRIX_NONE:
		for (int i=0; i<size; i++)
			resultvxM[i] = aVector[i];
		break;
	case MATRIX_DIAGONAL || MATRIX_LINEAR: //Random diagonal matrix
	for (int i=0; i<size; i++)
		resultvxM[i] = aVector[i]* rndMatrix[0][i][i];
	break;
	case MATRIX_RRM_EXP_MAP:{ //Random rotation matrix using exponential method
		for (int i = 0 ; i < size ; i ++) {
			resultvxM[i] = 0.0;
			for (int j = 0 ; j < size ; j ++) {
				resultvxM[i] += (aVector[i] * rndMatrix[0][i][j]);
			}
		}
	}
	break;
	case MATRIX_RRM_EUCLIDEAN_ONE:{ //Random rotation matrix using Euclidean rotation (ONLY ONE PLANE)
		bool exitFlag = false;
		//copy the vector
		for (int i=0; i<size; i++){
			resultvxM[i] = aVector[i];
		}
		//find the planes to rotate
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++){
				if (i != j && rndMatrix[0][i][j] != 0){
					resultvxM[i] = (aVector[i] * rndMatrix[0][i][i]) + (aVector[j] * rndMatrix[0][j][i]);
					resultvxM[j] = (aVector[j] * rndMatrix[0][j][j]) + (aVector[i] * rndMatrix[0][i][j]);
					exitFlag = true;
					break;
				}
			}
			if (exitFlag == false)
				continue;
			else
				break;
		}
	}
	break;
	case MATRIX_RRM_EUCLIDEAN_ALL:{ //Random rotation matrix using Euclidean rotation (ALL POSSIBLE PLANES)
		bool exitFlag = false;
		//copy the vector
		for (int i=0; i<size; i++){
			resultvxM[i] = aVector[i];
		}
		int matNum=0;
		//The next two loops are to traverse the (size*size-1)/2 random matrices
		//needed to rotate the vector in all possible planes
		for (int g=0; g<size-1; g++) {
			for (int h=g+1; h<size; h++) {
				//find the planes to rotate
				for (int i=0; i<size; i++){
					for (int j=0; j<size; j++){
						if (i != j && rndMatrix[matNum][i][j] != 0){
							resultvxM[i] = (aVector[i] * rndMatrix[matNum][i][i]) + (aVector[j] * rndMatrix[matNum][j][i]);
							resultvxM[j] = (aVector[j] * rndMatrix[matNum][j][j]) + (aVector[i] * rndMatrix[matNum][i][j]);
							exitFlag = true;
							break;
						}
					}
					if (exitFlag == false)
						continue;
					else
						break;
				}
				matNum++;
			}
		}
	}
	break;
	}
	//	cout << "\n Original vector: [" ;
	//	for (int i=0; i<size; i++){
	//		cout << " " << aVector[i] << " ";
	//	}
	//	cout << "]";
	//	cout << "\n Vector rotated: [" ;
	for (int i=0; i<size; i++){
		aVector[i] = resultvxM[i];
		//		cout << " " << aVector[i] << " ";
	}
	//	cout << "]" << endl;

	//return the aVector
	return(aVector);
}

//Compute the random matrix to employ
void Particle::computeRndMatrix(double *** rndMatrix, int RmatrixType){
	switch(RmatrixType){
	case MATRIX_DIAGONAL: //Random diagonal matrix
		for (int i=0; i<size; i++)
			rndMatrix[0][i][i] = problem->getRandom01();
		break;
	case MATRIX_LINEAR:{ //Random linear matrix
		double rnd = problem->getRandom01();
		for (int i=0; i<size; i++)
			rndMatrix[0][i][i] = rnd;
	}
	break;
	case MATRIX_RRM_EXP_MAP:{ //Random rotation matrix using exponential method
		//1.- Generate a random matrix
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++)
				rndMatrix[0][i][j] = problem->getRandomX(-0.5,0.5);
		}
		//2.- Generate the transpose of the rndMatrix
		double trans_rndMatrix[size][size];
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++)
				trans_rndMatrix[j][i]=rndMatrix[0][i][j];
		}
		//3.- Determine the rotation angle
		double angle = (problem->getRandomX(0.001,7)*PI)/180; //rotation between 0 and 10 degrees
		//4.- Subtract trans_rndMatrix to rndMatrix, multiply by the angle and add the identity matrix
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++){
				rndMatrix[0][i][j]=(rndMatrix[0][i][j]-trans_rndMatrix[i][j])*angle;
				if (i==j)
					rndMatrix[0][i][j]+=1.0; //Add the Identity matrix
			}
		}
	}
	break;
	case MATRIX_RRM_EUCLIDEAN_ONE:{ //Random rotation matrix using Euclidean rotation (ONLY ONE PLANE)
		//1.- Get the rotation angle
		double angle = (problem->getRandomX(0.001,7)*PI)/180; //rotation between 0 and 10 degrees
		//cout << " angle:  " << angle << endl;
		//2.- Randomly select two different planes to rotate
		int plane1 = (int)floor(RNG::randVal(0.0,(double)size-1));
		int plane2;
		for (unsigned int i=0; i<size; i++) {
			plane2 = (int)floor(RNG::randVal(0.0,(double)size-1));
			if (plane2 == plane1)
				continue;
			else
				break;
		}
		//3.- Generate the Euclidean RRM
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++){
				if (i == j){
					if (i == plane1 || i == plane2)
						rndMatrix[0][i][j] = cos(angle);
					else
						rndMatrix[0][i][j] = 1;
				}
				else {
					if (i==plane1 && j==plane2)
						rndMatrix[0][i][j] = -sin(angle);
					else if (i==plane2 && j==plane1)
						rndMatrix[0][i][j] = sin(angle);
					else
						rndMatrix[0][i][j] = 0;
				}
			}
		}
		//4.- Print the matrix
		//		cout << "\n Euclidean RRM one plane: " << endl;
		//		for (int i=0; i<size; i++){
		//			cout << "[" ;
		//			for (int j=0; j<size; j++){
		//				if (i == j){
		//					if (i == plane1 || i == plane2)
		//						cout << " " << rndMatrix[0][i][j] << " "; //cout << " cos(alpha) ";
		//					else
		//						cout << " " << rndMatrix[0][i][j] << " "; //cout << " 1 ";
		//				}
		//				else {
		//					if (i==plane1 && j==plane2)
		//						cout << " " << rndMatrix[0][i][j] << " "; //cout << " -sin(angle) ";
		//					else if (i==plane2 && j==plane1)
		//						cout << " " << rndMatrix[0][i][j] << " "; //cout << " sin(angle) ";
		//					else
		//						cout << " " << rndMatrix[0][i][j] << " "; //cout << " 0 ";
		//				}
		//			}
		//			cout << " ]" << endl;
		//		}
	}
	break;
	case MATRIX_RRM_EUCLIDEAN_ALL:{ //Random rotation matrix using Euclidean rotation (ALL PLANES)
		int matNum=0; //we need a hypercube for each vector (this method is really demanding)
		//The maximum number of RRM is (size*size-1)/2
		for (int g=0; g<size-1; g++) {
			for (int h=g+1; h<size; h++) {
				//1.- Randomly choose an angle between 0 and 7
				double angle = (problem->getRandomX(0.001,7)*PI)/180; //rotation between 0 and 10 degrees
				//2.- Generate the Euclidean RRM
				for (int i=0; i<size; i++){
					for (int j=0; j<size; j++){
						if (i == j){
							if (i == g || i == h)
								rndMatrix[matNum][i][j] = cos(angle);
							else
								rndMatrix[matNum][i][j] = 1;
						}
						else {
							if (i==g && j==h)
								rndMatrix[matNum][i][j] = -sin(angle);
							else if (i==h && j==g)
								rndMatrix[matNum][i][j] = sin(angle);
							else
								rndMatrix[matNum][i][j] = 0;
						}
					}
				}
				matNum++;
			}
		}
		//		//4.- Print the matrix
		//		cout << "\n Euclidean RRM all planes (last plane): " << endl;
		//		for (int i=0; i<size; i++){
		//			cout << "[" ;
		//			for (int j=0; j<size; j++){
		//				if (i == j){
		//					if (i == size-2 || i == size-1)
		//						cout << " cos(alpha) ";
		//					else
		//						cout << " 1 ";
		//				}
		//				else {
		//					if (i==size-2 && j==size-1)
		//						cout << " -sin(angle) ";
		//					else if (i==size-1 && j==size-2)
		//						cout << " sin(angle) ";
		//					else
		//						cout << " 0 ";
		//				}
		//			}
		//			cout << " ]" << endl;
		//		}
	}
	break;
	case MATRIX_NONE:{
		for (int i=0; i<size; i++)
			rndMatrix[0][i][i] = 1.0;
	}
	break;
	}
}

void Particle::computeEvaluation() {
	//cout << "Evaluation: ";
	current.eval = problem->getFunctionValue(current.x);
	//cout << current.eval << endl;
	if (current.eval < pbest.eval) {
		for (int i=0;i<size;i++) {
			pbest.x[i] = current.x[i];
		}
		pbest.eval=current.eval;
	}
}

/* update global best solution
   INPUT: * coefficient vector x
 * corresponding solution value eval = f(x)
 */
void Particle::updateGlobalBest(double* x, double eval){
	for(int j=0;j <size;j++){
		gbest.x[j]= x[j];
	}
	gbest.eval = eval;
}

double Particle::computeDistPbestGbest(){
	double sumDistSqr = 0.0;
	for(int i=0;i<size;i++){
		sumDistSqr += pow((this->gbest.x[i]-this->pbest.x[i]),2);
	}
	return sqrt(sumDistSqr);
}

double Particle::computeDistance(double * x, double * p){
	double sumDistSqr = 0.0;
	for(int i=0;i<size;i++){
		sumDistSqr += pow((x[i]-p[i]),2);
	}
	return sqrt(sumDistSqr);
}

double* Particle::getCurrentPosition() {
	return(current.x);
}

long double Particle::getCurrentEvaluation(){
	return(current.eval);
}

double* Particle::getCurrentVelocity(){
	return velocity;
}

int Particle::getgBestID(){
	return gBestID;
}

double* Particle::getPbestPosition() {
	return(pbest.x);
}

long double Particle::getPbestEvaluation(){
	return(pbest.eval);
}

void Particle::printPosition(){
	cout << "Solution: " << current.eval << " Id: " << id << endl;
	for(int i=0; i<size; i++){
		cout << current.x[i] << "  ";
	}
	cout << endl;

}

void Particle::printNeighborByID(int identifier){
	cout << "Solution: " << neighbours.at(identifier)->current.eval << " Id: " << neighbours.at(identifier)->id << endl;
	for(int i=0; i<size; i++){
		cout << neighbours.at(identifier)->current.x[i] << "  ";
	}
	cout << endl;
}

/*Check the neighborhood for the best particle */
int Particle::getBestOfNeibourhood(){
	double aux_eval=this->gbest.eval;;
	int best=-1;
	//int prev_gBestID = this->gBestID;

	//Check first personal best
	if(this->pbest.eval < this->gbest.eval){
		updateGlobalBest(this->pbest.x, this->pbest.eval);
		this->gBestID = this->id;
	}
	//Check after the rest of the neighborhood
	for(unsigned int i=0; i<this->neighbours.size();i++){
		if(aux_eval > this->neighbours[i]->getPbestEvaluation()){
			best = i;
		}
	}

	//New best particle in the neighborhood
	if(best!=-1){
		updateGlobalBest(this->neighbours[best]->getPbestPosition(), this->neighbours[best]->getPbestEvaluation());
		this->gBestID = this->neighbours[best]->getID();
	}
	//if (prev_gBestID != best)
	//cout << "New gBestID for particle -- " << this->id << endl;//remove

	return best;
}

int Particle::getRandomNeighbor(){
	int randomIndex = -1;
	for (unsigned int i=0; i<neighbours.size(); i++) {
		randomIndex = (int)floor(RNG::randVal(0.0,(double)neighbours.size()));

		//Distinct of itself
		if (neighbours.at(randomIndex)->id != id)
			break;
		else
			continue;
	}
	if (randomIndex == -1){
		cout << "Couldn't find any :/" << endl;
		exit (EXIT_FAILURE);
		return id;
	}
	else
		return randomIndex;
}


void Particle::addNeighbour(Particle* p){
	neighbours.push_back(p);
}

//From Frankesntein's PSO
unsigned int Particle::getNeighborhoodSize(){
	return neighbours.size();
}

int Particle::getID(){
	return id;
}

void Particle::eraseNeighborbyID(int nid){
	for(unsigned int i=0;i<neighbours.size();i++){
		if(nid == neighbours.at(i)->id){
			neighbours.erase(neighbours.begin()+i);
			break;
		}
	}
}

int Particle::getRandomNonAdjacentNeighborID(Configuration* config){
	int randomIndex,randomID;
	if(neighbours.size()>3){
		while(true){
			randomIndex = (int)floor(RNG::randVal(0.0,(double)neighbours.size()));
			//One random neighbor
			randomID = neighbours.at(randomIndex)->id;

			//First position of the vector
			if(id==0){
				//Distinct of itself, the right-side neighbor and the last neighbor
				if(randomID != id && randomID != id+1 && randomID != config->getSwarmSize()-1)
					return randomID;
				else
					continue;
			}
			else{
				//Last position of the vector
				if(id == config->getSwarmSize()-1){
					//Distinct of itself, the first neighbor and the left-side neighbor
					if(randomID != id && randomID != 0 && randomID != id-1)
						return randomID;
					else
						continue;
				}
				else{
					//Some position in within the bounds just check right- and left-side neighbors
					if(randomID != id && randomID != id+1 && randomID != id-1)
						return randomID;
					else
						continue;
				}
			}
		}
	}
	else
		return -1;
}

//int Particle::getParticle_gBestID(){
//	double aux_eval = this->gbest.eval;
//	int best= this->id;
//
//	//Check the neighborhood
//	for(unsigned int i=0; i<this->neighbours.size(); i++){
//		if(aux_eval > this->neighbours[i]->getPbestEvaluation()){
//			best = this->neighbours[i]->getID();
//		}
//	}
//	return best;
//}

///*Check the neighborhood for the best particle */
//void Particle::checkNeibourhood(){
//	double aux_eval;
//	int best=-1;
//
//	//Check first personal best
//	if(pbest.eval < gbest.eval){
//		updateGlobalBest(pbest.x, pbest.eval);
//	}
//	//Check after the rest of the neighborhood
//	aux_eval=gbest.eval;
//	for(unsigned int i=0; i<neighbours.size();i++){
//		if(aux_eval > neighbours[i]->getPbestEvaluation()){
//			best = i;
//		}
//	}
//
//	//New best particle in the neighborhood
//	if(best!=-1)
//		updateGlobalBest(neighbours[best]->getPbestPosition(), neighbours[best]->getPbestEvaluation());
//}

//double Particle::computeNewVelocity(Configuration* config, double vel, double u1, double u2, double perInf, double socInf, double pos, double additionalVal){
//	double new_vel = 0.0;
//
//	//Original PSO
//	if (config->getVelocityRule() == VEL_BASIC){
//		new_vel = vel +
//				(phi_1 * problem->getRandom01() * perInf) +
//				(phi_2 * problem->getRandom01() * socInf);
//	}
//	//Standard PSO
//	else if (config->getVelocityRule() == VEL_STANDARD){
//		new_vel = (inertia * vel) +
//				(phi_1 * problem->getRandom01() * perInf) +
//				(phi_2 * problem->getRandom01() * socInf);
//
//	}
//	//Linear PSO
//	else if (config->getVelocityRule() == VEL_LINEAR){
//		new_vel = (inertia * vel) +
//				(phi_1 * u1 * perInf) +
//				(phi_2 * u2 * socInf);
//	}
//	//Constriction coefficient PSO
//	else if (config->getVelocityRule() == VEL_CONSTRICTED){
//		new_vel = CONSTRICTION_COEFFICIENT * (vel +
//				(phi_1 * problem->getRandom01() * perInf) +
//				(phi_2 * problem->getRandom01() * socInf));
//	}
//	else if (config->getVelocityRule() == VEL_GUARAN_CONVERG){
//		//if (pbest.eval == gbest.eval)
//	}
//	//Standard 2011 PSO - Uses a Hypersphere distribution to generate new points
//	else if (config->getVelocityRule() == VEL_STANDARD2011){
//		//in this case additionalVal = (H(G||G-X||) - G-X) -> see StandarPSO Maurice Clerc for details
//		new_vel = (inertia * vel) + additionalVal;
//	}
//	//Use standard PSO
//	else
//		new_vel = (inertia * vel) +
//		(phi_1 * problem->getRandom01() * perInf) +
//		(phi_2 * problem->getRandom01() * socInf);
//
//
//	return new_vel;
//}


//// The computation of the radius and the random point in the HyperSphere
//// was taken from the publicly available from Maurice Clerc - Standard PSO 2011
//// https://www.particleswarm.info/Programs.html
//void Particle::getHypersphericalVector(double V2[], double V1[]){
//	double G[size];	//center of the sphere
//	double l[size]; //particle's Gbest
//	double radius = 0.0;	//radius G-X
//	double pw=1./(double)size;
//
//	//check that Pbest != Gbest
//	if (pbest.eval == gbest.eval){
//		//use a random neighbor as Gbest
//		int randNeighbor = getRandomNeighbor();
//		for (int i=0; i<size; i++)
//			l[i] = neighbours.at(randNeighbor)->current.x[i];
//	}
//	else
//		for (int i=0; i<size; i++)
//			l[i] = gbest.x[i];
//
//	//Compute G (center of the sphere) and V1 (radius of each dimension)
//	for (int i=0; i<size; i++){
//		double P = current.x[i] + (phi_1 * (pbest.x[i]-current.x[i]));
//		double L = current.x[i] + (phi_2 * (l[i]-current.x[i]));
//		G[i] = (current.x[i] + P + L)/3.0;
//		radius += pow(abs(current.x[i] - G[i]), 2);
//		V1[i] = G[i] - current.x[i];
//	}
//	radius = sqrt(radius); //this is the actual radius of the hyper-sphere
//
//	// Get a random vector in the hyper-sphere H(G||G-X||)
//	// ----------------------------------- Step 1.  Direction
//	double length=0.0;
//	for (int i=0;i<size;i++) {
//		V2[i] = RNG::randGauss(1.0);
//		length += pow(V2[i],2);
//	}
//	length=sqrt(length);
//	//----------------------------------- Step 2. Random radius
//	// Random uniform distribution on the sphere
//	double r = pow(RNG::randVal(0.0,1.0),pw);
//
//	//----------------------------------- Step 3. Random vector
//	for (int i=0;i<size;i++) {
//		V2[i] = radius*r*V2[i]/length;
//	}
//}


//// The computation of the radius and the random point in the HyperSphere
//// was taken from the publicly available code of Maurice Clerc - Standard PSO 2011
//// https://www.particleswarm.info/Programs.html
//void Particle::getHypersphericalVector(double V2[], double V1[], int numInformants, int *theInformants, double *** rndMatrix, int RmatrixType){
//	double G[size];	//center of the sphere
//	double l[size]; //particle's Gbest
//	double radius = 0.0;	//radius G-X
//	double pw=1./(double)size;
//	double * v_PosToInfPbest[numInformants];
//	bool pBestInformant = isPbestInformant(numInformants, theInformants); //flag to know if pBest is included in theInformants
//
//	//When the hierarchical topology is used, this verification might be irrelevant most of the time
//	//since the gBest particle can be located in a lower level in the hierarchy.
//	//Check if the particle is gBest
//	if (this->id == this->gBestID){
//		//use a random neighbor as Gbest
//		int randNeighbor = getRandomNeighbor();
//		if (randNeighbor != this->id){
//			cout << "\t\t ... using a random neighbor's pbest--" << endl;
//			for (int i=0; i<size; i++)
//				l[i] = neighbours.at(randNeighbor)->pbest.x[i]; //Use random neighbor's pBest
//		}
//		//initialization rule of Incremental PSO -- reinitialize position to model
//		else{
//			cout << "\t\t ... using reinitialization to model --" << endl;
//			for (int i=0; i<size; i++)
//				l[i] = problem->getRandomX() + problem->getRandom01()*(gbest.x[i]-current.x[i]);
//		}
//	}
//	else
//		for (int i=0; i<size; i++)
//			l[i] = gbest.x[i];
//
//	//Rotate if some type of rotation is chosen
//	//	if (RmatrixType != MATRIX_NONE){
//	//If gBest is not included in theInformant array, it will be added at the end
//	if (pBestInformant == false)
//		numInformants = numInformants+1;
//	//Create a temporary structure
//	for (int i=0; i<numInformants; i++)
//		v_PosToInfPbest[i] = new double [size];
//
//	//If pBest was not included
//	if (pBestInformant == false){
//		//Copy the informants in a temporary structure
//		//cout << "\t\t---Original vectors pk-xi --" << endl;
//		for (int j=0; j<numInformants-1; j++){
//			//cout << "\t p[" << theInformants[j] << "] -> { ";
//			if (this->id == this->gBestID){ //use the right gBest
//				//cout << " l[] ";
//				for (int i=0; i<size; i++){
//					//this is the vector we want to rotate
//					v_PosToInfPbest[j][i] = (l[i]-current.x[i]);
//					//cout << v_PosToInfPbest[j][i] << " ";
//				}
//			}
//			else
//				for (int i=0; i<size; i++){
//					//this is the vector we want to rotate
//					v_PosToInfPbest[j][i] = (neighbours.at(theInformants[j])->pbest.x[i]-current.x[i]);
//					//cout << v_PosToInfPbest[j][i] << " ";
//				}
//			//cout << "}" << endl;
//		}
//		for (int i=0; i<size; i++){
//			//this is the vector we want to rotate
//			v_PosToInfPbest[numInformants-1][i] = (pbest.x[i]-current.x[i]);
//			//cout << v_PosToInfPbest[j][i] << " ";
//		}
//	}
//	//pBest was already included
//	else {
//		//Copy the informants in a temporary structure
//		//cout << "\t\t---Original vectors pk-xi --" << endl;
//		for (int j=0; j<numInformants; j++){
//			//cout << "\t p[" << theInformants[j] << "] -> { ";
//			if (this->id == this->gBestID){ //use the right gBest
//				//cout << " l[] ";
//				for (int i=0; i<size; i++){
//					//this is the vector we want to rotate
//					v_PosToInfPbest[j][i] = (l[i]-current.x[i]);
//					//cout << v_PosToInfPbest[j][i] << " ";
//				}
//			}
//			else
//				for (int i=0; i<size; i++){
//					//this is the vector we want to rotate
//					v_PosToInfPbest[j][i] = (neighbours.at(theInformants[j])->pbest.x[i]-current.x[i]);
//					//cout << v_PosToInfPbest[j][i] << " ";
//				}
//			//cout << "}" << endl;
//		}
//	}
//	//Rotate each vector
//	//cout << "\t\t---Rotated vectors M * (pk-xi) --" << endl;
//	for (int j=0; j<numInformants; j++){ //the rest of informants
//		//cout << "\t p[" << theInformants[j] << "] -> { ";
//		multiplyVectorByRndMatrix(v_PosToInfPbest[j], rndMatrix, RmatrixType);
//		for (int i=0; i<size; i++){
//			//cout << v_PosToInfPbest[j][i] << " ";
//		}
//		//cout << "}" << endl;
//		if (j<numInformants)
//			computeRndMatrix(rndMatrix, RmatrixType); //compute a random matrix for the next informant
//	}
//	//	}
//
//	//Compute G (center of the sphere) and V1 (radius of each dimension)
//	for (int i=0; i<size; i++){
//		//		if (RmatrixType == MATRIX_NONE){
//		//			double R = 0.0;
//		//			double P = 0.0;
//		//			if (pBestInformant == true){ //pBest is included in theInformants
//		//				for (int j=0; j<numInformants; j++){
//		//					if (neighbours.at(theInformants[j])->getID() == this->gBestID)
//		//						R += current.x[i] + (phi_1 * (l[i]-current.x[i])); //always use l[]
//		//					else
//		//						R += current.x[i] + (phi_1 * (neighbours.at(theInformants[j])->pbest.x[i]-current.x[i]));
//		//				}
//		//				G[i] = (current.x[i] + R )/(numInformants + 1);
//		//			}
//		//			else{ //pBest is not included in theInformants
//		//				if (this->id == this->gBestID)
//		//					P = current.x[i] + (phi_1 * (l[i]-current.x[i])); //pBest == gBest: use l[]
//		//				else
//		//					P = current.x[i] + (phi_1 * (pbest.x[i]-current.x[i])); //pBest
//		//				//(Note that l[] will be P or R, but not both)
//		//				for (int j=0; j<numInformants; j++){
//		//					if (neighbours.at(theInformants[j])->getID() == this->gBestID) //pBest == gBest: use l[]
//		//						R += current.x[i] + (phi_1 * (l[i]-current.x[i])); //always use l[]
//		//					else
//		//						R += current.x[i] + (phi_1 * (neighbours.at(theInformants[j])->pbest.x[i]-current.x[i]));
//		//				}
//		//				G[i] = (current.x[i] + P + R )/(numInformants + 2);
//		//			}
//		//		}
//		//		else {
//		//Use the rotated vectors
//		double R = 0.0;
//		for (int j=0; j<numInformants; j++){
//			//cout << "\t\t  informant at [" << j << "] is " <<  theInformants[j] << endl; //remove
//			R += current.x[i] + (phi_1 * v_PosToInfPbest[j][i]); //v_PosToInfPbest is the vector already rotated
//		}
//		G[i] = (current.x[i] + R )/(numInformants + 1);
//		//		}
//		radius += pow(abs(current.x[i] - G[i]), 2);
//		V1[i] = G[i] - current.x[i];
//	}
//	radius = sqrt(radius); //this is the actual radius of the hyper-sphere
//
//	// Get a random vector in the hyper-sphere H(G||G-X||)
//	// ----------------------------------- Step 1.  Direction
//	double length=0.0;
//	for (int i=0;i<size;i++) {
//		V2[i] = RNG::randGauss(1.0);
//		length += pow(V2[i],2);
//	}
//	length=sqrt(length);
//	//----------------------------------- Step 2. Random radius
//	// Random uniform distribution on the sphere
//	double r = pow(RNG::randVal(0.0,1.0),pw);
//
//	//----------------------------------- Step 3. Random vector
//	for (int i=0;i<size;i++) {
//		V2[i] = radius*r*V2[i]/length;
//	}
//
//	//if (pBestInformant == false && RmatrixType != MATRIX_NONE)
//	if (pBestInformant == false)
//		numInformants-=1;
//}
