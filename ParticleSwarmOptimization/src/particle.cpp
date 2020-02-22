/*
 * particle.cpp
 *
 *  Created on: May 31, 2018
 *      Author: leonardo
 */
//cout << "\n So far, so good" << endl; //remove

#include "float.h"
#include "particle.h"
#include "iostream"
#include "rng.h"
#include <math.h>
#include <cmath>
#include <new>
#include "utils.h"

using namespace std;

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
Particle::Particle (Problem* problem, Configuration* config, int identifier, long int iteration){

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
	inertia = config->getOmega1();

	//cout << "Particle " << id << " : [ ";
	initializePosition(config, iteration);

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

	phi_1 = p.phi_1;
	phi_2 = p.phi_2;
	inertia=p.inertia;

	init = true;

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

		phi_1 = p.phi_1;
		phi_2 = p.phi_2;
		inertia=p.inertia;

		init = true;
	}
	return *this;
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


/* Initialize particle using uniformly random values */
void Particle::initUniform(Configuration* config){
	if (config->verboseMode()) cout << "\tvector::(initial)::[ " ;
	for(int i=0; i<size; i++){
		if(config->getCompetitionID() == CEC05 && config->getProblemID() == SHIFTED_ROTATED_GRIEWANK_CEC05)
			current.x[i] = problem->getRandomX(0.0,600.0);
		else if(config->getCompetitionID() == CEC05 && config->getProblemID() == ROTATED_HYBRIDCOMPOSITION4_NO_BOUNDS)
			current.x[i] = problem->getRandomX(2.0,5.0);
		else
			current.x[i] = problem->getRandomX(); //random values within the bounds of the function
		pbest.x[i]=current.x[i];
		velocity[i]=0;
		if (config->verboseMode()) cout << current.x[i] << " ";
	}
	if (config->verboseMode()) cout << "]" << endl;
	evaluateSolution();
}

/* Initialize particle using a model (as in Incremental PSO) */
void Particle::initToModel(){
	for (int i=0; i<size; i++){
		current.x[i] = current.x[i] + problem->getRandom01()*(gbest.x[i]-current.x[i]);
	}
	evaluateSolution();
}

void Particle::initializePosition(Configuration* config, long int iteration){
	switch (config->getParticleInitType()) {
	case PARTICLE_INIT_RANDOM:
		initUniform(config);
		break;
	case PARTICLE_INIT_MODEL:{ //Only available for dynamic PopCS and after the first iteration
		initUniform(config);
		if (config->getPopulationCS() != POP_CONSTANT && iteration > 0)
			initToModel();
	}
	break;
	default:
		initUniform(config);
		break;
	}
}

void Particle::setVelocityLimits(Configuration* config){
	minVelLimit=((config->getMaxInitRange()-config->getMinInitRange())/2.0)*-1;
	maxVelLimit=((config->getMaxInitRange()-config->getMinInitRange())/2.0);
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

void Particle::setPhi1(double new_phi_1){
	phi_1 = new_phi_1;
}
void Particle::setPhi2(double new_phi_2){
	phi_2 = new_phi_2;
}
double Particle::getPhi1(){
	return phi_1;
}
double Particle::getPhi2(){
	return phi_2;
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

bool Particle::ispBestIntheInformants(int numInformants){
	bool result = false;
	for (unsigned int j=0; j<InformantsPos.size(); j++)
		if (this->id == neighbours.at(InformantsPos[j])->getID() )
			result = true;

	return result;
}

/* Generate a new solution by updating the particle's position */
void Particle::move(Configuration* config, double minBound, double maxBound, long int iteration, double omega1, double omega2, double omega3,
		int numInformants, int lastLevelComplete, double alpha_t, double l, double delta){

	double vect_distribution[size];
	double vect_perturbation[size];

	bool pBestIntheInformants = ispBestIntheInformants(numInformants);
	//Since we are putting everything in one single structure we need to check if pBest is in the Informants or not
	//to avoid including it twice.
	if (!pBestIntheInformants)
		numInformants = numInformants+1;
	if (config->verboseMode()) cout << "\tvar::numInformants: " << numInformants << " ";
	vector<vector< double> > vect_PbestMinusPosition;
	vect_PbestMinusPosition.resize(numInformants, vector<double>(size));

	/*** PERTURBATION 1 (distribution-based)	--->	It has to be computed per Informant and/or per Dimension depending on the strategy
	 *** 						this component is applied directly in the DNNP members of DNPP  ***/
	/*** PERTURBATION 2	(additive)	--->	It has to be computed per dimension, but only once per particle  ***/
	for (int i=0;i<size;i++)
		vect_perturbation[i] = getPerturbationMagnitude(config->getPerturbation2Type(), alpha_t, delta); //Only additive perturbation


	/*** DISTRIBUTION VECTOR (DNNPs)***/
	switch (config->getDistributionNPP()) {
	case DIST_RECTANGULAR:{
		computeSubtractionPerturbationRotation(
				config,
				vect_PbestMinusPosition,
				numInformants,
				pBestIntheInformants,
				alpha_t,
				l);
		getRectangularDNPP(
				vect_distribution,
				numInformants,
				pBestIntheInformants,
				vect_PbestMinusPosition,
				config->getModelOfInfluence());
		break;
	}
	case DIST_SPHERICAL:{
		computeSubtractionPerturbationRotation(
				config,
				vect_PbestMinusPosition,
				numInformants,
				pBestIntheInformants,
				alpha_t,
				l);
		getSphericalDNPP(
				vect_distribution,
				numInformants,
				pBestIntheInformants,
				vect_PbestMinusPosition,
				config->getModelOfInfluence());
		break;
	}
	case DIST_ADD_STOCH:{
		computeSubtractionPerturbationRotation(
				config,
				vect_PbestMinusPosition,
				numInformants,
				pBestIntheInformants,
				alpha_t,
				l);
		getAdditiveStochasticDNPP(
				vect_distribution,
				numInformants,
				pBestIntheInformants,
				vect_PbestMinusPosition,
				config->getRandNeighbor(),
				config->getOperator_q());
		break;
	}
	}

	//Compute new position
	for (int i=0;i<size;i++) {
		velocity[i] = omega1 *  velocity[i] +
				omega2 * vect_distribution[i] +
				omega3 * vect_perturbation[i];

		//Clamp velocity
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

		//Clamp position
		if(current.x[i] < minBound)
			current.x[i] = minBound;
		if(current.x[i] > maxBound)
			current.x[i]= maxBound;
	}

	//Evaluate the objective function and update pbest if a new one has been found
	evaluateSolution();
	if (config->verboseMode()) cout << "\tParticle with ID:[" << this->id << "].status::MOVED" << endl;
	if (config->verboseMode()) cout << "\t------------------------------------------" << endl;

}

void Particle::computeSubtractionPerturbationRotation(
		Configuration* config,
		vector<vector< double> > &vect_PbestMinusPosition,
		int &numInformants,
		bool pBestIntheInformants,
		double alpha_t,
		double l_value) {

	double l[size]; //particle's Gbest
	bool increase_numInformants = false;
	int pertubType = config->getPerturbation1Type();

	if (config->verboseMode()){
		(pBestIntheInformants) ?
				cout << "\n\tvar::this->id: " << this->id << "\n\tvar::pBestIntheInformants: TRUE"
				<< " \n\tvar::DNPP: " << config->getDistributionNPP() << endl :
				cout << "\n\tvar::this->id: " << this->id << "\n\tvar::pBestIntheInformants: FALSE"
				<< " \n\tvar::DNPP: " << config->getDistributionNPP() << endl;
	}

	//1.- Check if the particle is pBest == gBest
	if (this->id == this->gBestID && pBestIntheInformants){
		if (config->verboseMode()) cout << "\t\tnotice::using reinitialization to model for particle "
				<< this->id << " -The new position (l[]) is somewhere around gBest)-" << endl;
		for (int i=0; i<size; i++)
			l[i] = problem->getRandomX() + problem->getRandom01()*(gbest.x[i]-current.x[i]);
	}

	if (config->getDistributionNPP() == DIST_ADD_STOCH) {
		//2.- Get the p^k of all Informants and then add pBest as the end of the Array
		if (!pBestIntheInformants){
			for (int j=0; j<numInformants-1; j++){ //Copy the informants in a temporary structure
				setPerturbationMagnitude(pertubType, current.x, neighbours[InformantsPos[j]]->pbest.x, alpha_t, l_value); //informant-wise
				for (int i=0; i<size; i++){
					vect_PbestMinusPosition.at(j).at(i) = (applyPerturbation(pertubType, neighbours.at(InformantsPos[j])->pbest.x[i])); //vector to be rotated
					setPerturbationMagnitude(pertubType, alpha_t); //dimensional-wise
				}
			}
			//Add pBest at the end of the Array
			setPerturbationMagnitude(pertubType, current.x, pbest.x, alpha_t, l_value); //informant-wise
			for (int i=0; i<size; i++){
				vect_PbestMinusPosition.at(numInformants-1).at(i) = (applyPerturbation(pertubType, pbest.x[i])); //vector to be rotated
				setPerturbationMagnitude(pertubType, alpha_t); //dimensional-wise
			}
		}
		else {
			//Note that in this case, when gBest = pBest, we would be using only one of these two informants.
			//However, for a model of influence such as best-of-neighborhood, this will actually be limiting
			//the number of informants to 1, which could provide only little information to create the random
			//vector in hyper-sphere. Therefore, we use l[], which is a position close to gBest created using
			//the method propose in Incremental PSO as gBest
			for (int j=0; j<numInformants; j++){ //Copy the informants in a temporary structure
				if ((neighbours.at(InformantsPos[j])->getID() == this->gBestID) && (this->id == this->gBestID)){
					//use gBest in the form of l[]
					setPerturbationMagnitude(pertubType, current.x, l, alpha_t, l_value); //informant-wise perturbation magnitude
					for (int i=0; i<size; i++){
						vect_PbestMinusPosition.at(j).at(i) = (applyPerturbation(pertubType, l[i]));  //vector to be rotated
						setPerturbationMagnitude(pertubType, alpha_t); //dimensional-wise perturbation magnitude
					}
					//create one extra space in vect_PbestMinusPosition for l[]
					vect_PbestMinusPosition.resize(numInformants+1, vector<double>(size));
					increase_numInformants = true;
					//add pBest at the end of the array
					setPerturbationMagnitude(pertubType, current.x, pbest.x, alpha_t, l_value); //informant-wise perturbation magnitude
					for (int i=0; i<size; i++){
						vect_PbestMinusPosition.at(numInformants).at(i) = (applyPerturbation(pertubType, pbest.x[i]));  //vector to be rotated
						setPerturbationMagnitude(pertubType, alpha_t); //dimensional-wise perturbation magnitude
					}
				}
				else{
					setPerturbationMagnitude(pertubType, current.x, neighbours[InformantsPos[j]]->pbest.x, alpha_t, l_value); //informant-wise
					for (int i=0; i<size; i++){
						vect_PbestMinusPosition.at(j).at(i) = (applyPerturbation(pertubType, neighbours.at(InformantsPos[j])->pbest.x[i]));  //vector to be rotated
						setPerturbationMagnitude(pertubType, alpha_t); //dimensional-wise
					}
				}
			}
			if (increase_numInformants)
				numInformants++;
		}
	}
	else {
		//2.- Get the p^k-x^i of all Informants and then add pBest as the end of the Array
		if (!pBestIntheInformants){
			for (int j=0; j<numInformants-1; j++){ //Copy the informants in a temporary structure
				setPerturbationMagnitude(pertubType, current.x, neighbours[InformantsPos[j]]->pbest.x, alpha_t, l_value); //informant-wise
				for (int i=0; i<size; i++){
					vect_PbestMinusPosition.at(j).at(i) = applyPerturbation(pertubType, neighbours.at(InformantsPos[j])->pbest.x[i])-current.x[i]; //vector to be rotated
					setPerturbationMagnitude(pertubType, alpha_t); //dimensional-wise
				}
			}
			//Add pBest at the end of the Array
			setPerturbationMagnitude(pertubType, current.x, pbest.x, alpha_t, l_value); //informant-wise
			for (int i=0; i<size; i++){
				vect_PbestMinusPosition.at(numInformants-1).at(i) = (applyPerturbation(pertubType, pbest.x[i])-current.x[i]); //vector to be rotated
				setPerturbationMagnitude(pertubType, alpha_t); //dimensional-wise
			}
		}
		else {
			//Note that in this case, when gBest = pBest, we would be using only one of these two informants. However, for a model of influence such as
			//best-of-neighborhood, this will actually be limiting the number of informants to 1, which can't be little information to create the random
			//vector in hypersphere. Therefore, we use l[], which is a position close to gBest created using the method propose in Incremental PSO.
			for (int j=0; j<numInformants; j++){ //Copy the informants in a temporary structure
				if ((neighbours.at(InformantsPos[j])->getID() == this->gBestID) && (this->id == this->gBestID)){
					//use gBest in the form of l[]
					setPerturbationMagnitude(pertubType, current.x, l, alpha_t, l_value); //informant-wise perturbation magnitude
					for (int i=0; i<size; i++){
						vect_PbestMinusPosition.at(j).at(i) = (applyPerturbation(pertubType, l[i])-current.x[i]);  //vector to be rotated
						setPerturbationMagnitude(pertubType, alpha_t); //dimensional-wise perturbation magnitude
					}
					//create one extra space in vect_PbestMinusPosition for l[]
					vect_PbestMinusPosition.resize(numInformants+1, vector<double>(size));
					increase_numInformants = true;
					//add pBest at the end of the array
					setPerturbationMagnitude(pertubType, current.x, pbest.x, alpha_t, l_value); //informant-wise perturbation magnitude
					for (int i=0; i<size; i++){
						vect_PbestMinusPosition.at(numInformants).at(i) = (applyPerturbation(pertubType, pbest.x[i])-current.x[i]);  //vector to be rotated
						setPerturbationMagnitude(pertubType, alpha_t); //dimensional-wise perturbation magnitude
					}
				}
				else{
					setPerturbationMagnitude(pertubType, current.x, neighbours[InformantsPos[j]]->pbest.x, alpha_t, l_value); //informant-wise
					for (int i=0; i<size; i++){
						vect_PbestMinusPosition.at(j).at(i) = (applyPerturbation(pertubType, neighbours.at(InformantsPos[j])->pbest.x[i])-current.x[i]);  //vector to be rotated
						setPerturbationMagnitude(pertubType, alpha_t); //dimensional-wise
					}
				}
			}
			if (increase_numInformants)
				numInformants++;
		}
		//3.- Rotate each vector
		/*** ROTATION MATRICES	--->	has to be computed per Informant ***/
		int numMatrices =(int)floor(((size*(size-1))/2)); //number of rotation matrices to rotate in all possible planes
		double ** rndMatrix[numMatrices];
		//Allocate memory to be used
		for(int i=0; i<numMatrices; i++){
			rndMatrix[i] = new double*[size];
			for (unsigned int j=0; j<size; j++)
				rndMatrix[i][j] = new double[size];
		}
		computeRndMatrix(rndMatrix, config->getRandomMatrix(), config->getRotationAgle()); //we need to compute this for each informant
		for (int j=0; j<numInformants; j++){ //the rest of informants
			multiplyVectorByRndMatrix(vect_PbestMinusPosition, j, rndMatrix, config->getRandomMatrix());
			if (j<numInformants) //avoid computing a new random matrix in the last iteration
				computeRndMatrix(rndMatrix, config->getRandomMatrix(), config->getRotationAgle()); //compute a random matrix for the next informant
		}
		//Deallocate memory used
		for(int i=0; i<numMatrices; i++){
			for (unsigned int j=0; j<size; j++)
				delete [] rndMatrix[i][j];
			delete [] rndMatrix[i];
		}
	}
}

//This function returns a random position in theInformants different from the pBest of the particle
int Particle::getRandomInformantPosition(int numInformants, bool pBestIntheInformants){
	int randomIndex = 0;

	//Note that numInformants (independently of the topology and model of influences) is at least 2,
	//that is, theInformants contains always pBest in some position and some other particle in other position
	if (pBestIntheInformants && (this->id != this->gBestID)){ //pBest somewhere in the Array we don't know
		for (int i=0; i<numInformants; i++) {
			randomIndex = (int)floor(RNG::randVal(0.0,(double)numInformants));

			//Distinct of itself
			if (neighbours.at(InformantsPos[randomIndex])->getID() != id){
				break;
			}
			else
				continue;
		}
	}
	else {
		for (int i=0; i<numInformants-1; i++) { //pBest is at the end of the Array
			randomIndex = (int)floor(RNG::randVal(0.0,(double)numInformants-1));
		}
	}
	return randomIndex;
}

//This function returns the position of the particle's pBest in theInformants array
int Particle::getPositionOfpBest(int numInformants, bool pBestIntheInformants){
	int pBestIndex = 0;

	//Note that numInformants (independently of the topology and model of influences) is at least 2,
	//that is, theInformants contains always pBest in some position and some other particle in other position
	if (pBestIntheInformants && (this->id != this->gBestID)){ //pBest somewhere in the Array we don't know
		for (int i=0; i<numInformants; i++) {
			//Distinct of itself
			if (neighbours.at(InformantsPos[i])->getID() == id){
				pBestIndex = i;
				break;
			}
			else
				continue;
		}
	}
	else
		pBestIndex = numInformants-1;

	return pBestIndex;
}

void Particle::getAdditiveStochasticDNPP(double vect_distribution[], int numInformants, bool pBestIntheInformants,
		vector<vector< double> > &vect_PbestMinusPosition, bool randNeighbor, int operatorQ){

	int p2Index;
	int p1Index = getPositionOfpBest(numInformants, pBestIntheInformants);

	if (randNeighbor)
		p2Index = getRandomInformantPosition(numInformants, pBestIntheInformants);
	else
		p2Index = 0;

	//cout << "\t\tnotice:simply dynamic PSO, p1 : [" << p1Index << "]" << "\t p2 : [" << p2Index << "]" <<endl;

	switch (operatorQ) {
	case Q_STANDARD:
		for (int i=0; i<size; i++){
			vect_distribution[i] = (
					(phi_1 * vect_PbestMinusPosition[p1Index][i])+
					(phi_2 * vect_PbestMinusPosition[p2Index][i]) ) / (phi_1+phi_2);
			vect_distribution[i] -= current.x[i];
		}
		break;
	case Q_GAUSSIAN:
		for (int i=0; i<size; i++){
			double center = (vect_PbestMinusPosition[p1Index][i]+ vect_PbestMinusPosition[p2Index][i]) / 2.0;
			double dispersion = abs (vect_PbestMinusPosition[p1Index][i] - vect_PbestMinusPosition[p2Index][i]);
			vect_distribution[i] = RNG::randGaussWithMean(dispersion, center);
			vect_distribution[i] -= current.x[i];
		}
		break;
	case Q_DISCRETE_2:
		for (int i=0; i<size; i++){
			int rndBernoulli = RNG::randBernoulli(0.5); //toss a coin
			vect_distribution[i] =(((1 + rndBernoulli ) * vect_PbestMinusPosition[p1Index][i]) +
					((1 - rndBernoulli) * vect_PbestMinusPosition[p2Index][i]))/2;
			vect_distribution[i] -= current.x[i];
		}
		break;
	case Q_CAUCHY_NORMAL:
		for (int i=0; i<size; i++){
			if (problem->getRandom01() <= 0.5){
				//to see how the gamma parameter is set see RNG::randCauchy()
				vect_distribution[i] = vect_PbestMinusPosition[p1Index][i] - current.x[i]; //we need to discount current.x[i] because it will be added later in the GVU formula
				vect_distribution[i] += (RNG::randCauchy(1.0) * abs(vect_PbestMinusPosition[p1Index][i] -
						vect_PbestMinusPosition[p2Index][i]));

			}
			else {
				vect_distribution[i] = vect_PbestMinusPosition[p2Index][i] - current.x[i]; //we need to discount current.x[i] because it will be added later in the GVU formula
				vect_distribution[i] += (RNG::randGaussWithMean(1.0,0) * abs(vect_PbestMinusPosition[p1Index][i] -
						vect_PbestMinusPosition[p2Index][i]));
			}
		}
		break;
	}
}

void Particle::getRectangularDNPP(double vect_distribution[], int numInformants, bool pBestIntheInformants,
		vector<vector< double> > &vect_PbestMinusPosition, int modelOfInflu){

	double varPhi2 = phi_2;
	//Compute vect_distribution
	for (int i=0; i<size; i++){
		//		if (i==0) cout << "\tPhi2 values: [ ";
		for (int j=0; j<numInformants; j++){
			if (pBestIntheInformants && (this->id != this->gBestID)) { //we look for pBest in the Array
				if (this->id == neighbours.at(InformantsPos[j])->getID())
					vect_distribution[i] += (phi_1 * vect_PbestMinusPosition[j][i]); //personal coefficient phi_1
				else {
					if (j==0 && modelOfInflu == MOI_RANKED_FI)
						varPhi2 = phi_2*((double)(numInformants-1)/numInformants);
					if (j>0 && modelOfInflu == MOI_RANKED_FI)
						varPhi2 = varPhi2/2.0;
					vect_distribution[i] += (varPhi2 * vect_PbestMinusPosition[j][i]); //social coefficient phi_1
					//					if (i==0) cout << varPhi2 << " ";
				}
			}
			else{
				if (j == numInformants-1) //pBest is at the end of the Array
					vect_distribution[i] += (phi_1 * vect_PbestMinusPosition[j][i]); //personal coefficient phi_1
				else {
					if (j==0 && modelOfInflu == MOI_RANKED_FI)
						varPhi2 = phi_2*((double)(numInformants-1)/numInformants);
					if (j>0 && modelOfInflu == MOI_RANKED_FI)
						varPhi2 = varPhi2/2.0;
					vect_distribution[i] += (varPhi2 * vect_PbestMinusPosition[j][i]); //social coefficient phi_1
					//					if (i==0) cout << varPhi2 << " ";
				}
			}
		}
		varPhi2 = phi_2;
		//		if (i==0) cout << "]" << endl;
	}
}

// The computation of the radius and the random point in the HyperSphere
// was taken from the publicly available code of Maurice Clerc - Standard PSO 2011
// https://www.particleswarm.info/Programs.html
void Particle::getSphericalDNPP(double vect_distribution[], int numInformants, bool pBestIntheInformants,
		vector<vector< double> > &vect_PbestMinusPosition, int modelOfInflu){
	double V2[size], V1[size]; //working space arrays
	double G[size];	//center of the sphere
	double radius = 0.0;	//radius G-X
	double pw=1./(double)size;

	//Compute G (center of the sphere) and V1 (radius of each dimension)
	double varPhi2 = phi_2;
	for (int i=0; i<size; i++){
		//		if (i==0) cout << "\tPhi2 values: [ ";
		//		if (i==1) cout << "\tPhi2 values: [ ";
		double R = 0.0;
		for (int j=0; j<numInformants; j++){
			if (pBestIntheInformants && (this->id != this->gBestID)){ //we look for pBest in the Array
				if (this->id == neighbours.at(InformantsPos[j])->getID())
					R += current.x[i] + (phi_1 * vect_PbestMinusPosition[j][i]); //personal coefficient phi_1
				else {
					if (j==0 && modelOfInflu == MOI_RANKED_FI)
						varPhi2 = phi_2*((double)(numInformants-1)/numInformants);
					if (j>0 && modelOfInflu == MOI_RANKED_FI)
						varPhi2 = varPhi2/2.0;
					R += current.x[i] + (varPhi2 * vect_PbestMinusPosition[j][i]); //social coefficient phi_2
					//					if (i==0) cout << varPhi2 << " ";
					//					if (i==1) cout << varPhi2 << " ";
				}
			}
			else {
				if (j == numInformants-1) //pBest is at the end of the Array
					R += current.x[i] + (phi_1 * vect_PbestMinusPosition[j][i]); //personal coefficient phi_1
				else {
					if (j==0 && modelOfInflu == MOI_RANKED_FI)
						varPhi2 = phi_2*((double)(numInformants-1)/numInformants);
					if (j>0 && modelOfInflu == MOI_RANKED_FI)
						varPhi2 = varPhi2/2.0;
					R += current.x[i] + (varPhi2 * vect_PbestMinusPosition[j][i]); //social coefficient phi_2
					//					if (i==0) cout << varPhi2 << " ";
					//					if (i==1) cout << varPhi2 << " ";
				}
			}
		}
		G[i] = (current.x[i] + R )/(numInformants + 1);
		radius += pow(abs(current.x[i] - G[i]), 2);
		V1[i] = G[i] - current.x[i];
		varPhi2 = phi_2;
		//		if (i==0) cout << "]" << endl;
		//		if (i==1) cout << "]" << endl;
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

	//Return a vector from a hyperspherical distribution with center in G
	for (int i=0;i<size;i++) {
		vect_distribution[i] = V2[i] + V1[i];
	}
}

double Particle::applyPerturbation(int pertubType, double pos_xi){
	double returnVal = 0.0;
	switch(pertubType){
	case PERT1_NONE:
		returnVal = pos_xi;
		break;
	case PERT1_NORMAL_DISTANCE:
		//Gaussian distribution
		returnVal = RNG::randGaussWithMean(perturbMagnitud1, pos_xi);
		break;
	case PERT1_NORMAL_SUCCESS:
		//Gaussian distribution
		returnVal = RNG::randGaussWithMean(perturbMagnitud1, pos_xi);
		break;
	case PERT1_CAUCHY_DISTANCE:
		//Cauchy distribution
		returnVal = RNG::randCauchy(perturbMagnitud1);
		break;
	case PERT1_CAUCHY_SUCCESS:
		//Cauchy distribution
		returnVal = RNG::randCauchy(perturbMagnitud1);
		break;
	}
	//cout << "\tvar::returnVal: " << returnVal << endl;
	return returnVal;
}

//Perturbation 1 -- informants-wise
void Particle::setPerturbationMagnitude(int pertubType, double * pos_x, double * pbest_x, double alpha_t, double l){
	double distance = 0.0;
	switch(pertubType){
	case PERT1_NONE: //Do not apply perturbation
		perturbMagnitud1 = 1.0;
		break;
	case PERT1_NORMAL_DISTANCE:
		distance = computeDistance(pos_x, pbest_x);
		if (distance != 0.0)
			perturbMagnitud1 = (l*distance);
		break;
	case PERT1_CAUCHY_DISTANCE: //Normally distributed (here, we only compute the std. deviation
		distance = computeDistance(pos_x, pbest_x);
		if (distance != 0.0)
			perturbMagnitud1 = (l*distance);
		break;
	case PERT1_NORMAL_SUCCESS:
		perturbMagnitud1 = alpha_t*(1-(2*problem->getRandomX(0,1)));
		break;
	case PERT1_CAUCHY_SUCCESS: //Additional rectangular
		perturbMagnitud1 = alpha_t*(1-(2*problem->getRandomX(0,1)));
		break;
	}
}
//Perturbation 1 -- dimension-wise
void Particle::setPerturbationMagnitude(int pertubType, double alpha_t){
	switch(pertubType){
	case PERT1_NONE: //Do not apply perturbation
		perturbMagnitud1 = 1.0;
		break;
	case PERT1_NORMAL_DISTANCE || PERT1_CAUCHY_DISTANCE:
	break;
	case PERT1_NORMAL_SUCCESS:
		perturbMagnitud1 = alpha_t*(1-(2*problem->getRandomX(0,1)));
		break;
	case PERT1_CAUCHY_SUCCESS: //Additional rectangular
		perturbMagnitud1 = alpha_t*(1-(2*problem->getRandomX(0,1)));
		break;
	default:
		break;
	}
}
//Perturbation 2 (additive perturbation) -- particle-wise
double Particle::getPerturbationMagnitude(int pertubType, double alpha_t, double delta){
	double returnVal=0.0;
	switch(pertubType){
	case PERT2_NONE: //Do not apply perturbation
		returnVal = 0.0;
		perturbMagnitud2 = returnVal;
		break;
	case PERT2_ADD_RECT : //Additional rectangular
		returnVal = alpha_t*(1-(2*problem->getRandomX(0,1)));
		perturbMagnitud2 = returnVal;
		break;
	case PERT2_ADD_NOISY: //Additional noisy
		returnVal = problem->getRandomX(-delta/2,delta/2);
		perturbMagnitud2 = returnVal;
		break;
	default:
		returnVal = 0.0;
		perturbMagnitud2 = returnVal;
		break;
	}
	return returnVal;
}

//Multiply by the random matrix
void Particle::multiplyVectorByRndMatrix(vector<vector< double> > &vect_PbestMinusPosition, int informant,
		double ** rndMatrix[], int RmatrixType){
	double resultvxM[size]; //working space variable

	switch(RmatrixType){
	case MATRIX_NONE:
		for (int i=0; i<size; i++)
			resultvxM[i] = vect_PbestMinusPosition[informant][i];
		break;
	case MATRIX_DIAGONAL || MATRIX_LINEAR: //Random diagonal matrix
	for (int i=0; i<size; i++)
		resultvxM[i] = vect_PbestMinusPosition[informant][i]* rndMatrix[0][i][i];
	break;
	case MATRIX_RRM_EXP_MAP:{ //Random rotation matrix using exponential method
		for (int i = 0 ; i < size ; i ++) {
			resultvxM[i] = 0.0;
			for (int j = 0 ; j < size ; j ++) {
				resultvxM[i] += (vect_PbestMinusPosition[informant][i] * rndMatrix[0][i][j]);
			}
		}
	}
	break;
	case MATRIX_RRM_EUCLIDEAN_ONE:{ //Random rotation matrix using Euclidean rotation (ONLY ONE PLANE)
		bool exitFlag = false;
		//copy the vector
		for (int i=0; i<size; i++){
			resultvxM[i] = vect_PbestMinusPosition[informant][i];
		}
		//find the planes to rotate
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++){
				if (i != j && rndMatrix[0][i][j] != 0){
					resultvxM[i] = (vect_PbestMinusPosition[informant][i] * rndMatrix[0][i][i]) + (vect_PbestMinusPosition[informant][j] * rndMatrix[0][j][i]);
					resultvxM[j] = (vect_PbestMinusPosition[informant][j] * rndMatrix[0][j][j]) + (vect_PbestMinusPosition[informant][i] * rndMatrix[0][i][j]);
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
			resultvxM[i] = vect_PbestMinusPosition[informant][i];
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
							resultvxM[i] = (vect_PbestMinusPosition[informant][i] * rndMatrix[matNum][i][i]) +
									(vect_PbestMinusPosition[informant][j] * rndMatrix[matNum][j][i]);
							resultvxM[j] = (vect_PbestMinusPosition[informant][j] * rndMatrix[matNum][j][j]) +
									(vect_PbestMinusPosition[informant][i] * rndMatrix[matNum][i][j]);
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
		vect_PbestMinusPosition[informant][i] = resultvxM[i];
		//		cout << " " << aVector[i] << " ";
	}
	//	cout << "]" << endl;

	//return the aVector
	//return(vect_PbestMinusPosition[informant]);
}

//Compute the random matrix to employ
void Particle::computeRndMatrix(double ** rndMatrix[], int RmatrixType, double angle){
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
		//double angle = (problem->getRandomX(0.001,7)*PI)/180; //rotation between 0 and 10 degrees
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
				//double angle = (problem->getRandomX(0.001,7)*PI)/180; //rotation between 0 and 10 degrees
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

void Particle::evaluateSolution() {
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
void Particle::updateGbestParticle(double* x, double eval){
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

void Particle::setgBestID(int gB_ID){
	gBestID = gB_ID;
}

/*Check the neighborhood for the best particle */
int Particle::getBestOfNeibourhood(){
	//double aux_eval=gbest.eval;
	double aux_eval=LDBL_MAX;
	int best=-1;
	bool pBestIsNeighbor = false;

	//Check if the particle is its own neighbor
	for(unsigned int i=0; i<neighbours.size(); i++){
		if (neighbours.at(i)->getID() == id)
			pBestIsNeighbor = true;
	}
	//If the particles is its own neighbor check personal best
	if(pbest.eval < gbest.eval && pBestIsNeighbor == true){
		updateGbestParticle(pbest.x, pbest.eval);
		gBestID = id;
		//aux_eval = LDBL_MAX;
		aux_eval = pbest.eval;
	}
	//Check after the rest of the neighborhood
	for(unsigned int i=0; i<neighbours.size();i++){
		if(aux_eval > neighbours.at(i)->getPbestEvaluation()){
			best = i;
		}
	}
	//New best particle in the neighborhood
	if(best!=-1){
		updateGbestParticle(this->neighbours[best]->getPbestPosition(), this->neighbours[best]->getPbestEvaluation());
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
		//cout << "Couldn't find any :/" << endl;
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
