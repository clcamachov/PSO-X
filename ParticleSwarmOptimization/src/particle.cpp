/*
 * particle.cpp
 *
 *  Created on: May 31, 2019
 *      Author: Christian L. Camacho Villalón
 */


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
	problem = NULL;
	size =-1;
	id = -1;
	ranking = -1;
	parent =-1;
	stereotype = -1;
	lbestID = -1;
	velocity = NULL;
	minVelLimit = 0;					//minimum velocity -v_max = (x_max-x_min)/2
	maxVelLimit = 0;					//maximum velocity  v_max
	phi_1 = 0;
	phi_2 = 0;
	inertia = 0;

	init = false;
}

/* Constructor*/
Particle::Particle (Problem* problem, Configuration* config, int identifier, long int iteration){

	this->problem = problem;
	size= problem->getProblemDimension();
	id = identifier;
	ranking = 0;
	parent = 0;
	stereotype = 0;
	lbestID = -1;

	current.x = new double[size];
	pbest.x = new double[size];
	lbest.x = new double[size];
	current.eval = LDBL_MAX;
	pbest.eval = LDBL_MAX;
	lbest.eval = LDBL_MAX;

	velocity = new double[size];
	setVelocityLimits(config); //Set velocity limits
	inertia = config->getOmega1();
	phi_1 = config->getPhi1();
	phi_2 = config->getPhi2();

	initializePosition(config, iteration, true);

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
	lbestID = p.lbestID;

	if(!init){
		current.x = new double[size];
		pbest.x = new double[size];
		lbest.x = new double[size];
		velocity = new double[size];
	}
	for(int i=0; i<size;i++){
		current.x[i] = p.current.x[i];
		pbest.x[i] = p.pbest.x[i];
		lbest.x[i] = p.lbest.x[i];
		velocity[i] = p.velocity[i];
	}
	current.eval = p.current.eval;
	pbest.eval = p.pbest.eval;
	lbest.eval = p.lbest.eval;

	minVelLimit = p.minVelLimit;
	maxVelLimit = p.maxVelLimit;
	inertia=p.inertia;
	phi_1 = p.phi_1;
	phi_2 = p.phi_2;

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
		lbestID = p.lbestID;

		if(!init){
			current.x = new double[size];
			pbest.x = new double[size];
			lbest.x = new double[size];
			velocity = new double[size];
		}
		for(int i=0; i<size;i++){
			current.x[i] = p.current.x[i];
			pbest.x[i] = p.pbest.x[i];
			lbest.x[i] = p.lbest.x[i];
			velocity[i] = p.velocity[i];
		}
		current.eval = p.current.eval;
		pbest.eval = p.pbest.eval;
		lbest.eval = p.lbest.eval;

		minVelLimit = p.minVelLimit;
		maxVelLimit = p.maxVelLimit;
		inertia=p.inertia;
		phi_1 = p.phi_1;
		phi_2 = p.phi_2;

		init = true;
	}
	return (*this);
}

/* destructor */
Particle::~Particle(){
	if(init){
		delete[] current.x;
		delete[] pbest.x;
		delete[] lbest.x;
		delete[] velocity;
	}
	init=false;
}

void Particle::setRandomPositionInBoundsWithProbability(Configuration* config){
	for (int j=0; j<size; j++){
		if (RNG::randVal(0,1) < (1/size)){
			if(config->getCompetitionID() == CEC05 && config->getProblemID() == SHIFTED_ROTATED_GRIEWANK_CEC05)
				current.x[j] = RNG::randVal(0.0,600.0);
			else if(config->getCompetitionID() == CEC05 && config->getProblemID() == ROTATED_HYBRIDCOMPOSITION4_NO_BOUNDS)
				current.x[j] = RNG::randVal(2.0,5.0);
			else
				current.x[j] = RNG::randVal(config->getMinInitBound(),config->getMaxInitBound());
		}
	}
	//	initializeVelocity(config);
}

void Particle::setRandomPositiongBestDerivative(Configuration* config, double* global_bestPos){
	for (int j=0; j<size; j++)
		current.x[j] = (global_bestPos[j]-current.x[j])/2;
	//	initializeVelocity(config);
}

void Particle::initializePosition(Configuration* config, long int iteration, bool updatePbest){
	if(config->getParticleInitType() == PARTICLE_INIT_RANDOM){
		initUniform(config);
	}
	if(config->getParticleInitType() == PARTICLE_INIT_MODEL){ //Only available for dynamic PopCS and after the first iteration
		initUniform(config);
		if ((config->getPopulationCS() != POP_CONSTANT) && (iteration > 0))
			initToModel();
	}
	if(updatePbest){
		for(int i=0; i<size; i++)
			pbest.x[i]=current.x[i];
		evaluateSolution();
	}
}
/* Initialize particle using uniformly random values */
void Particle::initUniform(Configuration* config){
	for(int i=0; i<size; i++){
		if(config->getCompetitionID() == CEC05 && config->getProblemID() == SHIFTED_ROTATED_GRIEWANK_CEC05)
			current.x[i] = RNG::randVal(0.0,600.0);
		else if(config->getCompetitionID() == CEC05 && config->getProblemID() == ROTATED_HYBRIDCOMPOSITION4_NO_BOUNDS)
			current.x[i] = RNG::randVal(2.0,5.0);
		else
			current.x[i] = RNG::randVal(config->getMinInitBound(),config->getMaxInitBound()); //random values within the bounds of the function
		velocity[i]=0;
	}
}
/* Initialize particle using a model (as in Incremental PSO) */
void Particle::initToModel(){
	for (int i=0; i<size; i++){
		current.x[i] = current.x[i] + RNG::randVal(0,1)*(lbest.x[i]-current.x[i]);
	}
}
void Particle::initializeVelocity(Configuration* config){
	for(int i=0; i<size; i++)
		velocity[i]=RNG::randVal(minVelLimit,maxVelLimit);
}

void Particle::setVelocityLimits(Configuration* config){
	if (config->useVelocityClamping()){
		if(config->getCompetitionID() == CEC05 &&
				(config->getProblemID() == SHIFTED_ROTATED_GRIEWANK_CEC05 || config->getProblemID() == ROTATED_HYBRIDCOMPOSITION4_NO_BOUNDS)){
			maxVelLimit = numeric_limits<double>::max();
			minVelLimit = -numeric_limits<double>::max();
		}
		else {
			maxVelLimit = ((config->getMaxInitRange()-config->getMinInitRange())/2.0);
			minVelLimit = (-maxVelLimit);
		}
	}
	else {
		if(config->getCompetitionID() == CEC05 &&
				(config->getProblemID() == SHIFTED_ROTATED_GRIEWANK_CEC05 || config->getProblemID() == ROTATED_HYBRIDCOMPOSITION4_NO_BOUNDS)){
			maxVelLimit = (config->getMaxInitRange());
			minVelLimit = (config->getMinInitRange());
		}
		else {
			maxVelLimit = (config->getMaxInitRange()*2);
			minVelLimit = (config->getMinInitRange()*2);
		}
	}
}

double Particle::getMinVelLimit(){
	return (minVelLimit);
}

double Particle::getMaxVelLimit(){
	return (maxVelLimit);
}

int Particle::getRanking(){
	return (ranking);
}

void Particle::setPhi1(double new_phi_1){
	phi_1 = new_phi_1;
}
void Particle::setPhi2(double new_phi_2){
	phi_2 = new_phi_2;
}
double Particle::getPhi1(){
	return (phi_1);
}
double Particle::getPhi2(){
	return (phi_2);
}
void Particle::setRanking(int rank){
	ranking = rank;
}

int Particle::getParent(){
	return (parent);
}

void Particle::setParent(int node){
	parent = node;
}

bool Particle::ispBestIntheInformants(int numInformants){
	bool result = false;
	for (unsigned int j=0; j<informants.size(); j++)
		if (this->id == neighbors.at(informants[j])->getID() )
			result = true;
	return (result);
}

/* Generate a new solution by updating the particle's position */
void Particle::move(Configuration* config, long int iteration, double omega1, double omega2, double omega3,
		int lastLevelComplete, int solImproved){

	double vect_distribution[size];
	double vect_perturbation[size];

	for (int i=0; i<size; i++){
		vect_perturbation[i]=0;
		vect_distribution[i]=0;
	}

	if (config->getPerturbation2CS() != PERT2_NONE && config->getOmega3CS() != O3_ZERO){
		/*** PERTURBATION 2	(additive)	--->	It has to be computed per dimension, but only once per particle  ***/
		getRandomAdditivePerturbation(config, vect_perturbation); //Random additive perturbation
	}

	if (config->getOmega2CS() != O2_ZERO){
		/* vect_PbestMinusPosition is a 2D vector for computing:
		 *		DNPP-rectangular: Mtx^k (pert(p^k)-x^i)
		 *		DNPP-spherical: vector P^i and L^i without the multiplication by varphi
		 *		DNPP-additive stochastic p'^k and p'^i in
		 */
		vector<vector< double> > vect_PbestMinusPosition;
		vect_PbestMinusPosition.resize( informants.size(), vector<double>(size));

		/*** DISTRIBUTION VECTOR (DNNPs)***/
		if (config->getDistributionNPP() == DIST_RECTANGULAR) {
			computeSubtractionPerturbationRotation(
					config,
					vect_PbestMinusPosition,
					iteration,
					solImproved);
			getRectangularDNPP(
					config,
					vect_distribution,
					vect_PbestMinusPosition );
		}
		if (config->getDistributionNPP() == DIST_SPHERICAL){
			computeSubtractionPerturbationRotation(
					config,
					vect_PbestMinusPosition,
					iteration,
					solImproved);
			getSphericalDNPP(
					config,
					vect_distribution,
					vect_PbestMinusPosition);
		}
		if (config->getDistributionNPP() == DIST_ADD_STOCH){
			computeSubtractionPerturbationRotation(
					config,
					vect_PbestMinusPosition,
					iteration,
					solImproved);
			getAdditiveStochasticDNPP(
					config,
					vect_distribution,
					vect_PbestMinusPosition);
		}
	}
	//Compute new position
	for (int i=0;i<size;i++) {
		velocity[i] = (omega1 *  velocity[i]) +
				(omega2 * vect_distribution[i]) +
				(omega3 * vect_perturbation[i]);

		//Clamp velocity
		if (config->useVelocityClamping()){
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
		if(current.x[i] < config->getMinInitBound())
			current.x[i] = config->getMinInitBound();
		if(current.x[i] > config->getMaxInitBound())
			current.x[i]= config->getMaxInitBound();
	}


	//Detect stagnation and introduce perturbation using Schmitt and Wanka technique
	if(config->detectParticleStagnated())
		detectStagnation(config);

	//Evaluate the objective function and update pbest if a new one has been found
	evaluateSolution();

	if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_COMPUTATIONS) ){
		cout << "\tvec::p[" << id << "].x(new): [";
		for(int i=0; i<size; i++){
			cout << fixed << current.x[i] << " ";
		} cout << "]" << endl;
	}

	if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE)) cout << "\n\tParticle with ID:[" << id << "].status()::MOVED" << endl;
	if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_VARIABLE)) cout << "\t------------------------------------------" << endl;
}

void Particle::detectStagnation(Configuration* config){
	//Compute norm of the velocity
	double vel_norm =0;
	for (int i=0;i<size;i++) {
		vel_norm += pow(velocity[i],2);
	}
	vel_norm = sqrt(vel_norm);

	double gBestpPos_norm = 0;
	for (int i=0;i<size;i++) {
		gBestpPos_norm += pow((lbest.x[i]-current.x[i]),2);
	}
	gBestpPos_norm = sqrt(gBestpPos_norm);

	if((vel_norm + gBestpPos_norm) < STAGNATION_PRECISION){
		//Compute new position
		for (int i=0;i<size;i++) {
			velocity[i] = ((2*RNG::randVal(0,1))-1)*STAGNATION_PRECISION;

			//Clamp velocity
			if (config->useVelocityClamping()){
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
			if(current.x[i] < config->getMinInitBound())
				current.x[i] = config->getMinInitBound();
			if(current.x[i] > config->getMaxInitBound())
				current.x[i]= config->getMaxInitBound();
		}
	}
}

void Particle::computeSubtractionPerturbationRotation(
		Configuration* config,
		vector<vector< double> > &vect_PbestMinusPosition, //This is the vector we aim to compute here
		long int iteration,
		int solImprov) {

	double pertMagnitude[size];		//vector containing the magnitude of the perturbation that will be applied to each dimension


	if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_COMPUTATIONS)){
		cout << "\n\tvec::p[" << id << "].x:\t  [";
		for(int i=0; i<size; i++){
			cout << fixed << current.x[i] << " ";
		} cout << "]" << endl;
	}
	//Get the p^k-x^i of all Informants
	for (unsigned int j=0; j<informants.size(); j++){

		if (config->verboseMode() && (config->verboseLevel() >= VERBOSE_LEVEL_COMPUTATIONS)){
			cout << "\tvec::inf.p[" << neighbors[informants[j]]->getID() << "].pb: [";
			for(int i=0;i<size;i++){
				cout << fixed << neighbors[informants[j]]->pbest.x[i] << " ";
			}
			cout << "]" << endl;
		}
		//informant-wise perturbation magnitude
		setPerturbation1Magnitude(config, pertMagnitude, current.x, neighbors.at(informants[j])->pbest.x);

		for (int i=0; i<size; i++){
			vect_PbestMinusPosition.at(j).at(i) = applyInformedPerturbation(config, pertMagnitude, neighbors.at(informants[j])->pbest.x[i], i, iteration);
			if (config->getDistributionNPP() != DIST_ADD_STOCH)
				vect_PbestMinusPosition.at(j).at(i) = vect_PbestMinusPosition.at(j).at(i) -current.x[i];
		}
	}
	if (config->getDistributionNPP() != DIST_ADD_STOCH) {
		//3.- Rotate each vector
		double ** rndMatrix = new double*[size];

		//Allocate memory to be used
		for(int i=0; i<size; i++){
			rndMatrix[i] = new double[size];
		}
		//Create the matrix and multiply vect_PbestMinusPosition by the matrix
		for (unsigned int j=0; j<informants.size(); j++){ //the rest of informants
			computeRndMatrix(config, rndMatrix, config->getRandomMatrix(), getAnAngle(config, solImprov, iteration));
			multiplyVectorByRndMatrix(config, vect_PbestMinusPosition, j, rndMatrix, config->getRandomMatrix(), solImprov, iteration);
		}
		//Deallocate memory used
		for(int i=0; i<size; i++)
			delete [] rndMatrix[i];
		delete [] rndMatrix;
	}
}

//This function returns a random position in informants that is different from the pBest of the particle
int Particle::getRandomInformantPosition(){
	int randomIndex = 0;

	//The size of informants is at least 1
	for (unsigned int i=0; i<informants.size(); i++) {
		randomIndex = (int)floor(RNG::randVal(0.0,(double)informants.size()));

		//Distinct of itself
		if (neighbors.at(informants[randomIndex])->getID() != id){
			break;
		}
		else
			continue;
	}
	return (randomIndex);
}

//This function returns the position of the particle's pBest in informants
int Particle::getPositionOfpBest(){
	int pBestIndex = 0;

	//The size of informants is at 1
	for (unsigned int i=0; i<informants.size(); i++) {
		if (neighbors.at(informants[i])->getID() == id){
			pBestIndex = i;
			break;
		}
		else
			continue;
	}
	return (pBestIndex);
}

void Particle::getAdditiveStochasticDNPP(Configuration* config, double vect_distribution[], vector<vector< double> > &vect_PbestMinusPosition){

	int p2Index;
	int p1Index = getPositionOfpBest();

	if (config->getRandNeighbor())
		p2Index = getRandomInformantPosition();
	else
		p2Index = 0;

	switch (config->getOperator_q()) {
	//Discount current.x[i] in all these operators; it will be added later in the GVU formula
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
			double center = (vect_PbestMinusPosition[p1Index][i] + vect_PbestMinusPosition[p2Index][i]) / 2.0;
			double dispersion = fabs (vect_PbestMinusPosition[p1Index][i] - vect_PbestMinusPosition[p2Index][i]);
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
			if (RNG::randVal(0,1) <= config->getOperatorCG_parm_r()){
				vect_distribution[i] = vect_PbestMinusPosition[p1Index][i] - current.x[i];
				vect_distribution[i] += (RNG::randCauchy(1.0) * fabs(vect_PbestMinusPosition[p1Index][i] -
						vect_PbestMinusPosition[p2Index][i]));
			}
			else {
				vect_distribution[i] = vect_PbestMinusPosition[p2Index][i] - current.x[i];
				vect_distribution[i] += (RNG::randGaussWithMean(1.0,0) * fabs(vect_PbestMinusPosition[p1Index][i] -
						vect_PbestMinusPosition[p2Index][i]));
			}
		}
		break;
	}
}

void Particle::computeAC(Configuration* config, double &c1, double &c2){

	//When the MoI is BoN there is not need to make any change in the values
	if (config->getModelOfInfluence() == MOI_BEST_OF_N){
		c1=phi_1;
		c2=phi_2;
	}
	if (config->getModelOfInfluence() == MOI_FI) {
		if (config->getAccelCoeffCS() == AC_CONSTANT){
			double varphi = phi_1 + phi_2;
			c1 = varphi/informants.size();
			c2 = c1;
		}
		else{
			c1 = phi_1;
			c2 = phi_2/(informants.size()-1);
		}
	}
	if (config->getModelOfInfluence() == MOI_RANKED_FI) {
		c1=phi_1;
		c2=phi_2;
	}
}

void Particle::getRectangularDNPP(Configuration* config, double vect_distribution[],
		vector<vector< double> > &vect_PbestMinusPosition){

	double varPhi1=0;
	double varPhi2=0;

	//Compute vect_distribution
	for (int i=0; i<size; i++){
		//Compute the value of varphi1 and varph2 according to the model of influence
		computeAC(config, varPhi1, varPhi2);

		for (unsigned int j=0; j<informants.size(); j++){
			if (id == neighbors.at(informants[j])->getID())
				vect_distribution[i] += (varPhi1 * vect_PbestMinusPosition[j][i]); //personal coefficient phi_1
			else {
				if (j==0 && config->getModelOfInfluence() == MOI_RANKED_FI)
					varPhi2 = phi_2*((double)(informants.size()-1)/informants.size());
				if (j>0 && config->getModelOfInfluence() == MOI_RANKED_FI)
					varPhi2 = varPhi2/2;
				vect_distribution[i] += (varPhi2 * vect_PbestMinusPosition[j][i]); //social coefficient phi_2
			}
		}
	}
}

// The computation of the radius and the random point in the HyperSphere
// was taken from the publicly available code of Maurice Clerc - Standard PSO 2011
// https://www.particleswarm.info/Programs.html
void Particle::getSphericalDNPP(Configuration* config, double vect_distribution[], vector<vector< double> > &vect_PbestMinusPosition){
	double V2[size];
	double V1[size]; 		//working space arrays
	double G[size];			//center of the sphere
	double radius = 0.0;	//radius G-X
	double pw=1./(double)size;

	//Compute G (center of the sphere) and V1 (radius of each dimension)
	double varPhi1=0;
	double varPhi2=0;

	for (int i=0; i<size; i++){
		//Compute the value of varphi1 and varph2 according to the model of influence
		computeAC(config, varPhi1, varPhi2);

		double R = 0.0;
		for (unsigned int j=0; j<informants.size(); j++){
			if (id == neighbors.at(informants[j])->getID())
				R += current.x[i] + (phi_1 * vect_PbestMinusPosition[j][i]); //personal coefficient phi_1
			else {
				if (j==0 && config->getModelOfInfluence() == MOI_RANKED_FI)
					varPhi2 = phi_2*((double)(informants.size()-1)/informants.size());
				if (j>0 && config->getModelOfInfluence() == MOI_RANKED_FI)
					varPhi2 = varPhi2/2.0;
				R += (varPhi2 * vect_PbestMinusPosition[j][i]); //social coefficient phi_2
			}
		}
		G[i] = (current.x[i] + current.x[i] + R )/(3);
		V1[i] = G[i] - current.x[i];

		radius += pow(fabs(current.x[i] - G[i]), 2);
		varPhi2 = phi_2;
	}
	radius = pow(radius, 1/2); //this is the actual radius of the hyper-sphere

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

double Particle::applyInformedPerturbation(Configuration* config, double pertMagnitude[], double pos_xi, int index, long int iteration){
	double returnVal = 0.0;

	if(config->getPerturbation1CS() == PERT1_GAUSSIAN){
		//Gaussian distribution
		iteration == 1 ? returnVal = pos_xi + RNG::randGauss(1) :
				returnVal = pos_xi + RNG::randGauss(1)*pertMagnitude[index];
	}
	else if(config->getPerturbation1CS() == PERT1_LEVY){

		double alpha = floor(RNG::randVal(10,20))/10;
		//alpha = 1 = Cauchy distribution
		//alpha = 2 = Gaussian distribution

		//Adjust the peak (parameter c in the function RNG::randLevy()) according to the values provided in Blackwell 2006
		//The peak of the distribution is obtained by multiplying c*σ and σ is already in the variable pert_vrand[index]
		if (alpha <= 1.2)
			pertMagnitude[index] = pertMagnitude[index] * 0.557;
		if (alpha > 1.2 && alpha <= 1.3)
			pertMagnitude[index] = pertMagnitude[index] * 0.562;
		if (alpha > 1.3 && alpha <= 1.4)
			pertMagnitude[index] = pertMagnitude[index] * 0.565;
		if (alpha > 1.4 && alpha <= 1.5)
			pertMagnitude[index] = pertMagnitude[index] * 0.568;
		if (alpha > 1.5 && alpha <= 1.6)
			pertMagnitude[index] = pertMagnitude[index] * 0.572;
		if (alpha > 1.6 && alpha <= 1.7)
			pertMagnitude[index] = pertMagnitude[index] * 0.576;
		if (alpha > 1.7 && alpha <= 1.8)
			pertMagnitude[index] = pertMagnitude[index] * 0.580;
		if (alpha > 1.8 && alpha <= 1.9)
			pertMagnitude[index] = pertMagnitude[index] * 0.585;
		if (alpha > 1.9 && alpha <= 2.0)
			pertMagnitude[index] = pertMagnitude[index] * 0.633;
		//Levy distribution
		returnVal = pos_xi + (RNG::randLevy(pertMagnitude[index],alpha));
	}
	else if(config->getPerturbation1CS() == PERT1_UNIFORM){
		//Cauchy distribution
		returnVal = pos_xi + (RNG::randVal(-1,1)*pertMagnitude[index]*pos_xi);
	}
	else
		returnVal = pos_xi;

	return (returnVal);
}

//Perturbation 1
void Particle::setPerturbation1Magnitude(Configuration* config, double pertMagnitude[], double * pos_x, double * pbest_x){

	switch (config->getMagnitude1CS()) {
	case MAGNITUDE_EUC_DISTANCE:{
		double distance = 0.0;
		//Compute the norm of pos_x - pbest_x
		for(int i=0;i<size;i++){
			distance += pow(pos_x[i]-pbest_x[i],2);
		}
		distance = sqrt(distance);

		if (distance == 0){ //use the last computed value
			for(int i=0;i<size;i++)
				pertMagnitude[i] = config->getMagnitude1();
		}
		else{
			for(int i=0;i<size;i++){
				pertMagnitude[i] = config->getMag1_parm_l()*distance;
			}
			config->setMagnitude1(pertMagnitude[0]);
		}
	}
	break;
	case MAGNITUDE_OBJ_F_DISTANCE:{
		double ob_distance = 0.0;

		//Compute the distance in terms of the objective function
		ob_distance = (lbest.eval-current.eval)/lbest.eval;

		if (ob_distance == 0){ //use the last computed value
			for(int i=0;i<size;i++)
				pertMagnitude[i] = config->getMagnitude1();
		}
		else{
			for(int i=0;i<size;i++)
				pertMagnitude[i] = config->getMag1_parm_m()*ob_distance;
		}
		config->setMagnitude1(pertMagnitude[0]);
	}
	break;
	case MAGNITUDE_SUCCESS:{
		for(int i=0;i<size;i++){
			pertMagnitude[i] = config->getMagnitude1();
		}
	}
	break;
	case MAGNITUDE_CONSTANT:{
		for(int i=0;i<size;i++){
			pertMagnitude[i] = config->getMagnitude1();
		}
	}
	break;
	}
}

//Perturbation 2 (additive perturbation) -- particle-wise
void Particle::getRandomAdditivePerturbation(Configuration* config, double vect_perturbation[]){
	double PM = 0;
	//1.0 - Set the perturbation magnitude
	switch (config->getMagnitude2CS()) {
	case MAGNITUDE_EUC_DISTANCE:{
		double distance = 0.0;
		int bestID = 0;

		//Get the position of gBestID in the neighbours array
		for(long unsigned int i=0;i<neighbors.size();i++){
			if (neighbors.at(i)->getID() == lbestID)
				bestID = i;
		}
		//Compute the norm of pos_x - lbest_x
		for(int i=0;i<size;i++){
			distance += pow((current.x[i] - neighbors.at(bestID)->pbest.x[i]),2);
		}
		distance = sqrt(distance);

		if (distance == 0){ //use the last computed value
			PM = config->getMagnitude2();
		}
		else
			PM = config->getMag2_parm_l()*distance;
		config->setMagnitude2(PM);
	}
	break;
	case MAGNITUDE_OBJ_F_DISTANCE:{
		double ob_distance = 0;

		//Compute the distance in terms of the objective function
		ob_distance = (lbest.eval-current.eval)/lbest.eval;

		if (ob_distance == 0)
			PM = config->getMagnitude2();
		else
			PM = config->getMag2_parm_m()*ob_distance;

		config->setMagnitude2(PM);
	}
	break;
	case MAGNITUDE_SUCCESS:
		PM = config->getMagnitude2();

		break;
	case MAGNITUDE_CONSTANT:
		PM = config->getMagnitude2();

		break;
	default:
		PM = 1.0;
		break;
	}

	//2.0 Compute the random vector
	if(config->getPerturbation2CS() == PERT2_RECTANGULAR){ //Additional rectangular
		for (int i=0; i<size; i++){
			vect_perturbation[i] = PM*(1-(2*RNG::randVal(0,1)));
		}
	}
	else if(config->getPerturbation2CS() == PERT2_NOISY){ //Additional noisy
		for (int i=0; i<size; i++){
			vect_perturbation[i] = RNG::randVal(-PM/2,PM/2);
		}
	}
	else {
		for (int i=0; i<size; i++){
			vect_perturbation[i] = 0;
		}
	}
}

//Multiply by the random matrix
void Particle::multiplyVectorByRndMatrix(Configuration* config, vector<vector< double> > &vect_PbestMinusPosition, int informant,
		double ** rndMatrix, int RmatrixType, int solImprov, long int iteration){
	double resultvxM[size]; //working space variable

	if (RmatrixType == MATRIX_NONE){
		for (int i=0; i<size; i++)
			resultvxM[i] = vect_PbestMinusPosition[informant][i];
	}
	if (RmatrixType == MATRIX_DIAGONAL || RmatrixType == MATRIX_LINEAR){
		for (int i=0; i<size; i++)
			resultvxM[i] = vect_PbestMinusPosition[informant][i] * rndMatrix[i][i];
	}
	if (RmatrixType == MATRIX_RRM_EXP_MAP){
		for (int i=0; i<size; i++){
			resultvxM[i] = 0.0;
			for (int j = 0 ; j < size ; j ++) {
				resultvxM[i] += (vect_PbestMinusPosition[informant][j] * rndMatrix[i][j]);
			}
		}
	}
	if (RmatrixType == MATRIX_RRM_EUCLIDEAN_ONE){
		//Copy the vector
		for (int i=0; i<size; i++){
			resultvxM[i] = vect_PbestMinusPosition[informant][i];
		}
		//Get two random planes
		int plane1 = (int)floor(RNG::randVal(0.0,(double)size-1));
		int plane2 = 0;
		for (int i=0; i<size; i++) {
			plane2 = (int)floor(RNG::randVal(0.0,(double)size-1));
			if (plane2 == plane1)
				continue;
			else
				break;
		}
		//compute the angle using any of the available strategies
		double angle = getAnAngle(config, solImprov, iteration);
		//Rotate the plane1 of the vector in direction of plane2
		//v_i = v_i * cos(α_k) - v_j * sin(α_k)
		resultvxM[plane1] = (vect_PbestMinusPosition[informant][plane1] * cos(angle)) - (vect_PbestMinusPosition[informant][plane2] * sin(angle));
		//v_j = v_j * cos(α_k) + v_i * sin(α_k)
		resultvxM[plane2] = (vect_PbestMinusPosition[informant][plane1] * sin(angle)) + (vect_PbestMinusPosition[informant][plane2] * cos(angle));
	}
	if (RmatrixType == MATRIX_RRM_EUCLIDEAN_ALL){
		//Copy the vector
		for (int i=0; i<size; i++){
			resultvxM[i] = vect_PbestMinusPosition[informant][i];
		}
		//Rotate the vector in all possible combination of planes
		for (int i=size-2; i>=0; i--){
			for (int j=size-1; j>i+1; j--){

				//compute the angle using any of the available strategies
				double angle = getAnAngle(config, solImprov, iteration);

				//v_i = v_i * cos(α_k) - v_j * sin(α_k)
				resultvxM[i] = (vect_PbestMinusPosition[informant][i] * cos(angle)) - (vect_PbestMinusPosition[informant][j] * sin(angle));
				//v_j = v_j * cos(α_k) + v_i * sin(α_k)
				resultvxM[j] = (vect_PbestMinusPosition[informant][i] * sin(angle)) + (vect_PbestMinusPosition[informant][j] * cos(angle));
			}
		}
	}

	//Copy the perturbed vector in the distribution vector
	for (int i=0; i<size; i++){
		vect_PbestMinusPosition[informant][i] = resultvxM[i];
	}
}

//Compute the random matrix to employ
void Particle::computeRndMatrix(Configuration* config, double ** rndMatrix, int RmatrixType, double angle){

	if (RmatrixType == MATRIX_DIAGONAL){//Random diagonal matrix
		for (int i=0; i<size; i++)
			rndMatrix[i][i] = RNG::randVal(0,1);
	}
	if (RmatrixType == MATRIX_LINEAR){ //Random linear matrix
		double rndVal = RNG::randVal(0,1);
		for (int i=0; i<size; i++)
			rndMatrix[i][i] = rndVal;
	}
	if (RmatrixType == MATRIX_RRM_EXP_MAP){ //Random rotation matrix using exponential method
		//1.- Generate a random matrix
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++)
				rndMatrix[i][j] = RNG::randVal(-0.5,0.5);
		}
		//2.- Generate the transpose of the rndMatrix
		double trans_rndMatrix[size][size];
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++)
				trans_rndMatrix[j][i]=rndMatrix[i][j];
		}
		//3.- Determine the rotation angle
		double expMapAngle = angle; //rotation between 0 and 10 degrees
		//4.- Subtract trans_rndMatrix to rndMatrix, multiply by the angle and add the identity matrix
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++){
				rndMatrix[i][j]=(rndMatrix[i][j]-trans_rndMatrix[i][j])*expMapAngle;
				if (i==j)
					rndMatrix[i][j]+=1.0; //Add the Identity matrix
			}
		}
	}
	if (RmatrixType == MATRIX_RRM_EUCLIDEAN_ONE){ //Random rotation matrix using Euclidean rotation (ONLY ONE PLANE)
		//No need to generate a matrix, the multiplication of the vector by the rotation matrix is done directly in multiplyVectorByRndMatrix()
		for (int i=0; i<size; i++)
			rndMatrix[i][i] = 1.0;
	}
	if (RmatrixType == MATRIX_RRM_EUCLIDEAN_ALL){ //Random rotation matrix using Euclidean rotation (ALL PLANES)
		//No need to generate a matrix, the multiplication of the vector by the rotation matrix is done directly in multiplyVectorByRndMatrix()
		for (int i=0; i<size; i++)
			rndMatrix[i][i] = 1.0;
	}
	if (RmatrixType == MATRIX_NONE){
		for (int i=0; i<size; i++)
			rndMatrix[i][i] = 1.0;
	}
	if (RmatrixType == MATRIX_GROUPED_INCREASING){//Random diagonal matrix using the increasing stochastic scaling grouping
		double * theMatrixGI;
		theMatrixGI = new double[size];
		theMatrixGI = config->getMatrixGI();
		for (int i=0; i<size; i++)
			rndMatrix[i][i] = theMatrixGI[i];
	}
}

//This function returns an angle in radian
double Particle::getAnAngle(Configuration* config, int solImprov, long int iteration){
	double angle = 0;

	if (config->getAngleCS() == ANGLE_CONSTANT) //use the same rotation angle in every iteration
		angle = (config->getRotationAgle()*PI)/180;

	if (config->getAngleCS() ==  ANGLE_NORMAL) //draw the angle from a Normal distribution with µ=0 and s.d. = angleSD
		angle = (RNG::randGauss(config->getAngleSD())*PI)/180;

	if (config->getAngleCS() ==  ANGLE_ADAPTIVE){ //draw the angle from a Normal distribution with µ=0 using an adaptive s.d.
		double alpha = config->get_angle_par_alpha();
		double beta = config->get_angle_par_beta();
		//use simpSwarm: the strategy is based on success.
		if (iteration == 1) { //define simpSwarm if it does not exist
			angle =  (RNG::randGauss(beta)*PI)/180;
		}
		else {
			double std_dev = ((alpha * ((double)solImprov/config->getSwarmSize()))/sqrt((double)config->getProblemDimension()))+beta;
			angle =  (RNG::randGauss(std_dev)*PI)/180;
		}
	}
	return (angle);
}


void Particle::evaluateSolution() {
	current.eval = problem->getFunctionValue(current.x);
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
void Particle::updatelBestParticle(double* x, double eval){
	for(int j=0; j<size; j++)
		lbest.x[j]= x[j];
	lbest.eval = eval;
}

double Particle::computeDistPbestGbest(){
	double sumDistSqr = 0.0;
	for(int i=0;i<size;i++){
		sumDistSqr += pow((this->lbest.x[i]-this->pbest.x[i]),2);
	}
	return (sqrt(sumDistSqr));
}

double* Particle::getCurrentPosition() {
	return(current.x);
}

long double Particle::getCurrentEvaluation(){
	return(current.eval);
}

double* Particle::getCurrentVelocity(){
	return (velocity);
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
	cout << "Solution: " << neighbors.at(identifier)->current.eval << " Id: " << neighbors.at(identifier)->id << endl;
	for(int i=0; i<size; i++){
		cout << neighbors.at(identifier)->current.x[i] << "  ";
	}
	cout << endl;
}

void Particle::setlBestID(int gB_ID){
	lbestID = gB_ID;
}
int Particle::getlBestID(){
	return (lbestID);
}

/*Check the neighborhood for the best particle */
int Particle::getBestOfNeibourhood(){
	double aux_eval=neighbors.at(0)->getPbestEvaluation();
	int pos = 0;

	//Look for best particle in the neighborhood
	for(unsigned int i=0; i<neighbors.size();i++){
		if(neighbors.at(i)->getPbestEvaluation() < aux_eval){
			aux_eval = neighbors.at(i)->getPbestEvaluation();
			pos = i;
		}
	}
	//New best particle in the neighborhood
	if(neighbors.at(pos)->getID() != lbestID ){
		updatelBestParticle(neighbors.at(pos)->getPbestPosition(), neighbors.at(pos)->getPbestEvaluation());
		lbestID = neighbors.at(pos)->getID();
		return (neighbors.at(pos)->getID());
	}
	else
		return (lbestID);
}

void Particle::addNeighbour(Particle* p){
	neighbors.push_back(p);
}

//From Frankesntein's PSO
unsigned int Particle::getNeighborhoodSize(){
	return (neighbors.size());
}

int Particle::getID(){
	return (id);
}

void Particle::setID(int newID){
	id = newID;
}

void Particle::eraseNeighborbyID(int nid){
	for(unsigned int i=0;i<neighbors.size();i++){
		if(nid == neighbors.at(i)->id){
			neighbors.erase(neighbors.begin()+i);
			break;
		}
	}
}

int Particle::getRandomNonAdjacentNeighborID(Configuration* config){
	int randomIndex;
	int randomID;
	if(neighbors.size()>3){
		while(true){
			randomIndex = (int)floor(RNG::randVal(0.0,(double)neighbors.size()));
			//One random neighbor
			randomID = neighbors.at(randomIndex)->id;

			//First position of the vector
			if(id==0){
				//Distinct of itself, the right-side neighbor and the last neighbor
				if(randomID != id && randomID != id+1 && randomID != config->getSwarmSize()-1){
					return (randomID);
				}
				else
					continue;
			}
			else{
				//Last position of the vector
				if(id == config->getSwarmSize()-1){
					//Distinct of itself, the first neighbor and the left-side neighbor
					if(randomID != id && randomID != 0 && randomID != id-1){
						return (randomID);
					}
					else
						continue;
				}
				else{
					//Some position in within the bounds just check right- and left-side neighbors
					if(randomID != id && randomID != id+1 && randomID != id-1){
						return (randomID);
					}
					else
						continue;
				}
			}
		}
	}
	else
		return (-1);
}
