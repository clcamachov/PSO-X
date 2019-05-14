/*
 * swarm.cpp
 *
 *  Created on: Jun 20, 2018
 *      Author: leonardo
 */
#include "iostream"
#include "float.h"
#include "utils.h"

#include "problem.h"
#include "config.h"
#include "swarm.h"

#include "rng.h"

using namespace std;

/* Variables to compute the inertia weight using the available control strategies */
double alpha = 1/pow(M_PI,2); 				//small positive constant 		IW_NONL_DEC - 4
double omega = 0.3;							//value between [0,1] 			IW_NONL_DEC_IMP - 5
double u = 1.0002;							//value between [1.0001,1.0005] IW_NONL_DEC_IMP - 5
double omegaXu = omega * u;					//								IW_NONL_DEC_IMP - 5
double zFunction;							//								IW_CHAOTIC_DEC - 7
int	   k = 7;								//positive integer constant		IW_OSCILLATING - 9
double simNumOfCos = 2.0*M_PI*((4*k)+6);	//								IW_OSCILLATING - 9
double a = 1;								//small positive constant		IW_LOG_DEC - 10
double omega_2 = 0;							//								IW_SELF_REGULATING - 11
double eta = 1;								//								IW_SELF_REGULATING - 11
double idealVelocity;						//								IW_VELOCITY_BASED - 12
double avVel;								//								IW_VELOCITY_BASED - 12
double deltaOmega = 0.1;					//small positive constant		IW_VELOCITY_BASED - 12
struct SimplifySwarm simpSwarm;				//								IW_RANKS_BASED - 14
double alpha_2 = 0.5;						//small constant in [0,1]		IW_CONVERGE_BASED
double beta_2 = 0.5;						//small constant in [0,1]		IW_CONVERGE_BASED


//Default constructor
Swarm::Swarm () {
	problem = 0;
	size = 0;
	best_particle = 0;
	ranked = false;
}

Swarm::~Swarm(){
	// Memory allocated dynamically
	for(long int i=0;i<size;i++)
		delete swarm.at(i);

	delete [] global_best.x;

	if (ranked){
		//free memory reserved to the use rankings
		delete [] simpSwarm.id;
		delete [] simpSwarm.eval;
	}
}

Swarm::Swarm (Problem* problem, Configuration* config){
	//cout << "Creating swarm.\n\n";

	this->problem = problem;
	size = config->getSwarmSize();

	/*Initialize global best*/
	global_best.x = new double[config->getProblemDimension()];
	for(unsigned int i=0;i<config->getProblemDimension();i++)
		global_best.x[i] = 0;
	global_best.eval=LDBL_MAX;

	for (long int i=0; i<size; i++) {
		Particle * aParticle = new Particle(problem, config, i);
		swarm.push_back(aParticle);

		if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
			updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
			best_particle = swarm.at(i);
		}
	}

	//Select one of the available topologies
	if (config->getTopology() == TOP_FULLYCONNECTED) {
		createFullyConnectedTopology();
	} else if (config->getTopology() == TOP_RING) {
		createRingTopology();
	} else if (config->getTopology() == TOP_STAR) {
		createStarTopology();
	} else if (config->getTopology() == TOP_RANDOM) {
		createRandomTopology();
	} else if (config->getTopology() == TOP_TIMEVARYING) {
		createFullyConnectedTopology();
		config->setTopologyUpdatePeriod(
				(int)floor((double)config->getTopologySchedule()/(config->getSwarmSize()-3)));
		RNG::initializePermutation(config->getSwarmSize());
	} else if (config->getTopology() == TOP_VONNEUMANN) {
		createVonNeumannTopology();
	}
	else {
		cerr << "Wrong topology" << endl;
		exit (-1);
	}

	//rankings
	if (config->getinertiaCS() == IW_RANKS_BASED || config->getinertiaCS() == IW_SUCCESS_BASED
			|| config->getinertiaCS() == IW_CONVERGE_BASED)
		ranked = true;

	//cout << "\nInitial best value: " << global_best.eval << "\n" << endl;

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
	//	for(unsigned int i=0; i< sizeof(getGlobalBest().x); i++){
	for(unsigned int i=0; i< dimensions; i++){
		cout << getGlobalBest().x[i] << "  ";
	}
	cout << " ]\n" << endl;
}


/*Move the swarm to new solutions */
void Swarm::moveSwarm(Configuration* config, long int iteration, const double minBound, const double maxBound) {
	//best value at iteration t-1
	double prev_eval=global_best.eval;

	//compute the inertia weight according to the strategy chosen
	computeInertia(config, iteration);

	cout << "iteration: " << iteration << endl; //remove

	//with this strategy the global_best particle has a different inertia weight than the rest of the swarm
	if (config->getinertiaCS() == IW_SELF_REGULATING){

		int best_particle_id = best_particle->getID(); //id of the best particle
		double inertia_bck = config->getInertia(); //backup the inertia weight

		//compute the inertia weight for the global_best particle
		config->setInertia(omega_2 + eta * ((config->getFinalIW() - config->getInitialIW()) / config->getMaxIterations()));
		swarm.at(best_particle->getID())->move(config, minBound,maxBound,iteration);

		//See if the best particle has a new better position
		if (swarm.at(best_particle_id)->getPbestEvaluation() < global_best.eval){
			updateGlobalBest(swarm.at(best_particle_id)->getPbestPosition(), swarm.at(best_particle_id)->getPbestEvaluation());
			best_particle = swarm[best_particle_id];
		}
		//restore inertia weight for the rest of the particles
		config->setInertia(inertia_bck);

		//Update the rest of the particles
		for (unsigned int i=0;i<swarm.size();i++){
			if (swarm.at(i)->getID() != best_particle_id) {

				swarm.at(i)->move(config, minBound,maxBound,iteration);

				if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
					updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
					best_particle = swarm[i];
				}
			}
		}
	}
	//with this strategy all particles has a different inertia weight
	else if (config->getinertiaCS() == IW_CONVERGE_BASED){
		if (iteration > 1) {
			for (unsigned int i=0; i<swarm.size(); i++){
				for (unsigned int j=0; j<sizeof(simpSwarm.id); j++){
					if (swarm.at(i)->getID() == simpSwarm.id[j]){
						//convergence factor
						long double c_i = abs(simpSwarm.eval[j] - swarm.at(i)->getCurrentEvaluation())/
								(simpSwarm.eval[j] + swarm.at(i)->getCurrentEvaluation());
						//diffusion factor
						long double d_i = abs(swarm.at(i)->getCurrentEvaluation() - prev_eval) /
								(swarm.at(i)->getCurrentEvaluation() + prev_eval);
						//set the inertia weight
						config->setInertia(1 - abs(alpha_2*(1-c_i)) / (1+d_i)*(1+beta_2));

						//particles move applying the rules for new velocity and new position
						swarm.at(i)->move(config, minBound,maxBound,iteration);

						if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
							updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
							best_particle = swarm[i];
						}
						break;
					}
				}
			}
		}
		else {
			// the inertia weight of the first iteration is computed in Swarm::computeInertia
			for (unsigned int i=0;i<swarm.size();i++){
				//particles move applying the rules for new velocity and new position
				swarm.at(i)->move(config, minBound,maxBound,iteration);

				if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
					updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
					best_particle = swarm[i];
				}
			}
		}

	}
	//Any other strategy
	else {
		for (unsigned int i=0;i<swarm.size();i++){

			//particles move applying the rules for new velocity and new position
			swarm.at(i)->move(config, minBound,maxBound,iteration);

			if (swarm.at(i)->getPbestEvaluation() < global_best.eval){
				updateGlobalBest(swarm.at(i)->getPbestPosition(), swarm.at(i)->getPbestEvaluation());
				best_particle = swarm[i];
			}
		}

	}

	if(prev_eval > global_best.eval) {
		//cout << "Iteration "<< iteration << " [ New global best: " << global_best.eval << " ]" << endl;
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

void Swarm::createRingTopology(){       // Every particle is neighbor of the adjacent particles
	int a,b;
	for(unsigned int i=0;i<swarm.size();i++){
		a=i-1;
		b=i+1;
		if(i==0)
			a=swarm.size()-1;
		if(i==(swarm.size()-1))
			b=0;

		swarm.at(i)->addNeighbour(swarm.at(a));
		swarm.at(i)->addNeighbour(swarm.at(b));
	}
}

void Swarm::createStarTopology(){        //A.k.a. Global best, all particles are neighbors of one central particle
	for(unsigned int i=1;i<swarm.size();i++){
		swarm.at(i)->addNeighbour(swarm.at(0));
		swarm.at(0)->addNeighbour(swarm.at(i));
	}
}

void Swarm::createRandomTopology(){
	int randomEdge;
	for(unsigned int i=0;i<swarm.size();i++){
		randomEdge = (int)floor(RNG::randVal(0.0,(double)swarm.size()));
		swarm.at(i)->addNeighbour(swarm.at(randomEdge));
	}
}

void Swarm::createVonNeumannTopology(){
	int a,b,c;
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

		swarm.at(i)->addNeighbour(swarm.at(a));
		swarm.at(i)->addNeighbour(swarm.at(b));
		swarm.at(i)->addNeighbour(swarm.at(c));
	}
}

void Swarm::updateTimeVaryingTopology(Configuration* config, long int iterations){
	//Topology update
	if( (iterations > 0) && (config->getEsteps() < swarm.size()-3) && (iterations%config->getTopologyUpdatePeriod() == 0)){
		unsigned int removals = 0;
		//cout << " -- esteps " << config->getEsteps() << endl;
		//cout << " -- Update topology at iteration: " << iterations << " Target: " << swarm.size()-(2+config->getEsteps()) << endl;
		RNG::shufflePermutation();
		while(removals < swarm.size()-(2+config->getEsteps())){
			for(unsigned int i=0;i<swarm.size();i++){
				int particleIndex = RNG::getPermutationElement(i);
				if( swarm.at(particleIndex)->getNeighborhoodSize() > 3 ){ //3 because a particle is a neighbor to itself
					int neighborID = swarm.at(particleIndex)->getRandomNonAdjacentNeighborID(config);

					//cout << " -- Erasing edge " << particleIndex << " <---> " << neighborID << endl;

					swarm.at(particleIndex)->eraseNeighborbyID(neighborID);
					swarm.at(neighborID)->eraseNeighborbyID(particleIndex);

					removals++;
				}
				if(removals == swarm.size()-(2+config->getEsteps()))
					break;
			}
		}
		config->setEsteps(config->getEsteps()+1);
		//cout << "Removals " << removals << endl;
		//for(unsigned int i=0;i<swarm.size();i++)
		//	cout << "particleIndex: " << i << " -- Neighbors: " << swarm.at(i)->getNeighborhoodSize() << endl;
	}
}

Solution Swarm::getGlobalBest(){
	return global_best;
}

double Swarm::computeAvgVelocity(Configuration* config){
	double sumVelocity = 0.0;
	for (unsigned int i=0;i<swarm.size();i++){
		for (unsigned int j=0;j<config->getProblemDimension();j++){
			sumVelocity = sumVelocity + abs(swarm.at(i)->getCurrentVelocity()[j]);
		}
	}
	return (sumVelocity/(swarm.size()*config->getProblemDimension()));
}

void Swarm::rankParticles(SimplifySwarm* simpSwarm){
	//The ranking are obtained via the mergesort algorithm implemented in Utils.cpp

	//Copy particle's id and evaluation in simpSwarm
	for (unsigned int i=0;i<swarm.size();i++){
		simpSwarm->id[i] = swarm.at(i)->getID();
		simpSwarm->eval[i] = swarm.at(i)->getCurrentEvaluation();
	}
	//Sort simpSwarm by the value of particles' evaluation
	Utils::mergeSort(simpSwarm, 0, swarm.size()-1); //mergeSort(array, left(LOWER) index, right (UPPER) index);

	//Set particles' rank
	for (unsigned int i=0;i<swarm.size();i++){
		swarm.at(simpSwarm->id[i])->setRanking(i+1);
	}

}

// This function computes the inertia weight according to the selected strategy
void Swarm::computeInertia(Configuration* config, long int iteration){
	/* Non-adaptive strategies */
	//IW_L_INC - 1 - Linear increasing
	if (config->getinertiaCS() == IW_L_INC) {
		//from Frankenstein's PSO
		if(iteration <= config->getIWSchedule()){
			config->setInertia(
					((double)(config->getIWSchedule() - iteration)/config->getIWSchedule())*
					(config->getFinalIW() - config->getInitialIW()) + config->getInitialIW()
			);
			//cout << config->getInertia() << endl;
		}
		else
			config->setInertia(config->getFinalIW());
	}
	//IW_L_DEC - 2 - Linear decreasing
	else if (config->getinertiaCS() == IW_L_DEC) {
		//from Frankenstein's PSO
		if(iteration <= config->getIWSchedule()){
			config->setInertia(
					((double)(config->getIWSchedule() - iteration)/config->getIWSchedule())*
					(config->getInitialIW() - config->getFinalIW()) + config->getFinalIW()
			);
			//cout << config->getInertia() << endl;
		}
		else
			config->setInertia(config->getFinalIW());
	}
	//IW_RANDOM - 3 - Random
	else if (config->getinertiaCS() == IW_RANDOM) {
		config->setInertia( 0.5 * (problem->getRandom01()/2.0));
		//cout << config->getInertia() << endl;
	}
	//IW_NONL_DEC - 4 - Nonlinear decreasing
	else if (config->getinertiaCS() == IW_NONL_DEC) {
		config->setInertia(
				config->getFinalIW() - (config->getFinalIW()-config->getInitialIW())*
				pow((double)(iteration)/config->getMaxIterations(),alpha)
		);
		//cout << iteration << endl;
		//cout << config->getInertia() << endl;
	}
	//IW_NONL_DEC_IMP - 5 - Nonlinear decreasing improved
	else if (config->getinertiaCS() == IW_NONL_DEC_IMP) {
		config->setInertia(
				pow(omegaXu, iteration)
		);
		//cout << iteration << endl;
		//cout << config->getInertia() << endl;
	}
	//IW_NONL_DEC_TIME - 6 - Nonlinear decreasing time-dependent
	else if (config->getinertiaCS() == IW_NONL_DEC_TIME) {
		config->setInertia(
				pow((2.0/iteration), 0.3)
		);
		//cout << iteration << endl;
		//cout << config->getInertia() << endl;
	}
	//IW_CHAOTIC_DEC - 7 Chaotic decreasing
	else if (config->getinertiaCS() == IW_CHAOTIC_DEC) {
		iteration == 1 ? zFunction = problem->getRandom01() : zFunction = 4*zFunction*(1-zFunction);
		config->setInertia(
				(zFunction*config->getInitialIW()) + (config->getFinalIW()-config->getInitialIW()) *
				((double)(config->getMaxIterations())-iteration)/config->getMaxIterations()
		);
		//cout << iteration << endl;
		//cout << config->getInertia() << endl;
	}
	//IW_EXP_DEC - 8 - Natural exponential decreasing
	else if (config->getinertiaCS() == IW_EXP_DEC) {
		config->setInertia(
				config->getInitialIW() + (config->getFinalIW()-config->getInitialIW())*
				exp((-10.0*iteration)/config->getMaxIterations())
		);
		//cout << iteration << endl;
		//cout << config->getInertia() << endl;
	}
	//IW_OSCILLATING - 9 - Oscillating
	else if (config->getinertiaCS() == IW_OSCILLATING) {
		if (iteration < (3*config->getMaxIterations())/4)
			config->setInertia(
					((config->getInitialIW() + config->getFinalIW()) /2.0) +
					((config->getFinalIW() + config->getInitialIW()) /2.0) *
					cos((simNumOfCos*iteration)/(config->getMaxIterations()*3.0))
			);
		else
			config->setInertia(config->getInitialIW());
		//cout << iteration << endl;
		//cout << config->getInertia() << endl;
	}
	//IW_LOG_DEC - 10 - Logarithm decreasing
	else if (config->getinertiaCS() == IW_LOG_DEC) {
		config->setInertia(
				config->getFinalIW() + (config->getInitialIW()-config->getFinalIW())*
				log10(((10.0*iteration)/config->getMaxIterations())+ a )
		);
		//cout << iteration << endl;
		//cout << config->getInertia() << endl;
	}


	/* Adaptive strategies */
	//IW_SELF_REGULATING - 11 - Self-regulating
	else if (config->getinertiaCS() == IW_SELF_REGULATING) {
		static double deltaOmega = (((double)config->getFinalIW()) - config->getInitialIW())/config->getMaxIterations();
		iteration == 1 ? omega_2 = config->getFinalIW() : omega_2 = omega_2-deltaOmega;
		config->setInertia(omega_2);
		//cout << iteration << endl;
		//cout << config->getInertia() << endl;
	}
	//IW_VELOCITY_BASED - 12 - Based on velocity information
	else if (config->getinertiaCS() == IW_VELOCITY_BASED) {

		static double T_0_95 = 95*config->getMaxIterations()/100; //iteration at which 95% of search process is completed

		if (iteration == 1)
			config->setInertia(config->getFinalIW());

		idealVelocity = swarm.at(0)->getMaxVelLimit() * ((1.0 + cos(M_PI*(iteration/T_0_95)))/2);
		avVel=computeAvgVelocity(config);	//average absolute velocity of the swarm

		if (avVel >= idealVelocity){
			(config->getInertia()-deltaOmega) >= config->getInitialIW() ?
					config->setInertia(config->getInertia()-deltaOmega) : config->setInertia(config->getInitialIW());
		}
		else{
			(config->getInertia()+deltaOmega) >= config->getFinalIW() ?
					config->setInertia(config->getFinalIW()) : config->setInertia(config->getInertia()+deltaOmega);
		}
		//cout << iteration << endl;
		//cout << avVel << " -- " << idealVelocity << endl;
		//cout << config->getInertia() << endl;
	}
	//IW_DOUBLE_EXP - 13 - Double exponential self-adaptive
	else if (config->getinertiaCS() == IW_DOUBLE_EXP) {
		//This strategy is implemented directly in Particle::move because the inertia weight is computed independently for each particle
	}
	//IW_RANKS_BASED - 14 - Rank-based
	else if (config->getinertiaCS() == IW_RANKS_BASED) {
		if (iteration == 1) {
			//Memory allocation to rank particles
			simpSwarm.eval = new long double [config->getSwarmSize()];
			simpSwarm.id = new int [config->getSwarmSize()];
		}
		//Rank particles
		rankParticles(&simpSwarm);
		//The computation of the inertia weight is implemented in Particle::move
	}
	//IW_SUCCESS_BASED - 15 Success-based
	else if (config->getinertiaCS() == IW_SUCCESS_BASED) {
		if (iteration == 1) {
			//simpSwarm contains a simplified copy of the swarm at t-1
			simpSwarm.eval = new long double [config->getSwarmSize()];
			simpSwarm.id = new int [config->getSwarmSize()];
			//Copy particle's id and evaluation in simpSwarm
			for (unsigned int i=0;i<swarm.size();i++){
				simpSwarm.id[i] = swarm.at(i)->getID();
				simpSwarm.eval[i] = swarm.at(i)->getCurrentEvaluation();
			}
			config->setInertia(config->getFinalIW()); //set inertia weight to its maximum value for the first iteration
		}
		else{
			int S_i = 0; //Number of solutions that improved after the last iteration
			for (unsigned int i=0; i<swarm.size(); i++){
				for (unsigned int j=0; j<sizeof(simpSwarm.id); j++){
					if (swarm.at(i)->getID() == simpSwarm.id[j]){
						//evaluate if the solution improved
						if (swarm.at(i)->getCurrentEvaluation() < simpSwarm.eval[j]){
							S_i++;
							break;
						}
						else
							break;
					}
				}
			}
			//set the inertia weight
			config->setInertia( config->getInitialIW() + ((config->getFinalIW()-config->getInitialIW()) *
					((double)S_i/swarm.size()))
			);
			//Copy particle's id and evaluation in simpSwarm for the next iteration
			for (unsigned int i=0;i<swarm.size();i++){
				simpSwarm.id[i] = swarm.at(i)->getID();
				simpSwarm.eval[i] = swarm.at(i)->getCurrentEvaluation();
			}
			//cout << iteration << endl;
			//cout << "P[" << id << "]" << "rank: "<< ranking << " eval:" << current.eval << endl;
			//cout << config->getInertia() << endl;
		}

	}
	//IW_CONVERGE_BASED - 16 Convergence-based
	else if (config->getinertiaCS() == IW_CONVERGE_BASED) {
		if (iteration == 1) {
			//simpSwarm contains a simplified copy of the swarm at t-1
			simpSwarm.eval = new long double [config->getSwarmSize()];
			simpSwarm.id = new int [config->getSwarmSize()];
			//Copy particle's id and evaluation in simpSwarm
			for (unsigned int i=0;i<swarm.size();i++){
				simpSwarm.id[i] = swarm.at(i)->getID();
				simpSwarm.eval[i] = swarm.at(i)->getCurrentEvaluation();
			}
			config->setInertia(1 - abs(alpha_2/(1+beta_2)));
		}
		else{
			//Copy particle's id and evaluation in simpSwarm for the next iteration
			for (unsigned int i=0;i<swarm.size();i++){
				simpSwarm.id[i] = swarm.at(i)->getID();
				simpSwarm.eval[i] = swarm.at(i)->getCurrentEvaluation();
			}
		}
	}
	//No strategy, inertia is constant during the execution
	else {
		config->setInertia(config->getInertia()); //kind of unnecessary, but ensures integrity
		//cout << inertia << endl;
		//return inertia;
	}
}
