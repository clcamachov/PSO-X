/*
 * config.h
 *
 *  Created on: May 31, 2018
 *      Author: leonardo
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <sstream> // for std::stringstream

//CEC2005 problems
#define	SHIFTED_SPHERE_CEC05									0
#define	SHIFTED_SCHWEFEL12_CEC05								1
#define	SHIFTED_ROTATED_HIGH_CONDITIONED_ELLIPTIC_CEC05			2
#define	NOISE_SHIFTED_SCHWEFEL12_CEC05							3
#define	SCHWEFEL26_GLOBAL_OPTIMUM_ON_BOUNDS_CEC05				4
#define	SHIFTED_ROSENBROCK_CEC05								5
#define	SHIFTED_ROTATED_GRIEWANK_CEC05							6
#define	SHIFTED_ROTATED_ACKLEY_GOOB_CEC05						7
#define	SHIFTED_RASTRIGIN_CEC05									8
#define	SHIFTED_ROTATED_RASTRIGIN_CEC05							9
#define	SHIFTED_ROTATED_WEIERSTRASS_CEC05						10
#define	BASIC_SCHWEFEL213_CEC05									11
#define	SHIFTED_EXPANDED_GRIEWANKROSENBROCK_CEC05				12
#define	SHIFTED_ROTATED_EXPANDED_SCHAFFERF6_CEC05				13
#define	BASIC_HYBRIDCOMPOSITION1								14
#define	ROTATED_HYBRIDCOMPOSITION1								15
#define	NOISE_ROTATED_HYBRIDCOMPOSITION1						16
#define	ROTATED_HYBRIDCOMPOSITION2								17
#define	ROTATED_HYBRIDCOMPOSITION2_NBGO							18
#define	ROTATED_HYBRIDCOMPOSITION2_GOOB							19
#define	ROTATED_HYBRIDCOMPOSITION3								20
#define	ROTATED_HYBRIDCOMPOSITION3_HCNM							21
#define	NONCONTINUOUS_ROTATED_HYBRIDCOMPOSITION3				22
#define	ROTATED_HYBRIDCOMPOSITION4								23
#define	ROTATED_HYBRIDCOMPOSITION4_NO_BOUNDS					24

//CEC2014 problems
#define SHIFTED_ROTATED_HIGH_CONDITIONED_ELLIPTIC_CEC14			0
#define ROTATED_BENT_CIGER_CEC14								1
#define ROTATED_DISCUS_CEC14									2
#define SHIFTED_ROTATED_ACKLEY_CEC14							3
#define SHIFTED_ROTATED_ROSENBROCK_CEC14						4
#define SHIFTED_ROTATED_GRIEWANK_CEC14							5
#define SHIFTED_ROTATED_RASTRIGIN_CEC14							6
#define SHIFTED_SCHWEFEL_CEC14									7
#define SHIFTED_ROTATED_SCHWEFEL_CEC14							8
#define SHIFTED_ROTATED_WEIERSTRASS_CEC14						9
#define SHIFTED_ROTATED_KATSUURA_CEC14							10
#define SHIFTED_ROTATED_HAPPYCAT_CEC14							11
#define SHIFTED_ROTATED_HGBAT_CEC14								12
#define H1_CEC14_CEC14											13
#define H2_CEC14_CEC14											14
#define H3_CEC14_CEC14											15
#define H4_CEC14_CEC14											16
#define H5_CEC14_CEC14											17
#define H6_CEC14_CEC14											18

//SOCO problems
#define	SHIFTED_SPHERE_SOCO 					0
#define SHIFTED_SCHWEFEL221_SOCO				1
#define SHIFTED_ROSENBROCK_SOCO 				2
#define SHIFTED_RASTRIGIN_SOCO					3
#define SHIFTED_GRIEWANK_SOCO 					4
#define SHIFTED_ACKLEY_SOCO 					5
#define SHIFTED_SCHWEFEL222_SOCO 				6
#define SHIFTED_SCHWEFEL12_SOCO					7
#define SHIFTED_EXTENDED_F10_SOCO 				8
#define SHIFTED_BOHACHEVSKY_SOCO 				9
#define SHIFTED_SCHAFFER_SOCO 					10
#define H1 										11
#define H2										12
#define H3										13
#define H4										14
#define H7										15
#define H8										16
#define H9										17
#define H10										18

//ABC-X PROBLEMS (MIXTURE)
//UNIMODAL FUNCTIONS
#define SHIFTED_SPHERE									0 //SOCO
#define SHIFTED_ROTATED_HIGH_CONDITIONED_ELLIPTIC		1 // CEC14
#define ROTATED_BENT_CIGER								2 // CEC14
#define ROTATED_DISCUS									3 // CEC14
#define SHIFTED_SCHWEFEL221								4 // SOCO
#define SHIFTED_SCHWEFEL12								5 // SOCO
#define SHIFTED_SCHWEFELS12_NOISE_IN_FITNESS			6 // CEC05
#define SHIFTED_SCHWEFEL222								7 // SOCO
#define SHIFTED_EXTENDED_F10 							8 // SOCO
#define SHIFTED_BOHACHEVSKY 							9 // SOCO
#define SHIFTED_SCHAFFER 								10 // SOCO
#define SCHWEFEL26_GLOBAL_OPTIMUM_ON_BOUNDS				11 // CEC05
//MULTIMODAL FUNCTIONS
#define SHIFTED_ACKLEY									12 // SOCO
#define SHIFTED_ROTATED_ACKLEY							13 // CEC14
#define	SHIFTED_ROSENBROCK								14 // SOCO
#define	SHIFTED_ROTATED_ROSENBROCK						15 // CEC14
#define SHIFTED_GRIEWANK								16 // SOCO
#define SHIFTED_ROTATED_GRIEWANK						17 // CEC14
#define SHIFTED_RASTRIGIN								18 // SOCO
#define	SHIFTED_ROTATED_RASTRIGIN						19 // CEC14
#define	SHIFTED_SCHWEFEL								20 // CEC14
#define SHIFTED_ROTATED_SCHWEFEL						21 // CEC14
#define	SHIFTED_ROTATED_WEIERSTRASS						22 // CEC14
#define SHIFTED_ROTATED_KATSUURA						23 // CEC14
#define	SHIFTED_ROTATED_HAPPYCAT						24 // CEC14
#define	SHIFTED_ROTATED_HGBAT							25 // CEC14
//HYBRID FUNCTIONS
#define	H1_SOCO								26 // SOCO
#define H2_SOCO								27 // SOCO
#define H3_SOCO								28 // SOCO
#define H4_SOCO								29 // SOCO
#define H7_SOCO								30 // SOCO
#define H8_SOCO								31 // SOCO
#define H9_SOCO								32 // SOCO
#define H10_SOCO							33 // SOCO
#define	H1_CEC14							34 //CEC14
#define H2_CEC14							35 //CEC14
#define H3_CEC14							36 //CEC14
#define H4_CEC14							37 //CEC14
#define H5_CEC14							38 //CEC14
#define H6_CEC14							39 //CEC14
//COMPOSITION FUNCTIONS
#define COMPOSITION_F1						40 //CEC05
#define COMPOSITION_F2						41 //CEC05
#define COMPOSITION_F3						42 //CEC05
#define COMPOSITION_F4						43 //CEC05
#define COMPOSITION_F5						44 //CEC05
#define COMPOSITION_F6						45 //CEC05
#define COMPOSITION_F7						46 //CEC05
#define COMPOSITION_F8						47 //CEC05
#define COMPOSITION_F9						48 //CEC05
#define COMPOSITION_F10						49 //CEC05

//BBOB
//#define SPHERE_FUNCTION_BBOB											0 //BBOB
//#define ELLIPSOIDAL_FUNCTION_BBOB										1 //BBOB
//#define RASTRIGIN_FUNCTION_BBOB											2 //BBOB
//#define BUCHE_RASTRIGIN_FUNCTION_BBOB									3 //BBOB
//#define LINEAR_SLOPE_BBOB												4 //BBOB
//#define ATTRACTIVE_SECTOR_FUNCTION_BBOB									5 //BBOB
//#define STEP_ELLIPSOIDAL_FUNCTION_BBOB									6 //BBOB
//#define ROSENBROCK_FUNCTION_ORIGINAL_BBOB								7 //BBOB
//#define ROSENBROCK_FUNCTION_ROTATED_BBOB								8 //BBOB
//#define ELLIPSOIDAL_FUNCTION_BBOB										9 //BBOB
//#define DISCUS_FUNCTION_BBOB											10 //BBOB
//#define BENT_CIGAR_FUNCTION_BBOB										11 //BBOB
//#define SHARP_RIDGE_FUNCTION_BBOB										12 //BBOB
//#define DIFFERENT_POWERS_FUNCTION_BBOB									13 //BBOB
//#define RASTRIGIN_FUNCTION_BBOB											14 //BBOB
//#define WEIERSTRASS_FUNCTION_BBOB										15 //BBOB
//#define SCHAFFERS_F7_FUNCTION_BBOB										16 //BBOB
//#define SCHAFFERS_F7_FUNCTION_MODERATELY_ILL_CONDITIONED_BBOB			17 //BBOB
//#define COMPOSITE_GRIEWANK_ROSENBROCK_FUNCTION_F8F2_BBOB				18 //BBOB
//#define SCHWEFEL_FUNCTION_BBOB											19 //BBOB
//#define GALLAGHER_GAUSSIAN_101_ME_PEAKS_FUNCTION_BBOB					20 //BBOB
//#define GALLAGHER_GAUSSIAN_21_HI_PEAKS_FUNCTION_BBOB					21 //BBOB
//#define KATSUURA_FUNCTION_BBOB											22 //BBOB
//#define LUNACEK_BI_RASTRIGIN_FUNCTION_BBOB								23 //BBOB

//competitions
#define CEC05	 				0
#define CEC14					1
#define SOFT_COMPUTING 			2
#define MIXTURE					3
#define BBOB					4

//Useful constants
#define PI acos(-1)
#define E  exp(1)
#define EPSILON 1E-15
#define DELTA 1E-30
#define CONSTRICTION_COEFFICIENT 0.7298
#define REINIT_PRECISION 0.00001
#define STAGNATION_PRECISION 0.0000001
#define MAX_DIMENSION 	1000
#define LINE_BUF_LEN    100
#define TRACE( x )
#define PRINTER(name) printer(#name, (name))

/*PSO definitions*/
// available topologies
#define TOP_FULLYCONNECTED 		0	//or STAR
#define TOP_RING 				1
#define TOP_WHEEL 				2
#define TOP_RANDOM 				3
#define TOP_VONNEUMANN			4
#define TOP_TIMEVARYING			5
#define TOP_HIERARCHICAL		6

// available models of influence
#define MOI_BEST_OF_N			0
#define MOI_FI					1
#define MOI_RANKED_FI			2

// available inertia control strategies (omega1)
// non adaptive
#define IW_CONSTANT						 0
#define IW_L_INC						 1
#define IW_L_DEC						 2
#define IW_RANDOM						 3
#define IW_NONL_DEC		 				 4
#define IW_NONL_DEC_IMP					 5
#define IW_NONL_DEC_TIME				 6
#define IW_CHAOTIC_DEC		 			 7
#define IW_EXP_DEC			 			 8
#define IW_OSCILLATING					 9
#define IW_LOG_DEC		 	 			10
// adaptive
#define IW_SELF_REGULATING 	 			11
#define IW_VELOCITY_BASED 				12
#define IW_DOUBLE_EXP 	 	 			13
#define IW_RANKS_BASED 		 			14
#define IW_SUCCESS_BASED 	 			15
#define IW_CONVERGE_BASED 	 			16

// available omega2 strategies
#define O2_EQUAL_TO_O1					0
#define O2_ZERO							1  //this makes the component Distribution = 0
#define O2_ONE							2
#define O2_RANDOM						3
#define O2_CONSTANT						4

// available omega3 strategies
#define O3_EQUAL_TO_O1					0
#define O3_ZERO							1  //this makes the component Perturbation = 0
#define O3_ONE							2
#define O3_RANDOM						3
#define O3_CONSTANT						4

// available acceleration coefficients strategies
#define AC_CONSTANT						0
#define AC_TIME_VARYING					1
#define AC_EXTRAPOLATED					2
#define AC_RANDOM						3

// available velocity update rules
#define VEL_BASIC						0
#define VEL_STANDARD					1
#define VEL_LINEAR						2
#define VEL_CONSTRICTED					3
#define VEL_GUARAN_CONVERG				4
#define VEL_FULLY_INFORMED				5
#define VEL_LOC_CON_TRANS_INV			6
#define VEL_STANDARD2011				7	//rotation invariant
#define VEL_ROTATION_INV				8	//random matrix

#define DIST_RECTANGULAR				0
#define DIST_SPHERICAL					1
#define DIST_ADD_STOCH 					2

#define PERT1_NONE						0
#define PERT1_NORMAL_DISTANCE			1
#define PERT1_NORMAL_SUCCESS			2
#define PERT1_CAUCHY_DISTANCE			3
#define PERT1_CAUCHY_SUCCESS			4

#define PERT1_DIST_L_INDEPENDENT		0
#define PERT1_DIST_L_USER_SUPPLIED		1

#define PERT2_NONE						0
#define PERT2_ADD_RECT					1
#define PERT2_ADD_NOISY					2

#define Q_STANDARD						0
#define Q_GAUSSIAN						1
#define Q_DISCRETE_2					2
#define Q_CAUCHY_NORMAL					3

#define MATRIX_NONE						0
#define MATRIX_DIAGONAL					1
#define MATRIX_LINEAR					2
#define MATRIX_RRM_EXP_MAP				3
#define MATRIX_RRM_EUCLIDEAN_ONE		4
#define MATRIX_RRM_EUCLIDEAN_ALL		5

#define ANGLE_CONSTANT					0
#define ANGLE_NORMAL					1
#define ANGLE_ADAPTIVE					2

#define POP_CONSTANT					0
#define POP_LADDERED					1
#define POP_INCREMENTAL					2

#define PARTICLE_INIT_RANDOM 			0
#define PARTICLE_INIT_MODEL 			1

class Configuration {

private:
	//general parameters
	unsigned long rngSeed;
	unsigned int maxFES;
	long int max_iterations;

	//Problem parameters
	unsigned int competitionID;
	unsigned int problemID;
	unsigned int problemDimension;
	double minInitRange;
	double maxInitRange;

	//Other parameters
	double startTime; // to calculate time elapsed
	bool useLogs;
	bool verbose;
	std::string outputPath;
	bool indStrategies;

	//Population
	long int particles;
	int populationCS;
	long int initialPopSize;
	long int finalPopSize;
	int particlesToAdd;
	int p_intitType;		//new particle initialization type

	//Acceleration coefficients
	double accelCoeffCS;
	double phi_1;
	double phi_2;
	double initialPhi1;
	double initialPhi2;
	double finalPhi1;
	double finalPhi2;

	//Topology parameters
	short topology;             // topology indicator
	unsigned int tSchedule;		//n, 2n, 3n, 4n, etc. (the lower the value the faster)
	unsigned int esteps;
	int topologyUpdatePeriod;
	int branching;

	//Model of influence
	short modelOfInfluence;

	//Inertia control parameters (omega1 in the GVU)
	short omega1CS;				// omega1 control strategy indicator (this is the inertiaCS)
	double inertia;				// actual variable of inertia. If no inertiaCS is given, this value is fixed during the
	double initialIW;
	double finalIW;
	unsigned int iwSchedule;	//n^2 , 2n^2 , 3n^2 , 4n^2, etc. (the lower the value the faster)
	short omega2CS;
	double omega2;
	short omega3CS;
	double omega3;
	//Additional values for inertia control parameter omega1CS (only adaptive)
	double iw_par_eta;				//from 0.1 to 1 in IW_SELF_REGULATING - 11
	double iw_par_lambda;		//from 0.1 to 1 small positive constant in IW_VELOCITY_BASED - 12
	double iw_par_alpha_2;			//from 0 to  1 in IW_CONVERGE_BASED - 16
	double iw_par_beta_2;			//from 0 to  1 in IW_CONVERGE_BASED - 16

	//Velocity
	bool useVelClamping;

	//Position
	bool reinitializePosition;

	//Perturbation
	short perturbation1; //distribution-based
	short perturbation2; //additive
	short pert1_par_l_CS;
	double pert1_par_l; //scaling factor for PERT1_NORMAL_DISTANCE
	int pert1_2_par_success; 		//success threshold of the additive rectangular perturbation
	int pert1_2_par_failure; 		//failure threshold of the additive rectangular perturbation

	//Random matrix
	short randomMatrix;
	short angleCS;
	double angleSD;
	double angle_par_alpha;
	double angle_par_beta;
	double rotation_angle; //self-explanatory

	//NPPDistribution
	short distributionNPP;
	short operator_q;
	bool randNeighbor;

	//Velocity rules
	int vRule;

public:

	Configuration();
	~Configuration();

	//Program
	bool getConfig(int argc, char *argv[]);
	void setStartTime(double stime);
	double getStartTime();

	void printParameters();
	void printUsage();
	void print();
	bool logOutput();
	bool verboseMode();
	std::string getOutputPath();

	//Problem
	unsigned long getRNGSeed();
	unsigned int getMaxFES();
	unsigned int getMaxIterations();

	void setDefaultParameters();
	unsigned int getCompetitionID();
	unsigned int getProblemID();
	unsigned int getProblemDimension();
	double getMinInitRange();
	double getMaxInitRange();

	void setMinInitRange(double lowerlimit);
	void setMaxInitRange(double upperlimit);

	double getMinInitBound();
	double getMaxInitBound();

	//PSO
	long int getSwarmSize();
	void setSwarmSize(long int new_size);
	int getPopulationCS();
	long int getInitialPopSize();
	long int getFinalPopSize();
	int getParticlesToAdd();
	void setParticlesToAdd(int new_pool_size);
	int getParticleInitType();
	bool useIndStrategies();

	//Velocity
	void setVelocityRule(int rule);
	int getVelocityRule();
	bool useVelocityClamping();

	//Inertia control (Omega1)
	short getOmega1CS();
	double getOmega1();
	void setOmega1(double new_inertia);
	double getInitialIW();
	double getFinalIW();
	unsigned int getIWSchedule();

	//Velocity
	bool isVelocityClamped();
	void setVelocityClamped(bool clamping);

	//Position
	bool useReinitialization();
	void setReinitialization(bool reinitilized);

	double get_iw_par_eta();
	double get_iw_par_deltaOmega();
	double get_iw_par_alpha_2();
	double get_iw_par_beta_2();

	//Omega2 and Omega3
	short getOmega2CS();
	short getOmega3CS();
	double getOmega2();
	void setOmega2(double new_omega2);
	double getOmega3();
	void setOmega3(double new_omega3);

	//Acceleration coefficients
	short getAccelCoeffCS();
	double getPhi1();
	double getPhi2();
	double getInitialPhi1();
	double getInitialPhi2();
	double getFinalPhi1();
	double getFinalPhi2();

	//Model of influence
	short getModelOfInfluence();

	//Perturbation
	short getPerturbation1Type();
	short getPerturbation2Type();
	short getPert1_par_l_CS();
	double getPert1_par_l();
	void setPert1_par_l(double l_val);

	short getRandomMatrix();
	short getAngleCS();
	double getRotationAgle();
	void setRotationAgle(double angle);
	double getAngleSD();
	void setAngleSD(double angle_sd);
	double get_angle_par_alpha();
	double get_angle_par_beta();

	int get_pert1_2_par_success();
	int get_pert1_2_par_failure();

	//Distribution
	short getDistributionNPP();
	short getOperator_q();
	bool getRandNeighbor();

	//Topology
	short getTopology();
	unsigned int getTopologySchedule();
	int getBranchingDegree();

	void setEsteps(unsigned int num_esteps);
	unsigned int getEsteps();

	void setTopologyUpdatePeriod(int period);
	int getTopologyUpdatePeriod();

};

#endif /* CONFIG_H_ */
