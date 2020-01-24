/*
 * utils.h
 *
 *  Created on: May 31, 2018
 *      Author: leonardo
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <string>
#include <vector>

using namespace std;

///* Simplified swarm structure */
//struct SimplifySwarm {
//	long double* eval;  /* value of the solution */
//	int* id;  			/* particle id */
//};

/* Simplified swarm structure */
struct SimplifySwarm {
	long double eval;  /* value of the solution */
	int id;  		   /* particle id */
};

///* Simplified swarm structure */
//struct SimplifySwarmArray {
//	long double* eval;  /* value of the solution */
//	int* id;  		   /* particle id */
//};

/* Data structure to organize particle hierarchically */
struct Node { // Represents a node of an n-ary tree
	int id;
	long double function_eval;
	int numSubNodos;
	struct Node *parent;
	struct Node **child;
};

class Utils{

public:
	static void loadRowVectorFromFile(string filename, int columns, double * row);
	static void loadRowVector(ifstream & fs, int columns, double * row);
	static void loadColumnVectorFromFile(string filename, int columns, double * column);
	static void loadColumnVector(ifstream & fs, int columns, double * column);
	static void loadMatrixFromFile(string filename, int N, int rows, int columns, double *** matrix3D);
	static void loadMatrixFromFile(string filename, int rows, int columns, double** matrix2D);
	static void loadMatrix(ifstream & fs, int rows, int columns, double ** matrix);
	static void mergeSort(SimplifySwarm* arr, int l, int r);
	static void merge(SimplifySwarm* arr, int l, int m, int r);
	static Node* newNode(Node* parent, int numSubNodos, int id, long double valor);
};



#endif /* UTILS_H_ */
