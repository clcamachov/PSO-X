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

using namespace std;

/* Solution structure */
struct SimplifySwarm {
	long double* eval;  /* value of the solution */
	int* id;  			/* particle id */
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
};



#endif /* UTILS_H_ */
