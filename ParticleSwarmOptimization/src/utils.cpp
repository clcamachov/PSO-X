/*
 * utils.cpp
 *
 *  Created on: Jun 4, 2018
 *      Author: leonardo
 */
#include "utils.h"

#include <string>
#include <sstream>
#include <fstream>

using namespace std;

void Utils::loadRowVector(ifstream & fs, int columns, double * row){

	//string to storage on line read from the file
	string aline, arow;

	/* The function std::getline takes a stream that it will extract a line
	 * from as its first argument. The second argument is a std::string
	 * that it will extract the line into. Without a third argument, the ending
	 *  of a "line" is considered to be where we see a newline character.
	 */
	getline(fs, aline);

	/* The std::stringstream type wraps a string with a stream interface.
	 * It means you can treat it like cout and cin, using << and >> to insert
	 *  and extract things from the string.*/
	istringstream lineStream(aline); //

	//cout << "  [lineStream]  ";
	for (int i = 0 ; i < columns ;i++) {
		lineStream >> row[i];
		//cout << row[i] << " ";
	}
	//cout << "  [lineStream]" << endl;

}

void Utils::loadColumnVector(ifstream & fs, int rows, double * column){

	string aline;
	getline(fs, aline);
	istringstream lineStream(aline);

	//cout << "  [lineStream]  ";
	for (int i = 0 ; i < rows ;i++) {
		lineStream >> column[i];
		//cout << column[i] << " ";
	}
	//cout << "  [lineStream]" << endl;
}

void Utils::loadRowVectorFromFile(string filename, int columns, double * row){
	ifstream file(filename.c_str());
	if (file){
		loadRowVector(file, columns, row);
		file.close();
	}
	else
		cout << "  " << filename << " can't be opened" << endl;

}


void Utils::loadColumnVectorFromFile(string filename, int rows, double * column){
	ifstream file(filename.c_str());
	if (file){
		loadColumnVector(file, rows, column);
		file.close();
	}
	else
		cout << "  " << filename << " can't be opened" << endl;
}

void Utils::loadMatrixFromFile(string filename, int N, int rows, int columns, double *** matrix3D){
	ifstream file(filename.c_str());

	if (file){
		for (int i = 0 ; i < N ; i++) {
			loadMatrix(file, rows, columns, matrix3D[i]);
		}
		file.close();
	}
	else
		cout << "  " << filename << " can't be opened" << endl;
}


void Utils::loadMatrixFromFile(string filename, int rows, int columns, double** matrix2D){
	ifstream file(filename.c_str());
	if (file){
		loadMatrix(file, rows, columns, matrix2D);
		file.close();
	}
	else
		cout << "  " << filename << " can't be opened" << endl;
}

void Utils::loadMatrix(ifstream & fs, int rows, int columns, double ** matrix){
	for (int i = 0 ; i < rows ; i ++) {
		loadRowVector(fs, columns, matrix[i]);
	}
}

//mergeSort(arr, 0, arr_size-1);
void Utils::mergeSort(SimplifySwarm* arr, int l, int r){
	if (l < r) {
		// Middle point of the array
		int m = l+(r-l)/2;

		// Sort first and second halves
		mergeSort(arr, l, m);
		mergeSort(arr, m+1, r);
		merge(arr, l, m, r);
	}
}

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void Utils::merge(SimplifySwarm* arr, int l, int m, int r){
	int i, j, k;
	int n1 = m - l + 1;
	int n2 =  r - m;

	/* create temp arrays */
	//struct SimplifySwarm L[n1], R[n2];
	struct SimplifySwarm L, R;
	L.eval = new long double [n1];
	L.id = new int [n1];
	R.eval = new long double [n2];
	R.id = new int [n2];

	/* Copy data to temp arrays L[] and R[] */
	for (i = 0; i < n1; i++){
		L.eval[i] = arr->eval[l + i];
		L.id[i] = arr->id[l + i];
	}
	for (j = 0; j < n2; j++){
		R.eval[j] = arr->eval[m + 1+ j];
		R.id[j] = arr->id[m + 1+ j];
	}

	/* Merge the temp arrays back into arr[l..r]*/
	i = 0; // Initial index of first subarray
	j = 0; // Initial index of second subarray
	k = l; // Initial index of merged subarray
	while (i < n1 && j < n2) {
		if (L.eval[i] <= R.eval[j]) {
			arr->eval[k] = L.eval[i];
			arr->id[k] = L.id[i];
			i++;
		}
		else {
			arr->eval[k] = R.eval[j];
			arr->id[k] = R.id[j];
			j++;
		}
		k++;
	}

	/* Copy the remaining elements of L[], if there
       are any */
	while (i < n1) {
		arr->eval[k] = L.eval[i];
		arr->id[k] = L.id[i];
		i++;
		k++;
	}

	/* Copy the remaining elements of R[], if there
       are any */
	while (j < n2) {
		arr->eval[k] = R.eval[j];
		arr->id[k] = R.id[j];
		j++;
		k++;
	}
	//delete [] L.id;
	//delete [] R.id;
	//delete [] L.eval;
	//delete [] R.eval;
}

// Utility function to create a new tree node
Node* Utils::newNode(Node* parent, int numSubNodos, int id, long double valor) {
	Node *node = new Node;
	node->id = valor;
	node->function_eval = valor;
	node->parent= parent;
	node->numSubNodos = numSubNodos; // <<--- AQUI
	if (numSubNodos > 0)
		//node->hijo =  malloc( numSubNodos * sizeof(Nodo*) );
		node->child = new Node*;
	else
		node->child = NULL;
	return node;
}
