/*
 * utils.cpp
 *
 *  Created on: May 31, 2019
 *      Author: christian
 */
#include "utils.h"

#include <string>
#include <sstream>
#include <fstream>

using namespace std;

void Utils::loadRowVector(ifstream & fs, int columns, double * row){

	//string to storage on line read from the file
	string aline;
	string arow;

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

