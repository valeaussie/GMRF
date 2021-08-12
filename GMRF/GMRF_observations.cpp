#define _CRT_SECURE_NO_WARNINGS

#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include "GMRF_header.h"
#include <cmath>

/********FUNCTIONS*********/
void F_print_matrix(std::vector < std::vector < int > > m);
void F_print_matrix(std::vector < std::vector < double > > M);
void F_print_vector(std::vector < double > v);
void F_print_vector(std::vector < int > v);



int max_steps = 10; //number of steps of the random walk
std::vector < std::vector < double > > mat_Z(max_steps, std::vector <double > (mu, 0)); //the matrix of the state of the observations
std::vector < int > O(max_steps, -1); //the vector containing the index of the observed elements of the grid


									  
									  //2D random walk
int GMRF_obs() {

	std::random_device rd; // create random device to seed generator
	std::mt19937 gen(rd()); // create generator with random seed
	std::uniform_real_distribution < double > uni(0., 1); // init uniform dist on (0,1]

	int x, y; //positions
	int current_step = 0;


	std::vector < double > A;
	std::vector < double > B;
	// Domain is x in [0,nu], y in [0,nu]
	x = 0;
	y = 0; // start at top left corner

	while (current_step < max_steps) {
		// Choose x or y direction randomly
		if (uni(gen) < 0.5) {
			// move in the x direction 
			if (uni(gen) < 0.5) {
				// step left
				if (x != 0) {
					x = x - 1;
				}
				else {
					x = x + 1; // reflection
				}
			}
			else {
				// step right
				if (x != nu-1) {
					x = x + 1;
				}
				else {
					x = x - 1; // reflection
				}
			}
		}
		else {
			// move in the y direction
			if (uni(gen) < 0.5) {
				// step down
				if (y != 0) {
					y = y - 1;
				}
				else {
					y = y + 1; // reflection
				}
			}
			else {
				// step up
				if (y != nu-1) {
					y = y + 1;
				}
				else {
					y = y - 1; // reflection
				}
			}
		}
		// Update
		current_step++;
		A.push_back(x);
		B.push_back(y);
	}


	//translate the x and y coordinates to the vector O of the indexes
	for (int k = 0; k < max_steps; k++) {
		O[k] = B[k] + nu * A[k];
	}


	//create the matrix of the state of the observations
	mat_Z[0][0] = X[0];
	for (int j = 1; j < max_steps; j++) {
		for (int i = 0; i < mu; i++) {
			mat_Z[j][i] = mat_Z[j-1][i];
			mat_Z[j][O[j]] = X[O[j]];
		}
	}



	//creates a dat file with the values of the matrix of the state of the observations mat_Z 
	//called "vector_Z.dat"
	std::ofstream outFile("./matrix_Z.dat");
	for (std::vector < double >  v : mat_Z) {
		for (double n : v) {
			outFile << n << " ";
		}
		outFile << std::endl;
	}
	outFile.close();


	return 0;
}