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



const int max_steps = 30; //number of steps of the random walk
int steps; //total number of element observed
std::vector < std::vector < double > > mat_Z(max_steps+1, std::vector < double > (mu, 0)); //the matrix of the state of the observations
std::vector < double > P(max_steps+1, -1); //the vector containing the index of the observed elements of the grid


									  
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
	std::vector < double > O(max_steps + 1, -1); //this vector contains the indexes of all observations at all times
	O[0] = 0;
	P[0] = 0; //this vector contains the indexes of the observed elements
	for (double k = 0; k < max_steps; k++) {
		P[k + 1] = ( B[k] + nu * A[k] );
		O[k + 1] = ( B[k] + nu * A[k] );
	}
	//sort and order this vector for future use
	std::cout << "O " << std::endl;
	for (double i : O) { std::cout << i << ","; }

	sort( P.begin(), P.end() );
	P.erase(unique(P.begin(), P.end()), P.end());

	std::cout << "P " << std::endl;
	for (double i : P) { std::cout << i << ","; }

	// set the total number of element observed
	steps = P.size();
	
	//create the matrix of the state of the observations
	mat_Z[0][0] = X[0];
	for (int j = 1; j < steps; j++) {
		for (int i = 0; i < mu; i++) {
			mat_Z[j][i] = mat_Z[j-1][i];
			mat_Z[j][P[j]] = X[P[j]];
		}
	}

	std::cout << "O" << std::endl;
	for (double i : O) { std::cout << i << std::endl; }


	
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

	//Creating a dat file with the values of the vector of the observed events "vect_obs_N"
	//at the current time N calling it "vector_Z.dat"
	std::ofstream outFile1("./vector_Z.dat");
	//outFile1 << endl;
	for (double n : mat_Z[steps - 1]) {
		outFile1 << n << std::endl;
	}
	outFile1.close();

	return 0;
}