#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "GMRF_header.h"
#include "numeric"
#include <Eigen/Sparse>



/********FUNCTIONS*********/
void F_print_matrix(std::vector < std::vector < int > > m);
void F_print_matrix(std::vector < std::vector < double > > M);
void F_print_vector(std::vector < double > v);
void F_print_vector(std::vector < int > v);

Eigen::SparseMatrix < double > Precision(int dim_grid, int param_beta);
std::vector < double > Chol(int dim_grid, Eigen::SparseMatrix < double > Prec_Q);

using namespace std;
using namespace Eigen;

int main() {

	std::random_device rd; // create random device to seed generator
	std::mt19937 generator(rd()); // create generator with random seed

	GMRF_func(); //calls the file with all the functions
	GMRF_model();  //simulates the Gaussian Markov Random Field
	GMRF_obs(); //simulate the state of the observations

	int n = 1000; //number of particles


	Eigen::SparseMatrix<double> Q = Precision(nu, adj);

	std::vector < double > X1 = Chol(nu, Q);





	//int GMRF_gold();

	return 0;

}