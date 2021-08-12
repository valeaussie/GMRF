#define _CRT_SECURE_NO_WARNINGS

#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "GMRF_header.h"
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>




const int nu = 4; //dimension of the square grid
const int mu = nu * nu; //number of columns (and rows) of the precision matrix
const int adj = 8; //the number of maximum adjacent elements in a grid

Eigen::SparseMatrix<double> Precision(int dim_grid, int param_beta);
std::vector < double > Chol(int dim_grid, Eigen::SparseMatrix < double > Prec_Q);

std::vector < double > X(mu); //the GMRF vector


int GMRF_model() {

	GMRF_func();//calls the file with all the functions
	
	//create the precision matrix
	Eigen::SparseMatrix<double> Q = Precision(nu, adj);

	//Cholesky decomposition Q = LL^T ad solution of L^T x = z 
	X = Chol(nu, Q);


	//creates a dat file with the values of X called "GMRF_vector_X.dat"
	std::ofstream outFile("./GMRF_vector_X.dat");
	for (double n : X) {
		outFile << n << std::endl;
	}
	outFile.close();

	
	return 0;
}