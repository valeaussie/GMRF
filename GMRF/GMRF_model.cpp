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
#include <Eigen/Eigenvalues>




const int nu = 4; //dimension of the square grid
const int mu = nu * nu; //number of columns (or rows) of the precision matrix
const int adj = 4; //the number of maximum adjacent elements in a grid (the degree)
Eigen::SparseMatrix < double > Q; //the precision matrix
std::vector < double > X; //the GMRF vector



Eigen::SparseMatrix < double > Precision( int dim_grid, int param_beta );
std::vector < double > Chol_and_LTsol( int dim_grid, Eigen::SparseMatrix < double > Prec_Q );


int GMRF_model() {

	GMRF_func();//calls the file with all the functions
	
	//create the precision matrix
	Q = Precision(nu, adj);

	//Cholesky decomposition Q = LL^T ad solution of L^T x = z 
	X = Chol_and_LTsol(nu, Q);


	//creates a dat file with the values of X called "GMRF_vector_X.dat"
	std::ofstream outFile("./GMRF_vector_X.dat");
	for (double n : X) {
		outFile << n << std::endl;
	}
	outFile.close();

	Eigen::MatrixXd denseQ = Eigen::MatrixXd(Q);
	Eigen::EigenSolver < Eigen::MatrixXd > ei;
	ei.compute(denseQ);
	std::cout << "The eigenvalues of Q are:" << std::endl;
	std::cout << ei.eigenvalues() << std::endl;

	std::cout << Q << std::endl;

	
	return 0;
}