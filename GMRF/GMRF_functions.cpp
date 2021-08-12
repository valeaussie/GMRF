#define _CRT_SECURE_NO_WARNINGS

#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "GMRF_header.h"
#include "numeric"
#include <Eigen/Dense>
#include <Eigen/Sparse>


using namespace std;

void F_print_matrix(vector < vector < int > > m);
void F_print_matrix(vector < vector < double > > M);
void F_print_vector(vector < double > v);
void F_print_vector(vector < int > v);

int GMRF_func() {

	return 0;
}


//*********** GENERIC FUNCTIONS ***************

//Prints a matrix of int
void F_print_matrix(vector < vector < int > > m) {
	for (const vector < int > v : m) {
		for (int x : v) cout << x << ' ';
		cout << endl;
	}
}
//Prints a matrix of signed doubles
void F_print_matrix(vector < vector < double > > M) {
	for (const vector < double > v : M) {
		for (double x : v) cout << x << ' ';
		cout << endl;
	}
}
//Prints a vector of doubles
void F_print_vector(vector < double > v) {
	for (const double x : v) cout << x << ' ';
	cout << endl;
}

//Prints a vector of doubles
void F_print_vector(vector < int > v) {
	for (const size_t x : v) cout << x << ' ';
	cout << endl;
}


//*********** GMRF SPECIFIC FUNCTIONS ************

//creates the precision matrix (sparse) for any value of the dimension of the grid
Eigen::SparseMatrix<double> Precision(int dim_grid, int adj_var) {

	float param_beta = 1.0 / adj_var;
	param_beta = floor(param_beta * 100.0) / 100.0; //this ensures that the precision matrix Q is diagonal dominant
	//this because Q must be positive definite

	int dim_prec = dim_grid * dim_grid;
	Eigen::SparseMatrix<double> Q(dim_prec, dim_prec);
	Q.reserve(Eigen::VectorXi::Constant(dim_prec, 9));
	for (int i = 0; i < dim_prec; i++) {
		for (int j = 0; j < dim_prec; j++) {
			Q.coeffRef(i, i) = 1;
		}
	}
	for (int i = 0; i < dim_prec; i++) {
		for (int j = 0; j < dim_prec; j++) {
			for (int k = 1; k < dim_grid - 1; k++) {
				for (int i = (dim_grid*k) + 1; i < (k + 1)*dim_grid - 1; i++) {
					for (int j = 0; j < dim_prec; j++) {
						Q.coeffRef(i, i - 1) = -param_beta;
						Q.coeffRef(i, i + 1) = -param_beta;
						Q.coeffRef(i, i - dim_grid) = -param_beta;
						Q.coeffRef(i, i - dim_grid - 1) = -param_beta;
						Q.coeffRef(i, i - dim_grid + 1) = -param_beta;
						Q.coeffRef(i, i + dim_grid) = -param_beta;
						Q.coeffRef(i, i + dim_grid - 1) = -param_beta;
						Q.coeffRef(i, i + dim_grid + 1) = -param_beta;
					}
				}
			}
		}
	}
	for (int i = 0; i < dim_prec; i++) {
		for (int j = 0; j < dim_prec; j++) {
			if (i == 0) {
				Q.coeffRef(i, i + 1) = -param_beta;
				Q.coeffRef(i, i + dim_grid) = -param_beta;
				Q.coeffRef(i, i + dim_grid + 1) = -param_beta;
			}
			if (i == dim_grid - 1) {
				Q.coeffRef(i, i - 1) = -param_beta;
				Q.coeffRef(i, i + dim_grid) = -param_beta;
				Q.coeffRef(i, i + dim_grid - 1) = -param_beta;
			}
			if (i == dim_prec - 1) {
				Q.coeffRef(i, i - 1) = -param_beta;
				Q.coeffRef(i, i - dim_grid) = -param_beta;
				Q.coeffRef(i, i - dim_grid - 1) = -param_beta;
			}
			if (i == dim_prec - dim_grid) {
				Q.coeffRef(i, i + 1) = -param_beta;
				Q.coeffRef(i, i - dim_grid) = -param_beta;
				Q.coeffRef(i, i - dim_grid + 1) = -param_beta;
			}
			for (int k = 1; k < dim_grid - 1; k++) {
				if (i == k) {
					Q.coeffRef(i, i - 1) = -param_beta;
					Q.coeffRef(i, i + 1) = -param_beta;
					Q.coeffRef(i, i + dim_grid) = -param_beta;
					Q.coeffRef(i, i + dim_grid - 1) = -param_beta;
					Q.coeffRef(i, i + dim_grid + 1) = -param_beta;
				}
				if (i == (k + 1)*dim_grid - 1) {
					Q.coeffRef(i, i - 1) = -param_beta;
					Q.coeffRef(i, i - dim_grid) = -param_beta;
					Q.coeffRef(i, i - dim_grid - 1) = -param_beta;
					Q.coeffRef(i, i + dim_grid) = -param_beta;
					Q.coeffRef(i, i + dim_grid - 1) = -param_beta;
				}
				if (i == dim_prec - k - 1) {
					Q.coeffRef(i, i - 1) = -param_beta;
					Q.coeffRef(i, i + 1) = -param_beta;
					Q.coeffRef(i, i - dim_grid) = -param_beta;
					Q.coeffRef(i, i - dim_grid - 1) = -param_beta;
					Q.coeffRef(i, i - dim_grid + 1) = -param_beta;
				}
				if (i == dim_prec - (k + 1)*dim_grid) {
					Q.coeffRef(i, i + 1) = -param_beta;
					Q.coeffRef(i, i - dim_grid) = -param_beta;
					Q.coeffRef(i, i - dim_grid + 1) = -param_beta;
					Q.coeffRef(i, i + dim_grid) = -param_beta;
					Q.coeffRef(i, i + dim_grid + 1) = -param_beta;
				}
			}
		}
	}
	return Q;
}


//Cholesky decomposition Q = LL^T ad solution of L^T x = z 
std::vector < double > Chol(int dim_grid, Eigen::SparseMatrix<double> Prec_Q) {

	int dim_prec = dim_grid * dim_grid;
	std::vector < double > GMRF_vec(dim_prec, 0);
	//class to sample from a multivariate normal distribution with mean zero
	struct normal_random_variable
	{
		normal_random_variable(Eigen::MatrixXd const& covar)
			: normal_random_variable(Eigen::VectorXd::Zero(covar.rows()), covar)
		{}

		normal_random_variable(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar)
			: mean(mean)
		{
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
			transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
		}

		Eigen::VectorXd mean;
		Eigen::MatrixXd transform;

		Eigen::VectorXd operator()() const
		{
			static std::mt19937 gen{ std::random_device{}() };
			static std::normal_distribution<> dist;

			return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](auto x) { return dist(gen); });
		}
	};

	//Sampling z from a standard multivariate normal distribution of size mu
	Eigen::MatrixXd covar(dim_prec, dim_prec);
	covar << Eigen::MatrixXd::Identity(dim_prec, dim_prec); //initialise my matris as an identity matrix
	normal_random_variable sample{ covar };
	Eigen::VectorXd z = sample();


	//Decomposition of Q = LL^T and solution of L^T x = z 
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
	solver.compute(Prec_Q).solve(z); //decomposition

	Eigen::VectorXd x = solver.solve(z); //Solution of L^T x = z 

	//creating an std::vector from Eigen::VectorXd
	Eigen::Map<Eigen::VectorXd>(GMRF_vec.data(), dim_prec) = x;

	return GMRF_vec;
}



