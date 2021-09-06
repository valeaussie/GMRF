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
std::vector < double > Chol_and_LTsol(int dim_grid, Eigen::SparseMatrix<double> Prec_Q) {

	int dim_prec = dim_grid * dim_grid;
	std::vector < double > GMRF_vec(dim_prec, 0);

	//Sampling z from a standard multivariate normal distribution of size mu
	Eigen::MatrixXd covar(dim_prec, dim_prec);
	covar << Eigen::MatrixXd::Identity(dim_prec, dim_prec); //initialise my matrix as an identity matrix
	normal_random_variable sample{ covar };
	Eigen::VectorXd z = sample();

	//Decomposition of Q = LL^T and solution of L^T x = z
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> choleskyQ;
	//Eigen::MatrixXd L = choleskyQ.compute(Prec_Q).matrixL();
	//Eigen::MatrixXd LT = L.transpose();
	Eigen::VectorXd x = Prec_Q.triangularView<Eigen::Upper>().solve(z); //Solution of L^T x = z 

	//creating an std::vector from Eigen::VectorXd
	Eigen::Map<Eigen::VectorXd>(GMRF_vec.data(), dim_prec) = x;

	return GMRF_vec;
}

//Sampling z from a standard multivariate normal distribution of size mu
Eigen::VectorXd SampleMND(int dim_of_matrix) {
	Eigen::MatrixXd covar(dim_of_matrix, dim_of_matrix);
	covar << Eigen::MatrixXd::Identity(dim_of_matrix, dim_of_matrix); //initialise my matrix as an identity matrix
	normal_random_variable sample{ covar };
	Eigen::VectorXd z = sample();

	return z;
}

//Cholesky decomposition Q = LL^T ad solution of L x = z 
std::vector < double > Chol_and_Lsol( int dim_grid, Eigen::SparseMatrix<double> Prec_Q ) {

	int dim_prec = dim_grid * dim_grid;
	std::vector < double > GMRF_vec( dim_prec, 0 );

	Eigen::VectorXd z = SampleMND(dim_prec);

	//Sampling z from a standard multivariate normal distribution of size mu
	//Eigen::MatrixXd covar(dim_prec, dim_prec);
	//covar << Eigen::MatrixXd::Identity(dim_prec, dim_prec); //initialise my matrix as an identity matrix
	//normal_random_variable sample{ covar };
	//Eigen::VectorXd z = sample();

	//Decomposition of Q = LL^T and solution of L x = z
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> choleskyQ;
	//Eigen::MatrixXd L = choleskyQ.compute(Prec_Q).matrixL();
	Eigen::VectorXd x = Prec_Q.triangularView<Eigen::Lower>().solve(z); //Solution of L x = z 

	//creating an std::vector from Eigen::VectorXd
	Eigen::Map<Eigen::VectorXd>(GMRF_vec.data(), dim_prec) = x;

	return GMRF_vec;
}

//creates a patition Q_AA of the (sparse) precision matrix Q for the elements in the set A
Eigen::SparseMatrix<double> Q_AA(std::vector < int > A, Eigen::SparseMatrix<double> Prec_Q) {
	int dim_A = A.size();
	Eigen::SparseMatrix<double> QAA(dim_A, dim_A);
	QAA.reserve(Eigen::VectorXi::Constant(dim_A, dim_A));
	for (int i = 0; i < dim_A; i++) {
		for (int j = 0; j < dim_A; j++) {
			QAA.coeffRef(i, j) = Prec_Q.coeffRef(A[i], A[j]);
		}
	}
	return QAA;
}

//creates a patition Q_AB of the (sparse) precision matrix Q for the elements in the set A and in the set B complment of A
Eigen::SparseMatrix<double> Q_AB(std::vector < int > B, Eigen::SparseMatrix<double> Prec_Q) {
	//creates a vector v of increasing integers from 1 to the square root of the precision matrix
	int dim_B = B.size();
	int max_val = sqrt(Prec_Q.size());
	int dim_A = max_val - dim_B;
	std::vector < int > A;
	std::vector < int > v(max_val, 0);
	for (int i = 0; i < max_val; i++) {
		v[i] = i;
	}
	//creates a vector B with the elements of v that are not in A
	std::remove_copy_if(v.begin(), v.end(), std::back_inserter(A),
		[&B](const int& arg)
	{ return (std::find(B.begin(), B.end(), arg) != B.end()); });
	//creates the partition Q_AB
	Eigen::SparseMatrix<double> QAB(dim_A, dim_B);
	QAB.reserve(Eigen::VectorXi::Constant(dim_B, 9));
	for (int i = 0; i < dim_B; i++) {
		for (int j = 0; j < dim_A; j++) {
			QAB.coeffRef(j, i) = Prec_Q.coeffRef(A[j], B[i]);
		}
	}
	return QAB;
}



