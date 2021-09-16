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



Eigen::SparseMatrix<double> Q_AA(std::vector < double > A, Eigen::SparseMatrix<double> Prec_Q);
Eigen::SparseMatrix<double> Q_AB(std::vector < double > A, Eigen::SparseMatrix<double> Prec_Q);
void F_outExp(std::vector <double> exp);

int GMRF_gold() {

	//calculate the mean of the conditional posterior

	Eigen::SparseMatrix < double > QAB = Q_AB(P, Q);

	//crete a vector of size mu of non negative integers
	std::vector < double > integers(mu); //vectors of size mu
	std::iota(std::begin(integers), std::end(integers), 0); //fill with 0, 1, ..., mu

	int size_of_obs = P.size();

	//creates a vector A with the elements of the vector "integers" that are not in P
	int size_of_A = (mu - size_of_obs);
	std::vector < double > A;
	std::remove_copy_if(integers.begin(), integers.end(), std::back_inserter(A),
		[](int const arg)
	{ return (std::find(P.begin(), P.end(), arg) != P.end()); });

	double* ptr2 = &A[0];
	Eigen::Map<Eigen::VectorXd> A_eigen(ptr2, size_of_A); //converting to eigen

	//calculate the matrice QAA
	Eigen::SparseMatrix < double > QAA = Q_AA(A, Q);

	//convert sparce matrices to dense matrices
	Eigen::MatrixXd denseQAA = Eigen::MatrixXd(QAA);
	Eigen::MatrixXd denseQAB = Eigen::MatrixXd(QAB);

	//typecast P into eigen
	std::vector < double > obs;
	for (int i : P) {
		obs.push_back(X[i]);
	}

	int obs_size = obs.size();
	double* ptr4 = &obs[0];
	Eigen::Map<Eigen::VectorXd> obs_eigen(ptr4, obs_size);
	
	
	//calculate the mean
	Eigen::VectorXd mean = -(denseQAA.inverse() * denseQAB) * obs_eigen;

	//typecast Eigen to std
	std::vector < double > mean_std(mean.size(), 0.0);
	Eigen::Map<Eigen::VectorXd>(mean_std.data(), mean.size()) = mean;


	//Create a dat files for expectations and variances
	F_outExp(mean_std);
	


	return 0;
}

//Creates a dat file with the values of the Expectations
void F_outExp(std::vector <double> exp) {
	std::ofstream outFile("./GMRF_gold_000.csv");
	outFile << std::endl;
	for (double n : exp) {
		outFile << n << std::endl;
	}
	outFile.close();
}