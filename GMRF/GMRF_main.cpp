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



/********FUNCTIONS*********/
void F_print_matrix(std::vector < std::vector < int > > m);
void F_print_matrix(std::vector < std::vector < double > > M);
void F_print_vector(std::vector < double > v);
void F_print_vector(std::vector < int > v);

Eigen::SparseMatrix < double > Precision(int dim_grid, int param_beta);
std::vector < double > Chol_and_LTsol(int dim_grid, Eigen::SparseMatrix < double > Prec_Q);
Eigen::SparseMatrix<double> Q_AA(std::vector < int > A, Eigen::SparseMatrix<double> Prec_Q);
Eigen::SparseMatrix<double> Q_AB(std::vector < int > A, Eigen::SparseMatrix<double> Prec_Q);
//struct normal_random_variable;

using namespace std;
using namespace Eigen;

int main() {

	std::random_device rd; // create random device to seed generator
	std::mt19937 generator(rd()); // create generator with random seed

	GMRF_func(); //calls the file with all the functions
	GMRF_model();  //simulates the Gaussian Markov Random Field
	GMRF_obs(); //simulate the state of the observations

	int n = 10; //number of particles

	//define the container for the sampled events and the sampled observations (0s and 1s)
	vector < vector < vector < double > > > sample(max_steps + 1, vector < vector < double > >(mu, vector < double >(n, 0.0)));
	//define the container for the new sampled events and the new sampled observations (0s and 1s)
	vector < vector < vector < double > > > corr_sample(max_steps + 1, vector < vector < double > >(mu, vector < double >(n, 0.0)));
	vector < vector < vector < double > > > resampled(max_steps + 1, vector < vector < double > >(mu, vector < double >(n, 0.0)));
	//define and initialise the containter for the unnormalised weights
	vector < vector < double > > un_weights(n, vector < double >(max_steps + 1, 0.0));
	//define and initialise the container for the normalised weights
	vector < vector < double > > weights(n, vector < double >(max_steps + 1, 0.0));

	
	//initialize:
	for (int i = 0; i < n; i++) {
		un_weights[i][0] = 1;
		weights[i][0] = 1.0 / n;
		for (int j = 0; j < max_steps + 1; j++) {
			sample[j][0][i] = X[0];
			resampled[j][0][i] = X[0];
			corr_sample[j][0][i] = X[0];
		}
	}

	
	//this needs to be done only once
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> choleskyQ;
	Eigen::SparseMatrix<double> L = choleskyQ.compute(Q).matrixL();
	Eigen::SparseMatrix<double> LT = L.transpose();

	
	//std::cout << QAB << std::endl;
	
	//Iterate	
	for (int j = 1; j < max_steps+1; j++) {
		for (int i = 0; i < n; i++) {
			vector < double > sum_vec;
			for (int k = 0; k < j; k++) {
				//sample[j][k][i] = resampled[j - 1][k][i];
				//corr_sample[j][k][i] = resampled[j - 1][k][i];
			}
			//vector < double > obs = mat_Z[j];
			Eigen::SparseMatrix< double > QAB = Q_AB(P, Q);
			//Eigen::SparseMatrix< double > QAB = Q_AB(O, Q);
			//Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
			//Eigen::VectorXi b = QAB * O;
			/*
			normal_distribution < double > normalDist1((phi * (sample[j][j - 1][i])), sigmasq);
			sample[j][j][i] = normalDist1(generator);
			corr_sample[j][j][i] = sample[j][j][i];*/
			//make the corrections
			/*
			for (int k = 0; k < j + 1; k++) {
				if (mat_Z[j][k] == 1) { corr_sample[j][k][i] = X[k]; }
			}
			//calculate the weights
			//this condition ensures that we are only calculating the partial weights that are not 1
			for (int k = 1; k < j + 1; k++) {
				if ((corr_sample[j][k][i] != sample[j][k][i]) || (corr_sample[j][k - 1][i] != sample[j][k - 1][i])) {
					double num_arg = pow((corr_sample[j][k][i] - phi * corr_sample[j][k - 1][i]), 2);
					double den_arg = pow((sample[j][k][i] - phi * sample[j][k - 1][i]), 2);
					double log_elem = (1 / (2 * sigmasq)) * (num_arg - den_arg);
					sum_vec.push_back(log_elem);
				}
				else { sum_vec.push_back(0); }
			}
			double sum_of_logs = accumulate(sum_vec.begin(), sum_vec.end(), 0.0);
			double W = exp(sum_of_logs);
			un_weights[i][j] = W;*/
		}
			/*
		//normalise the weights
		double sum_of_weights{ 0 };
		for (int i = 0; i < n; i++) {
			sum_of_weights = sum_of_weights + un_weights[i][j];
		}
		for (int i = 0; i < n; i++) {
			weights[i][j] = un_weights[i][j] / sum_of_weights;
		}
		//resampling
		vector < double > drawing_vector(n, 0.0);
		for (int i = 0; i < n; i++) {
			drawing_vector[i] = weights[i][j];
		}
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < j + 1; k++) {
				resampled[j][k][i] = corr_sample[j][k][i];
			}
		}
		for (int i = 0; i < n; i++) {
			discrete_distribution < int > discrete(drawing_vector.begin(), drawing_vector.end());
			resampled[j][j][discrete(generator)];
		}*/
	}

	











	//int GMRF_gold();

	return 0;

}