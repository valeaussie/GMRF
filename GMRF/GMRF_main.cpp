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
Eigen::VectorXd Chol_and_LTsol_eigen(int dim_grid, Eigen::SparseMatrix<double> Prec_Q);
std::vector < double > Chol_and_Lsol(int dim_grid, Eigen::SparseMatrix<double> Prec_Q);
Eigen::VectorXd Chol_and_Lsol_eigen(int dim_grid, Eigen::SparseMatrix<double> Prec_Q);
Eigen::SparseMatrix<double> Q_AA(std::vector < double > A, Eigen::SparseMatrix<double> Prec_Q);
Eigen::SparseMatrix<double> Q_AB(std::vector < double > A, Eigen::SparseMatrix<double> Prec_Q);
Eigen::VectorXd SampleMND(int dim_of_matrix);



int main() {

	GMRF_func(); //calls the file with all the functions
	GMRF_model();  //simulates the Gaussian Markov Random Field
	GMRF_obs(); //simulate the state of the observations

	int n = 10000; //number of particles


	std::vector < std::vector < std::vector < double > > > resampled(steps, 
		std::vector < std::vector < double > >(mu, std::vector < double >(n, 0.0)));
	//define and initialise the containter for the unnormalised weights
	std::vector < std::vector < double > > un_weights(n, std::vector < double >(steps, 0.0));
	//define and initialise the container for the normalised weights
	std::vector < std::vector < double > > weights(n, std::vector < double >(steps, 0.0));



	std::random_device rd;
	std::mt19937 generator(rd());

	//convert Q into a dense matrix
	Eigen::MatrixXd denseQ = Eigen::MatrixXd(Q);


	//draw the process for n particles
	Eigen::MatrixXd mat_sim(n, mu);
	for (int i = 0; i < n; i++) {
		Eigen::VectorXd sim = Chol_and_Lsol_eigen(nu, Q); 
		Eigen::RowVectorXd row_sim = sim.transpose();
		mat_sim.row(i) = row_sim;
	}


	//initialize:
	for (int i = 0; i < n; i++) {
		un_weights[i][0] = 1; //set equal weights at time 1
		weights[i][0] = 1.0 / n; //normalise the weights
		//for (int i = 0; i < steps; i++) {
			resampled[0][0][i] = X[0];
			for (int k = 1; k < mu; k++) {
				resampled[0][k][i] = mat_sim(i, k);
			}
		//}
	}

	//typecast X into eigen
	int X_size = X.size();
	double* ptr3 = &X[0];
	Eigen::Map<Eigen::VectorXd> X_eigen(ptr3, X_size);
	
	/*********Iterate********/
	for (int j = 1; j < steps; j++) {
		int o = P[j]; // index of observed element at time j
		//vector with the indexes of the observed elements up to now (j-1)
		std::vector < double > obs_now;
		for (int l = 0; l < j; l++) {
			obs_now.push_back(P[l]);
		}

		//populate the resample vector from previous time
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < mu; k++) {
				resampled[j][k][i] = resampled[j - 1][k][i];
			}
		}

		/*********make correction********/
		for (int i = 0; i < n; i++) {
			resampled[j][o][i] = X[o];
		}

		/*********calculate the weights********/
		for (int i = 0; i < n; i++) {
			//fill vector sim
			Eigen::VectorXd sim(mu);
			for (int k = 0; k < mu; k++) {
				sim(k) = resampled[j-1][k][i]; //populate Eigen::vector sim
			}
			
			//calculate weights for particle i at time j
			double w = exp( (denseQ.row(o) * sim - sim(o) ) * ( sim(o) - X_eigen(o) )
				+ (1 / 2) * ( (pow(sim(o), 2) - pow(X_eigen(o), 2) ) ) );
			un_weights[i][j] = w; //insert in matrix weigths
		}

		//normalise the weights
		double sum_of_weights{ 0 };
		for (int i = 0; i < n; i++) {
			sum_of_weights = sum_of_weights + un_weights[i][j];
		}
		for (int i = 0; i < n; i++) {
			weights[i][j] = un_weights[i][j] / sum_of_weights;
		}

		//resampling
		std::vector < double > drawing_vector(n, 0.0);
		for (int i = 0; i < n; i++) {
			drawing_vector[i] = weights[i][j];
		}

		double index_resampled;
		std::vector < int > vec_index;
		for (int i = 0; i < n; i++) {
			std::discrete_distribution < int > discrete(drawing_vector.begin(), drawing_vector.end());
			index_resampled = discrete(generator);
			vec_index.push_back(index_resampled);
		}
		std::vector < std::vector < double > > newmatrix(mu, std::vector < double >(n, 0.0));
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < mu; k++) {
				newmatrix[k][i] = resampled[j][k][vec_index[i]];
			}
		}
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < mu; k++) {
				resampled[j][k][i] = newmatrix[k][i];
			}
		}
	}

	std::ofstream outFile("./resampled_000.csv");
	outFile << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < mu; j++) {
			outFile << resampled[steps-1][j][i] << ",";
		}
		outFile << std::endl;
	}
	outFile.close();

	std::ofstream outFile3("./resampled_001.csv");
	outFile3 << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < mu; j++) {
			outFile3 << resampled[1][j][i] << ",";
		}
		outFile3 << std::endl;
	}
	outFile3.close();


	std::ofstream outFile2("./weights.csv");
	outFile2 << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < steps; j++) {
			outFile2 << weights[i][j] << ",";
		}
		outFile2 << std::endl;
	}
	outFile2.close();
	


	/*

			//this is the code to calculate the conditional probabilities


			Eigen::SparseMatrix < double > QAB = Q_AB(P, Q);

			obs.erase(
				std::remove(obs.begin(), obs.end(), 0),
				obs.end());

			//converting obs (std into eigen)
			int size_of_obs = obs.size();
			double* ptr1 = &obs[0];
			Eigen::Map<Eigen::VectorXd> obs_eigen(ptr1, size_of_obs);


			//calculate the mean of the canonical distribution
			Eigen::VectorXd QAB_times_XB = - (QAB * obs_eigen);

			//crete a vector of size mu of non negative integers
			std::vector < double > integers(mu); //vectors of size mu
			std::iota(std::begin(integers), std::end(integers), 0); //fill with 0, 1, ..., mu
			//creates a vector A with the elements of the vector "integers" that are not in P
			int size_of_A = (mu - size_of_obs);
			std::vector < double > A;
			std::remove_copy_if ( integers.begin(), integers.end(), std::back_inserter(A),
				[](int const arg )
			{ return (std::find(obs.begin(), obs.end(), arg) != obs.end()); });

			double* ptr2 = &A[0];
			Eigen::Map<Eigen::VectorXd> A_eigen(ptr2, size_of_A); //converting to eigen

			//calculate the matrix QAA
			Eigen::SparseMatrix < double > QAA = Q_AA( A, Q );

			Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower> choleskyQAA_Lower;
			Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> choleskyQAA_Upper;
			Eigen::VectorXd w1 = choleskyQAA_Lower.compute(QAA).solve(QAB_times_XB); //Solution of L w1 = b
			Eigen::VectorXd w2 = choleskyQAA_Upper.compute(QAA).solve(w1); //Solution of L^T w2 = w1
			Eigen::VectorXd z1 = SampleMND(size_of_A); //Sample from a normal distribution of size "size_of_a"
			Eigen::VectorXd w3 = choleskyQAA_Upper.compute(QAA).solve(z1); //Solution of L^T w3 = z1
			Eigen::VectorXd result = w2 + w3;//compute w2 + w3

			*/

	//create a .csv file containing the parameters
	std::ofstream outparam("./parameters.csv");
	outparam << "mu" << "," << "max_steps" << "," << "part" << std::endl;
	outparam << mu << "," << max_steps << "," << n << "," << std::endl;
	outparam.close();



	GMRF_gold();

	return 0;

}

