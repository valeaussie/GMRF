#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED



//DEFINITIONS

#include <Eigen/Dense>
#include <Eigen/Sparse>


extern const int nu; //dimension of the square grid
extern const int mu; //number of columns (and rows) of the precision matrix (nu*nu)
extern const int adj; //number of maximum adjacent elements in a grid
extern int steps; //number of steps of the random walk
extern const int max_steps; //number of max steps of the random walk
extern Eigen::SparseMatrix<double> Q; //Precision matrix
extern std::vector < double > X; //GMRF vector
extern std::vector < double > P; //the vector containing the index of the observed elements of the grid
extern std::vector < double > O; //the vector containing the index of the max observed elements of the grid
extern std::vector < std::vector < double > > mat_Z; //matrix of the state of the observations


int GMRF_func(); //contains all functions
int GMRF_model(); //draws from the GMRF
int GMRF_obs(); //draws the observations from a random walk
int GMRF_gold(); //calculates the analytical solutions

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

#endif
#pragma once
