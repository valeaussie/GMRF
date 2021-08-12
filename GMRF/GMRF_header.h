#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED



//DEFINITIONS


extern const int nu; //dimension of the square grid
extern const int mu; //number of columns (and rows) of the precision matrix (nu*nu)
extern const int adj; //number of maximum adjacent elements in a grid

extern std::vector < double > X; //GMRF vector
extern std::vector < std::vector < double > > mat_Z; //matrix of the state of the observations

int GMRF_func(); //contains all functions
int GMRF_model(); //draws from the GMRF
int GMRF_obs(); //draws the observations from a random walk
//int GMRF_gold(); //calculates the analytical solutions

#endif
#pragma once
