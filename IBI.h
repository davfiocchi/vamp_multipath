#ifndef IBI
#define IBI

#include <stdio.h>
#include <armadillo>
#include <math.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include "VectorUtils.h"

using namespace arma;

typedef struct {
	double obs;
	int period;
	int time;
} Pt;

/*
 * Returns the IBI vector from the ODF vector
 */
arma::rowvec getIBI_vector(arma::rowvec ODF);//, float t_dif, int fs, int w_dim = 512);

/*
 * Returns the matrix out of the comb filter
 * @startPeriod: beginning index of the period
 * @endPeriod: ending index of the period 
 * @frameLength: Length of the frame
 * @rayParam: parameter for creating the weigthing vector
*/
arma::mat getCombFilterMatrix(int startPeriod, int endPeriod, int frameLength, int rayParam);

/*
 * Returns the weighting vector
 * @dimension: resulting vector length
 * @rayParam: RayLeight parameter
*/
std::vector<double> getWeightCurve(int dimension, double rayParam);

/*
* Returns the vector smoothed by a rectangular window of dimension winDim
* @vec: input vector to be smoothed
* @winDim: window dimension
*/
std::vector<double> smooth(std::vector<double> vec, int winDim);

/*
* Implementation of utils.hwr(f - utils.smooth(f, 15, 'boxcar'))
* @vec: the frame f
*/
std::vector<double> hwr(std::vector<double> vec);

/*
 * Returns the the positive lag autocorrelation vector of vec
 * @vec: input vector
*/
std::vector<double> unbiasedAutocorrelation(std::vector<double> vec);

/*
* Compute the dot product between vector and a matrix as scipy does
* @vec: input vector
* @m: matrix
*/
arma::colvec dot(arma::colvec vec, arma::mat m);

typedef struct {
	std::vector<int> indeces;
	std::vector<double> tS;
} Tsi;

/*
* Get the candidate index and the transition score of a point given the set of all points
* @p: point to evaluate
* @pts: the whole set of points
* @sigma: parameter
*/
Tsi transitionScore(Pt p, std::vector<Pt> pts, float sigma, int threshold);

/*
* Returns the range of indexes where value-threshold and value should be inserted into vec
* @value: the value to be inserted
* @vec: the vector
* @threshold: threshold
*/
std::vector<int> searchSortedRange(std::vector<int> vec, int value, int threshold);

/*
* Returns the terminating index of the score
* @backlink: backlink
* @cumscore: cumscore
*/
int findTermination(std::vector<int> backlink, std::vector<double> cumscore);

/*
* Returns the interpolated version of the path vector interpolated
* @path: vector to be interpolated
* @periodsHop: hopsize between points
* @length: number of points in the returned vector
*/
arma::rowvec interpolatePeriods(rowvec path, int periodsHop, int length);

/*
 * Returns a new vector with elements from beginning to end
 * @beginning: start value
 * @end: end value
*/
std::vector<int> range(int beginning, int end);


#endif /* IBI */
#pragma once