//
//  BeatPeriod.hpp
//  BeatTracker
//
//  Created by Alessandro Ragano on 28/02/2017.
//  Copyright © 2017 Alessandro Ragano. All rights reserved.
//

#ifndef BeatPeriod_hpp
#define BeatPeriod_hpp
#include </usr/local/include/armadillo>
#include <stdio.h>

using namespace arma;

/*
* Returns the IBI vector given the ODF
* @odf: vector of ODF values
*/
rowvec getBeatPeriod(rowvec odf);

/*
* Partitioned odf vector into overlapped frames. Each column represents a frame
* @odf: vector of ODF values
* @frame_length: frame length dimension
* @hop_size: hopsize dimension
*/
mat overlapping_odf(rowvec odf, int frame_length, int hop_size);

/*
* Compute the adaptive moving mean on each odf overlapped frame
* @frames_odf: overlapped frame of the ODF
*/
mat adaptive_moving_mean(mat frames_odf);

/*
* Compute the half-way rectify given the ovelapped ODF and its moving mean
* @overlapped_frames_odf: overlapped frame of the ODF
* @adaptive_moving_mean: ODF moving mean
*/
mat hwr(mat overlapped_frames_odf, mat adaptive_moving_mean);

/*
* Computes the unbiased autocorrelation function
* @hwr_frames: half way rectified frames of the ODF
*/
mat unbiased_autocorrelation(mat hwr_frames);

/*
* Compute the output of a comb template, given its length and tau
* @frame_length: output dimension
* @tau: tau parameter
*/
colvec comb_template(int frame_length, int tau);

/*
* Compute the shift invariant comb filterbank
* @frame_length: length of the frame
*/
mat comb_filterbank(int frame_length);

/*
* Compute the product of the autocorrelation function with the comb filterbank(for each odf frame)
* @autocorrelation: autocorrelation function
* @frame_length: frame length dimension
*/
mat output(mat autocorrelation, int frame_length);

/*
* Compute the beat period as the tau_max for each odf frame
* @yg: output of the autocorrelation function and the combfilterbank
*/
rowvec beat_period(mat yg);

/*
* Interpolate the IBI values to obtain the same rate as the ODF
* @ibi: vector of IBI values
* @N: dimension of interpolation
*/
rowvec interpolate(rowvec ibi, int N);


#endif /* BeatPeriod_hpp */
