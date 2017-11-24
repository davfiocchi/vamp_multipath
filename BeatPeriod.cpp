//
//  BeatPeriod.cpp
//  BeatTracker
//
//  Created by Alessandro Ragano on 28/02/2017.
//  Copyright © 2017 Alessandro Ragano. All rights reserved.
//

#include "BeatPeriod.hpp"

using namespace arma;
using namespace std;

//params definition
int rgFrameLength = 512;
int rgHoplength = 128;

int Q = 16;	//samples for the adaptive moving mean

int startPeriod = 0;	//initial value for tau
int endPeriod = 128;	//final value for tau
int interval = endPeriod - startPeriod;

int beta = 43;	//rayleigh parameter


rowvec getBeatPeriod(rowvec odf)
{
	printf("---...GetIBI_vector begin!---\n");
	rowvec ibi;
	//not enough ODF samples to compute the IBI
	if (odf.n_elem <= 0)
		return rowvec((uword)0, fill::zeros);

	mat frames_odf = overlapping_odf(odf, rgFrameLength, rgHoplength);
	mat frames_odf_mean = adaptive_moving_mean(frames_odf);
	mat half_wave_odf = hwr(frames_odf, frames_odf_mean);
	mat autocorrelation = unbiased_autocorrelation(half_wave_odf);
	mat yg = output(autocorrelation, rgFrameLength);
	ibi = beat_period(yg);
	rowvec ibi_interpolated = interpolate(ibi, int(odf.n_elem));
	return ibi_interpolated;
}

// Partitioned odf function in overlapped frames. Each column represents a frame
mat overlapping_odf(rowvec odf, int frame_length, int hop_size)
{
	int n_frames = 1 + round((odf.n_elem - frame_length) / hop_size);
	colvec odf_t = odf.t();
	mat matrix_frames = mat(frame_length, n_frames);
	
	// Compute all the frames except the last one
	for (int i = 0; i<n_frames - 1; i++)
	{
		matrix_frames.col(i) = odf_t.subvec(i*(hop_size - 1), (frame_length - 1) + i*(hop_size - 1));
	}

	// Compute the last one with some samples plus zero padding
	if (frame_length + n_frames*hop_size>odf.n_elem)
	{
		for (int i = 0; i<frame_length; i++)
		{
			if ((n_frames - 1)*hop_size + i<odf.n_elem)
			{
				matrix_frames(i, n_frames - 1) = odf_t((n_frames - 1)*hop_size + i);
			}
			else
			{
				matrix_frames(i, n_frames - 1) = 0;
			}
		}
	}

	return matrix_frames;
}

// Compute the adaptive moving mean on each odf overlapped frame
mat adaptive_moving_mean(mat overlapped_frames_odf)
{
	size_t n_frames = overlapped_frames_odf.n_cols;
	size_t frame_length = overlapped_frames_odf.n_rows;
	mat matrix_frames = mat(frame_length, n_frames);

	for (int i = 0; i<n_frames; i++)
	{
		colvec frame = overlapped_frames_odf.col(i);
		for (int m = 0; m<frame_length; m++)
		{
			int start = m - Q / 2;
			int end = m + Q / 2;
			if (start <= 0)
				start = 0;
			if (end >= frame_length)
				end = frame_length - 1;
			matrix_frames(m, i) = mean(frame.subvec(start, end));
		}
	}

	return matrix_frames;
}

// Compute the half-way rectify
mat hwr(mat overlapped_frames_odf, mat adaptive_moving_mean)
{
	size_t n_frames = overlapped_frames_odf.n_cols;
	size_t frame_length = overlapped_frames_odf.n_rows;
	mat matrix_frames = mat(frame_length, n_frames);

	for (int i = 0; i<n_frames; i++)
	{
		colvec frame = overlapped_frames_odf.col(i) - adaptive_moving_mean.col(i);
		colvec tmp = frame + abs(frame);
		for (int m = 0; m<frame_length; m++)
		{
			if (tmp(m) != 0) {
				tmp(m) = tmp(m) / 2.0;
			}
		}
		matrix_frames.col(i) = tmp;
	}

	return matrix_frames;
}

// Compute the autocorrelation function
mat unbiased_autocorrelation(mat hwr_frames)
{
	size_t n_frames = hwr_frames.n_cols;
	size_t frame_length = hwr_frames.n_rows;
	mat autocorrelation = mat(frame_length, n_frames);

	for (int i = 0; i<n_frames; i++)
	{
		colvec frame = hwr_frames.col(i);
		for (int l = 0; l<frame_length; l++)
		{
			colvec numerator = colvec(frame_length);
			for (int m = 0; m<frame_length; m++)
			{
				if (l<m)
				{
					numerator(m) = frame(m)*frame(m - l);
				}
				else
				{
					numerator(m) = 0;
				}
			}

			autocorrelation(l, i) = sum(numerator) / double(abs(int(l - frame_length)));
		}
	}

	return autocorrelation;
}

// Compute the comb template
colvec comb_template(int frame_length, int tau)
{
	colvec lambda_tau = zeros(frame_length);
	int l = 0;

	for (l = 0; l<frame_length; l++)
	{
		for (int p = 1; p <= 4; p++)
		{
			for (int v = 1 - p; v <= p - 1; v++)
			{
				if ((l - tau*p + v >= 0) &&
					(p*tau + v >= 0) && (p*tau + v < frame_length))
				{
					//lambda_tau(l) += 1.0 / (2 * p - 1);
					lambda_tau(p*tau + v) = 1.0 / double(2 * p - 1);
				}
			}
		}
	}

	return lambda_tau;
}

// Compute the shift invariant comb filterbank
mat comb_filterbank(int frame_length)
{
	mat comb_filter = mat(frame_length, interval + 1);
	
	//calculate the comb filter output for each tau
	for (double tau = startPeriod; tau <= endPeriod; tau++)
	{
		//weigthing factor
		float weight_g = (((tau) / pow(beta, 2))*exp(-(pow((tau), 2) / (2 * pow(beta, 2)))));
		
		//comb template
		colvec cTemplate = comb_template(frame_length, tau);
		
		comb_filter.col(tau - startPeriod) = weight_g*cTemplate;
	}

	return comb_filter;
}

// Compute the product of the autocorrelation function with the comb filterbank(for each odf frame)
mat output(mat autocorrelation, int frame_length)
{
	int n_frames = autocorrelation.n_cols;
	mat yg = mat(interval, n_frames);
	mat fg = comb_filterbank(frame_length);
	
	//for each frame of the autocorrelation
	for (int i = 0; i<n_frames; i++)
	{
		colvec frame = autocorrelation.col(i);
		for (int tau = startPeriod; tau<endPeriod; tau++)
		{
			double product = 0;
			for (int l = 0; l<frame_length; l++)
			{
				product += frame(l)*fg(l, tau-startPeriod);
			}
			yg(tau-startPeriod, i) = product;
		}
	}

	return yg;
}

// Compute the beat period as the tau_max for each odf frame
rowvec beat_period(mat yg)
{
	int n_frames = yg.n_cols;
	colvec period = colvec(n_frames);

	//for each frame
	for (int i = 0; i<n_frames; i++)
	{
		colvec frame = yg.col(i);

		//extract the maximum index
		period(i) = frame.index_max();
	}

	return period.t();
}

rowvec interpolate(rowvec ibi, int N)
{
	uvec nonzero_idx = find(ibi);

	if (nonzero_idx.n_elem > 8) {	//enough samples different from zero: linear interpolation
		int L = nonzero_idx.n_elem;
		int R = round(N / L);
		int ibi_interp_size;

		if (L*R >= N)
			ibi_interp_size = L*R;
		else
			ibi_interp_size = N;

		rowvec in_loc = linspace<rowvec>(0, L - 1, L);
		rowvec out_loc = linspace<rowvec>(0, L - 1, ibi_interp_size);
		rowvec ibi_interp = rowvec(ibi_interp_size);

		interp1(in_loc, ibi.elem(nonzero_idx), out_loc, ibi_interp);

		return ibi_interp;
	}
	else {	//not enough samples different from zero: fill zero values with the mean of nonzero values

			//no element different from zero, no interpolation
		if (nonzero_idx.n_elem <= 0)
			return ibi;

		uvec zero_idx = find(ibi == 0);

		//set the zero values to the mean
		ibi.elem(zero_idx).fill(mean(ibi.elem(nonzero_idx)));

		return ibi;
	}
}
