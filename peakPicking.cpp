//
//  peakPicking.cpp
//  
//
//  Created by Mattia Gabbrielli on 08/10/16.
//
//

#include "peakPicking.h"
#include <armadillo>

using namespace std;
using namespace arma;

rowvec peakPicking(rowvec ODF, int Fs, int frameRate, double threshold)
{	
	/* Parameters setup */

	int hopsize = round(Fs/frameRate);
	int N = ODF.n_elem;
	double alpha = 0.5;

	double pre_max = 0.030;
	double post_max = 0.030;
	double pre_avg = 0.100;
	double post_avg = 0.070;

	int pre_max_int = round((pre_max*Fs)/hopsize);
	int post_max_int = round((post_max*Fs)/hopsize);
	int pre_avg_int = round((pre_avg*Fs)/hopsize);
	int post_avg_int = round((post_avg*Fs)/hopsize);

	int lower = max(pre_max_int,pre_avg_int);
	int upper = N - max(post_max_int,post_avg_int);
	
	rowvec g = zeros(N);
	rowvec peaks = zeros(N);
	
	/* ODF preprocessing */

	double mu = mean(ODF);
	double sd = stddev(ODF);
	ODF = (ODF - rowvec(N,mu))/sd;

	/* Peak detection algorithm */
	for (uword i = 1; i < N; ++i)
	{
		g(i) = max(ODF(i), alpha*g(i-1) + (1-alpha)*ODF(i));
		if (i > lower && i < upper)
		{	
			//ODF must be the local maxima in the interval [i - pre_max_int, i + post_max_int]
			if (ODF(i) == max(ODF.subvec(i - pre_max_int, i + post_max_int)))	
			{
				// ODF must be greater than the mean in the interval [i - pre_avg_int, i + post_avg_int]
				if (ODF(i) >= mean(ODF.subvec(i - pre_avg_int, i + post_avg_int)) + threshold)
				{
					// ODF must be thresholded
					if (ODF(i) >= g(i-1))
					{
						peaks(i) = 1.0;
					}
				}
			}
		}
	}

	return peaks;
}
