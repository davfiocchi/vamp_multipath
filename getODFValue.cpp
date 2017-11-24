//
//  getODFValue.cpp
//  
//
//  Created by Mattia Gabbrielli on 08/10/16.
//
//

#include "getODFValue.h"

using namespace std;
using namespace arma;

double getODFValue(cx_rowvec frameSpectrum, cx_rowvec pastFrameSpectrum)
{
	double envValue = 0.0;
	rowvec temp; 				// rowvec is a Row<double>

	if (pastFrameSpectrum.is_empty())		//if pastFrameSpectrum is defined -> spectralFlux
	{
		temp = norm(log10(frameSpectrum), 1) - norm(log10(pastFrameSpectrum), 1);
		temp = (temp + abs(temp)) / 2;
	}
	else 						//pastFrameSpectrum is NULL -> power spectrum envelope
	{
		//temp = abs(frameSpectrum) % abs(frameSpectrum);
		temp = pow(abs(frameSpectrum), 2);
	}

	envValue = sum(temp);

	return envValue;
}

rowvec normalizeODF(rowvec df, double maxIbi)
{
	rowvec normOdf = df / stddev(rowvec(df), 1);

	//---SMOOTH---
	int winDim = 2 * std::round(maxIbi)+1;
	//winDim = 11;
	if (normOdf.size() < winDim) {
		printf("ODF.normalize error: vector dimension (%d) < window dimension (%d)\n", normOdf.size(), winDim);
		return normOdf;
	}
	else if (winDim < 3) {
		printf("ODF.normalize error: window dimension (%d) too small\n", winDim);
		return normOdf;
	}

	std::vector<double> vNormOdf = conv_to<std::vector<double>>::from(normOdf);

	int beginIdx = 0;
	int endIdx = vNormOdf.end() - vNormOdf.begin() - 1;

	//mirror padding
	for (int i = 0; i < std::round(winDim / 2); i++) {
		vNormOdf.insert(vNormOdf.begin(), 1, vNormOdf[++beginIdx]);
		vNormOdf.insert(vNormOdf.end(), 1, vNormOdf[endIdx]);
	}

	//gaussian window
	rowvec win = rowvec(winDim);
	double maxValue = (1 / (2 * datum::pi))*exp(0);	//at index 0 we have the maximum value
	double sum = 0;
	for (int i = 0, n = -std::floor(winDim/2); i < winDim; i++, n++) {
		win(i) = (1 / (2 * datum::pi))*exp(-1 * pow(n, 2) / 2) / maxValue;
		sum += win(i);
	}

	//normalize window
	win = win / sum;
	
	//convolution
	//normOdf = conv(rowvec(vNormOdf), win, "same");
	std::vector<double> cnv = conv_to<std::vector<double>>::from(conv(rowvec(vNormOdf), win, "same"));
	cnv = std::vector<double>(cnv.begin() + std::round(winDim / 2), cnv.end() - std::round(winDim / 2));
	//printVec("Convolution: ", cnv);
	
	return rowvec(cnv);
}
