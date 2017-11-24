//
//  getODFValue.h
//  
//
//  Created by Mattia Gabbrielli on 08/10/16.
//
//

#ifndef getODFValue_h
#define getODFValue_h

#include <stdio.h>
#include <armadillo>
#include <math.h>


using namespace arma;

/*
*	getODFValue returns a value of the OnSetDetection function for the frame we pass.
*  There are two algorithm implemented: 1) Power Spectrum 2) Simple Flux
*	To use Power Spectrum, just pass an empty vector in pastframeSpectrum parameter e.g. getODFValue(frame, cx_rowvec())
*	To use simple flux, pass both the parameters.
*  Note that for simple flux is better if the pastFrameSpectrum is 2 frame before the current frame.
*/

double getODFValue(cx_rowvec frameSpectrum, cx_rowvec pastFrameSpectrum);

rowvec normalizeODF(rowvec df, double maxIbi);

#endif /* getODFValue_h */