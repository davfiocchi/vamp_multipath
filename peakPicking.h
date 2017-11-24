//
//  peakPicking.h
//  
//
//  Created by Mattia Gabbrielli on 08/10/16.
//
//

#ifndef peakPicking_h
#define peakPicking_h

#include <stdio.h>
#include <armadillo>

using namespace arma;



rowvec peakPicking(rowvec ODF, int Fs, int frameRate, double threshold=0.25);

#endif /* peakPicking_h */
