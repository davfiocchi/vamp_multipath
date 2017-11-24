//
//  IBI.cpp
//  
//
//  Created by Davide Fiocchi on 29/10/16.
//
//

#include "IBI.h"

using namespace arma;

rowvec getIBI_vector(rowvec ODF) 
{
	printf("---GetIBI_vector begin!---\n");

	//not enough ODF samples to compute the IBI
	if (ODF.n_elem <= 0)
		return rowvec((uword)0, fill::zeros);

	//params definition
	int startPeriod = 20;
	int endPeriod = 128;
	int rgFrameLength = 512;
	int rgHoplength = 128;
	int beta = 43;

	//zeropad the odf
	std::vector<double> odf = conv_to<std::vector<double>>::from(ODF);
	std::cout << "size of ODF original: " << odf.size() << std::endl;
	odf.insert(odf.begin(), std::round(rgFrameLength / 2), 0);
	odf.insert(odf.end(), std::round(rgFrameLength / 2), 0);
	std::cout << "size of ODF padded: " << odf.size() << std::endl;

	//retrieve rythmogram dimension
	int frameIdx;
	for (frameIdx = 0; (frameIdx*rgHoplength) + rgFrameLength < odf.size(); frameIdx++);
	int nFrames = frameIdx;
	int nPeriods = endPeriod - startPeriod + 1;

	//create rythmogram
	mat rythmogram = mat(nPeriods, nFrames, fill::zeros);
	std::cout << "Rythmogram dimension: " << arma::size(rythmogram) << std::endl;

	//create the matrix out of the comb filter
	mat F = getCombFilterMatrix(startPeriod, endPeriod, rgFrameLength, beta);
	std::cout << "CombFilter matrix dimension: " << arma::size(F) << std::endl;

	//analyze each odf window and estract the rhythmogram
	printf("Analyze odf windows\n");
	for (frameIdx = 0; (frameIdx*rgHoplength) + rgFrameLength < odf.size(); frameIdx++)
	{
		std::vector<double> frame = trimVec(odf, frameIdx*rgHoplength, (frameIdx*rgHoplength) + rgFrameLength - 1);

		//difference frame - smooth
		frame = hwr(frame - smooth(frame, 15));

		colvec ac = colvec(unbiasedAutocorrelation(frame));

		rythmogram.col(frameIdx) = dot(ac, F);
	}

	printf("Rhythmogram found, start Path Finding\n");


	//---PATH FINDING---
	std::vector<int> periods = range(startPeriod, endPeriod + 1);
	std::vector<int> times = range(0, nFrames);

	//build the vector of points
	std::vector<Pt> pts;
	Pt p;
	for (int iC = 0; iC < nFrames; iC++)
	{
		for (int iR = 0; iR < nPeriods; iR++) {
			p.obs = rythmogram(iR, iC);
			p.period = periods[iR];
			p.time = times[iC];

			pts.push_back(p);
		}
	}

	int nPts = pts.size();

	std::vector<int> backlink = std::vector<int>(nPts, -1);
	std::vector<double> cumscore = std::vector<double>(nPts, 0);
	std::vector<double> score;

	float sigma = 0.04;
	int threshold = 1;

	printf("Vector of points built, iterate over all points to build backlink and cumscore\n");
	for (int i = 0; i < nPts; i++)
	{
		p = pts[i];
		Tsi iCtS = transitionScore(p, pts, sigma, threshold);

		if (iCtS.indeces.size() != iCtS.tS.size())
			printf("We've got a situation here..two array dimension doesn't match\n");

		if (std::any_of(iCtS.indeces.begin(), iCtS.indeces.end(), [](int i) {return i < 0; })) {
			score = std::vector<double>(iCtS.tS);

			for (int i = 0; i < iCtS.indeces.size(); i++)
				if (iCtS.indeces[i] >= 0)
					score[i] += cumscore[i];
		}
		else {
			score = getValuesAtIndeces(cumscore,iCtS.indeces) + iCtS.tS;
		}

		std::vector<double>::iterator it = std::max_element(score.begin(), score.end());
		int ii = it - score.begin();

		backlink[i] = iCtS.indeces[ii];
		cumscore[i] = score[ii] + p.obs;
	}

	printf("Points analyzed, find termination\n");
	int terminationIndex = findTermination(backlink, cumscore);

	printf("Termination found! termination index: %d\n", terminationIndex);

	std::vector<int> path = std::vector<int>(1, terminationIndex);
	int indx = *(path.end() - 1);	//path's last value
	while (backlink[indx]>=0)
	{
		path.push_back(backlink[indx]);
		indx = backlink[indx];
	}

	std::vector<int> ibi_noInterp = std::vector<int>(path.size());
	for (int i = 0; i < path.size(); i++) {
		if (i != (path.size() - 1))
			ibi_noInterp[i] = path[i] - path[i + 1];
		else //last element
			ibi_noInterp[i] = path[i];
	}

	rowvec pathReverse = rowvec(path.size());
	for (int i = 0; i < path.size(); i++)
		pathReverse(i) = *(ibi_noInterp.end() - 1 - i);

	printVec("Ibi not interpolated: ", ibi_noInterp);

	rowvec ibi = interpolatePeriods(pathReverse, rgHoplength, ODF.n_elem);

	std::cout << "Ibi size: " << arma::size(ibi) << std::endl;
	
	return ibi;
}

mat getCombFilterMatrix(int startPeriod, int endPeriod, int frameLength, int rayParam)
{
	printf("Start getComFilterMatrix\n");
	int nImpulses = std::floor(frameLength / endPeriod);
	int nPeriods = endPeriod - startPeriod + 1;

	std::vector<double> wG = getWeightCurve(nPeriods, rayParam);

	mat F = mat(frameLength, nPeriods, fill::zeros);

	for (int period = startPeriod; period <= endPeriod; period++) {
		int col = period - startPeriod;

		for (int impulse = 1; impulse <= nImpulses; impulse++) {

			int w = 2 * impulse - 1;
			int hw = std::floor(w / 2);

			std::vector<int> rows = range(impulse*period-hw, 
				(int)fmin(impulse*period+hw+1,frameLength-1));

			for (auto itR = rows.begin(); itR < rows.end(); itR++) {
				F(*itR, col) = (1 / w) * wG[col];
			}
		}
	}

	printf("Matrix out of comb filter obtained!\n");

	return F;
}

std::vector<double> getWeightCurve(int dimension, double rayParam)
{
	printf("GetWeightCurve begin...");
	std::vector<double> w;

	//rayleigh weighting curve
	for (int n = 1; n <= dimension; n++)
		w.push_back((n / std::pow(rayParam, 2))*std::exp((-1 * std::pow(-n, 2))/ (2 * std::pow(rayParam, 2))));

	printf("weighting vector generated!\n");
	return w;
}

std::vector<double> smooth(std::vector<double> vec, int winDim)
{
	if (vec.size() < winDim) {
		printf("IBI.smooth error: vector dimension (%d) < window dimension (%d)\n", vec.size(), winDim);
		return vec;
	}
	else if (winDim < 3) {
		printf("IBI.smooth error: window dimension (%d) too small\n", winDim);
		return vec;
	}

	int beginIdx = 0;
	int endIdx = vec.end() - vec.begin() - 1;
	
	//mirror padding
	for (int i = 0; i < std::round(winDim / 2); i++) {
		vec.insert(vec.begin(), 1, vec[++beginIdx]);
		vec.insert(vec.end(), 1, vec[endIdx]);
	}

	//rectangular window normalized
	std::vector<double> win = std::vector<double>(winDim, 1 / double(winDim));
	
	//convolution
	std::vector<double> cnv = conv_to<std::vector<double>>::from(conv(rowvec(vec), rowvec(win), "same"));
	
	return std::vector<double>(cnv.begin()+std::round(winDim/2), cnv.end() - std::round(winDim / 2));
}

std::vector<double> hwr(std::vector<double> vec)
{
	for (int i = 0; i < vec.size(); i++)
		if (vec[i] < 0)
			vec[i] = 0;

	return vec;
}

std::vector<double> unbiasedAutocorrelation(std::vector<double> vec)
{
	int N = vec.size();	//equals to the framelength
	std::vector<int> lags = range(-N+1, N);

	//compute autocorrelation
	std::vector<double> ac;
	double sum;
	double max = 0;

	for (auto it = lags.begin(); it < lags.end(); it++) 
	{
		int lag = *it;
		sum = 0;
		
		for (int n = 0; n < N; n++) 
			if ((n + lag >= 0) && (n + lag < N)) 
				sum += vec[n + lag] * vec[n];

		ac.push_back(sum/(N-std::abs(lag)));

		//update maximum value
		if (sum > max)
			max = sum;
	}

	//normalize
	if (max != 0) {
		for (int i = 0; i < lags.size(); i++)
			ac[i] = ac[i] / max;
	}

	//return the positive lags
	if (N % 2 == 0)
		return trimVec(ac, std::floor(lags.size()/2), lags.size() - 1);
	else
		return trimVec(ac, std::ceil(lags.size()/2), lags.size() - 1);
}

colvec dot(colvec vec, mat m)
{
	if (vec.n_elem != m.n_rows) {
		printf("IBI.dot error: vector size (%d) and n_rows (%d) of matrix doesn't coincide\n", vec.size(), m.n_rows);
		return colvec(m.n_cols, fill::zeros);
	}

	std::vector<double> result;

	for (int i = 0; i < m.n_cols; i++) {
		result.push_back(arma::dot(vec,m.col(i)));
	}

	return colvec(result);
}

Tsi transitionScore(Pt p, std::vector<Pt> pts, float sigma, int threshold)
{
	Tsi retVal;

	std::vector<int> t_pts;
	for (int i = 0; i < pts.size(); i++)
		t_pts.push_back(pts[i].time);

	//get candidates
	std::vector<int> iC = searchSortedRange(t_pts, p.time, threshold);

	if (iC.size() == 0) {	//no candidates
		retVal.indeces = std::vector<int>(1, -1);
		retVal.tS = std::vector<double>(1, 1);

		return retVal;
	}

	retVal.indeces = std::vector<int>(iC);

	//calculate the transition score
	std::vector<Pt> p_pts = getValuesAtIndeces(pts, iC);

	double delta;
	for (int i = 0; i < p_pts.size(); i++) {
		delta = std::abs(std::log(p.period / p_pts[i].period));
		retVal.tS.push_back(std::exp(-1 * std::pow(delta, 2) / (2 * std::pow(sigma, 2))));
	}

	return retVal;
}

std::vector<int> searchSortedRange(std::vector<int> vec, int value, int threshold)
{
	std::vector<int>::iterator it1 = std::search_n(vec.begin(), vec.end(), 1, value-threshold);
	std::vector<int>::iterator it2 = std::search_n(vec.begin(), vec.end(), 1, value);

	int i1 = it1 - vec.begin();
	int i2 = it2 - vec.begin();
	
	return range(i1, i2);
}

int findTermination(std::vector<int> backlink, std::vector<double> cumscore)
{
	int beginIndxBackLink = *(backlink.end() - 1);
	if (beginIndxBackLink > (int)backlink.size()) {
		printf("IBI.findTermination error: the last element of backlink (%d) is bigger than the size(%d)\n", beginIndxBackLink, backlink.size());
		return 1;
	}

	int maxIndxCumScore;	//index maximum value in cum score
	std::vector<double>::iterator cmIt;
	if (beginIndxBackLink >= 0) {
		cmIt = std::max_element(cumscore.begin() + beginIndxBackLink, cumscore.begin() + backlink.size());
		maxIndxCumScore = cmIt - (cumscore.begin() + beginIndxBackLink);
	}
	else {
		cmIt = std::max_element(cumscore.begin(), cumscore.end());
		maxIndxCumScore = cmIt - cumscore.begin();
	}
	
	return beginIndxBackLink + maxIndxCumScore;
}

rowvec interpolatePeriods(rowvec path, int periodsHop, int length)
{
	rowvec periodsX = rowvec(path.size());
	rowvec periodsXnew = rowvec(length);
	rowvec ibi;

	for (int i = 0; i < path.size(); i++)
		periodsX(i) = i*periodsHop;

	for (int i = 0; i < length; i++)
		periodsXnew(i) = i;

	try{ interp1(periodsX, path, periodsXnew, ibi, "*linear");}
	catch (const std::exception&) 
	{ //in case there are no two values different from each other, fill the vector with the same value
		ibi = rowvec(length);
		ibi.fill(path(0));
	}

	for (int i = 0; i < length; i++) {
		if (periodsXnew(i) < periodsX(0))
			ibi(i) = path(0);
		if (periodsXnew(i) > periodsX(path.n_elem-1))
			ibi(i) = path(path.n_elem-1);
	}

	return ibi;
}

std::vector<int> range(int beginning, int end)
{
	if (end <= beginning) 
	{
		return std::vector<int>(0);
	}

	std::vector<int> vec;

	for (int i = beginning; i < end; i++)
		vec.push_back(i);

	return vec;
}