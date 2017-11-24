#ifndef V_UTILS
#define V_UTILS

#include <stdio.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
#include <cassert>

using namespace std;

/*
* Returns a slice of the input vector, from startIdx position to endIdx position (both included)
* @vec: input vector
* @startIdx: starting index
* @endIdx: ending index
*/
template <typename T>
vector<T> trimVec(vector<T> vec, int startIdx, int endIdx)
{
	if ((endIdx < startIdx) || (startIdx < 0) || (endIdx >= vec.size()))
		return vec;

	vector<T>::iterator begin = vec.begin();
	begin = begin + startIdx;

	vector<T>::iterator end = vec.begin();
	end = end + endIdx + 1;

	return vector<T>(begin, end);
}

/*
* Returns a vector vec[indeces]
* @vec: input vector
* @indeces: vector with ordered indeces
*/
template <typename T>
vector<T> getValuesAtIndeces(vector<T> vec, vector<int> indeces)
{
	if (vec.size() == 0 || indeces.size() == 0)
		return vector<T>(0);

	vector<T> out;

	for (int i = 0, j = 0; i < vec.size(); i++)
	{
		if (i == indeces[j]) {	//index present
			out.push_back(vec[i]);
			if (++j == indeces.size()) //no more index
				return out;
		}
	}

	return out;
}

/*
* Prints the values of a vector
* @vec: the vector to be printed
*/
template <typename T>
void printVec(char* header, vector<T> vec)
{
	cout << header << "[";
	for (vector<T>::iterator i = vec.begin(); i != vec.end(); i++)
		cout << ' ' << *i;

	cout << "] " << endl;
}


/*
* Implements the sum of two vectors (element by element)
* @a: first vector
* @b: second vector
*/
template <typename T>
vector<T> operator+(const vector<T>& a, const vector<T>& b)
{
	assert(a.size() == b.size());

	vector<T> result;
	result.reserve(a.size());

	transform(a.begin(), a.end(), b.begin(),
		back_inserter(result), plus<T>());
	return result;
}

/*
* Implements the difference of two vectors (element by element)
* @a: first vector
* @b: second vector
*/
template <typename T>
vector<T> operator-(const vector<T>& a, const vector<T>& b)
{
	assert(a.size() == b.size());

	vector<T> result;
	result.reserve(a.size());

	transform(a.begin(), a.end(), b.begin(),
		back_inserter(result), minus<T>());
	return result;
}

/*
* Implements the multiplication of two vectors (element by element)
* @a: first vector
* @b: second vector
*/
template <typename T>
vector<T> operator*(const vector<T>& a, const vector<T>& b)
{
	assert(a.size() == b.size());

	vector<T> result;
	result.reserve(a.size());

	transform(a.begin(), a.end(), b.begin(),
		back_inserter(result), multiplies<T>());
	return result;
}

#endif /* V_UTILS */
#pragma once