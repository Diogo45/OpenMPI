

#include "pch.h"
#include <cmath>

RandomNumberGenerator::RandomNumberGenerator(unsigned long seed)
{
	this->seed = seed;
	this->maxValueM = (long) pow(2, 32);
	this->multiplierA = 1664525;
	this->incrementC = 1013904223;
}

double RandomNumberGenerator::GetRandom()
{
	seed = (multiplierA * seed + incrementC) % maxValueM;
	return (double) seed / maxValueM;
}

double RandomNumberGenerator::GetRandom(double lo, double hi)
{
	seed = (multiplierA * seed + incrementC) % maxValueM;
	return ((double)seed / maxValueM) * (hi - lo) + lo;
}