#pragma once


class RandomNumberGenerator {

private:
	unsigned long maxValueM;
	unsigned long multiplierA;
	unsigned long incrementC;
public:
	unsigned long seed;

	RandomNumberGenerator(unsigned long seed);
	double GetRandom();
	double GetRandom(double lo, double hi);


};