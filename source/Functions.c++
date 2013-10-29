#ifndef ConvectionAlgorithms_CPP
#define ConvectionAlgorithms_CPP

#include "Functions.h"

static const double EPSILON = 1E-7;

bool isLessThan(double first, double second)
{
	if(first - second < -EPSILON)
		return true;
	else
		return false;
}

bool isLessOrEqual(double first, double second)
{
	if(first - second < EPSILON)
		return true;
	else
		return false;
}

bool isEqual(double first, double second)
{
	if(fabs(first - second) < EPSILON)
		return true;
	else
		return false;
}

bool isGreaterThan(double first, double second)
{
	if(first - second > EPSILON)
		return true;
	else
		return false;
}

bool isGreaterOrEqual(double first, double second)
{
	if(first - second > -EPSILON)
		return true;
	else
		return false;
}

bool isOdd(int N)
{
	return (N % 2 != 0);
}

bool isEven(int N)
{
	return (N % 2 == 0);
}

#endif



