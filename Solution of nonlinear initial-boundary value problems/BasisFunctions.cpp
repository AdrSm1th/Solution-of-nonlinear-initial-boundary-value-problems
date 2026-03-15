//BasisFunctions.cpp

#include "BasisFunctions.h"

std::vector<double> BasisFunction::evaluateBasis(double x, double x1, double x2, double x3) const
{
	double N1 = ((x - x3) * (x - x2)) / ((x1 - x3) * (x1 - x2));
	double N2 = ((x - x1) * (x - x3)) / ((x2 - x1) * (x2 - x3));
	double N3 = ((x - x1) * (x - x2)) / ((x3 - x1) * (x3 - x2));
	return { N1, N2, N3 };
}

std::vector<double> BasisFunction::evalueteDerivatives(double x, double x1, double x2, double x3) const
{
	double N1 = (2 * x - x2 - x3) / ((x1 - x3) * (x1 - x2));
	double N2 = (2 * x - x1 - x3) / ((x2 - x1) * (x2 - x3));
	double N3 = (2 * x - x1 - x2) / ((x3 - x1) * (x3 - x2));
	return { N1, N2, N3 };
}