//BasisFunctions.cpp

#include "BasisFunctions.h"

std::vector<double> BasisFunction::evaluateBasis(double xi) const
{
	double N1 = 0.5 * xi * (xi - 1);
	double N2 = 0.5 * xi * (xi + 1);
	double N3 = 1 - xi * xi;
	return { N1, N2, N3 };
}

std::vector<double> BasisFunction::evalueteDerivatives(double xi) const
{
	double N1 = xi - 0.5;
	double N2 = xi + 0.5;
	double N3 = -2 * xi;
	return { N1, N2, N3 };
}