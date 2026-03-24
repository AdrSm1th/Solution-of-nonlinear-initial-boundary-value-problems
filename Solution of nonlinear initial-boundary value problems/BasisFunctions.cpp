//BasisFunctions.cpp

#include "BasisFunctions.h"

std::vector<double> BasisFunction::evaluateBasis(double xi) const
{
	double N1 = 0.5 * xi * (xi - 1);
	double N2 = 1 - xi * xi;
	double N3 = 0.5 * xi * (xi + 1);
	return { N1, N2, N3 };
}

std::vector<double> BasisFunction::evalueteDerivatives(double xi) const
{
	double N1 = xi - 0.5;
	double N2 = -2 * xi;
	double N3 = xi + 0.5;
	return { N1, N2, N3 };
}