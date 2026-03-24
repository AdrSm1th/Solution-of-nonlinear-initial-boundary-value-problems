//BasisFunctions.cpp

#include "BasisFunctions.h"

//std::vector<double> BasisFunction::evaluateBasis(double x, double x1, double x2, double x3) const
//{
//	double N1 = (2 * (x - x2) * (x - x3)) / ((x1 - x2) * (x1 - x3));
//	double N2 = (2 * (x - x1) * (x - x3)) / ((x2 - x1) * (x2 - x3));
//	double N3 = (2 * (x - x1) * (x - x2)) / ((x3 - x1) * (x3 - x2));
//	return { N1, N2, N3 };
//}
//
//std::vector<double> BasisFunction::evalueteDerivatives(double x, double x1, double x2, double x3) const
//{
//	//double N1 = (2 * x - x2 - x3) / ((x1 - x3) * (x1 - x2));
//	//double N2 = (2 * x - x1 - x3) / ((x2 - x1) * (x2 - x3));
//	//double N3 = (2 * x - x1 - x2) / ((x3 - x1) * (x3 - x2));
//
//	//double N1 = (4 * x - 2 * x3 - 2 * x2) / ((x2 - x1) * x3 - x1 * x2 + x1 * x1);
//	//double N2 = (-4 * x + 2 * x3 + 2 * x1) / ((x2 - x1) * x3 - x2 * x2 + x1 * x2);
//	//double N3 = (4 * x - 2 * x3 - 2 * x2) / (x3 * x3 + (-x1 - x2) * x3 + x1 * x2);
//
//
//	return { N1, N2, N3 };
//}

std::vector<double> BasisFunction::evaluateBasis(double xi) const
{
	double N1 = 0.5 * xi * (xi - 1);
	double N2 = 1 - xi * xi;
	double N3 = 0.5 * xi * (xi + 1);
	return { N1, N2, N3 };
}

std::vector<double> BasisFunction::evalueteDerivatives(double xi) const
{
	//double N1 = (2 * x - x2 - x3) / ((x1 - x3) * (x1 - x2));
	//double N2 = (2 * x - x1 - x3) / ((x2 - x1) * (x2 - x3));
	//double N3 = (2 * x - x1 - x2) / ((x3 - x1) * (x3 - x2));

	//double N1 = (4 * x - 2 * x3 - 2 * x2) / ((x2 - x1) * x3 - x1 * x2 + x1 * x1);
	//double N2 = (-4 * x + 2 * x3 + 2 * x1) / ((x2 - x1) * x3 - x2 * x2 + x1 * x2);
	//double N3 = (4 * x - 2 * x3 - 2 * x2) / (x3 * x3 + (-x1 - x2) * x3 + x1 * x2);

	double N1 = xi - 0.5;
	double N2 = -2 * xi;
	double N3 = xi + 0.5;
	return { N1, N2, N3 };
}