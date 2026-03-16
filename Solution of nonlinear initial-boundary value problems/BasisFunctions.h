//BasisFunctions.h

#pragma once

#include <vector>

class BasisFunction
{
public:
	//std::vector<double> evaluateBasis(double x, double x1, double x2, double x3) const;
	//std::vector<double> evalueteDerivatives(double x, double x1, double x2, double x3) const;
	std::vector<double> evaluateBasis(double xi) const;
	std::vector<double> evalueteDerivatives(double xi) const;
};