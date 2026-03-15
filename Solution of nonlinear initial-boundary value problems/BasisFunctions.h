//BasisFunctions.h

#pragma once

#include <vector>

class BasisFunction
{
public:
	std::vector<double> evaluateBasis(double xi) const;
	std::vector<double> evalueteDerivatives(double xi) const;
};