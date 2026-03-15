//BasisFunctions.h

#pragma once

#include <vector>

class BasisFunction
{
public:
	std::vector<double> evaluateBasis(double x) const;
	std::vector<double> evalueteDerivatises(double x) const;
};