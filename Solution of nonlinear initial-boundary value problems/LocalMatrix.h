//LocalMatrix.h

#pragma once

#include <vector>
#include "Mesh.h"
#include "BasisFunctions.h"

struct LocalMatrices
{
	std::vector<double> A;
	std::vector<double> b;
};

class LocalAssembler
{
private:
	Mesh *mesh_;
	BasisFunction *basis_;

public:
	LocalAssembler(Mesh &mesh, BasisFunction &basis)
	{
		mesh_ = &mesh;
		basis_ = &basis;
	}

	LocalMatrices simpleIteration(int elem_id, const std::vector<double> &q);
};
