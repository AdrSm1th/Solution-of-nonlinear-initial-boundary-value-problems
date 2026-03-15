//LocalMatrix.cpp

#include "LocalMatrix.h"

LocalMatrices LocalAssembler::simpleIteration(int elem_id, const std::vector<double> &q_local)
{
	LocalMatrices matrices;
	matrices.A.resize(9);
	matrices.b.resize(3);

	auto elem_nodes = mesh_->getElementNodes(elem_id);

	double x1 = mesh_->getNodeCoord(elem_nodes[0]);
	double x2 = mesh_->getNodeCoord(elem_nodes[2]);
	double x3 = mesh_->getNodeCoord(elem_nodes[1]);

	std::vector<double> x(3), w(3);

	const std::vector<double> gauss_w = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };
	
	x[0] = (x1 + x2) / 2 - (x2 - x1) / 2 * 0.774596669241483;
	x[1] = (x1 + x2) / 2;
	x[2] = (x1 + x2) / 2 + (x2 - x1) / 2 * 0.774596669241483;

	w[0] = (x2 - x1) / 2 * 0.555555555555556;
	w[1] = (x2 - x1) / 2 * 0.888888888888889;
	w[2] = (x2 - x1) / 2 * 0.555555555555556;
	
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			//G component
			for (int ig = 0; ig < 3; ig++)
			{
				double u_h = 0;
				std::vector<double> psi = basis_->evaluateBasis(x[ig], x1, x2, x3);
				std::vector<double> dpsi_dx = basis_->evalueteDerivatives(x[ig], x1, x2, x3);
				for (int jg = 0; jg < 3; jg++)
				{
					u_h += q_local[jg] * psi[jg];
				}
				matrices.A[i * 3 + j] += mesh_->lambda(u_h) * dpsi_dx[i] * dpsi_dx[j] * w[ig];
			}

			//M component
			for (int im = 0; im < 3; im++)
			{
				std::vector<double> psi = basis_->evaluateBasis(x[im], x1, x2, x3);
				matrices.A[i * 3 + j] += mesh_->gamma() * psi[i] * psi[j] * w[im];
			}
		}

		for (int iF = 0; iF < 3; iF++)
		{
			std::vector<double> psi = basis_->evaluateBasis(x[iF], x1, x2, x3);
			matrices.b[i] += mesh_->f(x[iF]) * psi[i] * w[iF];
		}
	}

	return matrices;
}