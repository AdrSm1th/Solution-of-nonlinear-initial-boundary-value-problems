//LocalMatrix.cpp

#include "LocalMatrix.h"

LocalMatrices LocalAssembler::simpleIteration(int elem_id, const std::vector<double> &q_local)
{
	LocalMatrices matrices;
	matrices.A.resize(9);
	matrices.b.resize(3);

	auto elem_nodes = mesh_->getElementNodes(elem_id);

	double x1 = mesh_->getNodeCoord(elem_nodes[0]);
	double x2 = mesh_->getNodeCoord(elem_nodes[1]);
	double x3 = mesh_->getNodeCoord(elem_nodes[2]);
	double h = x3 - x1;

	std::vector<double> x(3), w(3);

	const std::vector<double> gauss_w = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };
	
	x[0] = (x1 + x3) / 2 - (x3 - x1) / 2 * 0.774596669241483;
	x[1] = (x1 + x3) / 2;
	x[2] = (x1 + x3) / 2 + (x3 - x1) / 2 * 0.774596669241483;

	//w[0] = (x3 - x1) / 2 * 0.555555555555556;
	//w[1] = (x3 - x1) / 2 * 0.888888888888889;
	//w[2] = (x3 - x1) / 2 * 0.555555555555556;

	w[0] = 0.555555555555556;
	w[1] = 0.888888888888889;
	w[2] = 0.555555555555556;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			//G component
			for (int ig = 0; ig < 3; ig++)
			{
				double u_h = 0;
				double xi = (2 * x[ig] - (x1 + x3)) / (x3 - x1);
				std::vector<double> psi = basis_->evaluateBasis(xi);
				std::vector<double> dpsi_dx = basis_->evalueteDerivatives(xi);

				for (int k = 0; k < 3; k++)
				{
					dpsi_dx[k] *= 2 / h;
				}

				for (int jg = 0; jg < 3; jg++)
				{
					u_h += q_local[elem_nodes[jg]] * psi[jg];
				}
				matrices.A[i * 3 + j] += mesh_->lambda(u_h) * dpsi_dx[i] * dpsi_dx[j] * w[ig] * (h / 2);
			}

			//M component
			for (int im = 0; im < 3; im++)
			{
				double xi = (2 * x[im] - (x1 + x3)) / (x3 - x1);
				std::vector<double> psi = basis_->evaluateBasis(xi);
				matrices.A[i * 3 + j] += mesh_->gamma() * psi[i] * psi[j] * w[im] * (h / 2);
			}
		}

		for (int iF = 0; iF < 3; iF++)
		{
			double xi = (2 * x[iF] - (x1 + x3)) / (x3 - x1);
			std::vector<double> psi = basis_->evaluateBasis(xi);
			matrices.b[i] += mesh_->f(x[iF]) * psi[i] * w[iF] * (h / 2);
		}
	}

	return matrices;
}