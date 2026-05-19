//LocalMatrix.cpp

#include "LocalMatrix.h"

LocalMatrices LocalAssembler::simpleIteration(int elem_id, const std::vector<double> &q)
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

	w[0] = 0.555555555555556;
	w[1] = 0.888888888888889;
	w[2] = 0.555555555555556;
	double J = 0;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				double u_h = 0;
				double xi = (2 * x[k] - (x1 + x3)) / (x3 - x1);
				std::vector<double> psi = basis_->evaluateBasis(xi);
				std::vector<double> dpsi_dx = basis_->evalueteDerivatives(xi);

				for (int jg = 0; jg < 3; jg++)
				{
					u_h += q[elem_nodes[jg]] * psi[jg];
				}
				//G component
				matrices.A[i * 3 + j] += mesh_->lambda(u_h) * dpsi_dx[i] * dpsi_dx[j] * w[k] * (2 / h);
				//M component
				matrices.A[i * 3 + j] += mesh_->gamma() * psi[i] * psi[j] * w[k] * (h / 2);
			}
		}

		//b component
		for (int iF = 0; iF < 3; iF++)
		{
			double xi = (2 * x[iF] - (x1 + x3)) / (x3 - x1);
			std::vector<double> psi = basis_->evaluateBasis(xi);
			std::vector<double> dpsi_dx = basis_->evalueteDerivatives(xi);
			matrices.b[i] += mesh_->f(x[iF]) * psi[i] * w[iF];
		}
		matrices.b[i] *= (h / 2);
	}

	return matrices;
}

LocalMatrices LocalAssembler::newtonIteration(int elem_id, const std::vector<double> &q)
{
   LocalMatrices matrices;
   matrices.A.resize(9, 0.0);
   matrices.b.resize(3, 0.0);

   auto elem_nodes = mesh_->getElementNodes(elem_id);
   double x1 = mesh_->getNodeCoord(elem_nodes[0]);
   double x2 = mesh_->getNodeCoord(elem_nodes[1]);
   double x3 = mesh_->getNodeCoord(elem_nodes[2]);
   double h = x3 - x1;

   const double gauss_pts[3] = { -0.774596669241483, 0.0, 0.774596669241483 };
   const double gauss_w[3] = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };

   for (int k = 0; k < 3; ++k)
   {
      double xi = gauss_pts[k];
      double x = 0.5 * (x1 + x3) + 0.5 * h * xi;

      std::vector<double> psi = basis_->evaluateBasis(xi);
      std::vector<double> dpsi_dxi = basis_->evalueteDerivatives(xi);

      std::vector<double> dpsi_dx(3);
      for (int i = 0; i < 3; ++i)
         dpsi_dx[i] = dpsi_dxi[i] * 2.0 / h;

      double u_h = 0.0, du_h = 0.0;
      for (int j = 0; j < 3; ++j)
      {
         u_h += q[elem_nodes[j]] * psi[j];
         du_h += q[elem_nodes[j]] * dpsi_dx[j];
      }

      double lam = mesh_->lambda(u_h);
      double lam_der = mesh_->lambda_derivative(u_h);
      double gamma = mesh_->gamma();
      double f_val = mesh_->f(x);

      double weight = gauss_w[k] * (h / 2.0);

      for (int i = 0; i < 3; ++i)
      {
         matrices.b[i] += (f_val * psi[i] + lam_der * du_h * u_h * dpsi_dx[i]) * weight;

         for (int j = 0; j < 3; ++j)
         {

            double A_ij = lam * dpsi_dx[i] * dpsi_dx[j] + gamma * psi[i] * psi[j];

            double S_ij = lam_der * du_h * dpsi_dx[i] * psi[j];
            matrices.A[i * 3 + j] += (A_ij + S_ij) * weight;
         }
      }
   }

   return matrices;
}