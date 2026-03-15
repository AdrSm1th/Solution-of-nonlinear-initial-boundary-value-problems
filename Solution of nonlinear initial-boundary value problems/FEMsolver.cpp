//FEMsolver.cpp

#include "FEMsolver.h"

int FEMsolver::index(int i, int j) const
{
	if (std::abs(i - j) > matrix_bandwidth_) return -1;

	int diag_offset = j - i;

	int diag_num = diag_offset + matrix_bandwidth_;

	return diag_num * mesh_->getNumNodes() + i;
}

void FEMsolver::assembleGlobalSystem(const std::vector<double> &q_current)
{
	std::fill(global_A_.begin(), global_A_.end(), 0.0);
	std::fill(global_b_.begin(), global_b_.end(), 0.0);

	for (int elem = 0; elem < mesh_->getNumElems(); elem++)
	{
		LocalMatrices matrices = assembler_->simpleIteration(elem, q_current);
		std::vector<int> elem_nodes = mesh_->getElementNodes(elem);
		
		for (int i = 0; i < 3; i++)
		{
			int i_global = elem_nodes[i];
			for (int j = 0; j < 3; j++)
			{
				int j_global = elem_nodes[j];
				(*this)(i_global, j_global) += matrices.A[i * 3 + j];
			}
			global_b_[i_global] += matrices.b[i];
		}
	}
}

void FEMsolver::applyBoundaryCondition()
{
	BoundaryCondition left = mesh_->getBoundaryCondition(true);
	BoundaryCondition right = mesh_->getBoundaryCondition(false);
	int n = mesh_->getNumNodes();
	for (int cond_id = 0; cond_id < 2; cond_id++)
	{
		bool left = cond_id == 0;
		BoundaryCondition cond = mesh_->getBoundaryCondition(left);
		int idx = left ? 0 : n - 1;

		switch (cond.type)
		{
			case 1:
			{
				for (int j = 0; j < n; j++) 
				{
					if (std::abs(j - idx) <= matrix_bandwidth_) (*this)(idx, j) = 0;
				}
				(*this)(idx, idx) = 1;
				global_b_[idx] = cond.u_g;

				break;
			}
			case 2:
			{
				auto elem_nodes = mesh_->getElementNodes(idx);

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

				int local_idx = (idx == elem_nodes[0]) ? 0 :
					(idx == elem_nodes[1]) ? 1 : 2;
				for (int i = 0; i < 3; i++)
				{
					std::vector<double> psi = basis_->evaluateBasis(x[i], x1, x2, x3);
					global_b_[idx] += cond.theta * psi[local_idx] * w[i];
				}

				break;
			}

			case 3:
			{
				auto elem_nodes = mesh_->getElementNodes(idx);

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

				int local_idx = (idx == elem_nodes[0]) ? 0 :
					(idx == elem_nodes[1]) ? 1 : 2;
				for (int i = 0; i < 3; i++)
				{
					std::vector<double> psi = basis_->evaluateBasis(x[i], x1, x2, x3);
					(*this) (idx, idx) += cond.beta * psi[i] * psi[local_idx] * w[i];
					global_b_[idx] += cond.u_beta * cond.beta * psi[local_idx] * w[i];
				}

				break;
			}
		}
	}
}

bool FEMsolver::solveLU()
{
	FEMsolver Acopy = (*this);

	int n = mesh_->getNumNodes();
	for (int k = 0; k < n; ++k) 
	{

		for (int j = k; j <= std::min(n - 1, k + matrix_bandwidth_); ++j) 
		{
			double sum = 0.0;
			for (int t = std::max(0, k - matrix_bandwidth_); t < k; ++t)
			{
				if (std::abs(t - k) <= matrix_bandwidth_ && std::abs(t - j) <= matrix_bandwidth_) 
				{
					sum += Acopy(k, t) * Acopy(t, j);
				}
			}
			Acopy(k, j) -= sum;
		}

		for (int i = k + 1; i <= std::min(n - 1, k + matrix_bandwidth_); ++i)
		{
			double sum = 0.0;
			for (int t = std::max(0, i - matrix_bandwidth_); t < k; ++t) 
			{
				if (std::abs(t - i) <= matrix_bandwidth_ && std::abs(t - k) <= matrix_bandwidth_)
				{
					sum += Acopy(i, t) * Acopy(t, k);
				}
			}

			if (std::abs(Acopy(k, k)) < 1e-15)
			{
				return false;
			}

			Acopy(i, k) = (Acopy(i, k) - sum) / Acopy(k, k);
		}
	}

	std::vector<double> y(n, 0.0);

	for (int i = 0; i < n; ++i) {
		double sum = 0.0;
		for (int j = std::max(0, i - matrix_bandwidth_); j < i; ++j)
		{
			sum += Acopy(i, j) * y[j];
		}
		y[i] = global_b_[i] - sum;
	}

	for (int i = n - 1; i >= 0; --i)
	{
		double sum = 0.0;
		for (int j = i + 1; j <= std::min(n - 1, i + matrix_bandwidth_); ++j)
		{
			sum += Acopy(i, j) * solution_[j];
		}

		if (std::abs(Acopy(i, i)) < 1e-15) throw std::runtime_error("Zero diagonal element in U");

		solution_[i] = (y[i] - sum) / Acopy(i, i);
	}

	return true;
}

std::vector<double> FEMsolver::getSolution() const { return solution_; }

double FEMsolver::computeResidualNorm()
{
	int n = mesh_->getNumNodes();
	std::vector<double> Aq(n);

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Aq[j] += (*this) (i, j) * solution_[j];
		}
	}

	double diff = 0, bNorm = 0;

	for (int i = 0; i < n; i++)
	{
		double r = Aq[i] - global_b_[i];
		diff += r * r;
	}
	diff = sqrt(diff);

	for (int i = 0; i < n; i++)
	{
		bNorm += global_b_[i] * global_b_[i];
	}
	bNorm = sqrt(bNorm);

	return diff / bNorm;
}