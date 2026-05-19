//FEMsolver.cpp

#include <fstream>
#include <iomanip>
#include "FEMsolver.h"

int FEMsolver::index(int i, int j) const
{
	if (std::abs(i - j) > matrix_bandwidth_) return -1;

	int diag_offset = i - j;

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

void FEMsolver::assembleNewtonSystem(const std::vector<double> &q_current)
{
	std::fill(global_A_.begin(), global_A_.end(), 0.0);
	std::fill(global_b_.begin(), global_b_.end(), 0.0);

	for (int elem = 0; elem < mesh_->getNumElems(); ++elem)
	{
		LocalMatrices matrices = assembler_->newtonIteration(elem, q_current);
		std::vector<int> elem_nodes = mesh_->getElementNodes(elem);

		for (int i = 0; i < 3; ++i)
		{
			int i_global = elem_nodes[i];
			for (int j = 0; j < 3; ++j)
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
				global_b_[idx] += cond.theta;
				break;
			}

			case 3:
			{
				(*this)(idx, idx) += cond.beta;
				global_b_[idx] += cond.beta * cond.u_beta;
				break;
			}
		}
	}
}

bool FEMsolver::solveLU()
{
	//std::vector<double> matrix_copy = global_A_;

	int n = mesh_->getNumNodes();
	for (int k = 0; k < n; ++k) 
	{
		for (int j = k; j <= std::min(k + matrix_bandwidth_, n - 1); ++j) 
		{
			double sum = 0.0;
			for (int t = std::max(0, k - matrix_bandwidth_); t < k; ++t)
			{
				if (std::abs(t - k) <= matrix_bandwidth_ && std::abs(t - j) <= matrix_bandwidth_) 
				{
					sum += (*this)(k, t) * (*this)(t, j);
				}
			}
			(*this)(k, j) -= sum;
		}

		for (int i = k + 1; i <= std::min(k + matrix_bandwidth_, n - 1); ++i)
		{
			double sum = 0.0;
			for (int t = std::max(0, i - matrix_bandwidth_); t < k; ++t) 
			{
				if (std::abs(t - i) <= matrix_bandwidth_ && std::abs(t - k) <= matrix_bandwidth_)
				{
					sum += (*this)(i, t) * (*this)(t, k);
				}
			}

			if (std::abs((*this)(k, k)) < 1e-15)
			{
				return false;
			}

			(*this)(i, k) = ((*this)(i, k) - sum) / (*this)(k, k);
		}
	}

	std::vector<double> y(n, 0.0);

	for (int i = 0; i < n; ++i) {
		double sum = 0.0;
		for (int j = std::max(0, i - matrix_bandwidth_); j < i; ++j)
		{
			sum += (*this)(i, j) * y[j];
		}
		y[i] = global_b_[i] - sum;
	}

	for (int i = n - 1; i >= 0; --i)
	{
		double sum = 0.0;
		for (int j = i + 1; j <= std::min(n - 1, i + matrix_bandwidth_); ++j)
		{
			sum += (*this)(i, j) * solution_[j];
		}

		if (std::abs((*this)(i, i)) < 1e-15) throw std::runtime_error("Zero diagonal element in U");

		solution_[i] = (y[i] - sum) / (*this)(i, i);
	}

	//global_A_ = matrix_copy;

	return true;
}

std::vector<double> FEMsolver::getSolution() const { return solution_; }

double FEMsolver::computeResidualNorm()
{
	int n = mesh_->getNumNodes();
	std::vector<double> Aq(n, 0);

	for (int i = 0; i < n; i++)
	{
		int start_j = i - 2;
		if (start_j < 0) start_j = 0;
		for (int j = start_j; j <= std::min(n - 1, i + 2); j++)
		{
			Aq[i] += (*this) (i, j) * solution_[j];
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

void FEMsolver::printGlobalA() const
{
	int n = mesh_->getNumNodes();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (abs(j - i) > matrix_bandwidth_) std::cout << 0 << " ";
			else std::cout << (*this)(i, j) << " ";
		}
		std::cout << std::endl;
	}
}

void FEMsolver::printGlobalb() const
{
	int n = mesh_->getNumNodes();
	for (int i = 0; i < n; i++)
	{
		std::cout << global_b_[i] << std::endl;
	}
}

void FEMsolver::printSolution() const
{
	std::ofstream output("output.txt");
	output << std::setprecision(16) << std::fixed;
	output << "x" << std::setw(29) << "q" << std::setw(30) << "q*" << std::setw(30) << "q - q*\n";
	for (int i = 0; i < mesh_->getNumNodes(); i++)
	{
		double q = mesh_->analyticalSolution(mesh_->getNodeCoord(i));
		output << mesh_->getNodeCoord(i) << std::setw(29) << q << std::setw(29) << solution_[i] << std::setw(25) << q - solution_[i] << std::endl;
	}
}