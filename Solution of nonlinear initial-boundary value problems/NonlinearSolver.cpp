//NonlinearSolver.cpp

#include "NonlinearSolver.h"

inline std::vector<double> CmultV(double c, std::vector<double> &vector)
{
	std::vector<double> res(vector.size());
	for (int i = 0; i < vector.size(); i++)
	{
		res[i] = vector[i] * c;
	}
	return res;
}

inline std::vector<double> VplusV(std::vector<double> &vector1, std::vector<double> &vector2)
{
	std::vector<double> res(vector1.size());
	for (int i = 0; i < vector1.size(); i++)
	{
		res[i] = vector1[i] + vector2[i];
	}
	return res;
}

std::vector<double> NonlinearSolver::solveNonlinear()
{
	for (int k = 0; k < params_.max_iterations; k++)
	{
		fem_solver_->assembleGlobalSystem(q_prev_);
		fem_solver_->solveLU();
		if(fem_solver_->computeResidualNorm() < params_.tolerance) return q_curr_;
		q_curr_ = fem_solver_->getSolution();
		std::vector<double> qcw = CmultV(params_.omega, q_curr_);
		std::vector<double> qpw = CmultV(1 - params_.omega, q_prev_);
		q_curr_ = VplusV(qcw, qpw);
		q_prev_ = q_curr_;
	}

	return q_curr_;
}

void NonlinearSolver::setInitialGuess(const std::vector<double> &q0)
{
	q_prev_ = q0;
}
