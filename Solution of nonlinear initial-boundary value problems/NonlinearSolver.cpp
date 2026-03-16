//NonlinearSolver.cpp

#include <iostream>
#include "NonlinearSolver.h"

std::vector<double> NonlinearSolver::solveNonlinear()
{
   int n = q_prev_.size();
   q_curr_.resize(n);

   for (int k = 0; k < params_.max_iterations; k++)
   {
      fem_solver_->assembleGlobalSystem(q_prev_);
      fem_solver_->applyBoundaryCondition();

      if (!fem_solver_->solveLU())
      {
         std::cerr << "LU decomposition failed at iteration " << k << std::endl;
         return q_prev_;
      }

      std::vector<double> q_new = fem_solver_->getSolution();

      for (int i = 0; i < n; i++)
      {
         q_curr_[i] = params_.omega * q_new[i] + (1.0 - params_.omega) * q_prev_[i];
      }

      fem_solver_->assembleGlobalSystem(q_curr_);
      fem_solver_->applyBoundaryCondition();
      double residual = fem_solver_->computeResidualNorm();

      std::cout << "Iteration " << k << ", residual = " << residual << std::endl;

      if (residual < params_.tolerance)
      {
         std::cout << "Converged after " << k + 1 << " iterations." << std::endl;
         return q_curr_;
      }

      q_prev_ = q_curr_;
   }
	return q_curr_;
}

void NonlinearSolver::setInitialGuess(const std::vector<double> &q0)
{
	q_prev_ = q0;
}
