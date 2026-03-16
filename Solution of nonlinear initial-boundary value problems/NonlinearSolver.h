//NonlinearSolver.h

#pragma once

#include <vector>
#include "FEMsolver.h"

struct SolverParams
{
	double tolerance = 1e-16;
	int max_iterations = 500;
	double omega = 1;
};

class NonlinearSolver
{
private:
	FEMsolver *fem_solver_;
	SolverParams params_;
	std::vector<double> q_prev_;
	std::vector<double> q_curr_;

public:
	NonlinearSolver(FEMsolver &solver)
	{
		fem_solver_ = &solver;
	}
	std::vector<double> solveNonlinear();
	void setInitialGuess(const std::vector<double> &q0);
};