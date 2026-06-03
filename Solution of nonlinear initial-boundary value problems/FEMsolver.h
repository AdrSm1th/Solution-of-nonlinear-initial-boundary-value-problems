//FEMsolver.h

#pragma once

#include <vector>
#include <iostream>
#include "Mesh.h"
#include "LocalMatrix.h"
#include "BasisFunctions.h"

class FEMsolver
{
private:
	Mesh *mesh_;
	LocalAssembler *assembler_;
	BasisFunction *basis_;
	std::vector<double> global_A_;
	std::vector<double> global_b_;
	std::vector<double> solution_;
	int matrix_bandwidth_;
	int index(int i, int j) const;
public:

	double &operator()(int i, int j)
	{
		int idx = index(i, j);
		if (idx == -1) throw std::out_of_range("Element outside band");
		return global_A_[idx];
	}
	const double &operator()(int i, int j) const
	{
		int idx = index(i, j);
		if (idx == -1) throw std::out_of_range("Element outside band");
		return global_A_[idx];
	}

	FEMsolver(Mesh &mesh, LocalAssembler &LA, BasisFunction &BS)
	{
		mesh_ = &mesh;
		assembler_ = &LA;
		basis_ = &BS;
		int n = mesh_->getNumNodes();
		matrix_bandwidth_ = 2;
		global_A_.resize(n * (2 * matrix_bandwidth_ + 1));
		global_b_.resize(n);
		solution_.resize(n);
	}

	void assembleGlobalSystem(const std::vector<double> &q_current);
	void assembleNewtonSystem(const std::vector<double> &q_current);

	void applyBoundaryCondition();

	bool solveLU();

	std::vector<double> getSolution() const;

	double computeResidualNorm();

	void assembleNewtonSystem(const std::vector<double> &q_current);

	void setSolution(const std::vector<double> &q);

	void printGlobalA() const;
	void printGlobalb() const;
	void printSolution() const;

};