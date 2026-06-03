	//Mesh.cpp

	#include <fstream>
	#include <iostream>
	#include "Mesh.h"

	bool Mesh::readMesh(const std::string& filename)
	{
		std::ifstream input(filename);
		input >> a_ >> b_ >> n_ >> gamma_;
		if (n_ < 3)
		{
			std::cout << "number of nodes must be more than 3\n";
			return false;
		}
		if (n_ % 2 == 0)
		{
			std::cout << "number of nodes must be odd\n";
			return false;
		}

		nodes_.resize(n_);

		int idx = 0;
		while(idx + 2 < n_)
		{
			elements_.push_back({ idx, idx + 1, idx + 2 });
			idx += 2;
		}

		input >> left_.type >> right_.type;

		switch (left_.type)
		{
			case 1:
			{
				input >> left_.u_g;
				break;
			}
			case 2:
			{
				input >> left_.theta;
				break;
			}
			case 3:
			{
				input >> left_.beta >> left_.u_beta;
				break;
			}
			default:
			{
				std::cout << "Invalid input: boundary condition must be 1, 2 or 3\n";
				return false;
			}
		}

		switch (right_.type)
		{
			case 1:
			{
				input >> right_.u_g;
				break;
			}
			case 2:
			{
				input >> right_.theta;
				break;
			}
			case 3:
			{
				input >> right_.beta >> right_.u_beta;
				break;
			}
			default:
			{
				std::cout << "Invalid input: boundary condition must be 1, 2 or 3\n";
				return false;
			}
		}

		return true;
	}

	void Mesh::generateUniformMesh()
	{
		double h = (b_ - a_) / (n_ - 1);

		for (int i = 0; i < n_ - 1; i++)
		{
			nodes_[i] = a_ + h * i;
		}
		nodes_[n_ - 1] = b_;
	}

	int Mesh::getNumNodes() const { return n_; }

	int Mesh::getNumElems() const { return elements_.size(); }

	double Mesh::getNodeCoord(int node_id) const { return nodes_[node_id]; }

	std::vector<int> Mesh::getElementNodes(int elem_id) const { return elements_[elem_id]; }

	bool Mesh::isBoundary(int node_id) const { return node_id == 0 || node_id == n_ - 1; }

	BoundaryCondition Mesh::getBoundaryCondition(bool left) const { return left ? left_ : right_; }

	double Mesh::lambda(double u) const { return log(u); }

	double Mesh::lambdaDerivative(double u) const { return 1 / u; }

	double Mesh::f(double x) const { return -(1 / x) + gamma_ * x; }

	double Mesh::gamma() { return gamma_; }

	double Mesh::analyticalSolution(double x) const { return x; }