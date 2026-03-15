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
	nodes_.resize(n_);
	elements_.resize(n_ - 2, std::vector<int>(0));

	for (int i = 0; i < n_ - 2; i++)
	{
		elements_[i].push_back(i);
		elements_[i].push_back(i + 1);
		elements_[i].push_back(i + 2);
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

int Mesh::getNumElems() const { return n_ - 2; }

double Mesh::getNodeCoord(int node_id) const { return nodes_[node_id]; }

std::vector<int> Mesh::getElementNodes(int elem_id) const { return elements_[elem_id]; }

bool Mesh::isBoundary(int node_id) const { return node_id == 0 || node_id == n_ - 1; }

BoundaryCondition Mesh::getBoundaryCondition(bool left) const { return left ? left_ : right_; }

double Mesh::lambda(double u) const { return u; }

double Mesh::f(double x) const { return x; }

double Mesh::gamma() { return gamma_; }