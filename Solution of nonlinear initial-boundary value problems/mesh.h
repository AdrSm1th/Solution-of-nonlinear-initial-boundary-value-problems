//Mesh.h

#pragma once

#include <string>
#include <vector>

class BoundaryCondition
{
public:
	int type;

	double u_g;
	double theta;
	double beta;
	double u_beta;
};

class Mesh
{
private:
	double a_, b_, gamma_;
	int n_;
	std::vector<double> nodes_;
	std::vector<std::vector<int>> elements_;
	BoundaryCondition left_, right_;

public:
	Mesh() : a_(0.0), b_(0.0), n_(0), gamma_(0.0), left_(), right_(){}
	bool readMesh(const std::string &filename);
	void generateUniformMesh();
	int genNumNodes() const;
	int genNumElems() const;
	double getNodeCoord(int node_id) const;
	std::vector<int> getElementNodes(int elem_id) const;
	bool isBoundary(int node_id) const;
	BoundaryCondition getBoundaryCondition(bool left) const;
	double f(double x) const;
	double lambda(double u) const;
	double gamma();
};