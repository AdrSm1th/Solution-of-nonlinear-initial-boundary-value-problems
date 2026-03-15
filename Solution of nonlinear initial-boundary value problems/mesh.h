//mesh.h

#pragma once

#include <string>
#include <vector>

class Mesh
{
public:
	std::vector<double> nodes_;
	std::vector<std::vector<int>> elements_;
	//std::map<int, BoundaryCondition> boundary_nodes_;

	void readMesh(const std::string &filename);
};