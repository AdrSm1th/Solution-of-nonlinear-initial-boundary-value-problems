//main.cpp

#include <iostream>
#include <vector>
#include "Mesh.h"
#include "BasisFunctions.h"
#include "LocalMatrix.h"
#include "FEMsolver.h"
#include "NonlinearSolver.h"

int main()
{
   Mesh mesh;
   mesh.readMesh("mesh.txt");
   mesh.generateUniformMesh();
   BasisFunction basis;
   LocalAssembler assembler(mesh, basis);
   FEMsolver FEMsolver(mesh, assembler, basis);
   NonlinearSolver nonlinearSolver(FEMsolver);
   std::vector<double> q_predict(mesh.getNumNodes(), 0.1);
   //std::vector<double> q_predict{0, 0.2, 0.5, 0.7, 1};
   //std::vector<double> q_predict{0, 0.25, 0.5, 0.75, 1};
   //std::vector<double> q_predict{ 0, 0.5, 1 };
   nonlinearSolver.setInitialGuess(q_predict);
   std::vector<double> q = nonlinearSolver.solveNonlinear();
   for (int i = 0; i < q.size(); i++)
   {
      std::cout << q[i] << std::endl;
   }

   return 0;
}