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
   if(mesh.readMesh("mesh.txt"))
   {
      mesh.generateUniformMesh();
      BasisFunction basis;
      LocalAssembler assembler(mesh, basis);
      FEMsolver FEMsolver(mesh, assembler, basis);
      //NonlinearSolver nonlinearSolver(FEMsolver);
      std::vector<double> q_predict(mesh.getNumNodes(), 2);
      //for (int i = 0; i < mesh.getNumNodes(); i++)
      //{
      //   q_predict[i] = mesh.analyticalSolution(mesh.getNodeCoord(i)) + 100;
      //}
      //std::vector<double> q_predict{0, 0.2, 0.5, 0.7, 1};
      //std::vector<double> q_predict{0, 0.25, 0.5, 0.75, 1};
      //std::vector<double> q_predict{ 0, 0.5, 1 };
      //nonlinearSolver.setInitialGuess(q_predict);
      //std::vector<double> q = nonlinearSolver.solveNonlinear();
      //std::vector<double> q = nonlinearSolver.solveNewton();
      //FEMsolver.printSolution();

      for (int i = 0; i <= 10; i++) {
         NonlinearSolver nonlinearSolver(FEMsolver, 0.1 + 0.1 * i);
         nonlinearSolver.setInitialGuess(q_predict);
         std::vector<double> q = nonlinearSolver.solveNewton();
      }
   }
   return 0;
}