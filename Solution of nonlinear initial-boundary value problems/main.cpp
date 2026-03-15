//main.cpp

#include <iostream>
#include "Mesh.h"

int main()
{
   Mesh mesh;
   mesh.readMesh("mesh.txt");
   mesh.generateUniformMesh();

   return 0;
}