#include <iostream>
#include <vector>
#include <assert.h>

#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../MatrixSolver/matrixsolver.h"
#include "./PDESolver.h"
#include "../BoundaryCondition/bval.h"


typedef std::vector<std::vector<double>> array2d;