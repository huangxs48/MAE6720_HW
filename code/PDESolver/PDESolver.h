#ifndef PDESOLVER_H
#define PDESOLVER_H

#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>

#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../MatrixSolver/matrixsolver.h"

class PDESolver{
  public:
    PDESolver();
    ~PDESolver();

    void FTCS(PartialDifferentialEquation *peq, Mesh *pmesh, std::vector<double> &u_in, std::vector<double> &u_out);

    void BTCS(PartialDifferentialEquation *peq, Mesh *pmesh, 
              MatrixSolver *pms, std::vector<double> &u_in, std::vector<double> &u_out);

    void CrankNicolson(PartialDifferentialEquation *peq, Mesh *pmesh, MatrixSolver *pms, 
    		      std::vector<double> &u_in, std::vector<double> &u_out);

};

#endif