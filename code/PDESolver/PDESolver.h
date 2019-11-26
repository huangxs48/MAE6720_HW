#ifndef PDESOLVER_H
#define PDESOLVER_H

#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>

#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../MatrixSolver/matrixsolver.h"
typedef std::vector<std::vector<double>> array2d;

class PDESolver{
  public:
    PDESolver();
    ~PDESolver();

    //For time marching parabolic equations
    std::vector<double> u_inner1d;

    void FTCS(PartialDifferentialEquation *peq, Mesh *pmesh, std::vector<double> &u_in, std::vector<double> &u_out);

    void BTCS(PartialDifferentialEquation *peq, Mesh *pmesh, 
              MatrixSolver *pms, std::vector<double> &u_in, std::vector<double> &u_out);

    void CrankNicolson(PartialDifferentialEquation *peq, Mesh *pmesh, MatrixSolver *pms, 
                  std::vector<double> &u_in, std::vector<double> &u_out);

    //For steady heat conduction, all are SOR method
    array2d u_inner2d; //inner points array at the end of every iteration, used in boundary function

    void pointJacobi(SteadyHeatConduction2d *peq, Mesh *pmesh, array2d &u_in, array2d &u_out, double omega);
    void pointGaussSiedel(SteadyHeatConduction2d *peq, Mesh *pmesh, array2d &u_in, array2d &u_out, double omega);

    void lineJacobi(SteadyHeatConduction2d *peq, Mesh *pmesh, MatrixSolver *pms, array2d &u_in, array2d &u_out, double omega);
    void lineGaussSiedel(SteadyHeatConduction2d *peq, Mesh *pmesh, MatrixSolver *pms, array2d &u_in, array2d &u_out, double omega);

};

#endif