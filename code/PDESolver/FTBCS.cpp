#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>

#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../MatrixSolver/matrixsolver.h"
#include "./PDESolver.h"

void PDESolver::FTBCS(BurgersEquation *peq, Mesh *pmesh, std::vector<double> &u){

  int is = pmesh->is, ie = pmesh->ie;
  int iu = pmesh->iu, il = pmesh->il;
  double Gamma = peq->Gamma;
  double C = peq->C;

  assert(peq != NULL);
  assert(u.size() == (iu - il +1));

  double res = 0.0;
  //calculate inner points
  for (int i=is; i<=ie; i++){
    double u_nx = u.at(i);
    double u_nxl = u.at(i-1);
    double u_nxr = u.at(i+1);
    double u_n1 = u_nx + Gamma*(u_nxl + u_nxr - 2*u_nx) + C*(u_nx - u_nxl);
    res += fabs(u_n1 - u_nx);
    u.at(i) = u_n1;
  }

  peq->Residual = res/(pmesh->NX1);

}