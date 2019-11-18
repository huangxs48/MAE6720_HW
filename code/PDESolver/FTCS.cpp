#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>

#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../MatrixSolver/matrixsolver.h"
#include "./PDESolver.h"

void PDESolver::FTCS(PartialDifferentialEquation *peq, Mesh *pmesh, std::vector<double> &u_in, std::vector<double> &u_out){
  //u_in includes boundaries points, u_out only has inner points, needs to add edge points manually, change this later
  int is = pmesh->is, ie = pmesh->ie;
  int iu = pmesh->iu, il = pmesh->il;
  double gamma = peq->gamma;

  assert(peq != NULL);
  assert(gamma <= 0.5);
  assert(u_in.size() == (iu - il +1));

  //apply left boundary
  u_out.push_back(0.0);

  //calculate inner points
  for (int i=is; i<=ie; i++){
    double u_nx = u_in.at(i);
    double u_nxl = u_in.at(i-1);
    double u_nxr = u_in.at(i+1);
    double u_n1 = u_nx + gamma*(u_nxl + u_nxr - 2*u_nx);
    u_out.push_back(u_n1);
  }

  //apply right boundary
  u_out.push_back(0.0);
 
  //check size of u_out
  assert(u_out.size() == u_in.size());
  //printf("ftcs: uout:%d uin:%d\n", u_out.size(), u_in.size());

  //return u_out;

}