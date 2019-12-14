#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>

#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../MatrixSolver/matrixsolver.h"
#include "./PDESolver.h"

void PDESolver::SecondUpwind(BurgersEquation *peq, Mesh *pmesh, std::vector<double> &u){

  int is = pmesh->is, ie = pmesh->ie;
  int iu = pmesh->iu, il = pmesh->il;
  double Gamma = peq->Gamma;
  double C = peq->C;

  assert(peq != NULL);
  assert(u.size() == (iu - il +1));

  double res = 0.0;
  std::vector<double> u_adv_; //vector to store intermediate step values
  u_adv_.push_back(0.0); //upwind, apply the left boundary first

  //calculate inner points, do advection intermediate step first
  for (int i=is; i<=iu; i++){ //here since it's upwind, we can update until the boundary point
    double u_nx = u.at(i);
    double u_nxl = u.at(i-1);
    double u_nx_adv_ = 0.0;

    //intermediate step for advection
    if (i==is){ //the first inner point can only do first order
      u_nx_adv_ = u_nx + C*(u_nx - u_nxl);
    }else{ //other points use second order
      double u_nxll = u.at(i-2);
      u_nx_adv_ = u_nx + (C/2.0) * (3.0*u_nx - 4.0*u_nxl + u_nxll);
    }
    u_adv_.push_back(u_nx_adv_);
  }
  
  assert(u_adv_.size()==u.size());

  //do the full step updating including advection and diffusion
  for (int i=is; i<=ie; i++){

    //points from last step
    double u_nx = u.at(i);
    double u_nxl = u.at(i-1);
    double u_nxr = u.at(i+1);
    //points from intermediate step
    double u_nx_ = u_adv_.at(i);
    double u_nxl_ = u_adv_.at(i-1);
    double adv_term = 0.0;

    //full step for both advection and diffusion
    if (i==is){
      adv_term = 0.5*(u_nx + u_nx_ + C*(u_nx_ - u_nxl_));
    } else{
      double u_nxll_ = u_adv_.at(i-2);
      adv_term = 0.5*(u_nx + u_nx_ + (C/2.0) * (3.0*u_nx_ - 4.0*u_nxl_ + u_nxll_));
    }
    double u_n1 = adv_term + Gamma*(u_nxl + u_nxr - 2*u_nx);
    res += fabs(u_n1 - u_nx);
    u.at(i) = u_n1;
  }
  //calculate residual
  peq->Residual = res/(pmesh->NX1);

}