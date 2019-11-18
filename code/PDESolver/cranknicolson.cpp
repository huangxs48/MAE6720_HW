#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>

#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../MatrixSolver/matrixsolver.h"
#include "./PDESolver.h"

/*
void printvector(std::vector<double> vec){
  std::cout<< "vec = ( ";
    for (int i=0; i<vec.size(); i++){
      std::cout<< vec.at(i) <<" ";
    }
  std::cout<<")\n";
}*/

void PDESolver::CrankNicolson(PartialDifferentialEquation *peq, Mesh *pmesh, MatrixSolver *pms, std::vector<double> &u_in, std::vector<double> &u_out){

  //u_in includes boundary points, u_out only has inner points, needs to add edge points manually, change this later
  int is = pmesh->is, ie = pmesh->ie;
  int iu = pmesh->iu, il = pmesh->il;
  double gamma = peq->gamma;

  assert(peq != NULL);
  assert(gamma <= 0.5);
  assert(u_in.size() == (iu - il +1));

  //prepare arrays, a,b,c,x,y have boundary points
  size_t N = ie-is+1;
  std::vector<double> a(N-1, -gamma/2.0);
  a.insert(a.begin(),0.0);
  std::vector<double> b(N, 1+gamma);
  std::vector<double> c(N-1, -gamma/2.0);
  c.push_back(0.0);

  std::vector<double> y(N, 0.0);
  std::vector<double> x(N, 0.0);

  //RHS of tridiagonal matrix equation
  for (int i=0; i<N; i++){
    x[i] = u_in.at(is+i);
    y[i] = gamma*0.5*u_in[is+i+1] + (1.0-gamma)*u_in[is+i] + gamma*0.5*u_in[is+i-1];
  }

  //call TDMA
  pms->TDMA(a, b, c, y, x);

  //prepare output, apply boundaries manually, need to change
  u_out = x;
  u_out.insert(u_out.begin(),0.0);
  u_out.push_back(0.0);


}

