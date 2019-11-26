#include <iostream>
#include <vector>
#include <assert.h>

#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../MatrixSolver/matrixsolver.h"
#include "./PDESolver.h"
#include "../BoundaryCondition/bval.h"


typedef std::vector<std::vector<double>> array2d;

void PDESolver::lineGaussSiedel(SteadyHeatConduction2d *peq, Mesh *pmesh, 
								                MatrixSolver *pms, array2d &u_in, array2d &u_out, double omega)
{
  BoundaryCondition *pbval;
  int is = pmesh->is, ie = pmesh->ie;
  int iu = pmesh->iu, il = pmesh->il;
  int js = pmesh->js, je = pmesh->je;
  int ju = pmesh->ju, jl = pmesh->jl;

  int NX1 = pmesh->NX1, NX2=pmesh->NX2;
  size_t N = NX1*NX2;

  //////////////////////////////////////////////////////////
  // i = is, including Neumann boundary
  // | 3/4   -1/4  0  ...|   |         |   |(1/4)*u(il,js)|
  // |-1/4  3/4  -1/4 ...|   |         |   |              |
  // | 0   -1/4   3/4 ...| * |  u_new  | = |   u_inner    |
  // |       ..      ....|   |         |   |              |
  // |                   |   |         |   |(1/4)*u(ie,ju)|

  // i = is+1->ie
  // |  1   -1/4  0 ...|   |         |   |(1/4)*u(il,js)|
  // |-1/4   1  -1/4...|   |         |   |              |
  // | 0   -1/4   1 ...| * |  u_new  | = |   u_inner    |
  // |       ..    ....|   |         |   |              |
  // |                 |   |         |   |(1/4)*u(ie,ju)|
  //////////////////////////////////////////////////////////

  //u_in has boundary points
  //prepares TDMA LHS array
  std::vector<double> a(NX1-1, -0.25*omega);
  a.insert(a.begin(), 0.0);
  std::vector<double> b(NX1, 1.0);
  std::vector<double> c(NX1-1, -0.25*omega);
  c.push_back(0.0);

  //for the first row
  std::vector<double> as(NX1-1, -0.25*omega);
  as.insert(as.begin(), 0.0);
  std::vector<double> bs(NX1, 1.0-(0.25*omega));
  std::vector<double> cs(NX1-1, -0.25*omega);
  cs.push_back(0.0);
  
  //sweep rows
  double res = 0.0;
  for (int j=js; j<=je; j++){
    std::vector<double> y(NX1, 0.0);
    std::vector<double> ys(NX1, 0.0);
    std::vector<double> x(NX1, 0.0);
    
    if (j == js){

      for (int i=0; i<NX1; i++){
        ys[i] = 0.25 * omega * (u_in.at(j+1).at(i+is)) + (1 - omega)*u_in.at(j).at(i+is);
      }
      pms->TDMA(as, bs, cs, ys, x);

    }else{

      for (int i=0; i<NX1; i++){
        y[i] = 0.25 * omega * (u_in.at(j-1).at(i+is) + u_in.at(j+1).at(i+is)) + (1 - omega)*u_in.at(j).at(i+is);
      }
      pms->TDMA(a, b, c, y, x);
    }

    //sanity check
    assert(x.size()==NX1);
    //add boundary points
    double x_inner = pbval->bvalfunc_x1_inner(il, j, u_in.at(j).at(is));
    double x_outter = pbval->bvalfunc_x1_outter(iu, j, u_in.at(j).at(ie));
    x.insert(x.begin(), x_inner);
    x.push_back(x_outter);

    //apply SOR and calculate residual 
    for(int i=0; i<x.size(); i++){
      
      res += fabs(x.at(i) - u_in.at(j).at(i));
    }

    //update this row
    u_in.at(j) = x;

  }

  pbval->Apply2DBval(pmesh, u_in);
  u_out = u_in;
  peq->Residual = res/(N);

}