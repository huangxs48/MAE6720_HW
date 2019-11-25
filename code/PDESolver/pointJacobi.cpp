#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>

#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../MatrixSolver/matrixsolver.h"
#include "./PDESolver.h"
#include "../BoundaryCondition/bval.h"


typedef std::vector<std::vector<double>> array2d;

void PDESolver::pointJacobi(SteadyHeatConduction2d *peq, Mesh *pmesh, array2d &u_in, array2d &u_out, double omega){
  BoundaryCondition *pbval;

  //u_in includes boundary points
  int is = pmesh->is, ie = pmesh->ie;
  int iu = pmesh->iu, il = pmesh->il;
  int js = pmesh->js, je = pmesh->je;
  int ju = pmesh->ju, jl = pmesh->jl;

  assert(peq != NULL);
  assert(u_in.size() == (ju - jl +1));
  assert(u_in.at(0).size() == (iu - il +1));

  array2d u_;//scratch array, only has inner points

  double res=0.0;
  for (int j=js; j<=je; j++){
    std::vector<double> row_;
    if(j==js){//incorporate BC into the equations
      for (int i=is; i<=ie; i++){
        double u_current = u_in.at(j).at(i);
        double sum_ = (1.0/3.0)*(u_in.at(j).at(i+1)+u_in.at(j).at(i-1)+u_in.at(j+1).at(i));
        double sum = omega*sum_ + (1.0-omega)*u_current;
        res += fabs(u_current-sum);
        row_.push_back(sum);
      }
      u_.push_back(row_);
    }else{
      for (int i=is; i<=ie; i++){
        double u_current = u_in.at(j).at(i);
        double sum_ = (1.0/4.0)*(u_in.at(j).at(i+1)+u_in.at(j).at(i-1)+u_in.at(j+1).at(i)+u_in.at(j-1).at(i));
        double sum = omega*sum_ + (1.0-omega)*u_current;
        res += fabs(u_current-sum);
        row_.push_back(sum);
      }
      u_.push_back(row_);
      
    }
  }

  pbval->Add2DBval(pmesh, u_, u_out);
  res = res/((je-is+1)*(ie-is+1));
  peq->Residual = res;

}