/* This is a routine for solving steady 2d laplace equation */

#include <cmath>
#include <iostream>
#include <vector>

#include "./MatrixSolver/matrixsolver.h"
#include "./PDE/PDE.h"
#include "./Mesh/mesh.h"
#include "./PDESolver/PDESolver.h"
#include "./Output/output.h"
#include "./BoundaryCondition/bval.h"

using namespace std;
typedef std::vector<std::vector<double>> array2d;

#define PI 3.14159265358979323846

//Must define boundary condition for four boundaries
double BoundaryCondition::bvalfunc_x1_inner(double xcoord, double ycoord, double input){
  return 0.0;
}
double BoundaryCondition::bvalfunc_x1_outter(double xcoord, double ycoord, double input){
  return 0.0;
}
double BoundaryCondition::bvalfunc_x2_inner(double xcoord, double ycoord, double input){
  return input;
}
double BoundaryCondition::bvalfunc_x2_outter(double xcoord, double ycoord, double input){
  return sin(PI*xcoord);
}

//Print time, cycle every step
void PrintScreen(Mesh *pmesh, int cycle){
  printf("step = %d, time = %g\n", cycle, pmesh->t_now);
}

void print2darray(array2d intput){
  printf("-------------\n");
  for (int i=0; i<intput.size(); i++){
    for (int j=0; j<intput.at(i).size(); j++){
      printf(" %g", intput.at(i).at(j));
    }
    printf("\n");
  }
}

array2d unitarray(int col, int row){
  array2d unit;
  for (int j=0; j<col; j++){
    std::vector<double> row_;
    for (int i=0; i<row; i++){
      if (i==j){
        row_.push_back(2.0);
      }else{
        row_.push_back(1.0);
      }
    }
    unit.push_back(row_);
  }
  return unit;
}

array2d initguess(Mesh *pmesh){
  array2d guess;
  double factor = 0.01;
  for (int j=pmesh->jl; j<=pmesh->ju; j++){
    std::vector<double> this_row;
    for (int i=pmesh->il; i<=pmesh->iu; i++){
      this_row.push_back((factor*(j))*sin(PI*pmesh->coord_NX1.at(i)));
    }
    guess.push_back(this_row);
  }
  return guess;
}

double checkerror(Mesh *pmesh, array2d u_next){

    int is=pmesh->is, ie=pmesh->ie;
    int js=pmesh->js, je=pmesh->je;
    int il=pmesh->il, iu=pmesh->iu;
    int jl=pmesh->jl, ju=pmesh->ju;

    //array2d alytclsol;
    double diff = 0;
    
    for (int j=js; j<=je; j++){
      //std::vector<double> this_row;
      for (int i=is; i<=ie; i++){
        double this_x = pmesh->MeshGrid.at(j).at(i)[0];
        double this_y = pmesh->MeshGrid.at(j).at(i)[1];
        double this_val = sin(PI*this_x)*cosh(PI*this_y)/cosh(PI);
        diff += fabs(this_val - u_next.at(j).at(i));
        //this_row.push_back(this_val);
      }
      //alytclsol.push_back(this_row);
    }

    diff /= (pmesh->NX1)*(pmesh->NX2);
    return diff;
    
}


int main(){

    BoundaryCondition *pbval;
    PDESolver *psol;
    Output *pout;
    MatrixSolver *pms;
    double omega = 1.0;

  //Initialize an Mesh object
    Mesh m2;
    Mesh *pmymesh = &m2;
    m2.NX1 = 39;
    m2.x1max = 1.0;
    m2.x1min = 0.0;

    m2.NX2 = 39;
    m2.x2max = 1.0;
    m2.x2min = 0.0;
    m2.MeshInitialize(); 

    int is=pmymesh->is, ie=pmymesh->ie;
    int js=pmymesh->js, je=pmymesh->je;
    int il=pmymesh->il, iu=pmymesh->iu;
    int jl=pmymesh->jl, ju=pmymesh->ju;

    //Make the mesh grid
    pmymesh->Uniform2dMeshConstructor(m2.x1min, m2.x1max, m2.x2min, m2.x2max, m2.NX1, m2.NX2);

    //Initialize equation object
    SteadyHeatConduction2d my_eq;
    SteadyHeatConduction2d *pmyeq = &my_eq;
    my_eq.Nx = pmymesh->NX1;
    my_eq.Ny = pmymesh->NX2;

    
    array2d u_init=initguess(pmymesh);
    pbval->Apply2DBval(pmymesh, u_init);

    array2d u_next;
    psol->lineGaussSiedel(pmyeq, pmymesh, pms, u_init, u_next, omega);
    u_init = u_next;
    u_next.clear();       
    int itr = 1;   

    printf("itr,res_before,res\n");   
    while (pmyeq->Residual>=1.e-6){
      double res_before = pmyeq->Residual;
      u_next.clear();
      psol->lineGaussSiedel(pmyeq, pmymesh, pms, u_init, u_next, omega);
      u_init = u_next;
      itr += 1; 
      double res_after = pmyeq->Residual;
      printf("%d,%g\n", itr,res_after);
    }

    printf("itr=%d, error=%g, residual=%g\n", itr, checkerror(pmymesh, u_next), pmyeq->Residual);  

    
    //User after works
    array2d error_arr;
    for (int j=jl; j<=ju; j++){
      std::vector<double> this_row;
      for (int i=il; i<=iu; i++){
        double diff;
        double this_x = pmymesh->MeshGrid.at(j).at(i)[0];
        double this_y = pmymesh->MeshGrid.at(j).at(i)[1];
        double this_val = sin(PI*this_x)*cosh(PI*this_y)/cosh(PI);
        diff = fabs(this_val - u_next.at(j).at(i));
        this_row.push_back(diff);
      }
      error_arr.push_back(this_row);
    }

    std::string name = "PJacobi";
    pout->Write2DArray(pmymesh, pmyeq, itr, error_arr, name);
    
    
    /*
    for(int itr =0; itr<=800; itr++){
      u_next.clear();
      psol->lineGaussSiedel(pmyeq, pmymesh, pms, u_init, u_next);
      u_init = u_next;
      //print2darray(u_init);
      printf("itr=%d, error=%g, residual=%g\n", itr, checkerror(pmymesh, u_next), pmyeq->Residual);
    }


    array2d error_arr;
    for (int j=jl; j<=ju; j++){
      std::vector<double> this_row;
      for (int i=il; i<=iu; i++){
        double diff;
        double this_x = pmymesh->MeshGrid.at(j).at(i)[0];
        double this_y = pmymesh->MeshGrid.at(j).at(i)[1];
        double this_val = sin(PI*this_x)*cosh(PI*this_y)/cosh(PI);
        diff = abs(this_val - u_next.at(j).at(i));
        this_row.push_back(diff);
      }
      error_arr.push_back(this_row);
    }
    
    //std::string name = "GS";
    //pout->Write2DArray(pmymesh, pmyeq, itr, error_arr, name);
    
    */
    

 return 0;
}