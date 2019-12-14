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

#define PI 3.14159265358979323846

//Initial Condition Function
double InitFunc(double x){
  return sin(x);
}

//Print time, cycle every step
void PrintScreen(Mesh *pmesh, int cycle){
  printf("step = %d, time = %g\n", cycle, pmesh->t_now);
}

//Must define boundary condition for four boundaries
double BoundaryCondition::bvalfunc_x1_inner(double xcoord, double ycoord, double input){
  return 0.0;
}
double BoundaryCondition::bvalfunc_x1_outter(double xcoord, double ycoord, double input){
  return 0.0;
}
double BoundaryCondition::bvalfunc_x2_inner(double xcoord, double ycoord, double input){
  return 0.0;
}
double BoundaryCondition::bvalfunc_x2_outter(double xcoord, double ycoord, double input){
  return 0.0;
}

int main(){

  //First, read in parameters
  //(TODO: change this I/O part)
    int numericalscheme; //integer number representing numberical scheme
    double gamma_input, c_input; //input gamma, c
    double timestep, deltax; //time step, space step

    char response; //let user choose input gamma or dt
    //int Nx_input = 50; //number of inner points
    double x1max = 1.0, x1min=0.0; //domain range

  //Initialize an Mesh object
    Mesh m1;
    Mesh *pmymesh = &m1;
    m1.NX1 = 89;
    m1.x1max = x1max;
    m1.x1min = x1min;
    //Mesh::MeshInitialize() calculates index for inner, outter points
    m1.MeshInitialize(); 
    //Mesh::UniformMeshConst returns vector of coordinates.
    vector<double> x1coord = m1.UniformMeshConstructor(m1.x1min, m1.x1max, m1.NX1);
    //define loop limits
    int is = m1.is, ie=m1.ie, il=m1.il, iu=m1.iu;
    double dx = m1.dx1;

  //Initialize a BurgersEquation object
  //Burgers Equation is an derived class of PartialDifferentialEquation 
    BurgersEquation u1;
    BurgersEquation *pmyeq = &u1;
    u1.N = m1.NX1;

    timestep = 0.001;
    deltax = m1.dx1;
    //double peclet = 
    c_input = 2.5;
    gamma_input = 0.05;

    u1.C = -c_input*timestep/deltax;
    u1.Gamma = gamma_input*timestep/(deltax*deltax);

    std::vector<double> u_init;
    for (int idx=0; idx< (iu-il); idx++){
      u_init.push_back(0.0);
    }
    u_init.push_back(1.0);
    std::vector<double> u_now = u_init;
    
    PDESolver *psol;
    Output *pout;

    int itr = 0;
    psol->SecondUpwind(pmyeq, pmymesh, u_now);
    
    /*
    for (int i=0; i<=100; i++){
      itr += 1;
      psol->FTCCS(pmyeq, pmymesh, u_now);
      printf("cycle = %d, res= %g\n", itr, pmyeq->Residual);
    }
    */
    

    while (pmyeq->Residual > 1.0e-6 ){
      psol->SecondUpwind(pmyeq, pmymesh, u_now);
      printf("cycle = %d, res= %g\n", itr, pmyeq->Residual);
      itr += 1;
    }  
    
  
    printf("peclet number: %g, mesh peclet number: %g\n", c_input/gamma_input, c_input*deltax/gamma_input);
    string namebase = "SUW89_";
    pout->Write1DArray(pmymesh, pmyeq, itr, x1coord, u_now, namebase);



 return 0;
}