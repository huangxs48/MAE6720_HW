#include <cmath>
#include <iostream>
#include <vector>

#include "./MatrixSolver/matrixsolver.h"
#include "./PDE/PDE.h"
#include "./Mesh/mesh.h"
#include "./PDESolver/PDESolver.h"
#include "./Output/output.h"

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

int main(){

  //First, read in parameters
  //(TODO: change this I/O part)
    int numericalscheme; //integer number representing numberical scheme
    double gamma_input; //input gamma
    double timestep, deltax; //time step, space step

    char response; //let user choose input gamma or dt
    int Nx_input; //number of inner points
    double x1max, x1min; //domain range

    //read in mesh related values 
    cout << "===== Define domain first. Need parameters to construct mesh =====\n";
    cout << "X1min = ? " <<endl;
    cin >> x1min;
    cout << "X1max = ?" <<endl;
    cin >> x1max;
    assert(x1max >= x1min && "x1max must be larger than or equal to x1min\n");
    cout << "number of inner points? \n";
    cin >> Nx_input;
    assert(Nx_input >= 0 && "Number of inner points must be larger than zero\n");

  //Initialize an Mesh object
    Mesh m1;
    Mesh *pmymesh = &m1;
    m1.NX1 = Nx_input;
    m1.x1max = x1max;
    m1.x1min = x1min;
    //Mesh::MeshInitialize() calculates index for inner, outter points
    m1.MeshInitialize(); 
    //Mesh::UniformMeshConst returns vector of coordinates.
    //(TODO: now only works for 1D, multidimension?)
    vector<double> x1coord = m1.UniformMeshConstructor(m1.x1min, m1.x1max, m1.NX1);
    //define loop limits
    int is = m1.is, ie=m1.ie, il=m1.il, iu=m1.iu;
    double dx = m1.dx1;

  //Initialize a HeatConduction object
  //HeatConduction is an derived class of PartialDifferentialEquation 
    HeatConduction u1;
    HeatConduction *pmyeq = &u1;
    u1.N = m1.NX1;

  //Read in FDE paramters
    cout << "==== Then specify FDE parameters ====" <<endl;
    cout << "Specify gamma or timestep, type G for specify gamma, type T for specify timestep"<<endl;
    cin >> response;
    if (response != 'G' && response != 'g' && response != 'T' && response != 't'){
      cout << "Specify gamma or timestep, type G for specify gamma, type T for specify timestep"<<endl;
      cin >> response;
    }
    if (response == 'G' || response == 'g'){
      cout << "what is gamma?" << endl;
      cin >> gamma_input;
      u1.gamma = gamma_input;
      timestep = u1.gamma * dx * dx;
    } else if (response == 'T' || response == 't'){
      cout << "what is timestep?" << endl;
      cin >> timestep;
      u1.gamma = timestep/(dx*dx);
    }
    double t_end = 2.0;

  //Get numerical scheme
    cout << "choose a scheme. 0 for FTCS, 1 for Crank-Nicolson, 2 for BTCS\n";
    cin >> numericalscheme;
    if (numericalscheme!=0 && numericalscheme!=1 && numericalscheme!=2){
      cout << "choose a scheme. 0 for FTCS, 1 for Crank-Nicolson, 2 for BTCS\n";
      cin >> numericalscheme;      
    }

  //print configuration for user
    printf("The code is configured as: \n");
    printf("dx = %g, gamma = %g, m1.NX1 = %d, timestep=%g\n", m1.dx1, u1.gamma, m1.NX1, timestep);


    //Set initial condition
    vector<double> u_init=u1.u_init;
    for (int i = il; i <= iu; i++){
      double u_0 = InitFunc(x1coord.at(i));
      u_init.push_back(u_0);
    }

    //Add PDE solver pointer, matrix solver pointer
    PDESolver *psol;
    MatrixSolver *pms;

    //vector scratches to store u in the integration loop
    vector<double> u_next;
    vector<double> u_now = u_init;

    //set integration time, iteration limits
    double t_now = 0.0;
    int itrlim = t_end/timestep;

    //set numerical scheme
    //numericalscheme = 0;
    Output *pout;


    //Choose a scheme outside of the integration loop, 
    //so "if (numericalscheme)" statement will not be processed every cycle
    if (numericalscheme == 0){
      string filename = "FTCS";
      printf("The scheme is FTCS, start evoling with time...\n");
      for (int itr=0; itr<=itrlim; itr++){
        PrintScreen(pmymesh, itr);
        psol->FTCS(pmyeq, pmymesh, u_now, u_next);
        u_now = u_next;
        u_next.clear();
        t_now += timestep;
        pmymesh->t_now = t_now;
      }
      pout->Write1DArray(pmymesh, pmyeq, itrlim, x1coord, u_now, filename);
    }else if (numericalscheme == 1){
      string filename = "CN";
      printf("The scheme is CrankNicolson, start evoling with time...\n");
      for (int itr=0; itr<=itrlim; itr++){
        PrintScreen(pmymesh, itr);
        psol->CrankNicolson(pmyeq, pmymesh, pms, u_now, u_next);
        u_now = u_next;
        u_next.clear();
        t_now += timestep;
        pmymesh->t_now = t_now;
      }
      pout->Write1DArray(pmymesh, pmyeq, itrlim, x1coord, u_now, filename);
    }else if (numericalscheme == 2){
      string filename = "BTCS";
      printf("The scheme is BTCS, start evoling with time...\n");
      for (int itr=0; itr<=itrlim; itr++){
        PrintScreen(pmymesh, itr);
        psol->BTCS(pmyeq, pmymesh, pms, u_now, u_next);
        u_now = u_next;
        u_next.clear();
        t_now += timestep;
        pmymesh->t_now = t_now;
      }
      pout->Write1DArray(pmymesh, pmyeq, itrlim, x1coord, u_now, filename);
    }
  


 return 0;
}