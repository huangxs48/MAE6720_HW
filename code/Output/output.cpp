#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#include "output.h"
#include "../Mesh/mesh.h" //use t_now
#include "../PDE/PDE.h" //use gamma

using namespace std;

Output::Output(){
};

Output::~Output(){

}

void Output::Write1DArray(Mesh *pmesh, PartialDifferentialEquation *peq, int itr, vector<double> &coord, vector<double> &u, string username){

  string name = username;
  name += to_string(itr);
  name += ".csv";

  ofstream output (name);

  assert(output.is_open());

  output << "time = "<< to_string(pmesh->t_now) << " , gamma = " << to_string(peq->gamma) << " , dx = " << to_string(pmesh->dx1) << "\n";
  output << "coord,u" <<"\n";

  assert(u.size()==coord.size());
  for (int it=0; it<u.size(); it++){
    output << coord.at(it) << "," << u.at(it) << "\n";
  }  

  output.close();

  
}