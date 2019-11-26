#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <assert.h>

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

void Output::Write2DArray(Mesh *pmesh, PartialDifferentialEquation *peq, int itr, array2d &u, std::string username){

  string name = username;
  name += to_string(itr);
  name += ".dat";

  ofstream output (name);
  assert(output.is_open());

  for (int j=0; j<u.size(); j++){
    for (int i=0; i<u.at(j).size()-1; i++){
      output << u.at(j).at(i) <<",";
    }
    output << u.at(j).back()<< "\n";
  }


}