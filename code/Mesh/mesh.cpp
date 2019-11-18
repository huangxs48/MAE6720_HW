#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>

#include "mesh.h"

/* Usage: specify NX1, NX2, NX3 and max and min in each direction, always Meshinitialize before use */


Mesh::Mesh(){
  //dim = 0;
  NX1 = 0, NX2 = 0, NX3 = 0;
  x1ratio = 1.0, x2ratio = 1.0, x3ratio = 1.0;
}

Mesh::~Mesh(){
}


void Mesh::MeshInitialize(){
  //assert(dim > 0);//must specify dimension
  t_now = 0.0;
  //set inner points index
  is = 1, js = 1, ks = 1;
  ie = 1, je = 1, ke = 1;
  //set coord index
  iu = 1, ju = 1, ju = 1;
  il = 0, jl = 0, kl = 0;
  if (NX1 > 0){
    ie = NX1;
    iu = NX1 + 1;
    assert(x1max > x1min);
    dx1 = (x1max - x1min)/(NX1+1);
  }
  if (NX2 > 0){
    je = NX2;
    ju = NX2 + 1;
    assert(x2max > x2min);
    dx2 = (x2max - x2min)/(NX2+1);
  }
  if (NX3 > 0){
    ke = NX3;
    ku = NX3 + 1;
    assert(x3max > x3min);
    dx3 = (x3max - x3min)/(NX3+1);
  }


}

//return a coordinate vector, inner points are is->ie, total points are il->iu
std::vector<double> Mesh::UniformMeshConstructor(double min, double max, int N){
  std::vector<double> coord;
  double dx = (max - min)/(N + 1);
  for (int index = 0; index <= N+1; index++){
    coord.push_back(min+dx*index);
  }
  return coord;
}





