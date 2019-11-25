#ifndef BVAL_H
#define BVAL_H

#include <cmath>
#include <iostream>
#include <vector>

#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../PDESolver/PDESolver.h"

typedef std::vector<std::vector<double>> array2d;

class BoundaryCondition{
public:
  BoundaryCondition();
  ~BoundaryCondition();

  void Apply1DBval(Mesh *pmesh, std::vector<double> &u_in, std::vector<double> &u_out);
  void Apply2DBval(Mesh *pmesh, array2d &u_in);
  void Add2DBval(Mesh *pmesh, array2d &u_in, array2d &u_out);

  //(TODO) need to rewrite boundaries conditions, maybe use pmesh and peq
  double bvalfunc_x1_inner(double xcoord, double ycoord, double input);
  double bvalfunc_x1_outter(double xcoord, double ycoord, double input);
  double bvalfunc_x2_inner(double xcoord, double ycoord, double input);
  double bvalfunc_x2_outter(double xcoord, double ycoord, double input);

};

#endif