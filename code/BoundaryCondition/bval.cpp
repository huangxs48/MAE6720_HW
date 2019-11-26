#include <cmath>
#include <iostream>
#include <vector>

#include "./bval.h"
#include "../PDE/PDE.h"
#include "../Mesh/mesh.h"
#include "../MatrixSolver/matrixsolver.h"

typedef std::vector<std::vector<double>> array2d;

//(TODO) need to re-think this part

void BoundaryCondition::Add2DBval(Mesh *pmesh, array2d &u_in, array2d &u_out){
  /* This function is used to add boundaries to an inner point array */
  /* u_in only has inner points, u_out has all the points */
	BoundaryCondition *pbval;
  int is=pmesh->is, ie=pmesh->ie, il=pmesh->il, iu=pmesh->iu;
  int js=pmesh->js, je=pmesh->je, jl=pmesh->jl, ju=pmesh->ju;
  int NX1=pmesh->NX1, NX2=pmesh->NX2;

	std::vector<double> x2_inner; // "top row" with only inner points
	std::vector<double> x2_outter; // "bottom row" with only inner points
  std::vector<double> row_coord = pmesh->coord_NX1; //each row has the same x coord
  std::vector<double> col_coord = pmesh->coord_NX2; //each col has the same y coord

  //add x2 boundaries
	for (int i=0; i<u_in.at(0).size(); i++){ //now only write inner points in x direction
    double xcoord = row_coord.at(is+i); //inner points
    x2_inner.push_back(pbval->bvalfunc_x2_inner(xcoord, jl, u_in.at(0).at(i)));
    x2_outter.push_back(pbval->bvalfunc_x2_outter(xcoord, ju, u_in.at(NX1-1).at(i)));
  }
  
  u_in.insert(u_in.begin(), x2_inner);
  u_in.push_back(x2_outter);

  for (int j=0; j<u_in.size(); j++){
    double ycoord = col_coord.at(j);
    u_in.at(j).insert(u_in.at(j).begin(), pbval->bvalfunc_x1_inner(il,ycoord,u_in.at(j).at(0)));
    u_in.at(j).push_back(pbval->bvalfunc_x1_outter(il,ycoord,u_in.at(j).at(NX2-1)));
  }
  u_out = u_in;
}


void BoundaryCondition::Apply2DBval(Mesh *pmesh, array2d &u_in){
  /* This function is used to apply B.C. to a full array */
  /* u_in, u_out both include boundary points */
  BoundaryCondition *pbval;
  int is=pmesh->is, ie=pmesh->ie, il=pmesh->il, iu=pmesh->iu;
  int js=pmesh->js, je=pmesh->je, jl=pmesh->jl, ju=pmesh->ju;

  std::vector<double> row_coord = pmesh->coord_NX1; //each row has the same x coord
  std::vector<double> col_coord = pmesh->coord_NX2; //each col has the same y coord

  //first, do the top and bottom rows
  for(int i=il; i<=iu; i++){
    double xcoord = row_coord.at(i);
    u_in.at(jl).at(i) = pbval->bvalfunc_x2_inner(xcoord, jl, u_in.at(js).at(i));
    u_in.at(ju).at(i) = pbval->bvalfunc_x2_outter(xcoord, jl, u_in.at(je).at(i));
  }

  //then, do the left and right cols
  for(int j=jl; j<=ju; j++){
    double ycoord = col_coord.at(j);
    u_in.at(j).at(il) = pbval->bvalfunc_x1_inner(il, ycoord, u_in.at(j).at(is));
    u_in.at(j).at(iu) = pbval->bvalfunc_x1_outter(iu, ycoord, u_in.at(j).at(ie));
  }


}


