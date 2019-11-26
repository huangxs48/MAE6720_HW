#ifndef OUTPUT_H
#define OUTPUT_H

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#include "../Mesh/mesh.h"
#include "../PDE/PDE.h" //use gamma

typedef std::vector<std::vector<double>> array2d;

class Output{
  public:
    Output();
    ~Output();

    void Write1DArray(Mesh *pmesh, PartialDifferentialEquation *peq, int itr, std::vector<double> &coord, std::vector<double> &u, std::string username);
    void Write2DArray(Mesh *pmesh, PartialDifferentialEquation *peq, int itr, array2d &u, std::string username);


};



#endif