#ifndef OUTPUT_H
#define OUTPUT_H

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#include "../Mesh/mesh.h"
#include "../PDE/PDE.h" //use gamma

class Output{
  public:
    Output();
    ~Output();

    void Write1DArray(Mesh *pmesh, PartialDifferentialEquation *peq, int itr, std::vector<double> &coord, std::vector<double> &u, std::string username);



};



#endif