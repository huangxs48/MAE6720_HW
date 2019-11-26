#ifndef MATRIXSOLVER_H
#define MATRIXSOLVER_H

#include <cmath>
#include <iostream>
#include <vector>

class MatrixSolver{ //probably add more matrix solver later

  public:
  	MatrixSolver();
  	~MatrixSolver();
  	
    void TDMA(const std::vector<double>& a, 
              const std::vector<double>& b,   
              const std::vector<double>& c,
              std::vector<double>& y,
              std::vector<double>& x);
    

};

#endif