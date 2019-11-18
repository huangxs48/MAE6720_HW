#include <cmath>
#include <iostream>
#include <vector>

#include "matrixsolver.h"

MatrixSolver::MatrixSolver(){
    
}

MatrixSolver::~MatrixSolver(){

}



#include <iostream>

using namespace std;



void MatrixSolver::TDMA(const std::vector<double>& a, 
                        const std::vector<double>& b,   
                        const std::vector<double>& c,
                        std::vector<double>& y,
                        std::vector<double>& x){
    size_t N = b.size();

    std::vector<double> c_(N, 0.0);
    std::vector<double> y_(N, 0.0);

    c_[0] = c[0] / b[0];
    y_[0] = y[0] / b[0];

    for (int i=1; i<N; i++){
        double m = 1.0/(b[i] - a[i]*c_[i-1]);
        c_[i] = c[i] * m;
        y_[i] = (y[i] - a[i]*y_[i-1]) * m;
    }

    for (int i=N; i-- > 0;){
        y[i] = y_[i] - c_[i]*y[i+1];
        x[i] = y[i];
        //printf("x[%d] = %g\n", i, x[i]);
    }

}