#ifndef PDE_H
#define PDE_H

#include <cmath>
#include <iostream>
#include <vector>

class PartialDifferentialEquation{

  public:
    PartialDifferentialEquation();
    ~PartialDifferentialEquation();

    double gamma; //gamma that sets courant condition
    int N; //number of inner points

};

class HeatConduction : public PartialDifferentialEquation{
  public:
    HeatConduction();
    ~HeatConduction();

    //initial condition
    std::vector<double> u_init;
    //double InitCondition(double x);
    //void SetInitCondition(int xs, int xe);

    //vectors used in implicit method
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> x;
    std::vector<double> y;


};

#endif