#ifndef PDE_H
#define PDE_H

#include <cmath>
#include <iostream>
#include <vector>

typedef std::vector<std::vector<double>> array2d;

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


class SteadyHeatConduction2d : public PartialDifferentialEquation{
  public:
    SteadyHeatConduction2d();
    ~SteadyHeatConduction2d();

    array2d u_init;
    int Nx, Ny;
    double Residual;
    
};

class BurgersEquation : public PartialDifferentialEquation{
  public:
    BurgersEquation();
    ~BurgersEquation();

    double Gamma, C;
    double P_e, P_me; //Peclet number and mesh Peclet number

    double Residual;

};

#endif