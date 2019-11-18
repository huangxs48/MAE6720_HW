#ifndef MESH_H
#define MESH_H

#include <cmath>
#include <iostream>
#include <vector>

class Mesh{
  public:
    Mesh();
    ~Mesh();

    //int dim; //dimension
    int NX1, NX2, NX3; //number of inner point in X,Y,Z direction
    double x1max, x1min, x2max, x2min, x3max, x3min; //
    int is, ie, js, je, ks, ke; //inner points index
    int il, iu, jl, ju, kl, ku; //index including boundaries
    double x1ratio, x2ratio, x3ratio; //later user? logarithm grid
    double dx1, dx2, dx3;

    double t_now;

    std::vector<double> x1v, x2v, x3v; //volume-centered coordinates in X,Y,Z direction

    //set up mesh indexs, dx
    void MeshInitialize();

    //return an array of coords
    //void UniformMeshConsructor(double min, double max, int N, std::vector<double> coord);
    std::vector<double> UniformMeshConstructor(double min, double max, int N);
    void LogMeshConstructor(std::vector<double> output, double min, double max, int N, double ratio);

  private:

    std::vector<double> x1f, x2f, x3f; 

};

#endif