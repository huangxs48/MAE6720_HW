#ifndef MESH_H
#define MESH_H

#include <cmath>
#include <iostream>
#include <vector>
#include <array>

typedef std::vector<std::vector<double>> array2d;
typedef std::vector<std::vector<std::array<double, 2>>> coord2d;

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

    //return an array of coords (TODO) rewrite 1D to void function
    std::vector<double> UniformMeshConstructor(double min, double max, int N);
  
    void Uniform2dMeshConstructor(double min_x1, double max_x1, double min_x2, double max_x2, int NX1, int NX2);
    coord2d MeshGrid; // to store the mesh grid
    std::vector<double> coord_NX1;
    std::vector<double> coord_NX2;

    void printcoord(coord2d coord);

  private:

    std::vector<double> x1f, x2f, x3f; 

};

#endif