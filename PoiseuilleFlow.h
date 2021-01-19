#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class LBM
{
private:
    int nx; // number of points in x
    int ny; // number of points in y
    double dx; // Grid spacing dx
    double dy; // Grid spacing dy
    int q; // Quantum number in DdQq model
    vector<vector<vector<double>>> f; // Distribution function f, pre-collision
    vector<vector<vector<double>>> f_a; // Distribution function f, after collision
    vector<double> f_eq; // Distribution function in equilibrium
    vector<double> wt; // Weight of distribution
    vector<int> ex; // x-direction of quantum
    vector<int> ey; // y-direction of quantum
    vector<int> rev; // Reverse direction of quantum
    vector<double> source; // Source term

    int time;
    double tau;
    double rho0;
    double dpdx;
    vector<vector<double>> ux;
    vector<vector<double>> uy;
    vector<vector<double>> rho;
    vector<vector<double>> indicator; // Indicator whether solid or fluid
    
public:
    LBM(int nx, int ny, double dx, double dy, int q, int time, double tau, double rho0, double dpdx);
    void LB_Method();
    void print();
};