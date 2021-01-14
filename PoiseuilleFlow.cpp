#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "VtkData.h"
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
    LBM(int nx, int ny, double dx, double dy, int q, int time, double tau, double rho0, double dpdx)
    : nx(nx), ny(ny), dx(dx), dy(dy), q(q), time(time), tau(tau), rho0(rho0), dpdx(dpdx)
    {
        this->wt.push_back(4.0/9.0);
        for(int i = 1; i < 5; i++)
            this->wt.push_back(1.0/9.0);
        for(int j = 5; j < 9; j++)
            this->wt.push_back(1.0/36.0);

        vector<vector<double>> f_temp;
        f_temp.assign(this->nx, vector<double>(this->ny, 0));
        this->f.assign(this->q, f_temp);

        vector<vector<double>> f_a_temp;
        f_a_temp.assign(this->nx, vector<double>(this->ny, 0));
        this->f_a.assign(this->q, f_a_temp);

        this->f_eq.assign(this->q, 0);

        this->ux.assign(this->nx, vector<double>(this->ny, 0.0));
        this->uy.assign(this->nx, vector<double>(this->ny, 0.0));
        this->rho.assign(this->nx, vector<double>(this->ny, 0.0));
        this->indicator.assign(this->nx, vector<double>(this->ny, 0.0));

        for(int i = 0; i < this->nx; i++)
        {
            for(int j = 0; j < this->ny; j++)
            {
                this->indicator[i][j] = 0; // 0 : fluid, 1 : solid
                if(j == 0 || j == ny - 1)
                    this->indicator[i][j] = 1;
                for(int alpha = 0; alpha < this->q; alpha++)
                    this->f[alpha][i][j] = this->wt[alpha]*rho0;
            }
        }

        this->ex.push_back(0); this->ey.push_back(0); // f0
        this->ex.push_back(1); this->ey.push_back(0); // f1
        this->ex.push_back(0); this->ey.push_back(1); // f2
        this->ex.push_back(-1); this->ey.push_back(0); // f3
        this->ex.push_back(0); this->ey.push_back(-1); // f4
        this->ex.push_back(1); this->ey.push_back(1); // f5
        this->ex.push_back(-1); this->ey.push_back(1); // f6
        this->ex.push_back(-1); this->ey.push_back(-1); // f7
        this->ex.push_back(1); this->ey.push_back(-1); // f8

        this->rev.assign(this->q, 0);

        for(int alpha = 0; alpha < this->q; alpha++)
        {
            for(int conv_a = alpha; conv_a < this->q; conv_a++) {
                if((this->ex[alpha] + this->ex[conv_a]) == 0 && (this->ey[alpha] + this->ey[conv_a]) == 0)
                {
                    this->rev[alpha] = conv_a;
                    this->rev[conv_a] = alpha;
                }
            }
        }

        source.assign(this->q, 0);
    }

    void LB_Method()
    {
        double rho_av; int kpor; double u2;
        double u_calpha; double u_calpha2;
        int i_a, j_a;
        for(int t = 1; t <= this->time; t++)
        {
            // Collision Step
            rho_av = 0.0; kpor = 0;
            for(int i = 0;i < this->nx; i++)
            {
                for(int j = 1; j < this->ny - 1; j++)
                {
                    // ux, uy, rho should always be initialized as 0 due to the accumulation of each time step.
                    this->ux[i][j] = 0.0; this->uy[i][j] = 0.0; this->rho[i][j] = 0.0;
                    for(int alpha = 0; alpha < this->q; alpha++)
                    {
                        this->rho[i][j] += this->f[alpha][i][j];
                        this->ux[i][j] += this->f[alpha][i][j]*ex[alpha];
                        this->uy[i][j] += this->f[alpha][i][j]*ey[alpha];
                    }
                    
                    this->ux[i][j] += this->dpdx/2.0;
                    this->ux[i][j] /= this->rho[i][j]; this->uy[i][j] /= this->rho[i][j];
                    u2 = this->ux[i][j]*this->ux[i][j] + this->uy[i][j]*this->uy[i][j];
                    for (int alpha = 0; alpha < this->q; alpha++)
                    {
                        u_calpha = this->ux[i][j]*this->ex[alpha] + this->uy[i][j]*this->ey[alpha];
                        u_calpha2 = u_calpha*u_calpha;
                        this->source[alpha] = (1.0-0.5/this->tau)*this->wt[alpha]*(3.0*(this->ex[alpha] - this->ux[i][j])
                                + 9.0*(this->ex[alpha]*this->ux[i][j] + this->ey[alpha]*this->uy[i][j])*this->ex[alpha])*this->dpdx;
                        this->f_eq[alpha] = this->wt[alpha]*this->rho[i][j]*(1.0 + 3.0*u_calpha + 4.5*u_calpha2 - 1.5*u2);
                        this->f_a[alpha][i][j] = this->f[alpha][i][j] - (this->f[alpha][i][j] - this->f_eq[alpha])/this->tau
                                + this->source[alpha]; //collision step
                    }
                }
            }

            // Streaming Step
            for(int i = 0;i < this->nx; i++)
            {
                for(int j = 1; j < this->ny - 1; j++)
                {
                    for(int alpha = 0; alpha < this->q; alpha++)
                    {
                        i_a = i + this->ex[alpha];
                        j_a = j + this->ey[alpha];
                        if(i_a < 0)
                            i_a = this->nx - 1;
                        if(i_a > this->nx - 1)
                            i_a = 0;
                        this->f[alpha][i_a][j_a] = this->f_a[alpha][i][j];                        
                    }
                }
            }

            // Bounceback Boundary Condition (No slip)
            for(int i = 0;i < this->nx; i++)
            {
                for(int j = 1; j < this->ny - 1; j++)
                {
                    if(this->indicator[i][j] == 0)
                    {
                        for(int alpha = 0; alpha < this->q; alpha++)
                        {
                            i_a = i + this->ex[alpha];
                            j_a = j + this->ey[alpha];
                            if(i_a < 0)
                                i_a = this->nx - 1;
                            if(i_a > this->nx - 1)
                                i_a = 0;

                            if(this->indicator[i_a][j_a] == 1)
                                this->f[alpha][i][j] = this->f[this->rev[alpha]][i_a][j_a];
                        }
                    }
                }
            }
        }
    }

    void print()
    {
        for(int j = 1; j < this->ny - 1; j++)
        {
            for(int i = 0;i < this->nx; i++)
            {
                cout << "[ " << this->ux[i][j] << ", " << this->uy[i][j] << " ]" << " ";                
            }
            cout << endl;
        }
        cout << endl;
    }
};

int main()
{
    LBM Poise(201, 21, 0.5, 0.5, 9, 20000, 0.8, 1.0, 1.0e-5);

    Poise.LB_Method();
    Poise.print();

    return 0;
}