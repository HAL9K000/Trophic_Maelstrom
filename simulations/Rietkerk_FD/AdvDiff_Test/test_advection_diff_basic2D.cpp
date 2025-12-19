#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <cstdlib>
#include <iomanip>
#include <sstream>
using namespace std;


const double Dx = 0.5; // Grid spacing
const int L = 256; // Domain size (assuming uniform grid spacing in x and y)
const int Nx = static_cast<int>(L / Dx) + 1; // Number of grid points in x and y directions

const double Dt = 0.25; // Time step
const double T_end = 100.0; // Simulation end time
const double vx = 1.0;    // Advection velocity in both x and y directions
const double vy = 1.0;    // Advection velocity in both x and y directions
const double v =1.0;
const double D = 0.2;   // Diffusion coefficient

const double a1 = 1; // Determines rate of advection (as opposed to diffusion)

const int numSteps = static_cast<int>(T_end / Dt);

// Function to initialize the concentration field
void initializeConcentration(std::vector<std::vector<double>>& C) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Nx; ++j) {
            double x = i * Dx;
            double y = j * Dx;
            C[i][j] = std::exp(-((x - 128.0)*(x - 128.0) + (y - 128.0)*(y - 128.0)) / 32.0);
        }
    }
}

// Function to apply periodic boundary conditions
void applyPeriodicBoundaryConditions(std::vector<std::vector<double>>& C) {
    // Periodic boundary conditions in x
    for (int j = 0; j < Nx; ++j) {
        C[0][j] = C[Nx - 2][j];
        C[Nx - 1][j] = C[1][j];
    }

    // Periodic boundary conditions in y
    for (int i = 0; i < Nx; ++i) {
        C[i][0] = C[i][Nx - 2];
        C[i][Nx - 1] = C[i][1];
    }
}

// Function to save the concentration field
void saveFrame(const std::vector<std::vector<double>>& C, double t) 
{
    //std::string fileName = "./frame_T_" + std::to_string(t) + "_adv_rate_" + std::to_string(a1) + ".csv";

    //Suggest a different filename where all trailing zeros in t and a1 are removed.
    stringstream tm, a, va, vb;
    tm << t; a << a1; va << vx; vb << vy;
    string fileName = "./frame_T_" + tm.str() + "_adv_rate_" + a.str() + "_vx_" + va.str() + "_vy_" +vb.str() + ".csv";
    // Save frame as 2D (256X256) CSV file.
    std::ofstream file(fileName);
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Nx; ++j) {
            file << C[i][j] << ",";
        }
        file << "\n";
    }
    file.close();
}

// Function to print the concentration field
void printConcentration(const std::vector<std::vector<double>>& C) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Nx; ++j) {
            std::cout << C[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// Function to simulate the advection-diffusion equation
void simulate(std::vector<std::vector<double>>& C) 
{
    std::vector<std::vector<double>> C_Copy(Nx, std::vector<double>(Nx, 0.0));
    for (int step = 0; step < numSteps; ++step) 
    {
        // Apply periodic boundary conditions
        applyPeriodicBoundaryConditions(C);
        for(int i = 0; i < Nx; ++i)
        {
            for(int j = 0; j < Nx; ++j)
            {
                C_Copy[i][j] = C[i][j];
            }
        }
        // Advection term
        if(step == 0)
        {
            saveFrame(C_Copy, 0);
        }
        for (int i = 1; i < Nx - 1; ++i) 
        {
            
            for (int j = 1; j < Nx - 1; ++j) 
            {
                //C[i][j] = C_Copy[i][j] + Dt*( - a1*vx/ Dx * (C_Copy[i][j] - C_Copy[i - 1][j]) - a1*vy/ Dx * (C_Copy[i][j+1] - C_Copy[i][j])
                C[i][j] = C_Copy[i][j] + Dt*( - a1*v/ Dx * (2*C_Copy[i][j] - C_Copy[i - 1][j] - C_Copy[i][j-1])
                            + (1 - a1)*D /(Dx * Dx) * (C_Copy[i + 1][j] + C_Copy[i][j + 1] + C_Copy[i][j - 1] 
                            + C_Copy[i - 1][j] - 4 * C_Copy[i][j] ));
            }
        }

        if ((step + 1) % static_cast<int>((T_end / 10) / Dt) == 0) {
            std::cout << "Saving frame " << (step + 1)*Dt << endl;
            saveFrame(C_Copy, (step + 1)*Dt);
        }
    }
}

int main() {
    // Initialize concentration field
    std::vector<std::vector<double>> C(Nx, std::vector<double>(Nx, 0.0));
    initializeConcentration(C);

    // Simulate the advection-diffusion equation
    simulate(C);

    // Print the final concentration field
    printConcentration(C);

    return 0;
}