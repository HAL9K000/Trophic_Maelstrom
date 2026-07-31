
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <list>
#include <thread>
using namespace std::chrono;
//static std::mt19937_64 rng(time(NULL));
#include<bits/stdc++.h>
using namespace std;
#include <sys/resource.h>
#include <omp.h>

#include <map>

#ifndef BERSERK_DYNAMICS_H

#define BERSERK_DYNAMICS_H



inline const int Sp = 3; //Number of species in the system.
inline const int Sp4_1 = 4*(Sp -2) +3; // Used for generating statistics on surviving runs
inline const int Sp2 = 2*Sp; // Used for generating statistics on surviving runs


//------------------------- TYPEDEFS --------------------------------------------------------------//

typedef std::vector<std::vector <double>> D2Vec_Double;
typedef std::vector<std::vector <int>> D2Vec_Int;

typedef std::vector<std::vector <std::vector <double>>> D3Vec_Double;
typedef std::vector<std::vector <std::vector <int>>> D3Vec_Int;

typedef std::vector<std::vector <std::vector <std::vector <double>>>> D4Vec_Double;
typedef std::vector<std::vector <std::vector <std::vector <int>>>> D4Vec_Int;



//---------------------Some Primitives ------------------------------------------------------ //



struct coordinates {
   int x;
   int y;
};

struct pde_val {
	int label;  //Stores grid location of PDE values
	std::vector<double> pde;   //PDE values at label index.
};

struct f_coordinates {
   float x;
   float y;
};

struct zd_coordinates {
   double x;
   double y;
   double z;
};



void increase_stack_limit(long long stack_size);
bool maxis(int a, int b);
void add_three(int a, int b, int c); //Test function.

void init_fullframe(D2Vec_Double &array, int Sp, int size);
void init_randconstframe(D2Vec_Double &array, int Sp, int size, double perc,   double c_high[], double c_low[]);
void init_constframe(D2Vec_Double &array,  int Sp, int size, double constant[]);
void init_randframe(D2Vec_Double &array, int Sp, int size, double mean[], double sd[]);
void init_quarterframe(D2Vec_Double &array, double c1, double c2, double c3, double c4, int L);
void init_solitarytear(D2Vec_Double &array, int length);

double mean_of_array(double array[],int size);
double standard_deviation_of_array(double array[],int size);
double mean_of_vector(vector<double> array,int size);
double standard_deviation_of_vector(vector<double> array,int size);
double occupied_sites_of_vector(vector<double> array,int size);
auto meansq_spread_of_vector(vector<double> array, int g, int c_x, int c_y);

void var_mean_incremental_surv_runs(double t_avg_var_rep_N[][Sp4_1], double X_curr[][Sp2], int size);

void var_mean_incremental_all_runs(double rep_avg_var[][3], double X_curr[][2], int size, int r);
void spreading_incremental_surv_runs(double t_N_R2_rep[][4], double X_curr[][3], int size);
void var_mean_incremental(double rep_avg_var[][2], vector<vector <double>> &X_curr, int size, int r);

// --------------------------- Partitioning functions to make life easier -----------------------------------//

template<typename T> std::vector<double> linspace(T start_in, T end_in, int num_in);
template<typename T>std::vector<double> lnspace(T start_in, T end_in, int log_points);
std::vector<double> logarithmic_time_bins(double t_max, double dt);
std::vector<double> logarithm10_time_bins(double t_max, double dt);

//----------------------------- Misc. Supporting Add-ons -------------------------------------------------

int theborderwall(D2Vec_Double &Rho_t, int g);
void determine_neighbours_R2( int g, int S, D3Vec_Int &neighbours_R2);
void determine_neighbours_Sq8(int g, int S, D3Vec_Int &neighbours_Sq8);

//----------------------------- Stochastic 2PAC Dornic Integration Machinery --------------------------------------------//

void f_2D(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, double alpha, double d, double rW, double W0,  
	double D[], double K[],  double t, double dt, double dx, int g);
void RK4_Integrate(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, double a, double b, double D[], double A[Sp][Sp], 
	double H[Sp][Sp], double E[], double M[], double t, double dt, double dx, int g);
void RK4_Wrapper_2D(D2Vec_Double &Rho, vector <double> &t_meas, D2Vec_Double &Rh0, double t_max, double a, double c, double gmax, 
	  double alpha, double d, double rW, double W0, double D[], double K[], double a_st, double a_end,  double dt, double dx, double dP, int r, int g);

void first_order_critical_exp_delta(D2Vec_Double &Rh0, int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double d, double rW, double W0, double D[], double K[], double dt, double dx, double dP, int r,  int g);

#endif