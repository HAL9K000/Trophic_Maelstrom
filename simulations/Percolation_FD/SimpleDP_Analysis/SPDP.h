
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>

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

#ifndef PERCOL_DYNAMICS_H
#define PERCOL_DYNAMICS_H

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

void increase_stack_limit(long stack_size);
bool maxis(int a, int b);

void init_fullframe(double array[], int size);
void init_randframe(double array[], double mean, double sd, int size);
void init_constframe(double array[], double constant, int size);
void init_quarterframe(double array[], double c1, double c2, double c3, double c4, int L);
void init_solitarytear(double array[], int length);

double mean_of_array(double array[],int size);
double standard_deviation_of_array(double array[],int size);
double mean_of_vector(vector<double> array,int size);
double standard_deviation_of_vector(vector<double> array,int size);
double occupied_sites_of_vector(vector<double> array,int size);
auto meansq_spread_of_vector(vector<double> array, int g, int c_x, int c_y);

void var_mean_incremental_surv_runs(double t_avg_var_rep_N[][5], double X_curr[][2], int size);
void var_mean_incremental_all_runs(double rep_avg_var[][3], double X_curr[][2], int size, int r);
void spreading_incremental_surv_runs(double t_N_R2_rep[][4], double X_curr[][3], int size);
void var_mean_incremental(double rep_avg_var[][2], vector<vector <double>> &X_curr, int size, int r);
template<typename T> std::vector<double> linspace(T start_in, T end_in, int num_in);
template<typename T>std::vector<double> lnspace(T start_in, T end_in, int log_points);
std::vector<double> logarithmic_time_bins(double t_max, double dt);
std::vector<double> logarithm10_time_bins(double t_max, double dt);

//----------------------------- Stochastic SPDP RK4 Integration Machinery --------------------------------------------//

void fP_2D(double f[], vector<double> &Rho_t, int n[][2][5], double a, double b, double D, double t, double dt, double dx, int g);
void percNoDiff_2D(vector<vector<double>> &Rho, vector<double> &t_stop, double Rh0[],
	 double t_max, double a, double b, double D,  double dt, double dx, int r,  int g);

void eulerf_2D(double f[], vector<double> &Rho_t, int n[][2][5], double p, double D, double t, double dt, double dx, int g);
void determine_neighbours_R2(int neighbours_R2[][2][5], int g);
void determine_neighbours_S8(int neighbours_R2[][3][3], int g);
void find_neighbours_R2(double neighbours_R2[2][5], vector<double> &Rho_t, int p, int q, int g);
int twilight_zoneS8(vector<double> &Rho_t, int neighbours_S8[][3][3], int i);
int twilight_zoneR2(vector<double> &Rho_t, int neighbours_R2[][2][5], int i);
void f_2D(double f[], vector<double> &Rho_t, int n[][2][5], double p, double D, double t, double dt, double dx, int g);
//void f_2D(double f[], vector<double> &Rho_t, double p, double D, double t, double dt, double dx, int g);
void percolation_2D(vector<vector<double>> &Rho, vector<double> &t_stop, double Rh0[],
  double t_max, double p, double D, double sigma, double dt, double dx, int r,  int g);
//void stoc_dem_percolation(double Rho_0[], vector<double> &t_stop, int div, double t_max,
//	double p_start, double p_end, double D, double sigma, double dt, double dx, int r,  int g);
//----------------------------- Dornic For DP --------------------------------------------//
void fDo_2D(double f[], vector<double> &Rho_t, int n[][2][5], double b, double t, double dt, int g);
void percolationDornic_2D(vector<vector<double>> &Rho, double Rh0[], 
	 double t_max, double a, double b, double D, double sigma, double dt, double dx, int r,  int g);
void spreading_exp_Dornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, double Rh0[],
	 double t_max, double a, double b, double D, double sigma, double dt, double dx, int r, int g , int c_x, int c_y);
void percolationHalfDornic_2D(vector<vector<double>> &Rho, vector<double> &t_stop, double Rh0[],
	 double t_max, double p, double D, double sigma, double dt, double dx, int r,  int g);


//----------------------------- Crt Exp Determination Machinery --------------------------------------------//
void critical_exp_delta(double Rho_0[], int div, double t_max, double a_start, double a_end, 
double b, double D, double sigma, double dt, double dx, int r,  int g);
void critical_exp_finite_scaling_collapse(int div, int g_st, int g_end, double t_max,  double a_c, 
	double b, double D, double sigma, double dt, double dx, int r);
void critical_exp_finite_scaling_stationary(int div, int g_st, int g_end, double t_lim,  double t_max, double a_c, 
	double b, double D, double sigma, double dt, double dx, int r);
void critical_exp_theta_delta(double Rho_0[], int div, double t_max, double a_c, double b, 
		double D, double sigma, double dt, double dx, int r,  int g);

//----------------------------- Dornic For Expanded DP --------------------------------------------//
void f_exp2D(vector<double> &f, vector<double> &Rho_t, double a, double b, double c, double D, double t, double dt, double dx, int g);
void RK4_Integrate(vector<double> &Rho_t, double a, double b, double c, double D, double t, double dt, double dx, int g);
void expanded_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, double Rh0[],
	 double t_max, double a, double b, double c, double D, double sigma, double dt, double dx, int r,  int g);

void first_order_critical_exp_delta(double Rho_0[], int div, double t_max, double a_start, double a_end, 
	double b, double c, double D, double sigma, double dt, double dx, int r,  int g);

void terminal_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, double Rh0[],
	 double t_max, double a, double b, double c, double D, double sigma, double dt, double dx, int r,  int g);
void first_order_terminal_diff(int div, double t_event, double p_start, double p_end, double a,
	double b, double c, double D, double sigma, double dt, double dx, int r,  int g);


void frames_expanded_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, double Rh0[],
	 double t_max, double t_event, double a, double b, double c, double D, double sigma, double dt, double dx, int g, int r);
void capture_frame_decay(double Rho_0[], int div, double t_max, double t_event, double a_start, double a_end, 
	double b, double c, double D, double sigma, double dt, double dx, int r,  int g);

void test_gamma_poisson();
#endif



/** void f_2D(double f[], vector<double> &Rho_t, double p, double D, double t, double dt, double dx, int g)
{
  //Vector function that evaluates drho/dt = f(p,t) at given p and t.

	for(int i=0; i < g*g; i++)
	{
		int a = int(i/g); int b = i%g; //Getting 2D coordinates.
		double n_R2[2][5] ={0.0}; //Stores neighbours in a ball of radius 2 around a,b.
		//First row stores horizontal neighbours, second row stores vertical neighbours.

		find_neighbours_R2(n_R2, Rho_t, a, b, g); //Finds these neighbours.

		f[i] = -Rho_t[i]*Rho_t[i] - D*(1/(12.0*dx*dx))*(1*(n_R2[0][0] + n_R2[0][4] + n_R2[1][0] + n_R2[1][4])
		-16*(n_R2[0][1] + n_R2[0][3] + n_R2[1][1] + n_R2[1][3]) +60*n_R2[0][2]);
		// RK4 update of Non-linear parts plus diffusion term
		// 4th order Laplacian from 2SHOC.
	}
} */

/** void stoc_dem_percolation(double Rho_0[], vector<double> &t_stop, int div, double t_max,
	double p_start, double p_end, double D, double sigma, double dt, double dx, int r,  int g)
{
	vector<double> p_space = linspace(p_start, p_end, div);
  // The pspace to iterate over.
	t_stop = linspace(t_max/3.0, t_max, 2); // 2 Time-stamps at which screenshots will be captured.
	for( int i =0 ; i< t_stop.size(); i++)
		{ cout << t_stop[i] <<endl; }

	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

	auto start = high_resolution_clock::now();

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */
	  /**
      #pragma omp for nowait schedule(static)
      for (int i=0; i < p_space.size(); i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on p Value:\t" << p_space[i] <<endl;
        cout << message.str();
        std::vector<vector <double>> Rho_p;
		double Rh0[g*g]={0.0}; //Stores initial density conditions
		for(int j =0; j < g*g; j++)
		{	Rh0[j] = Rho_0[j]; }
		percolation_2D(Rho_p, t_stop, Rh0, t_max, p_space[i], D, sigma, dt, dx, r, g);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), Rho_p.begin(), Rho_p.end());

      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting critical exponent data for each grid size in order.
          stringstream message3;
          message3 << "Is this happening?\n";
          cout << message3.str();

      }
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "RK4 Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, tm ,p1, p2, rini, Dm, Sm;

  L << g; tm << t_max; p1 << setprecision(4) << p_start; p2 << setprecision(4) << p_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_dp;
  // Creating a file instance called output to store output data as CSV.
	output_dp.open("../Data/SPDP/P_c_DP_G_" + L.str() + "_T_" + tm.str() + "_p1_"+ p1.str() +
	"_p2_"+ p2.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	// Output = | p (Order Parameter) | Replicate (r) | Time(t)  | a*L + b | Cell Density[a][b] |
	output_dp << " p , # Tr No , t ,  i*L + j ,  Rho_ij(t) \n";
	cout << "The vector elements are: "<< endl;
  cout << " p , # Tr No , t ,  i*L + j ,  Rho_ij(t) \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_dp << setprecision(8) << vec[i][0] << "," << static_cast<int>(vec[i][1]) << "," << setprecision(7)
		<< vec[i][2] << "," << static_cast<int>(vec[i][3]) << "," << setprecision(16) << vec[i][4] << endl;
		if( i%(5000) ==1)
    {
			cout << setprecision(5) << vec[i][0] << "," << static_cast<int>(vec[i][1]) << "," << setprecision(7)
			<< vec[i][2] << "," << vec[i][3] << "," << setprecision(16) << vec[i][4] << endl;
    }
	}
	output_dp.close();


} */


///-----------------------CODE SNIPPET FOR DEBUGGING DORNIC----------------------------------------//


/**stringstream m10;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m10 << "In Replicate:\t" << j << "\t For Time t: " <<t << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< "Rhox_AVG CALCULATION SUCCESSFUL" << endl; cout << m10.str();
			
stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "CHECKERED PAST AND STOC UPDATE FOR Replicate:\t" << j << "\t At Time t: " << t << "\t for Thread Rank:\t " 
			<< omp_get_thread_num() <<endl; cout << m7.str();**/

///FIND NAN PROPAGATION THROUGHOUT GRID.
/**vector <double> Err_Rho;
			if(isnan(rhox_avg) && s == 1)
			{	
				stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m0 << "In Replicate:\t" << j << " AT t:\t" << t << " <Rho(t)>x:\t " << rhox_avg <<  
				" TURNS UNSTABLE!!!!: \n" << endl; 
				cout << m0.str();
				//Rhox_avg turns nan for the first time.
				//Now identify culprit.
				for(int x=0; x <g*g; x++)
				{
					if(isnan(Rho_dt[x]))
					{	
						Err_Rho.push_back(x);
							//cout << "Located at index: \t" << x << " | " << static_cast<int>(x/g) << " , " << x%g <<endl;
						for(int y=0; y < 5; y++)
						{
							if(isnan(Rho_dt[nR2[x][0][y]]) || isnan(Rho_dt[nR2[x][1][y]]))
							{
								//cout << "BAD NEIGHBOUR AS WELL!!!" <<endl;
							}
						}
					}
					
				}
				cout << "THE CULPRITS ARE: " << Err_Rho.size() << "  " << g*g <<endl;
				if(po == 1)
				{
					int c = Err_Rho.size();
					if(c>20)
					{	c=20;}
					for(int z=0; z<c; z++)
					{	cout<< Err_Rho[z] << " "; } cout<< "THAT'S ALL!" << endl; po =-1;
					for(int z=c; z>0; z--)
					{	cout<< Err_Rho[Err_Rho.size()-z] << " "; } cout<<endl;
				}
				if(Err_Rho.size() == g*g)
				{
					cout << "At TIME t:\t" << t <<  "ALL CELLS ARE NAN. \n" <<endl; s =-1;
				}
			}**/


// CHECK FOR NANS

			/**if(DRho[i] == 0.0)
			{
				//First NAN encountered.
				cout << "DRHO IS ZERO!"<<endl;
				stringstream m2;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m2 << "In Replicate:\t" << j <<  " at t:\t" << t << " in INDEX:\t" << i 
				<<  " we have 0 \n"  << endl; 
				cout << m2.str(); co=-1;
			}
			if(isnan(Rho_dt[i]) && eo ==1)
			{
				//First NAN encountered.
				cout << "FIRST NAN IN RHO_dt!"<<endl;
				stringstream m2;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m2 << "In Replicate:\t" << j << "\t for <Rho(t)>x:\t " << rhox_avg << " at t:\t" << t << " in INDEX:\t" << i << " we have: \n" 
				<< "Rho_dt[i] (after Stoc Update): " << Rho_dt[i] << " with Gamma: " << gamma(rng) << " & Poisson: " << poisson(rng)
				<< "\n with Alpha[i]: " << alpha_i << " and mu: " << mu << " and lambda: " << lambda << endl; 
				cout << m2.str(); eo=-1;
			}
			if(isnan(alpha_i) && so==1)
			{
				cout << "EXAMINING which of the terms of alpha_i are at fault:\t" <<endl;
				cout << "INDEX:\t" << i << "AT TIME:\t" << t << "\n" << endl;
				cout << "C: " << DRho[nR2[i][0][2]] << " " << DRho[nR2[i][1][2]] << endl;
				cout << "T1: " << DRho[nR2[i][0][1]] << " T2: " << DRho[nR2[i][0][0]] << endl;
				cout << "B1: " << DRho[nR2[i][0][3]] << " B2: " << DRho[nR2[i][0][4]] << endl;
				cout << "L1: " << DRho[nR2[i][1][1]] << " L2: " << DRho[nR2[i][1][0]] << endl;
				cout << "R1: " << DRho[nR2[i][1][3]] << " R2: " << DRho[nR2[i][1][4]] << "\n\n" << endl;
				so =-1;
				
			}
	    	}**/

			/**
			 * if( rho_rep_avg_var[i][0] < 26 && isnan(rho_rep_avg_var[i][1]))
				{	stringstream m2;     //To make cout thread-safe as well as non-garbled due to race conditions.
					m2 << "In Replicate:\t" << j << "\t for AVG <<Rho(t)>x>r:\t " << rho_rep_avg_var[i][0] << "at t:\t" << rho_rep_avg_var[i][0] << " YOU ARE NICKED: \n" 
					<< endl; cout << m2.str();	
				}


				for(int i =0; i< tot_iter; i++)
				{
				if( rho_rep_avg_var[i][0] < 26 && isnan(rho_rep_avg_var[i][1]))
				{	stringstream m2;     //To make cout thread-safe as well as non-garbled due to race conditions.
					m2 << "In Replicate:\t" << j << "\t for AVG <<Rho(t)>x>r:\t " << rho_rep_avg_var[i][1] << "at t:\t" << rho_rep_avg_var[i][0] << " YOU ARE NICKED: \n" 
					<< endl; cout << m2.str();	
				}
			}
			 * 
			 */


//-----------------------RK4 UPDATE FOR DORNIC INTEGRATION ---------------------------------------------------------//


/** fDo_2D(K1, Rho_dt, nR2, b, t, dt, g); //K1 updated. 
			cout<< "K1" << endl;
			for(int i = 0; i<g*g; dx, i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K1[i];  } //Recall Runge Kutta-Update Schema for K2 = f(t+dt/2, y+dt/2*K1)
			cout << endl;
			fDo_2D(K2, Rho_M, nR2, b, t +dt/2.0, dt, g); //K2 updated.
			cout<< "K2" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K2[i]; } //Recall Runge Kutta-Update Schema for K3
		  	cout << endl;
			fDo_2D(K3, Rho_M, nR2, b, t +dt/2.0, dt, g); //K3 updated.
			cout<< "K3" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt)*K3[i]; } //Recall Runge Kutta-Update Schema for K4
			cout << endl;
			fDo_2D(K4, Rho_M, nR2, b, t + dt, dt, g); //K4 updated.

			for(int i=0;i<Rho_dt.size();i++)
			{ Rho_dt[i]+= (dt/6.0)*( K1[i] + 2*K2[i] + 2*K3[i] + K4[i]); } **/
