//#pragma once
#include <string>
//#include <format>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

#if defined(BARRACUDA) || defined(__CUDACC__)
#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftXt.h>
#endif


#include <regex>
#include <limits>
#include <random>
#include <vector>
#include <algorithm>
#include <ranges>
#include <span>
#include <chrono>
#include <cmath>
#include <list>
#include <thread>
#include <mutex>
/** NOTE: The following include is for the exprtk library. Download and add to your project/include path.
// if not already present. The library is used for parsing mathematical expressions. */
#include <exprtk.hpp>
/** NOTE: The following include is for the FFTW3 library. Download, configure, make, install and add to your project/include path.
// if not already present. The library is used for performing discrete FFTs. */
#include <fftw3.h>

using namespace std::chrono;
//static std::mt19937_64 rng(time(NULL));
#include<bits/stdc++.h>
using namespace std;
#include <sys/resource.h>
#include <omp.h>

#include <map>

#if defined(__GNUC__) && (__GNUC__ >= 9)
#include <filesystem>
namespace fs = std::filesystem;
#endif

#ifndef BERSERK_DYNAMICS_H

#define BERSERK_DYNAMICS_H

#ifndef SPB
	#define SPB 3 //Number of biota species in the system (default value used if NOT defined).
#endif

#ifndef INIT
	#define INIT 0 //Initial condition for the system (0 = Homogeneous, 1 = Random Heterogenous Speckles, 2 = Burn-in)
#endif

#if SPB == 3
	#include "MultiSPDP_constants_3Sp.h"
#elif SPB == 2
	#include "MultiSPDP_constants_2Sp.h"
#elif SPB == 1
	#include "MultiSPDP_constants_1Sp.h"
#else
	#error "Number of species not supported. Pass valid compiler flag as -DSPB=2 or -DSPB=3"
#endif

//------------------------- STRUCTS & TYPEDEFS --------------------------------------------------------------//

// Plain POD for device memory
struct DoubleIntPair {
    double first;
    int    second;
};

typedef std::vector<std::vector <double>> D2Vec_Double;
typedef std::vector<std::vector <int>> D2Vec_Int;

typedef std::vector<std::vector <std::vector <double>>> D3Vec_Double;
typedef std::vector<std::vector <std::vector <int>>> D3Vec_Int;

typedef std::vector<std::vector <std::vector <std::vector <double>>>> D4Vec_Double;
typedef std::vector<std::vector <std::vector <std::vector <int>>>> D4Vec_Int;

typedef std::vector<std::vector<std::pair<double, double>>> D2Vec_Pair_Double;
typedef std::vector<std::vector<std::pair<int, int>>> D2Vec_Pair_Int;

typedef std::vector<double> Vec_Double;
typedef std::vector<int> Vec_Int;

//------------------------- GLOBAL CONSTANTS --------------------------------------------------------------//

// CUDA constants (actually defined in the .cu kernel files)
#if defined(__CUDACC__) || defined(BARRACUDA)
#ifdef DEFINE_CUDA_CONSTANTS
__constant__ double CuA[CuSpNV] /** ={aij, ajm} */;
__constant__ double CuH[CuSpNV] /** ={hij, hjm} */;
__constant__ double CuE[CuSpNV] /** ={ej, em} */;
__constant__ double CuM[CuSpNV] /** ={mj, mm} */;
__constant__ double CuD[Sp] /** ={D0, D1, D2, D3, D4} */;
__constant__ double CuK[3] /** ={K0, K1, K2} */;
__constant__ double CuV[CuSpNV] /** ={v0, v1} */;
__constant__ DoubleIntPair d_r_frac[CuSpB]; /** ={ {fr_j, idx_j}, {fr_m, idx_m} } **/;
#endif
#endif


inline const int SpB = Sp; //Number of biota species in the system.
inline const int Sp_NV = Sp-1; //Number of grazer and predator species in system.
inline const int Sp4_1 = 4*(Sp) +1; // Used for generating statistics on surviving runs
inline const int Sp2 = 2*Sp; // Used for generating statistics on surviving runs
inline const int SpB2 = 2*SpB; // Used for generating statistics on surviving runs
inline const int SpB3 = 3*SpB; // Used for generating statistics on surviving runs
inline const int SpB4 = 4*SpB; // Used for generating statistics on surviving runs

// Global user-defined parameters for the Rietkerk model.
inline double dt2 /** = dt/2.0*/; inline double dt6 /** = dt/6.0*/;
inline double dx2 /** = dx*dx*/; inline double dx1_2 /** = 1.0/(dx*dx)*/; inline double dxb2 /** = dx/2.0/ */; inline double dx1 /** =  1/dx */;

// Global user-defined parameters for the predator-prey trophic chain.
inline double K_P; /** Carrying capacity of predator */
inline double K_P1; /** Inverse of carrying capacity of predator */

// Rietkerk Model Parameters, set in set_global_user_Rietkerk_params() function.
//inline double cgmax /** = c*gmax */; inline double K2W0 /**= K[2]*W0 */;inline double A01H01 /** = A[0][1]*H[0][1] */;inline double A12H12 /**= A[1][2]*H[1][2] **/;

// Strings that are used for I/O operations and file naming in the simulations.
inline string prefix = "DiC-NREF-NI";
inline string frame_folder = "../Data/DP/Frames/Stochastic/"+ std::to_string(SpB) +"Sp/" + prefix + "_"; //Folder to store frames.
inline string prelim_folder = "../Data/DP/Prelims/Stochastic/"+ std::to_string(SpB) +"Sp/"+ prefix +"_"; //Folder to store preliminary data.
inline const string frame_prefix = "/FRAME_P_c_DP_G_"; //Prefix for frame files.
inline const string gamma_prefix = "/GAMMA_G_"; //Prefix for gamma files.
inline const string prelim_prefix = "/PRELIM_AGGRAND_P_c_ID_"; //Prefix for preliminary data files.
inline const string replicate_prefix = "/PRELIM_TSERIES_P_c_DP_G_"; //Prefix for replicate time-series data files.
inline const string movement_prefix = "/PRELIM_MOVSERIES_P_c_DP_G_"; //Prefix for movement time-series data files.

inline const string input_prefix = "/FRAME_T_"; //Prefix for input files.
inline string stat_prefix = "../Data/DP/Stochastic/"+ std::to_string(SpB) +"Sp/1stCC_"  + prefix + "_P_c_G_";

inline string input_frame_parenfolder = "../Input/DP/"+ std::to_string(SpB) +"Sp"; //Parendirectory for input data.
inline string input_frame_subfolder =""; //Subfolder for input data.
// NOTE: Complete path to input data is input_frame_parenfolder + input_frame_subfolder + input_prefix +"*.csv"
inline std::map<string, string> input_keys = { {"outPRE", ""}, {"a", "0"}, {"a_c", "0"}, {"ascaled", "0"}, {"T", "0"}, {"dP", "0"}, {"Geq", "0"}, {"R", ""}, 
{"MAMP", "0"}, {"SEP", "0"}, {"WID", "0"}, {"MSF", "0"}};

//inline const string prelimheader = " a , r, L, t , <<P(x; t)>_x>_r, Var[<P(x; t)>_x]_r, # Surviving Runs P(x; t),"
//	" # Active Sites P(x; t), <<G(x; t)>_x>_r, Var[<G(x; t)>_x]_r, # Surviving Runs G(x; t), # Active Sites G(x; t), \n";

// Strings that are used for the Expertk library to parse mathematical expressions.
// These represent MFT-versions of the species Coexistance equilibria wrt to the order parameter (Rainfall, given by a).

inline vector<string> MFT_Vec_CoexExpr(2*Sp, ""); //Vector of MFT Coexistance expressions for all species.

//inline exprtk::symbol_table<double> global_symbol_table; //Symbol table for the Expertk library.
inline vector<double> frame_tmeas; //Vector to store time measurements for the frames.

// Cache to store FFTW3 K-Hat kernels for fixed FFT transforms of central neighbourhoods for reuse.
inline std::map<int, fftw_complex*> kHat_kernel_FFT_cache;
// Map that stores neighbourhood counts for different perception ranges R as {R, effective neighbour count}
inline std::map<int, int> nR_Perp_R; 


//---------------------Some Primitives ------------------------------------------------------ //

//Stolen from: https://web.archive.org/web/20200629174939/http://timmurphy.org/2013/06/27/template-metaprogramming-in-c/
// Generic exponent compile-time calculations. Pow<2,3>::result == 2^3
template < long long B, long long E>
struct Pow
{
    static const long long result = B * Pow<B, E-1>::result;
};

template <long long B>
struct Pow<B, 0>
{
    static const long long result = 1;
};

// Bailey-Borwein-Plouffe formula for calculating pi
// http://en.wikipedia.org/wiki/Bailey-Borwein-Plouffe_formula
template <long long int N>
struct CalculatePi
{
    static const constexpr double pi =
    (
        1.0/Pow<16,N>::result *
        (
            4.0/(8*N + 1.0) - 2.0/(8*N + 4.0) -
            1.0/(8*N + 5.0) - 1.0/(8*N + 6.0)
        )
    ) + CalculatePi<N-1>::pi;
};

template <>
struct CalculatePi<-1>
{
    static const constexpr  double pi = 0;
};

// Generic Permutation compile-time calculations. Perm<N,K>::result == N!/(N-K)!
template <long long N, long long K>
struct Perm
{
	static const long long result = N * Perm<N-1, K-1>::result;
};

template <long long N>
struct Perm<N, 0>
{
	static const long long result = 1;
};

// Generic Combination compile-time calculations. Comb<N,K>::result == N!/(K!(N-K)!)
template <long long N, long long K>
struct Comb
{
	static const long long result = Perm<N,K>::result / Perm<K, K>::result;
};

template <long long N>
struct Comb<N, 0>
{
	static const long long result = 1;
}; 

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

// Struct for FFTW3 plan management
struct FFTW3_CentralPlanner {
    int L;
    fftw_complex* rho_fft = nullptr; // Stores FFT of input (density) field
    fftw_complex* kernel_fft = nullptr; // Stores FFT of kernel (neighbourhood) input field.
    fftw_complex*  prod_fft = nullptr; // Stores products of frequencies in Fourier space
    double* rho_spatial; // Stores input (density) field in spatial domain
    double* kernel_spatial; // Stores neighourhodd kernel in spatial domain
    double* convolution_spatial; // Stores output (convolution result) field in spatial domain
    fftw_plan forward_rho  = nullptr; // Plan for forward FFT of density field
    fftw_plan forward_kernel  = nullptr; // Plan for forward FFT of kernel field
    fftw_plan backward  = nullptr;
    
	FFTW3_CentralPlanner(int L_) : L(L_) {
        int L2 = L * L;
        rho_fft = fftw_alloc_complex(L2);
        kernel_fft = fftw_alloc_complex(L2);
        prod_fft = fftw_alloc_complex(L2);
        rho_spatial = fftw_alloc_real(L2);
        kernel_spatial = fftw_alloc_real(L2);
        convolution_spatial = fftw_alloc_real(L2);
        forward_rho = fftw_plan_dft_r2c_2d(L, L, rho_spatial, rho_fft, FFTW_ESTIMATE);
        //forward_kernel = fftw_plan_dft_r2c_2d(L, L, kernel_spatial, kernel_fft, FFTW_ESTIMATE);
        backward = fftw_plan_dft_c2r_2d(L, L, prod_fft, convolution_spatial, FFTW_ESTIMATE);
    }
    
    ~FFTW3_CentralPlanner() {
        if (forward_rho) fftw_destroy_plan(forward_rho);
        if (forward_kernel) fftw_destroy_plan(forward_kernel);
        if (backward) fftw_destroy_plan(backward);
        if (rho_fft) fftw_free(rho_fft);
        if (kernel_fft) fftw_free(kernel_fft);
        if (prod_fft) fftw_free(prod_fft);
        fftw_free(rho_spatial);
        fftw_free(kernel_spatial);
        fftw_free(convolution_spatial);
    }
};

inline const double PI = CalculatePi<14>::pi;
inline const int SpB_Mov = 2*SpB + Comb<SpB, 2>::result; // Used for generating statistics on movement of surviving runs

//void set_global_user_Rietkerk_params(double c,double gmax,double alpha, double rW, double W0, double (&D)[Sp], 
//	double (&K)[3], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB])

template <typename T> int sgn(T val);
void increase_stack_limit(long long stack_size);
bool maxis(int a, int b);
string format_str(const std::string& s, const std::map<string, string>& values);
void precompute_FFT_KHat_kernels(FFTW3_CentralPlanner& fft_centralplanner, int L, int Rmax, double nVeg_frac_min = 0.2, double nVeg_frac_max = 1.0);
void add_three(int a, int b, int c); //Test function.
void set_Prefix(string& user_prefix, double mG = 0, double mP = 0);
void set_LocustPrefix(string& user_prefix);
void set_input_Prefix(string& user_inputframepath, string& user_prefix, double user_a_c, double user_dP, double user_Geq = -1, double user_input_T = -1);
void set_global_system_params(double dt, double dx);
void set_global_predator_params(double Km);

void display_symbol_table(const exprtk::symbol_table<double>& symbol_table);
// Custom comparator for sorting pairs based on the squared distance
bool comparePairs(const pair<int, int>& a, const pair<int, int>& b); 

void init_fullframe(D2Vec_Double &array, int Sp, int size);

void init_randbistableframe(D2Vec_Double &array, int size, double R, double R_c, double perc,  double c_high[], double c_low[]);
// NOTE: The following function REQUIRES the exprtk library to be included in the project.
void init_exprtk_homogenousMFTframe(D2Vec_Double &array, int size, double R, double R_c, double c_spread[]);
void init_exprtk_randbiMFTframe(D2Vec_Double &array, int size, double R, double R_c, double dP, double perc,  double c_spread[]);
void init_exprtk_readCSVcolumns_frame(D2Vec_Double &array, vector<int> &const_index, const string &parendir, const string &filenamePattern, 
					const std::vector<std::string> &read_cols, double a, double a_c,  int L, int r, double c_spread[]);

void init_csvconstframe(D2Vec_Double &array, D2Vec_Double &const_ind_val, const std::string& filename, const vector<int> &columns, int size);
void init_burnin_wrapper(D2Vec_Double &Rho_dt, double a, double a_c, int a_scalingfactor, double dP, double perc,  int L, int r, double c_spread[] );
void init_randconstframe(D2Vec_Double &array, int Sp, int size, double perc,   double c_high[], double c_low[]);
void init_constframe(D2Vec_Double &array,  int Sp, int size, double constant[]);
void init_randframe(D2Vec_Double &array, int Sp, int size, double mean[], double sd[]);
void init_quarterframe(D2Vec_Double &array, double c1, double c2, double c3, double c4, int L);
void init_gaussframe(D2Vec_Double &array, int size, vector<double> &sd, vector<double> &amp, int s_gauss = SpB);
void init_gaussiantear(D2Vec_Double &array, int Sp_lim, int g, double cx[], double cy[], const vector<double> &sd,const vector<double> &amp);
void init_solitarytear(D2Vec_Double &array, int length);

double mean_of_array(double array[],int size);
double standard_deviation_of_array(double array[],int size);
double mean_of_vector(vector<double> array,int size);
double variance_of_vector(vector<double> array,int size);
double standard_deviation_of_vector(vector<double> array,int size);
double occupied_sites_of_vector(vector<double> array,int size);
auto meansq_spread_of_vector(vector<double> array, int g, int c_x, int c_y);

void generic_SPBmean_surv_runs(D2Vec_Double &t_avg_var_rep_N, const D2Vec_Double &X_curr, int size, int spcolnum, int j);
void var_mean_incremental_surv_runs(D2Vec_Double &t_avg_var_rep_N, const D2Vec_Double &X_curr, int size, int j);

void var_mean_incremental_all_runs(double rep_avg_var[][3], double X_curr[][2], int size, int r);
void spreading_incremental_surv_runs(double t_N_R2_rep[][4], double X_curr[][3], int size);
void var_mean_incremental(double rep_avg_var[][2], vector<vector <double>> &X_curr, int size, int r);

// --------------------------- Partitioning functions to make life easier -----------------------------------//

template<typename T> std::vector<double> linspace(T start_in, T end_in, int num_in);
template<typename T>std::vector<double> lnspace(T start_in, T end_in, int log_points);
// By default, switchsort_and_bait() will simply sort the vector in ascending order.
template<typename T>std::vector<T> switchsort_and_bait(std::vector<T> vals, T low = 0, T high = 0, 
	int num_points =50, string insert_type = "None",  bool erase = false);
std::vector<double> logarithmic_time_bins(double t_max, double dt);
std::vector<double> logarithm10_time_bins(double t_max, double dt);

//----------------------------- Misc. Supporting Add-ons -------------------------------------------------

void recursive_dir_create(const string& path_to_dir);
int findMaxRepNo(const string& parendir, const string& filenamePattern);
int theborderwall(D2Vec_Double &Rho_t, int g);
void determine_neighbours_R2( int g, int S, D3Vec_Int &neighbours_R2);
void determine_neighbours_Sq4(int g, D3Vec_Int &neighbours_Sq4);
void determine_neighbours_Sq8(int g, int S, D3Vec_Int &neighbours_Sq8);
void set_density_dependentmortality(double m);
std::vector<std::pair<int, int>> computeNeighboringSitesCentral(int R);
std::vector<int> getNeighborsInBall(const vector<int>& sortedNeighboringSites, double f);
void generateNeighboringSitesFromCentral(const auto range_CentralNeighboringSites, 
										 vector<int>& neighboringSites, int i, int j, int L);
//void generateNeighboringSitesFromCentral(const vector<pair<int, int>>& centralNeighboringSites, 
//                                         vector<int>& neighboringSites, int i, int j, int L); 

//----------------------------- Regular DP Integration Machinery --------------------------------------------//

//----------------------------- Stochastic DP-Dornic Integration Machinery --------------------------------------------//

void save_prelimframe(D2Vec_Double &Rho_t, const string &parendir, const string &filenamePattern, double a, double a_st, double a_end, 
		double t, double dt, double dx, double dP, int r, int g, string header ="", bool overwrite = false, bool delete_previous = false);
void save_frame(D2Vec_Double &Rho_t, const string &parendir, const string &filenamePattern, double a, double a_st, double a_end, 
		double t, double dt, double dx, double dP, int r, int g, string header ="", bool overwrite = false);

void save_framefileswrapper(int index, int tot_iter, int j, int thrID,  double t, double dt, double dx, D2Vec_Double &Rho_dt, D2Vec_Double &DRho,
	D2Vec_Double &rho_rep_avg_var, vector<double> &t_meas, D2Vec_Double &gamma, D2Vec_Pair_Double &v_eff, double t_max, 
	double a, double b, double c, int a_scalingfactor, double (&D)[Sp], double (&v)[SpB], double sigma[], double a_st, double a_end, double a_c,
	double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double dP, int r, int g, double Gstar =-1., double Vstar = -1.);
int save_prelimfileswrapper(int index, int tot_iter, int j, int thrID,  double t, double dt, double dx, D2Vec_Double &Rho_dt, D2Vec_Double &DRho, D2Vec_Double &Rho_M, D2Vec_Double &Rho_Mov,
	D2Vec_Double &rho_rep_avg_var, vector<double> &t_meas, D2Vec_Double &gamma, D2Vec_Pair_Double &v_eff, double t_max, 
	double a, double b, double c, int a_scalingfactor, double (&D)[Sp], double (&v)[SpB], double sigma[], double a_st, double a_end, double a_c, 
	double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double dP, int r, int g, double Gstar =-1., double Vstar = -1.);

//------------------- Only Vegetation -------------------//

void f_DP_Dor_1Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double b, double c, double t, double dt, double g);

//------------------- Vegetation + Grazer-------------------//

// FFTW3 accelerated version [ O(N^2LOG(N)) ].
void calc_gamma_2Sp_NonRefugia_FFT(D2Vec_Double& Rho_t, D2Vec_Double& gamma, double (&Rho_avg)[Sp], vector<std::pair<double, int>>& rfrac, 
    vector<std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L, FFTW3_CentralPlanner& fft_central_plan);
// Conventional O(N^2*K^2) version.
void calc_gamma_2Sp_NonRefugia(const vector<pair<int, int>>& centralNeighboringSites, D2Vec_Double &Rho_t, D2Vec_Double &gamma,  
	double (&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, double nVeg_frac, int r_max, int L);
void f_DP_Dor_2Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double b, double c, 
	double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t, double dt, double g);


//------------------- Vegetation + Grazer + Predator -------------------//


//Calculates gamma for 3 Sp Rietkerk model (details in PDF, 
// ASSUMING PREDATORS AREN'T DETERRED BY HIGH LOCAL VEGETATION DENSITY)
// FFTW3 accelerated version [ O(N^2LOG(N)) ].
void calc_gamma_3Sp_NonRefugia_FFT(D2Vec_Double& Rho_t, D2Vec_Double& gamma, double (&Rho_avg)[Sp], vector<std::pair<double, int>>& rfrac, 
    vector<std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L, FFTW3_CentralPlanner& fft_central_plan);
// Conventional O(N^2*K^2) version.
void calc_gamma_3Sp_NonRefugia(const vector<pair<int, int>>& centralNeighboringSites, D2Vec_Double &Rho_t, D2Vec_Double &gamma,  
			double (&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, vector <std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L);
//Calculates gamma for 3 Sp Rietkerk model (details in PDF)
void calc_gamma_3Sp(const vector<pair<int, int>>& centralNeighboringSites, D2Vec_Double &Rho_t, D2Vec_Double &gamma,  
			double (&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, vector <std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L);
void f_DP_Dor_3Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double b, double c, 
	double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t, double dt, double g);
void RK4_Integrate_Stochastic_MultiSp(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D2Vec_Double &K1, D2Vec_Double &K2, D2Vec_Double &K3, D2Vec_Double &K4, D3Vec_Int &nR2,
	double b, double c, double (&Dxd2)[Sp], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t,double dt,double dx, int g);
void dP_Dornic_2D_MultiSp(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double b, double c, double (&D)[Sp], double (&v)[SpB],
    double sigma[], double a_st, double a_end, double a_c, double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double pR[], 
	int (&dtV)[SpB], double clow[], double dt, double dx, double dP, int r, int g, FFTW3_CentralPlanner& fft_central_plan, double Gstar = -1, double Vstar = -1);
void first_order_critical_exp_delta_stochastic_MultiSp(int div, double t_max, double a_start, double a_end, double a_c, double b,  double c, 
	double (&D)[Sp], double (&v)[SpB], double sigma[], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double pR[], int (&dtV)[SpB], double clo[],
	double dt, double dx, double dP, int r,  int g, double Gstar = -1, double Vstar = -1);

//============================ CUDA Kernels for DP-Dornic Model ===================================//

#if defined(__CUDACC__) || defined(BARRACUDA)

cudaError_t copyToDeviceConstantMemory_GammaSweepTerms(const DoubleIntPair* dr_frac);
void reportPerceptionArrays_Interface(int2* d_origin_Neighbourhood, double *d_rho_inverse, 
    int* d_eff_sizes, int nR_Perp_length, int r_max, int L);
__global__ void reportPerceptionArrays(int2* d_origin_Neighbourhood, double *d_rho_inverse, 
    int* d_eff_sizes, int nR_Perp_length, int r_max, int L);

__global__ void calc_gamma_3Sp_NonRefugia_kernel( const double* __restrict__ d_Rho_dt, double* __restrict__ d_gamma, std::pair<double, double>* d_v_eff,
    int2* d_origin_Neighbourhood, const double* __restrict__ d_rho_inverse, const int* __restrict__ d_eff_sizes, int L, int L2, int max_neighbors, double eps, double inv_eps);
void calc_gamma_vel_NonRefugia_CUDA(int2* d_origin_Neighbourhood, double* d_Rho_dt, double* d_gamma, std::pair<double, double>* d_v_eff, double *d_rho_inverse,
	 int* d_eff_sizes, double(&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, vector <std::pair<int, int>>& dtV_counter, int nR_Perp_length, int r_max, int L);
#endif

#endif