
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
#include <cassert>

#if defined(BARRACUDA) || defined(__CUDACC__)
#include <cuda_runtime.h>
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


/** Defining a compiler macro SPB if not already defined.
 * This is used to determine the number of species in the system (default =3)
 * This is done by including the appropriate header file based on the value of SPB, which in turn
 * defines the appropriate (inline) constants for the system.
 * This is done by passing the compiler flag -DSPB=2 or -DSPB=3 to the compiler.
*/
#ifndef SPB
	#define SPB 3 //Number of biota species in the system (default value used if NOT defined).
#endif

#ifndef INIT
	#define INIT 1 //Initial condition for the system (1 = Random Heterogenous Speckles, 0 = Homogeneous, 2 = Burn-in)
#endif

#ifndef MV_INVARIANCE
	#define MV_INVARIANCE 0 //Flag to determine movement invariance of the system [ 0 = Distance invariant, 1 = Velocity invariant]
#endif

#if SPB == 3
#include "rietkerk_bjork_constants_3Sp.h"
#elif SPB == 2
#include "rietkerk_bjork_constants_2Sp.h"
#elif SPB == 1
#include "rietkerk_bjork_constants_1Sp.h"
#else
#error "Number of species not supported. Pass valid compiler flag as  -DSPB=1/2/3"
#endif


//------------------------- GLOBAL CONSTANTS --------------------------------------------------------------//


// CUDA constants (actually defined in the .cu kernel files)
#if defined(__CUDACC__) || defined(BARRACUDA)
inline const int MAX_GAUSS_D_STENCIL_SIZE = 200; // Maximum size of Gaussian stencil for distance FKE.
inline const int MAX_GAUSS_VXY_STENCIL_SIZE = 200; // Maximum size of Gaussian stencil for velocity FKE.
inline bool use_constmem_gaussstencils = true; // Flag to indicate whether to use constant memory for Gaussian stencils.
#ifdef DEFINE_CUDA_CONSTANTS
//__constant__ double CuA[CuSpNV] /** ={aij, ajm} */;
//__constant__ double CuH[CuSpNV] /** ={hij, hjm} */;
//__constant__ double CuE[CuSpNV] /** ={ej, em} */;
//__constant__ double CuM[CuSpNV] /** ={mj, mm} */;
//__constant__ double CuD[Sp] /** ={D0, D1, D2, D3, D4} */;
//__constant__ double CuK[3] /** ={K0, K1, K2} */;
//__constant__ double CuV[CuSpNV] /** ={v0, v1} */;

__constant__ int d_sigD_ScaleBounds[CuSpNV] /** ={sigD_ScaleBounds[0], sigD_ScaleBounds[1]} */;
__constant__ int2 d_sigvD_ScaleBounds[CuSpNV] /** ={sigvD_ScaleBounds[0], sigvD_ScaleBounds[1]} */;
__constant__ double d_mu_vel_prefactor[CuSpNV] /** ={mu_vel_prefactor[0], mu_vel_prefactor[1]} */;

__constant__ int d_size_gauss_D[CuSpNV]; __constant__ int d_size_gauss_VXY[CuSpNV];
/** Stores the cumulative sizes of the Gaussian stencil for each species for the distance FKE
* This is used to determine the size of the stencil (and thereby access it) for each species in the device kernel. */;
__constant__ float d_gaussian_stencil_D[MAX_GAUSS_D_STENCIL_SIZE]; 
__constant__ float d_gaussian_stencil_VXY[MAX_GAUSS_VXY_STENCIL_SIZE];
// Gaussian stencils for distance and velocity FKEs respectively.
#endif
#endif


//inline const int Sp = 5; //Total number of species in the system.
inline const int SpB = Sp - 2; //Number of biota species in the system.
inline const int Sp_NV = Sp - 3; //Number of grazer and predator species in system.
inline const int Sp4_1 = 4 * (Sp - 2) + 3; // Used for generating statistics on surviving runs
inline const int Sp2 = 2 * Sp; // Used for generating statistics on surviving runs
inline const int SpB2 = 2 * SpB; // Used for generating statistics on surviving runs
inline const int SpB3 = 3 * SpB; // Used for generating statistics on surviving runs
inline const int SpB4 = 4 * SpB; // Used for generating statistics on surviving runs

// Global user-defined parameters for the Rietkerk model.
inline double dt2 /** = dt/2.0*/; inline double dt6 /** = dt/6.0*/; inline double L2 /** =g*g */;
inline double dx2 /** = dx*dx*/; inline double dx1_2 /** = 1.0/(dx*dx)*/; inline double dxb2 /** = dx/2.0/ */; inline double dx1 /** =  1/dx */;
inline const double epsilon = 1.0e-10; //Small value used for numerical stability.

// Global user-defined parameters used to set movement parameter terms (SEE NOTES and Johnson et al. 2008) for more details).
inline double Dt_ay; /** Analytical time-step */
/** Velocity parameter, ~ 3.0/Dt_analytical */
inline vector<double> beta_vel(SpB, 0.0); //Vector of velocity parameters for velocity FKE.
// Gaussian terms for analytic FKE distributions.
inline vector<std::pair<double, double>> sigma_v(SpB, { 0.0, 0.0 }); //Vector of Gaussian terms for velocity FKE.
inline vector<double> sigma_D(SpB, 0.0); //Vector of Gaussian terms for distance given by Diffusion process
inline vector<std::pair<double, double>> sigma_vD(SpB, { 0.0, 0.0 }); //Vector of Gaussian terms for distance FKE given by velocity OU.
inline vector<int> sig_D_ScaledBounds(SpB, 0);
//The bounds for the scaled distance OU terms, generally SET to { int(4*sigma_D/dx)}.
inline vector<std::pair<int, int>> sig_vD_ScaledBounds(SpB, { 0, 0 });
//The bounds for the scaled velocity OU terms for the distance FPE, scaled by dx,
// generally SET to { int(4*sigma_vD.first/dx), int(4*sigma_vD.second/dx}).
inline std::vector<std::vector <double>> gaussian_stencil_D(SpB, vector<double>(1, 0.0));
inline std::vector<std::vector <double>> gaussian_stencil_VXY(SpB, vector<double>(1, 0.0));
// Gaussian stencils for distance and velocity FPEs respectively.

inline vector<double> mu_vel_prefactor(SpB, 0.0); // = (1/dx)*(Dt_ay - ((1 - exp_beta_Dt[s])/beta_vel[s])), useful to save computation time.

inline vector<double> exp_beta_Dt(SpB, 0.0); // = exp(-beta_vel*Dt_analytical) for velocity FKE.
inline vector<double> exp_2beta_Dt(SpB, 0.0); // = exp(-2*beta_vel*Dt_analytical) for velocity FKE.


// Global user-defined parameters for the predator-prey trophic chain.
inline double K_P; /** Carrying capacity of predator */
inline double K_P1; /** Inverse of carrying capacity of predator */

inline double input_T = 0; /** Time for input data */

// Rietkerk Model Parameters, set in set_global_user_Rietkerk_params() function.
//inline double cgmax /** = c*gmax */; inline double K2W0 /**= K[2]*W0 */;inline double A01H01 /** = A[0][1]*H[0][1] */;inline double A12H12 /**= A[1][2]*H[1][2] **/;

// Strings that are used for I/O operations and file naming in the simulations.
inline string prefix = "DiC-NREF-LI";
inline string frame_folder = "../Data/Rietkerk/Frames/Stochastic/" + std::to_string(SpB) + "Sp/" + prefix + "_"; //Folder to store frames.
inline string prelim_folder = "../Data/Rietkerk/Prelims/Stochastic/" + std::to_string(SpB) + "Sp/" + prefix + "_"; //Folder to store preliminary data.
inline const string frame_prefix = "/FRAME_P_c_DP_G_"; //Prefix for frame files.
//inline const string errorframe_prefix = "/ERROR_"; //Prefix for frame files.
inline const string gamma_prefix = "/GAMMA_G_"; //Prefix for gamma files.
inline const string prelim_prefix = "/PRELIM_AGGRAND_P_c_ID_"; //Prefix for preliminary data files.
inline const string replicate_prefix = "/PRELIM_TSERIES_P_c_DP_G_"; //Prefix for replicate time-series data files.
inline const string movement_prefix = "/PRELIM_MOVSERIES_P_c_DP_G_"; //Prefix for movement time-series data files.
//inline const string frame_header = "a_c,  x,  P(x; t), G(x; t), Pr(x; t), W(x; t), O(x; t) \n"; //Header for frame files.
inline const string input_prefix = "/FRAME_T_"; //Prefix for input files.
inline string stat_prefix = "../Data/Rietkerk/Stochastic/" + std::to_string(SpB) + "Sp/1stCC_Rietkerk_" + prefix + "_P_c_G_";

inline string input_frame_parenfolder = "../Input/Rietkerk/" + std::to_string(SpB) + "Sp"; //Parendirectory for input data.
inline string input_frame_subfolder = ""; //Subfolder for input data.
// NOTE: Complete path to input data is input_frame_parenfolder + input_frame_subfolder + input_prefix +"*.csv"
inline std::map<string, string> input_keys = { {"outPRE", ""}, {"a", "0"}, {"a_c", "0"}, {"T", "0"}, {"dP", "0"}, {"Geq", "0"}, {"R", ""},
{"MAMP", "0"}, {"SEP", "0"}, {"WID", "0"}, {"MSF", "0"} };

// Strings that are used for the Expertk library to parse mathematical expressions.
// These represent MFT-versions of the species Coexistance equilibria wrt to the order parameter (Rainfall, given by a).
inline vector<string> MFT_Vec_CoexExpr(2 * Sp, ""); //Vector of MFT Coexistance expressions for all species.
//inline exprtk::symbol_table<double> global_symbol_table; //Symbol table for the Expertk library.
inline vector<double> frame_tmeas; //Vector to store time measurements for the frames.


//------------------------- TYPEDEFS --------------------------------------------------------------//

typedef std::vector<std::vector <double>> D2Vec_Double;
typedef std::vector<std::vector <float>> D2Vec_Float;
typedef std::vector<std::vector <int>> D2Vec_Int;

typedef std::vector<std::vector <std::vector <double>>> D3Vec_Double;
typedef std::vector<std::vector <std::vector <float>>> D3Vec_Float;
typedef std::vector<std::vector <std::vector <int>>> D3Vec_Int;

typedef std::vector<std::vector <std::vector <std::vector <double>>>> D4Vec_Double;
typedef std::vector<std::vector <std::vector <std::vector <float>>>> D4Vec_Float;
typedef std::vector<std::vector <std::vector <std::vector <int>>>> D4Vec_Int;

typedef std::vector<std::vector<std::pair<double, double>>> D2Vec_Pair_Double;
typedef std::vector<std::vector<std::pair<float, float>>> D2Vec_Pair_Float;
typedef std::vector<std::vector<std::pair<int, int>>> D2Vec_Pair_Int;

typedef std::vector<double> Vec_Double;
typedef std::vector<float> Vec_Float;
typedef std::vector<int> Vec_Int;





//---------------------Some Primitives ------------------------------------------------------ //

//Stolen from: https://web.archive.org/web/20200629174939/http://timmurphy.org/2013/06/27/template-metaprogramming-in-c/
// Generic exponent compile-time calculations. Pow<2,3>::result == 2^3
template < long long B, long long E>
struct Pow
{
static const long long result = B * Pow<B, E - 1>::result;
};

template <long long B>
struct Pow<B, 0>
{
static const long long result = 1;
};

// Generic Permutation compile-time calculations. Perm<N,K>::result == N!/(N-K)!
template <long long N, long long K>
struct Perm
{
static const long long result = N * Perm<N - 1, K - 1>::result;
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
static const long long result = Perm<N, K>::result / Perm<K, K>::result;
};

template <long long N>
struct Comb<N, 0>
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
		1.0 / Pow<16, N>::result *
		(
			4.0 / (8 * N + 1.0) - 2.0 / (8 * N + 4.0) -
			1.0 / (8 * N + 5.0) - 1.0 / (8 * N + 6.0)
			)
		) + CalculatePi<N - 1>::pi;
};

template <>
struct CalculatePi<-1>
{
static const constexpr  double pi = 0;
};

// Generic template to 

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

inline const double PI = CalculatePi<14>::pi;

inline const int SpB_Mov = 2 * SpB + Comb<SpB, 2>::result; // Used for generating statistics on movement of surviving runs

//void set_global_user_Rietkerk_params(double c,double gmax,double alpha, double rW, double W0, double (&D)[Sp], 
//	double (&K)[3], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB])

template <typename T> int sgn(T val);
template <typename T> int sgn_index(T val);

template <typename T> T min_positive(const std::vector<T>& vec); // Template function to return smallest positive value in a vector
template <typename T> T min_positive_arr(T* arr, int size); // Template function to return smallest positive value in an array
void increase_stack_limit(long long stack_size);
bool maxis(int a, int b);
string format_str(const std::string& s, const std::map<string, string>& values);
void add_three(int a, int b, int c); //Test function.
void set_Prefix(string& user_prefix, double mG = 0.0, double mP = 0.0);
void set_input_Prefix(string& user_inputframepath, string& user_prefix, double user_a_c, double user_dP, double user_Geq = -1, double user_input_T = -1);
void set_global_system_params(double dt, double dx);
void set_global_FKE_params(double(&D)[Sp], double(&v)[SpB], double(&dtAdv)[SpB], double Dt_ay, double dx, int L, int CI = 4);
void set_global_predator_params(double Km);

void display_symbol_table(const exprtk::symbol_table<double>& symbol_table);
// Custom comparator for sorting pairs based on the squared distance
bool comparePairs(const pair<int, int>& a, const pair<int, int>& b);

void init_fullframe(D2Vec_Double& array, int Sp, int size);

void init_randbistableframe(D2Vec_Double& array, int size, double R, double R_c, double perc, double c_high[], double c_low[]);
// NOTE: The following function REQUIRES the exprtk library to be included in the project.
void init_exprtk_homogenousMFTframe(D2Vec_Double& array, int size, double R, double R_c, double c_spread[]);
void init_exprtk_randbiMFTframe(D2Vec_Double& array, int size, double R, double R_c, double dP, double perc, double c_spread[]);
void init_exprtk_readCSVcolumns_frame(D2Vec_Double& array, vector<int>& const_index, const string& parendir, const string& filenamePattern,
const std::vector<std::string>& read_cols, double a, double a_c, int L, int r, double c_spread[]);
void init_exprtk_randbiMFTframe_OLD(D2Vec_Double& array, int size, double R, double R_c, double dP, double perc, double c_spread[]);


void init_csvconstframe(D2Vec_Double& array, D2Vec_Double& const_ind_val, const std::string& filename, const vector<int>& columns, int size);
void init_burnin_wrapper(D2Vec_Double& Rho_dt, double a, double a_c, double dP, double perc, int L, int r, double c_spread[]);
void init_randconstframe(D2Vec_Double& array, int Sp, int size, double perc, double c_high[], double c_low[]);
void init_constframe(D2Vec_Double& array, int Sp, int size, double constant[]);
void init_randframe(D2Vec_Double& array, int Sp, int size, double mean[], double sd[]);
void init_quarterframe(D2Vec_Double& array, double c1, double c2, double c3, double c4, int L);
void init_gaussframe(D2Vec_Double& array, int size, vector<double>& sd, vector<double>& amp, int s_gauss = SpB);
void init_gaussiantear(D2Vec_Double& array, int Sp_lim, int g, double cx[], double cy[], const vector<double>& sd, const vector<double>& amp);
void init_solitarytear(D2Vec_Double& array, int length);

double mean_of_array(double array[], int size);
double standard_deviation_of_array(double array[], int size);
double mean_of_vector(vector<double> array, int size);
double variance_of_vector(vector<double> array, int size);
double standard_deviation_of_vector(vector<double> array, int size);
double occupied_sites_of_vector(vector<double> array, int size);
auto meansq_spread_of_vector(vector<double> array, int g, int c_x, int c_y);

void generic_SPBmean_surv_runs(D2Vec_Double& t_avg_var_rep_N, const D2Vec_Double& X_curr, int size, int spcolnum, int j);
void var_mean_incremental_surv_runs(D2Vec_Double& t_avg_var_rep_N, const D2Vec_Double& X_curr, int size, int j);

void var_mean_incremental_all_runs(double rep_avg_var[][3], double X_curr[][2], int size, int r);
void spreading_incremental_surv_runs(double t_N_R2_rep[][4], double X_curr[][3], int size);
void var_mean_incremental(double rep_avg_var[][2], vector<vector <double>>& X_curr, int size, int r);

// --------------------------- Partitioning functions to make life easier -----------------------------------//

template<typename T> std::vector<double> linspace(T start_in, T end_in, int num_in);
template<typename T>std::vector<double> lnspace(T start_in, T end_in, int log_points);
template<typename T>std::vector<T> switchsort_and_bait(std::vector<T> vals, T low = 0, T high = 0, int num_points = 50,
string insert_type = "None", bool erase = false);
std::vector<double> logarithmic_time_bins(double t_max, double dt);
std::vector<double> logarithm10_time_bins(double t_max, double dt);

//----------------------------- Misc. Supporting Add-ons -------------------------------------------------

void recursive_dir_create(const string& path_to_dir);
int findMaxRepNo(const string& parendir, const string& filenamePattern);
int theborderwall(D2Vec_Double& Rho_t, int g);
void determine_neighbours_R2(int g, int S, D3Vec_Int& neighbours_R2);
void determine_neighbours_Sq4(int g, D3Vec_Int& neighbours_Sq4);
void determine_neighbours_Sq8(int g, int S, D3Vec_Int& neighbours_Sq8);
void set_density_dependentmortality(double m);
std::vector<std::pair<int, int>> computeNeighboringSitesCentral(int R);
std::vector<int> getNeighborsInBall(const vector<int>& sortedNeighboringSites, double f);
void generateNeighboringSitesFromCentral(const auto range_CentralNeighboringSites,
vector<int>& neighboringSites, int i, int j, int L);
//void generateNeighboringSitesFromCentral(const vector<pair<int, int>>& centralNeighboringSites, 
//                                         vector<int>& neighboringSites, int i, int j, int L); 
//Create Gaussian stencil of size Rbound + 1 to store precomputed values for FKE calculations.
std::vector<double> precompute_gaussian(double sigma, int Rbound);

//----------------------------- Regular Rietkerk Integration Machinery --------------------------------------------//

void f_2D(D2Vec_Double& f, D2Vec_Double& Rho_M, D3Vec_Int& nR2, double a, double c, double gmax, double alpha, double d, double rW, double W0,
double D[], double K[], double t, double dt, double dx, int g);
void RK4_Integrate(D2Vec_Double& Rho_t, D2Vec_Double& Rho_tsar, double a, double b, double D[], double A[Sp][Sp],
double H[Sp][Sp], double E[], double M[], double t, double dt, double dx, int g);
void RK4_Wrapper_2D(D2Vec_Double& Rho, vector <double>& t_meas, double t_max, double a, double c, double gmax,
double alpha, double d, double rW, double W0, double D[], double K[], double a_st, double a_end, double dt, double dx, double dP, int r, int g);

void first_order_critical_exp_delta(int div, double t_max, double a_start, double a_end, double c, double gmax, double alpha,
double d, double rW, double W0, double D[], double K[], double dt, double dx, double dP, int r, int g);


//----------------------------- Stochastic Rietkerk-Dornic Integration Machinery --------------------------------------------//

//------------------- Only Vegetation (+ Soil Water + Surface Water) -------------------//

void f_2Dor(D2Vec_Double& f, D2Vec_Double& Rho_M, D3Vec_Int& nR2, double a, double c, double gmax,
double alpha, double rW, double W0, double(&Dxd2)[Sp], double(&K)[3], double t, double dt, double g);
void RK4_Integrate_Stochastic(D2Vec_Double& Rho_dt, D2Vec_Double& Rho_tsar, D3Vec_Int& nR2, double a, double c, double gmax, double alpha,
double d, double rW, double W0, double D[], double K[], double t, double dt, double dx, int g);
void rietkerk_Dornic_2D(D2Vec_Double& Rho, vector <double>& t_meas, double t_max, double a, double c, double gmax, double alpha, double d,
double rW, double W0, double D[], double K[], double sigma[], double a_st, double a_end, double dt, double dx, double dP, int r, int g);
void first_order_critical_exp_delta_stochastic(int div, double t_max, double a_start, double a_end, double c, double gmax, double alpha,
double d, double rW, double W0, double D[], double K[], double sigma[], double dt, double dx, double dP, int r, int g);

//------------------- Vegetation + Grazer (+ Soil Water + Surface Water) -------------------//

void f_DorFKE_2Sp(D2Vec_Double& f, D2Vec_Double& Rho_M, D3Vec_Int& nR2, double a, double c, double gmax, double alpha, double rW,
double W0, double(&Dxd2)[Sp], double(&K)[3], double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB],
double t, double dt, double dx1_2, double g);
void RK4_Integrate_Stochastic_2Sp(D2Vec_Double& Rho_t, D2Vec_Double& Rho_tsar, D2Vec_Double& K1, D2Vec_Double& K2, D2Vec_Double& K3, D2Vec_Double& K4, D3Vec_Int& nR2,
double a, double c, double gmax, double alpha, double rW, double W0, double(&Dxd2)[Sp], double(&K)[3], double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB],
double t, double dt, double dx, int g);
void calc_gamma_vel_2Sp_NonRefugia(const vector<pair<int, int>>& centralNeighboringSites, std::mt19937_64& rng, D2Vec_Double& Rho_t, D2Vec_Double& gamma,
D2Vec_Pair_Double& v_eff, double(&v)[SpB], double(&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, double nVeg_frac, int r_max, int L);

//------------------- Vegetation + Grazer + Predator (+ Soil Water + Surface Water) -------------------//
void save_prelimframe(D2Vec_Double& Rho_t, const string& parendir, const string& filenamePattern, double a, double a_st, double a_end,
double t, double dt, double dx, double dP, int r, int g, string header = "", bool overwrite = false, bool delete_previous = false);
void save_frame(D2Vec_Double& Rho_t, const string& parendir, const string& filenamePattern, double a, double a_st, double a_end,
double t, double dt, double dx, double dP, int r, int g, string header = "", bool overwrite = false);

void save_framefileswrapper(int index, int tot_iter, int j, int thrID, double t, double dt, double dx, D2Vec_Double& Rho_dt, D2Vec_Double& DRho,
D2Vec_Double& rho_rep_avg_var, vector<double>& t_meas, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, double t_max,
double a, double c, double gmax, double alpha, double rW, double W0, double(&D)[Sp], double(&v)[SpB],
double(&K)[3], double sigma[], double a_st, double a_end, double a_c, double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB],
double dP, int r, int g, double Gstar /* =-1.*/, double Vstar /* = -1.*/);

int save_prelimfileswrapper(int index, int tot_iter, int j, int thrID, double t, double dt, double dx, D2Vec_Double& Rho_dt, D2Vec_Double& DRho, D2Vec_Double& Rho_M, D2Vec_Double& Rho_Mov,
D2Vec_Double& rho_rep_avg_var, vector<double>& t_meas, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, double t_max,
double a, double c, double gmax, double alpha, double rW, double W0, double(&D)[Sp], double(&v)[SpB],
double(&K)[3], double sigma[], double a_st, double a_end, double a_c, double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB],
double dP, int r, int g, double Gstar /* =-1.*/, double Vstar /* = -1.*/);


// A wrapper function to save the preliminary data for the Rietkerk model. Returns 1 if detects no vegetation left, else 0.
int save_fileswrapper(int index, int tot_iter, int j, int thrID, double t, double dt, double dx, D2Vec_Double& Rho_dt, D2Vec_Double& DRho, D2Vec_Double& Rho_M, D2Vec_Double& Rho_Mov,
D2Vec_Double& rho_rep_avg_var, vector<double>& t_meas, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, double t_max,
double a, double c, double gmax, double alpha, double rW, double W0, double(&D)[Sp], double(&v)[SpB], double(&K)[3], double sigma[], double a_st, double a_end, double a_c,
double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double dP, int r, int g, double Gstar = -1., double Vstar = -1.);

//Calculates gamma and velocity for 3 Sp Rietkerk model (details in PDF, 
// ASSUMING PREDATORS AREN'T DETERRED BY HIGH LOCAL VEGETATION DENSITY)
void calc_gamma_vel_3Sp_NonRefugia(const vector<pair<int, int>>& centralNeighboringSites, std::mt19937_64& rng, D2Vec_Double& Rho_t, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, double(&v)[SpB],
double(&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, vector <std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L);
//Calculates gamma for 3 Sp Rietkerk model (details in PDF)
void calc_gamma_vel_3Sp(const vector<pair<int, int>>& centralNeighboringSites, std::mt19937_64& rng, D2Vec_Double& Rho_t, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, double(&v)[SpB],
double(&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, vector <std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L);

void f_DorFKE_3Sp(D2Vec_Double& f, D2Vec_Double& Rho_M, D3Vec_Int& nR2, double a, double c, double gmax, double alpha, double rW, double W0,
double(&Dxd2)[Sp], double(&K)[3], double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double t, double dt, double dx1_2, double g);
void RK4_Integrate_Stochastic_MultiSp(D2Vec_Double& Rho_t, D2Vec_Double& Rho_tsar, D3Vec_Int& nR2, double a, double c, double gmax, double alpha,
double rW, double W0, double(&D)[Sp], double(&K)[3], double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double t, double dt, double dx, int g);

void advdiff_FKE_MultiSp(D2Vec_Double& Rho_t, D2Vec_Double& Rho_tsar, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, int g/**D2Vec_Double& gaussian_stencil_D, D2Vec_Double& gaussian_stencil_VXY, */);

void rietkerk_DorFPE_2D_MultiSp(D2Vec_Double& Rho, vector <double>& t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0,
double(&D)[Sp], double(&v)[SpB], double(&K)[3], double sigma[], double a_st, double a_end, double a_c, double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double pR[],
int(&dtV)[SpB], double clow[], double dt, double dx, double dP, int r, int g, double Gstar = -1.0, double Vstar = -1.0);
void first_order_critical_exp_delta_stochastic_MultiSp(int div, double t_max, double a_start, double a_end, double a_c, double c, double gmax, double alpha,
double rW, double W0, double(&D)[Sp], double(&v)[SpB], double(&K)[3], double sigma[], double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double pR[],
int(&dtV)[SpB], double clo[], double dt, double dx, double dP, int r, int g, double Gstar = -1.0, double Vstar = -1.0);

//============================ CUDA Kernels for Rietkerk-Dornic Model ===================================//

void reportConstantMemory();
// CUDA kernel to launch the FKE- Advection Diffusion CUDA kernel for the Rietkerk-Dornic model for Multi-Species systems
// NOTE: DEFAULT ARGUMENTS ASSUME USE OF CONSTANT MEMORY FOR GAUSSIAN STENCILS!
void launchAdvDiffKernel_MultiSp(double* Rho_t, double* Rho_tsar, double* gamma, std::pair<double, double>* v_eff, int L,
double* gaussian_stencil_D = nullptr, double* gaussian_stencil_VXY = nullptr);

#if defined(__CUDACC__) || defined(BARRACUDA)
// Copy host data for advection terms to device constant memory
cudaError_t copyToDeviceConstantMemory_AdvTerms(const int* sig_D_ScaledBounds, const int2* sigma_vD_ScBounds, const double* mu_vel_prefactor, const int* d_size_gauss_D, const int* d_size_gauss_VXY);
cudaError_t copyToDeviceConstantMemory_GaussianStencils(size_t total_size_gauss_D, size_t total_size_gauss_VXY);
template <typename T> void reportGaussianStencilValues(T* d_gaussian_stencil_D, T* d_gaussian_stencil_VXY);
#endif


#if defined(__CUDACC__)
// AtomicAdd function analogue for double precision floating point numbers for devices with compute capability < 6.0
__device__ double atomicAdd_Double(double* address, double val);
__global__ void reportConstantMemory_AdvTerms();
template <typename T, bool useConstantGaussMemory> __global__ void reportGaussianStencils(T* gaussian_stencil_D = nullptr, T* gaussian_stencil_VXY = nullptr);
// CUDA kernel to perform the FKE- Advection Diffusion step for the Rietkerk-Dornic model for Multi-Species systems.
// Template argument to switch between using constant memory (useConstantGaussMemory= true) or global memory for Gaussian stencils.
template<bool useConstantGaussMemory> __global__ void updateMovement(double* Rho_t, double* Rho_tsar, double* gamma, std::pair<double, double>* v_eff, int L,
double* gaussian_stencil_D = nullptr, double* gaussian_stencil_VXY = nullptr);
#endif


#endif