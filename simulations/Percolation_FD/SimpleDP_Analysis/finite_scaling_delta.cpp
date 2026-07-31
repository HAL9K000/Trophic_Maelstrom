#include "SPDP.h"

int main()
{
  increase_stack_limit(2048L); //Increase stack limit to 512 MB.

  double a_c; double b=1; int div; double t_max; double t_lim=0;
  double D; double sigma; double dt; double dx; int r;  int g_st; int g_end;
  

  cout << "Enter desired time-step (dt): ";
  cin >> dt;

  /**cout << "Enter desired step-distance (dx): ";
  cin >> dx;**/

  cout << "Enter maximum duration of simulation (end of time-window) (in hr): ";
  cin >> t_max;

  //cout << "Enter starting time-window of capture (in hr): ";
  //cin >> t_lim;

  cout << "Enter starting grid-size (L1): ";
  cin >> g_st;

  cout << "Enter ending grid-size (L2): ";
  cin >> g_end;

  cout << "Enter number of divisions of grid-sizes (n+1): ";
  cin >> div;

  cout << "Enter the desired number of infra-replicates (r): ";
  cin >> r;

  cout << "Enter critical threshold: ";
  cin >> a_c;

  /**cout << "Enter the desired stochastic coefficient (Sigma): ";
  cin >> sigma;

  cout << "Enter the desired diffusion coefficient (D): ";
  cin >> D;**/ 
  D=0.25; sigma =sqrt(2.0); dx=1.0;

  //D=0.25; sigma=0.25; div=11; a_start= 0.1; a_end = 0.11;
  //r=5; g=128; t_max= 100; dt=0.1; dx=0.5;

  //Left half of random grid has initial concentration 1/3rd that of the right half.
  //critical_exp_delta(Rho_0, div, t_max, a_start, a_end, b, D, sigma, dt, dx, r, g);
  //critical_exp_finite_scaling_collapse(div, g_st, g_end, t_max,  a_c, b, D, sigma, dt, dx, r);
  critical_exp_finite_scaling_stationary(div, g_st, g_end, t_lim, t_max,  a_c, b, D, sigma, dt, dx, r);
  return 0;
  
}