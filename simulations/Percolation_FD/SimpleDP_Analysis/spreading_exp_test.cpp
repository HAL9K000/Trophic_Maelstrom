#include "SPDP.h"

int main()
{
  increase_stack_limit(2048L); //Increase stack limit to 512 MB.

  double a_c; double b=1; int div; double t_max;
  double D; double sigma; double dt; double dx; int r;  int g;
  

  cout << "Enter desired time-step (dt): ";
  cin >> dt;

  /**cout << "Enter desired step-distance (dx): ";
  cin >> dx;**/

  cout << "Enter maximum duration of simulation (in hr): ";
  cin >> t_max;

  cout << "Enter the desired grid-size (L): ";
  cin >> g;

  cout << "Enter the desired number of infra-replicates (r): ";
  cin >> r;

  cout << "Enter number of divisions of supra-replicates (n+1): ";
  cin >> div;

  cout << "Enter critical threshold: ";
  cin >> a_c;

  /**cout << "Enter the desired stochastic coefficient (Sigma): ";
  cin >> sigma;

  cout << "Enter the desired diffusion coefficient (D): ";
  cin >> D;**/ 
  D=0.25; sigma =sqrt(2.0); dx=1.0;

  //D=0.25; sigma=0.25; div=11; a_start= 0.1; a_end = 0.11;
  //r=5; g=128; t_max= 100; dt=0.1; dx=0.5;

  double Rho_0[g*g]; init_solitarytear(Rho_0, g); //Returns Rho_0 with a full initial frame filled with ones.

  //Left half of random grid has initial concentration 1/3rd that of the right half.
  //percolationHalfDornic_2D(Rho, t_stop, Rho_0, t_max,  p_start, D, sigma,  dt, dx, r,  g);
  //critical_exp_delta(Rho_0, div, t_max, a_start, a_end, b, D, sigma, dt, dx, r, g);
  critical_exp_theta_delta(Rho_0, div, t_max, a_c, b, D, sigma, dt, dx, r, g);

  //test_gamma_poisson();
  
  
  cout<< "INITIAL CONDITIONS:" <<endl;

  for(int i=0; i< g; i++)
  {   //
    for(int j=0; j< g; j++)
      {
        //cout << "(" << i <<"," << j << ") " << Rho_0[g*i +j] << " ";
      }
  }
  cout << endl;
  return 0;
  
}