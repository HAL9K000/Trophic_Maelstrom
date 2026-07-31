#include "SPDP.h"
//#include "Debug.h"

int main()
{
  increase_stack_limit(2048L); //Increase stack limit to 512 MB.

  double a_start; double a_end; double b= -2; double c= 1; int div; double t_max; double t_event;
  double D; double sigma; double dt; double dx; int r =5;  int g;
  

  cout << "Enter desired time-step (dt): ";
  cin >> dt;

  /**cout << "Enter desired step-distance (dx): ";
  cin >> dx;**/

  cout << "Enter maximum duration of simulation (in hr): ";
  cin >> t_max;

  cout << "Enter time of frame capture (in hr): ";
  cin >> t_event;

  cout << "Enter the desired grid-size (L): ";
  cin >> g;

  cout << "Enter starting a-value (should be >= -b^{2}/4c): ";
  cin >> a_start;

  cout << "Enter ending a-value: ";
  cin >> a_end;

  cout << "Enter number of divisions of a (n+1): ";
  cin >> div;

  cout << "Enter number of replicates per value of measured a (R=5 is a good choice): ";
  cin >> r;

  /**cout << "Enter the desired stochastic coefficient (Sigma): ";
  cin >> sigma;

  cout << "Enter the desired diffusion coefficient (D): ";
  cin >> D;**/ 
  D=6.25; sigma =sqrt(1.0); dx=2.5;

  //D=0.25; sigma=0.25; div=11; a_start= 0.1; a_end = 0.11;
  //r=5; g=128; t_max= 100; dt=0.1; dx=0.5;

  double Rho_0[g*g]; init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.

  //critical_exp_delta(Rho_0, div, t_max, a_start, a_end, b, D, sigma, dt, dx, r, g);

  //test_gamma_poisson();

  
  capture_frame_decay(Rho_0, div, t_max, t_event, a_start, a_end, b, c, D, sigma, dt, dx, r,  g);
  
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