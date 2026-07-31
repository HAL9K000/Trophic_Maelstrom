#include "SPDP.h"
//#include "Debug.h"

int main()
{
  increase_stack_limit(2048L); //Increase stack limit to 512 MB.

  double p_start; double p_end; double a; double b= -2; double c= 1; int div; double t_max; double t_event;
  double D; double sigma; double dt; double dx; int r;  int g;
  

  cout << "Enter desired time-step (dt): ";
  cin >> dt;

  /**cout << "Enter desired step-distance (dx): ";
  cin >> dx;**/

  cout << "Enter time of frame capture (in hr): ";
  cin >> t_event;

  cout << "Enter the desired grid-size (L): ";
  cin >> g;

  cout << "Enter a-value: ";
  cin >> a;

  cout << "Enter starting initial-value (should be >= 0): ";
  cin >> p_start;

  cout << "Enter ending initial-value (consider setting less < 2): ";
  cin >> p_end;

  cout << "Enter number of divisions of p (n+1): ";
  cin >> div;

  cout << "Enter number of replicates ";
  cin >> r;

  /**cout << "Enter the desired stochastic coefficient (Sigma): ";
  cin >> sigma;

  cout << "Enter the desired diffusion coefficient (D): ";
  cin >> D;**/ 
  D=6.25; sigma =sqrt(1.0); dx=2.5;

  //D=0.25; sigma=0.25; div=11; a_start= 0.1; a_end = 0.11;
  //r=5; g=128; t_max= 100; dt=0.1; dx=0.5;

  //double Rho_0[g*g]; init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.

  //critical_exp_delta(Rho_0, div, t_max, a_start, a_end, b, D, sigma, dt, dx, r, g);

  //test_gamma_poisson();
  first_order_terminal_diff(div, t_event, p_start, p_end, a, b, c, D, sigma, dt, dx, r, g);
  
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