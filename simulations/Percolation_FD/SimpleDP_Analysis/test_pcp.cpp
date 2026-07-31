#include "SPDP.h"

int main()
{
  increase_stack_limit(128L); //Increase stack limit to 128 MB.

  double p_start; double p_end; int div; double t_max;
  double D; double sigma; double dt; double dx; int r;  int g;
  

  cout << "Enter desired time-step (dt): ";
  cin >> dt;

  cout << "Enter desired step-distance (dx): ";
  cin >> dx;

  cout << "Enter maximum duration of simulation (in hr): ";
  cin >> t_max;

  /** cout << "Enter the desired grid-size (L): ";
  cin >> g;

  cout << "Enter the desired number of replicates (r): ";
  cin >> r; **/

  cout << "Enter starting p-value: ";
  cin >> p_start;

  /** cout << "Enter ending p-value: ";
  cin >> p_end;

  cout << "Enter number of divisions of p (n+1): ";
  cin >> div; **/

  cout << "Enter the desired stochastic coefficient (Sigma): ";
  cin >> sigma;

  cout << "Enter the desired diffusion coefficient (D): ";
  cin >> D;

  r= 3; g=5;
  double Rho_0[g*g] ={0.75};
  for(int i=0; i< g; i++)
  {   //
    for(int j=0; j< g; j++)
    {
      Rho_0[g*i + j] = 0.1;
      
    }
  }
  for(int i=0; i< g; i++)
  {   //
    for(int j=0; j< g/2.0; j++)
      Rho_0[g*i + j] = 0.1;
  }
  
  vector<double> t_stop = {t_max/3.0, t_max}; vector <vector <double>> Rho;
  //Left half of random grid has initial concentration 1/3rd that of the right half.
  percolation_2D(Rho, t_stop, Rho_0, t_max,  p_start, D, sigma,  dt, dx, r,  g);
  //stoc_dem_percolation(Rho_0, t_stop, div, t_max, p_start, p_end, D, sigma, dt, dx, r, g);


  cout << " p , # Tr No , t ,  i*L + j ,  Rho_ij(t) \n";
  for(int i=0; i< Rho.size(); i++)
  {
    cout << setprecision(5) << Rho[i][0] << "," << static_cast<int>(Rho[i][1]) << "," << setprecision(7)
	    << Rho[i][2] << "," << Rho[i][3] << "," << setprecision(16) << Rho[i][4] << endl;
  }
  
  cout<< "INITIAL CONDITIONS:" <<endl;

  for(int i=0; i< g; i++)
  {   //
    for(int j=0; j< g; j++)
      { int ten = g*i +j;
        cout << "(" << i <<"," << j << ") " << Rho_0[g*i +j] << " ";}
  }
  cout << endl;
  return 0;
  
}