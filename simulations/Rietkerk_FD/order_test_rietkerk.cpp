#include "rietkerk_bjork_basic.h"

int main(int argc, char *argv[])
{
  increase_stack_limit(2048L); //Increase stack limit to 512 MB.

  string preFIX; // Prefix for the output files
  double a =0.575;
  double c, gmax, alpha, d, rW, W0, t_max, dt, dx; //int g, div;
  double a_start, a_end; double r; double dP; // Kick for high initial state
  int g, div;

  c= 10; d =0.25; gmax= 0.05; alpha =0.2; W0= 0.2; rW = 0.2; // From Bonachela et al 2015

  double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2

  double k0, k1, k2; double d0, d1, d2; //Slightly 
  
  dx= 0.2 ; //From Bonachela et al 2015 (in m)
  d0 = 0.001; d1 = 0.001; d2= 0.1; //From Bonachela et al 2015 (in m^2/day)
  k0= 0; k1 = 5; k2 =5;

  double D[Sp] ={d0, d1, d2}; //Diffusion coefficients for species.
  double K[3] ={k0, k1, k2}; //Diffusion coefficients for species.

  /**

  double beta = c*gmax - d;
  double p0istar, p0jstar, p0mstar; // Analytic steady state values.

  p0jstar = d*K[1]/beta; 
  p0istar = (c/d)*(R - rW*p0jstar);
  p0mstar = (R/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
  */
  
  if(argc == 1)
  {
    cout << "Enter desired time-step: ";
    cin >> dt;

    cout << "Enter maximum duration of simulation (in days): ";
    cin >> t_max;

    cout << "Enter the desired grid-size (L) (dx = 0.2 m): ";
    cin >> g;

    cout << "Enter the desired number of replicates (r): ";
    cin >> r;

    cout << "Note that a relevant value of 'a' may be construed as:\t" << setprecision(10) << a <<endl;

    cout << "Enter starting a-value: ";
    cin >> a_start;

    cout << "Enter ending a-value: ";
    cin >> a_end;

    cout << "Enter number of divisions of p (n+1): ";
    cin >> div;

    cout << "Enter kick for high state: ";
    cin >> dP;

    cout << "Enter Prefix (common choices include 'DiC-REF', 'DDM-DiC-BURNIN', 'nDiC' etc.) If not provided, default value: " 
        + prefix + " will be used: ";
    cin >> preFIX;
  }
  // If the user has entered the correct number of arguments, then the arguments are read from the command line.
  else if(argc == 10)
  {
    dt = atof(argv[1]);
    t_max = atof(argv[2]);
    g = atoi(argv[3]);
    r = atoi(argv[4]);
    a_start = atof(argv[5]);
    a_end = atof(argv[6]);
    div = atoi(argv[7]);
    dP = atof(argv[8]);
    preFIX = argv[9];
  }
  else
  {
    cout << "Please enter the correct number of arguments.\n";
    cout << "The correct number of arguments is 9.\n";
    cout << "The arguments are as follows: dt, t_max, g, r, a_start, a_end, div, dP, Prefix\n";
    exit(1);
  }

  //cout << "Enter 1 to read from csv or 0 to input standardised initial condions"

  /* cout << "Enter initial density of Primary Producer (kg/m^2): ";
  cin >> p0i;

  cout << "Enter initial density of Primary Consumer (kg/m^2): ";
  cin >> p0j;

	cout << "Enter initial density of Secondary Producer (kg/m^2): ";
  cin >> p0m; */

  //add_two(3.14, 1.78);
  cout << maxis(3.14, 1.78) <<endl;

  int go =8;
  add_three(1,2,go);

  set_global_system_params(dt, dx); //Set the global parameters for the simulation.
  cout << "Global parameters set, with dt/2.0 = " << dt2 << " and dx*dx = " << dx2 <<  " and 1/(dx*dx) = " << dx1_2 << "\n";

  cout << "\n=============================================================\n";

  if(preFIX == "N/A" || preFIX == "n/a" || preFIX == "NA" || preFIX == "na" || preFIX == "N" || preFIX == "n" || preFIX == "")
  {
    cout << "No prefix provided. Using default prefix: " << prefix << endl;
    preFIX = prefix;
  }

  set_Prefix(preFIX, 0.0, 0.0, false); //Set the prefix (and thus save folders) to the user-defined prefix.
  cout << "Save directory for frames: " << frame_folder << "\n";
  cout << "Save directory for preliminary data: " << prelim_folder << "\n";
  cout << "Save directory for final data: " << stat_prefix << "\n";

  //D2Vec_Double Rho_0(Sp, vector<double> (g*g, 1.0));
  // Rho_0 is 2D [Spx(L*L)] vector initialised to 1.00.

  //init_fullframe(Rho_0, Sp, g*g); //Returns Rho_0 with a full initial frame filled with ones.
  //double p0i = 1.0; double p0j= 0.05;`
  //double mean[Sp] = {p0i, p0j}; double sd[Sp] = {p0i/4.0, p0j/4.0};
	//init_randframe(Rho_0, Sp,  g*g, mean, sd); //Returns Rho_0 with a full initial frame filled with 0.2.
  //double mean[Sp] = {p0i, p0j, p0m}; double sd[Sp] = {p0i/2.0, p0j/2.0, p0m/2.0};
  //init_randframe(Rho_0, Sp,  g*g, mean, sd); //Returns Rho_0 with a full random frame.

  //double perc = 0.015; double c_high[Sp] ={p0i, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

  //init_randconstframe(Rho_0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.

  //init_constframe(Rho_0, Sp,  g*g, mean); //Returns Rho_0 with a full initial frame filled with 0.2.

  stringstream ast, est, dPo; ast << a_start; est  << a_end; dPo << dP;

  //string path_to_folder = frame_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Geq_" + geq.str();
  //stringstream foldername;
	//foldername << "../Data/Rietkerk/Frames/Non-Stochastic/" << ast.str() << "-" << est.str() << "_dP_" << dPo.str() << "/";
  //string path_to_folder = "../Data/Rietkerk/Frames/Non-Stochastic/" + ast.str() + "-" + est.str() + "_dP_" + dPo.str();
  string path_to_folder = frame_folder +  ast.str() + "-" + est.str() + "_dP_" + dPo.str();
  recursive_dir_create(path_to_folder);

  //stringstream foldername2;
	//foldername2 << "../Data/Rietkerk/Prelims/Non-Stochastic/" << ast.str() << "-" << est.str() << "_dP_" << dPo.str() << "/";
  //path_to_folder = "../Data/Rietkerk/Prelims/Non-Stochastic/" + ast.str() + "-" + est.str() + "_dP_" + dPo.str();
  path_to_folder = prelim_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "/TimeSeries";
  recursive_dir_create(path_to_folder);
	// Creating a string stream instance to store the values of the parameters in the file name.
	// This is done to avoid overwriting of files.
  /** 
	struct stat info1;
	if( stat( foldername2.str().c_str(), &info1 ) != 0 )
	{
		cout << "Cannot access " << foldername2.str() << ". Creating directory." << endl;
		const int dir_err = system(("mkdir " + foldername2.str()).c_str());
		if (-1 == dir_err)
		{
			printf("Error creating directory!n");
			exit(1);
		}
	}
  */
  path_to_folder = "../Data/Rietkerk/Non-Stochastic";
  recursive_dir_create(path_to_folder);
  recursive_dir_create("../Data/Rietkerk/Non-Stochastic/"+ std::to_string(SpB) +"Sp");
	// Creating a string stream instance to store the values of the parameters in the file name.
	// This is done to avoid overwriting of files.

  first_order_critical_exp_delta(div, t_max, a_start, a_end, c, gmax, alpha, d, rW, W0, D, K, dt, dx, dP, r, g);

  //tupac_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, auto &Rh0,
	 //double t_max, double a[], double b[], double c[], double D[], double sigma[], double dt, double dx, int r,  int g)

  return 0;
}