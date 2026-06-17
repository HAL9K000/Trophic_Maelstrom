#include "rietkerk_bjork_basic.h"

int main(int argc, char *argv[])
{
  string preFIX; // Prefix for the output files
  double a_c =1/24.0;
  double c, gmax, alpha, d, rW, W0, t_max, dt, dx; //int g, div;
  double a_start, a_end; double r; double dP; // Kick for high initial state
  int g, div;

  #if defined(INIT) && INIT > 1
    string error_string = "ERROR: The initial condition MUST be set to 1 (Random Speckles) or 0 (Homogeneous).\n";
    cout << error_string; cerr << error_string;
    exit(1);
  #endif

  //c= 10; d =0.25; gmax= 0.05; alpha =0.2; W0= 0.2; rW = 0.2; // From Bonachela et al 2015
  
  // Using units of mass = kg, length = km, time = hr.
  c = 10000; gmax = 0.05*pow(10, -3.0)/24.0; d = 0.25/24.0; alpha =0.2/24.0; W0 = 0.2; rW = 0.2/24.0; // From Bonachela et al 2015

  //double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2

  double k0, k1, k2; double d0, d1, d2; double s0; //Slightly 


  /** 
  dx= 0.2 ; //From Bonachela et al 2015 (in m)
  d0 = 0.001; d1 = 0.001; d2= 0.1; //From Bonachela et al 2015 (in m^2/day)
  k0= 0; k1 = 5; k2 =5;
  s0 = sqrt(0.025); // ~ D/(dx)^2 (in g^{0.5}/(m day))
  */
  /* */
  dx= 0.1 ; //From Bonachela et al 2015 (in km)
  //d0 = pow(10.0, -9.0)/24.0; d1 = pow(10.0, -9.0)/24.0; d2= pow(10.0, -7.0)/24.0; //From Bonachela et al 2015 (in km^2/hr)
  d0 = 0.00025/24.0; d1 = 0.00025/24.0; d2= 0.025/24.0; //From Bonachela et al 2015 (in km^2/hr)
  k0= 0; k1 = 5; k2 =5000;
  s0 = 1*sqrt(d0/(dx*dx)); // ~ D0/(dx)^2 (in kg^{0.5}/(km hr))
  //s1 = 1; // ~ D/(dx)^2 (in kg^{0.5}/(km hr))

  double D[Sp] ={d0, d1, d2}; //Diffusion coefficients for species.
  double K[3] ={k0, k1, k2}; //Diffusion coefficients for species.
  double sigma[Sp] ={s0, 0, 0}; //Demographic stochasticity coefficients for species.
  double v[SpB] ={0}; //Velocity of species.

  double kappa = c*gmax - d;
  double p0istar, Wstar, p0mstar; // Analytic steady state values.

  Wstar = d*K[1]/kappa; 
  //p0istar = (c/d)*(R - rW*Wstar);
  //p0mstar = (R/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);

  double A[SpB][SpB]; double H[SpB][SpB]; 
  double E[SpB] ={1.0}; double M[SpB] ={d}; double pR[SpB] ={0.0};
  

  cout << "Presets are as follows:\n";
  cout << "For the vegetation:\n";
  cout << "c = " << c << "\t gmax = " << gmax << "\t d = " << d << "\t alpha = " << alpha << "\t W0 = " << W0 << "\t rW = " << rW << endl;
  cout << "\nMFT biomass density of Soil Water (W*) = " << Wstar << " kg/km^2\n" << endl;

  cout << "Diffusion Constants For Species:\n";
  cout << "D0 = " << setprecision(16) << D[0] << "\t D1 = " << D[1] << "\t D2 = " << setprecision(16) << D[2] <<  endl;

  cout << "This is a 1Sp Stochastic Rietkerk Vegetation Model Script With UNIT/RANDOM Initial Conditions\n";
  cout << "Header Species Check: " << std::to_string(SpB) << "\n";

  cout << "\n============================================================\n";
  
  if(argc == 1)
  {
    cout << "Please enter all the arguments within the script.\n";
    cout << "Enter desired time-step (~0.25 is a good choice): ";
    cin >> dt;

    cout << "Enter maximum duration of simulation (in hrs): ";
    cin >> t_max;

    cout << "Enter the desired grid-size (L) (dx = 100 m): ";
    cin >> g;

    cout << "Enter the desired number of replicates (r): ";
    cin >> r;

    cout << "Note that a relevant value of 'a' may be construed as:\t" << setprecision(10) << a_c <<endl;

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
  // If there are otherwise an arbitrary number of arguements, terminate the program.
  else
  {
    cout << "Please enter the correct number of arguments.\n";
    cout << "The correct number of arguments is 9.\n";
    cout << "The arguments are as follows: dt, t_max, g, r, a_start, a_end, div, dP, PREFIX.\n";
    exit(1);
  }

  cout << "\n System-level details:\n";
  cout << "dx = " << dx << "\t dt = " << dt << "\t t_max = " << t_max << "\t g = " << g << "\t r = " << r << endl;

  cout << "\n=============================================================\n";

  if(preFIX == "N/A" || preFIX == "n/a" || preFIX == "NA" || preFIX == "na" || preFIX == "N" || preFIX == "n" || preFIX == "")
  {
    cout << "No prefix provided. Using default prefix: " << prefix << endl;
    preFIX = prefix;
  }

  set_Prefix(preFIX); //Set the prefix (and thus save folders) to the user-defined prefix.
  cout << "Save directory for frames: " << frame_folder << "\n";
  cout << "Save directory for preliminary data: " << prelim_folder << "\n";
  cout << "Save directory for final data: " << stat_prefix << "\n";

  set_global_system_params(dt, dx); //Set the global parameters for the simulation.
  cout << "Global parameters set, with dt/2.0 = " << dt2 << " and dx*dx = " << dx2 <<  " and 1/(dx*dx) = " << dx1_2 << "\n";

  // Next, create the integer array dtV[SpB], which is initialised to {0}.
  int dtV[SpB] = {0}; // Vegetation does not advect.
  
  
  // Equations for MFT E Eqilibrium values  as functions of a (Rainfall).
  double mV = c/d; double cV =  -rW*Wstar*c/d; //V* = (c/d)*(R - rW*Wstar);
  double mW = Wstar; double cW = 0;
  double mW_Prev = 1/rW; double cW_Prev = 0.0;
  double mO = 120.4744181; double cO = 0.58243061; // O* = R*(V* + K2)/[(V* + K2*W0)*alpha]
  double mO_Prev = 1/(alpha*W0); double cO_Prev = 0.0; 
  //double mO = 596.37; double cO = 0.0;
  #if defined(INIT) && INIT == 0
    // In the homogenous case, the vegetation is set initially to 1 (i.e. appropriate scaling_factor values)
    string MFT_V = std::to_string(0.0); string MFT_V_Prev = std::to_string(0.0);
    string MFT_W = std::to_string(0.0); string MFT_W_Prev = std::to_string(0.0);
    string MFT_O = std::to_string(0.0); string MFT_O_Prev = std::to_string(0.0);
  #else
    string MFT_V = std::to_string(mV) + " * a + " + std::to_string(cV);
    string MFT_V_Prev = std::to_string(0);
    string MFT_W = std::to_string(mW) + " * a - " + std::to_string(cW);
    string MFT_W_Prev = std::to_string(mW_Prev) + " * a";
    string MFT_O = std::to_string(mO) + " * a + " + std::to_string(cO);
    string MFT_O_Prev = std::to_string(mO_Prev) + " * a";
  #endif
  

  MFT_Vec_CoexExpr.assign({MFT_V_Prev, MFT_W_Prev, MFT_O_Prev, MFT_V, MFT_W, MFT_O});
  
  
  // NOTE: The following clow is used for analytic MFT based initial conditions.
  #if defined(INIT) && INIT == 0
    double scaling_factor[2*Sp] = {1, 1, 1,  1, 1, 1};
    // USED FOR HOMOGENEOUS MFT BASED INITIAL CONDITIONS.
  #else
    double scaling_factor[2*Sp] = {0, 1, 1,  5, 1, 1};
    // USED FOR PERIODIC MFT BASED INITIAL CONDITIONS.
  #endif
  
  // NOTE: The following clow is used for analytic MFT based initial conditions.
  //double scaling_factor[2*Sp] = {0, 1, 1,  5, 1, 1};
  double chigh[2*Sp]; 

  stringstream ast, est, dPo, geq; ast << a_start; est  << a_end; dPo << dP, geq << setprecision(5) << Wstar;

  
	// This is done to avoid overwriting of files.
  string path_to_folder = frame_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Geq_" + geq.str();
  recursive_dir_create(path_to_folder);

  path_to_folder = prelim_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Geq_" + geq.str() +"/TimeSeries";
  recursive_dir_create(path_to_folder);

  recursive_dir_create("../Data/Rietkerk/Stochastic/"+ std::to_string(SpB) +"Sp");

  first_order_critical_exp_delta_stochastic_MultiSp(div, t_max, a_start, a_end, a_c, c, gmax, alpha, rW, W0, D, v, K, sigma, A, H, E, M, pR, dtV, scaling_factor, dt, dx, dP, r, g, Wstar);
  
  //first_order_critical_exp_delta_stochastic(div, t_max, a_start, a_end, c, gmax, alpha, d, rW, W0, D, K, sigma, dt, dx, dP, r, g);
  return 0;
}