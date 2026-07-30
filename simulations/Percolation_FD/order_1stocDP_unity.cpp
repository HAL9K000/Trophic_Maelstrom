#include "MultiSPDP.h"

int main(int argc, char* argv[])
{
  increase_stack_limit(1024L); //Increase stack limit to 1024 MB.

  string preFIX; // Prefix for the output files
  //double a_c =1/24.0;
  double a_c, b, c, gmax, alpha, d, rW, W0, t_max, dt, dx; //int g, div;
  double a_start, a_end; double r; double dP; // Kick for high initial state
  int g, div;

  //c= 10; d =0.25; gmax= 0.05; alpha =0.2; W0= 0.2; rW = 0.2; // From Bonachela et al 2015
  /**
  dx= 0.2 ; //From Bonachela et al 2015 (in m)
  d0 = 0.001; d1=0; d2 = 0.001; d3= 0.1; //From Bonachela et al 2015 (in m^2/day)
  k0= 0; k1 = 5; k2 =5;
  s0 = sqrt(0.025); // ~ D/(dx)^2 (in g^{0.5}/(m day))
  */

  // Using units of mass = kg, length = km, time = hr.
  // ASSUMING MASS OF  GRAZER = 20 kg, MASS OF PREDATOR = 100 kg.
  b = pow(10.0, -6.00); // ~ aij in km^2/(hr kg)

  //double p0i = 0.5; double p0j= p0i/200; double p0j= 2.25; double p0m= 8; // In g/m^2

  double k0, k1, k2; double d0, d1, d2, d3, d4; double s0;

  dx = 0.01; //From Bonachela et al 2015 (in km)
  d0 = 0.00025 / 24.0; d1 = 0.0298; d2 = 0.05221; d3 = 0.00025 / 24.0; d4 = 0.025 / 24.0; //From Bonachela et al 2015 (in km^2/hr)
  s0 = 1 * sqrt(d0 / (dx * dx)); // ~ D0/(dx)^2 (in kg^{0.5}/(km hr))

  // Allometric scaling for velocity: 
  // v (m/hr) = 43.706*M^1.772*(1 - e^{-14.27*(M)^(-1.865)}), where M is mass in kg.

  double D[Sp] = { d0 }; //Diffusion coefficients for species.
  double sigma[Sp] = { s0 }; //Demographic stochasticity coefficients for species.
  double v[Sp] = { 0 }; //Velocity of species.

  /**

  double beta = c*gmax - d;
  double p0istar, p0jstar, p0mstar; // Analytic steady state values.

  p0jstar = d*K[1]/beta;
  p0istar = (c/d)*(R - rW*p0jstar);
  p0mstar = (R/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
  */
  double aij, hij, ej, mj; double ajm, hjm, em, mm; // Parameters for the grazer and predator.

  aij = 3.6 * pow(10.0, -6.08) * pow(20.0, -0.37); // in km^2/(hr kg)
  hij = 1; // Handling time in hrs
  ej = 0.45; mj = 0.061609 * pow(20.0, -0.25) / 8760.0; // Mortality rate in hr^{-1}



  double Km = pow(10.0, 1.22) * pow(100.0, -0.31); // Carrying capacity of predator in kg/km^2s

  //double Gstar = mm/((em -mm*hjm)*ajm); //Steady state for grazer.
  double H[SpB][SpB] = { {0} };    // Handling time matrix [T]. No canabilism, symmetric effects.
  double A[SpB][SpB] = { {0} };  // Attack rate matrix []. No canabilism, symmetric effects. 
  double M[SpB] = { d };     // Mortality Rate of Species. (in hr^{-1})
  double E[SpB] = { 1.0 }; //Efficiency of consumption.
  //double D[Sp] ={d0, d1}; //Diffusion coefficients for species. (in km^2/hr)
  double pR[Sp] = { 0.0 }; //Perception radius of species (in km)

  double Vstar = 0;
  double init_frac_grazpred = 1.0; //Initial fraction of grazers and predators in the high state.

  cout << "Presets are as follows:\n";
  cout << "For the vegetation:\n";
  cout << "c = " << c << "\t b = " << b << endl;
  cout << "\n System-level details:\n";
  cout << "pi = " << PI << endl;

  cout << "Diffusion Constants For Species:\n";
  cout << "D0 = " << setprecision(16) << D[0] << endl;

  cout << "This is a 1Sp (One) Stochastic DP Model Script WITHOUT DDM\n";
  cout << "Header Species Check: " << std::to_string(SpB) << "\n";

  cout << "\n============================================================\n";

  // There are eight user-defined parameters, as follow below: dt, t_max, g, r, a_start, a_end, div, dP.

  //Check if the user has entered the correct number of arguments. If not prompt the user to enter all the arguments within the script.
  // First check if no arguments are entered. In this case, the user will be prompted to enter all the arguments using cin.

  if (argc == 1)
  {
    cout << "Please enter all the arguments within the script.\n";
    cout << "Enter desired time-step (~0.1 hr is a good choice): ";
    cin >> dt;

    cout << "Enter maximum duration of simulation (in hrs): ";
    cin >> t_max;

    cout << "Enter the desired grid-size (L) (dx = 100 m): ";
    cin >> g;

    cout << "Enter the desired number of replicates (r): ";
    cin >> r;

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
  else if (argc == 10)
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
    cout << "The arguments are as follows: dt, t_max, g, r,a_start, a_end, div, dP, PREFIX.\n";
    exit(1);
  }

  cout << "\n System-level details:\n";
  cout << "dx = " << dx << "\t dt = " << dt << "\t t_max = " << t_max << "\t g = " << g << "\t r = " << r << endl;

  cout << "\n=============================================================\n";


  if (preFIX == "N/A" || preFIX == "n/a" || preFIX == "NA" || preFIX == "na" || preFIX == "N" || preFIX == "n" || preFIX == "")
  {
    cout << "No prefix provided. Using default prefix: " << prefix << endl;
    preFIX = prefix;
  }
  //Defining the save folders.

  set_Prefix(preFIX); //Set the prefix (and thus save folders) to the user-defined prefix.
  cout << "Save directory for frames: " << frame_folder << "\n";
  cout << "Save directory for preliminary data: " << prelim_folder << "\n";
  cout << "Save directory for final data: " << stat_prefix << "\n";

  set_global_system_params(dt, dx); //Set the global parameters for the simulation.
  cout << "Global parameters set, with dt/2.0 = " << dt2 << " and dx*dx = " << dx2 << " and 1/(dx*dx) = " << dx1_2 << "\n";
  int dtV[SpB] = { 0 };
  //INITIAL CONDITIONS:

  // Equations for MFT E Eqilibrium values  as functions of a (Rainfall).


  //double mG = E[1]/(kappa); double cG = E[1]*b*M[1]/(kappa*kappa); 
  // Coexistance equilibrium density of grazer as a function of a. (Geq = mG*a + cG)
  double mV = 1 / b; // Vegetation only equilibrium density as a function of a. (Veq = a/b =mV*a)

  string MFT_PreV = std::to_string(0.0);
  string MFT_V = std::to_string(0.0);


  MFT_Vec_CoexExpr.assign({ MFT_PreV, MFT_V });
  // NOTE: The following clow is used for analytic MFT based initial conditions. 
  double scaling_factor[2 * Sp] = { 1, 1 };

  // NORMAL INIT CLOSE TO MFT
  double chigh[2 * Sp] = { dP, 10000.0 + dP };

  //add_two(3.14, 1.78);
  cout << maxis(3.14, 1.78) << endl;
  int go = 8;
  add_three(1, 2, go);

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

  stringstream ast, est, dPo, veq; ast << a_start; est << a_end; dPo << dP; veq << setprecision(5) << Vstar;


  // This is done to avoid overwriting of files.
  string path_to_folder = frame_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Veq_" + veq.str();
  recursive_dir_create(path_to_folder);

  path_to_folder = prelim_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str()
    + "_Veq_" + veq.str() + "/TimeSeries";
  recursive_dir_create(path_to_folder);


  recursive_dir_create("../Data/DP/Stochastic/" + std::to_string(SpB) + "Sp");

  first_order_critical_exp_delta_stochastic_MultiSp(div, t_max, a_start, a_end, a_c, b, c, D, v, sigma, A, H, E, M, pR, dtV, scaling_factor, dt, dx, dP, r, g, -1.0, Vstar);


  return 0;
}