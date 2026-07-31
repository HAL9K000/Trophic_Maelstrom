#include "SPDP.h"
#include "Debug.h"

//--------------------Defining Primitives------------------------------------//

void increase_stack_limit(long stack_size){

	//Source: https://stackoverflow.com/questions/2279052/increase-stack-size-in-linux-with-setrlimit
	//Credit: https://stackoverflow.com/users/253056/paul-r
	//Used with modification

	// Necessary because of complexity of integrations.

	const rlim_t kStackSize = stack_size * 1024L * 1024L;   // min stack size = 128 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }
}

bool maxis(int a, int b) 
{
	//Binary function necessary for determining max in a list.
 	return (a < b); 
}

void init_fullframe(double array[], int size)
{
	//Returns array with all elements initialised to 1
	for(int i=0; i< size; i++)
	{	array[i] =1 ;}
}

void init_constframe(double array[], double constant, int size)
{
	//Returns array with all elements initialised to a constant value: "const"
	for(int i=0; i< size; i++)
	{	array[i] =constant ;}
}

void init_randframe(double array[], double mean, double sd, int size)
{
	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);
	//Returns array with all elements initialised to random values drawn from a Gaussian distribution with mean mu and SD sigma,
	// with all values 
	for(int i=0; i< size; i++)
	{	
		normal_distribution<double> norm;
		norm = normal_distribution<double>(mean, sd);
		array[i] =norm(rng);
		if(array[i] < 0)
			array[i] = 0.0; //Values less than 0 are unrealistic and not permissible.
	}
}

void init_quarterframe(double array[], double c1, double c2, double c3, double c4, int L)
{
	//Returns array with elements initialised to a constant values: "const1", "const2", "const3", "const4" in clockwise quadrants
	// starting from the top left.
	int l2 = static_cast<int>(L/2);
	for(int i=0; i< l2; i++)
	{
		int il2 = i + l2;	
		for(int j=0; j< l2; j++)
		{
			array[i*L + j] =c1; //Filling up top-left quadrant
			array[i*L + l2 + j] =c2; //Filling up top-right quadrant
			array[il2*L + l2 + j] =c3; //Filling up bottom-right quadrant
			array[il2*L  + j] =c4; //Filling up bottom-right quadrant
		}
			
		
	}
}

void init_solitarytear(double array[], int length)
{
	//Returns array with only central element initalised to 1, all others 0.
	for(int i=0; i< length*length; i++)
		array[i] =0;
	int lhalf = static_cast<int>(length/2.0);
	array[length*lhalf +lhalf] =1;
}

double mean_of_array(double array[],int size){

	double sum = 0.0;

	for (int i =0; i<size; ++i){
		sum += array[i];
	}
	return sum/(double)size;
}

double standard_deviation_of_array(double array[],int size){

	double mean = mean_of_array(array,size);
	double sum = 0.0;

	for (int i = 0; i < size; ++i){
		sum += (array[i]-mean)*(array[i]-mean);
	}

	double variance = sum/(double)size;

	return sqrt(variance);
}

double mean_of_vector(vector<double> array,int size){

	double sum = 0.0;

	for (int i =0; i<size; ++i){
		sum += array[i];
	}
	return sum/(double)size;
}

double standard_deviation_of_vector(vector<double> array,int size){

  double mean = mean_of_vector(array,size);
	float sum = 0.0;

	for (int i =0; i<size; ++i){
		sum += (array[i]-mean)*(array[i]-mean);
	}
  double variance = sum/(double)size;
	return sqrt(variance);
}

double occupied_sites_of_vector(vector<double> array,int size){
//Calculates number of occupied sites in a vector.
	double sum = 0.0;

	for (int i =0; i<size; ++i){
		if (array[i] > 0)
			sum+=1;
	}
	return sum;
}

auto meansq_spread_of_vector(vector<double> array, int g, int c_x, int c_y){
	//Calculates number of occupied sites in a vector. c_x and c_y represent the Cartesian coordinates of the central site.
	struct coord {        // Declare a local structure 
    	double i1; int i2; };
	double sq_dist = 0.0; int sum=0;
	for (int i =0; i<g*g; ++i){
		if (array[i] > 0)
			{	int p = int(i/g); int q = i%g; //Indices of occupied site.
				sum+=1;
				sq_dist += (p-c_x)*(p-c_x) + (q-c_y)*(q-c_y); //Calcuates norm-2 distance.
			}
	}
	return coord {sq_dist, sum};
}



void var_mean_incremental_surv_runs(double t_avg_var_rep_N[][5], double X_curr[][2], int size)
{
	/** X_curr contains time-series data (indexed as (N(t), <rho(t)>x)) from the current surviving replicate (given by r) with a length "size".
	 *  rep_avg_var stores the time(t), running density average and variance of X_curr, as well as N(t) and number of survivng replicates 
	 * for each time-step from previous surviving replicates ( <= r-1) respectively.
	 * In this function, X_curr is used to update the running avg and variance of rep_avg_var to account for current surviving replicate r(t).
	 * Note rep_avg_Var[t][3] stores number of surviving runs at t.
	*/
	double mean_prev[size];
	for(int i=0; i<size; i++)
	{
		int r= t_avg_var_rep_N[i][3]+1; //Number of surviving replicates
		if(X_curr[i][0] == 0.0)
		{	//Non-surviving run. No-updates required. Also <rho(t)>x == 0, it remains 0 for t' > t.
			break;
		}
		else
		{	//Number of surviving run at time t.
			mean_prev[i] = t_avg_var_rep_N[i][1]; //Stores X_n-1 mean data as temp variable.
			if(t_avg_var_rep_N[i][3] == 0.0)
			{	//First surviving run encountered at time t, encountered here.
				t_avg_var_rep_N[i][1] = X_curr[i][1]; t_avg_var_rep_N[i][2] = 0.0; t_avg_var_rep_N[i][4] = X_curr[i][0];
			}
			else
			{	//Not the first surviving run, which means r = rep_avg_var[i][3]+1 >=2
				t_avg_var_rep_N[i][1] = (X_curr[i][1] + (r-1)*t_avg_var_rep_N[i][1])/r; //Note r>= 1.
				t_avg_var_rep_N[i][4] = (X_curr[i][0] + (r-1)*t_avg_var_rep_N[i][4])/r; //Note r>= 1.
				//Calculating and storing incremental variance (V^2_n) [For formula derivation, refer to notes.]
				t_avg_var_rep_N[i][2] = (t_avg_var_rep_N[i][1] - mean_prev[i])*(t_avg_var_rep_N[i][1] - mean_prev[i]) +
				(1/(r-1))*( (r-2)*t_avg_var_rep_N[i][2] + (X_curr[i][1] - t_avg_var_rep_N[i][1])*(X_curr[i][1] - t_avg_var_rep_N[i][1]));
			}
			t_avg_var_rep_N[i][3] = r; //At the end, update number of succesful runs.
		}
	}

}

void var_mean_incremental_all_runs(double rep_avg_var[][3], double X_curr[][2], int size, int r)
{
	/** X_curr contains time-series data (indexed as (t, data)) from the current replicate (given by r) with a length "size".
	 *  rep_avg_var stores the running average and variance of X_curr for each time-step from previous replicates ( <= r-1).
	 * In this function, X_curr is used to update the running avg and variance of rep_avg_var to account for current surviving r(t).
	*/
	double mean_prev[size];
	for(int i=0; i<size; i++)
	{
		rep_avg_var[i][0] = X_curr[i][0]; //Store time points.
		mean_prev[i] = rep_avg_var[i][1]; //Stores X_n-1 mean data as temp variable.
		//Next calculate and store incremental mean (X_n).
		rep_avg_var[i][1] = (X_curr[i][1] + (r-1)*rep_avg_var[i][1])/r; //Note r>= 2.
		//Next calculate and store incremental variance (V^2_n) [For formula derivation, refer to notes.]
		rep_avg_var[i][2] = (rep_avg_var[i][1] - mean_prev[i])*(rep_avg_var[i][1] - mean_prev[i]) +
			(1/(r-1))*( (r-2)*rep_avg_var[i][2] + (X_curr[i][1] - rep_avg_var[i][1])*(X_curr[i][1] - rep_avg_var[i][1]));
	}

}

void spreading_incremental_surv_runs(double t_N_R2_rep[][4], double X_curr[][3], int size)
{
	//ONLY APPLICABLE IN THE CASE WHEN UPDATING ACROSS REPLICATE MEANS FOR SPREADING EXPERIMENTS.
	
	/** X_curr contains time-series data (indexed as (N(t), MSD, # Surviving Runs)) from the current surviving replicate (given by r) with a length "size".
	 *  rep_avg_var stores the time(t), running density average and variance of X_curr, as well as N(t) and number of survivng replicates 
	 * Note rep_avg_Var[t][3] stores number of surviving runs at t.
	*/
	for(int i=0; i<size; i++)
	{
		int r= t_N_R2_rep[i][3]+1; //Number of surviving replicates
		if(X_curr[i][0] == 0.0)
		{	//Non-surviving run. No-updates required. Also <rho(t)>x == 0, it remains 0 for t' > t.
			break;
		}
		else
		{	//Number of surviving run at time t.
			t_N_R2_rep[i][1] = (X_curr[i][0] + (r-1)*t_N_R2_rep[i][1])/r; //Note r>= 1.
			t_N_R2_rep[i][2] = (X_curr[i][1] + (r-1)*t_N_R2_rep[i][2])/r; //Note r>= 1.
			t_N_R2_rep[i][3] = r; //At the end, update number of succesful runs.
		}
	}
}

template<typename T>std::vector<double> linspace(T start_in, T end_in, int num_in)
{
	//Equivalent of numpy's linspace method.

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end);
  return linspaced;
}

template<typename T>std::vector<double> lnspace(T start_in, T end_in, int log_points)
{
       double logarithmicBase = exp(1);
	   std::vector<double> logspaced;
       double logMin = log(static_cast<double>(start_in));
       double logMax = log(static_cast<double>(end_in));
       double delta = (logMax - logMin) / (log_points-1);
	   if (log_points == 0) { return logspaced; }
  	   if (log_points == 1)
    	{
      		logspaced.push_back(static_cast<double>(start_in));
      		return logspaced;
    	}
       for (int i = 0; i < log_points-1; ++i)
       {
           double v = pow(logarithmicBase, logMin + delta*i);
		   logspaced.push_back(v);
       }
	   logspaced.push_back(static_cast<double>(end_in));
  	   return logspaced;
}

std::vector<double> logarithmic_time_bins(double t_max, double dt)
{	// Computes and returns ln-distributed points from t= e^2 to ln(t_max) (the latter rounded down to 1 decimal place)
	std::vector<double> logspaced;
	if(t_max <= exp(2.0))
	{
		cout << "GOUT!!"<<endl;
		return logspaced; // No points recorded if t_max < e^{2}
	}

	double tend_one = static_cast<double>(floor(log(t_max)* 10.)) / 10.0;
	// Finds the closest power of e (rounded down to one decimal place) from t_max.
	double pow =2.0;
	while(pow <= tend_one)
	{
		double t_raw = exp(pow); // Getting raw time-stamp.
		//Rounding down raw time-stamp to dt.
		double t_fin = static_cast<int>(t_raw/dt)*dt;

		if(t_fin < exp(2))
			t_fin = t_fin+dt; //Ensure t_fin is strictly above e^{2}

		if(logspaced.size() > 0  && t_fin != logspaced[logspaced.size() -1])
		{ 	logspaced.push_back(t_fin);	} // This is to prevent duplicate values from arising in logspaced.
		else if(logspaced.size() == 0)
			logspaced.push_back(t_fin);
		pow += 0.04; //25 time-stamps per decade.
	}
	cout << "TROUT! " << logspaced.size() <<endl;
	return logspaced;


}

std::vector<double> logarithm10_time_bins(double t_max, double dt)
{	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place)
	std::vector<double> logspaced;
	if(t_max <= 1)
	{
		cout << "GOUT!!"<<endl;
		return logspaced; // No points recorded if t_max < e^{2}
	}

	double tend_one = static_cast<double>(floor(log10(t_max)* 10.)) / 10.0;
	// Finds the closest power of e (rounded down to one decimal place) from t_max.
	double power =0.0;
	while(power <= tend_one)
	{
		double t_raw = pow(10, power); // Getting raw time-stamp.
		//Rounding down raw time-stamp to dt.
		double t_fin = static_cast<int>(t_raw/dt)*dt;

		if(logspaced.size() == 1)
				cout << "Feq " << logspaced[0] << " HODOR " << t_fin << " modor " << t_fin - logspaced[0] << endl;

		if(t_fin <= 1.0)
			t_fin = t_fin+dt; //Ensure t_fin is strictly above 1.0

		if(logspaced.size() > 0  && t_fin != logspaced[logspaced.size() -1])
		{ 
			logspaced.push_back(t_fin);
		} // This is to prevent duplicate values from arising in logspaced.
		else if(logspaced.size() == 0)
			logspaced.push_back(t_fin);
		power += 0.04; //25 time-stamps per decade.
	}
	int n = logspaced.size();
	for(int i=n-1; i >0; i--)
	{	//Rooting out duplicates.
		if(logspaced[i] < logspaced[i-1] +dt/2.0)
			logspaced.erase(logspaced.begin() + i);
	}
	cout << "TROUT! " << logspaced.size() <<endl;
	return logspaced;


}

int theborderwall(vector<double> &Rho_t, int g)
{
	//Checks if the border sites are occupied.
	int g2 = g*(g-1);
	for(int i =0; i< g; i++)
	{
		if(Rho_t[i] > 0 || Rho_t[g*i] > 0  || Rho_t[g2 + i] > 0 || Rho_t[g*i + g - 1 ] > 0)
			return 1; //Top, left, bottom and right borders are occupied by active sites.
	}
	return 0; //No active sites detected at the edges.
}

//----------------------------- Stochastic SPDP RK4 Integration Machinery --------------------------------------------//

void determine_neighbours_R2(int neighbours_R2[][2][5], int g)
{
	//Stores neighbouring sites (Van-Neumann Radius of 2) of each cell site in a 
	// g^2*2*5 matrix.
	int a=0; int b= 0; //Useful for accounting purposes
	for (int i=0; i<g*g; i++)
	{
		int p = int(i/g); int q = i%g;
		//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
		for(int j=0; j < 5; j++)
		{
			int k= j-2;
			a= ((p+k)%g + g)%g; // Eq to (p+k)mod g. -2 <= k <= 2
			b= ((q+k)%g + g)%g; // Eq to (q+k)mod g. -2 <= k <= 2

			//Recall, first row stores vertical neighbours, second row stores horizontal neighbours.
			neighbours_R2[i][0][j] = a*g + q; neighbours_R2[i][1][j] = p*g + b;
		}
	}

}

void determine_neighbours_S8(int neighbours_R2[][3][3], int g)
{
	//Stores neighbouring sites in a square (Norm 1) radius of 1
	int a=0; int b= 0; //Useful for accounting purposes
	for (int i=0; i<g*g; i++)
	{
		int p = int(i/g); int q = i%g;
		//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
		for(int y=0; y < 3; y++)
		{	int ky= y-1; a= ((p+ky)%g + g)%g; // Eq to (p+ky)mod g. -1 <= ky <= 1
			for(int x=0; x < 3; x++)
			{
				int kx = x-1;
				b= ((q+kx)%g + g)%g; // Eq to (q+kx)mod g. -1 <= kx <= 1
				//Recall, first row stores vertical neighbours, second row stores horizontal neighbours.
				neighbours_R2[i][y][x] = a*g + b; 
			}	
		}
	}

}

int twilight_zoneS8(vector<double> &Rho_t, int neighbours_S8[][3][3], int i)
{
	//Checks if Rho_dt[i] is 0 and is surrounded by empty patches (radius of 1). Returns 0 if so, else 1.
	//double epsi = 0.000000000001; //Define Epsilon as 10^{-12}
	double sum=0.0;
	for(int x =0; x<3; x++)
	{
		for(int y=0; y<3; y++) //Recall, S8 resembles Laplacian schematic.
		{		
			if(Rho_t[neighbours_S8[i][x][y]] > 0.0) 
				return 1;
			sum += Rho_t[neighbours_S8[i][x][y]];
		}
	}
	if(sum > 0.0)
		return 1;
	return 0; //Only occurs all neighbours in a radius of 1 are empty.

}

int twilight_zoneR2(vector<double> &Rho_t, int neighbours_R2[][2][5], int i)
{
	//Checks if Rho_dt[i] is 0 and is surrounded by empty patches (radius of 2). Returns 0 if so, else 1.

	for(int j =0; j<5; j++)
	{
		//Recall, first row stores horizontal neighbours, second row stores vertical neighbours.
		if(Rho_t[neighbours_R2[i][0][j]] != 0)
			return 1;
		else if(Rho_t[neighbours_R2[i][1][j]] != 0)
			return 1;
	}
	return 0; //Only occurs all neighbours in a radius of 2 are empty.

}



void f_2D(double f[], vector<double> &Rho_t, int n[][2][5], double p, double D, double t, double dt, double dx, int g)
{
  //Vector function that evaluates drho/dt = f(p,t) at given p and t.

	for(int i=0; i < g*g; i++)
	{
		f[i] = -Rho_t[i]*Rho_t[i] - D*(1/(12.0*dx*dx))*(1*(Rho_t[n[i][0][0]] + Rho_t[n[i][0][4]] 
		+ Rho_t[n[i][1][0]] + Rho_t[n[i][1][4]]) -16*(Rho_t[n[i][0][1]] + Rho_t[n[i][0][3]] 
		+ Rho_t[n[i][1][1]] + Rho_t[n[i][1][3]]) +60*Rho_t[n[i][0][2]]);

		// RK4 update of Non-linear parts plus diffusion term
		// 4th order Laplacian from 2SHOC.
	}
}

void fDo_2D(double f[], vector<double> &Rho_t, int n[][2][5], double b, double t, double dt, int g)
{
  //RK4 Integration For Old Dornic Methode, where, as you may recall the Diffusion term is treated analytically.

	for(int i=0; i < g*g; i++)
	{
		f[i] = -b*Rho_t[i]*Rho_t[i];
		// RK4 update of Non-linear parts
	}
}

void eulerf_2D(double f[], vector<double> &Rho_t, int n[][2][5], double p, double D, double t, double dt, double dx, int g)
{
  //Performs Euler update of non-linear and diffusion terms and stores it in f
	for(int i=0; i < g*g; i++)
	{
		f[i] = Rho_t[i] + (-Rho_t[i]*Rho_t[i] - D*(1/(12.0*dx*dx))*(1*(Rho_t[n[i][0][0]] + Rho_t[n[i][0][4]] 
		+ Rho_t[n[i][1][0]] + Rho_t[n[i][1][4]]) -16*(Rho_t[n[i][0][1]] + Rho_t[n[i][0][3]] 
		+ Rho_t[n[i][1][1]] + Rho_t[n[i][1][3]]) +60*Rho_t[n[i][0][2]]))*dt;

		// RK4 update of Non-linear parts plus diffusion term
		// 4th order Laplacian from 2SHOC.
	}
}

void fP_2D(double f[], vector<double> &Rho_t, int n[][2][5], double a, double b, double D, double t, double dt, double dx, int g)
{
  //Vector function that evaluates drho/dt = f(p,t) at given p and t.

	for(int i=0; i < g*g; i++)
	{
		f[i] = a*Rho_t[i] -b*Rho_t[i]*Rho_t[i] - D*(1/(12.0*dx*dx))*(1*(Rho_t[n[i][0][0]] + Rho_t[n[i][0][4]] 
		+ Rho_t[n[i][1][0]] + Rho_t[n[i][1][4]]) -16*(Rho_t[n[i][0][1]] + Rho_t[n[i][0][3]] 
		+ Rho_t[n[i][1][1]] + Rho_t[n[i][1][3]]) +60*Rho_t[n[i][0][2]]);

		// RK4 update of Non-linear parts plus diffusion term
		// 4th order Laplacian from 2SHOC.
	}
}


void find_neighbours_R2(double neighbours_R2[2][5], vector<double> &Rho_t, int p, int q, int g)
{
	//Stores neighbour density in a ball of radius 2 around p,q.
	int a=0; int b= 0; //Useful for accounting purposes

	//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
	for(int j=0; j < 5; j++)
	{
			int k= j-2;
			a= ((p+k)%g + g)%g; // Eq to (p+k)mod g. -2 <= k <= 2
			b= ((q+k)%g + g)%g; // Eq to (q+k)mod g. -2 <= k <= 2

			//Recall, first row stores horizontal neighbours, second row stores vertical neighbours.
			neighbours_R2[0][j] = Rho_t[a*g + q]; neighbours_R2[1][j] = Rho_t[p*g + b];
	}
}

void percNoDiff_2D(vector<vector<double>> &Rho, vector<double> &t_stop, double Rh0[],
	 double t_max, double a, double b, double D, double dt, double dx, int r,  int g)
{
	int n_R2[g*g][2][5] ={0};
	//n_R2 stores neighbours of all sites in a ball of radius 2.
	//Assigning neigbours.
	determine_neighbours_R2(n_R2, g);
	
	for(int j=0; j< 1; j++)
	{
		// Iterating over replicates.
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		thread_local std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd);

		vector <double> Rho_dt(g*g, 0.0); //Kept constant, only updated at end of turn
		double K1[g*g] ={0.0}; double K2[g*g] ={0.0}; double K3[g*g] ={0.0}; double K4[g*g] ={0.0};	
		//Initialise RK4 terms for future use.
		vector <double> Rho_M(g*g, 0.0); //Stores intermediate values of Rho at dt/2 and dt for RK4 calculations.
		cout<< "Initially: "<<endl;
	  for(int i=0; i< g*g; i++)
	  {
				Rho_dt[i] = Rh0[i]; //Assigning initial conditions
				cout<< Rho_dt[i] << " ";
	  } cout<<endl;

	   	stringstream m2;     //To make cout thread-safe as well as non-garbled due to race conditions.
     	m2 << "POCKET!" <<endl; cout << m2.str(); 

		//int t_max = static_cast<int>(t_stop[t_stop.size() -1]); //Storing last element of the same here.
		double t=0; //Initialise t
		int s=0;
		
		while (t < t_max + dt)
		{
			//End at the end of all time.
			//cout << "SOCKET" << endl;

			stringstream m3;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m3 << "B4RK4" <<endl; cout << m3.str();

			fP_2D(K1, Rho_dt, n_R2, a, b, D, t, dt, dx, g); //K1 updated. 
			cout<< "K1" << endl;
			for(int i = 0; i<g*g; dx, i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K1[i];  
			  //cout << K1[i] << " ";
			} //Recall Runge Kutta-Update Schema for K2 = f(t+dt/2, y+dt/2*K1)
			cout << endl;
			fP_2D(K2, Rho_M, n_R2, a, b, D, t +dt/2.0, dt, dx, g); //K2 updated.
			cout<< "K2" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K2[i]; cout << K2[i] << " ";  } //Recall Runge Kutta-Update Schema for K3
		  	cout << endl;
			fP_2D(K3, Rho_M, n_R2, a, b, D, t +dt/2.0, dt, dx, g); //K3 updated.
			cout<< "K3" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt)*K3[i]; cout << K3[i] << " ";  } //Recall Runge Kutta-Update Schema for K4
			cout << endl;
			fP_2D(K4, Rho_M, n_R2, a, b, D, t + dt, dt, dx, g); //K4 updated.
			cout<< "K4" << endl;
			for(int i = 0; i<g*g; i++)
			{
				cout << K4[i] << " ";
			}
			cout << endl;
			 stringstream m5;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m5 << "A7RK4" <<endl; cout << m5.str();
			//Now that K1, K2, K3, K4 updated for all sites, do the stochastic update.
			for(int i = 0; i<g*g; i++)
			{
			Rho_dt[i]+= (dt/6.0)*( K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);
			}
			//Impproved Dornic  Integration of Linear & Stochastic Term
			//CREDITS: PAULA VILLA MARTIN.
			t+=dt; //Update timestep
			cout << "MOULIN ROUGE Final Rho(t):" <<endl;
			for(int i=0;i<Rho_dt.size();i++)
			{ cout<< Rho_dt[i] << " "; } cout<<endl;

			

			stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m4 << "SOCKET!" <<endl; cout << m4.str();
			

			//Updating Pijm at end of sweep.

			//cout << t << endl;
			if(t > t_stop[s] -dt/2.0 && t <= t_stop[s] + dt/2.0)
			{
				stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
      			m4 << "LOCKET!:\t" << t << " " << a  <<endl; cout << m4.str();
				for ( int k = 0; k< Rho_dt.size(); k++)
			  {	// Recall Rho = | p (Order Parameter) | Replicate (r) | Time(t)  | a*L + b | Cell Density[a][b] |
					Rho.push_back({a, static_cast<double>(j), t, static_cast<double>(k), Rho_dt[k]});
			  }
				s = (s+1)%t_stop.size(); // To prevent out of bounds exception.
			}

		} // End of t loops



		vector<double>().swap(Rho_M); vector<double>().swap(Rho_dt);
		// Remove dynamic vectors from memory and free up allocated space.

		

	}



}


void percolation_2D(vector<vector<double>> &Rho, vector<double> &t_stop, double Rh0[],
	 double t_max, double p, double D, double sigma, double dt, double dx, int r,  int g)
{
	int n_roll =10000; int n_star = 1000; int n_R2[g*g][2][5] ={0};
	//n_R2 stores neighbours of all sites in a ball of radius 2.
	//Assigning neigbours.
	determine_neighbours_R2(n_R2, g);

	//For debugging.
	int i1 =0 ; int i2 = int(g*g/2); int i3 = g -1;
	vector <vector <double>> KRi1; vector <vector <double>> KRi2;
	vector <vector <double>> KRi3; double diff[2][g*g]; double Rhog[g*g];
	//Stores K1, K2, K3, K4 plus Rho_dt and Diffusion details for each index.

	int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j <<endl; cout << m1.str();
		
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		//thread_local std::default_random_engine rng(rd); //engine
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number


		/** Recall in Mean Field Limit (MFT), DP is given by:
		d(rho)/dt = (2p -1)rho - 1*rho^(2) **/
		//Defining the Dornic variables
		double b = 2*p -1.0;
		double lambda_exp = exp( (b)*dt); double lambda = 2*(b)/(sigma*sigma*(lambda_exp -1.0));

		vector <double> Rho_dt(g*g, 0.0); //Kept constant, only updated at end of turn
		double K1[g*g] ={0.0}; double K2[g*g] ={0.0}; double K3[g*g] ={0.0}; double K4[g*g] ={0.0};	
		//Initialise RK4 terms for future use.
		vector <double> Rho_M(g*g, 0.0); //Stores intermediate values of Rho at dt/2 and dt for RK4 calculations.
		cout<< "Initially: "<<endl;
	  for(int i=0; i< g*g; i++)
	  {
				Rho_dt[i] = Rh0[i]; //Assigning initial conditions
				cout<< Rho_dt[i] << " ";
	  } cout<<endl;

	   	stringstream m2;     //To make cout thread-safe as well as non-garbled due to race conditions.
     	m2 << "POCKET!" <<endl; cout << m2.str(); 
		double t=0; //Initialise t
		int s=0;
		
		while (t < t_max + dt)
		{
			//End at the end of all time.
			//cout << "SOCKET" << endl;

			stringstream m3;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m3 << "B4RK4" <<endl; cout << m3.str();

			f_2D(K1, Rho_dt, n_R2, p, D, t, dt, dx, g); //K1 updated. 
			cout<< "K1" << endl;
			for(int i = 0; i<g*g; dx, i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K1[i];  
			  cout << K1[i] << " ";
			} //Recall Runge Kutta-Update Schema for K2 = f(t+dt/2, y+dt/2*K1)
			cout << endl;
			f_2D(K2, Rho_M, n_R2, p, D, t +dt/2.0, dt, dx, g); //K2 updated.
			cout<< "K2" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K2[i]; cout << K2[i] << " ";  } //Recall Runge Kutta-Update Schema for K3
		  	cout << endl;
			f_2D(K3, Rho_M, n_R2, p, D, t +dt/2.0, dt, dx, g); //K3 updated.
			cout<< "K3" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt)*K3[i]; cout << K3[i] << " ";  } //Recall Runge Kutta-Update Schema for K4
			cout << endl;
			f_2D(K4, Rho_M, n_R2, p, D, t + dt, dt, dx, g); //K4 updated.
			cout<< "K4" << endl;
			/**for(int i = 0; i<g*g; i++)
			{
				cout << K4[i] << " ";
			}*/
			cout << endl;

			//Debugging section.
			for(int i = 0; i<g*g; i++)
			{
				diff[0][i] = - D*(1/(12.0*dx*dx))*(1*(Rho_dt[n_R2[i][0][0]] + Rho_dt[n_R2[i][0][4]] 
				+ Rho_dt[n_R2[i][1][0]] + Rho_dt[n_R2[i][1][4]]) -16*(Rho_dt[n_R2[i][0][1]] + Rho_dt[n_R2[i][0][3]] 
				+ Rho_dt[n_R2[i][1][1]] + Rho_dt[n_R2[i][1][3]]) +60*Rho_dt[n_R2[i][0][2]]);
				diff[1][i] = - D*(1/(12.0*dx*dx))*(1*(Rho_M[n_R2[i][0][0]] + Rho_M[n_R2[i][0][4]] 
				+ Rho_M[n_R2[i][1][0]] + Rho_M[n_R2[i][1][4]]) -16*(Rho_M[n_R2[i][0][1]] + Rho_M[n_R2[i][0][3]] 
				+ Rho_M[n_R2[i][1][1]] + Rho_M[n_R2[i][1][3]]) +60*Rho_M[n_R2[i][0][2]]);
			}

			//eulerf_2D(K1,Rho_dt, n_R2, p, D, t, dt, dx, g); //Euler update, stored in K1

			 stringstream m5;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m5 << "A7RK4 Rho at t" << t <<endl; cout << m5.str();
			//Now that K1, K2, K3, K4 updated for all sites, do the stochastic update.

			poisson_distribution<int> poisson; gamma_distribution <double> gamma;

			//Impproved Dornic  Integration of Linear & Stochastic Term
			//CREDITS: PAULA VILLA MARTIN.
			

			for(int i=0;i<Rho_dt.size();i++)
			{
	        Rho_dt[i]+= (dt/6.0)*( K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);
			// Deterministic RK4 complete.
			//Rho_dt[i] = K1[i];
			Rhog[i] = Rho_dt[i]; //Stores Rho_dt RK4 update values (which serves as initial conditions for non-linear terms).
			double mu = lambda*lambda_exp*Rho_dt[i];
			cout << Rho_dt[i] << " ";

	        poisson = poisson_distribution<int>(lambda*lambda_exp*Rho_dt[i]);
			//cout << "Poisson Pill:\t" << poisson(rng) <<endl;
	        gamma = gamma_distribution<double>(poisson(rng), 1.0);
					//Defining shapes of poisson and gamma distributions.
	        Rho_dt[i]= gamma(rng)/lambda;

	    	} cout<<endl;
			cout << "MOULIN ROUGE Final Rho(t) at t:" << t <<endl;
			for(int i=0;i<Rho_dt.size();i++)
			{ cout<< Rho_dt[i] << " "; } cout<<endl;

			

			stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m4 << "SOCKET!" <<endl; cout << m4.str();
			
			t+=dt; //Update timestep

			//Updating Pijm at end of sweep.

			//cout << t << endl;
			if(t > t_stop[s] -dt/2.0 && t <= t_stop[s] + dt/2.0)
			{
				stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
      			m4 << "LOCKET!:\t" << t << " " << p  <<endl; cout << m4.str();
				for ( int k = 0; k< Rho_dt.size(); k++)
			  {	// Recall Rho = | p (Order Parameter) | Replicate (r) | Time(t)  | a*L + b | Cell Density[a][b] |
					Rho.push_back({p, static_cast<double>(j), t, static_cast<double>(k), Rho_dt[k]});
			  }
				s = (s+1)%t_stop.size(); // To prevent out of bounds exception.
			}

			KRi1.push_back({static_cast<double>(j), t, K1[i1], K2[i1], K3[i1], K4[i1], 
			static_cast<double>(Rho_dt[i1]-Rhog[i1]), Rho_dt[i1], diff[0][i1], diff[1][i1]});

			KRi2.push_back({static_cast<double>(j), t, K1[i2], K2[i2], K3[i2], K4[i2], 
			static_cast<double>(Rho_dt[i2]-Rhog[i2]), Rho_dt[i2], diff[0][i2], diff[1][i2]});
			
			KRi3.push_back({static_cast<double>(j), t, K1[i3], K2[i3], K3[i3], K4[i3], 
			static_cast<double>(Rho_dt[i3]-Rhog[i3]), Rho_dt[i3], diff[0][i3], diff[1][i3]});

		} // End of t loops



		vector<double>().swap(Rho_M); vector<double>().swap(Rho_dt);
		// Remove dynamic vectors from memory and free up allocated space.

		

	}

	stringstream L, dT ,p1, dX, rini, Dm, Sm;

  	L << g; dT << dt; p1 << setprecision(4) << p; dX << setprecision(4) << dx;
  	rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma;
  	// setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_dp1; ofstream output_dp2; ofstream output_dp3;
  	// Creating a file instance called output to store output data as CSV.
	output_dp1.open("New/cDebug_Modulei1_G_" + L.str() + "_dt_" + dT.str() + "_dx_"+ dX.str() +
	"_p_"+ p1.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");
	output_dp2.open("New/cDebug_Modulei2_G_" + L.str() + "_dt_" + dT.str() + "_dx_"+ dX.str() +
	"_p_"+ p1.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");
	output_dp3.open("New/cDebug_Modulei3_G_" + L.str() + "_dt_" + dT.str() + "_dx_"+ dX.str() +
	"_p_"+ p1.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	// Output = | Replicate (r) | Time(t)  | K1 | ... | K4 | Stoc Rho_dt | Stoc - Dem Rho_dt |
	output_dp1 << " # Tr No , t ,  K1 , K2 , K3 , K4 , Stoc - Dem Rho_ij(t), Stoc Rho_ij(t), Diff_K1, Diff_K2, \n";
	output_dp2 << " # Tr No , t ,  K1 , K2 , K3 , K4 , Stoc - Dem Rho_ij(t), Stoc Rho_ij(t), Diff_K1, Diff_K2, \n";
	output_dp3 << " # Tr No , t ,  K1 , K2 , K3 , K4 , Stoc - Dem Rho_ij(t), Stoc Rho_ij(t), Diff_K1, Diff_K2, \n";
	cout << "The vector elements are: "<< endl;
  	cout << " p , # Tr No , t ,  i*L + j ,  Rho_ij(t) \n";

	for(int i=0; i< KRi1.size(); i++)
	{
		output_dp1 << static_cast<int> (KRi1[i][0])  << "," << setprecision(6) << KRi1[i][1] << "," << setprecision(9)
		<< KRi1[i][2] << "," << setprecision(9) << KRi1[i][3] << "," << setprecision(9) << KRi1[i][4] << "," 
		<< setprecision(9) << KRi1[i][5] << ","<< setprecision(16) << KRi1[i][6] << "," << setprecision(16) << KRi1[i][7] 
		<< "," << setprecision(16) << KRi1[i][8] << ","<< setprecision(16) << KRi1[i][9] << endl;

		output_dp2 << static_cast<int> (KRi2[i][0])  << "," << setprecision(6) << KRi2[i][1] << "," << setprecision(9)
		<< KRi2[i][2] << "," << setprecision(9) << KRi2[i][3] << "," << setprecision(9) << KRi2[i][4] << "," 
		<< setprecision(9) << KRi2[i][5] << ","<< setprecision(16) << KRi2[i][6] << ","<< setprecision(16) << KRi2[i][7] 
		<< ","  << setprecision(16) << KRi2[i][8] << ","<< setprecision(16) << KRi2[i][9] << endl;

		output_dp3 << static_cast<int> (KRi3[i][0])  << "," << setprecision(6) << KRi3[i][1] << "," << setprecision(9)
		<< KRi3[i][2] << "," << setprecision(9) << KRi3[i][3] << "," << setprecision(9) << KRi3[i][4] << "," 
		<< setprecision(9) << KRi3[i][5] << ","<< setprecision(16) << KRi3[i][6] << ","<< setprecision(16) << KRi3[i][7] 
		<< "," << setprecision(16) << KRi3[i][8] << "," << setprecision(16) << KRi3[i][9] <<  endl;

	}
	output_dp1.close(); output_dp2.close(); output_dp3.close();



}

void test_gamma_poisson()
{
	int nrolls=10000;  // number of experiments
  	int nstars=600;    // maximum number of stars to distribute
	for(int j=0; j< 4; j++)
	{
		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j <<endl; cout << m1.str();
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		thread_local std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd);

		poisson_distribution<int> poisson(3.7); gamma_distribution <double> gamma(2.0,2.0);

		int g[24]={};

  		for (int i=0; i<nrolls; ++i) 
		{
    	double number = gamma(rng);
    	if (number<12) ++g[int(2*number)];
  		}

  		std::cout << "gamma_distribution (2.0,2.0):" << std::endl;

  		for (int i=0; i<24; ++i) 
		{
			double vi = i/2.0;
			if(vi < 9.5)
			{
				std::cout << setprecision(3) << vi << "-" << setprecision(3) << (vi+ 0.5) << ": ";
			}
			else if(vi == 9.5)
			{
				std::cout << setprecision(3)  << vi << "-" << setprecision(2) << (vi+ 0.5) << ": ";
			}
			else
			{ std::cout << setprecision(3)  << vi << "-" << setprecision(3) << (vi+ 0.5) << ": "; }
    		
    		std::cout << std::string(g[i]*nstars/nrolls,'*') << std::endl;
  		}
		

		int p[12]={};
		for (int i=0; i<nrolls; ++i) 
		{
    		int number = poisson(rng);
    		if (number<12) ++p[number];
  		}

  		std::cout << "poisson_distribution (mean=3.75):" << std::endl;
  		for (int i=0; i<12; ++i)
    		std::cout << i << ": " << std::string(p[i]*nstars/nrolls,'*') << std::endl;
		
	}
}


//----------------------------- Stochastic SPDP Dornic Integration Machinery --------------------------------------------//



void percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, double Rh0[],
	 double t_max, double a, double b, double D, double sigma, double dt, double dx, int r,  int g)
{
	//int nR2[g*g][2][5] ={0}; //n_R2 stores neighbours of all sites in a ball of radius 2.
	int nR2[g*g][3][3] ={0}; //nR2 stores neighbours of all sites in a ball of NORM 1 (SQUARE) radius 1.
	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_S8(nR2, g); //Assigning neigbours.

	//int tot_iter = static_cast<int>(t_max/dt); // Stores total number of iterations (i.e. time steps logged) per replicate.
	int tot_iter = static_cast<int>(1.0/dt) + 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	double rho_rep_avg_var[tot_iter][5] ={0.0}; 
	//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
	double diff_coefficient = D/(dx*dx); //= D*(1.0/(6.0*dx*dx)); //Saves us some trouble later.
	double beta = a - 4*D/(dx*dx); //a - (10.0/3.0)*D/(dx*dx); //double beta = a - 4*D/(dx*dx);
	double lambda_exp = exp( (beta)*dt); double lambda = 2*(beta)/(sigma*sigma*(lambda_exp -1.0));

	double surv_rep[t_meas.size() -1] ={0.0}; //Measures number of surviving replicates at each time point.

	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	int ind=0;
	for(double t =0; t <= 1; t+=dt)
	{ rho_rep_avg_var[ind][0] = t; ind+=1; }
	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[ind+i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep
	
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
		<< "\t with Total Iterations:\t " << tot_iter  <<endl; cout << m1.str();
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		//thread_local std::default_random_engine rng(rd); //engine
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number
		//Defining the Dornic variables
		
		//double K1[g*g] ={0.0}; double K2[g*g] ={0.0}; double K3[g*g] ={0.0}; double K4[g*g] ={0.0};	
		//Initialise RK4 terms for future use.
		vector <double> Rho_dt(g*g, 0.0); //Updated through the turn, in turn used to update DRho
		vector <double> DRho(g*g, 0.0); //Dummy variable, Kept constant, only updated at end of turn.
		double Rho_M[tot_iter][2] ={0.0}; //Stores <rho>x values per time step for a single replicate.

		//cout<< "Initially: "<<endl;
	  	for(int i=0; i< g*g; i++)
	  	{
				Rho_dt[i] = Rh0[i]; //Assigning initial conditions
				DRho[i] = Rh0[i];
	  	} //cout<<endl;

		poisson_distribution<int> poisson; gamma_distribution <double> gamma;

		//int t_max = static_cast<int>(t_stop[t_stop.size() -1]); //Storing last element of the same here.
		double t=0; int index = 0;  //Initialise t
		int s=0; int po=1; int so =1; 
		//int co=1; int eo=1; 

		//stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        //m6 << "INITIAL CONDITIONS ASSIGNED FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() <<endl; cout << m6.str();
		while (t < t_max + dt)
		{
			double rhox_avg, rhox_num;
			//Start by updating the Rho vector.

			if(t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{
				rhox_avg = mean_of_vector(Rho_dt, g*g); //Finds spatial average of densities at given t.
				if(rhox_avg > 0)
					surv_rep[index] += 1;
			}


			if(t <= 1.0 || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{
				rhox_avg = mean_of_vector(Rho_dt, g*g); //Finds spatial average of densities at given t.
				rhox_num = occupied_sites_of_vector(Rho_dt, g*g); //Finds number of occupied at given t.
				Rho_M[s][0] = rhox_num; Rho_M[s][1] = rhox_avg;  s+=1;//Rho_M.push_back({t, rhox_avg});
				if(rhox_num == 0.0)
					break; //Finish negotiations.
				if( t > 1.0)
					index+=1;
			}
			

			// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC.
			for(int i=0;i<Rho_dt.size();i++)
			{
			//int check= twilight_zoneS8(DRho, nR2, i);
			//Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.
			
			double alpha_i = diff_coefficient*(1*(DRho[nR2[i][0][1]] + DRho[nR2[i][2][1]] +
			DRho[nR2[i][1][0]] + DRho[nR2[i][1][2]]));
			/**double alpha_i = - diff_coefficient*(1*(Rho_dt[nR2[i][0][0]] + Rho_dt[nR2[i][0][4]] 
				+ Rho_dt[nR2[i][1][0]] + Rho_dt[nR2[i][1][4]]) -16*(Rho_dt[nR2[i][0][1]] + Rho_dt[nR2[i][0][3]] 
				+ Rho_dt[nR2[i][1][1]] + Rho_dt[nR2[i][1][3]])); **/
			/**double alpha_i = - diff_coefficient*(1*(DRho[nR2[i][0][0]] + DRho[nR2[i][0][4]] 
				+ DRho[nR2[i][1][0]] + DRho[nR2[i][1][4]]) -16*(DRho[nR2[i][0][1]] + DRho[nR2[i][0][3]] 
				+ DRho[nR2[i][1][1]] + DRho[nR2[i][1][3]]));*/
			/**double alpha_i = diff_coefficient*(1*(DRho[nR2[i][0][0]] + DRho[nR2[i][0][2]] 
				+ DRho[nR2[i][2][0]] + DRho[nR2[i][2][2]]) +4*(DRho[nR2[i][0][1]] + DRho[nR2[i][2][1]] 
				+ DRho[nR2[i][1][0]] + DRho[nR2[i][1][2]])); **/

			if(alpha_i == 0 && DRho[nR2[i][1][1]] == 0) 
			{
				continue; 
			} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.
			
			double mu = -1.0 + 2.0*alpha_i/(sigma*sigma);
	        poisson = poisson_distribution<int>(lambda*lambda_exp*Rho_dt[i]);
			
	        gamma = gamma_distribution<double>(mu + 1.0 + poisson(rng), 1.0);
					
	        Rho_dt[i]= gamma(rng)/lambda;

			}

			//Euler Integration of remaining terms
			for(int i=0;i<Rho_dt.size();i++)
			{
			Rho_dt[i] = (Rho_dt[i])/( 1 + b*Rho_dt[i]*dt); 
			DRho[i] =Rho_dt[i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
			}

			t+=dt; //Update timestep

		} // End of t loops

		if(j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "INITIAL UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << tot_iter << endl; cout << m7.str();
			for(int i =0; i< tot_iter; i++)
			{	rho_rep_avg_var[i][1] = Rho_M[i][1]; rho_rep_avg_var[i][2] = 0.0; rho_rep_avg_var[i][4] = Rho_M[i][0];
				if(Rho_M[i][1] > 0) //Ensuring only surviving runs are considered
				{ rho_rep_avg_var[i][3] = 1;  } //One surviving run as of now
				else if(Rho_M[i][1] == 0) //Ensuring only surviving runs are considered
				{ rho_rep_avg_var[i][3] = 0;  } //No surviving run as of now

			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			/**stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "LATER UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << Rho_M.size() << endl; cout << m7.str();**/
			// Second or higher replicate, use incremental advances.
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		vector<double>().swap(Rho_dt); vector<double>().swap(DRho); //vector<vector <double>>().swap(Rho_M);
		// Remove dynamic vectors from memory and free up allocated space.

		if(j == 2)
		{
			stringstream L, tm ,d3, p1, rini, Dm, Sm;

  			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a;
  			rini << j; Sm << setprecision(3) << sigma;
			// Three replicates are over.
			ofstream output_dp;
  			// Creating a file instance called output to store output data as CSV.
			output_dp.open("../Data/SPDP/PRELIMINARY_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a1_"+ p1.str() +
			"_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

			// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
			output_dp << " a , L, t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";
			for(int i=0; i< tot_iter; i++)
			{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |
				if(i < ind)
				{ output_dp << a << "," << g << "," << rho_rep_avg_var[i][0] << "," << rho_rep_avg_var[i][1] << "," <<
				rho_rep_avg_var[i][2] << "," << rho_rep_avg_var[i][3] << "," << rho_rep_avg_var[i][4] << "," << j <<endl; }
				else
				{
					output_dp << a << "," << g << "," << rho_rep_avg_var[i][0] << "," << rho_rep_avg_var[i][1] << "," <<
					rho_rep_avg_var[i][2] << "," << rho_rep_avg_var[i][3] << "," << rho_rep_avg_var[i][4] << "," << surv_rep[i -ind] <<endl;
				}
			}

			output_dp.close();
		} //End of logging

	} // End of r loop.

	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|    L 		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |
		if(i < ind)
		{ Rho.push_back({a, static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], 
				rho_rep_avg_var[i][2], rho_rep_avg_var[i][3], rho_rep_avg_var[i][4], static_cast<double>(r)}); }
		else
		{
			Rho.push_back({a, static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], 
				rho_rep_avg_var[i][2], rho_rep_avg_var[i][3], rho_rep_avg_var[i][4], surv_rep[i -ind]});
		}
	}
}

void spreading_exp_Dornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, double Rh0[],
	 double t_max, double a, double b, double D, double sigma, double dt, double dx, int r, int g , int c_x, int c_y)
{
	//Used to determine the critical exponents derived from spreading experiments at a =a_c (central point, (c_x, c_y))
	int nR2[g*g][3][3] ={0}; //nR2 stores neighbours of all sites in a ball of NORM 1 (SQUARE) radius 1.
	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_S8(nR2, g); //Assigning neigbours.

	int tot_iter = static_cast<int>(1.0/dt) + 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	
	double diff_coefficient = D/(dx*dx); // D*(1/(6.0*dx*dx)); //Saves us some trouble later.
	double beta = a - 4*D/(dx*dx); //a - (10.0/3.0)*D/(dx*dx); //double beta = a - 4*D/(dx*dx);
	double lambda_exp = exp( (beta)*dt); double lambda = 2*(beta)/(sigma*sigma*(lambda_exp -1.0));

	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events
	double rho_t_N_R2_rep[tot_iter][4] ={0.0}; 
	//Stores time, <N(t)>r, <R2(t)>r, r(t) where r is the number of surviving runs (at t) respectively.
	int ind=0;
	for(double t =0; t <= 1; t+=dt)
	{ rho_t_N_R2_rep[ind][0] = t; ind+=1; }
	for(int i =0; i< t_meas.size(); i++)
	{ rho_t_N_R2_rep[ind+i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep

	int jo=0;
	
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.
		if( j%100 == 5)
		{
			stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m1 << "Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< "\t with Total Iterations:\t " << tot_iter  <<endl; cout << m1.str();
		}
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		poisson_distribution<int> poisson; gamma_distribution <double> gamma;
		//thread_local std::default_random_engine rng(rd); //engine
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number
		//Defining the Dornic variables
		vector <double> Rho_dt(g*g, 0.0); //Updated through the turn, in turn used to update DRho
		vector <double> DRho(g*g, 0.0); //Dummy variable, Kept constant, only updated at end of turn.
		double Rho_M[tot_iter][3] ={0.0}; //Stores N(t), R2(t) and Surviving Rep # (t) per time step for a single replicate.
	  	for(int i=0; i< g*g; i++)
	  	{
				Rho_dt[i] = Rh0[i]; //Assigning initial conditions
				DRho[i] = Rh0[i];
	  	} //cout<<endl;
		double t=0; //Initialise t
		// int po=1; int so =1; 
		int index = 0; int s=0;
		//int co=1; int eo=1;
		while (t < t_max + dt)
		{
			if( theborderwall(Rho_dt, g) == 1)
			{
				if( jo%10 == 0)
				{
					stringstream m2;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m2 << "Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() << "\t at time t:\t " << setprecision(6) << t
				<< "\t IS DOWN FOR THE COUNT!\t "  <<endl; cout << m2.str(); jo+=1;
				}
				jo+=1;
				break; //Terminate rest of time-based simulation if the metastatic spread reaches your extremeties. Cut out the tumour before it gets too large!
			}	
			//Start by updating the Rho vector.
			if(t <= 1.0 || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{	//Credits: https://www.educative.io/answers/how-to-return-multiple-values-from-a-function-in-cpp17
				auto [sq_dist, N] = meansq_spread_of_vector(Rho_dt, g, c_x, c_y); //Finds spatial average of densities at given t.
				if(N == 0)
					break; //No living cells in grid, will remain this way.
				Rho_M[s][0] = N; Rho_M[s][1] = sq_dist/N; Rho_M[s][2] = 1.0;  s+=1;//Rho_M.push_back({t, rhox_avg});
				if( t > 1.0)
					index+=1;
			}		
			// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC.
			for(int i=0;i<Rho_dt.size();i++)
			{
			double alpha_i = diff_coefficient*(1*(DRho[nR2[i][0][1]] + DRho[nR2[i][2][1]] +
			DRho[nR2[i][1][0]] + DRho[nR2[i][1][2]]));
			/** double alpha_i = D/(dx*dx)*(1*(Rho_dt[nR2[i][0][1]] + Rho_dt[nR2[i][0][3]] +
			Rho_dt[nR2[i][1][1]] + Rho_dt[nR2[i][1][3]]));*/
			/**double alpha_i = - diff_coefficient*(1*(Rho_dt[nR2[i][0][0]] + Rho_dt[nR2[i][0][4]] 
				+ Rho_dt[nR2[i][1][0]] + Rho_dt[nR2[i][1][4]]) -16*(Rho_dt[nR2[i][0][1]] + Rho_dt[nR2[i][0][3]] 
				+ Rho_dt[nR2[i][1][1]] + Rho_dt[nR2[i][1][3]])); **/
			/**double alpha_i = - diff_coefficient*(1*(DRho[nR2[i][0][0]] + DRho[nR2[i][0][4]] 
				+ DRho[nR2[i][1][0]] + DRho[nR2[i][1][4]]) -16*(DRho[nR2[i][0][1]] + DRho[nR2[i][0][3]] 
				+ DRho[nR2[i][1][1]] + DRho[nR2[i][1][3]]));*/
			/**double alpha_i = diff_coefficient*(1*(DRho[nR2[i][0][0]] + DRho[nR2[i][0][2]] 
				+ DRho[nR2[i][2][0]] + DRho[nR2[i][2][2]]) +4*(DRho[nR2[i][0][1]] + DRho[nR2[i][2][1]] 
				+ DRho[nR2[i][1][0]] + DRho[nR2[i][1][2]])); **/
			
			if(alpha_i == 0 && DRho[nR2[i][1][1]] == 0) 
			{	continue; } //Rho_dt[i] will remain 0 as neighbours are also at 0.
			double mu = -1 + 2.0*alpha_i/(sigma*sigma);
	        poisson = poisson_distribution<int>(lambda*lambda_exp*Rho_dt[i]);
			
	        gamma = gamma_distribution<double>(mu + 1 + poisson(rng), 1.0);
					
	        Rho_dt[i]= gamma(rng)/lambda;

			}

			//Euler Integration of remaining terms
			for(int i=0;i<Rho_dt.size();i++)
			{
			Rho_dt[i] = (Rho_dt[i])/( 1 + b*Rho_dt[i]*dt); 
			DRho[i] =Rho_dt[i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
			}

			t+=dt; //Update timestep

		} // End of t loops

		if(j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "INITIAL UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << tot_iter << endl; cout << m7.str();
			for(int i =0; i< tot_iter; i++)
			{	rho_t_N_R2_rep[i][1] = Rho_M[i][0]; rho_t_N_R2_rep[i][2] = Rho_M[i][1]; 
				if(Rho_M[i][0] > 0) //Ensuring only surviving runs are considered
				{ rho_t_N_R2_rep[i][3] = 1;  } //One surviving run as of now
				else if(Rho_M[i][0] == 0) //Ensuring only surviving runs are considered
				{ rho_t_N_R2_rep[i][3] = 0;  } //No surviving run as of now

			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{

			/**stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "LATER UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << Rho_M.size() << endl; cout << m7.str();**/
			// Second or higher replicate, use incremental advances.
			spreading_incremental_surv_runs(rho_t_N_R2_rep, Rho_M, tot_iter); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
			continue;
		}

		vector<double>().swap(Rho_dt); //vector<vector <double>>().swap(Rho_M);
		// Remove dynamic vectors from memory and free up allocated space.

	} // End of r loop.

	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|    t 		|     <N(t)>r			|    <MSD(N(t))>r    |    #Surviving Runs    |
		int tr= omp_get_thread_num();
		Rho.push_back({static_cast<double>(tr), a, rho_t_N_R2_rep[i][0], rho_t_N_R2_rep[i][1], rho_t_N_R2_rep[i][2], rho_t_N_R2_rep[i][3]});
	}
	stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m7 << "IN SUMMATION WE HAVE A BREACH COUNT OF:\t" << jo << "\t for Thread Rank:\t " << omp_get_thread_num() 
	<< endl; cout << m7.str();
			// Second or higher replicate, use incremental advances.

}



void critical_exp_delta(double Rho_0[], int div, double t_max, double a_start, double a_end, 
	double b, double D, double sigma, double dt, double dx, int r,  int g)
{

	init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "FUCK" <<endl;
  // The pspace to iterate over.
    vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  for (int i=0; i< t_measure.size(); i++)
  {
	cout << t_measure[i] << " ";
  } cout << endl;

  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = 5 * microToSeconds;     //5 seconds
    
    cout<<"Delay 1 in progress... ("<<delay1/microToSeconds<<"s)"<<endl;
	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 25)
	{
		omp_set_num_threads(25); //Limiting use on Chunk.
	}

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
		**/
		percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_space[i], b, D, sigma, dt, dx, r, g);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting critical exponent data for each grid size in order.
          stringstream message3;
          message3 << "Is this happening?\n";
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, tm ,d3, p1, p2, rini, Dm, Sm;

  L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a_start; p2 << setprecision(4) << a_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_dp;
  // Creating a file instance called output to store output data as CSV.
	output_dp.open("../Data/SPDP/DEBUG_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a1_"+ p1.str() +
	"_a2_"+ p2.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	cout << "Save file name: " <<endl;
	cout << "../Data/SPDP/DEBUG_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a1_"+ p1.str() +
	"_a2_"+ p2.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv" << endl;

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
	output_dp << " a , L,  t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";
	cout << "The vector elements are: "<< endl;
  	cout << "  a , L, t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_dp << setprecision(8) << vec[i][0] << "," << vec[i][1] << "," << vec[i][2] << "," << setprecision(16)
		<< vec[i][3] << "," << setprecision(16) << vec[i][4] << "," << vec[i][5] << "," <<
		setprecision(16) << vec[i][6] << "," << vec[i][7] << endl;
		if( i%(1000) ==1)
    {
			cout << setprecision(5) << vec[i][0] << "," << vec[i][1] << "," << vec[i][2] << "," << setprecision(9)
			<< vec[i][3] << "," <<  setprecision(9) << vec[i][4] << "," << vec[i][5] << "," << vec[i][6] << "," << vec[i][7] << endl;
    }
	}
	output_dp.close();
}

void critical_exp_finite_scaling_collapse(int div, int g_st, int g_end, double t_max, double a_c, 
	double b, double D, double sigma, double dt, double dx, int r)
{
	
	vector<double> gv_space = linspace(g_st, g_end, div);
	vector<int> g_space(gv_space.size()); //Convert to INT.
	for (int i=0; i < g_space.size(); i++)
		g_space[i] = static_cast<int>(gv_space[i]);
	vector<double>().swap(gv_space);
  // The grid sizes to iterate over.
    vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  for (int i=0; i< t_measure.size(); i++)
  {
	cout << t_measure[i] << " ";
  } cout << endl;

  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = 5 * microToSeconds;     //5 seconds
    
    cout <<"Delay 1 in progress... ("<<delay1/microToSeconds<<"s)"<<endl;
	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 25)
	{
		omp_set_num_threads(25); //Limiting use on Chunk.
	}

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */
	  
      #pragma omp for nowait
      for (int i=0; i < g_space.size(); i++)
      {
        //type="Gam";

        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on Grid Size:\t" << g_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
		double Rho_0[g_space[i]*g_space[i]];
		init_fullframe(Rho_0, g_space[i]*g_space[i]); //Returns Rho_0 with a full initial frame filled with ones.
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
		**/
		percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_c, b, D, sigma, dt, dx, r, g_space[i]);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      
    	#pragma omp critical
    	vec.insert(vec.end(), vec_private.begin(), vec_private.end());
    	// Inserting critical exponent data for each grid size in order.
    	stringstream message3;
    	message3 << "Is this happening?\n";
    	cout << message3.str();

      
		vector<vector <double>>().swap(vec_private);
		// Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream AC, tm ,d3, g1, g2, rini, Dm, Sm;

  tm << t_max; d3 << setprecision(3) << dt; g1 << setprecision(4) << g_st; g2 << setprecision(4) << g_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma; AC << a_c;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_dp;
  // Creating a file instance called output to store output data as CSV.
	output_dp.open("../Data/SPDP/DEBUG_P_c_Finite_Scaling_Collapse_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + AC.str() + 
	 "_L1_"+ g1.str() + "_L2_"+ g2.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
	output_dp << " a , L, t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";
	cout << "The vector elements are: "<< endl;
  	cout << "  a , L, t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_dp << setprecision(8) << vec[i][0] << "," << vec[i][1] << "," << setprecision(16)
		<< vec[i][2] << "," << setprecision(16) << vec[i][3] << "," << vec[i][4] << "," <<
		setprecision(16) << vec[i][5] << "," << vec[i][6] << endl;
		if( i%(1000) ==1)
    {
		cout << setprecision(5) << vec[i][0] << "," << vec[i][1] << "," << setprecision(9)
		<< vec[i][2] << "," <<  setprecision(9) << vec[i][3] << "," << vec[i][4] << "," << vec[i][5] << "," << vec[i][6] << endl;
    }
	}
	output_dp.close();
}

void critical_exp_finite_scaling_stationary(int div, int g_st, int g_end, double t_lim, double t_max, double a_c, 
	double b, double D, double sigma, double dt, double dx, int r)
{
	
	vector<double> gv_space = lnspace(g_st, g_end, div);
	vector<int> g_space(gv_space.size()); //Convert to INT.
	for (int i=0; i < g_space.size(); i++)
		g_space[i] = static_cast<int>(gv_space[i]);
	vector<double>().swap(gv_space);
  // The grid sizes to iterate over.
    //vector <double> t_measure = linspace(t_lim, t_max, 30);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  /**cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  for (int i=0; i< t_measure.size(); i++)
  {
	cout << t_measure[i] << " ";
  } cout << endl; **/

  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = 5 * microToSeconds;     //5 seconds
    
    cout <<"Delay 1 in progress... ("<<delay1/microToSeconds<<"s)"<<endl;
	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 25)
	{
		omp_set_num_threads(25); //Limiting use on Chunk.
	}

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

	  double t_start; double t_end; vector <double> t_measure;
	  
      #pragma omp for nowait
      for (int i=0; i < g_space.size(); i++)
      {
        //type="Gam";

        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on Grid Size:\t" << g_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
		double Rho_0[g_space[i]*g_space[i]];
		init_fullframe(Rho_0, g_space[i]*g_space[i]); //Returns Rho_0 with a full initial frame filled with ones.
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
		**/

		if(g_space[i] > 30 && g_space[i] <= 50)
		{	t_start = 300; t_end = 900;}
		else if(g_space[i] > 50 && g_space[i] <= 100)
		{	t_start = 1000; t_end = 2500;}
		else if(g_space[i] > 100 && g_space[i] <= 150)
		{	t_start = 1250; t_end = 8000;}
		else if(g_space[i] > 150 && g_space[i] <= 200)
		{	t_start = 1400; t_end = 10000;}
		else if(g_space[i] > 200)
		{	t_start = 2000; t_end = 15000;}
		t_measure = linspace(t_start, t_end, 30);
		percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_end+10, a_c, b, D, sigma, dt, dx, r, g_space[i]);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      
    	#pragma omp critical
    	vec.insert(vec.end(), vec_private.begin(), vec_private.end());
    	// Inserting critical exponent data for each grid size in order.
    	stringstream message3;
    	message3 << "Is this happening?\n";
    	cout << message3.str();

      
		vector<vector <double>>().swap(vec_private);
		// Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream AC, tm ,d3, g1, g2, rini, Dm, Sm;

  tm << t_max; d3 << setprecision(3) << dt; g1 << setprecision(4) << g_st; g2 << setprecision(4) << g_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma; AC << a_c;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_dp;
  // Creating a file instance called output to store output data as CSV.
	output_dp.open("../Data/SPDP/DEBUG_P_c_Finite_Scaling_Stationary_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + AC.str() + 
	 "_L1_"+ g1.str() + "_L2_"+ g2.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
	output_dp << " a , L, t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";
	cout << "The vector elements are: "<< endl;
  	cout << "  a , L, t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_dp << setprecision(8) << vec[i][0] << "," << vec[i][1] << "," << setprecision(16)
		<< vec[i][2] << "," << setprecision(16) << vec[i][3] << "," << vec[i][4] << "," <<
		setprecision(16) << vec[i][5] << "," << vec[i][6] << endl;
		if( i%(1000) ==1)
    {
		cout << setprecision(5) << vec[i][0] << "," << vec[i][1] << "," << setprecision(9)
		<< vec[i][2] << "," <<  setprecision(9) << vec[i][3] << "," << vec[i][4] << "," << vec[i][5] << "," << vec[i][6] << endl;
    }
	}
	output_dp.close();
}


void critical_exp_theta_delta(double Rho_0[], int div, double t_max, double a_c, double b, 
		double D, double sigma, double dt, double dx, int r,  int g)
{

	init_solitarytear(Rho_0, g); //Returns Rho_0 with a single central seed admist a sea of drowned despair.
	int c_x; int c_y; c_x = c_y = static_cast<int>(g/2);

    vector <double> t_measure = logarithm10_time_bins(t_max, dt); 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  for (int i=0; i< t_measure.size(); i++)
  {
	cout << t_measure[i] << " ";
  } cout << endl;

	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 22)
	{
		omp_set_num_threads(22); //Limiting use on Chunk.
	}

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */
	  
      #pragma omp for nowait
      for (int i=0; i < div; i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Supra-Replicate:\t" << i*r << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
		**/
		spreading_exp_Dornic_2D(CExpRho_a, t_measure, Rho_0, t_max, a_c, b, D, sigma, dt, dx, r, g , c_x, c_y);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      #pragma omp critical
      
      vec.insert(vec.end(), vec_private.begin(), vec_private.end());
      // Inserting critical exponent data for each grid size in order.
      stringstream message3;
      message3 << "Is this happening?\n";
      cout << message3.str();

      
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, tm ,d3, p1, p2, rini, Dm, Sm;

  L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(5) << a_c; //p2 << setprecision(4) << a_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_dp;
  // Creating a file instance called output to store output data as CSV.
	output_dp.open("../Data/SPDP/SPREADING_P_c_ThDel_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a1_"+ p1.str() +
	"_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	// Output =   | 	a		|    t 		|     <N(t)>r			|    <MSD(N(t))>r    |    #Surviving Runs    |
	output_dp << " #, a , t ,  <N(t)>_r, <MSD(N(t))>_r, <P(t)>_r \n";
	cout << "The vector elements are: "<< endl;
  	cout << " #, a , t ,  <N(t)>_r, <MSD(N(t))>_r, <P(t)>_r \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_dp << vec[i][0] << "," << vec[i][1] << "," << setprecision(16)
		<< vec[i][2] << "," << setprecision(16) << vec[i][3] << "," << vec[i][4] << "," << vec[i][5] << endl;
		if( i%(1000) ==1)
    {
			cout << vec[i][0] << "," << vec[i][1] << "," << setprecision(9)
			<< vec[i][2] << "," << setprecision(9) << vec[i][3] << "," << vec[i][4] << "," << vec[i][5] << endl;
    }
	}
	output_dp.close();
}

void f_exp2D(vector<double> &f, vector<double> &Rho_t, double a, double b, double c, double D, double t, double dt, double dx, int g)
{
	for(int i=0; i < g*g; i++)
	{
		f[i] = -b*Rho_t[i]*Rho_t[i] - c*Rho_t[i]*Rho_t[i]*Rho_t[i];
		// RK4 update of Non-linear parts
	}
}

void RK4_Integrate(vector<double> &Rho_t, double a, double b, double c, double D, double t, double dt, double dx, int g)
{

	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	double dt6 = dt/6.0; double dt2 = dt/2.0; 
	for(int i=0; i < g*g; i++)
	{
		double K1, K2, K3, K4, Rho_M;
		K1 = -b*Rho_t[i]*Rho_t[i] - c*Rho_t[i]*Rho_t[i]*Rho_t[i];
		Rho_M = Rho_t[i] + (dt2)*K1;
		K2 = -b*Rho_M*Rho_M - c*Rho_M*Rho_M*Rho_M; Rho_M = Rho_t[i] + (dt2)*K2;
		K3 = -b*Rho_M*Rho_M - c*Rho_M*Rho_M*Rho_M; Rho_M = Rho_t[i] + (dt)*K3;
		K4 = -b*Rho_M*Rho_M - c*Rho_M*Rho_M*Rho_M;

		Rho_t[i]+= (dt6)*( K1 + 2.0*K2 + 2.0*K3 + K4);

		if(K1 < 0 || K2 < 0 || K3 < 0 || K4 < 0 || isfinite(Rho_t[i]) == false)
		{
			if( Rho_t[i] < 0 || isfinite(Rho_t[i]) == false || Rho_t[i] > 5)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RK4 WAS KO'ED WITH:\t" << Rho_t[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with K1[i]:   " << K1 << "\t, K2[i]:\t" << K2 
				<< "\t, K2[i]:\t" << K3 << "\t AND K4[i]:\t" << K4 << endl; //cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			
		}
	}
}



void expanded_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, double Rh0[],
	 double t_max, double a, double b, double c, double D, double sigma, double dt, double dx, int r,  int g)
{
	//Integrating first order PDE: arho - brho^2 - crho^3 +Diff + Stocj
	//int nR2[g*g][2][5] ={0}; //n_R2 stores neighbours of all sites in a ball of radius 2.
	int nR2[g*g][3][3] ={0}; //nR2 stores neighbours of all sites in a ball of NORM 1 (SQUARE) radius 1.
	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_S8(nR2, g); //Assigning neigbours.

	//int tot_iter = static_cast<int>(t_max/dt); // Stores total number of iterations (i.e. time steps logged) per replicate.
	int tot_iter = static_cast<int>(1.0/dt) + 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	double rho_rep_avg_var[tot_iter][5] ={0.0}; 
	//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
	double diff_coefficient = D/(dx*dx); //= D*(1.0/(6.0*dx*dx)); //Saves us some trouble later.
	double beta = a - 4*D/(dx*dx); //a - (10.0/3.0)*D/(dx*dx); //double beta = a - 4*D/(dx*dx);
	double lambda_exp = exp( (beta)*dt); double lambda = 2*(beta)/(sigma*sigma*(lambda_exp -1.0));
	double bdt = b*dt; double cdt = c*dt;

	double surv_rep[t_meas.size() -1] ={0.0}; //Measures number of surviving replicates at each time point.

	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	int ind=0;
	for(double t =0; t <= 1; t+=dt)
	{ rho_rep_avg_var[ind][0] = t; ind+=1; }
	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[ind+i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep

	double alpha_prime = diff_coefficient*4; //Assume each Rho in neighbourhood for alpha estimation averages 1.
	double poiss_ru = lambda*lambda_exp*2.5; //Mean
	double mu_nought = 2.0*alpha_prime/(sigma*sigma) + poiss_ru;
	poiss_ru += 4*sqrt(poiss_ru); //Mean of Poisson Sampler + 4 SD.
	mu_nought += 6*sqrt(mu_nought); //Mean of  Mean of Gamma Sampler + 6 SD.


	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t Lambda:\t" << lambda << "\t Beta:\t " << beta << "\t dt,dx:\t " << dt << "," << dx << "\t Poisson Cut-Off:\t " << poiss_ru 
	<< "\t and  Gamma Cut-Off:\t " << mu_nought  <<endl; //cout << m1.str();
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	
	
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
		<< "\t with Total Iterations:\t " << tot_iter  <<endl; cout << m1.str();
		errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		//thread_local std::default_random_engine rng(rd); //engine
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number
		//Defining the Dornic variables
		
		//double K1[g*g] ={0.0}; double K2[g*g] ={0.0}; double K3[g*g] ={0.0}; double K4[g*g] ={0.0};	
		//Initialise RK4 terms for future use.
		vector <double> Rho_dt(g*g, 0.0); //Updated through the turn, in turn used to update DRho
		vector <double> DRho(g*g, 0.0); //Dummy variable, Kept constant, only updated at end of turn.
		vector <double> DummyRho(g*g, 0.0); //Dummy variable, Kept constant, only updated at end of turn.
		
		double Rho_M[tot_iter][2] ={0.0}; //Stores <rho>x values per time step for a single replicate.

		//double K1[g*g] ={0.0}; double K2[g*g] ={0.0}; double K3[g*g] ={0.0}; double K4[g*g] ={0.0}; //double Rho_Mid[g*g] ={0.0};

		//cout<< "Initially: "<<endl;
	  	for(int i=0; i< g*g; i++)
	  	{
				Rho_dt[i] = Rh0[i]; //Assigning initial conditions
				DRho[i] = Rh0[i];
	  	} //cout<<endl;

		poisson_distribution<int> poisson; gamma_distribution <double> gamma;

		//int t_max = static_cast<int>(t_stop[t_stop.size() -1]); //Storing last element of the same here.
		double t=0; int index = 0;  //Initialise t
		int s=0; int po=1; int so =1; int lo=1;
		//int co=1; int eo=1; 

		//stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        //m6 << "INITIAL CONDITIONS ASSIGNED FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() <<endl; cout << m6.str();
		while (t < t_max + dt)
		{

			double rhox_avg, rhox_num;
			//Start by updating the Rho vector.

			if(t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{
				rhox_avg = mean_of_vector(Rho_dt, g*g); //Finds spatial average of densities at given t.
				if(rhox_avg > 0)
					surv_rep[index] += 1;
			}


			if(t <= 1.0 || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{
				rhox_avg = mean_of_vector(Rho_dt, g*g); //Finds spatial average of densities at given t.
				rhox_num = occupied_sites_of_vector(Rho_dt, g*g); //Finds number of occupied at given t.
				Rho_M[s][0] = rhox_num; Rho_M[s][1] = rhox_avg;  s+=1;//Rho_M.push_back({t, rhox_avg});
				if(rhox_num == 0.0)
					break; //Finish negotiations.
				if( t > 1.0)
					index+=1;
			}
			

			// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC.
			for(int i=0;i<Rho_dt.size();i++)
			{
			double alpha_i = diff_coefficient*(1*(DRho[nR2[i][0][1]] + DRho[nR2[i][2][1]] +
			DRho[nR2[i][1][0]] + DRho[nR2[i][1][2]]));

			
			

			if(alpha_i == 0 && DRho[nR2[i][1][1]] == 0) 
			{
				continue; 
			} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

			if(alpha_i < 0 && po ==1)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "We are NOT OKAY:\t" << Rho_dt[i] << "\t at index:  " << i << " at time:\t:" << t <<
				" with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i << endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				//po=0;
			}
			double maxtek = max({DRho[nR2[i][0][1]],DRho[nR2[i][2][1]],DRho[nR2[i][1][0]], DRho[nR2[i][1][2]]}, maxis);
			if(maxtek - DRho[nR2[i][1][1]] > 5 || maxtek - DRho[nR2[i][1][1]] < -5)
			{
				//Checking if Laplacian is causing our troubles.
				
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "ALPHA GO-GO-TRON Rho[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " at time:\t:" << t <<
				" with Max(Neighbour(Rho)):   " << maxtek << "\t and alpha_i:\t" << alpha_i << endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				//po=0;
			}
			
			double mu = -1.0 + 2.0*alpha_i/(sigma*sigma);
			double ziggy = lambda*lambda_exp*Rho_dt[i]; double gru;
			if(ziggy == 0.0)
			{	//stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		//m6 << "SICK!!" <<endl; cout << m6.str(); 
				gru = mu + 1.0;  }
			else if(ziggy < 0.0 || isnan(ziggy) == true || isinf(ziggy) == true)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "POISSON IS MESSED UP:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i << endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			else
			{
				poisson = poisson_distribution<int>(lambda*lambda_exp*Rho_dt[i]);
				gru = mu + 1.0 + poisson(rng);
			}

			gamma = gamma_distribution<double>(gru, 1.0);
					
	        Rho_dt[i]= gamma(rng)/lambda;

			if(gru < 0 || isnan(gru) == true || isinf(gru) == true)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "GAMMA MESSED UP BIG TIME Rho*[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i <<  
				"\t and GRU:\t" << gru <<endl; //cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}

			if(gru > mu_nought || Rho_dt[i] > 5 || alpha_i > 20 || ziggy > poiss_ru)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "BLOWING UP, Rho*[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i <<  
				"\t and GRU:\t" << gru << "\t and Poisson RU:\t" << ziggy << endl; //cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
	        
	        

			if(isnan(Rho_dt[i]) == true || isinf(Rho_dt[i]) == true)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "YOU HIT ROCK BOTTOM WITH: Rho*[t,i]\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i << 
				"\t and GAMMA ENTRY:\t" << gru << endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			DummyRho[i] =Rho_dt[i];
			}

			//RK4 Integration of remaining terms

			RK4_Integrate(Rho_dt, a, b, c, D, t, dt, dx, g);

			for(int i=0;i<Rho_dt.size();i++)
			{
			//Euler Integration of remaining terms
			//Rho_dt[i] = (1 - (bdt + cdt*Rho_dt[i])*Rho_dt[i])*Rho_dt[i];
			//Rho_dt[i] = Rho_dt[i]/(1 + (dt*Rho_dt[i]*(c*Rho_dt[i] + b)));

			if(isnan(Rho_dt[i]) == true )
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RK4 BLEW IT Rho[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and Analytic Integration Term:\t"
				<< DummyRho[i]<< endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				lo = -1;
			}

			if( Rho_dt[i] < 0 || isinf(Rho_dt[i]) == true and so==1)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RHO FALLS BELOW O:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and Analytic Integration Term:\t" 
				<< DummyRho[i]<<  endl; cout << m6.str(); //so=0;
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			DRho[i] =Rho_dt[i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
			}
			/**
			if( int(t)%200== 0)
			{
				stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
				int gs = int(g/2); m4 << "RHO at 200 in Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() << " at time:\t:" << t << "\n";
				for(int i=int(g/4); i<int(3*g/4); i++)
				{
					m4 << Rho_dt[gs*i] << " ";
				}
        		m4 << endl; cout << m4.str();
			}
			**/

			t+=dt; //Update timestep

			if(lo == -1)
			{
				exit(3);
			}

		} // End of t loops

		if(int(omp_get_thread_num()) == 5)
		{
			size_t currentSize = getCurrentRSS( ); //Check Memory usage.
			size_t peakSize    = getPeakRSS( );
			stringstream m7;
			m7 << "In replicate " << j << ", current size in MB: " << currentSize/(1024.0*1024.0) 
			<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m7.str();
			errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();
		}

		if(j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "INITIAL UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << tot_iter << endl; cout << m7.str();
			for(int i =0; i< tot_iter; i++)
			{	rho_rep_avg_var[i][1] = Rho_M[i][1]; rho_rep_avg_var[i][2] = 0.0; rho_rep_avg_var[i][4] = Rho_M[i][0];
				if(Rho_M[i][1] > 0) //Ensuring only surviving runs are considered
				{ rho_rep_avg_var[i][3] = 1;  } //One surviving run as of now
				else if(Rho_M[i][1] == 0) //Ensuring only surviving runs are considered
				{ rho_rep_avg_var[i][3] = 0;  } //No surviving run as of now

			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			/**stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "LATER UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << Rho_M.size() << endl; cout << m7.str();**/
			// Second or higher replicate, use incremental advances.
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		vector<double>().swap(Rho_dt); vector<double>().swap(DRho); //vector<vector <double>>().swap(Rho_M);
		// Remove dynamic vectors from memory and free up allocated space.

		if(j == 2)
		{
			stringstream L, tm ,d3, p1, rini, Dm, Sm;

  			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a;
  			rini << j; Sm << setprecision(3) << sigma;
			// Three replicates are over.
			ofstream output_dp;
  			// Creating a file instance called output to store output data as CSV.
			output_dp.open("../Data/SPDP/PRELIMINARY_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a1_"+ p1.str() +
			"_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

			// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
			output_dp << " a , c, L, t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";
			for(int i=0; i< tot_iter; i++)
			{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |
				if(i < ind)
				{ output_dp << a << ","<< c << "," << g << "," << rho_rep_avg_var[i][0] << "," << rho_rep_avg_var[i][1] << "," <<
				rho_rep_avg_var[i][2] << "," << rho_rep_avg_var[i][3] << "," << rho_rep_avg_var[i][4] << "," << j <<endl; }
				else
				{
					output_dp << a << "," << c << "," << g << "," << rho_rep_avg_var[i][0] << "," << rho_rep_avg_var[i][1] << "," <<
					rho_rep_avg_var[i][2] << "," << rho_rep_avg_var[i][3] << "," << rho_rep_avg_var[i][4] << "," << surv_rep[i -ind] <<endl;
				}
			}

			output_dp.close();
		} //End of logging

	} // End of r loop.

	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|    L 		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |
		if(i < ind)
		{ Rho.push_back({a, c, static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], 
				rho_rep_avg_var[i][2], rho_rep_avg_var[i][3], rho_rep_avg_var[i][4], static_cast<double>(r)}); }
		else
		{
			Rho.push_back({a,c, static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], 
				rho_rep_avg_var[i][2], rho_rep_avg_var[i][3], rho_rep_avg_var[i][4], surv_rep[i -ind]});
		}
	}

	stringstream m10;
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m10.str();
	errout.open(thr, std::ios_base::app); errout << m10.str(); errout.close();
}


void terminal_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, double Rh0[],
	 double t_max, double a, double b, double c, double D, double sigma, double dt, double dx, int r,  int g)
{
	//Integrating first order PDE: arho - brho^2 - crho^3 +Diff + Stocj
	//int nR2[g*g][2][5] ={0}; //n_R2 stores neighbours of all sites in a ball of radius 2.
	int nR2[g*g][3][3] ={0}; //nR2 stores neighbours of all sites in a ball of NORM 1 (SQUARE) radius 1.
	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_S8(nR2, g); //Assigning neigbours.

	//int tot_iter = static_cast<int>(t_max/dt); // Stores total number of iterations (i.e. time steps logged) per replicate.
	int tot_iter = 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	double rho_rep_avg_var[r][5] ={0.0}; 
	//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
	double diff_coefficient = D/(dx*dx); //= D*(1.0/(6.0*dx*dx)); //Saves us some trouble later.
	double beta = a - 4*D/(dx*dx); //a - (10.0/3.0)*D/(dx*dx); //double beta = a - 4*D/(dx*dx);
	double lambda_exp = exp( (beta)*dt); double lambda = 2*(beta)/(sigma*sigma*(lambda_exp -1.0));
	double bdt = b*dt; double cdt = c*dt;

	//double surv_rep[t_meas.size() -1] ={0.0}; //Measures number of surviving replicates at each time point.

	t_max = t_meas[t_meas.size() -1] + 100; //Only compute until end of sampling events

	for(int i =0; i< r; i++)
	{ rho_rep_avg_var[i][0] = t_meas[t_meas.size() -1];}  //Setting up time index of rho_t_N_R2_rep

	double alpha_prime = diff_coefficient*4; //Assume each Rho in neighbourhood for alpha estimation averages 1.
	double poiss_ru = lambda*lambda_exp*2.5; //Mean
	double mu_nought = 2.0*alpha_prime/(sigma*sigma) + poiss_ru;
	poiss_ru += 4*sqrt(poiss_ru); //Mean of Poisson Sampler + 4 SD.
	mu_nought += 6*sqrt(mu_nought); //Mean of  Mean of Gamma Sampler + 6 SD.


	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t Lambda:\t" << lambda << "\t Beta:\t " << beta << "\t dt,dx:\t " << dt << "," << dx << "\t Poisson Cut-Off:\t " << poiss_ru 
	<< "\t and  Gamma Cut-Off:\t " << mu_nought  <<endl; //cout << m1.str();
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	
	
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
		<< "\t with Total Iterations:\t " << tot_iter  <<endl; cout << m1.str();
		if(j%5 == 0)
			errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		//thread_local std::default_random_engine rng(rd); //engine
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number
		//Defining the Dornic variables
		
		//double K1[g*g] ={0.0}; double K2[g*g] ={0.0}; double K3[g*g] ={0.0}; double K4[g*g] ={0.0};	
		//Initialise RK4 terms for future use.
		vector <double> Rho_dt(g*g, 0.0); //Updated through the turn, in turn used to update DRho
		vector <double> DRho(g*g, 0.0); //Dummy variable, Kept constant, only updated at end of turn.
		vector <double> DummyRho(g*g, 0.0); //Dummy variable, Kept constant, only updated at end of turn.
		
		double Rho_M[2][2] ={0.0}; //Stores <rho>x values at t = 0 and t_meas[0] for a single replicate.

		//double K1[g*g] ={0.0}; double K2[g*g] ={0.0}; double K3[g*g] ={0.0}; double K4[g*g] ={0.0}; //double Rho_Mid[g*g] ={0.0};

		//cout<< "Initially: "<<endl;
	  	for(int i=0; i< g*g; i++)
	  	{
				Rho_dt[i] = Rh0[i]; //Assigning initial conditions
				DRho[i] = Rh0[i];
	  	} //cout<<endl;

		poisson_distribution<int> poisson; gamma_distribution <double> gamma;

		//int t_max = static_cast<int>(t_stop[t_stop.size() -1]); //Storing last element of the same here.
		double t=0; int index = 0;  //Initialise t
		int s=0; int po=1; int so =1; int lo=1;
		//int co=1; int eo=1; 

		//stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        //m6 << "INITIAL CONDITIONS ASSIGNED FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() <<endl; cout << m6.str();
		while (t < t_max + dt)
		{

			double rhox_avg, rhox_num;
			//Start by updating the Rho vector.


			if(t == 0.0 || t >= t_meas[0] -dt/2.0 && t < t_meas[0] +dt/2.0)
			{
				rhox_avg = mean_of_vector(Rho_dt, g*g); //Finds spatial average of densities at given t.
				rhox_num = occupied_sites_of_vector(Rho_dt, g*g); //Finds number of occupied at given t.
				Rho_M[s][0] = rhox_num; Rho_M[s][1] = rhox_avg;  s+=1;//Rho_M.push_back({t, rhox_avg});
				if(rhox_num == 0.0)
					break; //Finish negotiations.
			}

			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC.
			for(int i=0;i<Rho_dt.size();i++)
			{
			double alpha_i = diff_coefficient*(1*(DRho[nR2[i][0][1]] + DRho[nR2[i][2][1]] +
			DRho[nR2[i][1][0]] + DRho[nR2[i][1][2]]));

			
			

			if(alpha_i == 0 && DRho[nR2[i][1][1]] == 0) 
			{
				continue; 
			} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

			if(alpha_i < 0 && po ==1)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "We are NOT OKAY:\t" << Rho_dt[i] << "\t at index:  " << i << " at time:\t:" << t <<
				" with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i << endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				//po=0;
			}
			double maxtek = max({DRho[nR2[i][0][1]],DRho[nR2[i][2][1]],DRho[nR2[i][1][0]], DRho[nR2[i][1][2]]}, maxis);
			if(maxtek - DRho[nR2[i][1][1]] > 5 || maxtek - DRho[nR2[i][1][1]] < -5)
			{
				//Checking if Laplacian is causing our troubles.
				
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "ALPHA GO-GO-TRON Rho[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " at time:\t:" << t <<
				" with Max(Neighbour(Rho)):   " << maxtek << "\t and alpha_i:\t" << alpha_i << endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				//po=0;
			}
			
			double mu = -1.0 + 2.0*alpha_i/(sigma*sigma);
			double ziggy = lambda*lambda_exp*Rho_dt[i]; double gru;
			if(ziggy == 0.0)
			{	//stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		//m6 << "SICK!!" <<endl; cout << m6.str(); 
				gru = mu + 1.0;  }
			else if(ziggy < 0.0 || isnan(ziggy) == true || isinf(ziggy) == true)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "POISSON IS MESSED UP:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i << endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			else
			{
				poisson = poisson_distribution<int>(lambda*lambda_exp*Rho_dt[i]);
				gru = mu + 1.0 + poisson(rng);
			}

			gamma = gamma_distribution<double>(gru, 1.0);
					
	        Rho_dt[i]= gamma(rng)/lambda;

			if(gru < 0 || isnan(gru) == true || isinf(gru) == true)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "GAMMA MESSED UP BIG TIME Rho*[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i <<  
				"\t and GRU:\t" << gru <<endl; //cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}

			if(gru > mu_nought || Rho_dt[i] > 5 || alpha_i > 20 || ziggy > poiss_ru)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "BLOWING UP, Rho*[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i <<  
				"\t and GRU:\t" << gru << "\t and Poisson RU:\t" << ziggy << endl; //cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}

			if(isnan(Rho_dt[i]) == true || isinf(Rho_dt[i]) == true)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "YOU HIT ROCK BOTTOM WITH: Rho*[t,i]\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i << 
				"\t and GAMMA ENTRY:\t" << gru << endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			DummyRho[i] =Rho_dt[i];
			}

			//RK4 Integration of remaining terms

			RK4_Integrate(Rho_dt, a, b, c, D, t, dt, dx, g);

			for(int i=0;i<Rho_dt.size();i++)
			{

			if(isnan(Rho_dt[i]) == true )
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RK4 BLEW IT Rho[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and Analytic Integration Term:\t"
				<< DummyRho[i]<< endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				lo = -1;
			}

			if( Rho_dt[i] < 0 || isinf(Rho_dt[i]) == true and so==1)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RHO FALLS BELOW O:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and Analytic Integration Term:\t" 
				<< DummyRho[i]<<  endl; cout << m6.str(); //so=0;
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			DRho[i] =Rho_dt[i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
			}

			t+=dt; //Update timestep

			if(lo == -1)
			{
				exit(3);
			}

		} // End of t loops

		if(int(omp_get_thread_num()) == 5 && j%5 == 0)
		{
			size_t currentSize = getCurrentRSS( ); //Check Memory usage.
			size_t peakSize    = getPeakRSS( );
			stringstream m7;
			m7 << "In replicate " << j << ", current size in MB: " << currentSize/(1024.0*1024.0) 
			<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m7.str();
			errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();
		}

		rho_rep_avg_var[j][1] = Rho_M[1][1]; rho_rep_avg_var[j][2] = Rho_M[0][1]; rho_rep_avg_var[j][3] = Rho_M[1][0];

		vector<double>().swap(Rho_dt); vector<double>().swap(DRho); //vector<vector <double>>().swap(Rho_M);
		// Remove dynamic vectors from memory and free up allocated space.

	} // End of r loop.

	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	for(int j=0; j < r; j++)
	{	// Recall Rho is: | 	a	 |	c	|    L 		|     r		|    t 		|     <Rho(0)>x			|    <Rho(t)>x	   |       #Active Sites (at t) |
		
		Rho.push_back({a,c, static_cast<double>(g), static_cast<double>(j), rho_rep_avg_var[j][0], rho_rep_avg_var[j][2], 
			rho_rep_avg_var[j][1], rho_rep_avg_var[j][3]});
		
	}

	stringstream m10;
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m10.str();
	errout.open(thr, std::ios_base::app); errout << m10.str(); errout.close();
}



void first_order_critical_exp_delta(double Rho_0[], int div, double t_max, double a_start, double a_end, 
	double b, double c, double D, double sigma, double dt, double dx, int r,  int g)
{

	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	init_constframe(Rho_0, 0.25, g*g); //Returns Rho_0 with a full initial frame filled with 0.2.
	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "NOSTRA" <<endl;
  // The pspace to iterate over.
    vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  for (int i=0; i< t_measure.size(); i++)
  {
	cout << t_measure[i] << " ";
  } cout << endl;

  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = 5 * microToSeconds;     //5 seconds
    
    cout<<"Delay 1 in progress... ("<<delay1/microToSeconds<<"s)"<<endl;
	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

  	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );

	cout << "On initialisation, current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl;

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 25)
	{
		omp_set_num_threads(25); //Limiting use on Chunk.
	}

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
		**/
		expanded_percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_space[i], b, c, D, sigma, dt, dx, r, g);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting critical exponent data for each grid size in order.
          stringstream message3;
		  size_t currSize = getCurrentRSS( ); //Check Memory usage.
		  size_t peaksSize    = getPeakRSS( );
			 
          message3 << "Is this happening?\n" << "In replicate " << i << ", current size in MB: " << currSize/(1024.0*1024.0) 
		  << " and Peak Size (in MB): " << peaksSize/(1024.0*1024.0) << endl;
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, coco, tm ,d3, p1, p2, rini, Dm, Sm, dix, dimitri;

  L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a_start; p2 << setprecision(4) << a_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma; coco << setprecision(4) << c; dix << setprecision(2) << dx;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_1stdp;
  // Creating a file instance called output to store output data as CSV.
	output_1stdp.open("../Data/SPDP/1stOrder_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_c_"+ coco.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	cout << "Save file name: " <<endl;
	cout << "../Data/SPDP/1stOrder_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_c_"+ coco.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv";

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
	output_1stdp << " a , c , L,  t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";
	cout << "The vector elements are: "<< endl;
  	cout << "  a , c , L, t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_1stdp << setprecision(8) << vec[i][0] << "," << vec[i][1] << "," << vec[i][2] << "," << vec[i][3] 
		<< "," << setprecision(16) << vec[i][4] << "," << setprecision(16) << vec[i][5] << "," <<
		 vec[i][6] << "," << setprecision(16) << vec[i][7] << "," << vec[i][8] <<  endl;
		if( i%(1000) ==1)
    {
			cout << setprecision(5) << vec[i][0] << "," << vec[i][1] << "," << vec[i][2] << "," << vec[i][3] << "," 
			<<  setprecision(9) << vec[i][4] << "," << setprecision(9) << vec[i][5] << "," << vec[i][6] << "," << vec[i][7] 
			<< "," << vec[i][8] << endl;
    }
	}
	output_1stdp.close();
}

void first_order_terminal_diff(int div, double t_event, double p_start, double p_end, double a,
	double b, double c, double D, double sigma, double dt, double dx, int r,  int g)
{

	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//init_constframe(Rho_0, 0.25, g*g); //Returns Rho_0 with a full initial frame filled with 0.2.
	vector<double> p_space = linspace(p_start, p_end, div);
	cout << "NOSTRA" <<endl;
  // The pspace to iterate over.
    vector <double> t_measure; t_measure.push_back(t_event);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  for (int i=0; i< t_measure.size(); i++)
  {
	cout << t_measure[i] << " ";
  } cout << endl;

  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = 5 * microToSeconds;     //5 seconds
    
    cout<<"Delay 1 in progress... ("<<delay1/microToSeconds<<"s)"<<endl;
	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

  	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );

	cout << "On initialisation, current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl;

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 33)
	{
		omp_set_num_threads(33); //Limiting use on Chunk.
	}

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;
	  double Rho_0[g*g];

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < p_space.size(); i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on initial Rho Value:\t" << p_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
		
		init_constframe(Rho_0, p_space[i], g*g); //Returns Rho_0 with a full initial frame filled with p_space[i].

        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
		**/
		terminal_percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_event, a, b, c, D, sigma, dt, dx, r, g);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting critical exponent data for each grid size in order.
          stringstream message3;
		  size_t currSize = getCurrentRSS( ); //Check Memory usage.
		  size_t peaksSize    = getPeakRSS( );
			 
          message3 << "Is this happening?\n" << "In replicate " << i << ", current size in MB: " << currSize/(1024.0*1024.0) 
		  << " and Peak Size (in MB): " << peaksSize/(1024.0*1024.0) << endl;
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, coco, tm ,d3, p1, p2, rini, Dm, Sm, dix, aoili;

  L << g; tm << t_event; d3 << setprecision(3) << dt; p1 << setprecision(4) << p_start; p2 << setprecision(4) << p_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma; coco << setprecision(4) << c; dix << setprecision(2) << dx;
  aoili << setprecision(4) << a;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_1stdp;
  // Creating a file instance called output to store output data as CSV.
	output_1stdp.open("../Data/SPDP/1stOrder_Rep_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ aoili.str() + "_c_"+ coco.str() +
	"_p1_"+ p1.str() + "_p2_"+ p2.str() + "_dx_"+ dix.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	cout << "Save file name: " <<endl;
	cout << "../Data/SPDP/1stOrder_Rep_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ aoili.str() + "_c_"+ coco.str() +
	"_p1_"+ p1.str() + "_p2_"+ p2.str() + "_dx_"+ dix.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv";

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
	output_1stdp << "  a , c , L, r, t ,  <<Rho(x; 0)>_x>, <<Rho(x; t)>_x>, # Active Sites \n";
	cout << "The vector elements are: "<< endl;
  	cout << "  a , c , L, r, t ,  <<Rho(x; 0)>_x>, <<Rho(x; t)>_x>, # Active Sites \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_1stdp << setprecision(8) << vec[i][0] << "," << vec[i][1] << "," << vec[i][2] << "," << vec[i][3] 
		<< "," <<  vec[i][4] << "," << setprecision(16) << vec[i][5] << "," << setprecision(16) <<
		 vec[i][6] << "," << setprecision(16) << vec[i][7] <<  endl;
		if( i%(1000) ==1)
    {
			cout << setprecision(5) << vec[i][0] << "," << vec[i][1] << "," << vec[i][2] << "," << vec[i][3] << "," 
			<<  vec[i][4] << "," << setprecision(9) << vec[i][5] << "," << setprecision(9) <<  vec[i][6] << "," << vec[i][7] << endl;
    }
	}
	output_1stdp.close();
}

void frames_expanded_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, double Rh0[],
	 double t_max, double t_event, double a, double b, double c, double D, double sigma, double dt, double dx, int g, int r)
{
	// Find frame ay a given time t_event
	//Integrating first order PDE: arho - brho^2 - crho^3 +Diff + Stocj
	
	int nR2[g*g][3][3] ={0}; //nR2 stores neighbours of all sites in a ball of NORM 1 (SQUARE) radius 1.
	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_S8(nR2, g); //Assigning neigbours.
	int tot_iter = static_cast<int>(1.0/dt) + 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	double rho_rep_avg_var[tot_iter][5] ={0.0}; 
	//Stores time.
	double diff_coefficient = D/(dx*dx); //= D*(1.0/(6.0*dx*dx)); //Saves us some trouble later.
	double beta = a - 4*D/(dx*dx); //a - (10.0/3.0)*D/(dx*dx); //double beta = a - 4*D/(dx*dx);
	double lambda_exp = exp( (beta)*dt); double lambda = 2*(beta)/(sigma*sigma*(lambda_exp -1.0));
	double bdt = b*dt; double cdt = c*dt;
	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	int ind=0;
	for(double t =0; t <= 1; t+=dt)
	{ rho_rep_avg_var[ind][0] = t; ind+=1; }
	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[ind+i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep

	double alpha_prime = diff_coefficient*4; //Assume each Rho in neighbourhood for alpha estimation averages 1.
	double poiss_ru = lambda*lambda_exp*2.5; //Mean
	double mu_nought = 2.0*alpha_prime/(sigma*sigma) + poiss_ru;
	poiss_ru += 4*sqrt(poiss_ru); //Mean of Poisson Sampler + 4 SD.
	mu_nought += 6*sqrt(mu_nought); //Mean of  Mean of Gamma Sampler + 6 SD.

	double rho_cut_off =0.55*g*g; //Used as cut_off to capture frames.

	t_max = t_event + 50;


	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t Lambda:\t" << lambda << "\t Beta:\t " << beta << "\t dt,dx:\t " << dt << "," << dx << "\t Poisson Cut-Off:\t " << poiss_ru 
	<< "\t and  Gamma Cut-Off:\t " << mu_nought  <<endl; //cout << m1.str();
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	for(int j=0; j < r; j++)
	{
		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1 << "For Thread Rank:\t " << omp_get_thread_num() 
		<< "\t with Total Iterations:\t " << tot_iter  << "  in Replicate:\t" << j << endl; cout << m1.str();
		errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);

		//Defining the Dornic variables
		vector <double> Rho_dt(g*g, 0.0); //Updated through the turn, in turn used to update DRho
		vector <double> DRho(g*g, 0.0); //Dummy variable, Kept constant, only updated at end of turn.
		vector <double> DummyRho(g*g, 0.0); //Dummy variable, Kept constant, only updated at end of turn.
		
		double Rho_M[tot_iter][2] ={0.0}; //Stores # Active Sites, <rho>x  values per time step for a single replicate.

		//cout<< "Initially: "<<endl;
		for(int i=0; i< g*g; i++)
		{
			Rho_dt[i] = Rh0[i]; //Assigning initial conditions
			DRho[i] = Rh0[i];
		} //cout<<endl;

		poisson_distribution<int> poisson; gamma_distribution <double> gamma;

		double t=0; int index = 0;  //Initialise t
		int s=0; int po=1; int so =1; int lo=1; int glo =1;
		//int co=1; int eo=1; 
	
		

		while (t < t_max + dt)
		{

			double rhox_avg, rhox_num;
			//Start by updating the Rho vector.
			

			if(t <= 1.0 || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{
				rhox_num = occupied_sites_of_vector(Rho_dt, g*g);
				rhox_avg = mean_of_vector(Rho_dt, g*g); //Finds spatial average of densities at given t.
				//rhox_num = occupied_sites_of_vector(Rho_dt, g*g); //Finds number of occupied at given t.
				Rho_M[s][0] = rhox_num; Rho_M[s][1] = rhox_avg;  s+=1;//Rho_M.push_back({t, rhox_avg});
				//if(rhox_num == 0.0)
					//break; //Finish negotiations.
				if( t > 1.0)
					index+=1;
			}
			/**
			if(glo == 1 && rhox_num < rho_cut_off )
			{
				glo =0;
				//Dump frame.
				stringstream L, tm ,d3, p1, dX, Dm, Sm;

  				L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a;
  				dX << setprecision(3) << dx; Sm << setprecision(3) << sigma; Dm << setprecision(3) << D;
				ofstream output_frame;
  				// Creating a file instance called output to store output data as CSV.
				output_frame.open("../Data/SPDP/Frames/HFRAME_P_c_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a1_"+ p1.str() +
				"_Sig_"+ Sm.str() + "_dx_" + dX.str() + ".csv");
				output_frame << " a, c, t , x,  Rho(x; t)\n";
				for(int i=0; i< g*g; i++)
				{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |
					
					output_frame << a << ","<< c << "," << t << "," << i << "," << Rho_dt[i] <<endl;
				}
				output_frame.close();
			} */
			else if(t >= t_event -dt/2.0 &&  t < t_event +dt/2.0) //&& glo == 1 )
			{
				glo =0;
				//Dump frame.
				stringstream L, tm ,d3, p1, dX, Dm, Sm , repli;

  				L << g; tm << t_event; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; repli << r;
  				dX << setprecision(3) << dx; Sm << setprecision(3) << sigma; Dm << setprecision(3) << D;
				ofstream output_frame;
				std::string filename;
				filename = "../Data/SPDP/Frames/1RBOKCHOYFRAME_P_c_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a1_"+ p1.str() +
					"_Sig_"+ Sm.str() + "_dx_" + dX.str() + "_R_" +repli.str() + ".csv";
				if(j == 0)
					output_frame.open(filename, std::ios_base::out); //If first replicate initiate fresh CSV
				else
					output_frame.open(filename, std::ios_base::app); //Else append.
  				// Creating a file instance called output to store output data as CSV.
				output_frame << " a, c, t , x,  Rho(x; t)\n";
				for(int i=0; i< g*g; i++)
				{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |
					
					output_frame << a << ","<< c << "," << t << "," << i << "," << Rho_dt[i] <<endl;
				}
				output_frame.close();

			}

			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC.
			for(int i=0;i<Rho_dt.size();i++)
			{
				double alpha_i = diff_coefficient*(1*(DRho[nR2[i][0][1]] + DRho[nR2[i][2][1]] +
				DRho[nR2[i][1][0]] + DRho[nR2[i][1][2]]));


				if(alpha_i == 0 && DRho[nR2[i][1][1]] == 0) 
				{
					continue; 
				} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

				if(alpha_i < 0)
				{
					stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        			m6 << "In Replicate:\t" << j << " \tWe are NOT OKAY:\t" << Rho_dt[i] << "\t at index:  " << i << " at time:\t:" << t <<
					" with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i << endl; cout << m6.str();
					errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					//po=0;
				}
				/**
				double maxtek = max({DRho[nR2[i][0][1]],DRho[nR2[i][2][1]],DRho[nR2[i][1][0]], DRho[nR2[i][1][2]]}, maxis);
				if(maxtek - DRho[nR2[i][1][1]] > 5 || maxtek - DRho[nR2[i][1][1]] < -5)
				{
				//Checking if Laplacian is causing our troubles.
				
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "ALPHA GO-GO-TRON Rho[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " at time:\t:" << t <<
				" with Max(Neighbour(Rho)):   " << maxtek << "\t and alpha_i:\t" << alpha_i << endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				//po=0;
				}
				*/
				double mu = -1.0 + 2.0*alpha_i/(sigma*sigma);
				double ziggy = lambda*lambda_exp*Rho_dt[i]; double gru;
				if(ziggy == 0.0)
				{	//stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		//m6 << "SICK!!" <<endl; cout << m6.str(); 
				gru = mu + 1.0;  }
				/**
				else if(ziggy < 0.0 || isnan(ziggy) == true || isinf(ziggy) == true)
				{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "POISSON IS MESSED UP:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i << endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				} */
				else
				{
					poisson = poisson_distribution<int>(lambda*lambda_exp*Rho_dt[i]);
					gru = mu + 1.0 + poisson(rng);
				}

				gamma = gamma_distribution<double>(gru, 1.0);
					
	        	Rho_dt[i]= gamma(rng)/lambda;

				if(gru < 0 || isnan(gru) == true || isinf(gru) == true)
				{
					stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        			m6 << "GAMMA MESSED UP BIG TIME Rho*[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
					<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i <<  
					"\t and GRU:\t" << gru <<  "  in Replicate:\t" << j <<  endl; //cout << m6.str();
					errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				}

				if(gru > mu_nought || Rho_dt[i] > 5 || alpha_i > 20 || ziggy > poiss_ru)
				{
					stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        			m6 << "BLOWING UP, Rho*[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
					<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i <<  
					"\t and GRU:\t" << gru << "\t and Poisson RU:\t" << ziggy << "  in Replicate:\t" << j <<  endl; //cout << m6.str();
					errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				}
	        
	        
				/**
				if(isnan(Rho_dt[i]) == true || isinf(Rho_dt[i]) == true)
				{
					stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        			m6 << "YOU HIT ROCK BOTTOM WITH: Rho*[t,i]\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
					<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and alpha_i:\t" << alpha_i << 
					"\t and GAMMA ENTRY:\t" << gru << endl; cout << m6.str();
					errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				} */
				DummyRho[i] =Rho_dt[i];
			}

			//RK4 Integration of remaining terms

			RK4_Integrate(Rho_dt, a, b, c, D, t, dt, dx, g);

			for(int i=0;i<Rho_dt.size();i++)
			{
			//Euler Integration of remaining terms
			//Rho_dt[i] = (1 - (bdt + cdt*Rho_dt[i])*Rho_dt[i])*Rho_dt[i];
			//Rho_dt[i] = Rho_dt[i]/(1 + (dt*Rho_dt[i]*(c*Rho_dt[i] + b)));

			if(isnan(Rho_dt[i]) == true )
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RK4 BLEW IT Rho[t,i]:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and Analytic Integration Term:\t"
				<< DummyRho[i]<< "  in Replicate:\t" << j <<  endl; cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
				lo = -1;
			}

			if( Rho_dt[i] < 0 || isinf(Rho_dt[i]) == true and so==1)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RHO FALLS BELOW O:\t" << Rho_dt[i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[i] << "\t and Analytic Integration Term:\t" 
				<< DummyRho[i]<< "  in Replicate:\t" << j <<   endl; cout << m6.str(); //so=0;
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			DRho[i] =Rho_dt[i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
			}

			t+=dt; //Update timestep

			if(lo == -1)
			{
				exit(3); //Exit as soon as Rho_dt[i] becomes NaN to prevent infinite simulations.
			}

		} // End of t loops

		/**

		if(j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			for(int i =0; i< tot_iter; i++)
			{	rho_rep_avg_var[i][1] = Rho_M[i][1]; rho_rep_avg_var[i][2] = 0.0; rho_rep_avg_var[i][4] = Rho_M[i][0];
				if(Rho_M[i][1] > 0) //Ensuring only surviving runs are considered
				{ rho_rep_avg_var[i][3] = 1;  } //One surviving run as of now
				else if(Rho_M[i][1] == 0) //Ensuring only surviving runs are considered
				{ rho_rep_avg_var[i][3] = 0;  } //No surviving run as of now

			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			// Second or higher replicate, use incremental advances.
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		**/
		
		vector<double>().swap(Rho_dt); vector<double>().swap(DRho); //vector<vector <double>>().swap(Rho_M);
		// Remove dynamic vectors from memory and free up allocated space.

	} //End of r loop.

	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	/**
	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|	 c		|    L 		|    t 		|     <<Rho(t)>>x,r			|   #Active Sites | #Surviving Runs
		Rho.push_back({a,c, static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], rho_rep_avg_var[i][4], rho_rep_avg_var[i][3]});	
	}
	**/

	stringstream m10;
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m10.str();
	errout.open(thr, std::ios_base::app); errout << m10.str(); errout.close();
}



void capture_frame_decay(double Rho_0[], int div, double t_max, double t_event, double a_start, double a_end, 
  double b, double c, double D, double sigma, double dt, double dx, int r,  int g)
{

  //Get histograms for rho(x,t) for all x at a given t in the stationary phase. Use this to construct Prob density.

  //init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
  init_randframe(Rho_0, 1.00, 0.25, g*g); //Returns Rho_0 with a full initial frame filled with 0.2.
  //init_quarterframe(Rho_0, 1.0, 1.0, 1.0 , 0.0, g); //Returns Rho_0 with a initial frame filled with right side set to 1 and left to 0.

  cout << " A subset of the initial frame:\t" << endl;
  for(int i = 0; i < 10; i++)
  {
	for(int j =0; j < 10; j++)
		cout << setprecision(3) << Rho_0[i*g + j] << " ";
	cout << endl;
  }

  vector<double> a_space = linspace(a_start, a_end, div);
  cout << "NOSTRA" <<endl;
  // The pspace to iterate over.
  vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  for (int i=0; i< t_measure.size(); i++)
  {
	cout << t_measure[i] << " ";
  } cout << endl;

  std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

  size_t currentSize = getCurrentRSS( ); //Check Memory usage.
  size_t peakSize    = getPeakRSS( );

	cout << "On initialisation, current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl;

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 25)
	{
		omp_set_num_threads(25); //Limiting use on Chunk.
	}

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
		**/
		frames_expanded_percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, t_event, a_space[i], b, c, D, sigma, dt, dx, g, r);
        //Spit out frames at time "t_event"

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting critical exponent data for each grid size in order.
          stringstream message3;
		  size_t currSize = getCurrentRSS( ); //Check Memory usage.
		  size_t peaksSize    = getPeakRSS( );
			 
          message3 << "Is this happening?\n" << "In replicate " << i << ", current size in MB: " << currSize/(1024.0*1024.0) 
		  << " and Peak Size (in MB): " << peaksSize/(1024.0*1024.0) << endl;
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, coco, tm ,d3, p1, p2, rini, Dm, Sm, dix, dimitri;

  L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a_start; p2 << setprecision(4) << a_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma; coco << setprecision(4) << c; dix << setprecision(2) << dx;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.

  /**
	ofstream output_1stdp;
  // Creating a file instance called output to store output data as CSV.
	output_1stdp.open("../Data/SPDP/1stOrder_FRAME_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_c_"+ coco.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	cout << "Save file name: " <<endl;
	cout << "../Data/SPDP/1stOrderBOK_FRAME_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_c_"+ coco.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv";

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
	output_1stdp << " a , c , L,  t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";
	cout << "The vector elements are: "<< endl;
  	cout << "  a , c , L, t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r, # Surviving Runs, # Active Sites, # Surviving Runs Alt \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_1stdp << setprecision(8) << vec[i][0] << "," << vec[i][1] << "," << vec[i][2] << "," << vec[i][3] 
		<< "," << setprecision(16) << vec[i][4] << "," << setprecision(16) << vec[i][5] << "," <<
		 vec[i][6] << "," << setprecision(16) << vec[i][7] << "," << vec[i][8] <<  endl;
		if( i%(1000) ==1)
    {
			cout << setprecision(5) << vec[i][0] << "," << vec[i][1] << "," << vec[i][2] << "," << vec[i][3] << "," 
			<<  setprecision(9) << vec[i][4] << "," << setprecision(9) << vec[i][5] << "," << vec[i][6] << "," << vec[i][7] 
			<< "," << vec[i][8] << endl;
    }
	}
	output_1stdp.close(); */
}