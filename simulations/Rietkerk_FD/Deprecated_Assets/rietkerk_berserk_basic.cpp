#include "rietkerk_berserk_basic.h"
#include "Debug.h"


//-------------- Template definitions and declarations must occur in each individual file.

/**
template<int D, typename T> struct createVec : public vector<createVec<D - 1, T>> 
{
  static_assert(D >= 1, "Vector dimension must be > 0!");
  template<typename... Args> createVec(int n = 0, Args... args) : vector<createVec<D - 1, T>>(n, createVec<D - 1, T>(args...)) 
  {
  }
};
template<typename T> struct createVec<1, T> : public vector<T> 
{
  createVec(int n = 0, const T& val = T()) : vector<T>(n, val) 
  {
  }
};

**/


//--------------------Defining Primitives------------------------------------//

void increase_stack_limit(long long stack_size){

	//Source: https://stackoverflow.com/questions/2279052/increase-stack-size-in-linux-with-setrlimit
	//Credit: https://stackoverflow.com/users/253056/paul-r
	//Used with modification

	// Necessary because of complexity of integrations.

	const rlim_t kStackSize = stack_size * 1024LL * 1024LL;   // min stack size = 128 MB
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


void add_three(int a, int b, int c)
{
	int d = a + b +c; cout << d << endl;
}

void init_fullframe(D2Vec_Double &array, int Sp, int size)
{
	//Returns array with all elements initialised to 1
	for(int s=0; s< Sp; s++) 
	{
		for(int j=0; j< size; j++)
		{	array[s][j] =1 ;}
	}
}

void init_constframe(D2Vec_Double &array, int Sp, int size, double constant[])
{
	//Returns array with all elements initialised to a constant value: "const"
	for(int s=0; s< Sp; s++) 
	{
		for(int i=0; i< size; i++)
		{	array[s][i] =constant[s] ;}
	}
}


//Initialises a frame (D2Vec_Double vector) with variable percent of the elements randomly set to const c_high, others set to const val c_low.
void init_randconstframe(D2Vec_Double &array, int Sp, int size, double perc,  double c_high[], double c_low[] )
{
	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);
	//Returns array with all elements initialised to a constant value: "const"
	
	for(int i=0; i< size; i++)
	{	
		uniform_real_distribution<double> unif;
		unif = uniform_real_distribution<double>(0.0, 1.0);
		double r = unif(rng);
		for(int s=0; s< Sp; s++) 
		{
			if(r <= perc)
				array[s][i] =c_high[s];
			else
				array[s][i] =c_low[s];
		}
	}
}




void init_randframe(D2Vec_Double &array, int Sp, int size, double mean[], double sd[])
{
	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);
	//Returns array with all elements initialised to random values drawn from a Gaussian distribution with mean mu and SD sigma,
	// with all values 
	for(int s=0; s< Sp; s++) 
	{
		for(int i=0; i< size; i++)
		{	
		normal_distribution<double> norm;
		norm = normal_distribution<double>(mean[s], sd[s]);
		array[s][i] =norm(rng);
		if(array[s][i] < 0)
			array[s][i] = 0.0; //Values less than 0 are unrealistic and not permissible.
		}
	}
}

void init_quarterframe(D2Vec_Double &array, int Sp, int g, double c1[], double c2[], double c3[], double c4[])
{
	//Returns array with elements initialised to a constant values: "const1", "const2", "const3", "const4" in clockwise quadrants
	// starting from the top left.
	int l2 = static_cast<int>(g/2);
	for(int s=0; s< Sp; s++) 
	{
	for(int i=0; i< l2; i++)
	{
		int il2 = i + l2;	
		for(int j=0; j< l2; j++)
		{
			array[s][i*g + j] =c1[s]; //Filling up top-left quadrant
			array[s][i*g + l2 + j] =c2[s]; //Filling up top-right quadrant
			array[s][il2*g + l2 + j] =c3[s]; //Filling up bottom-right quadrant
			array[s][il2*g  + j] =c4[s]; //Filling up bottom-right quadrant
		}
			
	}	
	}
}

void init_solitarytear(D2Vec_Double &array, int Sp, int length)
{
	//Returns array with only central element initalised to 1, all others 0.
	for(int s=0; s< Sp; s++) 
	{
	for(int i=0; i< length*length; i++)
		array[s][i] =0;
	int lhalf = static_cast<int>(length/2.0);
	array[s][length*lhalf +lhalf] =1;
	}
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

void var_mean_incremental_surv_runs(double t_avg_var_rep_N[][Sp4_1], double X_curr[][Sp2], int size)
{
	/** X_curr contains time-series data (indexed as (N(t), <rho(t)>x)) from the current surviving replicate (given by r) with a length "size".
	 *  rep_avg_var stores the time(t), running density average and variance of X_curr, as well as N(t) and number of survivng replicates 
	 * for each time-step from previous surviving replicates ( <= r-1) respectively.
	 * In this function, X_curr is used to update the running avg and variance of rep_avg_var to account for current surviving replicate r(t).
	 * Note rep_avg_Var[t][3] stores number of surviving runs at t.
	*/
	double mean_prev[size];
	for(int s=0; s < Sp -2; s++)
	{
	for(int i=0; i<size; i++)
	{
		int r= t_avg_var_rep_N[i][4*s+ 5]+1; //Number of surviving replicates
		if(X_curr[i][2*s+ 0] == 0.0)
		{	//Non-surviving run. No-updates required. Also <rho(t)>x == 0, it remains 0 for t' > t.
			return;
		}
		else
		{	//Number of surviving run at time t.
			mean_prev[i] = t_avg_var_rep_N[i][4*s+ 3]; //Stores X_n-1 mean data as temp variable.
			if(t_avg_var_rep_N[i][4*s + 5] == 0.0)
			{	//First surviving run encountered at time t, encountered here.
				t_avg_var_rep_N[i][4*s+ 3] = X_curr[i][2*s+ 1]; t_avg_var_rep_N[i][4*s+ 4] = 0.0; 
				t_avg_var_rep_N[i][4*s+ 6] = X_curr[i][2*s+ 0];
			}
			else
			{	//Not the first surviving run, which means r = rep_avg_var[i][3]+1 >=2
				t_avg_var_rep_N[i][4*s+ 3] = (X_curr[i][2*s+1] + (r-1)*t_avg_var_rep_N[i][4*s+ 3])/r; //Note r>= 1.
				t_avg_var_rep_N[i][4*s+ 6] = (X_curr[i][2*s+0] + (r-1)*t_avg_var_rep_N[i][4*s+6])/r; //Note r>= 1.
				//Calculating and storing incremental variance (V^2_n) [For formula derivation, refer to notes.]
				t_avg_var_rep_N[i][4*s+4] = (t_avg_var_rep_N[i][4*s+3] - mean_prev[i])*(t_avg_var_rep_N[i][4*s+3] - mean_prev[i]) +
				(1/(r-1))*( (r-2)*t_avg_var_rep_N[i][4*s+4] + 
				(X_curr[i][2*s+1] - t_avg_var_rep_N[i][4*s+3])*(X_curr[i][2*s+1] - t_avg_var_rep_N[i][4*s+3]));
			}
			t_avg_var_rep_N[i][4*s+5] = r; //At the end, update number of succesful runs.
		}
	}
	}
	for(int s=Sp -2; s < Sp; s++)
	{
		for(int i=0; i<size; i++)
		{
		//For soil water and surface water.
		int r= t_avg_var_rep_N[i][4*s+ 5]+1; //Number of surviving replicates
		if(t_avg_var_rep_N[i][4*s + 5] == 0.0)
				t_avg_var_rep_N[i][s] = X_curr[i][2*s+ 1]; //First surviving run encountered at time t, encountered here.
		else
		{
			t_avg_var_rep_N[i][s] = (X_curr[i][2*s+1] + (r-1)*t_avg_var_rep_N[i][s])/r; //Note r>= 1.
		}
		
		t_avg_var_rep_N[i][4*s+5] = r; //At the end, update number of succesful runs.
		}
	}

}

// --------------------------- Partitioning functions to make life easier -----------------------------------//

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

//----------------------------- Misc. Supporting Add-ons ------------------------------------------------- //

int theborderwall(D2Vec_Double &Rho_t, int g)
{
	//Checks if the border sites are occupied.
	for(int s=0; s <Sp; s++)
	{
		int g2 = g*(g-1);
		for(int i =0; i< g; i++)
		{
			if(Rho_t[s][i] > 0 || Rho_t[s][g*i] > 0  || Rho_t[s][g2 + i] > 0 || Rho_t[s][g*i + g - 1 ] > 0)
				return 1; //Top, left, bottom and right borders are occupied by active sites.
		}
		
	}
	return 0; //No active sites detected at the edges.
}

void determine_neighbours_R2( int g, D3Vec_Int &neighbours_R2)
{
	//Stores neighbouring sites (Van-Neumann Radius of 2) of each cell site in a 
	// S*(g^2)*2*5 matrix.
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

void determine_neighbours_Sq8(int g, D3Vec_Int &neighbours_Sq8)
{
    
	//Stores neighbouring sites of species s and at site i in a square (Norm 1) radius of 1

    // Neighbours_Sq8 has dimensions S*(L*L)*3*3 where S is the number of species (=2) in this case
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
				//Recall, for given s and i first (third) row stores vertical neighbours, second (fourth) row stores horizontal neighbours.
                neighbours_Sq8[i][y][x] = a*g + b; 
			}	
		}
	}

}


//----------------------------- Stochastic 2PAC Dornic Integration Machinery --------------------------------------------//


// ALT INTEGRATION MACHINERY WHERE MIXING OF TIMESTEPS DUE TO LAPLACIAN MAY BE AN ISSUE.



void f_2D(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, double alpha, double d, double rW, double W0,  
	double D[], double K[],  double t, double dt, double dx2, int g)
{
	//Vector function that updates an array containing (dP/dt, dW/dt, DO/dt) for each site in the lattice.
    //Based on the Rietkerk model for plant vegetation dynamics.

	for(int i=0; i < g*g; i++)
	{
        //Equations for the density of plants, soil water and surface water at each site.
        //Note that the Laplacian is calculated using reflective boundary conditions.
        //The Laplacian is calculated using the 5-point stencil method.
        //Equivalent to dP/dt = c*g_max*P*W/(W+K1) - d*P + D*(Laplacian of P)
        f[0][i] = c*gmax*Rho_M[1][i]*Rho_M[0][i]/(Rho_M[1][i] +K[1]) -d*Rho_M[0][i] + (D[0]/(dx2))*(Rho_M[0][nR2[i][0][1]]  + Rho_M[0][nR2[i][2][1]]  
        + Rho_M[0][nR2[i][1][0]]  + Rho_M[0][nR2[i][1][2]]  - 4*Rho_M[0][i]);
        //Equivalent to dW/dt = alpha*(P+K2*W)/(P+K2) - rW*W + D*(Laplacian of W)
        f[1][i] = alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[2][i] -gmax*Rho_M[1][i]*Rho_M[0][i]/(Rho_M[1][i] +K[1]) - rW*Rho_M[1][i] 
        + (D[1]/(dx2))*(Rho_M[1][nR2[i][0][1]]  + Rho_M[1][nR2[i][2][1]]  + Rho_M[1][nR2[i][1][0]]  + Rho_M[1][nR2[i][1][2]]  - 4*Rho_M[1][i]);
        //Equivalent to dO/dt = a - alpha*(P+K*W)/(P+K)*O + D*(Laplacian of O)
		f[2][i] = a - alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[2][i] + (D[2]/(dx2))*(Rho_M[2][nR2[i][0][1]]  + Rho_M[2][nR2[i][2][1]] 
        + Rho_M[2][nR2[i][1][0]]  + Rho_M[2][nR2[i][1][2]]  - 4*Rho_M[2][i]);
		
	}
	
	

}


//Modify the function below to perform RK4 integration of the 2D Vector Rho_t, with a laplacian and time step of dt and spatial step dx.
//The function should also take in the parameters of the model, and the time t at which the integration is being performed.
//The function should also take in the number of grid points g, and the number of species Sp.
//The function should also take in the number of threads being used, and the thread number.
//The function should also take in the number of iterations performed so far, and the number of iterations to be performed.


void RK4_Integrate(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D3Vec_Int &nR2, double a, double c, double gmax, double alpha, double d, double rW,
	  double W0, double D[], double K[],  double t, double dt, double dx, int g)
{
    

    // Rho_t stores value of rho at time t.

	// Rho_tstar stores value of rho*, provided by sampling Gaussian.
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	double dt6 = dt/6.0; double dt2 = dt/2.0; double dx2 = dx*dx;

    D2Vec_Double K1(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K2(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double K3(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K4(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double Rho_M(Sp, vector<double> (g*g, 0.0));

	//double K1[Sp][g*g] ={0.0}; double K2[Sp][g*g] ={0.0}; double K3[Sp][g*g] ={0.0}; double K4[Sp][g*g] ={0.0};
	//double Rho_M[Sp] ={0.0, 0.0};


	f_2D(K1, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t, dt, dx2, g); //K1 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K1[s][i];
	}

	f_2D(K2, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt2, dt, dx2, g); //K2 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K2[s][i];
	}
	f_2D(K3, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt2, dt, dx2, g); //K3 updated.
	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt)*K3[s][i];
	}
	f_2D(K4, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt, dt, dx2, g); //K4 updated.
    
	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
		{
			Rho_t[s][i]+= (dt6)*( K1[s][i] + 2.0*K2[s][i] + 2.0*K3[s][i] + K4[s][i]);
	

			if( Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || isnan(Rho_t[s][i]) == true)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RK4 WAS KO'ED WITH:\t" << Rho_t[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " For Species:\t:" << s << " with K1[s][i]:   " << K1[s][i] << "\t, K2[s][i]:\t" << K2[s][i] 
				<< "\t, K3[s][i]:\t" << K3[s][i] << "\t AND K4[s][i]:\t" << K4[s][i] << endl; //cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			
		}
	}
}

void RK4_Wrapper_2D(D2Vec_Double &Rho, vector <double> &t_meas, D2Vec_Double &Rh0, double t_max, double a, double c, double gmax, 
	  double alpha, double d, double rW, double W0, double D[], double K[], double a_st, double a_end, double dt, double dx, double dP, int r, int g)
{
	// Rho_0 stores value of rho at time t=0.
	// Rho_t stores value of rho at time t.
	// Rho_tstar stores value of rho*, provided by sampling Gaussian.
	// Rho_M stores value of rho at intermediate time steps.

	//D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0));
	//double Rho_tsar[Sp][g*g] ={0.0};
	//double Rho_M[Sp][g*g] ={0.0, 0.0};

	D3Vec_Int nR2 (g*g, D2Vec_Int (3, vector<int> (3, 0)));
    //nR2 is 3D [(L*L)x3x3] vector initialised to 0, which stores neighbours of all sites i (second dim) in a ball of NORM 1 (SQUARE) radius 1 (2ND & 3RD dim).
	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_Sq8(g, nR2); //Assigning neigbours.
	int tot_iter = static_cast<int>(1.0/dt) + 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0}; 

	int ind=0;
	for(double t =0; t <= 1; t+=dt)
	{ rho_rep_avg_var[ind][0] = t; ind+=1; }
	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[ind+i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep

	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t a:\t" << a << "\t c:\t " << c << "\t dt,dx:\t " << dt << "," << dx << "\t gmax:\t " << gmax 
	<< "\t and  Alpha:\t " << alpha  <<endl; //cout << m1.str();
	cout << m0.str();
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	
	
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		D2Vec_Double Rho_dt(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		double Rho_M[tot_iter][Sp2] ={0.0}; //Stores <rho>x values per time step for a single replicate.

		for(int i=0; i < g*g; i++)
		{
			for(int s= 0; s <Sp; s++)
			{
				Rho_dt[s][i] = Rh0[s][i]; //Assigning initial conditions
				Rho_tsar[s][i] = Rh0[s][i];
			}
			
		}

		double t=0; int index = 0;  //Initialise t
		int iter_index=0; int po=1; int so =1; int lo=1;

		while (t < t_max + dt)
		{


			if(t <= 1.0 || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{

				double rhox_avg, rhox_num; //Stores spatial average and number of occupied sites at given t.
                for(int s=0; s<Sp; s++)
                {
                    vector <double> temp= {Rho_dt[s].begin(),Rho_dt[s].end()}; //Rho_dt for species 's'
                    rhox_avg = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
				    rhox_num = occupied_sites_of_vector(temp, g*g); //Finds number of occupied at given t.
				    Rho_M[iter_index][2*s] = rhox_num; Rho_M[iter_index][2*s +1] = rhox_avg; 

					vector<double>().swap(temp); //Flush temp out of memory.0

                    if(rhox_num == 0.0 and s== 0)
					    break; //Finish negotiations as all basal species are dead.
                }
				 iter_index+=1;//Rho_M.push_back({t, rhox_avg});

				if( t >= 2000 && t <= 20000)
				{
					//Basically indicates we are at t_max. SAVE FRAME AT T_MAX FOR FUTURE ANALYSIS.
					
					stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, alph, w0t, dix, dimitri;

  					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  					rini << j; Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; gm << setprecision(3) << gmax;
					w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; dimitri << setprecision(3) << dP;
					a1 << a_st; a2  << a_end;
					// Three replicates are over.

					ofstream frame_dp;
  					// Creating a file instance called output to store output data as CSV
					frame_dp.open("../Data/Rietkerk/Frames/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_Basic_P_c_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() 
					+ "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() + "_dx_"+ dix.str() + "_alph_"+ alph.str() + "_gmax_"+ gm.str() + "_R_"+ rini.str() + ".csv");

					// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|
					frame_dp << "a_c,  x,  P(x; tmax), W(x; tmax), O(x; tmax) \n";
					for(int i=0; i< g*g; i++)
					{
						frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] <<endl;
				
					}

					frame_dp.close();

					
				}
				
				if( t > 1.0)
					index+=1;

			}

			for(int s= 0; s <Sp; s++)
				Rho_tsar[s] = Rho_dt[s]; //Assigning initial conditions

			RK4_Integrate(Rho_dt, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t, dt, dx, g);

			for(int s= 0; s <Sp; s++)
			{
				for(int i=0; i<g*g; i++)
				{
					if( Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true and so==1)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        				m6 << "RHO FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << Rho_tsar[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}

				}
				
			}
			

			t += dt;

			

		} //End of Time Loop


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
			{	
				for(int s =0; s< 1; s++)
				{
					rho_rep_avg_var[i][4*s+ 3] = Rho_M[i][2*s + 1]; // Avg density of frame at given time point i.
					rho_rep_avg_var[i][4*s+ 4] = 0.0; // Variance in avg frame density.
					rho_rep_avg_var[i][4*s+ 6] = Rho_M[i][2*s]; // Number of active sites in frame at given time point i.

					if(Rho_M[i][2*s] > 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s+5 ] = 1;  } //One surviving run as of now
					else if(Rho_M[i][2*s] == 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s + 5] = 0;  } //No surviving run as of now
				}
				rho_rep_avg_var[i][1] = Rho_M[i][2*1 + 1]; // Avg density of frame at given time point i.
				rho_rep_avg_var[i][2] = Rho_M[i][2*2 + 1]; // Avg density of frame at given time point i.
				
				

			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			// Second or higher replicate, use incremental advances.
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		if((j+1)%3 == 1)
		{	//Start logging at every multiple of 3
			stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, Dm, gm;

  			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dimitri << dP;
  			rini << j; Dm << setprecision(3) << D[0]; gm << gmax; 
			// Three replicates are over.

			ofstream output_dp;
  			// Creating a file instance called output to store output data as CSV.
			output_dp.open("../Data/Rietkerk/Prelims/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/PRELIM_AGGRAND_Basic_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
			"_D0_"+ Dm.str() + "_gmax_"+ gm.str() + "_R_"+ rini.str() + ".csv");

			// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
			output_dp << " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs, # Active Sites, \n";
			for(int i=0; i< tot_iter; i++)
			{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |

				output_dp << a << ","<< j << "," << g << "," << rho_rep_avg_var[i][0] << "," << rho_rep_avg_var[i][1] << "," << rho_rep_avg_var[i][2] << ",";
				for(int s=0; s <1; s++)
				{ output_dp << rho_rep_avg_var[i][4*s + 3] << "," <<
				rho_rep_avg_var[i][4*s +4] << "," << rho_rep_avg_var[i][4*s +5] << "," << rho_rep_avg_var[i][4*s +6] << ","; 
				}
				output_dp <<  endl;
				
			}

			output_dp.close();
		}




	} //End of Replicate Loop


	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	
	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|    L 		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |

		Rho.push_back({a, static_cast<double>(r), static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], 
					rho_rep_avg_var[i][2]});
		for(int s=0;s < Sp-2; s++)
		{

			
			Rho[i].insert(Rho[i].end(), { rho_rep_avg_var[i][4*s +3], rho_rep_avg_var[i][4*s +4], 
				rho_rep_avg_var[i][4*s +5], rho_rep_avg_var[i][4*s +6]}); 		
			/**if(s == 0)
			{	//Basal species
				Rho.push_back({a, c, static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][4*s + 1], 
					rho_rep_avg_var[i][4*s +2], rho_rep_avg_var[i][4*s + 3], rho_rep_avg_var[i][4*s + 4]});
			}
			else
			{ 	//Extend rows, don't create them.
				Rho[i].insert(Rho[i].end(), { rho_rep_avg_var[i][4*s +1], rho_rep_avg_var[i][4*s +2], 
				rho_rep_avg_var[i][4*s +3], rho_rep_avg_var[i][4*s +4]}); 
			}**/
		}
	}
	

	stringstream m10;
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m10.str();
	errout.open(thr, std::ios_base::app); errout << m10.str(); errout.close();
}



void first_order_critical_exp_delta(D2Vec_Double &Rh0, int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double d, double rW, double W0,  double D[], double K[], double dt, double dx, double dP, int r,  int g)
{

	
	double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2
	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//double mean[Sp] = {1.0, 0.05}; double sd[Sp] = {0.25, 0.0125};
	//init_randframe(Rho_0, Sp,  g*g); //Returns Rho_0 with a full initial frame filled with 0.2.
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
  	} 	cout << endl;

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

	double perc = 0.015; double c_high[Sp] ={dP, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

  	  
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {

		init_randconstframe(Rh0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho1(t)>>x,r			|    Var[<Rho1(t)>x],r    |
		**/

		RK4_Wrapper_2D(CExpRho_a, t_measure, Rh0, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , a_start, a_end, dt, dx, dP, r, g);
		//expanded_percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_space[i], b, c, D, sigma, dt, dx, r, g);
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

	stringstream L, coco, tm ,d3, p1, p2, rini, Dm0, Dm1, gm, alph, dix, dimitri;

	// Creating a string stream instance to store the values of the parameters in the file name.
	// This is done to avoid overwriting of files.

	//Check recursively if a folder exists, if not create one.
	//https://stackoverflow.com/questions/12510874/how-to-check-if-a-directory-exists
	/**
	struct stat info;
	if( stat( "../Data/Rietkerk/Frames", &info ) != 0 )
	{
		cout << "Cannot access ../Data/Rietkerk/Frames. Creating directory." << endl;
		const int dir_err = system("mkdir ../Data/Rietkerk/Frames");
		if (-1 == dir_err)
		{
			printf("Error creating directory!n");
			exit(1);
		}
	}*/

	

  	L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << a_start; p2  << a_end; rini << r; alph << alpha;
  	Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; gm << gmax; coco << setprecision(4) << c; dix << setprecision(2) << dx;
  	// setprecision() is a stream manipulator that sets the decimal precision of a variable.

	//Input a folder name with user-defined numbers a_st and a_end. Check recursively if this folder exists, if not create one.
	//https://stackoverflow.com/questions/12510874/how-to-check-if-a-directory-exists

	/**
	stringstream foldername;
	foldername << "../Data/Rietkerk/Frames/" << p1.str() << "-" << p2.str() << "/";
	// Creating a string stream instance to store the values of the parameters in the file name.
	// This is done to avoid overwriting of files.

	struct stat info2;
	if( stat( foldername.str().c_str(), &info2 ) != 0 )
	{
		cout << "Cannot access " << foldername.str() << ". Creating directory." << endl;
		const int dir_err = system(("mkdir " + foldername.str()).c_str());
		if (-1 == dir_err)
		{
			printf("Error creating directory!n");
			exit(1);
		}
	}
	*/

	ofstream output_1stdp;
  // Creating a file instance called output to store output data as CSV.
	output_1stdp.open("../Data/Rietkerk/1stOrderCC_Rietkerk_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + "_gmax_"+ gm.str() + "_alpha_"+ alph.str() + "_R_"+ rini.str() + ".csv");

	cout << "Save file name: " <<endl;
	cout << "../Data/Rietkerk/1stOrderCC_Rietkerk_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + "_gmax_"+ gm.str() + "_alpha_"+ alph.str() + rini.str() + ".csv";

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

	string header = " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs, # Active Sites";;
	output_1stdp << header << "\n";
	cout << "The vector elements are: "<< endl;
  	cout << header << "\n";

	for(int i=0; i< vec.size(); i++)
	{
		for(int j=0; j< vec[i].size() -1; j++)
			output_1stdp << setprecision(16) << vec[i][j] << ",";
		output_1stdp << setprecision(16) << vec[i][vec[i].size() -1] <<  endl;
		if( i%(1000) ==1)
    	{
			for(int j=0; j< vec[i].size() -1; j++)
				cout << setprecision(16) << vec[i][j] << ",";
			cout << setprecision(16) << vec[i][vec[i].size() -1] <<  endl;
    	}
	}
	output_1stdp.close();
}