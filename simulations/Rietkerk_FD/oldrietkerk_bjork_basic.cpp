#include "rietkerk_bjork_basic.h"
#include "Debug.h"

double Mm; // Density dependent mortality of predator.

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

template <typename T> int sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}

template <typename T> int sgn_index(T val) 
{
    if(T(0) < val)
		return 0;
	else
		return 1;
	//Template designed specifically for use in determining the index for the UFDM Scheme.
	//Specifically, the notation follows entirely from determine_neighbours_Sq4().
	//Left (or up) neighbour (index 0) if val > 0, right (or down) neighbour (index 1) if val < 0.
}


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

bool comparePairs(const pair<int, int>& a, const pair<int, int>& b) 
{
    return (a.first * a.first + a.second * a.second) < (b.first * b.first + b.second * b.second);
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

//Initialises a frame (D2Vec_Double vector) with variable percent of the elements randomly set to const c_high, others set to const val c_low.
// The first Sp elements of c_high and c_low are used if R < R_c, otherwise the last Sp elements are used.
void init_randbistableframe(D2Vec_Double &array, int size, double R, double R_c, double perc,  double c_high[], double c_low[])
{

	mutex errorMutex; // Mutex to make error messages thread-safe
	if(Sp != 5)
	{
		std::lock_guard<std::mutex> lock(errorMutex);
		std::cerr << "Error: Number of species must be 3 for this function." << std::endl;
	}
	//c_high and c_low are arrays of size 2*Sp, with the values of the constants for each species.
	// If R < R_c, then the first Sp elements of c_high and c_low are passed to init_randconstframe() to initialise the frame.
	// If R >= R_c, then the last Sp elements of c_high and c_low are passed to init_randconstframe() to initialise the frame.
	if(R < R_c)
	{
		init_randconstframe(array, Sp, size, perc, c_high, c_low);
	}
	else
	{
		init_randconstframe(array, Sp, size, perc, c_high + Sp, c_low + Sp);
	}
	
}

void init_csvconstframe(D2Vec_Double &array, D2Vec_Double &const_ind_val, const std::string& filename, const vector<int> &columns, int size) 
{	
	/**
	Initialises a frame (stored finally in D2Vec_Double array) with constant values from a CSV file.
	The columns derived from the CSV file are stored in the vector<int> columns, and generally include 
	V(x,t), W(x,t) and O(x,t). The inital values of other higher order trophic species are stored in 
	the const_ind_val vector (which stores tuples of species index and initial value).
	NOTE:
	The CSV file is expected to have a header line starting with "a_c", and the columns are expected to be comma-separated.
	*/

    // Create an empty 2D vector with Sp rows (and no columns) called result storing V(x,t), W(x,t) and O(x,t) in its respective columns.
    D2Vec_Double result(3, vector<double>(1, 0.0));

	mutex errorMutex; // Mutex to make error messages thread-safe

    // Open the CSV file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) 
    {
		std::lock_guard<std::mutex> lock(errorMutex);
        std::cerr << "Error opening file: " << filename << std::endl;
    }

    std::string line;
    int columnIndex = columns[columns.size()-1];
    while (std::getline(inputFile, line)) 
    {
        std::stringstream ss(line);
        std::string cell;
		
		//Check if ss starts with "a_c", if so, skip the line.
		if(ss.str().find("a_c") != std::string::npos)
			continue; //Removing header line, if it exists, from CSV file.

        // Extract the desired columns
        for (int i = 0; i <= columnIndex; ++i) 
        {
            if (!std::getline(ss, cell, ',')) 
            {
				lock_guard<mutex> lock(errorMutex);
                std::cerr << "Error: Column index " << columnIndex << " out of bounds." << std::endl;
                inputFile.close();  
            }
            //Check if i is in columns vector
            if (std::find(columns.begin(), columns.end(), i) != columns.end()) 
            {
                // Store the cell value in the result vector as double
                double value;
                try 
                {
                    value = std::stod(cell);
                    result[i-2].push_back(value);
                } 
                catch (const std::invalid_argument& e) 
                {
					std::lock_guard<std::mutex> lock(errorMutex);
                    std::cerr << "Error: Unable to convert cell value to double. Column index: " << i
					<< " , for file: " << filename << " on thread #" << omp_get_thread_num() 
					<< "\n with the offending line: \n" << ss.str() << std::endl;
					result[i-2].push_back(0);
                } 
                catch (const std::out_of_range& e) 
                {
					std::lock_guard<std::mutex> lock(errorMutex);
                    std::cerr << "Error: Cell value out of range for double. Column index: " << i 
					<< " , for file: " << filename << " on thread #" << omp_get_thread_num() 
					<< "\n with the offending line: \n" << ss.str() << std::endl;
					result[i-2].push_back(0);
                }
            }
        }

    }

    // First store result in array.

    for(int s=0; s< 1; s++)
    {
        for(int i=0; i< size; i++)
        {
            array[s][i] = result[s][i+1];
        }
    }

    for(int i=0; i< size; i++)
    {
            array[Sp-2][i] = result[1][i+1];    // Storing W in the second last row of array
            array[Sp-1][i] = result[2][i+1];    // Storing O in the last row of array
    }

    // const_ind_val is a 2D vector of size (Sp - columns.size()) x 2, which stores the constant index 
    //for each column not being extracted (and thus not in result) and the constant value to be substituted in the array at that index value.
    // These values are pre-assigned in the main function.

    for(int i=0; i< const_ind_val.size(); i++)
    {
        int index = const_ind_val[i][0];
        double value = const_ind_val[i][1];

        for(int j=0; j< size; j++)
        {	array[index][j] = value;	}
    }
    // Close the file
    inputFile.close();

    //Free up memory allocated to result
    vector <vector<double>>().swap(result);
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

void init_gaussframe(D2Vec_Double &array, int size, vector<double> &sd, vector<double> &amp, int s_gauss /* = SpB*/)
{
	/**
	 * Accepts a 2D vector array of size Sp x g^2, and fills part of it with a 2D Gaussian distribution with mean (cx, cy) and standard deviation sd (and amplitude amp).
	 * The other part of the array is filled with constant values provided in the array amp.
	 * The split occurs at the index s_gauss (optional argument, default is SpB).
	 * The way the array is filled is as follows:
	 * 1. The first s_gauss species are filled with a 2D Gaussian distribution with mean (cx, cy) and standard deviation sd (and amplitude amp).
	 *    This occurs through the function init_gaussiantear(). cx. cy and sd MUST be arrays of size s_gauss.
	 * 2. cx and cy are arrays of size s_gauss, and these values are filled with the Cartesian coordinates of the centre of the Gaussian distribution for each species.
	 *    These values are chosen in a counter-clockwise manner starting from the top-left corner of the array (then top-right, bottom-right, bottom-left and so on depending on s_gauss).
	 * 3. sd[] is an optional arguement that defaults to an array of size SpB with elements set to  [L/8.0, L/16.0, L/16.0 ....].
	 *    If the size of sd[] is less than s_gauss, the remaining elements are set to the last value of sd[]. 
	 * 	  If the size of sd[] is greater than s_gauss, the extra elements are ignored.
	 * 4. amp[] is generally an array of size Sp, with the values of the constants for each species. The first s_gauss elements are used to fill the Gaussian distribution via init_gaussiantear().
	 *    If amp[] is not of size Sp, the function will provide a warning and fill/ignore the elements accordingly (if size < Sp, the remaining elements are set to the last value of amp[]).
	 * 5. The remaining Sp-s_gauss species are filled with the constant values provided in the array amp.
	 * 6. The function returns the filled array.
	 * 7. It does so by designating two temporary D2Vec_Double arrays, gauss_frame and const_frame, to store the Gaussian and constant values respectively.
	 */
	mutex errorMutex; // Mutex to make error messages thread-safe
	int g = static_cast<int>(sqrt(size));
	//Initialise the Gaussian frame
	D2Vec_Double gauss_frame(s_gauss, vector<double>(g*g, 0.0));
	// Make cx, cy as dynamic arrays of size s_gauss.
	double *cx = new double[s_gauss]; double *cy = new double[s_gauss];
	// Fill cx and cy arrays with Cartesian coordinates of the centre of the Gaussian distribution for each species, starting from the top-left corner.
	for(int i=0; i< s_gauss; i++)
	{
		// cx and cy reprsent Cartesian coordinates in a counter-clockwise manner starting from the top-left corner of the array.
		// i.e cx =[ 0, 1, 1, 0, 0, 1 ....] and cy = [0, 0, 1, 1, 0, 0 ....]
		cx[i] = (i%4)%3 == 0 ? int(g/8) : int(7*g/8);
		cy[i] = (i%4) < 2 ? int(g/8) : int(7*g/8);
	}
	//If size of sd[] is NOT s_gauss, fill the remaining elements with the last value of sd[], else ignore the extra elements.

	if(sd.size() < s_gauss)
	{
		double last_sd = sd[sd.size() - 1];
		for(int i=sd.size(); i< s_gauss; i++)
			sd.push_back(last_sd);
	}
	else if(sd.size() > s_gauss)
	{
		sd.resize(s_gauss);
	}
	
	init_gaussiantear(gauss_frame, s_gauss, g, cx, cy, sd, amp);

	//Initialise the constant frame.
	D2Vec_Double const_frame(Sp - s_gauss, vector<double>(g*g, 0.0));
	for(int i=0; i< Sp - s_gauss; i++)
	{
		for(int j=0; j< g*g; j++)
		{
			const_frame[i][j] = amp[s_gauss + i];
		}
	}

	//Combine the two frames to get the final frame.
	for(int s=0; s< s_gauss; s++) {
		for(int i=0; i< g*g; i++)
			array[s][i] = gauss_frame[s][i];
		
	}
	for(int s=0; s< Sp - s_gauss; s++) {
		for(int i=0; i< g*g; i++)
			array[s + s_gauss][i] = const_frame[s][i];	
	}

	//Debugging.
	std::lock_guard<std::mutex> lock(errorMutex);
	ostringstream oss;
	oss << "Gaussian frame initialised with " << s_gauss << " species.\n";
	oss << "Amp values:\n";
	for(int i=0; i< Sp; i++)
		oss << amp[i] << " ";
	oss << "\n";
	oss << "SD values:\n";
	for(int i=0; i< s_gauss; i++)
		oss << " ( " << i << " , " << sd[i] << " ) \t";
	oss << "\n";
	oss << "cx values & cy values:\n";
	for(int i=0; i< s_gauss; i++)
		oss << i << " " << cx[i] << " " << cy[i] << "\n";

	cout << oss.str();
	cerr << oss.str() << endl;

	// Swap memory to free up space.
	vector <vector<double>>().swap(gauss_frame);
	vector <vector<double>>().swap(const_frame);
	vector <double>().swap(amp);
	// Delete cx and cy arrays.
	delete[] cx; delete[] cy;
}

void init_gaussiantear(D2Vec_Double &array, int Sp_lim, int g, double cx[], double cy[], const vector<double> &sd,const vector<double> &amp)
{
	/**
	 * Returns array with elements initialised to a 2D Gaussian tear drop with mean (cx, cy) and standard deviation sd.
	 * array is a 2D vector of size Sp x g^2. cx, cy, amp & sd are arrays of size Sp, with the values of the constants for each species.
	 */
	for(int s=0; s< Sp_lim; s++) 
	{
		double A = amp[s]/(2*PI*sd[s]*sd[s]);
		for(int i=0; i< g; i++)
		{
			for(int j=0; j< g; j++)
				array[s][i*g + j] = amp[s]*exp(-((i-cx[s])*(i-cx[s]) + (j-cy[s])*(j-cy[s]))/(2*sd[s]*sd[s]));
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

void var_mean_incremental_surv_runs(D2Vec_Double &t_avg_var_rep_N, const D2Vec_Double &X_curr, int size, int j)
{
	/** X_curr contains time-series data (indexed as (N(t), <rho(t)>x)) from the current surviving replicate (given by r) with a length "size".
	 *  rep_avg_var stores the time(t), running density average and variance of X_curr, as well as N(t) and number of survivng replicates 
	 * for each time-step from previous surviving replicates ( <= r-1) respectively.
	 * In this function, X_curr is used to update the running avg and variance of rep_avg_var to account for current surviving replicate r(t).
	 * Note rep_avg_Var[t][3] stores number of surviving runs at t.
	*/
	vector<double> mean_prev(size);
	for(int s=0; s < Sp -2; s++)
	{
	for(int i=0; i<size; i++)
	{
		int r= t_avg_var_rep_N[i][4*s+ 5]+1; //Number of surviving replicates
		if(X_curr[i][2*s+ 0] == 0.0)
		{	//Non-surviving run. No-updates required. Also <rho(t)>x == 0, it remains 0 for t' > t.
			break;
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
	for(int i=0; i<size; i++)
	{	//For soil water (W) and surface water (0) respectively
		t_avg_var_rep_N[i][1] = (X_curr[i][2*(Sp-2)+1] + (j)*t_avg_var_rep_N[i][1])/(j+1);
		t_avg_var_rep_N[i][2] = (X_curr[i][2*(Sp-1)+1] + (j)*t_avg_var_rep_N[i][2])/(j+1); 
		//Note r>= 1 as all replicates of water are surviving as per definition.
	}

	//Swap memory to free up space.
	vector <double>().swap(mean_prev);
	

}
/**
void var_mean_incremental_surv_runs_MultiSp(double t_avg_var_rep_N[][Sp4_1], double X_curr[][Sp2], int size, int j)
{
	/** X_curr contains time-series data (indexed as (N(t), <rho(t)>x)) from the current surviving replicate (given by r) with a length "size".
	 *  rep_avg_var stores the time(t), running density average and variance of X_curr, as well as N(t) and number of survivng replicates 
	 * for each time-step from previous surviving replicates ( <= r-1) respectively.
	 * In this function, X_curr is used to update the running avg and variance of rep_avg_var to account for current surviving replicate r(t).
	 * Note rep_avg_Var[t][3] stores number of surviving runs at t.
	*/
/**
	double mean_prev[size];
	for(int s=0; s < Sp; s++)
	{
		int s_eff = s-2;
		if(s ==0)
		{  s_eff = 0; }
		else if(s_eff <= 0)
		{	continue; }
		for(int i=0; i<size; i++)
		{
			int r= t_avg_var_rep_N[i][4*s_eff+ 5]+1; //Number of surviving replicates
			//Number of surviving run at time t.
			if(X_curr[i][2*s+ 0] != 0.0)
			{
				// Only surviving runs are considered.
				mean_prev[i] = t_avg_var_rep_N[i][4*s_eff+ 3]; //Stores X_n-1 mean data as temp variable.
				if(t_avg_var_rep_N[i][4*s_eff + 5] == 0.0)
				{	//First surviving run encountered at time t, encountered here.
					t_avg_var_rep_N[i][4*s_eff+ 3] = X_curr[i][2*s+ 1]; t_avg_var_rep_N[i][4*s+ 4] = 0.0; 
					t_avg_var_rep_N[i][4*s_eff+ 6] = X_curr[i][2*s+ 0];
				}
				else
				{	//Not the first surviving run, which means r = rep_avg_var[i][3]+1 >=2
					t_avg_var_rep_N[i][4*s_eff+ 3] = (X_curr[i][2*s+1] + (r-1)*t_avg_var_rep_N[i][4*s_eff+ 3])/r; //Note r>= 1.
					t_avg_var_rep_N[i][4*s_eff+ 6] = (X_curr[i][2*s+0] + (r-1)*t_avg_var_rep_N[i][4*s_eff+6])/r; //Note r>= 1.
					//Calculating and storing incremental variance (V^2_n) [For formula derivation, refer to notes.]
					t_avg_var_rep_N[i][4*s_eff+4] = (t_avg_var_rep_N[i][4*s_eff+3] - mean_prev[i])*(t_avg_var_rep_N[i][4*s_eff+3] - mean_prev[i]) +
					(1/(r-1))*( (r-2)*t_avg_var_rep_N[i][4*s_eff+4] + 
					(X_curr[i][2*s+1] - t_avg_var_rep_N[i][4*s_eff+3])*(X_curr[i][2*s+1] - t_avg_var_rep_N[i][4*s_eff+3]));
				}
				t_avg_var_rep_N[i][4*s_eff+5] = r; //At the end, update number of succesful runs.
			}
		}
	}
	for(int s=1; s < 2; s++)
	{
		for(int i=0; i<size; i++)
		{
		//For soil water and surface water.
			t_avg_var_rep_N[i][s] = (X_curr[i][2*s+1] + (j-1)*t_avg_var_rep_N[i][s])/j; //Note r>= 1.
		}
	}

}
**/

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

// This function finds the maximum replicate number in the directory, given the filename pattern, so to avoid overwriting files.
int findMaxRepNo(string& parendir, const string& filenamePattern)
{
	 int maxRepNo = 0;
    //Filenamepattern is of the form: "/FRAME_RAND_ThreeSp_P_c_DP_G_128_T_57544_dt_0.12_a_0.049_D1_0.0298_D2_0.0522_dx_0.1_R_"
	// The pattern is used to match the filenames in the directory.

    std::regex pattern(filenamePattern + "(\\d+).csv");

    for (const auto& entry : std::filesystem::directory_iterator(parendir)) 
	{
        std::smatch match;
        string matchStr = "/" + entry.path().filename().string();
        //cout << "Checking: " << matchStr << " against " << filenamePattern << endl;
        // matchStr stores the filename
        //cout << regex_match (matchStr, match, pattern) << endl;
        if (regex_match (matchStr, match, pattern)) 
        {   
            //stringstream m0;
            int repNo = std::stoi(match[1]);
            //m0 << "For pattern: " << filenamePattern + "(\\d+).csv" << "  Matched: " << match[1] << endl; cout << m0.str();
            if (repNo > maxRepNo) 
            	maxRepNo = repNo;
        }
    }
    return maxRepNo;
}

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

void determine_neighbours_Sq4(int g, D3Vec_Int &neighbours_Sq4)
{
    
	//Stores neighbouring sites of species s and at site i in a square (Norm 1) radius of 1

    // Neighbours_Sq8 has dimensions S*(L*L)*2*2 where S is the number of species (=2) in this case
	// First row stores vertical neighbours, second row stores horizontal neighbours.
	int a=0; int b= 0; //Useful for accounting purposes
    
	for (int i=0; i<g*g; i++)
	{
		int p = int(i/g); int q = i%g;
		//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.

		int y_up = ((p-1)%g + g)%g; int y_down = ((p+1)%g + g)%g; // Eq to (p+ky)mod g. -1 <= ky <= 1
		int x_left = ((q-1)%g + g)%g; int x_right = ((q+1)%g + g)%g; // Eq to (q+kx)mod g. -1 <= kx <= 1
		neighbours_Sq4[i][0][0] = y_up*g +q; neighbours_Sq4[i][0][1] = y_down*g +q; 
		neighbours_Sq4[i][1][0] = p*g + x_left; neighbours_Sq4[i][1][1] = p*g + x_right;
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

//Sets the density dependent mortality rate.
void set_density_dependentmortality(double m)
{	Mm = m;	}

std::vector<std::pair<int, int>> computeNeighboringSitesCentral(int R)
{
    vector<pair<int, int>> centralNeighboringSites;
    // Compute neighboring sites for the central lattice site (0, 0) within a ball of radius R
    for (int x = -R; x <= R; ++x) 
    {
        for (int y = -R; y <= R; ++y) 
        {
            int distanceSquared = x * x + y * y;
            if (distanceSquared <= R * R) 
            {
                centralNeighboringSites.emplace_back(x, y);
            }
        }
    }
    // Sort the vector of neighboring sites based on squared distance
    // Sort the centralNeighboringSites based on squared distance
    sort(centralNeighboringSites.begin(), centralNeighboringSites.end(), comparePairs);

    return centralNeighboringSites;
}
std::vector<int> getNeighborsInBall(const vector<int>& sortedNeighboringSites, double f) 
{
	int numNeighbors = sortedNeighboringSites.size();
    int maxDistanceIndex = int(f*f*numNeighbors);
    return vector<int>(sortedNeighboringSites.begin(), sortedNeighboringSites.begin() + maxDistanceIndex);
}
void generateNeighboringSitesFromCentral(const vector<pair<int, int>>& centralNeighboringSites, 
                                         vector<int>& neighboringSites, int i, int j, int L)
{
    // Resize neighboringSites if necessary
    neighboringSites.resize(centralNeighboringSites.size());

    // Generate neighboring sites for the lattice site (i, j) based on the central neighboring sites for the lattice site (0, 0)
    for (size_t k = 0; k < centralNeighboringSites.size(); ++k) {
        const auto& offset = centralNeighboringSites[k];
        int dx = offset.first;
        int dy = offset.second;

        int nx = (i + dx + L) % L;  // Equivalent to (i + dx) mod L with reflective boundary conditions.
        int ny = (j + dy + L) % L;  // Equivalent to (j + dy) mod L with reflective boundary conditions.
        neighboringSites[k] = nx * L + ny;
    }
}

/**
std::vector<int> generateNeighboringSitesFromCentral(const vector<pair<int, int>>& centralNeighboringSites, int i, int j, int L)
{
	vector<int> neighboringSites;
    // Generate neighboring sites for the lattice site (i, j) based on the central neighboring sites for the lattice site (0, 0)
    for (const auto& offset : centralNeighboringSites) {
        int dx = offset.first;
        int dy = offset.second;

        int nx = (i + dx + L) % L;  // Equivalent to (i + dx) mod L with reflective boundary conditions.
        int ny = (j + dy + L) % L;  // Equivalent to (j + dy) mod L with reflective boundary conditions.
        neighboringSites.push_back(nx * L + ny);
    }

    return neighboringSites;
} 
*/


//----------------------------- Regular Rietkerk Integration Machinery --------------------------------------------//


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

void RK4_Wrapper_2D(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, 
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

	//t_meas[0] = 100; t_meas[1] = 200; t_meas[2] = 300; 

	//double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0};
	D2Vec_Double rho_rep_avg_var(tot_iter, vector<double> (Sp4_1, 0.0));

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

	double beta = c*gmax - d;
  	double p0istar, p0jstar, p0mstar; // Analytic steady state values.

	if(a < 1)
	{
		p0istar = 0.0; p0jstar = a/rW; p0mstar = a/(alpha*W0);
	}
	else
	{
		p0jstar = d*K[1]/beta; 
 		p0istar = (c/d)*(a - rW*p0jstar);
  		p0mstar = (a/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
	}
	
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.
		//double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2


		/**

		stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, alph, w0t, dix, dimitri;

  		L << g; tm << 0; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
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
			frame_dp << a << "," << i << ","<< Rh0[0][i] << "," << Rh0[1][i] << "," << Rh0[2][i] <<endl;
				
		}

		frame_dp.close();

		*/


		D2Vec_Double Rho_dt(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		//double Rho_M[tot_iter][Sp2] ={0.0}; //Stores <rho>x values per time step for a single replicate.
		D2Vec_Double Rho_M(tot_iter, vector<double> (Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		double perc = 0.015; double c_high[Sp] ={p0istar + dP, p0jstar, p0mstar}; double c_low[Sp] ={p0istar, p0jstar, p0mstar};
		
		init_randconstframe(Rho_dt, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
		

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1 << "Initial Conditions:\t Per:\t" << perc  << "\t C_High[0], C_High[1], C_High[2]:\t " << c_high[0] << "," << c_high[1] << "," << c_high[2]   
		<< "\t C_Lo[0], C_Lo[1], C_Lo[2]:\t " << c_low[0] << "," << c_low[1] << "," << c_low[2] << ",\t R: " << a << endl; //cout << m1.str();
		cout << m1.str();
		//std::ofstream errout; //Save Error Logs
		//std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
		errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();

		for(int i=0; i < g*g; i++)
		{
			for(int s= 0; s <Sp; s++)
			{
				Rho_tsar[s][i] = Rho_dt[s][i];
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

				if( t >= 100 && t <= 2500 || t== 0 ||  iter_index >= tot_iter -3 &&  iter_index <= tot_iter)
				{
					//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.
					
					stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, alph, w0t, dix, dimitri;

  					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  					rini << j; Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; gm << setprecision(3) << gmax;
					w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; dimitri << setprecision(3) << dP;
					a1 << a_st; a2  << a_end;
					// Three replicates are over.

					ofstream frame_dp;
  					// Creating a file instance called output to store output data as CSV
					frame_dp.open("../Data/Rietkerk/Frames/Non-Stochastic/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_Basic_P_c_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() 
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
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter, j); 
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
			output_dp.open("../Data/Rietkerk/Prelims/Non-Stochastic/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/PRELIM_AGGRAND_Basic_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
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







void first_order_critical_exp_delta(int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double d, double rW, double W0,  double D[], double K[], double dt, double dx, double dP, int r,  int g)
{

	
	//double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2
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

	//double perc = 0.015; double c_high[Sp] ={dP, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

  	  
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {

		//init_randconstframe(Rh0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho1(t)>>x,r			|    Var[<Rho1(t)>x],r    |
		**/

		RK4_Wrapper_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , a_start, a_end, dt, dx, dP, r, g);
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
	output_1stdp.open("../Data/Rietkerk/Non-Stochastic/1stOrderCC_Rietkerk_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + "_gmax_"+ gm.str() + "_alpha_"+ alph.str() + "_R_"+ rini.str() + ".csv");

	cout << "Save file name: " <<endl;
	cout << "../Data/Rietkerk/Non-Stochastic/1stOrderCC_Rietkerk_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() +
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


// ------------------------------------- Stochastic Integration Machinery ------------------------------------- //


// ALT INTEGRATION MACHINERY WHERE MIXING OF TIMESTEPS DUE TO LAPLACIAN MAY BE AN ISSUE.

void f_2Dor(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, 
	double alpha, double d, double rW, double W0, double D[], double K[], double t, double dt, double dx1_2, double g)
{
	//Vector function that updates an array containing (dP/dt, dW/dt, DO/dt) for each site in the lattice.
    //Based on the Rietkerk model for plant vegetation dynamics with the Dornic twist where linear and stoch term for vegetation are already taken care of.

	for(int i=0; i < g*g; i++)
	{
        //Equations for the density of plants, soil water and surface water at each site.
        //Note that the Laplacian is calculated using reflective boundary conditions.
        //The Laplacian is calculated using the 5-point stencil method.
        /**
		 * NOTE!!!!!:  Equivalent to dP/dt = c*g_max*P*W/(W+K1)
		 *  Linear and stochastic terms taken care of by Dornic integration routine previously.
		**/
		
        f[0][i] = c*gmax*Rho_M[1][i]*Rho_M[0][i]/(Rho_M[1][i] +K[1]);
        //Equivalent to dW/dt = alpha*(P+K2*W)/(P+K2) - rW*W + D*(Laplacian of W)
        f[1][i] = alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[2][i] -gmax*Rho_M[1][i]*Rho_M[0][i]/(Rho_M[1][i] +K[1]) - rW*Rho_M[1][i] 
        + (D[1]*dx1_2)*(Rho_M[1][nR2[i][0][0]]  + Rho_M[1][nR2[i][0][1]]  + Rho_M[1][nR2[i][1][0]]  + Rho_M[1][nR2[i][1][1]]  - 4*Rho_M[1][i]);
        //Equivalent to dO/dt = a - alpha*(P+K*W)/(P+K)*O + D*(Laplacian of O)
		f[2][i] = a - alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[2][i] + (D[2]*dx1_2)*(Rho_M[2][nR2[i][0][0]]  + Rho_M[2][nR2[i][0][1]] 
        + Rho_M[2][nR2[i][1][0]]  + Rho_M[2][nR2[i][1][1]]  - 4*Rho_M[2][i]);
		
	}

}

void RK4_Integrate_Stochastic(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D3Vec_Int &nR2,double a,double c,double gmax,double alpha,
		double d, double rW, double W0, double D[], double K[],double t,double dt,double dx, int g)
{
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	double dt6 = dt/6.0; double dt2 = dt/2.0; double dx1_2 = 1/(dx*dx);

	D2Vec_Double K1(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K2(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double K3(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K4(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double Rho_M(Sp, vector<double> (g*g, 0.0));

	f_2Dor(K1, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t, dt, dx1_2, g); //K1 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K1[s][i];
	}

	f_2Dor(K2, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt2, dt, dx1_2, g); //K2 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K2[s][i];
	}
	f_2Dor(K3, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt2, dt, dx1_2, g); //K3 updated.
	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt)*K3[s][i];
	}
	f_2Dor(K4, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt, dt, dx1_2, g); //K4 updated.
    
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

void rietkerk_Dornic_2D(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double d, 
	  double rW, double W0, double D[], double K[], double sigma[], double a_st, double a_end,  double dt, double dx, double dP, int r, int g)
{

	D3Vec_Int nR2 (g*g, D2Vec_Int (2, vector<int> (2, 0)));
    //nR2 is 3D [(L*L)x2x2] vector initialised to 0, which stores neighbours of all sites i (second dim) in a ball of NORM 1 (SQUARE) radius 1 (2ND & 3RD dim).
	// nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 

	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_Sq4(g, nR2); //Assigning neigbours.
	int tot_iter = 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	//t_meas[0] = 100; t_meas[1] = 200; t_meas[2] = 300; 

	//Defining Dornic variables for integration schema.
	//double bdt = b*dt; double cdt = c*dt; 
	double dx2 = dx*dx;
	double diff_coefficient[Sp]={0.0}; double beta[Sp]={0.0}; double lambda_exp[Sp]={0.0}; double lambda[Sp]={0.0};

	beta[0] = -d - 4*D[0]/(dx*dx);
	lambda_exp[0]= exp( (beta[0])*dt); lambda[0] = 2*(beta[0])/(sigma[0]*sigma[0]*(lambda_exp[0] -1.0));
	for(int s=0; s < Sp; s++)
	{
		diff_coefficient[s] = D[s]/(dx2);
	}
	
	//double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0}; 
	D2Vec_Double rho_rep_avg_var(tot_iter, vector<double> (Sp4_1, 0.0));
	//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
	rho_rep_avg_var[0][0] = 0;
	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[1+i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep

	double alpha_prime = diff_coefficient[0]*4; //Assume each Rho in neighbourhood for alpha estimation averages 1.
	double poiss_ru = lambda[0]*lambda_exp[0]*2.5; //Mean
	double mu_nought = 2.0*alpha_prime/(sigma[0]*sigma[0]) + poiss_ru;  // For Gamma, Beta =1. So Mean = 2*alpha/(sigma^2) + lambda = alpha/beta
	poiss_ru += 5*sqrt(poiss_ru); //Mean of Poisson Sampler + 5 SD.
	mu_nought += 6*sqrt(mu_nought); //Mean of  Mean of Gamma Sampler + 6 SD.

	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t a:\t" << a << "\t c:\t " << c << "\t dt,dx:\t " << dt << "," << dx << "\t gmax:\t " << gmax 
	<< "\n Lambda:\t" << lambda[0] << "\t Beta:\t " << beta[0] << "\t dt,dx:\t " << dt << "," << dx << "\t Poisson Cut-Off  For Veg:\t " << poiss_ru 
	<< "\t and  Gamma Cut-Off For Veg:\t " << mu_nought  << "\t and  Alpha  For Veg:\t " << alpha  <<endl; //cout << m1.str();
	cout << m0.str();
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	double kappa = c*gmax - d;
  	double p0istar, p0jstar, p0mstar; // Analytic steady state values.
	double astar = rW*d*K[1]/kappa; //Analytic critical point for a.

	if(a < astar)
	{
		p0istar = 0.0; p0jstar = a/rW; p0mstar = a/(alpha*W0);
	}
	else
	{
		p0jstar = d*K[1]/kappa; 
 		p0istar = (c/d)*(a - rW*p0jstar);
  		p0mstar = (a/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
	}

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
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number



		D2Vec_Double Rho_dt(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double DRho(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.

		D2Vec_Double Rho_M(tot_iter, vector<double> (Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		double perc = 0.015; double c_high[Sp] ={p0istar + dP, p0jstar, p0mstar}; double c_low[Sp] ={p0istar, p0jstar, p0mstar};
		
		init_randconstframe(Rho_dt, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
		
		for(int i=0; i < g*g; i++)
		{
			for(int s=0; s < Sp; s++)
			{
				DRho[s][i] = Rho_dt[s][i]; //Assigning initial conditions
				Rho_tsar[s][i] =Rho_dt[s][i];
			}
		}
		stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_1 << "Initial Conditions:\t Per:\t" << perc  << "\t C_High[0], C_High[1], C_High[2]:\t " << c_high[0] << "," << c_high[1] << "," << c_high[2]   
		<< "\t C_Lo[0], C_Lo[1], C_Lo[2]:\t " << c_low[0] << "," << c_low[1] << "," << c_low[2] << ",\t R: " << a << endl; //cout << m1.str();
		cout << m1.str();
		//std::ofstream errout; //Save Error Logs
		//std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
		errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();

		poisson_distribution<int> poisson; gamma_distribution <double> gamma; normal_distribution<double> norm(0.0, 1.0);

		double t=0; int index = 0;  //Initialise t
		int iter_index=0; int po=1; int so =1; int lo=1;

		while( t < t_max + dt)
		{
			// Basic Book-keeping below. Prints frames out to file at given time points. 
			//Also updates number and avg density of active sites at given time points.
			if(t == 0 || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
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
				iter_index+=1;
				if( t >= 2000 && t <= 5000 || t== 0 ||  index >= tot_iter -3 &&  index <= tot_iter-1)
				{
					//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.
					
					stringstream L, tm ,d3, p1, rini, cgm, a1, a2, Dm0, Dm1, alph, w0t, dix, dimitri, sig0, sig1;

  					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  					rini << j; Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; cgm <<  c*gmax;
					w0t << setprecision(3) << W0; alph <<  alpha; dimitri <<  dP;
					a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1];
					// Three replicates are over.

					ofstream frame_dp;
  					// Creating a file instance called output to store output data as CSV
					frame_dp.open("../Data/Rietkerk/Frames/Stochastic/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_Basic_P_c_DP_G_" 
					+ L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_alph_"+ alph.str() + "_cgmax_"+ cgm.str() 
					+ "_a_"+ p1.str() 
					+ "_dx_"+ dix.str() +  "_R_"+ rini.str() + ".csv"); //+ "_DV_"+ Dm0.str() + "_DW_"+ Dm1.str() 

					// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|
					frame_dp << "a_c,  x,  P(x; tmax), W(x; tmax), O(x; tmax) \n";
					for(int i=0; i< g*g; i++)
					{
						frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] <<endl;
				
					}

					frame_dp.close();
				}
				if( t > 0.0)
					index+=1;


			}

			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC. ONLY PLANTS EXPERIENCING DORNIC INTEGRATION.
			for(int s=0; s< 1; s++)
            {
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					double alpha_i = diff_coefficient[s]*(1*(DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
					DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]]));

					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					double mu = -1.0 + 2.0*alpha_i/(sigma[s]*sigma[s]);
					double ziggy = lambda[s]*lambda_exp[s]*Rho_dt[s][i]; 
					double gru; // Stores the Gamma random variable.
					if(ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if(ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng)*sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if(gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng)*sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i]= gauss/lambda[s];
						/**
						stringstream m6_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6_2 << "LARGE GRU: " << gru << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
						<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and Gamma (Gaussian) Approx Value:\t" << gauss << "\t and VEG RHO*:\t" << Rho_dt[s][i] << endl; 
						cout << m6_2.str(); errout.open(thr, std::ios_base::app); errout << m6_2.str(); errout.close(); */
					}
					else
					{
						gamma = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i]= gamma(rng)/lambda[s];
					}
						if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
						{
							stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
							m6 << "YOU HIT ROCK BOTTOM WITH: Rho*[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
							<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i << 
							"\t and GAMMA ENTRY:\t" << gru << endl; cout << m6.str();
							errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
						}
						Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
				}
			}

			// Finally RK4 integration of the remaining terms.

			RK4_Integrate_Stochastic(Rho_dt, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t, dt, dx, g);

			for(int s=0; s< Sp; s++)
			{
					
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					if( Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true and so==1)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}

					if( isnan(Rho_dt[s][i] == true))
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO NAAN PARATHA:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
						lo == -1;
					}
					DRho[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
					Rho_tsar[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
				}
			}

			t+=dt; //Update timestep

			if(lo == -1)
			{
				exit(3);
			}
		} // End of time loop.

		// Freeing up memory.
		vector<vector <double>>().swap(Rho_dt); vector<vector <double>>().swap(DRho); vector<vector <double>>().swap(Rho_tsar);

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
				for(int s =0; s< Sp-2; s++)
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
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter, j); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		// Freeing up memory.
		vector<vector <double>>().swap(Rho_M);

		if((j+1)%3 == 1)
		{	//Start logging at every multiple of 3
			stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, Dm, cgm,  sig0;

  			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dimitri << dP;
  			rini << j; Dm << setprecision(3) << D[0]; cgm << c*gmax; sig0 << sigma[0]; a1 << a_st; a2  << a_end;
			// Three replicates are over.

			ofstream output_dp;
  			// Creating a file instance called output to store output data as CSV.
			output_dp.open("../Data/Rietkerk/Prelims/Stochastic/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/PRELIM_AGGRAND_Basic_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
			"_cgmax_"+ cgm.str() + "_R_"+ rini.str() + ".csv"); //+ "_Sig0_" + sig0.str() //"_D0_"+ Dm.str() + 

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

	} // End of r loop.


	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	
	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|	r	|    L 		|    t 		|		<<W(x,t)>>x,r	|		<<O(x,t)>>x,r	|
		//     <<Rho(t)>>x,r		|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |

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



void first_order_critical_exp_delta_stochastic(int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double d, double rW, double W0,  double D[], double K[], double sigma[], double dt, double dx, double dP, int r,  int g)
{

	
	//double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2
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
	if(nProcessors > 32)
	{
		omp_set_num_threads(32); //Limiting use on Chunk. Don't be greedy.
	}

	//double perc = 0.015; double c_high[Sp] ={dP, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

  	  
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {

		//init_randconstframe(Rh0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho1(t)>>x,r			|    Var[<Rho1(t)>x],r    |
		**/
		rietkerk_Dornic_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , sigma, a_start, a_end, dt, dx, dP, r, g);
		//RK4_Wrapper_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , a_start, a_end, dt, dx, dP, r, g);
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
			 
          message3 << "Is this happening?\n" << "In thread: " << i << ", current size in MB: " << currSize/(1024.0*1024.0) 
		  << " and Peak Size (in MB): " << peaksSize/(1024.0*1024.0) << endl;
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, coco, tm ,d3, p1, p2, rini, Dm0, Dm1, cgm, alph, dix, dimitri, Sig0;

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

	

  	L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << a_start; p2  << a_end; rini << r; alph << alpha; cgm << c*gmax; coco << setprecision(4) << c;
  	Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; Sig0 << sigma[0]; dix << setprecision(2) << dx;
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
	output_1stdp.open("../Data/Rietkerk/Stochastic/1stOrderCC_Rietkerk_STOC_P_c_Delta_DP_G_" + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + "_cgmax_"+ cgm.str() + "_alpha_"+ alph.str() + "_R_"+ rini.str() + ".csv");
	 //+ "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() + "_S0_" + Sig0.str() +

	cout << "Save file name: " <<endl;
	cout << "../Data/Rietkerk/Stochastic/1stOrderCC_Rietkerk_STOC_P_c_Delta_DP_G_" + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() + "_S0_" + Sig0.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + "_cgmax_"+ cgm.str() + "_alpha_"+ alph.str() + rini.str() + ".csv";

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

	string header = " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs, # Active Sites";
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


//------------------- Vegetation + Grazer (+ Soil Water + Surface Water) -------------------//

void f_2Dor_2Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, 
	double alpha, double rW, double W0, double D[], double K[], double A[Sp][Sp], double H[Sp][Sp], double E[], double t, double dt, double dx1_2, double g)
{
	//cout << "Placeholder" << endl;
	//Vector function that updates an array containing ( dP/dt, dW/dt, DO/dt,  dG/dt) for each site in the lattice.
    //Based on the Rietkerk model for plant vegetation dynamics with the Dornic twist where linear and stoch term for vegetation (and grazer) are already taken care of.

	for(int i=0; i < g*g; i++)
	{
        //Equations for the density of plants, soil water, surface water and grazers at each site.
        //Note that the Laplacian is calculated using reflective boundary conditions.
        //The Laplacian is calculated using the 5-point stencil method.
        /**
		 * NOTE!!!!!:  Equivalent to dP/dt = c*g_max*P*W/(W+K1) - aij*hij*P*G/(1+aij*hij*V)
		 *  Linear and stochastic terms taken care of by Dornic integration routine previously.
		**/
        f[0][i] = c*gmax*Rho_M[2][i]*Rho_M[0][i]/(Rho_M[2][i] +K[1]) -(A[0][1]*Rho_M[0][i]/(1 + A[0][1]*H[0][1]*Rho_M[0][i]))*Rho_M[1][i];
		//Equivalent to dG/dt = e*(aij*V*P)/(1+aij*hij*V)  [NO PREDATION].
		f[1][i] = (E[1]*A[0][1]*Rho_M[0][i]/(1 + A[0][1]*H[0][1]*Rho_M[0][i]))*Rho_M[1][i];
        //Equivalent to dW/dt = alpha*(P+K2*W0)/(P+K2)*O - rW*W + D*(Laplacian of W)
        f[2][i] = alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[3][i] -gmax*Rho_M[2][i]*Rho_M[0][i]/(Rho_M[2][i] +K[1]) - rW*Rho_M[2][i] 
        + (D[2]*dx1_2)*(Rho_M[2][nR2[i][0][0]]  + Rho_M[2][nR2[i][0][1]]  + Rho_M[2][nR2[i][1][0]]  + Rho_M[2][nR2[i][1][1]]  - 4*Rho_M[2][i]);
        //Equivalent to dO/dt = a - alpha*(P+K*W)/(P+K)*O + D*(Laplacian of O)
		f[3][i] = a - alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[3][i] + (D[3]*dx1_2)*(Rho_M[3][nR2[i][0][0]]  + Rho_M[3][nR2[i][0][1]] 
        + Rho_M[3][nR2[i][1][0]]  + Rho_M[3][nR2[i][1][1]]  - 4*Rho_M[3][i]);
		
		
	}
}
void RK4_Integrate_Stochastic_2Sp(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D3Vec_Int &nR2,double a,double c,double gmax,double alpha,
		double rW, double W0, double D[], double K[], double A[Sp][Sp], double H[Sp][Sp], double E[], double t,double dt,double dx, int g)
{
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	double dt6 = dt/6.0; double dt2 = dt/2.0; double dx1_2 = 1/(dx*dx);

	D2Vec_Double K1(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K2(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double K3(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K4(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double Rho_M(Sp, vector<double> (g*g, 0.0));

	f_2Dor_2Sp(K1, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t, dt, dx1_2, g); //K1 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K1[s][i];
	}

	f_2Dor_2Sp(K2, Rho_M, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t + dt2, dt, dx1_2, g); //K2 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K2[s][i];
	}
	f_2Dor_2Sp(K3, Rho_M, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t + dt2, dt, dx1_2, g); //K3 updated.
	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt)*K3[s][i];
	}
	f_2Dor_2Sp(K4, Rho_M, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t + dt, dt, dx1_2, g); //K4 updated.
    
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
void rietkerk_Dornic_2D_2Sp(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0, 
		double D[], double v[], double K[], double sigma[], double a_st, double a_end, double A[Sp][Sp], double H[Sp][Sp], double E[], double M[], double pR[], double dt, double dx, double dP, int r, int g)
{
	//Store max value in pR[] in a variable r_max.
	double r_max_effective = (*max_element(pR, pR + Sp))/dx;
	std::vector<std::pair<int, int>> origin_Neighbourhood = computeNeighboringSitesCentral(r_max_effective); 

	D3Vec_Int nR2 (g*g, D2Vec_Int (2, vector<int> (2, 0)));
    //nR2 is 3D [(L*L)x2x2] vector initialised to 0, which stores neighbours of all sites i (second dim) in a ball of NORM 1 (SQUARE) radius 1 (2ND & 3RD dim).
	// nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 

	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_Sq4(g, nR2); //Assigning neigbours.
	//Append the value 0.0 to the start of t_meas.
	//t_meas.insert(t_meas.begin(), 0.0); //Inserting 0.0 at the start of t_meas.
	int tot_iter = t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";

	/**
	stringstream m0_1;	//To make cout thread-safe as well as non-garbled due to race conditions.
    m0_1 << "Values in t_measure (which is of total size " << t_meas.size() << ") are :" <<endl;
  	for (int i=0; i< t_meas.size(); i++)
  	{
		m0_1 << t_meas[i] << " ";
  	} 	m0_1 << endl;//cout << m1.str();
	cout << m0_1.str();
	errout.open(thr, std::ios_base::out); errout << m0_1.str(); errout.close();
	**/

	//Defining Dornic variables for integration schema.
	//double bdt = b*dt; double cdt = c*dt; 
	double dx2 = dx*dx; double dx1 = 1/dx;
	double diff_coefficient[Sp]={0.0}; double beta[Sp]={0.0}; double lambda_exp[Sp]={0.0}; double lambda[Sp]={0.0};

	beta[0] = -M[0] - 4*D[0]/(dx2);
	lambda_exp[0]= exp( (beta[0])*dt); lambda[0] = 2*(beta[0])/(sigma[0]*sigma[0]*(lambda_exp[0] -1.0));
	diff_coefficient[0] = D[0]/(dx2);
	for(int s=1; s < Sp-2; s++)
	{
		diff_coefficient[s] = D[s]/(dx2);
		//diff_coefficient[s] = (D[s] -(dx/2.0)*(v[s]*(1 - dt*vdx[s])))/(dx2);
		// Advection leads to excess diffusion, hence the correction term.
		if(diff_coefficient[s] <= 0)
		{
			stringstream m0;	//To make cout thread-safe as well as non-garbled due to race conditions.
    		m0 << "You Fucked Up With Species:\t" << s << "\t with D_eff[s]:\t " << diff_coefficient[s] 
			<< "\t dt,dx:\t " << dt << "," << dx << "\t v:\t " << v[s] <<endl; //cout << m1.str();
			cout << m0.str();
			errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();
		}
	}
	
	//double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0}; 
	D2Vec_Double rho_rep_avg_var(tot_iter, vector<double> (Sp4_1, 0.0));
	//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively. rho_rep_avg_var[0][0] = 0;
	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep

	double alpha_prime = diff_coefficient[0]*4; //Assume each Rho in neighbourhood for alpha estimation averages 1.
	double poiss_ru = lambda[0]*lambda_exp[0]*2.5; //Mean
	double mu_nought = 2.0*alpha_prime/(sigma[0]*sigma[0]) + poiss_ru;  // For Gamma, Beta =1. So Mean = 2*alpha/(sigma^2) + lambda = alpha/beta
	poiss_ru += 5*sqrt(poiss_ru); //Mean of Poisson Sampler + 5 SD.
	mu_nought += 6*sqrt(mu_nought); //Mean of  Mean of Gamma Sampler + 6 SD.

	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t a:\t" << a << "\t c:\t " << c << "\t dt,dx:\t " << dt << "," << dx << "\t gmax:\t " << gmax 
	<< "\n Lambda:\t" << lambda[0] << "\t Beta:\t " << beta[0] << "\t dt,dx:\t " << dt << "," << dx << "\t Poisson Cut-Off  For Veg:\t " << poiss_ru 
	<< "\t and  Gamma Cut-Off For Veg:\t " << mu_nought  << "\t and  Alpha  For Veg:\t " << alpha  <<endl; //cout << m1.str();
	cout << m0.str();
	errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	double kappa = c*gmax - M[0];
  	double p0istar, p0jstar, p0mstar, p0nstar; // Analytic steady state values.

	//Represent soil water, surface water, basal vegetation and grazer respectively.
	
	double astar = rW*M[0]*K[1]/kappa; //Analytic critical point for a.
	if(a < 1/24.0)
	{
		p0istar = 0.0; p0mstar = a/rW; p0nstar = a/(alpha*W0); p0jstar =  dP/50000.0; //= 0.4
		//p0istar = a/rW; p0jstar = a/(alpha*W0); p0mstar = 0; p0nstar = dP/50.0;
	}
	else
	{
		//p0jstar = M[0]*K[1]/kappa; p0istar = (c/M[0])*(a - rW*p0jstar); p0mstar = (a/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
		p0mstar = M[0]*K[1]/kappa; //rhoi*= 8000; pojstar = 160;
 		p0istar = (c/M[0])*(a - rW*p0mstar);  p0jstar = p0istar/50.0;
  		p0nstar = (a/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
	}
	
	//const vector<int> initcsv_columns = {2, 3, 4}; //Columns to be read from csv file.
	//D2Vec_Double const_ind_val = {{1, dP/50000.0 }}; 
	// Const value of dp/50000 for Rho indexed 1 (i.e. GRAZER). Other initial values determined by burn-in csv files through init_csvconstframe(...)

	auto start_t = high_resolution_clock::now();
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
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number



		D2Vec_Double Rho_dt(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double DRho(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		//double Rho_M[tot_iter][Sp2] ={0.0}; //Stores <rho>x values per time step for a single replicate.
		D2Vec_Double Rho_M(tot_iter, vector<double> (Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		double perc = 0.015; double c_high[Sp] ={p0istar +dP, p0jstar, p0mstar, p0nstar}; double c_low[Sp] ={p0istar, p0jstar, p0mstar, p0nstar};
		init_randconstframe(Rho_dt, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.


		/**

		stringstream  rain, rep_j, grid; 
		rain  << a; rep_j << int(j%2); grid << g;
    	string initcsv_filename= "Burning_In_Frames/FRAME_RAND_Basic_P_c_DP_G_" + grid.str() +"_T_91201_dt_0.2_alph_0.00833333_cgmax_0.0208333_a_" 
		+ rain.str() +"_dx_0.1_R_"  + rep_j.str() + ".csv";
		const vector<int> initcsv_columns = {2, 3, 4}; //Columns to be read from csv file.
		init_csvconstframe(Rho_dt, const_ind_val, initcsv_filename, initcsv_columns, g*g); //Initialise Rho_dt with burn-in frame.
		
		vector<double> temp_vec= {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species 's'
		double init_veg_num = occupied_sites_of_vector(temp_vec, g*g); //Finds number of occupied at given t.
		vector<double>().swap(temp_vec); //Flush temp_vec out of memory.

		

		if(init_veg_num == 0.0)
		{
			stringstream m1_1;
			m1_1 << "RUN-TIME WARNING: ZERO Active BURN-IN veg sites for Thread Rank:\t " << omp_get_thread_num() 
			<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << " CONSIDER RE-RUNNING OVER A DIFFERENT RANGE." <<endl; 
			cout << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
		}

		stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_1 << "Initial Conditions: BURN-IN # Active Veg sites:\t"<< init_veg_num << " for Thread Rank:\t " << omp_get_thread_num() 
		<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << endl; 
		cout << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();

		*/

		
		for(int i=0; i < g*g; i++)
		{
			for(int s=0; s < Sp; s++)
			{
				DRho[s][i] = Rho_dt[s][i]; //Assigning initial conditions
				Rho_tsar[s][i] =Rho_dt[s][i];
			}
		}
		
		stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_1 << "Initial Conditions:\t Per:\t" << perc  << "\t C_High[0], C_High[1], C_High[2], C_High[3] :\t " << c_high[0] << "," << c_high[1] << "," << c_high[2] << "," << c_high[3]   
		<< "\t C_Lo[0], C_Lo[1], C_Lo[2], C_Lo[3]:\t " << c_low[0] << "," << c_low[1] << "," << c_low[2] << "," << c_low[3] << ",\t R: " << a << endl; //cout << m1.str();
		cout << m1_1.str();
		//std::ofstream errout; //Save Error Logs
		//std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
		errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
		

		poisson_distribution<int> poisson; gamma_distribution <double> gamma; 
		uniform_real_distribution<double> unif; normal_distribution<double> norm(0.0, 1.0);
		

		double t=0; int index = 0;  //Initialise t
		int iter_index=0; int po=1; int so =1; int lo=1; int counter =0;

		while( t < t_max + dt)
		{

			// Basic Book-keeping below. Prints frames out to file at given time points. 
			//Also updates number and avg density of active sites at given time points.
			if(t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{

				double rhox_avg, rhox_num; //Stores spatial average and number of occupied sites at given t.
				for(int s=0; s<Sp; s++)
				{
					vector <double> temp= {Rho_dt[s].begin(),Rho_dt[s].end()}; //Rho_dt for species 's'
                    rhox_avg = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
				    rhox_num = occupied_sites_of_vector(temp, g*g); //Finds number of occupied at given t.
				    Rho_M[index][2*s] = rhox_num; Rho_M[index][2*s +1] = rhox_avg; 

					vector<double>().swap(temp); //Flush temp out of memory.0

                    if(rhox_num == 0.0 and s== 0)
					    break; //Finish negotiations as all basal species are dead.
				}

				vector <double> temp_alt= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
				rhox_num = occupied_sites_of_vector(temp_alt, g*g); //Finds number of occupied at given t.
				vector<double>().swap(temp_alt); //Flush temp out of memory.0

				/**
				stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m4 << "STATUS UPDATE AT TIME [t,i]\t" << t_meas[index] << " with # non-zero sites in VEG Rho_dt_0:  " << Rho_M[0][0] 
				<< " AND # Non-Zero in VEG DRho_0:  " << rhox_num <<   "  and # Non-Zero in Grazer Rho_dt_1:  " <<  Rho_M[0][2] << endl; cout << m4.str();
				errout.open(thr, std::ios_base::app); errout << m4.str(); errout.close(); */
				if( t >= 99 && t <= 10000 || t== 0 ||  index >= tot_iter -6 &&  index <= tot_iter-1)
				{
					//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.
					
					stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, alph, w0t, aij, hij, dix, dimitri, sig0, sig1;

  					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  					rini << j; Dm0 << D[0]*pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; gm << setprecision(3) << gmax;
					w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; dimitri  << dP;
					a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1];
					aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1];
					// Three replicates are over.

					ofstream frame_dp;

					if(t < 760)
					{
						//Only save one in three frames here.
						if(index%3 ==2 || t== 0)
						{
							frame_dp.open("../Data/Rietkerk/Frames/Stochastic/2Sp/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_TwoSp_P_c_DP_G_" 
							+ L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str()  
							+ "_D1_"+ Dm1.str() + "_dx_"+ dix.str() + "_R_"+ rini.str() + ".csv");

							// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|   W(x, tmax) 		|    O(x, tmax) 		|
							frame_dp << "a_c,  x,  P(x; t), G(x; t), W(x; t), O(x; t) \n";
							for(int i=0; i< g*g; i++)
							{	frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] << "," << Rho_dt[3][i] << endl;	}

							frame_dp.close();

							stringstream m3;
							m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
							cout << m3.str(); errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
						}
					}
					else if(t >= 760 && t < 6000)
					{
						//Only save one in two frames here.
						if(index%2 ==0)
						{
							frame_dp.open("../Data/Rietkerk/Frames/Stochastic/2Sp/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_TwoSp_P_c_DP_G_" 
							+ L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str()  
							+ "_D1_"+ Dm1.str() + "_dx_"+ dix.str() + "_R_"+ rini.str() + ".csv");

							// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|   W(x, tmax) 		|    O(x, tmax) 		|
							frame_dp << "a_c,  x,  P(x; t), G(x; t), W(x; t), O(x; t) \n";
							for(int i=0; i< g*g; i++)
							{	frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] << "," << Rho_dt[3][i] << endl;	}

							frame_dp.close();
							stringstream m3;
							m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
							cout << m3.str(); errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
						}
					}
					else
					{
						//Save all frames here.
						frame_dp.open("../Data/Rietkerk/Frames/Stochastic/2Sp/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_TwoSp_P_c_DP_G_"
						+ L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str()
						+ "_D1_"+ Dm1.str() + "_dx_"+ dix.str() + "_R_"+ rini.str() + ".csv");

						// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|   W(x, tmax) 		|    O(x, tmax) 		|
						frame_dp << "a_c,  x,  P(x; t), G(x; t), W(x; t), O(x; t) \n";
						for(int i=0; i< g*g; i++)
						{	frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] << "," << Rho_dt[3][i] << endl;	}

						frame_dp.close();
						stringstream m3;
						m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
					}
  					// Creating a file instance called output to store output data as CSV
					/**
					frame_dp.open("../Data/Rietkerk/Frames/Stochastic/2Sp/BURNIN_" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RANDBURNIN_TwoSp_P_c_DP_G_" 
					+ L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str()  
					+ "_D1_"+ Dm1.str() + "_dx_"+ dix.str() + "_R_"+ rini.str() + ".csv");

					// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|   W(x, tmax) 		|    O(x, tmax) 		|
					frame_dp << "a_c,  x,  P(x; t), G(x; t), W(x; t), O(x; t) \n";
					for(int i=0; i< g*g; i++)
					{	frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] << "," << Rho_dt[3][i] << endl;	}

					frame_dp.close();
					


					if(index%2 ==0)
					{
						stringstream m3;
						m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
						cout << m3.str(); errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
					} */
					
				}
				index+=1;
				/**
				stringstream m4_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m4_1 << "NEW INDEX VALUE\t" << index << " , WITH NEW T_MEAS:  " << t_meas[index] << endl; cout << m4_1.str();
				errout.open(thr, std::ios_base::app); errout << m4_1.str(); errout.close(); */
			}
			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC. ONLY LIVING MATTER EXPERIENCES DORNIC INTEGRATION.
			//First up vegetation.

			vector <double> temp_1= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
			double rhox_DR = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.
			vector<double>().swap(temp_1); //Flush temp out of memory.0

			temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
			double  rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
			vector<double>().swap(temp_1); //Flush temp out of memory.0
			/**
			stringstream m4_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
			m4_1 << "AFTER SAVING FRAME, STATUS UPDATE AT TIME [t,i]\t" << t_meas[index] << " with # non-zero sites in VEG Rho_dt_0:  " << rhox_DR 
			<< " AND # Non-Zero in VEG DRho_0:  " << rhox_dt <<  endl; cout << m4_1.str();
			errout.open(thr, std::ios_base::app); errout << m4_1.str(); errout.close(); **/
			//int counter_nzero_alpha=0;
			for(int s=0; s< 1; s++)
            {
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					double alpha_i = diff_coefficient[s]*((DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
					DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]]));

					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						/**
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "Much of VEG TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i << " \t with alpha_i, Rho(i), DRho(i):  "
						<< alpha_i << " , " << Rho_dt[s][i] << " , " << DRho[s][i] << "\t and Diff Const (D0/dx2): " << setprecision(16) << diff_coefficient[s] << endl; cout << m6.str();
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close(); */
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					//counter_nzero_alpha +=1;

					double mu = -1.0 + 2.0*alpha_i/(sigma[s]*sigma[s]);
					double ziggy = lambda[s]*lambda_exp[s]*Rho_dt[s][i]; 
					double gru; // Stores the Gamma random variable.
					if(ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if(ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng)*sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;

						/**
						stringstream m6_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6_2 << "LARGE ZIGGY: " << ziggy << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
						<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and Gaussian Approx Value:\t" << gauss << "\t and GRU:\t" << gru << endl; 
						cout << m6_2.str(); errout.open(thr, std::ios_base::app); errout << m6_2.str(); errout.close(); */
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if(gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng)*sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i]= gauss/lambda[s];
						/**
						stringstream m6_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6_2 << "LARGE GRU: " << gru << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
						<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and Gamma (Gaussian) Approx Value:\t" << gauss << "\t and VEG RHO*:\t" << Rho_dt[s][i] << endl; 
						cout << m6_2.str(); errout.open(thr, std::ios_base::app); errout << m6_2.str(); errout.close(); */
					}
					else
					{
						gamma = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i]= gamma(rng)/lambda[s];
					}
					/**
					stringstream m6_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
					m6_1 << "DORNIC VEG RHO* " << Rho_dt[s][i] << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
					<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i  << "\t and POISSON ENTRY:\t" << ziggy << "\t and GAMMA ENTRY:\t" << gru << endl; 
					cout << m6_1.str(); errout.open(thr, std::ios_base::app); errout << m6_1.str(); errout.close(); **/
					if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "VEG HIT ROCK BOTTOM WITH: Rho*"<< s <<"[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i << 
						"\t and GAMMA ENTRY:\t" << gru << endl; cout << m6.str();
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}
					Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
				}
			}
			// Now for higher order species.
			vector <double> temp_veg= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
			double rhox_num_veg = occupied_sites_of_vector(temp_veg, g*g); //Finds number of occupied sites at given t.
			double rhox_avg_veg = mean_of_vector(temp_veg, g*g); //Finds spatial average of densities at given t.
			vector<double>().swap(temp_veg); //Flush temp out of memory. 
			double nR_fac = 1 - rhox_num_veg/(g*g); //Factor to reduce the number of neighbours for gamma_i estimation
			if (nR_fac < 0.35)
			{	nR_fac = 0.35; }
			/**
			if(counter%50000 == 0)
			{

				temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
				rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, thr]\t" << t << " , " << omp_get_thread_num()  << "\t VEG DORNIC DONE with non-zero elements in DRho_0:  " 
				<< rhox_num_veg << "\t and non-zero elements in Rho_dt_0:  " << rhox_dt << "\t and INDEX VAL, T_MEAS:\t" << index <<  " , " << t_meas[index] << endl; cout << m5_1.str();
				errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			*/
			for(int s=1; s< Sp-2; s++)
            {
				for(int i=0;i<Rho_dt[0].size();i++)
				{

					//Firstly gamma needs to be estimated for each site.
					double gamma_i = 0.0;
					//To do this, first get neighbouring sites (i.e. perceptual range of species)
					int c_i = int(i/g); int c_j = i%g; //Current x and y coordinates of site i.
					std::vector<int> nR_Perp(origin_Neighbourhood.size(), 0); //Nearest Neighbours at site i
					generateNeighboringSitesFromCentral(origin_Neighbourhood, nR_Perp, c_i, c_j, g);
					// If nR_fac < 1, then reduce the number of neighbours to estimate gamma_i.
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i  
						<< "\t NEIGHBOURS FOUND WITH SIZE:" << nR_Perp.size() << endl; cout << m5_2.str();
						errout.open(thr, std::ios_base::app); errout << m5_2.str(); errout.close();
					}
					if(nR_fac < 1)
					{	nR_Perp.resize(int(nR_fac*nR_fac*nR_Perp.size())); } //= getNeighborsInBall(nR_Perp, nR_fac); }

					// Summing over all values in the perceptual range of species s (i.e. nR_Perp)
					
					for(int k=0; k< nR_Perp.size(); k++)
					{
						gamma_i += DRho[s-1][nR_Perp[k]];
					}
					//Finally normalising gamma_i by the number of sites in the perceptual range.
					//vector <double> temp= {Rho_dt[s].begin(),Rho_dt[s].end()}; //Rho_dt for species 's'
                    //double rhox_avg = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
					if(rhox_avg_veg == 0)
						gamma_i = 0;
					else
						gamma_i = gamma_i/(nR_Perp.size()*rhox_avg_veg); //Normalising gamma_i by the number of sites in the perceptual range.

					vector <int>().swap(nR_Perp); // vector <double>().swap(temp);  Remove dynamic vectors from memory and free up allocated space.

					//Finally check if gamma_i > 1, if so set it to 1.
					if(gamma_i > 1)
					{	gamma_i = 1; }
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i  
						<< "\t GAMMA FOUND WITH GAMMA_I:" << gamma_i << endl; cout << m5_2.str();
						errout.open(thr, std::ios_base::app); errout << m5_2.str(); errout.close();
					}

					//Next sampling the advection vector from v[s].
					unif = uniform_real_distribution<double>(0.0, 1.0);
					double ran = unif(rng);
					double vx = v[s]*cos(2*PI*ran); double vy = v[s]*sin(2*PI*ran); //Advection vector.
					double vx_abs = vx*sgn(vx); double vy_abs = vy*sgn(vy); //Absolute value of advection vector.

					// RECALL:	nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 

					double alpha_i = gamma_i*diff_coefficient[s]*(1*(DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
									DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]])) 
									+dx1*(1-gamma_i)*(vx_abs*DRho[s][nR2[i][1][sgn_index(vx)]]+ vy_abs*DRho[s][nR2[i][0][sgn_index(vy)]]);
					// NOTE: dx1 is the inverse of dx. sgn_index(vx) returns 0 if vx is positive, else 1.
					// Recall nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1;
					// ALSO NOTE: diff_coefficient[s] = D[s]/(dx*dx)
					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					//Beta depends on the gamma_i value at each lattice site for each species, hence needs to be updated at each iteration.

					beta[s] = -M[s] - 4*gamma_i*diff_coefficient[s] - (1-gamma_i)*(vx_abs + vy_abs)*dx1; // NOTE: diff_coefficient[s] = D[s]/(dx*dx)
					lambda_exp[s]= exp( (beta[s])*dt); lambda[s] = 2*(beta[s])/(sigma[s]*sigma[s]*(lambda_exp[s] -1.0));

					double mu = -1.0 + 2.0*alpha_i/(sigma[s]*sigma[s]);
					double ziggy = lambda[s]*lambda_exp[s]*Rho_dt[s][i]; 
					double gru; // Stores the Gamma random variable.
					if(ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if(ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng)*sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if(gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng)*sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i]= gauss/lambda[s];
					}
					else
					{
						gamma = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i]= gamma(rng)/lambda[s];
					}
					if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "YOU HIT ROCK BOTTOM WITH: Rho*"<< s <<"[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i << 
						"\t and GAMMA ENTRY:\t" << gru << endl; cout << m6.str();
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
						exit(3);
					}
					Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
					/**
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i  
						<< "\t DORNIC VALUE OF GRAZER UPDATED WITH RHO*(I,t): " << Rho_tsar[s][i] << endl; cout << m5_2.str();
						errout.open(thr, std::ios_base::app); errout << m5_2.str(); errout.close();
					}
					*/
				}
			}
			/**
			if(counter%50000 == 0)
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, thr]\t" << t << " , " << omp_get_thread_num()  << "\t GRAZER DORNIC DONE." << endl; cout << m5_1.str();
				errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			*/

			// Finally RK4 integration of the remaining terms.

			RK4_Integrate_Stochastic_2Sp(Rho_dt, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t, dt, dx, g);
			/**
			if(counter%50000 == 0)
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, thr]\t" << t << " , " << omp_get_thread_num()  << "\t RK4 INTEGRATION OF REMAINDER TERMS DONE." << endl; cout << m5_1.str();
				errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			*/


			for(int s=0; s< Sp; s++)
			{
					
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					if( Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true and so==1)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO "<<s<<"  FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}

					if( isnan(Rho_dt[s][i] == true))
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO "<<s<<"  NAAN PARATHA:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
						lo == -1;
					}
					DRho[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
					Rho_tsar[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
				}
			}

			t+=dt; //Update timestep
			counter+=1; //Update counter

			if(counter%90000 ==1)
			{
				//Find memory usage and time taken for every 50000 iterations.
				//size_t currentSize = getCurrentRSS( ); //Check Memory usage.
				//size_t peakSize    = getPeakRSS( );
				auto end_t = high_resolution_clock::now();
				auto elapsed_min = duration_cast<minutes>(end_t-start_t);
				// Reset the clock
				start_t = high_resolution_clock::now();
				stringstream m7;
				m7 << "STATUS UPDATE AT [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j 
				//<< " ,with current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) 
				<< " and time taken for " << counter << " iterations: " << elapsed_min.count() << " min." <<  endl; cout << m7.str();
				errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();

				//Report syst every 50000 iterations.
			

			}

			if(lo == -1)
			{
				exit(3);
			}
		} // End of time loop.
		
		vector<vector <double>>().swap(Rho_dt); vector<vector <double>>().swap(DRho); vector<vector <double>>().swap(Rho_tsar);
		// Freeing up memory.
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
				for(int s =0; s< Sp-2; s++)
				{
					rho_rep_avg_var[i][4*s+ 3] = Rho_M[i][2*s + 1]; // Avg density of frame at given time point i.
					rho_rep_avg_var[i][4*s+ 4] = 0.0; // Variance in avg frame density.
					rho_rep_avg_var[i][4*s+ 6] = Rho_M[i][2*s]; // Number of active sites in frame at given time point i.

					if(Rho_M[i][2*s] > 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s+5 ] = 1;  } //One surviving run as of now
					else if(Rho_M[i][2*s] == 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s + 5] = 0;  } //No surviving run as of now
				}
				rho_rep_avg_var[i][1] = Rho_M[i][2*(Sp-2) + 1]; // Avg density of frame at given time point i.
				rho_rep_avg_var[i][2] = Rho_M[i][2*(Sp-1) + 1]; // Avg density of frame at given time point i.
				
			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			// Second or higher replicate, use incremental advances.
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter, j); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		// Freeing up memory.
		vector<vector <double>>().swap(Rho_M); 

		if((j+1)%3 == 1)
		{	//Start logging at every multiple of 3
			stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, Dm, cgm,  sig0;

  			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dimitri << dP;
  			rini << j; Dm << D[0]; cgm << c*gmax; sig0 << sigma[0];
			// Three replicates are over.

			ofstream output_dp;
  			// Creating a file instance called output to store output data as CSV.
			output_dp.open("../Data/Rietkerk/Prelims/Stochastic/2Sp/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/PRELIM_AGGRAND_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
			"_D0_"+ Dm.str() + "_cgmax_"+ cgm.str() + "_R_"+ rini.str() + ".csv");

			// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
			output_dp << " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs Rho0, # Active Sites Rho0, <<Rho1(x; t)>_x>_r, Var[<Rho1(t)>_x]_r, # Surviving Runs Rho1, # Active Sites Rho1, \n";
			for(int i=0; i< tot_iter; i++)
			{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |

				output_dp << a << ","<< j << "," << g << "," << rho_rep_avg_var[i][0] << "," << rho_rep_avg_var[i][1] << "," << rho_rep_avg_var[i][2] << ",";
				for(int s=0; s <Sp-2; s++)
				{
					output_dp << rho_rep_avg_var[i][4*s + 3] << "," <<
					rho_rep_avg_var[i][4*s +4] << "," << rho_rep_avg_var[i][4*s +5] << "," << rho_rep_avg_var[i][4*s +6] << ","; 
				}
				output_dp <<  endl;
				
			}

			output_dp.close();
		}
		

	} // End of r loop.

	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|    r 		|     L 		|    t 		|     <<W(t)>>x,r			|      <<O(t)>>x,r		|
		//..     <<Rho0(t)>>x,r			|    Var0[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |....

		Rho.push_back({a, static_cast<double>(r), static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], 
					rho_rep_avg_var[i][2]});
		for(int s=0;s < Sp-2; s++)
		{
			Rho[i].insert(Rho[i].end(), { rho_rep_avg_var[i][4*s +3], rho_rep_avg_var[i][4*s +4], 
			rho_rep_avg_var[i][4*s +5], rho_rep_avg_var[i][4*s +6]});		
		}
	}
	stringstream m10;
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m10.str();
	errout.open(thr, std::ios_base::app); errout << m10.str(); errout.close();
}


void first_order_critical_exp_delta_stochastic_2Sp(int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double rW, double W0,  double D[], double v[], double K[], double sigma[], double A[Sp][Sp], double H[Sp][Sp], double E[], double M[], double pR[],
	double dt, double dx, double dP, int r,  int g)
{

	//double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2
	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//double mean[Sp] = {1.0, 0.05}; double sd[Sp] = {0.25, 0.0125};
	//init_randframe(Rho_0, Sp,  g*g); //Returns Rho_0 with a full initial frame filled with 0.2.
	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "NOSTRA" <<endl;
  	// The pspace to iterate over.
    vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

	t_measure.insert(t_measure.begin(), 0.0);

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
	if(nProcessors > 32)
	{
		omp_set_num_threads(32); //Limiting use on Chunk. Don't be greedy.
	}

	//double perc = 0.015; double c_high[Sp] ={dP, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

  	  
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {

		//init_randconstframe(Rh0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho1(t)>>x,r			|    Var[<Rho1(t)>x],r    |
		**/
		rietkerk_Dornic_2D_2Sp(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, rW, W0, D, v, K , sigma, a_start, a_end, A,H,E,M, pR, dt, dx, dP, r, g);
		//RK4_Wrapper_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , a_start, a_end, dt, dx, dP, r, g);
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
			 
          message3 << "Is this happening?\n" << "In thread: " << i << ", current size in MB: " << currSize/(1024.0*1024.0) 
		  << " and Peak Size (in MB): " << peaksSize/(1024.0*1024.0) << endl;
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, coco, tm ,d3, p1, p2, rini, Dm0, Dm1, aij, hij, gm, alph, dix, dimitri, Sig0;

	

  	L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << a_start; p2  << a_end; rini << r; alph << alpha; gm << gmax; coco << setprecision(4) << c;
  	Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; Sig0 << sigma[0]; dix << setprecision(2) << dx; aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1];

	ofstream output_1stdp;
  // Creating a file instance called output to store output data as CSV.
	output_1stdp.open("../Data/Rietkerk/Stochastic/2Sp/1stOrderCC_Rietkerk_STOC_P_c_Delta_DP_G_" + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D1_"+ Dm1.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + 
	"_alpha_"+ alph.str() + "_R_"+ rini.str() + ".csv");

	cout << "Save file name: " <<endl;
	cout << "../Data/Rietkerk/Stochastic/2Sp/1stOrderCC_Rietkerk_STOC_P_c_Delta_DP_G_" + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D1_"+ Dm1.str() + "_S0_" + Sig0.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_aij_"+ aij.str() + "_hij_"+ hij.str() + "_dx_"+ dix.str() + 
	"_gmax_"+ gm.str() + "_alpha_"+ alph.str() + "_R_"+ rini.str() + ".csv" <<endl;

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

	string header = " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs Rho0, # Active Sites Rho0, <<Rho1(x; t)>_x>_r, Var[<Rho1(t)>_x]_r, # Surviving Runs Rho1, # Active Sites Rho1";
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



//------------------- Vegetation + Grazer + Predator (+ Soil Water + Surface Water) -------------------//

void f_2Dor_3Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, 
	double alpha, double rW, double W0, double (&D)[Sp], double (&K)[3], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t, double dt, double dx1_2, double g)
{
	//Vector function that updates an array containing ( dP/dt, dW/dt, DO/dt,  dG/dt, dPr/dt) for each site in the lattice.
    //Based on the Rietkerk model for plant vegetation dynamics with the Dornic twist where linear and stoch term for vegetation (and grazer) are already taken care of.

	for(int i=0; i < g*g; i++)
	{
        //Equations for the density of plants, soil water, surface water and grazers at each site.
        //Note that the Laplacian is calculated using reflective boundary conditions.
        //The Laplacian is calculated using the 5-point stencil method.
        /**
		 * NOTE!!!!!:  Equivalent to dP/dt = c*g_max*P*W/(W+K1) - aij*hij*P*G/(1+aij*hij*V)
		 *  Linear and stochastic terms taken care of by Dornic integration routine previously.
		**/
        f[0][i] = c*gmax*Rho_M[3][i]*Rho_M[0][i]/(Rho_M[3][i] +K[1]) -(A[0][1]*Rho_M[0][i]/(1 + A[0][1]*H[0][1]*Rho_M[0][i]))*Rho_M[1][i];
		//Equivalent to dG/dt = ej*(aij*V*P)/(1+aij*hij*V)  -ajm*G*P/(1+ajm*hjm*G).
		f[1][i] = (E[1]*A[0][1]*Rho_M[0][i]/(1 + A[0][1]*H[0][1]*Rho_M[0][i]))*Rho_M[1][i]  
					-A[1][2]*Rho_M[1][i]/(1 + A[1][2]*H[1][2]*Rho_M[1][i])*Rho_M[2][i];
        //Equivalent to dPr/dt = em*(ajm*G*Pr)/(1+ajm*hjm*G)
		f[2][i] = (E[2]*A[1][2]*Rho_M[1][i]/(1 + A[1][2]*H[1][2]*Rho_M[1][i]))*Rho_M[2][i]; //-Mm*Rho_M[2][i]*Rho_M[2][i];
		//Equivalent to dW/dt = alpha*(P+K2*W0)/(P+K2)*O - g_max*P*W/(W+K1) - rW*W + D*(Laplacian of W)
        f[3][i] = alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[4][i] -gmax*Rho_M[3][i]*Rho_M[0][i]/(Rho_M[3][i] +K[1]) - rW*Rho_M[3][i] 
        + (D[3]*dx1_2)*(Rho_M[3][nR2[i][0][0]]  + Rho_M[3][nR2[i][0][1]]  + Rho_M[3][nR2[i][1][0]]  + Rho_M[3][nR2[i][1][1]]  - 4*Rho_M[3][i]);
        //Equivalent to dO/dt = a - alpha*(P+K2*W)/(P+K2)*O + D*(Laplacian of O)
		f[4][i] = a - alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[4][i] + (D[4]*dx1_2)*(Rho_M[4][nR2[i][0][0]]  + Rho_M[4][nR2[i][0][1]] 
        + Rho_M[4][nR2[i][1][0]]  + Rho_M[4][nR2[i][1][1]]  - 4*Rho_M[4][i]);
		
	}
}
void RK4_Integrate_Stochastic_MultiSp(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D3Vec_Int &nR2,double a,double c,double gmax,double alpha,
		double rW, double W0, double (&D)[Sp], double (&K)[3], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t,double dt,double dx, int g)
{
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	double dt6 = dt/6.0; double dt2 = dt/2.0; double dx1_2 = 1/(dx*dx);

	D2Vec_Double K1(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K2(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double K3(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K4(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double Rho_M(Sp, vector<double> (g*g, 0.0));

	f_2Dor_3Sp(K1, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t, dt, dx1_2, g); //K1 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K1[s][i];
	}

	f_2Dor_3Sp(K2, Rho_M, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t + dt2, dt, dx1_2, g); //K2 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K2[s][i];
	}
	f_2Dor_3Sp(K3, Rho_M, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t + dt2, dt, dx1_2, g); //K3 updated.
	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt)*K3[s][i];
	}
	f_2Dor_3Sp(K4, Rho_M, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t + dt, dt, dx1_2, g); //K4 updated.
    
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
				<< "\t, K3[s][i]:\t" << K3[s][i] << "\t, K4[s][i]:\t" << K4[s][i] << "\t AND Rho(t-dt)" << Rho_tsar[s][i] << endl; //cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			
		}
	}

}

void save_frame(D2Vec_Double &Rho_t, double a, double a_st, double a_end, double c,double gmax,double alpha, double rW, double W0, double (&D)[Sp], double (&K)[3], 
	double sigma[], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t, double dt, double dx, double dP, int r, int g, double Gstar /*=-1*/)
{
	stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, Dm2, alph, w0t, aij, hij, dix, dimitri, sig0, sig1, geq;
	//Save the frame to a file.
	L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
	rini << r; Dm0 << D[0]*pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2]; gm << setprecision(3) << gmax;
	w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; dimitri  << dP;
	a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1]; geq << setprecision(5) << Gstar;
	aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1];

	string filenamePattern = frame_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
    + "_a_" + p1.str()  + "_D1_"+ Dm1.str() + "_D2_"+ Dm2.str() + "_dx_"+ dix.str() + "_R_";
	//Next, designate the naive filename
	string filename = filenamePattern + rini.str() + ".csv";
	string parendir = "";

	
	// Creating a file instance called output to store output data as CSV.
	if(Gstar != -1)
		parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str();
	else
		parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();
	

	filename = parendir + filename;
	// If the file already exists, increment the replicate number by 1 of the maximum replicate number corresponding to the filename pattern.
	// This is to avoid overwriting files.
	if (fs::exists((filename))) 
	{
		int maxRepNo = findMaxRepNo(parendir, filenamePattern);
		//maxRepNo++; //Increment the maximum replicate number by 1.
		filename = parendir + filenamePattern + std::to_string(maxRepNo + 1) + ".csv";
		stringstream m1;

		m1 << "File exists. New filename: " << filename << " \n and maxRepNo: " << maxRepNo << endl;
		cout << m1.str();
	}

	std::ostringstream oss;
	oss << frame_header;
	for(int i=0; i< g*g; i++)
	{	
		oss << a << "," << i;
		for(int s=0; s < Sp; s++)
			oss << "," << Rho_t[s][i];
		oss << "\n";
	}

	FILE *fp = fopen(filename.c_str(), "w");
	if(fp)
	{
		fprintf(fp, "%s", oss.str().c_str());
		fclose(fp);
	}
	else
	{
		std::cerr << "Error: Could not open file " << filename << std::endl;
	}
	
}

//Calculates gamma for 3 Sp Rietkerk model (details in PDF) [ 0 < ]
void calc_gamma_3Sp(const vector<pair<int, int>>& centralNeighboringSites, D2Vec_Double &Rho_t, vector <double> &gamma_i,  
			double (&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, double nVeg_frac, int r_max, int i, int L)
{
	mutex errorMutex; // Mutex to make error messages thread-safe
	double eps = 1.0e-50; //Small number to avoid division by zero.
	std::vector<int> nR_Perp;
	int c_i = int(i/L); int c_j = i%L; //Current x and y coordinates of site i.
	generateNeighboringSitesFromCentral(centralNeighboringSites, nR_Perp, c_i, c_j, L);

	if(nVeg_frac < 0.35)
		nVeg_frac = 0.35; //Minimum fraction of vegetation cover.

	if(nVeg_frac < 1)
		nR_Perp.resize(int(nVeg_frac*nVeg_frac*nR_Perp.size()));
	// In this case, r_effective is reduced to nVeg_frac*r_max, hence the reduction in the number of nearest neighbours.
	// rfrac is sorted in descending order.
	double fr_prev = 1.0; //Fraction of perception radius to max radius.
	for(const auto& frac: rfrac)
	{
		int s = frac.second;
		double fr = frac.first; //Fraction of perception radius to max radius.

		if( fr == 0.0)
		{	gamma_i[s] = 0.0; continue;	}
		else if(fr < fr_prev)
		{
			nR_Perp.resize(int(fr*fr*nR_Perp.size()/(fr_prev*fr_prev)));
			fr_prev = fr;
		}
		
		// Gamma for grazers (indexed 1)
		if (s == 1)
		{	
			
			if(Rho_avg[s] < eps)
			{	gamma_i[s] = 1.0; continue;	} // No grazers left
			else if( Rho_avg[2] < eps)
			{	
				for(int k=0; k< nR_Perp.size(); k++)
				{		gamma_i[s] += Rho_t[s-1][nR_Perp[k]];	}
			} // No predators left
			else
			{
				if(Rho_avg[0] == 0)
				{	gamma_i[s] = 0.0;	} // No vegetation
				else
				{	
					double rho_pred1 = 1/Rho_avg[2];
					for(int k=0; k< nR_Perp.size(); k++)
					{		
						gamma_i[s] += Rho_t[s-1][nR_Perp[k]]*(1 -Rho_t[2][nR_Perp[k]]*rho_pred1);	
					}
					
					gamma_i[s] = gamma_i[s]/(nR_Perp.size()*Rho_avg[s-1]);	
				}
			}

		}
		if(s == 2)
		{
			if(Rho_avg[s] < eps)
			{	gamma_i[s] = 1.0; continue;	} // No predators left
			else if( Rho_avg[0] < eps)
			{	
				for(int k=0; k< nR_Perp.size(); k++)
				{		gamma_i[s] += Rho_t[s-1][nR_Perp[k]];	}
			} // No vegetation left
			else
			{
				if(Rho_avg[1] < eps)
				{	gamma_i[s] = 0.0;	} // No grazers left
				else
				{	
					double rho_veg1 = 1/Rho_avg[0];
					for(int k=0; k< nR_Perp.size(); k++)
					{		
						gamma_i[s] += Rho_t[s-1][nR_Perp[k]]*(1 -Rho_t[0][nR_Perp[k]]*rho_veg1);	
					}
					gamma_i[s] = gamma_i[s]/(nR_Perp.size()*Rho_avg[s-1]);	
				}
			}
		}
		if(gamma_i[s] > 1)
		{	gamma_i[s] = 1.0;	} //Ensuring gamma is not greater than 1.
		else if(gamma_i[s] < 0)
		{	gamma_i[s] = 0.0;	} //Ensuring gamma is not less than 0.
		else if(isnan(gamma_i[s]) == true || isinf(gamma_i[s]) == true)
		{
			errorMutex.lock();
			std::cerr << "Gamma is NaN with value: " << gamma_i[s] << " For [s, thr, i]\t" << s << " , " << omp_get_thread_num() 
					  << " , " << i << " ] with Rho_avg[0]: " << Rho_avg[0] << " and Rho_avg[1]: " << Rho_avg[1] << " and Rho_avg[2]: " << Rho_avg[2] 
					  << " \n and Rho_t[0][i]: " << Rho_t[0][i] << " and Rho_t[1][i]: " << Rho_t[1][i] << " and Rho_t[2][i]: " << Rho_t[2][i] 
					  << " \n and Nveg_frac: " << nVeg_frac << " and r_max: " << r_max << " and rfrac.size(): " << rfrac.size() << " and nR_Perp.size(): " << nR_Perp.size()
					  << std::endl;
			errorMutex.unlock();
		}
	}
	vector <int>().swap(nR_Perp);
}

void rietkerk_Dornic_2D_MultiSp(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0, 
	double (&D)[Sp], double v[], double (&K)[3], double sigma[], double a_st, double a_end, double a_c, double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double pR[], 
	double chigh[], double clow[], double dt, double dx, double dP, int r, int g)
{
	double perc =0.015; //Percentage of high density patches.
	//Store max value in pR[] in a variable r_max.
	//NOTE: SpB = Sp - 2; //Number of species excluding water terms.
	double r_max_effective = (*max_element(pR, pR + SpB))/dx;
	// Store pR[]/dx as a fraction of r_max_effective in r_frac[].
	//Initialise vector<std::pair<int, int>> r_frac with initial size SpB and set all values to 0.0.
	vector <std::pair<double, int>> r_frac(SpB, {0.0, 0});
 
	for(int s=0; s< Sp-2; s++)
	{	r_frac[s] = {pR[s]/(dx*r_max_effective), s};	}

	// Sort the r_frac[] array in descending order and store the indices in r_sort[].
	sort(r_frac.begin(), r_frac.end(), [](std::pair<double, int> &a, std::pair<double, int> &b) 
				{ return a.first > b.first; });
	std::vector<std::pair<int, int>> origin_Neighbourhood = computeNeighboringSitesCentral(r_max_effective);

	if(omp_get_thread_num()== 1)
	{
		for(int i=0; i< r_frac.size(); i++)
		{
			stringstream m0;	//To make cout thread-safe as well as non-garbled
			m0 << "r_frac[" << i << "]:\t" << r_frac[i].first << "\t and r_sort[" << i << "]:\t" << r_frac[i].second << endl;
			cout << m0.str();
		}
	}

	//Initialise variables for the RK4 integration.
	D3Vec_Int nR2 (g*g, D2Vec_Int (2, vector<int> (2, 0)));
	//nR2 is 3D [(L*L)x2x2] vector initialised to 0, which stores neighbours of all sites i (second dim) in a ball of NORM 1 (SQUARE) radius 1 (2ND & 3RD dim).
	// nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 

	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_Sq4(g, nR2); //Assigning neigbours.
	//Append the value 0.0 to the start of t_meas.
	//t_meas.insert(t_meas.begin(), 0.0); //Inserting 0.0 at the start of t_meas.
	int tot_iter = t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	//double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0}; 
	D2Vec_Double rho_rep_avg_var(tot_iter, vector<double> (Sp4_1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively. rho_rep_avg_var[0][0] = 0;
	//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively. rho_rep_avg_var[0][0] = 0;
	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep


	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";

	//Defining Dornic variables for integration schema.
	//double bdt = b*dt; double cdt = c*dt; 
	double dx2 = dx*dx; double dx1 = 1/dx; double dx1_2 = 1/(dx2); double dxb2 = dx/2.0;
	double diff_coefficient[Sp]={0.0}; double beta[Sp]={0.0}; double lambda_exp[Sp]={0.0}; double lambda[Sp]={0.0}; double sigma2_1[Sp]={0.0};
	vector<std::pair<double, double>> diff_eff(SpB, {0.0, 0.0}); //Stores effective diffusion coefficient for each species after correcting for advection.
	vector <double> gamma_i(SpB, 0.0); //Stores gamma for each species.

	beta[0] = -M[0] - 4*D[0]/(dx2);
	lambda_exp[0]= exp( (beta[0])*dt); lambda[0] = 2*(beta[0])/(sigma[0]*sigma[0]*(lambda_exp[0] -1.0));
	diff_coefficient[0] = D[0]/(dx2); sigma2_1[0] = 1/(sigma[0]*sigma[0]);
	for(int s=1; s < Sp-2; s++)
	{
		diff_coefficient[s] = D[s]/(dx2);
		sigma2_1[s] = 1/(sigma[s]*sigma[s]);
		//diff_coefficient[s] = (D[s] -(dx/2.0)*(v[s]*(1 - dt*vdx[s])))/(dx2);
		// Advection leads to excess diffusion, hence the correction term.
		if(diff_coefficient[s] <= 0)
		{
			stringstream m0;	//To make cout thread-safe as well as non-garbled due to race conditions.
    		m0 << "You Fucked Up With Species:\t" << s << "\t with D_eff[s]:\t " << diff_coefficient[s] 
			<< "\t dt,dx:\t " << dt << "," << dx << "\t v:\t " << v[s] <<endl; //cout << m1.str();
			cout << m0.str();
			errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();
		}
	}
	double Gstar = M[2]/((E[2] -M[2]*H[1][2])*A[1][2]); // MFT Coexistance Steady state for grazer.
	double alpha_prime = diff_coefficient[0]*4; //Assume each Rho in neighbourhood for alpha estimation averages 1.
	double poiss_ru = lambda[0]*lambda_exp[0]*2.5; //Mean
	double mu_nought = 2.0*alpha_prime/(sigma[0]*sigma[0]) + poiss_ru;  // For Gamma, Beta =1. So Mean = 2*alpha/(sigma^2) + lambda = alpha/beta
	poiss_ru += 5*sqrt(poiss_ru); //Mean of Poisson Sampler + 5 SD.
	mu_nought += 6*sqrt(mu_nought); //Mean of  Mean of Gamma Sampler + 6 SD.

	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t a:\t" << a << "\t c:\t " << c << "\t dt,dx:\t " << dt << "," << dx << "\t gmax:\t " << gmax 
	<< "\n Lambda:\t" << lambda[0] << "\t Beta:\t " << beta[0] << "\t dt,dx:\t " << dt << "," << dx << "\t Poisson Cut-Off  For Veg:\t " << poiss_ru 
	<< "\t and  Gamma Cut-Off For Veg:\t " << mu_nought  << "\t and  Alpha  For Veg:\t " << alpha  
	<< "\n with MFT biomass density of Grazer at Coexistance = " << Gstar << " kg/km^2\n" << endl; //cout << m1.str();
	omp_get_thread_num() == 1 ? cout << m0.str() : cout << "";
	errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	//const vector<int> initcsv_columns = {2, 3, 4}; //Columns to be read from csv file.
	
	// Const value of dp/50000 for Rho indexed 1 (i.e. GRAZER). Other initial values determined by burn-in csv files through init_csvconstframe(...)

	auto start_t = high_resolution_clock::now(); //Start the clock
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
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number

		D2Vec_Double Rho_dt(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double DRho(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.

		D2Vec_Double Rho_M(tot_iter, vector<double> (Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		//double Rho_M[tot_iter][Sp2] ={0.0}; //Stores <rho>x values per time step for a single replicate.

		double Rhox_avg[Sp] = {0.0}; double Rhox_num[Sp] = {0.0}; //Stores spatial average and number of occupied sites at given t.

		//init_randbistableframe(Rho_dt, g*g, a, a_c,  perc, chigh, clow); // Returns a frame with random speckles of high and low density.

		// Initialise vector <double> amp to elements of clow[].
		vector <double> amp(Sp, 500.0);  //Setting amplitude of gaussian distributions of vegetation to 500.
		for (int s=1; s< Sp; s++)
			amp[s] = clow[s+ Sp];
		
		vector <double> sd{g/8.0, g/16.0}; //Setting standard deviation of gaussian distributions of vegetation to g/8 and g/16.

		init_gaussframe(Rho_dt, g*g, sd, amp);
		



		// The first SpB species stored in Rho_dt are to be initialised on frame as Gaussian distributions.
		// Species 0 is the vegetation, Species 1 is the grazer and Species 2 is the predator.
		// Species 0 should be centered near the top right corner of the grid, Species 1 near the bottom right corner and Species 2 near the bottom left corner.

		

		/**
		
		stringstream  rain, rep_j, grid; 
		rain  << a; rep_j << int(j%2); grid << g;
    	string initcsv_filename= "Burning_In_Frames/FRAME_RAND_Basic_P_c_DP_G_" + grid.str() +"_T_91201_dt_0.2_alph_0.00833333_cgmax_0.0208333_a_" 
		+ rain.str() +"_dx_0.1_R_"  + rep_j.str() + ".csv";
		const vector<int> initcsv_columns = {2, 3, 4}; //Columns to be read from csv file.
		D2Vec_Double const_ind_val = {{1, dP/1000.0 }, {2, dP/5000.0}}; 
		init_csvconstframe(Rho_dt, const_ind_val, initcsv_filename, initcsv_columns, g*g); //Initialise Rho_dt with burn-in frame.
		
		vector<double> temp_vec= {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species 's'
		double init_veg_num = occupied_sites_of_vector(temp_vec, g*g); //Finds number of occupied at given t.
		vector<double>().swap(temp_vec); //Flush temp_vec out of memory.

		

		if(init_veg_num == 0.0)
		{
			stringstream m1_1;
			m1_1 << "RUN-TIME WARNING: ZERO Active BURN-IN veg sites for Thread Rank:\t " << omp_get_thread_num() 
			<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << " CONSIDER RE-RUNNING OVER A DIFFERENT RANGE." <<endl; 
			cout << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
		}

		stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_1 << "Initial Conditions: BURN-IN # Active Veg sites:\t"<< init_veg_num << " for Thread Rank:\t " << omp_get_thread_num() 
		<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << endl; 
		cout << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
		
		*/

		for(int i=0; i < g*g; i++)
		{
			for(int s=0; s < Sp; s++)
			{
				DRho[s][i] = Rho_dt[s][i]; //Assigning initial conditions
				Rho_tsar[s][i] =Rho_dt[s][i];
			}
		}
		
		stringstream m1_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_2 << "Initial Conditions:\t Per:\t" << perc  << "\t C_High[0], C_High[1], C_High[2]:\t " << chigh[0] << "," << chigh[1] << "," << chigh[2]   
		<< "\t C_Lo[0], C_Lo[1], C_Lo[2]:\t " << clow[0] << "," << clow[1] << "," << clow[2] << ",\t R: " << a << endl; //cout << m1.str();
		if(omp_get_thread_num()== 1)
			cout << m1_2.str();  
		errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close(); 
		

		poisson_distribution<int> poisson; gamma_distribution <double> gamma_distr; 
		uniform_real_distribution<double> unif(0.0, 1.0); normal_distribution<double> norm(0.0, 1.0);

		double t=0; int index = 0;  //Initialise t
		int iter_index=0; int po=1; int so =1; int lo=1; int counter =0;

		while( t < t_max + dt)
		{
			// Basic Book-keeping below. Prints frames out to file at given time points. 
			//Also updates number and avg density of active sites at given time points.
			if(t == 0 || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{
				double rhox_avg, rhox_num; //Stores spatial average and number of occupied sites at given t.
				for(int s=0; s<Sp; s++)
				{
					vector <double> temp= {Rho_dt[s].begin(),Rho_dt[s].end()}; //Rho_dt for species 's'
                    rhox_avg = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
				    rhox_num = occupied_sites_of_vector(temp, g*g); //Finds number of occupied at given t.
				    Rho_M[index][2*s] = rhox_num; Rho_M[index][2*s +1] = rhox_avg; 

					vector<double>().swap(temp); //Flush temp out of memory.0

                    //if(rhox_num == 0.0 and s== 0)
					//    break; //Finish negotiations as all basal species are dead.
				}

				vector <double> temp_alt= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
				rhox_num = occupied_sites_of_vector(temp_alt, g*g); //Finds number of occupied at given t.
				vector<double>().swap(temp_alt); //Flush temp out of memory.0

				if( t >= 15 && t <= 10000 || t== 0 ||  index >= tot_iter -8 &&  index <= tot_iter-1)
				{
					//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.
					
					stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, Dm2, alph, w0t, aij, hij, dix, dimitri, sig0, sig1, geq;

  					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  					rini << j; Dm0 << D[0]*pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2]; gm << setprecision(3) << gmax;
					w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; dimitri  << dP;
					a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1]; geq << setprecision(5) << Gstar;
					aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1];
					// Three replicates are over.

					//ofstream frame_dp;


					save_frame(Rho_dt, a, a_st, a_end, c, gmax, alpha, rW, W0, D, K, sigma, A, H, E, t, dt, dx, dP, j, g, Gstar);

					stringstream m3;
					m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
					cout << m3.str(); errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
					/**

					if(t < 760)
					{
						//Only save one in three frames here.
						if(index%3 ==2 || t== 0)
						{
							save_frame(Rho_dt, a, a_st, a_end, c, gmax, alpha, rW, W0, D, K, sigma, A, H, E, t, dt, dx, dP, j, g, Gstar);

							stringstream m3;
							m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
							cout << m3.str(); errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
						}
					}
					else if(t >= 760 && t < 6000)
					{	//Only save one in two frames here.
						if(index%2 ==0)
						{
							save_frame(Rho_dt, a, a_st, a_end, c, gmax, alpha, rW, W0, D, K, sigma, A, H, E, t, dt, dx, dP, j, g, Gstar);
							stringstream m3;
							m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
							cout << m3.str(); errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
						}
					}
					else
					{
						//Save all frames here.
						save_frame(Rho_dt, a, a_st, a_end, c, gmax, alpha, rW, W0, D, K, sigma, A, H, E, t, dt, dx, dP, j, g, Gstar);
						stringstream m3;
						m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
					}

					*/
				}

				vector <double> temp_1= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
				double rhox_DR = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
				double  rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				index+=1;
				/**
				stringstream m4_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m4_1 << "AFTER SAVING FRAME, STATUS UPDATE AT TIME [t,i]\t" << t_meas[index] << " with # non-zero sites in VEG Rho_dt_0:  " << rhox_DR 
				<< " AND # Non-Zero in VEG DRho_0:  " << rhox_dt <<  endl; cout << m4_1.str();
				errout.open(thr, std::ios_base::app); errout << m4_1.str(); errout.close(); **/
				//int counter_nzero_alpha=0;
			}

			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC. ONLY LIVING MATTER EXPERIENCES DORNIC INTEGRATION.
			//First up vegetation.

			for(int s=0; s< 1; s++)
            {
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					double alpha_i = diff_coefficient[s]*((DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
					DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]]));

					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						/**
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "MUCH VEG TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i << " \t with alpha_i, Rho(i), DRho(i):  "
						<< alpha_i << " , " << Rho_dt[s][i] << " , " << DRho[s][i] << "\t and Diff Const (D0/dx2): " << setprecision(16) << diff_coefficient[s] << endl; cout << m6.str();
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close(); */
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					//counter_nzero_alpha +=1;

					double mu = -1.0 + 2.0*alpha_i*sigma2_1[s]; // NOTE: sigma2_1 = 1/(sigma^2)
					double ziggy = lambda[s]*lambda_exp[s]*Rho_dt[s][i]; 
					double gru; // Stores the Gamma random variable.
					if(ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if(ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng)*sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;

						/**
						stringstream m6_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6_2 << "LARGE ZIGGY: " << ziggy << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
						<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and Gaussian Approx Value:\t" << gauss << "\t and GRU:\t" << gru << endl; 
						cout << m6_2.str(); errout.open(thr, std::ios_base::app); errout << m6_2.str(); errout.close(); */
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if(gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng)*sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i]= gauss/lambda[s];
						/**
						stringstream m6_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6_2 << "LARGE GRU: " << gru << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
						<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and Gamma (Gaussian) Approx Value:\t" << gauss << "\t and VEG RHO*:\t" << Rho_dt[s][i] << endl; 
						cout << m6_2.str(); errout.open(thr, std::ios_base::app); errout << m6_2.str(); errout.close(); */
					}
					else
					{
						gamma_distr = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i]= gamma_distr(rng)/lambda[s];
					}
					/**
					stringstream m6_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
					m6_1 << "DORNIC VEG RHO* " << Rho_dt[s][i] << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
					<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i  << "\t and POISSON ENTRY:\t" << ziggy << "\t and GAMMA ENTRY:\t" << gru << endl; 
					cout << m6_1.str(); errout.open(thr, std::ios_base::app); errout << m6_1.str(); errout.close(); **/
					if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "VEG HIT ROCK BOTTOM WITH: Rho*"<< s <<"[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i << 
						"\t and GAMMA ENTRY:\t" << gru << endl; cout << m6.str();
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}
					Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
				}
			} // End of Vegetation Integration
			
			// Book-keeping for determining gamma_i for higher order species. Calculating Rho averages per species at each time step.
			for(int s=0; s < Sp -2; s++)
			{
				vector <double> temp= {DRho[s].begin(),DRho[s].end()}; //Rho_dt for species '0'
				Rhox_avg[s] = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
				vector<double>().swap(temp); //Flush temp out of memory.0
			}
			vector <double> temp_veg= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
			double rhox_num_veg = occupied_sites_of_vector(temp_veg, g*g); //Finds number of occupied sites at given t.
			vector<double>().swap(temp_veg); //Flush temp out of memory. 
			double nR_fac = 1 - rhox_num_veg/(g*g); //Factor to reduce the number of neighbours for gamma_i estimation
			if (nR_fac < 0.35)
			{	nR_fac = 0.35; }

			/**
			if(counter%50000 == 0)
			{

				temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
				rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, thr]\t" << t << " , " << omp_get_thread_num()  << "\t VEG DORNIC DONE with non-zero elements in DRho_0:  " 
				<< rhox_num_veg << "\t and non-zero elements in Rho_dt_0:  " << rhox_dt << "\t and INDEX VAL, T_MEAS:\t" << index <<  " , " << t_meas[index] << endl; cout << m5_1.str();
				errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			*/

			D2Vec_Double gamma(SpB, vector<double> (g*g, 0.0)); //Stores gamma for each species at each site.

			// Now for higher order species.
			for(int s=1; s< Sp-2; s++)
            {
				if(Rhox_avg[s] == 0.0)
					continue;
				//If the average density of a species is 0, then skip the Dornic integration for that species.

				
				for(int i=0;i<Rho_dt[0].size();i++)
				{

					//Firstly gamma needs to be estimated for each site.
					if(s ==1)
						calc_gamma_3Sp(origin_Neighbourhood, DRho, gamma_i, Rhox_avg, r_frac, nR_fac, r_max_effective, i, g);
					// This function calculates gamma for each species at each site i in an array indexed by species.
					
					//Finally check if gamma_i > 1, if so set it to 1.
					if(gamma_i[s] > 1)
					{	gamma_i[s] = 1; }
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j 
						<< "\t GAMMA FOUND WITH GAMMA_I:" << gamma_i[s] << endl; cout << m5_2.str();
						errout.open(thr, std::ios_base::app); errout << m5_2.str(); errout.close();
					}
					if( gamma_i[s] < 0.0)
					{
						stringstream m5_3;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_3 << "GAMMA_I FALLS BELOW 0,  AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j  
						<< "\t WITH GAMMA_I:" << gamma_i[s] << " for species: " << s << endl; cout << m5_3.str();
						errout.open(thr, std::ios_base::app); errout << m5_3.str(); errout.close();
					}
					if(isnan(gamma_i[s]) == true || isinf(gamma_i[s]) == true)
					{
						stringstream m5_3;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_3 << "GAMMA_I BLOWS THE LID,  AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j  
						<< "\t WITH GAMMA_I:" << gamma_i[s] << " for species: " << s << " with Rho_avg[0]: " << Rhox_avg[0] << " and Rho_avg[1]: " << Rhox_avg[1] << " and Rho_avg[2]: " << Rhox_avg[2] 
					  	<< " \n and Rho_t[0][i]: " << Rho_dt[0][i] << " and Rho_t[1][i]: " << Rho_dt[1][i] << " and Rho_t[2][i]: " << Rho_dt[2][i]<< endl; cout << m5_3.str();
						errout.open(thr, std::ios_base::app); errout << m5_3.str(); errout.close();
						save_frame(Rho_dt, a, a_st, a_end, c, gmax, alpha, rW, W0, D, K, sigma, A, H, E, t, dt, dx, dP, j, g, Gstar);
					}

					//Next sampling the advection vector from v[s].
					//unif = uniform_real_distribution<double>(0.0, 1.0);
					double ran = unif(rng);
					double vx = v[s]*cos(2*PI*ran); double vy = v[s]*sin(2*PI*ran); //Advection vector.
					double vx_abs = vx*sgn(vx); double vy_abs = vy*sgn(vy); //Absolute value of advection vector.

					//RECALL: dx1_2 = 1/(dx2); dxb2 = dx/2.0;
					diff_eff[s].first = (D[s] -(dxb2)*(vx_abs*(1 - dt*dx1*vx_abs)))*dx1_2;
					diff_eff[s].second = (D[s] -(dxb2)*(vy_abs*(1 - dt*dx1*vy_abs)))*dx1_2;
					//diff_coefficient[s] = (D[s] -(dx/2.0)*(v[s]*(1 - dt*vdx[s])))/(dx2);
					// Advection leads to excess diffusion, hence the correction term.

					if(diff_eff[s].first < diff_coefficient[0])
						diff_eff[s].first = diff_coefficient[0];
					
					if(diff_eff[s].second < diff_coefficient[0])
						diff_eff[s].second = diff_coefficient[0];
					// If the effective diffusion coefficient is less than D[0]/dx2, set it to D[0]/dx2.
					

					// RECALL:	nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 
					// NOTE: dx1 is the inverse of dx. sgn_index(vx) returns 0 if vx is positive, else 1.
					// Recall nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1;
					// ALSO NOTE: diff_coefficient[s] = D[s]/(dx*dx); diff_eff[s].first = D_X_eff[s]/(dx*dx); diff_eff[s].second = D_Y_eff[s]/(dx*dx);
					
					// Standard diffusion term.
					//double alpha_i = gamma_i[s]*diff_coefficient[s]*(1*(DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
					//				DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]])) 
					//				+dx1*(1-gamma_i[s])*(vx_abs*DRho[s][nR2[i][1][sgn_index(vx)]]+ vy_abs*DRho[s][nR2[i][0][sgn_index(vy)]]);

					// Advection corrected diffusion term.
					double alpha_i = gamma_i[s]*(diff_eff[s].second*(DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]]) 
												+ diff_eff[s].first*(DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]])) 
									+dx1*(1-gamma_i[s])*(vx_abs*DRho[s][nR2[i][1][sgn_index(vx)]]+ vy_abs*DRho[s][nR2[i][0][sgn_index(vy)]]);
					
					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					//Beta depends on the gamma_i value at each lattice site for each species, hence needs to be updated at each iteration.
					// Standard diffusion term.
					//beta[s] = -M[s] - 4*gamma_i[s]*diff_coefficient[s] - (1-gamma_i[s])*(vx_abs + vy_abs)*dx1; // NOTE: diff_coefficient[s] = D[s]/(dx*dx)
					// Advection corrected diffusion term.
					beta[s] = -M[s] - 2*gamma_i[s]*(diff_eff[s].first + diff_eff[s].second) - (1-gamma_i[s])*(vx_abs + vy_abs)*dx1;

					lambda_exp[s]= exp( (beta[s])*dt); lambda[s] = 2*(beta[s])*(sigma2_1[s])/(lambda_exp[s] -1.0); // NOTE: sigma2_1 = 1/(sigma^2);

					double mu = -1.0 + 2.0*alpha_i*sigma2_1[s]; // NOTE: sigma2_1 = 1/(sigma^2)
					double ziggy = lambda[s]*lambda_exp[s]*DRho[s][i]; 
					double gru; // Stores the Gamma random variable.
					if(ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if(ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng)*sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if(gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng)*sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i]= gauss/lambda[s];
					}
					else
					{
						gamma_distr = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i]= gamma_distr(rng)/lambda[s];
					}
					if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "YOU HIT ROCK BOTTOM WITH: Rho*"<< s <<"[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i << 
						"\t ,GAMMA ENTRY:\t" << gru << "\t and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" <<endl; cout << m6.str();
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
						exit(3);
					}
					Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
					
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j 
						<< "\t DORNIC VALUE OF SP " << s << " UPDATED WITH RHO*(I,t): " << Rho_tsar[s][i] 
						<< " and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" << endl; cout << m5_2.str();
						errout.open(thr, std::ios_base::app); errout << m5_2.str(); errout.close();
					}

					gamma[s][i] = gamma_i[s]; //Store gamma for each species at each site.
					
				}
			}

			if(t==0 || t >= t_meas[index-1] -dt/2.0 && t < t_meas[index-1] +dt/2.0)
			{
				//Save the gamma values for each species at each site at each time point.
				stringstream L, tm ,d3, p1, rini, a1, a2, Dm1, Dm2, dix, dimitri, geq;

  				L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  				rini << j; Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2];
				a1 << a_st; a2  << a_end; dimitri  << dP; geq << setprecision(5) << Gstar;

				string filenamePattern = gamma_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str()  + "_D1_"+ Dm1.str() + "_D2_"+ Dm2.str() + "_dx_"+ dix.str() + "_R_";
				
				string parendir = "";
				// Creating a file instance called output to store output data as CSV.
				if(Gstar != -1)
					parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str();
				else
					parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();

				string gammaheader = "a_c,  x,  Gam0(x; t), Gam1(x; t), Gam2(x; t) \n"; //Header for output frame.

				string filename = filenamePattern + rini.str() + ".csv";
				string filepath = parendir + "/" + filename;
				ostringstream oss;
				oss << gammaheader;
				for(int i=0; i<g*g; i++)
				{
					oss << a << ", " << i << ", " << gamma[0][i] << ", " << gamma[1][i] << ", " << gamma[2][i] << "\n";
				}

				FILE *gammafile = fopen(filepath.c_str(), "w");
				if(gammafile)
				{
					fprintf(gammafile, "%s", oss.str().c_str());
					fclose(gammafile);
				}
				else
				{
					std::cerr << "Error: Could not open file " << filepath << " for writing." << std::endl;
				}

			}

			//Swap the gamma vector to free up memory.
			vector<vector<double>>().swap(gamma);
			
			if(counter%50000 == 0)
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j << "\t GRAZER + PRED DORNIC DONE." << endl; 
				cout << m5_1.str(); errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			

			// Finally RK4 integration of the remaining terms.

			RK4_Integrate_Stochastic_MultiSp(Rho_dt, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t, dt, dx, g);
			
			if(counter%50000 == 0)
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num()  << " , " << j << "\t RK4 INTEGRATION OF REMAINDER TERMS DONE." 
				<< endl; cout << m5_1.str(); errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			


			for(int s=0; s< Sp; s++)
			{
					
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					if( Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true and so==1)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO "<<s<<"  FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i] << "\t and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" <<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}

					if( isnan(Rho_dt[s][i] == true))
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO "<<s<<"  NAAN PARATHA:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<< "\t and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" <<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
						lo == -1;
					}
					DRho[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
					Rho_tsar[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
				}
			}

			t+=dt; //Update timestep
			counter+=1; //Update counter

			if(counter%40000 ==1)
			{
				//Find memory usage and time taken for every 50000 iterations.
				//size_t currentSize = getCurrentRSS( ); //Check Memory usage.
				//size_t peakSize    = getPeakRSS( );
				auto end_t = high_resolution_clock::now();
				auto elapsed_min = duration_cast<minutes>(end_t-start_t);
				// Reset the clock
				start_t = high_resolution_clock::now();
				stringstream m7;
				m7 << "STATUS UPDATE AT [t, a, thr, j]\t" << t << " , "  << a << " , " << omp_get_thread_num() << " , " << j 
				//<< " ,with current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) 
				<< " and time taken for " << counter << " iterations: " << elapsed_min.count() << " min." <<  endl; cout << m7.str();
				errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();

				//Report syst every 50000 iterations.
			

			}

			if(lo == -1)
			{
				exit(3);
			}
		} // End of time loop.
		
		vector<vector <double>>().swap(Rho_dt); vector<vector <double>>().swap(DRho); vector<vector <double>>().swap(Rho_tsar);
		// Freeing up memory.
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
				for(int s =0; s< Sp-2; s++)
				{
					rho_rep_avg_var[i][4*s+ 3] = Rho_M[i][2*s + 1]; // Avg density of frame at given time point i.
					rho_rep_avg_var[i][4*s+ 4] = 0.0; // Variance in avg frame density.
					rho_rep_avg_var[i][4*s+ 6] = Rho_M[i][2*s]; // Number of active sites in frame at given time point i.

					if(Rho_M[i][2*s] > 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s+5 ] = 1;  } //One surviving run as of now
					else if(Rho_M[i][2*s] == 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s + 5] = 0;  } //No surviving run as of now
				}
				rho_rep_avg_var[i][1] = Rho_M[i][2*(Sp-2) + 1]; // Avg density of frame at given time point i.
				rho_rep_avg_var[i][2] = Rho_M[i][2*(Sp-1) + 1]; // Avg density of frame at given time point i.
				
			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			// Second or higher replicate, use incremental advances.
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter, j); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		// Freeing up memory.
		vector<vector <double>>().swap(Rho_M); 

		if((j+1)%2 == 1 || (j+1)%2 == 0)
		{	//Start logging at every multiple of 2
			stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, rini_prev, Dm, cgm,  sig0, geq, jID;
			
			a1 << a_st; a2 << a_end;
  			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(5) << a; dimitri << dP;
  			rini << j; Dm << setprecision(4) << D[2]; cgm << c*gmax; sig0 << sigma[0]; geq << setprecision(5) << Gstar;
			// Three replicates are over.

			unif = uniform_real_distribution<double>(0.0, 1.0);
			double ran_jid = (unif(rng)*1000.0)/1000.0; // Random number between 0 and 1.

			jID << ran_jid;
			
			string path_to_dir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str();
			string filename = "/PRELIM_AGGRAND_P_c_ID_"+ jID.str() +"_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
			"_D2_"+ Dm.str() + "_R_"+ rini.str() + ".csv";

			std::ostringstream oss;
			oss << " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs Rho0,"
			" # Active Sites Rho0, <<Rho1(x; t)>_x>_r, Var[<Rho1(t)>_x]_r, # Surviving Runs Rho1, # Active Sites Rho1," 
			"<<Rho2(x; t)>_x>_r, Var[<Rho2(t)>_x]_r, # Surviving Runs Rho2, # Active Sites Rho2, \n";

			for(int i=0; i< tot_iter; i++)
			{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |
				oss << a << ","<< j << "," << g << "," << rho_rep_avg_var[i][0] << "," << rho_rep_avg_var[i][1] << "," << rho_rep_avg_var[i][2] << ",";
				for(int s=0; s <Sp-2; s++)
				{
					oss << rho_rep_avg_var[i][4*s + 3] << "," << rho_rep_avg_var[i][4*s +4] << "," << rho_rep_avg_var[i][4*s +5] << "," << rho_rep_avg_var[i][4*s +6] << ","; 
				}
				oss <<  "\n";
			}

			FILE *fp = fopen((path_to_dir + filename).c_str(), "w");
			if(fp)
			{
				fprintf(fp, "%s", oss.str().c_str());
				fclose(fp);
				stringstream m7;
				m7 << "PRELIM File " << filename << " written successfully." << std::endl; cout << m7.str();
				errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();
			}
			else
			{
				stringstream m7;
				m7 << "Error writing to file: " << path_to_dir + filename << " with error: " << strerror(errno) << std::endl; 
				cout << m7.str(); cerr << m7.str(); errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();
			}

			// Delete the previous file (if it exists) to save space.
			#if defined(__GNUC__) && (__GNUC__ >= 9)

			if(j > 2)
			{
				// Delete the previous file (if it exists) to save space.
				// Reset rini value.
				rini_prev << j-2;

				// Using regex to delete the file.
				string filenamePattern = "PRELIM_AGGRAND_P_c_ID_.*_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
				"_D2_"+ Dm.str() + "_R_"+ rini_prev.str() + ".csv";

				std::regex filePattern(filenamePattern);

				// Iterate over the files in the parent directory (path_to_dir) and delete the file if it matches the pattern.
				for (const auto & entry : fs::directory_iterator(path_to_dir))
				{
					if (std::regex_match(entry.path().filename().string(), filePattern))
					{
						if(fs::exists(entry.path()) && fs::is_regular_file(entry.path()))
						{
							try
							{	
								fs::remove(entry.path());
								stringstream m8;
								m8 << "File " << entry.path().string() << " deleted successfully." << std::endl; cout << m8.str();
								errout.open(thr, std::ios_base::app); errout << m8.str(); errout.close();
							}
							catch (const fs::filesystem_error& e)
							{
								stringstream m8;
								m8 << "Error deleting file: " << entry.path().string() << " with error: " << e.what() << std::endl; cout << m8.str();
								cerr << m8.str(); errout.open(thr, std::ios_base::app); errout << m8.str(); errout.close();
							}
						}
						else
						{
							stringstream m8;
							m8 << "File: " << entry.path().string() << " does not exist." << std::endl; cout << m8.str();
							cerr << m8.str(); errout.open(thr, std::ios_base::app); errout << m8.str(); errout.close();
						}
						
					}
				}
					
			}
			#endif

    
		}

	} // End of r loop.

	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|    r 		|     L 		|    t 		|     <<W(t)>>x,r			|      <<O(t)>>x,r		|
		//..     <<Rho0(t)>>x,r			|    Var0[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |....

		Rho.push_back({a, static_cast<double>(r), static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], 
					rho_rep_avg_var[i][2]});
		for(int s=0;s < Sp-2; s++)
		{
			Rho[i].insert(Rho[i].end(), { rho_rep_avg_var[i][4*s +3], rho_rep_avg_var[i][4*s +4], 
			rho_rep_avg_var[i][4*s +5], rho_rep_avg_var[i][4*s +6]});		
		}
	}
	stringstream m10;
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m10.str();
	errout.open(thr, std::ios_base::app); errout << m10.str(); errout.close();
}


void first_order_critical_exp_delta_stochastic_3Sp(int div, double t_max, double a_start, double a_end, double a_c,  double c, double gmax, double alpha,
	double rW, double W0,  double (&D)[Sp], double v[], double (&K)[3], double sigma[], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double pR[], double ch[], double clo[],
	double dt, double dx, double dP, int r,  int g)
{
	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//double mean[Sp] = {1.0, 0.05}; double sd[Sp] = {0.25, 0.0125};
	//init_randframe(Rho_0, Sp,  g*g); //Returns Rho_0 with a full initial frame filled with 0.2.

	set_density_dependentmortality(M[2]); //Sets the density dependent mortality rate.

	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "NOSTRA" <<endl;
  	// The pspace to iterate over.
    //vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  	// Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

	vector <double> t_measure = linspace(40.0, 800.0, 20); // Computes and returns linearly distributed points from t= 20 to 240 (in steps of 20) 
  	// Returns time-points measured on a linear scale from 20 to 240 in steps of 20.
	t_measure.insert(t_measure.begin(), 0.0);

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
	if(nProcessors > 32)
	{
		//omp_set_num_threads(32); //Limiting use on Chunk. Don't be greedy.
		cout << "WARNING: Number of Processors is: " << nProcessors << ". Don't be greedy." << endl;
	}

	//double perc = 0.015; double c_high[Sp] ={dP, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

  	  
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {

		//init_randconstframe(Rh0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho1(t)>>x,r			|    Var[<Rho1(t)>x],r    |
		**/
		rietkerk_Dornic_2D_MultiSp(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, rW, W0, D, v, K , sigma, a_start, a_end, a_c, A,H,E,M, pR, ch, clo, dt, dx, dP, r, g);
		//RK4_Wrapper_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , a_start, a_end, dt, dx, dP, r, g);
		//expanded_percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_space[i], b, c, D, sigma, dt, dx, r, g);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);
		//void rietkerk_Dornic_2D_MultiSp(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0,  double D[], double v[], double K[], double sigma[], double a_st, double a_end, double a_c, double A[SpB][SpB], double H[SpB][SpB], double E[], double M[], double pR[], chigh[],  clow[],  dt, dx, dP,  r, g)

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
			 
          message3 << "Is this happening?\n" << "In thread: " << i << ", current size in MB: " << currSize/(1024.0*1024.0) 
		  << " and Peak Size (in MB): " << peaksSize/(1024.0*1024.0) << endl;
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	//Store a uniformly distributed random number between 0 and 1 (rounded to 3 decimal places).
	double id = round(unif(rng)* 1000.0) / 1000.0; //Random number between 0 and 1.
	//Round to 3 decimal places.
	stringstream L, coco, tm ,d3, p1, p2, rini, Dm0, Dm1, aij, hij, cgm, alph, dix, dimitri, Sig0, geq, ID;

	double Gstar = M[2]/((E[2] -M[2]*H[1][2])*A[1][2]); // MFT estimate of Grazer density at coexistence.

  	L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << a_start; p2  << a_end; rini << r; 
	alph << alpha; cgm << c*gmax; coco << setprecision(4) << c;
  	Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; Sig0 << sigma[0]; dix << setprecision(2) << dx; 
	aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; geq << setprecision(5) << Gstar;
	ID << id;

	ofstream output_1stdp;
  // Creating a file instance called output to store output data as CSV.
	output_1stdp.open(stat_prefix + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D1_"+ Dm1.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + 
	"_Geq_"+ geq.str() /* + "_ID_" + ID.str() */ + "_R_"+ rini.str() + ".csv");

	cout << "Save file name: " <<endl;
	cout << stat_prefix + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D1_"+ Dm1.str() + "_S0_" + Sig0.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + 
	"_cgmax_"+ cgm.str() + "_Geq_"+ geq.str() /* + "_ID_" + ID.str() */ + "_R_"+ rini.str() + ".csv" <<endl;

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

	string header = " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs Rho0, # Active Sites Rho0,"
	" <<Rho1(x; t)>_x>_r, Var[<Rho1(t)>_x]_r, # Surviving Runs Rho1, # Active Sites Rho1, <<Rho2(x; t)>_x>_r, Var[<Rho2(t)>_x]_r, # Surviving Runs Rho2, # Active Sites Rho2";
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


	



