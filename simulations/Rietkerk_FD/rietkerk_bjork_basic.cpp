#include "rietkerk_bjork_basic.h"
#include "Debug.h"

//double Mm; // Density dependent mortality of predator.

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


#if defined(BARRACUDA)
//Helper functions for type casting
// Scalar specialization
template <typename Src, typename Dst>
__host__ __device__ inline void cast_value(const Src &src, Dst &dst) {
    dst = static_cast<Dst>(src);
}

// Pair specialization
template <typename SrcT1, typename SrcT2, typename DstT1, typename DstT2>
__host__ __device__ inline void cast_value(const std::pair<SrcT1, SrcT2> &src,
                                           std::pair<DstT1, DstT2> &dst) {
    dst.first  = static_cast<DstT1>(src.first);
    dst.second = static_cast<DstT2>(src.second);
}

// Copy a vector of type Tsrc with N elements to device array of type Tdst
template <typename Tsrc, typename Tdst>
void transfer_castedVec_to_device(const std::vector<Tsrc> &host_src, Tdst* d_dst_out, size_t N) {
    //static_assert(std::is_arithmetic<Tsrc>::value && std::is_arithmetic<Tdst>::value,
    //              "upload_casted_to_device: only numeric types are supported");

    //const size_t N = host_src.size();
    std::vector<Tdst> tmp(N);

    // Explicit cast (avoids narrowing warnings)
    for (size_t i = 0; i < N; ++i)
        cast_value(host_src[i], tmp[i]);
		//tmp[i] = static_cast<Tdst>(host_src[i]);

    cudaError_t err = cudaMemcpy(d_dst_out, tmp.data(), sizeof(Tdst) * N, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        throw std::runtime_error("cudaMemcpy (H-->D) failed: " + std::string(cudaGetErrorString(err)));
    }
}

template <typename Tsrc, typename Tdst>
void transfer_castedVec_to_host(const Tsrc *d_src, std::vector<Tdst> &host_dst, size_t N) 
{
    //static_assert(std::is_arithmetic<Tsrc>::value && std::is_arithmetic<Tdst>::value,
    //              "download_casted_to_host: only numeric types are supported");
    //const size_t N = host_dst.size();
    std::vector<Tsrc> tmp(N);

    cudaError_t err = cudaMemcpy(tmp.data(), d_src, sizeof(Tsrc) * N, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        throw std::runtime_error("cudaMemcpy (D-->H) failed: " + std::string(cudaGetErrorString(err)));
    }

    for (size_t i = 0; i < N; ++i)
        cast_value(host_dst[i], tmp[i]);
		//host_dst[i] = static_cast<Tdst>(tmp[i]);
}

#endif


//--------------------Defining CUDA Context Kernel------------------------------------//



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

// Function to replace placeholders in the string with corresponding values from the map
// Placeholders are in the format {key}, where key is the key in the map
string format_str(const std::string& s, const std::map<string, string>& values) 
{
    string result = s;
    std::regex placeholder_pattern(R"(\{([a-zA-Z0-9_]+)\})");  // Regex to match {key}
    std::smatch match;

    // Find all placeholders in the string
    while (std::regex_search(result, match, placeholder_pattern)) {
        std::string placeholder = match[0];  // Full match, e.g., "{Pre}"
        std::string key = match[1];          // Capture group, e.g., "Pre"

        // Replace with corresponding value from the map if key exists
        auto it = values.find(key);
        if (it != values.end()) {
            result.replace(match.position(0), placeholder.length(), it->second);
        } else {
            throw std::runtime_error("Key not found: " + key);
        }
    }

    return result;
}

// Function to precompute FFTs of neighborhood kernels for different radii R.
// These FFTs are stored in a cache indexed by R for efficient reuse during simulations.
// 
void precompute_FFT_KHat_kernels(FFTW3_CentralPlanner& fft_centralplanner, int L, int Rmax, double nVeg_frac_min /*= 0.2 */, double nVeg_frac_max /*= 1.0 */)
{
    int L2 = L * L;
    //int new_length_max = int(nVeg_frac_max * nVeg_frac_max * centralNeighboringSites.size());
	//int new_length_min = int(nVeg_frac_min * nVeg_frac_min * centralNeighboringSites.size());
    //if (new_length_max < 1) new_length_max = 1; if (new_length_min < 1) new_length_min = 1;
    
	for(int R = 0; R <= Rmax; R++)	
    //for (int eff_nR_Perp_size = new_length_min; eff_nR_Perp_size <= new_length_max; eff_nR_Perp_size++)	
	{
		double frac = double(R) / double(Rmax);
		// Compute sorted neighboring site indices (centralized around (0,0)) for current R.
		std::vector<std::pair<int, int>> originK_R = computeNeighboringSitesCentral(R);
		int eff_nR_Perp_size = originK_R.size();
		//int eff_nR_Perp_size = int(frac * frac * centralNeighboringSites.size());
		if (eff_nR_Perp_size < 1) eff_nR_Perp_size = 1;
		// Check if kernel for this radius is already computed
		//if (kHat_kernel_FFT_cache.find(R) != kHat_kernel_FFT_cache.end()) {
		//	continue; // Kernel already computed, skip
		//}
		
		kHat_kernel_FFT_cache[R] = nullptr; // Initialize cache entries to nullptr
		// Create kernel in spatial domain which is 1 at each (periodic) neighbor location, 0 elsewhere.
        std::fill(fft_centralplanner.kernel_spatial, fft_centralplanner.kernel_spatial + L2, 0.0);
        // Build kernel: 1/N at each neighbor location, accounting for periodic boundaries
        double kernel_weight = 1.0/ eff_nR_Perp_size;
        for (int k = 0; k < eff_nR_Perp_size; k++) 
		{
            int dx = originK_R[k].first;
            int dy = originK_R[k].second;
            // Convert to array index with periodic wrapping
            int kx = (dx + L) % L;
            int ky = (dy + L) % L;
            fft_centralplanner.kernel_spatial[kx * L + ky] = kernel_weight;
        }
        // FFT of the kernel, and the forward plan creation.
        fftw_complex* kernel_fft = fftw_alloc_complex(L2);
		fftw_plan forward_kernel = fftw_plan_dft_r2c_2d(L, L, fft_centralplanner.kernel_spatial, kernel_fft, FFTW_ESTIMATE);
		fftw_execute(forward_kernel);
        //fftw_execute_dft_r2c(fft_centralplanner.forward_kernel, fft_centralplanner.kernel_spatial, kernel_fft);
        kHat_kernel_FFT_cache[R] = kernel_fft;
		// Map R to eff_nR_Perp_size for later use.
		nR_Perp_R[R] = eff_nR_Perp_size;
		fftw_destroy_plan(forward_kernel); // Clean up plan after execution
    }
}

bool comparePairs(const pair<int, int>& a, const pair<int, int>& b) 
{
    return (a.first * a.first + a.second * a.second) < (b.first * b.first + b.second * b.second);
}

void set_Prefix(string& user_prefix, double mG /** = 0.0 */, double mP /** = 0.0 */, bool stochastic /** = true */ )
{
	//Sets the global prefix for the output files.
	prefix = user_prefix;
	// If mG and mP are integers, convert them to integers, else round them to 2 decimal places.
	mG = (mG == int(mG)) ? int(mG) : round(mG * 100) / 100; 
	mP = (mP == int(mP)) ? int(mP) : round(mP * 100) / 100;
	// This is done to ensure that the string representation of mG and mP are consistent across different functions.
	ostringstream mG_str, mP_str;
	// Format the string representation of mG and mP to 2 decimal places at most and remove trailing zeros.
	mG_str << std::defaultfloat << mG; mP_str <<  std::defaultfloat << mP;
	string SpB_str ="";
	// If non-zero values of mG and mP are provided, include them in the directory structure.
	if(mG > 0.0 && mP > 0.0)
	  SpB_str = std::to_string(SpB) +"Sp_mG_" + mG_str.str() + "_mP_" + mP_str.str();
	else if(mG > 0.0 && mP == 0.0)
	  SpB_str = std::to_string(SpB) +"Sp_mG_" + mG_str.str();
	else
	  SpB_str = std::to_string(SpB) +"Sp";
	  if(stochastic)
	  {
		frame_folder = "../Data/Rietkerk/Frames/Stochastic/"+ SpB_str +"/" + prefix + "_";  //Folder to store frames.
	  	prelim_folder  = "../Data/Rietkerk/Prelims/Stochastic/"+ SpB_str +"/"+ prefix +"_"; //Folder to store preliminary data.
	  	stat_prefix =  "../Data/Rietkerk/Stochastic/1stOrderCC_Rietkerk_" + prefix + "_P_c_G_"; //Header for frame files.
	  }
	  else
	  {
		frame_folder = "../Data/Rietkerk/Frames/Non-Stochastic/"+ SpB_str +"/" + prefix + "_";  //Folder to store frames.
	  	prelim_folder  = "../Data/Rietkerk/Prelims/Non-Stochastic/"+ SpB_str +"/"+ prefix +"_"; //Folder to store preliminary data.
	  	stat_prefix =  "../Data/Rietkerk/Non-Stochastic/1stOrderCC_Rietkerk_" + prefix + "_P_c_G_"; //Header for frame files.
	  }

}


void set_LocustPrefix(string& user_prefix)
{
	//Sets the global prefix for the output files.
	prefix = user_prefix; 
	frame_folder = "../Data/LocRiet/Frames/Stochastic/"+ std::to_string(SpB) +"Sp/" + prefix + "_";  //Folder to store frames.
	prelim_folder  = "../Data/LocRiet/Prelims/Stochastic/"+ std::to_string(SpB) +"Sp/"+ prefix +"_"; //Folder to store preliminary data.
	stat_prefix =  "../Data/LocRiet/Stochastic/"+ std::to_string(SpB) +"Sp/1stOrderCC_Rietkerk_" + prefix + "_STOC_P_c_G_"; //Header for frame files.
}

void set_input_Prefix(string& user_inputframepath, string& user_prefix, double user_a_c, double user_dP, double user_Geq /* = -1*/, double user_input_T /* = -1*/)
{
	//Sets the global prefix for the INPUT files (for Burn-In Frames)
	input_frame_subfolder = user_inputframepath; //Folder to store input data.
	input_frame_parenfolder = "../Input/Rietkerk/"+ std::to_string(SpB) +"Sp"; //Folder to store input data.

	// Update map input_keys with the new values
	input_keys["outPRE"] = user_prefix;
	input_keys["a_c"] = std::to_string(user_a_c);
	input_keys["dP"] = std::to_string(user_dP);
	if(user_Geq >= 0)
		input_keys["Geq"] = std::to_string(user_Geq);

	if(user_input_T >= 0)
		input_keys["T"] = std::to_string(user_input_T);
	
}

void set_global_system_params(double dt, double dx)
{
	//Sets the global parameters for the simulation.
	//These are the parameters that are user-defined and are used in the simulation.
	//These are set by the user in the main function

	dt2 = dt/2.0; dt6= dt/6.0;
	dx2 = dx*dx; dx1 = 1.0/dx; dx1_2 = 1.0/(dx*dx); double dxb2 = dx/2.0;
}

void set_global_predator_params(double Km)
{
	/** Sets the global parameters for the predator.
	// Km is the carrying capacity of the predator.
	//These are set by the user in the main function **/
	K_P = Km; // Carrying Capacity of the predator
	K_P1 = 1.0/Km; // Inverse of the carrying capacity of the predator
}

/** // set_global_user_Rietkerk_params(...) Commented out for now, as it is not used in the current implementation.
void set_global_user_Rietkerk_params(double c,double gmax,double alpha, double rW, double W0, double (&D)[Sp], 
	double (&K)[3], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB])
{
	// Sets the global parameters for the Rietkerk model, as defined by the user.
	// This helps save computational time by avoiding repeated calculations.

	cgmax = c*gmax; K2W0 = K[2]*W0; 
	// For 2 Sp Rietkerk model, set only A01H01
	if (SpB >= 2)
	{
		A01H01 = A[0][1]*H[0][1];
	}
	// For 3 Sp Rietkerk model, set A12H12 as well.
	if (SpB >= 3)
	{
		A12H12 = A[1][2]*H[1][2];
	}
}
// */


void display_symbol_table(const exprtk::symbol_table<double>& symbol_table) {
    // Get variable names
    std::vector<std::string> variable_names;
    symbol_table.get_variable_list(variable_names);
    std::cout << "Variables:\n";
    for (const auto& var_name : variable_names) {
        if (symbol_table.symbol_exists(var_name)) {
			//double* var = symbol_table.get_variable(var_name);
			//auto& var = symbol_table.variable_ref(var_name);
            std::cout << "  " << var_name << " with value " << symbol_table.get_variable(var_name) << "\n";
        }
    }
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

void init_gradientframe(D2Vec_Double &array, int size, std::map<int, std::pair<double, double>> gradient, string orientation /*= "L"*/)
{
	//Initialises a frame (D2Vec_Double vector) with a gradient of values.
	// The gradient is defined by the map gradient, which contains as keys the indices of the species and as values the gradient range (min, max).
	// The orientation of the gradient is defined by the string orientation, which can be:
	// "T" (top to bottom) "B" (bottom to top) "L" (left to right) [Default] or "R" (right to left).
	int L = sqrt(size);

	//Iterating over keys in the map gradient
	for(auto const& [key, value] : gradient)
	{
		//Iterating over the range of the gradient
		for(int i=0; i< size; i++)
		{
			//Setting the value of the species at index i to the gradient value at index i.
			if(orientation == "T")
				array[key][i] = value.first + (value.second - value.first) * int(i/L) / L;
			else if(orientation == "B")
				array[key][i] = value.second - (value.second - value.first) * int(i/L) / L;
			else if(orientation == "L")
				array[key][i] = value.first + (value.second - value.first) * ( i % L) / L;
			else if(orientation == "R")
				array[key][i] = value.second - (value.second - value.first) * (i % L) / L;
			else
				array[key][i] = value.first + (value.second - value.first) * ( i % L) / L;
		}
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
		init_randconstframe(array, Sp, size, perc, c_high, c_low);
	
	else
		init_randconstframe(array, Sp, size, perc, c_high + Sp, c_low + Sp);
	
}

void init_exprtk_homogenousMFTframe(D2Vec_Double &array, int size, double R, double R_c, double c_spread[])
{
	/** SUMMARY: Same as init_exprtk_randbiMFTframe(), but with the following differences:
	 * dP and perc are not considered and all elements of the frame are set to the same value.
	 * These values are given by cspread[s]*MFT_Vec_CoexExpr[s] for each species s.
	 */

	mutex errorMutex; // Mutex to make error messages thread-safe
	// Create expression objects for each species
    std::vector<exprtk::expression<double>> expressions(2*Sp);
    exprtk::parser<double> parser;

	// Create a local symbol table, copying the global symbol table
	exprtk::symbol_table<double> local_symbol_table;// = global_symbol_table;

	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);

	uniform_real_distribution<double> unif;
	unif = uniform_real_distribution<double>(0.0, 1.0);
	int thr = omp_get_thread_num();
	//errorMutex.lock();
	local_symbol_table.add_variable("a", R);
	//errorMutex.unlock();

	// Compile the expressions
	for(int s=0; s< 2*Sp; s++) 
	{	
		std::lock_guard<std::mutex> lock(errorMutex);

		std::string expression_string = MFT_Vec_CoexExpr[s];
		expressions[s].register_symbol_table(local_symbol_table);
		if (!parser.compile(expression_string, expressions[s]))
		{
			
			std::cerr << "Error: Unable to compile expression for species " << s << " for [a, thr, expr]: " 
			<< R << " , " << thr << " , " << expression_string << "\n";
		}

		std::cerr << "Value of expression for species " << s << " for [a, thr, expr]: " 
			<< R << " , " << thr << " , " << expressions[s].value() << "\n";
	}

	//display_symbol_table(local_symbol_table);

	if (R < R_c)
	{
		for(size_t i=0; i< size; i++)
		{
			for(size_t s=0; s< Sp; s++)
			{
				// If the MFT condition before the critical point is 0, set the species to c_spread[s]
				// Otherwise, set the species to the MFT value multiplied by the constant values in c_spread
				if(expressions[s].value() == 0)
					array[s][i] = c_spread[s];
				else
					array[s][i] = c_spread[s]*expressions[s].value();
			}
		}
	}
	else if ( R >= R_c)
	{
		for(size_t i=0; i< size; i++)
		{
			// Set the other species to the constant values
			for(size_t s=0; s< Sp; s++)
			{
				if(expressions[s + Sp].value() != 0)
					array[s][i] = c_spread[s + Sp]*expressions[s+Sp].value();
				else
					array[s][i] = c_spread[s + Sp];
			}	
			// Set the other species to their MFT values multiplied by the constant values in c_spread
		}
	}
}


void init_exprtk_randbiMFTframe(D2Vec_Double &array, int size, double R, double R_c, double dP, double perc,  double c_spread[])
{

	/** Summary of init_exprtk_randbiMFTframe() function.
	 * @brief Initialises a frame (D2Vec_Double vector) as follows:
	 * If R < R_c, then the first Sp elements of c_spread are used to initialise the frame.
	 * In this case, the first row of the frame (corresponding to the first species) has a random distribution of values.
	 * (1 -perc)% of the values are set to 0 (c_spread[0]), and the remaining perc% are set to dP.
	 * The remaining rows of the frame (corresponding to the other species) are set to the constant values given by c_spread + 1 to c_spread + Sp-1.
	 * If R >= R_c, then the last Sp elements of c_spread are used to initialise the frame.
	 * In this case, the first row of the frame (corresponding to the first species) has a random distribution of values.
	 * (1 -perc)% of the values are set to MFT_Vec_CoexExpr[0], and the remaining perc% are set to c_spread[Sp]*MFT_Vec_CoexExpr[0] or dP, whichever is greater.
	 * The other "i" rows of the frame are set to the constant values given by MFT_Vec_CoexExpr[i]*c_spread[i+Sp]  where i ranges from 1 to Sp-1.
	 */

	mutex errorMutex; // Mutex to make error messages thread-safe
	// Create expression objects for each species
    std::vector<exprtk::expression<double>> expressions(2*Sp);
    exprtk::parser<double> parser;

	// Create a local symbol table, copying the global symbol table
	exprtk::symbol_table<double> local_symbol_table;// = global_symbol_table;

	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);

	uniform_real_distribution<double> unif;
	unif = uniform_real_distribution<double>(0.0, 1.0);
	int thr = omp_get_thread_num();
	//errorMutex.lock();
	local_symbol_table.add_variable("a", R);
	//errorMutex.unlock();

	// Compile the expressions
	for(int s=0; s< 2*Sp; s++) 
	{	
		std::lock_guard<std::mutex> lock(errorMutex);

		std::string expression_string = MFT_Vec_CoexExpr[s];
		expressions[s].register_symbol_table(local_symbol_table);
		if (!parser.compile(expression_string, expressions[s]))
		{
			
			std::cerr << "Error: Unable to compile expression for species " << s << " for [a, thr, expr]: " 
			<< R << " , " << thr << " , " << expression_string << "\n";
		}

		std::cerr << "Value of expression for species " << s << " for [a, thr, expr]: " 
			<< R << " , " << thr << " , " << expressions[s].value() << "\n";
	}

	//display_symbol_table(local_symbol_table);

	if (R < R_c)
	{
		double highval_veg = std::max(dP, c_spread[0]*(expressions[0].value()));
		/** Use the max function to calculate the maximum of the two values, dP and c_spread[Sp]*expressions[0].value()
		// and assign it to the first species. */
		for(size_t i=0; i< size; i++)
		{
			double r = unif(rng);
			if(r <= perc)
				array[0][i] = highval_veg; // perc fraction of the vegetation species are set to dP
			else
				array[0][i] = c_spread[0]*expressions[0].value(); // Set the other sites of the first species to the MFT value.
			
			for(size_t s=1; s< Sp; s++)
			{
				// If the MFT condition before the critical point is 0, set the species to c_spread[s]
				// Otherwise, set the species to the MFT value multiplied by the constant values in c_spread
				if(expressions[s].value() == 0)
					array[s][i] = c_spread[s];
				else
					array[s][i] = c_spread[s]*expressions[s].value();
			}
		}
	}
	else if ( R >= R_c)
	{
		double highval_veg = std::max(dP, c_spread[Sp]*(expressions[Sp].value()));
		/** Use the max function to calculate the maximum of the two values, dP and c_spread[Sp]*expressions[0].value()
		// and assign it to the first species. */
		for(size_t i=0; i< size; i++)
		{
			double r = unif(rng);
			if(r <= perc)
				array[0][i] = highval_veg;
			else
				array[0][i] = expressions[Sp].value(); // Set the other sites of the first species to the MFT value.
			
			// Set the other species to the constant values
			for(size_t s=1; s< Sp; s++)
			{
				if(expressions[s + Sp].value() != 0)
					array[s][i] = c_spread[s + Sp]*expressions[s+Sp].value();
				else
					array[s][i] = c_spread[s + Sp];
			}	
			// Set the other species to their MFT values multiplied by the constant values in c_spread
		}
	}

	
	local_symbol_table.remove_variable("a");
	errorMutex.lock();
	// Unlock the mutex
	std::cerr << "Expressions filled for [R, thr]: " << R << " , " << thr << "\n";
	errorMutex.unlock();
}

void init_exprtk_randbiMFTframe_OLD(D2Vec_Double &array, int size, double R, double R_c, double dP, double perc,  double c_spread[])
{

	/** Summary of init_exprtk_randbiMFTframe() function.
	 * @brief Initialises a frame (D2Vec_Double vector) as follows:
	 * If R < R_c, then the first Sp elements of c_spread are used to initialise the frame.
	 * In this case, the first row of the frame (corresponding to the first species) has a random distribution of values.
	 * (1 -perc)% of the values are set to 0 (c_spread[0]), and the remaining perc% are set to dP.
	 * The remaining rows of the frame (corresponding to the other species) are set to the constant values given by c_spread + 1 to c_spread + Sp-1.
	 * If R >= R_c, then the last Sp elements of c_spread are used to initialise the frame.
	 * In this case, the first row of the frame (corresponding to the first species) has a random distribution of values.
	 * (1 -perc)% of the values are set to MFT_Vec_CoexExpr[0], and the remaining perc% are set to c_spread[Sp]*MFT_Vec_CoexExpr[0] or dP, whichever is greater.
	 * The other "i" rows of the frame are set to the constant values given by MFT_Vec_CoexExpr[i]*c_spread[i+Sp]  where i ranges from 1 to Sp-1.
	 */

	mutex errorMutex; // Mutex to make error messages thread-safe
	// Create expression objects for each species
    std::vector<exprtk::expression<double>> expressions(Sp);
    exprtk::parser<double> parser;

	// Create a local symbol table, copying the global symbol table
	exprtk::symbol_table<double> local_symbol_table;// = global_symbol_table;

	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);

	uniform_real_distribution<double> unif;
	unif = uniform_real_distribution<double>(0.0, 1.0);
	int thr = omp_get_thread_num();
	//errorMutex.lock();
	local_symbol_table.add_variable("a", R);
	//errorMutex.unlock();

	// Compile the expressions
	for(int s=0; s< Sp; s++) 
	{	
		std::lock_guard<std::mutex> lock(errorMutex);

		std::string expression_string = MFT_Vec_CoexExpr[s];
		expressions[s].register_symbol_table(local_symbol_table);
		if (!parser.compile(expression_string, expressions[s]))
		{
			
			std::cerr << "Error: Unable to compile expression for species " << s << " for [a, thr, expr]: " 
			<< R << " , " << thr << " , " << expression_string << "\n";
		}

		//std::cerr << "Value of expression for species " << s << " for [a, thr, expr]: " 
		//	<< R << " , " << thr << " , " << expressions[s].value() << "\n";
	}

	//display_symbol_table(local_symbol_table);

	if (R < R_c)
	{
		
		for(size_t i=0; i< size; i++)
		{
			double r = unif(rng);
			if(r <= perc)
				array[0][i] = dP; // perc fraction of the vegetation species are set to dP
			else
				array[0][i] = c_spread[0]; // Set the other species to the constant values
			
			for(size_t s=1; s< SpB; s++) 
				array[s][i] = c_spread[s]; // Set the other species to the constant values
			for(size_t s=SpB; s< Sp; s++) 
			{
				// Set the non-biota species to their MFT values multiplied by the constant values in c_spread
				array[s][i] = c_spread[s]*expressions[s].value();
			}
		}
	}
	else if ( R >= R_c)
	{
		double highval_veg = std::max(dP, c_spread[Sp]*(expressions[0].value()));
		/** Use the max function to calculate the maximum of the two values, dP and c_spread[Sp]*expressions[0].value()
		// and assign it to the first species. */
		for(size_t i=0; i< size; i++)
		{
			double r = unif(rng);
			if(r <= perc)
				array[0][i] = highval_veg;
			
			else
				array[0][i] = expressions[0].value(); // Set the other sites of the first species to the MFT value.
			
			// Set the other species to the constant values
			for(size_t s=1; s< Sp; s++) 
				array[s][i] = c_spread[s + Sp]*expressions[s].value();
			// Set the other species to their MFT values multiplied by the constant values in c_spread
		}
	}

	//errorMutex.lock();
	local_symbol_table.remove_variable("a");
	// Unlock the mutex
	//errorMutex.unlock();
}


void init_exprtk_readCSVcolumns_frame(D2Vec_Double &array, vector<int> &const_index, const string &parendir, const string &filenamePattern, 
					const std::vector<std::string> &read_cols, double a, double a_c,  int L, int r, double c_spread[])
{	
	/** SUMMARY: 
	 * @brief Initialises a frame (D2Vec_Double vector called "array") as follows:
	 * The frame is initialised by reading the columns matching entries in read_cols in the CSV files in the directory "parendir" 
	 * with the filename pattern "filenamePattern". The columns are read from the file corresponding to the value of "r", which is
	 * the replicate number. Maximum replicate number is determined by findMaxRepNo(parendir, filenamePattern), so that r values are within bounds.
	 * Column values that do not match the entries in read_cols are ignored, and are generated using expertk expressions (present in MFT_Vec_CoexExpr).
	 * 
	 * The values of the columns "W(x;t)" and "O(x;t)" are stored in the last two rows of the array.
	 * For all other columns, the values are stored in the rows corresponding to the column index - 2.
	*/
    std::mutex errorMutex;

	int maxRepNo = findMaxRepNo(parendir, filenamePattern);

	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);

	uniform_int_distribution<int> unif_int;
	unif_int = uniform_int_distribution<int>(0,  maxRepNo);

	r = r + unif_int(rng); // Randomly select a replicate number
	r = r % (maxRepNo + 1); // Ensure that r is within bounds
	// File name template: "/FRAME_T_{Tval}_a_{a_c}_R_{Rval}.csv"
	string basefilename = filenamePattern + std::to_string(r) + ".csv";
    
    // Loop through all files in the parent directory matching the pattern.
    for (const auto &entry : fs::directory_iterator(parendir)) 
	{
        std::string filename = "/" + entry.path().filename().string();

		

        // Check if filename matches the pattern.
        if (filename.find(basefilename) == std::string::npos) 
		{
            continue; // Skip non-matching files.
        }
		else
		{
			std::lock_guard<std::mutex> lock(errorMutex);
			stringstream m0;
			m0 << "Filename: " << filename << " matches pattern: " << basefilename << endl;
			cerr << m0.str(); cout << m0.str();
		}

        std::ifstream inputFile(entry.path());
        if (!inputFile.is_open()) {
            std::lock_guard<std::mutex> lock(errorMutex);
            std::cerr << "Error opening file: " << filename << std::endl;
            continue;
        }
		

        std::string line;
        std::getline(inputFile, line); // Read header

        // Parse the header line to map column names to indices
        std::unordered_map<std::string, int> header_map;
        std::stringstream ss(line);
        std::string col_name;
        int col_idx = 0;
        while (std::getline(ss, col_name, ',')) {
            //col_name.erase(std::remove_if(col_name.begin(), col_name.end(), isspace), col_name.end());
			col_name.erase(std::remove_if(col_name.begin(), col_name.end(), [](unsigned char c) { return std::isspace(c); }), col_name.end());
            header_map[col_name] = col_idx++;
        }

		/**  Print the header map for debugging
		stringstream m1;
		m1 << "Header map for file: " << filename;
		for (const auto &pair : header_map)
			m1 << "\n" << pair.first << " : " << pair.second;
		cerr << m1.str(); cout << m1.str();
		*/

        // Verify that all columns in read_cols are present in the CSV header
        for (const auto &col : read_cols) 
		{
            if (header_map.find(col) == header_map.end()) 
			{
				stringstream m0;
				m0 << "Column missing in CSV file: " << col << " in file " << filename 
				<< "with [a, thr, r]: " << a  << " , " << omp_get_thread_num() << " , " << r << std::endl;
                inputFile.close();
				cerr << m0.str(); cout << m0.str();
                throw std::invalid_argument("Column missing in CSV file");
				exit(5);
            }
        }

        // Read the data and store the selected columns in the array
        int row = 0;
        while (std::getline(inputFile, line) && row < L * L) 
		{
            std::stringstream line_ss(line);
            std::string cell;
            int csv_col = 0;
            while (std::getline(line_ss, cell, ',')) 
			{
                // Find corresponding column in read_cols
                for (size_t i = 0; i < read_cols.size(); ++i) 
				{
                    if (header_map[read_cols[i]] == csv_col) 
					{
                        try 
						{
							//Check if read_cols[i] == "W(x;t)" or "O(x;t)"
							if (read_cols[i] == "W(x;t)")
							  array[Sp-2][row] = std::stod(cell); // Store value in array
							else if (read_cols[i] == "O(x;t)")
							  array[Sp-1][row] = std::stod(cell); // Store value in array
							else
                            	array[csv_col-2][row] = std::stod(cell); // Store value in array
                        }
						catch (const std::invalid_argument& e) 
                		{
							std::lock_guard<std::mutex> lock(errorMutex);
                    		std::cerr << "Error: Unable to convert cell value to double. Column index: " << i
							<< " , for file: " << filename << " on thread #" << omp_get_thread_num() 
							<< "\n with the offending line: \n" << line_ss.str()<< std::endl;
							//Check if read_cols[i] == "W(x;t)" or "O(x;t)"
							if (read_cols[i] == "W(x;t)")
							  array[Sp-2][row] = 0; // Store value in array
							else if (read_cols[i] == "O(x;t)")
							  array[Sp-1][row] = 0; // Store value in array
							else
                            	array[csv_col-2][row] = 0; // Store value in array
                		} 
                		catch (const std::out_of_range& e) 
                		{
							std::lock_guard<std::mutex> lock(errorMutex);
                    		std::cerr << "Error: Cell value out of range for double. Column index: " << i 
							<< " , for file: " << filename << " on thread #" << omp_get_thread_num() 
							<< "\n with the offending line: \n" << line_ss.str() << std::endl;
							if (read_cols[i] == "W(x;t)")
							  array[Sp-2][row] = 0; // Store value in array
							else if (read_cols[i] == "O(x;t)")
							  array[Sp-1][row] = 0; // Store value in array
							else
                            	array[csv_col-2][row] = 0; // Store value in array
                		}
						catch (const std::exception &e) 
						{
                            std::lock_guard<std::mutex> lock(errorMutex);
                            std::cerr << "Error: Unable to convert cell value to double. Column index: " << i
							<< " , for file: " << filename << " on thread #" << omp_get_thread_num() 
							<< "\n with the offending line: \n" << line_ss.str() << std::endl;

                        }
                    }
                }
                csv_col++;
            }
            row++;
        } // End of while loop for reading the CSV file

		// Check if sizeof const_ind_val is not 0, if so, use expertk to set the values in the array.
		if(const_index.size() != 0)
		{
			
			// Create expression objects for each species
			std::vector<exprtk::expression<double>> expressions(2*Sp);
    		exprtk::parser<double> parser;

			// Create a local symbol table, copying the global symbol table
			exprtk::symbol_table<double> local_symbol_table;// = global_symbol_table;

			uniform_real_distribution<double> unif;
			unif = uniform_real_distribution<double>(0.0, 1.0);
			int thr = omp_get_thread_num();

			local_symbol_table.add_variable("a", a);
			local_symbol_table.add_variable("a_c", a_c);

			// Compile the expressions
			for(int s=0; s< 2*Sp; s++)
			{
				std::lock_guard<std::mutex> lock(errorMutex);

				std::string expression_string = MFT_Vec_CoexExpr[s];
				expressions[s].register_symbol_table(local_symbol_table);
				if (!parser.compile(expression_string, expressions[s]))
				{
					std::cerr << "Error: Unable to compile expression for species " << s << " for [a, thr, expr]: "
					<< a << " , " << thr << " , " << expression_string << "\n";
				}
			}

			// Now iterate over the const_ind_val vector and set the values in the array.
			for(auto& s : const_index)
			{
				if (a < a_c)
				{
					for(size_t i=0; i< L*L; i++)
					{
						// If the MFT condition before the critical point is 0, set the species to c_spread[s]
						// Otherwise, set the species to the MFT value multiplied by the constant values in c_spread
						if(expressions[s].value() == 0)
							array[s][i] = c_spread[s];
						else
							array[s][i] = c_spread[s]*expressions[s].value();
					}
				}
				else if ( a >= a_c)
				{
					// Set the other species to the constant values
					for(size_t i=0; i< L*L; i++)
					{
						if(expressions[s + Sp].value() != 0)
							array[s][i] = c_spread[s + Sp]*expressions[s+Sp].value();
						else
							array[s][i] = c_spread[s + Sp];
					}	
				}
			} // End of for loop for const_index
        
    	} // End of if(const_index.size() != 0)

		
		inputFile.close();
	} // End of for loop for directory_iterator
}

void init_burnin_wrapper(D2Vec_Double &Rho_dt, double a, double a_c, double dP, double perc,  int L, int r, double c_spread[] )
{
	stringstream  rain, dPO, Lgrid, t_val; rain  << a; Lgrid << L; dPO << dP; t_val << 76000; 
	map<string, string> local_input_keys = input_keys; // Copy input_keys to local_input_keys.
	local_input_keys["T"] = t_val.str(); local_input_keys["a"] = rain.str(); local_input_keys["dP"] = dPO.str();
	local_input_keys["L"] = Lgrid.str();
	//Update the values of T, a, dP and L in local_input_keys.
	string local_input_frame_subfolder = input_frame_subfolder; //Copy input_frame_subfolder to local_input_frame_subfolder.
	// Formatted string for input_frame_subfolder.
	string local_input_folder = format_str(local_input_frame_subfolder, local_input_keys); //Format input_folder_pattern with local_input_keys.

	string initcsv_parendir= input_frame_parenfolder + local_input_folder; //Parent directory for csv files.
	vector<int> const_species; //Vector to store species to be initialised with constant values.
	string csv_filename_pattern; //Filename pattern for csv files.
	vector<string> initcsv_columns_temp; //Columns to be read from csv file.
	// Check if local_input_folder contains "HEXBLADE"
	if(local_input_folder.find("HEXBLADE") != std::string::npos || local_input_folder.find("TRUEHEX") != std::string::npos)
	{
		///** INPUT WHEN READING ARTIFICIAL HEX FILES (Filename patterns: FRAME_T_0_a_0_R_{}.csv) 
		csv_filename_pattern = input_prefix + "0_a_0_R_";
		initcsv_columns_temp = {"P(x;t)"};
		
		for(int s=1; s< Sp; s++)
			const_species.push_back(s); //Species to be initialised with constant values (ALL EXCEPT VEGETATION).
		//*/
	}
	else
	{
		///**INPUT Filename Patterns are of the form: FRAME_T_{}_a_{}_R_{}.csv (WHEN READING RIETKERK FRAMES)
		csv_filename_pattern = input_prefix + t_val.str() + "_a_" + rain.str() + "_R_";
		initcsv_columns_temp = {"P(x;t)", "G(x;t)", "Pr(x;t)", "W(x;t)", "O(x;t)"}; //Columns to be read from csv file.
		//for(int s=1; s< SpB; s++)
		//	const_species.push_back(s); //Species to be initialised with constant values.
		//*/
	}

	//Copy the contents of initcsv_columns_temp to initcsv_columns.
	const vector<string> initcsv_columns = initcsv_columns_temp;

	// Remove initcsv_columns_temp from memory.
	vector<string>().swap(initcsv_columns_temp);

	//Check first if initcsv_parendir exists, if not throw an error and exit.
	if(!std::filesystem::exists(initcsv_parendir))
	{
		stringstream m1_1;
		m1_1 << "RUN-TIME ERROR: Parent Directory " << initcsv_parendir << 
		"\nfor Initialisation CSV Files does not exist for Thread Rank:\t " << omp_get_thread_num() 
		<< "  with a_value:\t" << a << "\t and Replicate:\t" << r << " EXITING." <<endl; 
		cout << m1_1.str(); cerr << m1_1.str();
		exit(5);
	}

	init_exprtk_readCSVcolumns_frame(Rho_dt, const_species, initcsv_parendir, csv_filename_pattern, initcsv_columns, a, a_c, L, r, c_spread);

	vector<double> temp_vec= {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species 's'
	double init_veg_num = occupied_sites_of_vector(temp_vec, L*L); //Finds number of occupied at given t.
	vector<double>().swap(temp_vec); //Flush temp_vec out of memory.

	if(init_veg_num == 0.0)
	{
		stringstream m1_1;
		m1_1 << "RUN-TIME WARNING: ZERO Active BURN-IN veg sites for Thread Rank:\t " << omp_get_thread_num() 
		<< "  with a_value:\t" << a << "\t and Replicate:\t" << r << " CONSIDER RE-RUNNING OVER A DIFFERENT RANGE.\n"; 
		cout << m1_1.str(); 
	}

	stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
	m1_1 << "Initial Conditions: BURN-IN # Active Veg sites:\t"<< init_veg_num << " for Thread Rank:\t " << omp_get_thread_num() 
	<< "  with a_value:\t" << a << "\t and Replicate:\t" << r << "\n"; 
	cout << m1_1.str();
}

void init_csvconstframe(D2Vec_Double &array, D2Vec_Double &const_ind_val, const string &parendir, const std::string& filename, const vector<int> &columns, int size) 
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

double variance_of_vector(vector<double> array,int size){

  double mean = mean_of_vector(array,size);
	float sum = 0.0;

	for (int i =0; i<size; ++i){
		sum += (array[i]-mean)*(array[i]-mean);
	}
  double variance = sum/(double)size;
	return variance;
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

void generic_SPBmean_surv_runs(D2Vec_Double &t_avg_var_rep_N, const D2Vec_Double &X_curr, int size, int spcolnum, int j)
{
	/** X_curr contains time-series data (indexed as (N(t), <rho(t)>x)) from the current surviving replicate (given by r) 
	 * with a length (number of rows) "size", and spcolnum is the columns per species in X_curr (so X_curr will have spcolnum*SpB columns).
	 *  rep_avg_var stores the time(t), running density average and variance of X_curr, as well as N(t) and number of survivng replicates 
	 * for each time-step from previous surviving replicates ( <= r-1) respectively.
	 * In this function, X_curr is used to update the running avg and variance of rep_avg_var to account for current surviving replicate r(t).
	 * Note rep_avg_Var[t][3] stores number of surviving runs at t.
	*/
	for(int s=0; s < SpB; s++)
	{
		for(int i=0; i<size; i++)
		{
			// For each species s, save the values of X_curr in the corresponding columns of t_avg_var_rep_N, iterating using spcolnum.
			for(int k = 0; k< spcolnum; k++)
			{
				t_avg_var_rep_N[i][spcolnum*s + k + 1] = X_curr[i][spcolnum*s + k];
			}
		}
	}
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
	for(int s=0; s < SpB; s++)
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

/** Create a function switch_and_bait() that accepts a vector of any type, along with the following arguments:
// 1. low: Lower bound of the range to be switched (default = 0).
// 2. high: Upper bound of the range to be switched (default = 0).
// 3. insert_type: Type of insertion (None, linspace, lnspace) (default = None).
// 4. num_points: Number of points to be inserted (default = 50).
// 5. erase: Boolean to determine if the original vals in the range corresponding to {low ,high} should be erased (default = false).
// The function accepts a vector of any type, and if erase is true, it erases the elements in the range {low, high} and inserts the new elements based on insert_type.
// If insert_type is None, the function only sorts and returns the vector. 
// If insert_type is linspace, it inserts num_points linearly spaced points between low and high.
// If insert_type is lnspace, it inserts num_points logarithmically spaced points between low and high.
// The function returns the sorted vector. */ 

template<typename T>std::vector<T> switchsort_and_bait(std::vector<T> vals, T low, T high, 
	int num_points /**= 50 */, string insert_type /** = "None" */,  bool erase /**= false*/)
{
	// Function to switch and bait elements in a vector.
	if(erase)
		vals.erase(std::remove_if(vals.begin(), vals.end(), [low, high](T x){return x >= low && x <= high;}), vals.end());
	
	// Next check if low < high, if not, swap the values.
	if(low > high){
		T temp = low; low = high; high = temp;
	}
	if(insert_type == "linspace")
	{
		std::vector<double> lin_spaced = linspace(low, high, num_points);
		vals.insert(vals.end(), lin_spaced.begin(), lin_spaced.end());
	}
	else if(insert_type == "lnspace")
	{
		std::vector<double> log_spaced = lnspace(low, high, num_points);
		vals.insert(vals.end(), log_spaced.begin(), log_spaced.end());
	}
	else
		std::cout << "Invalid insert_type. Sorting and returning original vector. \n";
	
	std::sort(vals.begin(), vals.end());
	vals.erase(unique( vals.begin(), vals.end() ), vals.end()); //Remove duplicates if any.
	return vals;
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

// Recursive directory creation function
void recursive_dir_create(const string& path_to_dir)
{	
	#if defined(__GNUC__) && (__GNUC__ >= 9)
	//Creates a directory and all its parent directories if they do not exist.
	if (!std::filesystem::exists(path_to_dir))
	{
		recursive_dir_create(path_to_dir.substr(0, path_to_dir.find_last_of("/\\")));
		std::filesystem::create_directory(path_to_dir);
	}
	#else
		//Creates a directory and all its parent directories if they do not exist (using system commands for older compilers).
		#if defined(_WIN32) || defined(_WIN64)
			//Windows
			string command = "mkdir " + path_to_dir;
			system(command.c_str());
		#else
			//Linux
			string command = "mkdir -p " + path_to_dir;
			system(command.c_str());
		#endif
	#endif

}

// This function finds the maximum replicate number in the directory, given the filename pattern, so to avoid overwriting files.
int findMaxRepNo(const string& parendir, const string& filenamePattern)
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
//void set_density_dependentmortality(double m)
//{	Mm = m;	}

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

void generateNeighboringSitesFromCentral(const auto range_CentralNeighboringSites, 
										 vector<int>& neighboringSites, int i, int j, int L)
{
	// range_CentralNeighboringSites is a subrange of pairs of integers representing the central neighboring sites.

	// Resize neighboringSites if necessary
	if (neighboringSites.size() != range_CentralNeighboringSites.size())
		neighboringSites.resize(range_CentralNeighboringSites.size());
	// Generate neighboring sites for the lattice site (i, j) based on the central neighboring sites for the lattice site (0, 0)
	for (size_t k = 0; k < range_CentralNeighboringSites.size(); ++k) {
		const auto& offset = range_CentralNeighboringSites[k];
		int dx = offset.first;
		int dy = offset.second;

		int nx = (i + dx + L) % L;  // Equivalent to (i + dx) mod L with reflective boundary conditions.
		int ny = (j + dy + L) % L;  // Equivalent to (j + dy) mod L with reflective boundary conditions.
		neighboringSites[k] = nx * L + ny;
	}
}



//----------------------------- Regular Rietkerk Integration Machinery --------------------------------------------//


// ALT INTEGRATION MACHINERY WHERE MIXING OF TIMESTEPS DUE TO LAPLACIAN MAY BE AN ISSUE.



void f_2D(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, double alpha, double d, double rW, double W0,  
	double (&Dxd2)[Sp], double (&K)[3],  double t, double dt, int g)
{
	//Vector function that updates an array containing (dP/dt, dW/dt, DO/dt) for each site in the lattice.
    //Based on the Rietkerk model for plant vegetation dynamics.
	#pragma omp simd // SIMD optimization for better performance.	
	for(int i=0; i < g*g; i++)
	{
        //Equations for the density of plants, soil water and surface water at each site.
        //Note that the Laplacian is calculated using reflective boundary conditions.
        //The Laplacian is calculated using the 5-point stencil method.
        //Equivalent to dP/dt = c*g_max*P*W/(W+K1) - d*P + D*(Laplacian of P)
        f[0][i] = c*gmax*Rho_M[1][i]*Rho_M[0][i]/(Rho_M[1][i] +K[1]) -d*Rho_M[0][i] + (Dxd2[0])*(Rho_M[0][nR2[i][0][1]]  + Rho_M[0][nR2[i][2][1]]  
        + Rho_M[0][nR2[i][1][0]]  + Rho_M[0][nR2[i][1][2]]  - 4*Rho_M[0][i]);
        //Equivalent to dW/dt = alpha*(P+K2*W)/(P+K2) - rW*W + D*(Laplacian of W)
        f[1][i] = alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[2][i] -gmax*Rho_M[1][i]*Rho_M[0][i]/(Rho_M[1][i] +K[1]) - rW*Rho_M[1][i] 
        + (Dxd2[1])*(Rho_M[1][nR2[i][0][1]]  + Rho_M[1][nR2[i][2][1]]  + Rho_M[1][nR2[i][1][0]]  + Rho_M[1][nR2[i][1][2]]  - 4*Rho_M[1][i]);
        //Equivalent to dO/dt = a - alpha*(P+K*W)/(P+K)*O + D*(Laplacian of O)
		f[2][i] = a - alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[2][i] + (Dxd2[2])*(Rho_M[2][nR2[i][0][1]]  + Rho_M[2][nR2[i][2][1]] 
        + Rho_M[2][nR2[i][1][0]]  + Rho_M[2][nR2[i][1][2]]  - 4*Rho_M[2][i]);
		
	}
	
	

}


//Modify the function below to perform RK4 integration of the 2D Vector Rho_t, with a laplacian and time step of dt and spatial step dx.
//The function should also take in the parameters of the model, and the time t at which the integration is being performed.
//The function should also take in the number of grid points g, and the number of species Sp.
//The function should also take in the number of threads being used, and the thread number.
//The function should also take in the number of iterations performed so far, and the number of iterations to be performed.


void RK4_Integrate(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D2Vec_Double &K1, D2Vec_Double &K2,  D2Vec_Double &K3, D2Vec_Double &K4, D3Vec_Int &nR2, 
	double a, double c, double gmax, double alpha, double d, double rW, double W0, double (&Dxd2)[Sp], double (&K)[3],  double t, double dt, double dx, int g)
{
    

    // Rho_t stores value of rho at time t.

	// Rho_tstar stores value of rho*, provided by sampling Gaussian.
	//std::ofstream errout; //Save Error Logs
	//std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	double dt6 = dt/6.0; double dt2 = dt/2.0; //double dx2 = dx*dx;

    
    //D2Vec_Double Rho_M(Sp, vector<double> (g*g, 0.0));

	//double K1[Sp][g*g] ={0.0}; double K2[Sp][g*g] ={0.0}; double K3[Sp][g*g] ={0.0}; double K4[Sp][g*g] ={0.0};
	//double Rho_M[Sp] ={0.0, 0.0};


	f_2D(K1, Rho_t, nR2, a, c, gmax, alpha, d, rW, W0, Dxd2, K, t, dt, g); //K1 updated.


	for(int s= 0; s <Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K1[s][i];
	}

	f_2D(K2, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, Dxd2, K, t + dt2, dt,  g); //K2 updated.

	for(int s= 0; s <Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K2[s][i];
	}
	f_2D(K3, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, Dxd2, K, t + dt2, dt,  g); //K3 updated.
	for(int s= 0; s <Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt)*K3[s][i];
	}
	f_2D(K4, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, Dxd2, K, t + dt, dt, g); //K4 updated.
    
	for(int s= 0; s <Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for(int i=0; i < g*g; i++)
		{
			Rho_t[s][i]+= (dt6)*( K1[s][i] + 2.0*K2[s][i] + 2.0*K3[s][i] + K4[s][i]);
		}
		#ifdef DEBUG
		bool isKoed = false; // Flag to indicate if the RK4 integration was "KO'ed" due to invalid values.
		for(int i=0; i < g*g; i++)
		{
			if( Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || isnan(Rho_t[s][i]) == true)
			{
				if(isKoed == false)
				{
					stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
					m6 << "RK4 WAS KO'ED WITH:\t" << Rho_t[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
					<< " at time:\t:" << t << " For Species:\t:" << s << " with K1[s][i]:   " << K1[s][i] << "\t, K2[s][i]:\t" << K2[s][i] 
					<< "\t, K3[s][i]:\t" << K3[s][i] << "\t AND K4[s][i]:\t" << K4[s][i] << "\n"; cerr << m6.str();
						//cout << m6.str();
					//errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					isKoed = true;
				}
			}
		}
		#endif //DEBUG
			
	}
}

void RK4_Wrapper_2D(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, 
	  double alpha, double d, double rW, double W0, double (&D)[Sp], double (&K)[3], double a_st, double a_end, double dt, double dx, double dP, int r, int g)
{
	// Rho_0 stores value of rho at time t=0.
	// Rho_t stores value of rho at time t.
	// Rho_tstar stores value of rho*, provided by sampling Gaussian.
	// Rho_M stores value of rho at intermediate time steps.

	//D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0));
	//double Rho_tsar[Sp][g*g] ={0.0};
	//double Rho_M[Sp][g*g] ={0.0, 0.0};

	int thrID; //Unique ID for each thread, used for saving files.
	
	D3Vec_Int nR2 (g*g, D2Vec_Int (3, vector<int> (3, 0)));
    //nR2 is 3D [(L*L)x3x3] vector initialised to 0, which stores neighbours of all sites i (second dim) in a ball of NORM 1 (SQUARE) radius 1 (2ND & 3RD dim).
	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_Sq8(g, nR2); //Assigning neigbours.
	//int tot_iter = static_cast<int>(1.0/dt) + 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	int tot_iter = t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	//t_meas[0] = 100; t_meas[1] = 200; t_meas[2] = 300; 

	//double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0};
	D2Vec_Double rho_rep_avg_var(tot_iter, vector<double> (Sp4_1, 0.0));

	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep
	//int ind=0;
	//for(double t =0; t <= 1; t+=dt)
	//{ rho_rep_avg_var[ind][0] = t; ind+=1; }
	//for(int i =0; i< t_meas.size(); i++)
	//{ rho_rep_avg_var[ind+i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep

	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t a:\t" << a << "\t c:\t " << c << "\t dt,dx:\t " << dt << "," << dx << "\t gmax:\t " << gmax 
	<< "\t and  Alpha:\t " << alpha  << "\n"; //cout << m1.str();
	cout << m0.str();
	//std::ofstream errout; //Save Error Logs
	cerr << m0.str(); // Output to standard error stream for debugging purposes.
	//std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	//errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	double beta = c*gmax - d; double dx2 = dx*dx; // dx2 is the square of the spatial step.
	double diff_coefficient[Sp]={0.0};
	for(int s= 0; s < Sp; s++)
		diff_coefficient[s] = D[s]/(dx2);
  	double p0istar, p0jstar, p0mstar; // Analytic steady state values.

	if(a <= 1)
	{
		p0istar = 0.0; p0jstar = a/rW; p0mstar = a/(alpha*W0);
	}
	else
	{
		p0jstar = d*K[1]/beta; 
 		p0istar = (c/d)*(a - rW*p0jstar);
  		p0mstar = (a/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
	}

	//double def_sigma[Sp] = {0, 0, 0}; // Non-stochastic default.
	//double v[SpB] ={0}; //Velocity of species.
	//double A[SpB][SpB]; double H[SpB][SpB]; 
  	//double E[SpB] ={1.0}; double M[SpB] ={d}; double pR[SpB] ={0.0};

	auto start_t = high_resolution_clock::now(); //Start the clock
	for(int j=0; j< r; j++)
	{

		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		uniform_real_distribution<double> unif(0.0, 1.0); normal_distribution<double> norm(0.0, 1.0);
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

		//Define thrID as a random number between 0 and 10000, iff j =0.
		thrID = (j == 0) ? int(round(unif(rng)*10000)) : thrID; //Random number between 0 and 1000.


		D2Vec_Double Rho_dt(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		//double Rho_M[tot_iter][Sp2] ={0.0}; //Stores <rho>x values per time step for a single replicate.

		D2Vec_Double gamma(SpB, vector<double> (g*g, 0.0)); //Stores gamma for each species at each site.
		// 2D vector of pairs of dim: SpBx(L*l), called v. Stores the velocity of each species at each site.
		//D2Vec_Pair_Double v_eff(SpB, std::vector<std::pair<double, double>>(g*g, {0.0, 0.0}));

		D2Vec_Double K1(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K2(Sp, vector<double> (g*g, 0.0));
    	D2Vec_Double K3(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K4(Sp, vector<double> (g*g, 0.0));
		D2Vec_Double Rho_M(tot_iter, vector<double> (Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		double perc = 0.015; double c_high[Sp] ={p0istar + dP, p0jstar, p0mstar}; double c_low[Sp] ={p0istar, p0jstar, p0mstar};
		
		init_randconstframe(Rho_dt, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
		

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1 << "Initial Conditions:\t Per:\t" << perc  << "\t C_High[0], C_High[1], C_High[2]:\t " << c_high[0] << "," << c_high[1] << "," << c_high[2]   
		<< "\t C_Lo[0], C_Lo[1], C_Lo[2]:\t " << c_low[0] << "," << c_low[1] << "," << c_low[2] << ",\t R: " << a << endl; //cout << m1.str();
		cout << m1.str(); cerr << m1.str(); // Output to standard error stream for debugging purposes.
		//std::ofstream errout; //Save Error Logs
		//std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
		//errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();

		for(int i=0; i < g*g; i++)
		{
			for(int s= 0; s <Sp; s++)
			{
				Rho_tsar[s][i] = Rho_dt[s][i];
			}	
		}

		double t=0; int index = 0;  //Initialise t
		int iter_index=0; int po=1; int so =1; int lo=1;
		int counter=0;

		while (t < t_max + dt)
		{
			if(t == 0.0 || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
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
				//iter_index+=1;//Rho_M.push_back({t, rhox_avg});

				// Saving frames to file at given time points.
				if( index == int(tot_iter*0.6) || index == int(tot_iter*0.85) || index == int(tot_iter*0.9) || index == int(tot_iter*0.95) || index == int(tot_iter-1))
				{

					D2Vec_Double rho_rep_avg_var_temp(index+1, vector<double> (Sp4_1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
					//std::copy(rho_rep_avg_var.begin(), rho_rep_avg_var.begin() + index, rho_rep_avg_var_temp.begin());
					var_mean_incremental_surv_runs(rho_rep_avg_var_temp, Rho_M, index+1, 0);

					//Time is stored in the first column.
					for(int i=0; i < index+1; i++)
						rho_rep_avg_var_temp[i][0] = t_meas[i];

					stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, Dm, gm, jID;

					L << g; tm << t_meas[index]; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dimitri << dP;
					rini << j; Dm << setprecision(3) << D[0]; gm << gmax; jID << thrID; a1 << a_st; a2  << a_end;
					// Three replicates are over.
					
					string parendir ="";
					parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/TimeSeries"; // Parent directory for timeseries files.
					string filenamePattern = replicate_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
					"_D0_"+ Dm.str() + "_jID_"+ jID.str() + "_R_";

					save_prelimframe(rho_rep_avg_var_temp, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);
				}

				if( t >= 1000 && t <= 2000 || t== 0 ||  index >= tot_iter -15 &&  index <= tot_iter)
				{
					//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.
					
					stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, alph, w0t, dix, dimitri,jID;

  					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  					rini << j; Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; gm << setprecision(3) << gmax;
					w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; dimitri << setprecision(3) << dP;
					a1 << a_st; a2  << a_end; jID << thrID;
					// Three replicates are over.

					//string parendir = "../Data/Rietkerk/Frames/Non-Stochastic/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();
					string parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();
					string filenamePattern = "/FRAME_RAND_Basic_P_c_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() 
					+ "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() + "_jID_"+ jID.str() + "_dx_"+ dix.str() + "_alph_"+ alph.str() + "_gmax_"+ gm.str() + "_R_";
					string filename = filenamePattern + rini.str() + ".csv"; // Naive filename
					save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, "a_c,  x,  P(x; tmax), W(x; tmax), O(x; tmax) \n");
					/** OLD WAY OF SAVING FRAME FILES 
					 ofstream frame_dp;
  					// Creating a file instance called output to store output data as CSV
					frame_dp.open("../Data/Rietkerk/Frames/Non-Stochastic/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_Basic_P_c_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() 
					+ "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() + "_dx_"+ dix.str() + "_alph_"+ alph.str() + "_gmax_"+ gm.str() + "_R_"+ rini.str() + ".csv");

					// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|
					frame_dp << "a_c,  x,  P(x; tmax), W(x; tmax), O(x; tmax) \n";
					for(int i=0; i< g*g; i++)
					{
						frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] <<"\n";
				
					}

					frame_dp.close();
					*/
					
				}
				
				index+=1;

			}

			for(int s= 0; s <Sp; s++)
				Rho_tsar[s] = Rho_dt[s]; //Assigning initial conditions

			RK4_Integrate(Rho_dt, Rho_tsar, K1, K2, K3, K4, nR2, a, c, gmax, alpha, d, rW, W0, diff_coefficient, K, t, dt, dx, g);

			for(int s= 0; s <Sp; s++)
			{
				for(int i=0; i<g*g; i++)
				{
					if( Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true and so==1)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        				m6 << "RHO FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << Rho_tsar[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<<  "\n"; cout << m6.str(); cerr << m6.str(); //so=0;
						exit(1);
						//rrout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}

				}
				
			}
			

			t += dt;
			counter+=1; //Update counter

			if(counter%40000 ==1)
			{
				auto end_t = high_resolution_clock::now();
				auto elapsed_min = duration_cast<minutes>(end_t-start_t);
				// Reset the clock
				start_t = high_resolution_clock::now();
				stringstream m7;
				m7 << "STATUS UPDATE AT [t, a, thr, j]\t " << t << " , "  << a << " , " << omp_get_thread_num() << " , " << j 
				//<< " ,with current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) 
				<< " and time taken for " << counter << " iterations: " << elapsed_min.count() << " min." <<  "\n"; cout << m7.str();
				cerr << m7.str(); // Output to standard error stream for debugging purposes.
				//errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();

				//Report syst every 50000 iterations.
			

			}

			

		} //End of Time Loop


		if(int(omp_get_thread_num()) == 5)
		{
			size_t currentSize = getCurrentRSS( ); //Check Memory usage.
			size_t peakSize    = getPeakRSS( );
			stringstream m7;
			m7 << "In replicate " << j << ", current size in MB: " << currentSize/(1024.0*1024.0) 
			<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << "\n"; cout << m7.str();
			cerr << m7.str(); // Output to standard error stream for debugging purposes.
			//errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();
		}

		if(j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "INITIAL UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << tot_iter << "\n"; cout << m7.str(); cerr << m7.str();
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
  			rini << j; Dm << setprecision(3) << D[0]; gm << gmax; a1 << a_st; a2  << a_end;
			// Three replicates are over.
			
			string parendir ="";
			parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/"; // Parent directory for prelim files.
			string filenamePattern = prelim_prefix_basic + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
			"_D0_"+ Dm.str() + "_gmax_"+ gm.str() + "_R_";

			save_prelimframe(rho_rep_avg_var, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, true);
			/** OLD WAY OF SAVING PRELIM FILES 
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
				output_dp <<  "\n";
				
			}

			output_dp.close();
			//*/
		}




	} //End of Replicate Loop


	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	cerr << m9.str(); // Output to standard error stream for debugging purposes.
	//errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	
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
	cerr << m10.str(); // Output to standard error stream for debugging purposes.
	//errout.open(thr, std::ios_base::app); errout << m10.str(); errout.close();
}



void first_order_critical_exp_delta(int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double d, double rW, double W0,  double (&D)[Sp], double (&K)[3], double dt, double dx, double dP, int r,  int g)
{

	
	//double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2
	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//double mean[Sp] = {1.0, 0.05}; double sd[Sp] = {0.25, 0.0125};
	//init_randframe(Rho_0, Sp,  g*g); //Returns Rho_0 with a full initial frame filled with 0.2.
	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "NOSTRA" <<"\n";
  	// The pspace to iterate over.
    vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	t_measure.insert(t_measure.begin(), 0.0);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  	cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<"\n";
  	for (int i=0; i< t_measure.size(); i++)
  	{
		cout << t_measure[i] << " ";
  	} 	cout << "\n";

	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd); // Seed the generator
	std::uniform_int_distribution<int> unif_int(0, 15); //Use this to generate random numbers between 0 and 15.
	
  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = unif_int(rng)* microToSeconds;     //5 seconds
	// Random delay to stagger the start of the parallel threads.
	usleep(static_cast<int>(delay1)); //Sleep for delay1 microseconds
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
	double alpha, double rW, double W0, double (&Dxd2)[Sp], double (&K)[3], double t, double dt, double g)
{
	//Vector function that updates an array containing (dP/dt, dW/dt, DO/dt) for each site in the lattice.
    //Based on the Rietkerk model for plant vegetation dynamics with the Dornic twist where linear and stoch term for vegetation are already taken care of.
	// The following values are precomputed to reduce time complexity.
	double cgmax = c*gmax; double K2W0 = K[2]*W0;
	for(int i=0; i < g*g; i++)
	{
        //Equations for the density of plants, soil water and surface water at each site.
        //Note that the Laplacian is calculated using reflective boundary conditions.
        //The Laplacian is calculated using the 5-point stencil method.
        /**
		 * NOTE!!!!!:  Equivalent to dP/dt = c*g_max*P*W/(W+K1)
		 *  Linear and stochastic terms taken care of by Dornic integration routine previously.
		**/
		
        f[0][i] = cgmax*Rho_M[1][i]*Rho_M[0][i]/(Rho_M[1][i] +K[1]);
        //Equivalent to dW/dt = alpha*(P+K2*W)/(P+K2) - rW*W + D*(Laplacian of W)
        f[1][i] = alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[2][i] -gmax*Rho_M[1][i]*Rho_M[0][i]/(Rho_M[1][i] +K[1]) - rW*Rho_M[1][i] 
        + (Dxd2[1])*(Rho_M[1][nR2[i][0][0]]  + Rho_M[1][nR2[i][0][1]]  + Rho_M[1][nR2[i][1][0]]  + Rho_M[1][nR2[i][1][1]]  - 4*Rho_M[1][i]);
        //Equivalent to dO/dt = a - alpha*(P+K*W)/(P+K)*O + D*(Laplacian of O)
		f[2][i] = a - alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[2][i] + (Dxd2[2])*(Rho_M[2][nR2[i][0][0]]  + Rho_M[2][nR2[i][0][1]] 
        + Rho_M[2][nR2[i][1][0]]  + Rho_M[2][nR2[i][1][1]]  - 4*Rho_M[2][i]);
		
	}

}


void RK4_Integrate_Stochastic_1Sp(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D2Vec_Double &K1, D2Vec_Double &K2,  D2Vec_Double &K3, D2Vec_Double &K4, 
	D3Vec_Int &nR2, double a, double c, double gmax, double alpha, double rW, double W0, double (&Dxd2)[Sp], double (&K)[3], double t,double dt,double dx, int g)
{
	
	/** // OLD DECLARATIONS of K1, K2, K3, K4 and Rho_M. Skipped as declaration and assignment is done in rietkerk_Dornic_2D_MultiSp
	 * and doing it here takes O(n^2) time each time this function is called.
	*/

	// RECALL: Dxd2[s] = D[s]/dx2
	f_2Dor(K1, Rho_t, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, t, dt,  g); //K1 updated.

	for(int s= 0; s <Sp; s++)
	{
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K1[s][i];
	}

	f_2Dor(K2, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, t + dt2, dt,  g); //K2 updated.

	for(int s= 0; s <Sp; s++)
	{
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K2[s][i];
	}
	f_2Dor(K3, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, t + dt2, dt, g); //K3 updated.
	for(int s= 0; s <Sp; s++)
	{
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt)*K3[s][i];
	}
	f_2Dor(K4, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, t + dt, dt, g); //K4 updated.
    
	for(int s= 0; s <Sp; s++)
	{
		for(int i=0; i < g*g; i++)
		{
			Rho_t[s][i]+= (dt6)*( K1[s][i] + 2.0*K2[s][i] + 2.0*K3[s][i] + K4[s][i]);

			if( Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || isnan(Rho_t[s][i]) == true)
			{
				std::ofstream errout; //Save Error Logs
				std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RK4 WAS KO'ED WITH:\t" << Rho_t[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " For Species:\t:" << s << " with K1[s][i]:   " << K1[s][i] << "\t, K2[s][i]:\t" << K2[s][i] 
				<< "\t, K3[s][i]:\t" << K3[s][i] << "\t, K4[s][i]:\t" << K4[s][i] << "\t AND Rho(t-dt)" << Rho_tsar[s][i] << endl; //cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
		}
	}

}


//------------------- Vegetation + Grazer (+ Soil Water + Surface Water) -------------------//

void f_2Dor_2Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, 
	double alpha, double rW, double W0, double (&Dxd2)[Sp], double (&K)[3], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t, double dt, double dx1_2, double g)
{
	//cout << "Placeholder" << endl;
	//Vector function that updates an array containing ( dP/dt, dW/dt, DO/dt,  dG/dt) for each site in the lattice.
    //Based on the Rietkerk model for plant vegetation dynamics with the Dornic twist where linear and stoch term for vegetation (and grazer) are already taken care of.
	// The following values are precomputed to reduce time complexity.
	double cgmax = c*gmax; double K2W0 = K[2]*W0; double A01H01 = A[0][1]*H[0][1]; double EA1 = E[1]*A[0][1];

	for(int i=0; i < g*g; i++)
	{
        //Equations for the density of plants, soil water, surface water and grazers at each site.
        /** NOTE!!!!!:  Equivalent to dP/dt = c*g_max*P*W/(W+K1) - aij*hij*P*G/(1+aij*hij*V)
		 *  Linear and stochastic terms taken care of by Dornic integration routine previously.
		**/
        f[0][i] = cgmax*Rho_M[2][i]*Rho_M[0][i]/(Rho_M[2][i] +K[1]) -(A[0][1]*Rho_M[0][i]/(1 + A01H01*Rho_M[0][i]))*Rho_M[1][i];
		//Equivalent to dG/dt = e*(aij*V*P)/(1+aij*hij*V)  [NO PREDATION].
		f[1][i] = (EA1*Rho_M[0][i]/(1 + A01H01*Rho_M[0][i]))*Rho_M[1][i];
        //Equivalent to dW/dt = alpha*(P+K2*W0)/(P+K2)*O - rW*W + D*(Laplacian of W)
        f[2][i] = alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[3][i] -gmax*Rho_M[2][i]*Rho_M[0][i]/(Rho_M[2][i] +K[1]) - rW*Rho_M[2][i] 
        + (Dxd2[2])*(Rho_M[2][nR2[i][0][0]]  + Rho_M[2][nR2[i][0][1]]  + Rho_M[2][nR2[i][1][0]]  + Rho_M[2][nR2[i][1][1]]  - 4*Rho_M[2][i]);
        //Equivalent to dO/dt = a - alpha*(P+K*W)/(P+K)*O + D*(Laplacian of O)
		f[3][i] = a - alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[3][i] + (Dxd2[3])*(Rho_M[3][nR2[i][0][0]]  + Rho_M[3][nR2[i][0][1]] 
        + Rho_M[3][nR2[i][1][0]]  + Rho_M[3][nR2[i][1][1]]  - 4*Rho_M[3][i]);
	}
}

void RK4_Integrate_Stochastic_2Sp(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D2Vec_Double &K1, D2Vec_Double &K2, D2Vec_Double &K3, D2Vec_Double &K4, D3Vec_Int &nR2,
		double a,double c,double gmax,double alpha, double rW, double W0, double (&Dxd2)[Sp], double (&K)[3], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t,double dt,double dx, int g)
{
	/** // OLD DECLARATIONS of K1, K2, K3, K4 and Rho_M. Skipped as declaration and assignment is done in rietkerk_Dornic_2D_MultiSp
	 * and doing it here takes O(n^2) time each time this function is called.
	*/
	
	// RECALL: Dxd2[s] = D[s]/dx2
	f_2Dor_2Sp(K1, Rho_t, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A,H,E, t, dt, dx1_2, g); //K1 updated.

	for(int s= 0; s <Sp; s++)
	{
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K1[s][i];
	}

	f_2Dor_2Sp(K2, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A,H,E, t + dt2, dt, dx1_2, g); //K2 updated.

	for(int s= 0; s <Sp; s++)
	{
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K2[s][i];
	}
	f_2Dor_2Sp(K3, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A,H,E, t + dt2, dt, dx1_2, g); //K3 updated.
	for(int s= 0; s <Sp; s++)
	{
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt)*K3[s][i];
	}
	f_2Dor_2Sp(K4, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A,H,E, t + dt, dt, dx1_2, g); //K4 updated.
    
	for(int s= 0; s <Sp; s++)
	{
		bool isKoed = false;
		for(int i=0; i < g*g; i++)
		{
			Rho_t[s][i]+= (dt6)*( K1[s][i] + 2.0*K2[s][i] + 2.0*K3[s][i] + K4[s][i]);

			if( Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || isnan(Rho_t[s][i]) == true)
			{
				if(isKoed == false)
				{
					std::ofstream errout; //Save Error Logs
					std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
					stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
					m6 << "RK4 WAS KO'ED WITH VAL:\t" << Rho_t[s][i] << " For [s, thr, a, i, Pre]\t" << s << " , " << omp_get_thread_num() << " , " << a 
					<< " , " << i << " , " << prefix << " ] with K1[s][i]:   " << K1[s][i] << "\t, K2[s][i]:\t" << K2[s][i] 
					<< "\t, K3[s][i]:\t" << K3[s][i] << "\t, K4[s][i]:\t" << K4[s][i] << "\t AND Rho(t-dt)" << Rho_tsar[s][i] << "\n"; //cout << m6.str();
					errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					cerr << m6.str();
					isKoed = true;
				}
			}
		}
	}
}

// FFTW3-based gamma sweep calculation for 3-Species Non-Refugia model.
void calc_gamma_2Sp_NonRefugia_FFT(D2Vec_Double& Rho_t, D2Vec_Double& gamma, double (&Rho_avg)[Sp], vector<std::pair<double, int>>& rfrac, 
    vector<std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L, FFTW3_CentralPlanner& fft_central_plan)    
{
	mutex errorMutex; // Mutex to make error messages thread-safe
    int L2 = L * L;
    double eps = 1.0e-12;

	int rmax_eff = int(nVeg_frac * r_max); //Effective r_max based on vegetation fraction.
	

    double rho_inverse[SpB];
    for (int s = 0; s < SpB; s++)
        Rho_avg[s] >= eps ? rho_inverse[s] = 1 / Rho_avg[s] : rho_inverse[s] = 1 / eps; //Avoid division by zero.
    
	for (const auto& frac : rfrac) {
        int s = frac.second;
        double fr = frac.first;

		int r_eff = int(fr * rmax_eff); //Effective r based on species fraction and vegetation cover.
		int eff_nR_Perp_size = nR_Perp_R[r_eff]; //Get precomputed effective nR_Perp size.
        //int eff_nR_Perp_size = int(fr * fr * new_length);
        if (fr == 0.0) continue; //Skip if fraction is zero.
		if (s == 0) continue; //Skip vegetation.
        if (eff_nR_Perp_size < 1) eff_nR_Perp_size = 1;
        if (Rho_avg[s] < eps) 
		{
			#pragma omp simd
            for (int i = 0; i < L2; i++) gamma[s][i] = 0.0;
            continue;
        }
        if (Rho_avg[s-1] < eps) 
		{
			#pragma omp simd
            for (int i = 0; i < L2; i++) gamma[s][i] = 0.0;
            continue;
        }
        
        // Prepare spatial field based on species
        if (s == 1) 
		{
            // Copy vegetation density directly
            std::copy(Rho_t[s-1].begin(), Rho_t[s-1].end(), fft_central_plan.rho_spatial);

        }
        
        // Forward FFT of density (Rho) field
        fftw_execute_dft_r2c(fft_central_plan.forward_rho, fft_central_plan.rho_spatial, fft_central_plan.rho_fft);
        
        // Get precomputed kernel FFT
        fftw_complex* kernel_fft = kHat_kernel_FFT_cache[r_eff];
        
		double real_part, imag_part;
        // Multiply in Fourier space (complex multiplication)
		#pragma omp simd
        for (int i = 0; i < L2; i++) 
		{
			// RECALL: Z = (a + ib)*(c + id) = (ac - bd) + i(ad + bc)
            real_part = fft_central_plan.rho_fft[i][0] * kernel_fft[i][0] - 
                fft_central_plan.rho_fft[i][1] * kernel_fft[i][1];
            imag_part = fft_central_plan.rho_fft[i][0] * kernel_fft[i][1] + 
                fft_central_plan.rho_fft[i][1] * kernel_fft[i][0];
            fft_central_plan.prod_fft[i][0] = real_part; 
			fft_central_plan.prod_fft[i][1] = imag_part;
        }
        
        // Inverse FFT to get convolution result
        fftw_execute_dft_c2r(fft_central_plan.backward, fft_central_plan.prod_fft, fft_central_plan.convolution_spatial);
        // Normalize and apply species-specific scaling.
        double normalisation = 1.0 / L2;  // FFTW3 doesn't normalise - multiplication by 1/(L*L) is needed.
		#pragma omp simd
        for (int i = 0; i < L2; i++) 
		{
            gamma[s][i] = fft_central_plan.convolution_spatial[i] * normalisation * rho_inverse[s-1];
		}
		// Clamp to [0, 1]
        for (int i = 0; i < L2; i++) 
		{   
			if (gamma[s][i] > 1.0) gamma[s][i] = 1.0;
            else if (gamma[s][i] < 0.0) gamma[s][i] = 0.0;
			else if (isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
			{
				errorMutex.lock();
				std::cerr << "Gamma is NaN with value: " << gamma[s][i] << " For [s, thr, i]\t" << s << " , " << omp_get_thread_num()
					<< " , " << i << " ] with Rho_avg[0]: " << Rho_avg[0] << " and Rho_avg[1]: " << Rho_avg[1] << " and Rho_avg[2]: " << Rho_avg[2]
					<< " \n and Rho_t[0][i]: " << Rho_t[0][i] << " and Rho_t[1][i]: " << Rho_t[1][i] << " and Rho_t[2][i]: " << Rho_t[2][i]
					<< " \n and Nveg_frac: " << nVeg_frac << " and r_eff: " << r_eff << " and rfrac.size(): " << rfrac.size() << " and nR_Perp_R[r_eff]: " << nR_Perp_R[r_eff]
					<< std::endl;
				errorMutex.unlock();
			}
        }
    }
}

// Calculate the gamma values for 2Species model with Vegetation and Grazers.
void calc_gamma_2Sp_NonRefugia(const vector<pair<int, int>>& centralNeighboringSites, D2Vec_Double &Rho_t, D2Vec_Double &gamma, double (&Rho_avg)[Sp],   
	vector <std::pair<double, int>>& rfrac, vector<std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L)
{	
	if(nVeg_frac < 0.35)
		nVeg_frac = 0.35; //Minimum fraction of vegetation cover.

	int L2 = L*L;
	double eps = 1.0e-12; //Small number to avoid division by zero.

	int new_length = int(nVeg_frac*nVeg_frac*centralNeighboringSites.size());
	new_length < 1 ? new_length = 1 : new_length = new_length; // At least one nearest neighbor (the site itself) is considered.
	// Make a copy of the centralNeighboringSites, resized to the fraction nVeg_frac*nVeg_frac.
	auto range_CentralNeighboringSites = std::ranges::subrange(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);
	//vector<pair<int, int>> centralNeighboringSitesCopy(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);

	std::vector<int> nR_Perp(new_length);

	double rho_inverse[SpB]; //Inverse of average density.
	for(int s=0; s< SpB; s++)
		Rho_avg[s] >= eps ? rho_inverse[s] = 1/Rho_avg[s] : rho_inverse[s] = 1/eps; //Avoid division by zero.
	// Avoid division by zero.

	for( int i=0; i < L2; i++)
	{
		mutex errorMutex; // Mutex to make error messages thread-safe
		
		int c_i = int(i/L); int c_j = i%L; //Current x and y coordinates of site i.

		generateNeighboringSitesFromCentral(range_CentralNeighboringSites, nR_Perp, c_i, c_j, L);
		//generateNeighboringSitesFromCentral(centralNeighboringSites, nR_Perp, c_i, c_j, L);

		//if(nVeg_frac < 1)
		//	nR_Perp.resize(int(nVeg_frac*nVeg_frac*nR_Perp.size()));
		// In this case, r_effective is reduced to nVeg_frac*r_max, hence the reduction in the number of nearest neighbours.
		// rfrac is sorted in descending order.
		double fr_prev = 1.0; //Fraction of perception radius to max radius.
		for(const auto& frac: rfrac)
		{
			int s = frac.second;
			double fr = frac.first; //Fraction of perception radius to max radius.

			int eff_nR_Perp_size = int(fr*fr*nR_Perp.size());
			
			if( fr == 0.0)
			{	continue;	} // No perception radius for this species (indicates vegetation species)
			if(Rho_avg[s] < eps)
				{	gamma[s][i] = 0.0; continue;	} // No species left, value assigned doesn't matter.
			else if(Rho_avg[s-1] < eps)
				{	gamma[s][i] = 0.0;	continue; } // No resource left, consumers will advect to extinction.
			
			// Gamma for grazers (indexed 1)
			for(int k=0; k< eff_nR_Perp_size; k++)
			{	gamma[s][i] += Rho_t[s-1][nR_Perp[k]];	}

			gamma[s][i] = (gamma[s][i]*rho_inverse[s-1])/(eff_nR_Perp_size);

			if(gamma[s][i] > 1)
			{	gamma[s][i] = 1.0;	} //Ensuring gamma is not greater than 1.
			else if(gamma[s][i] < 0)
			{	gamma[s][i] = 0.0;	} //Ensuring gamma is not less than 0.
			else if(isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
			{
				errorMutex.lock();
				std::cerr << "Gamma is NaN with value: " << gamma[s][i] << " For [s, thr, i]\t" << s << " , " << omp_get_thread_num() 
						<< " , " << i << " ] with Rho_avg[0]: " << Rho_avg[0] << " and Rho_avg[1]: " << Rho_avg[1] << " and Rho_avg[2]: " << Rho_avg[2] 
						<< " \n and Rho_t[0][i]: " << Rho_t[0][i] << " and Rho_t[1][i]: " << Rho_t[1][i]
						<< " \n and Nveg_frac: " << nVeg_frac << " and r_max: " << r_max << " and rfrac.size(): " << rfrac.size() << " and nR_Perp.size(): " << nR_Perp.size()
						<< std::endl;
				errorMutex.unlock();
			}
		}
		vector <int>().swap(nR_Perp);
	}
}

//------------------- Vegetation + Grazer + Predator (+ Soil Water + Surface Water) -------------------//

void f_2Dor_3Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, 
	double alpha, double rW, double W0, double (&Dxd2)[Sp], double (&K)[3], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t, double dt, double dx1_2, double g)
{
	//Vector function that updates an array containing ( dP/dt, dW/dt, DO/dt,  dG/dt, dPr/dt) for each site in the lattice.
    //Based on the Rietkerk model for plant vegetation dynamics with the Dornic twist where linear and stoch term for vegetation (and grazer) are already taken care of.
	
	// The following values are precomputed to reduce time complexity.
	double cgmax = c*gmax; double K2W0 = K[2]*W0; double A01H01 = A[0][1]*H[0][1]; double A12H12 = A[1][2]*H[1][2];
	// The following are the values of E*A and A*H for each species, which reduces time complexity by avoiding repeated calculations.
	double EA1 = E[1]*A[0][1]; double EA2 = E[2]*A[1][2]; int L2 = g*g; //L2 is the total number of sites in the lattice.
	//vector<double> EA(SpB);
	//for(int s=1; s< SpB; s++)
	//	EA[s] = E[s]*A[s-1][s];
	#pragma omp simd // SIMD optimization for better performance.	
	for(int i=0; i < L2; i++)
	{
        //Equations for the density of plants, soil water, surface water and grazers at each site.
        //Note that the Laplacian is calculated using reflective boundary conditions.
        //The Laplacian is calculated using the 5-point stencil method.
        /**
		 * NOTE!!!!!:  Equivalent to dP/dt = c*g_max*P*W/(W+K1) - aij*hij*P*G/(1+aij*hij*V)
		 *  Linear and stochastic terms taken care of by Dornic integration routine previously.
		**/

		// RECALL: Dxd2[s] = D[s]/dx2 AND K_P1 = Inverse of Carrying Capacity of Predator (defined as Km in main())
        f[0][i] = cgmax*Rho_M[3][i]*Rho_M[0][i]/(Rho_M[3][i] +K[1]) -(A[0][1]*Rho_M[0][i]/(1 + A01H01*Rho_M[0][i]))*Rho_M[1][i];
		//Equivalent to dG/dt = ej*(aij*V*P)/(1+aij*hij*V)  -ajm*G*P/(1+ajm*hjm*G).
		f[1][i] = (EA1*Rho_M[0][i]/(1 + A01H01*Rho_M[0][i]))*Rho_M[1][i]  
					-A[1][2]*Rho_M[1][i]/(1 + A12H12*Rho_M[1][i])*Rho_M[2][i];
        //Equivalent to dPr/dt = em*(ajm*G*Pr)/(1+ajm*hjm*G)
		f[2][i] = (EA2*Rho_M[1][i]/(1 + A12H12*Rho_M[1][i]))*Rho_M[2][i]; //(1 - Rho_M[2][i]*K_P1)*  -Mm*Rho_M[2][i]*Rho_M[2][i];
		//Equivalent to dW/dt = alpha*(P+K2*W0)/(P+K2)*O - g_max*P*W/(W+K1) - rW*W + D*(Laplacian of W)
        f[3][i] = alpha*(Rho_M[0][i]+ K2W0)/(Rho_M[0][i] +K[2])*Rho_M[4][i] -gmax*Rho_M[3][i]*Rho_M[0][i]/(Rho_M[3][i] +K[1]) - rW*Rho_M[3][i] 
        + (Dxd2[3])*(Rho_M[3][nR2[i][0][0]]  + Rho_M[3][nR2[i][0][1]]  + Rho_M[3][nR2[i][1][0]]  + Rho_M[3][nR2[i][1][1]]  - 4*Rho_M[3][i]);
        //Equivalent to dO/dt = a - alpha*(P+K2*W)/(P+K2)*O + D*(Laplacian of O)
		f[4][i] = a - alpha*(Rho_M[0][i]+ K2W0)/(Rho_M[0][i] +K[2])*Rho_M[4][i] + (Dxd2[4])*(Rho_M[4][nR2[i][0][0]]  + Rho_M[4][nR2[i][0][1]] 
        + Rho_M[4][nR2[i][1][0]]  + Rho_M[4][nR2[i][1][1]]  - 4*Rho_M[4][i]);
		
	}
}

void RK4_Integrate_Stochastic_MultiSp(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D2Vec_Double &K1, D2Vec_Double &K2, D2Vec_Double &K3, D2Vec_Double &K4, D3Vec_Int &nR2,
		double a,double c,double gmax,double alpha, double rW, double W0, double (&Dxd2)[Sp], double (&K)[3], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t,double dt,double dx, int g)
{
	//double dt6 = dt/6.0; double dt2 = dt/2.0; double dx1_2 = 1/(dx*dx);

	/** // OLD DECLARATIONS of K1, K2, K3, K4 and Rho_M. Skipped as declaration and assignment is done in rietkerk_Dornic_2D_MultiSp
	 * and doing it here takes O(n^2) time each time this function is called.
	D2Vec_Double K1(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K2(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double K3(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K4(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double Rho_M(Sp, vector<double> (g*g, 0.0));
	*/

	// RECALL: Dxd2[s] = D[s]/dx2
	f_2Dor_3Sp(K1, Rho_t, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A,H,E, t, dt, dx1_2, g); //K1 updated.

	for(int s= 0; s <Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K1[s][i];
	}

	f_2Dor_3Sp(K2, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A,H,E, t + dt2, dt, dx1_2, g); //K2 updated.

	for(int s= 0; s <Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K2[s][i];
	}
	f_2Dor_3Sp(K3, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A,H,E, t + dt2, dt, dx1_2, g); //K3 updated.
	for(int s= 0; s <Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for(int i=0; i < g*g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt)*K3[s][i];
	}
	f_2Dor_3Sp(K4, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A,H,E, t + dt, dt, dx1_2, g); //K4 updated.
    
	for(int s= 0; s <Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for(int i=0; i < g*g; i++)
		{
			Rho_t[s][i]+= (dt6)*( K1[s][i] + 2.0*K2[s][i] + 2.0*K3[s][i] + K4[s][i]);
		}
		#ifdef DEBUG
		bool isKoed = false; // Flag to indicate if the RK4 integration was "KO'ed" due to invalid values.
		for(int i=0; i < g*g; i++)
		{
			if( Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || isnan(Rho_t[s][i]) == true)
			{
				if(isKoed == false)
				{
					std::ofstream errout; //Save Error Logs
					std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
					stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
					m6 << "RK4 WAS KO'ED WITH VAL:\t" << Rho_t[s][i] << " For [s, thr, a, i, Pre]\t" << s << " , " << omp_get_thread_num() << " , " << a 
					<< " , " << i << " , " << prefix << " ] with K1[s][i]:   " << K1[s][i] << "\t, K2[s][i]:\t" << K2[s][i] 
					<< "\t, K3[s][i]:\t" << K3[s][i] << "\t, K4[s][i]:\t" << K4[s][i] << "\t AND Rho(t-dt)" << Rho_tsar[s][i] << "\n"; //cout << m6.str();
					//errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					cerr << m6.str();
					isKoed = true;
				}
			}
		}
		#endif //DEBUG
	}
}
// Function to save the preliminary frame data to disk.
// NOTE: Rho is assumed to be a matrix with the following structure:
//  < Sp Independent Cols (SpB_indepnt_cols in number) >  < Sp Dependent Cols which are SpB*(spcols) in number >
void save_prelimframe(D2Vec_Double &Rho_t, const string &parendir, const string &filenamePattern, double a, double a_st, double a_end, double t, 
	double dt, double dx, double dP, int r, int g, string header /* =""*/, bool overwrite /* = false*/, bool delete_previous /* = false*/)
{
	stringstream rini; rini << r;
	//Next, designate the naive filename
	string basefilename = filenamePattern + rini.str() + ".csv";

	string filename = parendir + basefilename;
	// If the file already exists, increment the replicate number by 1 of the maximum replicate number corresponding to the filename pattern.
	// This is to avoid overwriting files.
	if (fs::exists((filename)) && !overwrite) 
	{
		int maxRepNo = findMaxRepNo(parendir, filenamePattern);
		//maxRepNo++; //Increment the maximum replicate number by 1.
		filename = parendir + filenamePattern + std::to_string(maxRepNo + 1) + ".csv";
		stringstream m1;

		m1 << "File exists. New filename: " << filename << " \n and maxRepNo: " << maxRepNo << "\n";
		cout << m1.str();
		basefilename = filenamePattern + std::to_string(maxRepNo + 1) + ".csv";
	}

	

	std::ostringstream oss;
	if(header != "")
		oss << header;
	else
		oss << frame_header;
		// Using default header (defined in header file) if no header is provided.

	//tot_iter is the number of rows in the Rho_t matrix.
	double tot_rows = Rho_t.size();
	for(int i=0; i< tot_rows; i++)
	{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |
		oss << a << ","<< r+1 << "," << g << ","; // << Rho_t[i][0] << "," << Rho_t[i][1] << "," << Rho_t[i][2] << ",";
		
		// Iterate over size of Rho_t[i] and write the dependent columns to the file.
		for(int k=0; k < Rho_t[i].size(); k++)
			oss << Rho_t[i][k] << ",";

		oss <<  "\n";
		/** OLD CODE: Write the dependent columns to the file. 
		oss << Rho_t[i][0] << "," << Rho_t[i][1] << "," << Rho_t[i][2] << ","

		for(int s=0; s <SpB; s++)
		{
			oss << Rho_t[i][4*s + 3] << "," << Rho_t[i][4*s +4] << "," << Rho_t[i][4*s +5] << "," << Rho_t[i][4*s +6] << ","; 
		}
		//*/
		
	}

	FILE *fp = fopen(filename.c_str(), "w");
	if(fp)
	{
		fprintf(fp, "%s", oss.str().c_str());
		fclose(fp);
		cout << "File: " << basefilename << " written to disk successfully. \n";
	}
	else
	{
		std::cerr << "Error: Could not open PRELIM file " << filename << std::endl;
	}

	// Delete the previous file (if it exists) to save space.
	#if defined(__GNUC__) && (__GNUC__ >= 9)

	if(r > 2 && delete_previous)
	{
		// Delete the previous file (if it exists) to save space.
		// Reset rini value.
		stringstream rini_prev; rini_prev << r-2;
		// Get substring from the string filenamePattern starting after the sequence "_DP_G_" to the end of the string.
		std::smatch m;
		std::regex_search(filenamePattern, m, std::regex("_DP_G_"));
		string filenamePattern_prev = prelim_prefix + ".*_DP_G_" +  m.suffix().str() + rini_prev.str() + ".csv";
		//stringstream m7; m7 << "For (r,a,t)" << r << " , " << a << " , " << t 
		//<< " Previous filename pattern: " << filenamePattern_prev << "with m: " << m.str() << "\n"; cout << m7.str();

		// Using regex to delete the file with the previous replicate number.
		std::regex filePattern(filenamePattern_prev);

		// Iterate over the files in the parent directory (path_to_dir) and delete the file if it matches the pattern.
		for (const auto & entry : fs::directory_iterator(parendir))
		{
			string matchStr = "/" + entry.path().filename().string();
			if (std::regex_search(matchStr, filePattern))
			{
				if(fs::exists(entry.path()) && fs::is_regular_file(entry.path()))
				{
					try
					{	
						fs::remove(entry.path());
						stringstream m8;
						m8 << "File " << entry.path().filename().string() << " deleted successfully. \n"; cout << m8.str();
					}
					catch (const fs::filesystem_error& e)
					{
						stringstream m8;
						m8 << "Error deleting file: " << entry.path().string() << " with error: " << e.what() << "\n"; cout << m8.str();
						cerr << m8.str();
					}
				}
				else
				{
					stringstream m8;
					m8 << "File: " << entry.path().string() << " does not exist. \n"; cout << m8.str();
					cerr << m8.str();
				}	
			}
		}	
	}
	#endif
}

void save_frame(D2Vec_Double &Rho_t, const string &parendir, const string &filenamePattern, double a, double a_st, double a_end, double t, double dt, double dx, double dP, int r, int g, string header /* =""*/, bool overwrite /* = false*/)
{
	stringstream rini; rini << r;
	//Next, designate the naive filename
	string basefilename = filenamePattern + rini.str() + ".csv";

	string filename = parendir + basefilename;
	// If the file already exists, increment the replicate number by 1 of the maximum replicate number corresponding to the filename pattern.
	// This is to avoid overwriting files.
	if (fs::exists((filename)) && !overwrite) 
	{
		int maxRepNo = findMaxRepNo(parendir, filenamePattern);
		//maxRepNo++; //Increment the maximum replicate number by 1.
		filename = parendir + filenamePattern + std::to_string(maxRepNo + 1) + ".csv";
		stringstream m1;

		m1 << "File exists. New filename: " << filename << " \n and maxRepNo: " << maxRepNo << endl;
		cout << m1.str();
	}

	std::ostringstream oss;
	if(header != "")
		oss << header;
	else
		oss << frame_header;
		// Using default header (defined in header file) if no header is provided.

	double Sp = Rho_t.size();
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


void save_framefileswrapper(int index, int tot_iter, int j, int thrID,  double t, double dt, double dx, D2Vec_Double &Rho_dt, D2Vec_Double &DRho,
	D2Vec_Double &rho_rep_avg_var, vector<double> &t_meas, D2Vec_Double &gamma, D2Vec_Pair_Double &v_eff, double t_max, 
	double a, double c, double gmax, double alpha, double rW, double W0, double (&D)[Sp], double (&v)[SpB], 
	double (&K)[3], double sigma[], double a_st, double a_end, double a_c, double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], 
    double dP, int r, int g, double Gstar /* =-1.*/, double Vstar /* = -1.*/)
{
	// FRAME SAVING

	//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.
	
	stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, Dm2, alph, w0t, aij, hij, dix, dimitri, sig0, sig1, geq, veq, thr_ID;

	L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
	rini << j; Dm0 << D[0]*pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2]; thr_ID << thrID; 
	a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1]; dimitri  << dP; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;
	//gm << setprecision(3) << gmax; w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; 
	// Three replicates are over.

	string parendir = "";

	// Creating a file instance called output to store output data as CSV.
	if(Gstar != -1)
		parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str();
	else if(Vstar != -1)
		parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str();
	else
		parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();
	
	string filenamePattern = frame_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
	+ "_a_" + p1.str()  + "_D1_"+ Dm1.str() + "_jID_"+ thr_ID.str() + "_dx_"+ dix.str() + "_R_";
	//Next, designate the naive filename
	string filename = filenamePattern + rini.str() + ".csv";

	#if SPB > 1
		string gammafilenamePattern = gamma_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() 
		+ "_a_" + p1.str()  + "_D1_"+ Dm1.str() + "_jID_"+ thr_ID.str() + "_dx_"+ dix.str() + "_R_";
	#endif
	
	// SAVE ALL FRAMES
	save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
	#if SPB > 1 // Save gamma frames.
		save_frame(gamma, parendir, gammafilenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
	#endif
	stringstream m3;
	m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
	cout << m3.str();
	// */
	

	// SAVE SELECTED FRAMES
	/** FRAME SAVING STRATEGY FOR LOW DT (DT ~ 0.01)
	if(t < 250)
	{
		//Only save one in 4 frames here.
		if(index%4 ==0 || t== 0)
		{
			save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
			#if SPB > 1 // Save gamma frames.
				save_frame(gamma, parendir, gammafilenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
			#endif
			stringstream m3;
			m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
			cout << m3.str(); 
		}
	}
	else if(t >= 250 && t < 600)
	{	//Only save one in two frames here.
		if(index%2 ==0)
		{
			save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
			#if SPB > 1 // Save gamma frames.
				save_frame(gamma, parendir, gammafilenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
			#endif
			stringstream m3;
			m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
			cout << m3.str(); 
		}
	} //*/

	/** PREVIOUS FRAME SAVING STRATEGY (DT ~ 0.1)
	if(t < 760)
	{
		//Only save one in three frames here.
		if(index%3 ==2 || t== 0)
		{
			save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
			#if SPB > 1 // Save gamma frames.
				save_frame(gamma, parendir, gammafilenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
			#endif
			stringstream m3;
			m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
			cout << m3.str();
		}
	}
	else if(t >= 760 && t < 6000)
	{	//Only save one in two frames here.
		if(index%2 ==0)
		{
			save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
			#if SPB > 1 // Save gamma frames.
				save_frame(gamma, parendir, gammafilenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
			#endif

			stringstream m3;
			m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
			cout << m3.str(); 
		}
	}//*/

	/** FRAME SAVING STRATEGY FOR DT WITH FRAME_MEAS WITH LINEAR WINDOWS
	if(t > 1500 && t < 6000)
	{	//Only save one in two frames here.
		if(index%2 ==0)
		{
			save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
			#if SPB > 1 // Save gamma frames.
				save_frame(gamma, parendir, gammafilenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
			#endif

			stringstream m3;
			m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
			cout << m3.str(); 
		}
	}
	else
	{
		//Save all frames here.
		save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
		#if SPB > 1 // Save gamma frames.
			save_frame(gamma, parendir, gammafilenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
		#endif
		stringstream m3;
		m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
		cout << m3.str(); 
	}
	//*/
	

	//
}

// A wrapper function to save the prelim files for each thread, called from the main simulation loop in rietkerk_Dornic_2D_MultiSp()
// Also saves future prelim files if vegetation extinction is detected.
int save_prelimfileswrapper(int index, int tot_iter, int j, int thrID,  double t, double dt, double dx, D2Vec_Double &Rho_dt, D2Vec_Double &DRho, D2Vec_Double &Rho_M, D2Vec_Double &Rho_Mov,
	D2Vec_Double &rho_rep_avg_var, vector<double> &t_meas, D2Vec_Double &gamma, D2Vec_Pair_Double &v_eff, double t_max, 
	double a, double c, double gmax, double alpha, double rW, double W0, double (&D)[Sp], double (&v)[SpB], 
	double (&K)[3], double sigma[], double a_st, double a_end, double a_c, double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], 
    double dP, int r, int g, double Gstar /* =-1.*/, double Vstar /* = -1.*/)
{

	// BLOCK FOR CALCULATING AND TEMP PRELIMINARY FRAMES
	
	// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
	// and save to file.
	D2Vec_Double rho_rep_avg_var_temp(index+1, vector<double> (Sp4_1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
	//std::copy(rho_rep_avg_var.begin(), rho_rep_avg_var.begin() + index, rho_rep_avg_var_temp.begin());
	var_mean_incremental_surv_runs(rho_rep_avg_var_temp, Rho_M, index+1, 0);

	#if SPB > 1 // Save movement time-series for each species in Rho_Mov.
		D2Vec_Double gamma_mov_avg_temp(index+1, vector<double> (int(Rho_Mov[0].size()) + 1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
		// Get cols_per_sp by dividing the number of columns in Rho_Mov by the number of species (SpB).
		int cols_per_sp = int(Rho_Mov[0].size()/SpB);
		generic_SPBmean_surv_runs(gamma_mov_avg_temp, Rho_Mov, index+1, cols_per_sp, 0);
		//Time is stored in the first column.
		for(int i=0; i < index+1; i++)
			gamma_mov_avg_temp[i][0] = t_meas[i];

	#endif


	//Time is stored in the first column.
	for(int i=0; i < index+1; i++)
		rho_rep_avg_var_temp[i][0] = t_meas[i];

	//Finally save to file.
	stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, Dm, geq, veq, jID; // cgm, sig0;
	a1 << a_st; a2 << a_end;
	L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(5) << a; dimitri << dP; jID << thrID;
	rini << j; Dm << setprecision(4) << D[2]; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;

	string parendir= "";
	if(Gstar != -1)
		parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str() + "/TimeSeries";
	else if(Vstar != -1)
		parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str() + "/TimeSeries";
	else
		parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/TimeSeries";

	//double ran_jid = (unif(rng)*1000.0)/1000.0; jID << ran_jid; // Random number between 0 and 1.
	
	string filenamePattern = replicate_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
	"_jID_"+ jID.str() + "_R_";

	save_prelimframe(rho_rep_avg_var_temp, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);

	#if SPB > 1 // Save movement time-series for each species in Rho_Mov.
	string filenamePattern_mov = movement_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
	"_jID_"+ jID.str() + "_R_";
	save_prelimframe(gamma_mov_avg_temp, parendir, filenamePattern_mov, a, a_st, a_end, t, dt, dx, dP, j, g, movmheader, false, false);
	#endif

	

	// Check if no vegetation is left at this current index. If so, break out of the while loop.
	if( Rho_M[index][0] == 0.0)
	{
		stringstream m3_1;
		m3_1 << "RUN-TIME WARNING: ZERO Active VEG sites for TIME [t, a, thr, j]\t" 
		<< t << " , " << a << " , " << omp_get_thread_num()  << " , " << j << " Skipping to next replicate .... \n"; 
		cout << m3_1.str(); cerr << m3_1.str();

		//Next save future preliminary frames with 0.0 values for population densities after this time point.
		int index_prelim_range[5] = {int(tot_iter*0.6), int(tot_iter*0.85), int(tot_iter*0.9), int(tot_iter*0.95), int(tot_iter-1)};
		int iter_index =0;
		for(int iter_index =0; iter_index < 5; iter_index++)
		{	
			int new_index = index_prelim_range[iter_index];
			if( index > new_index)
				continue; //Skip to next iteration if index is greater than the current range.
			// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
			// and save to file.
			D2Vec_Double rho_rep_avg_var_future_temp(new_index+1, vector<double> (Sp4_1, 0.0)); 
			//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
			var_mean_incremental_surv_runs(rho_rep_avg_var_future_temp, Rho_M, new_index+1, 0);
			stringstream tm_future; tm_future << t_meas[new_index]; // Time at which the frame is saved.

			//Time is stored in the first column.
			for(int i=0; i < new_index+1; i++)
				rho_rep_avg_var_future_temp[i][0] = t_meas[i];

			string filenamePattern_future = replicate_prefix + L.str() + "_T_" + tm_future.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
			"_jID_"+ jID.str() + "_R_";

			//Finally save to file.
			save_prelimframe(rho_rep_avg_var_future_temp, parendir, filenamePattern_future, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);

			#if SPB > 1 // Save movement time-series for each species in Rho_Mov.
				D2Vec_Double gamma_mov_avg_future_temp(new_index+1, vector<double> (int(Rho_Mov[0].size()) + 1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
				// Get cols_per_sp by dividing the number of columns in Rho_Mov by the number of species (SpB).
				int cols_per_sp = int(Rho_Mov[0].size()/SpB);
				generic_SPBmean_surv_runs(gamma_mov_avg_future_temp, Rho_Mov, new_index+1, cols_per_sp, 0);
				//Time is stored in the first column.
				for(int i=0; i < new_index+1; i++)
					gamma_mov_avg_future_temp[i][0] = t_meas[i];

				string filenamePattern_mov_future = movement_prefix + L.str() + "_T_" + tm_future.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
				"_jID_"+ jID.str() + "_R_";
				save_prelimframe(gamma_mov_avg_future_temp, parendir, filenamePattern_mov_future, a, a_st, a_end, t, dt, dx, dP, j, g, movmheader, false, false);
				vector<vector<double>>().swap(gamma_mov_avg_future_temp); //Flush temp out of memory.
			#endif

			vector<vector<double>>().swap(rho_rep_avg_var_future_temp); //Flush temp out of memory.
		}
		return 1; //Indicate that the current replicate has no vegetation left at this time point (and break out of the while loop).
	}

	vector<vector<double>>().swap(rho_rep_avg_var_temp); //Flush temp out of memory.
	 
	// */

	return 0;


}

// FFTW3-based gamma sweep calculation for 3-Species Non-Refugia model.
void calc_gamma_3Sp_NonRefugia_FFT(D2Vec_Double& Rho_t, D2Vec_Double& gamma, double (&Rho_avg)[Sp], vector<std::pair<double, int>>& rfrac, 
    vector<std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L, FFTW3_CentralPlanner& fft_central_plan)
    
{
	mutex errorMutex; // Mutex to make error messages thread-safe
    int L2 = L * L;
    double eps = 1.0e-12;
    //int new_length = int(nVeg_frac * nVeg_frac * centralNeighboringSites.size());
    //if (new_length < 1) new_length = 1;

	int rmax_eff = int(nVeg_frac * r_max); //Effective r_max based on vegetation fraction.
	

    double rho_inverse[SpB];
    for (int s = 0; s < SpB; s++)
        Rho_avg[s] >= eps ? rho_inverse[s] = 1 / Rho_avg[s] : rho_inverse[s] = 1 / eps; //Avoid division by zero.
    
	for (const auto& frac : rfrac) {
        int s = frac.second;
        double fr = frac.first;

		int r_eff = int(fr * rmax_eff); //Effective r based on species fraction and vegetation cover.
		int eff_nR_Perp_size = nR_Perp_R[r_eff]; //Get precomputed effective nR_Perp size.
        //int eff_nR_Perp_size = int(fr * fr * new_length);
        if (fr == 0.0) continue; //Skip if fraction is zero.
		if (s == 0) continue; //Skip vegetation.
        if (eff_nR_Perp_size < 1) eff_nR_Perp_size = 1;
        if (Rho_avg[s] < eps) 
		{
			#pragma omp simd
            for (int i = 0; i < L2; i++) gamma[s][i] = 0.0;
            continue;
        }
        if (Rho_avg[s-1] < eps) 
		{
			#pragma omp simd
            for (int i = 0; i < L2; i++) gamma[s][i] = 0.0;
            continue;
        }
        
        // Prepare spatial field based on species
        if (s == 1) 
		{
            if (Rho_avg[2] < eps) 
			{	// Copy vegetation density directly
                std::copy(Rho_t[s-1].begin(), Rho_t[s-1].end(), fft_central_plan.rho_spatial);
            } else {
                // Apply predator avoidance term
				#pragma omp simd
                for (int i = 0; i < L2; i++)
				{
            		fft_central_plan.rho_spatial[i] = Rho_t[s-1][i] * (1.0 - Rho_t[2][i] * rho_inverse[2]);
                }
            }
        } else if (s == 2) {
            // Copy grazer density
            std::copy(Rho_t[s-1].begin(), Rho_t[s-1].end(), fft_central_plan.rho_spatial);
        }
        
        // Forward FFT of density (Rho) field
        fftw_execute_dft_r2c(fft_central_plan.forward_rho, fft_central_plan.rho_spatial, fft_central_plan.rho_fft);
        
        // Get precomputed kernel FFT
        fftw_complex* kernel_fft = kHat_kernel_FFT_cache[r_eff];
        
		double real_part, imag_part;
        // Multiply in Fourier space (complex multiplication)
		#pragma omp simd
        for (int i = 0; i < L2; i++) 
		{
			// RECALL: Z = (a + ib)*(c + id) = (ac - bd) + i(ad + bc)
            real_part = fft_central_plan.rho_fft[i][0] * kernel_fft[i][0] - 
                fft_central_plan.rho_fft[i][1] * kernel_fft[i][1];
            imag_part = fft_central_plan.rho_fft[i][0] * kernel_fft[i][1] + 
                fft_central_plan.rho_fft[i][1] * kernel_fft[i][0];
            fft_central_plan.prod_fft[i][0] = real_part; 
			fft_central_plan.prod_fft[i][1] = imag_part;
        }
        
        // Inverse FFT to get convolution result
        fftw_execute_dft_c2r(fft_central_plan.backward, fft_central_plan.prod_fft, fft_central_plan.convolution_spatial);
        // Normalize and apply species-specific scaling.
        double normalisation = 1.0 / L2;  // FFTW3 doesn't normalise - multiplication by 1/(L*L) is needed.
		#pragma omp simd
        for (int i = 0; i < L2; i++) 
		{
            gamma[s][i] = fft_central_plan.convolution_spatial[i] * normalisation * rho_inverse[s-1];
		}
		// Clamp to [0, 1]
        for (int i = 0; i < L2; i++) 
		{   
			if (gamma[s][i] > 1.0) gamma[s][i] = 1.0;
            else if (gamma[s][i] < 0.0) gamma[s][i] = 0.0;
			else if (isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
			{
				errorMutex.lock();
				std::cerr << "Gamma is NaN with value: " << gamma[s][i] << " For [s, thr, i]\t" << s << " , " << omp_get_thread_num()
					<< " , " << i << " ] with Rho_avg[0]: " << Rho_avg[0] << " and Rho_avg[1]: " << Rho_avg[1] << " and Rho_avg[2]: " << Rho_avg[2]
					<< " \n and Rho_t[0][i]: " << Rho_t[0][i] << " and Rho_t[1][i]: " << Rho_t[1][i] << " and Rho_t[2][i]: " << Rho_t[2][i]
					<< " \n and Nveg_frac: " << nVeg_frac << " and r_eff: " << r_eff << " and rfrac.size(): " << rfrac.size() << " and nR_Perp_R[r_eff]: " << nR_Perp_R[r_eff]
					<< std::endl;
				errorMutex.unlock();
			}
        }
    }
}

//Calculates gamma for 3 Sp Rietkerk model (details in PDF, 
// ASSUMING PREDATORS AREN'T DETERRED BY HIGH LOCAL VEGETATION DENSITY)
void calc_gamma_3Sp_NonRefugia(const vector<pair<int, int>>& centralNeighboringSites, D2Vec_Double &Rho_t, D2Vec_Double &gamma,  
			double (&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, vector <std::pair<int, int>>& dtV_counter,  double nVeg_frac, int r_max, int L)
{	
	//if(nVeg_frac < 0.35)
	//	nVeg_frac = 0.35; //Minimum fraction of vegetation cover.

	int L2 = L*L;
	double eps = 1.0e-12; //Small number to avoid division by zero.

	int new_length = int(nVeg_frac*nVeg_frac*centralNeighboringSites.size());
	new_length < 1 ? new_length = 1 : new_length = new_length; // At least one nearest neighbor (the site itself) is considered.
	// Make a copy of the centralNeighboringSites, resized to the fraction nVeg_frac*nVeg_frac.
	auto range_CentralNeighboringSites = std::ranges::subrange(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);
	//vector<pair<int, int>> centralNeighboringSitesCopy(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);

	std::vector<int> nR_Perp(new_length);

	double rho_inverse[SpB];// ={1/eps, 1/eps, 1/eps}; //Inverse of average density.
	for(int s=0; s< SpB; s++)
		Rho_avg[s] >= eps ? rho_inverse[s] = 1/Rho_avg[s] : rho_inverse[s] = 1/eps; //Avoid division by zero.
	// Avoid division by zero.

	#if defined(DEBUG)
	stringstream m1;
	m1 << " calc_gamma_3sp UPDATE, for thr: " << omp_get_thread_num() << " WITH nVeg_frac: " << nVeg_frac << " , new_length: " << new_length 
	<< " and AVG( Rho[0], Rho[1], Rho[2]): " << Rho_avg[0] << " , " << Rho_avg[1] << " , " << Rho_avg[2] << " and nR_Perp.size(): " << nR_Perp.size() << "\n";
	cout << m1.str(); cerr << m1.str();
	#endif

	for( int i=0; i < L2; i++)
	{
		mutex errorMutex; // Mutex to make error messages thread-safe
		
		int c_i = int(i/L); int c_j = i%L; //Current x and y coordinates of site i.

		generateNeighboringSitesFromCentral(range_CentralNeighboringSites, nR_Perp, c_i, c_j, L);

		//if(nVeg_frac < 1)
		//	nR_Perp.resize(int(nVeg_frac*nVeg_frac*nR_Perp.size()));
		// In this case, r_effective is reduced to nVeg_frac*r_max, hence the reduction in the number of nearest neighbours.
		// rfrac is sorted in descending order.
		double fr_prev = 1.0; //Fraction of perception radius to max radius.
		for(const auto& frac: rfrac)
		{
			int s = frac.second;
			double fr = frac.first; //Fraction of perception radius to max radius.

			int eff_nR_Perp_size = int(fr*fr*nR_Perp.size());
			
			if( fr == 0.0)
			{	continue;	} // No perception radius for this species (indicates vegetation species)
			if( eff_nR_Perp_size < 1 && fr > 0)
			{	eff_nR_Perp_size = 1;	} // At least one nearest neighbor (the site itself) is considered.
			if(Rho_avg[s] < eps)
				{	gamma[s][i] = 0.0; continue;	} // No species left, value assigned doesn't matter.
			else if(Rho_avg[s-1] < eps)
				{	gamma[s][i] = 0.0;	continue; } // No resource left, consumers will advect to extinction.
			
			/** 
			if(dtV_counter[s].first != 1 && dtV_counter[s].second != 0)
			{	continue;	} 
			//*/
			//Skip this species if it is not updated at this time step (depending on advection timescale)
			double gamma_sum = 0.0; //Sum of gamma values for this species at site i (needed for SIMD reduction).
			// Gamma for grazers (indexed 1)
			if (s == 1)
			{	
				
				if( Rho_avg[2] < eps)
				{
					#pragma omp simd reduction(+:gamma_sum )	
					for(int k=0; k< eff_nR_Perp_size; k++)
					{		gamma_sum  += Rho_t[s-1][nR_Perp[k]];	}
				} // No predators left
				else
				{
					//double rho_pred1 = 1/Rho_avg[2];
					#pragma omp simd reduction(+:gamma_sum )
					for(int k=0; k< eff_nR_Perp_size; k++)
					{		
						gamma_sum += Rho_t[s-1][nR_Perp[k]]*(1 -Rho_t[2][nR_Perp[k]]*rho_inverse[2]);	
					}	
				}
				// RECALL : rho_inverse[s] = 1/Rho_avg[s]
				gamma[s][i] = (gamma_sum *rho_inverse[s-1])/(eff_nR_Perp_size);	

			}
			if(s == 2)
			{
				#pragma omp simd reduction(+:gamma_sum )
				for(int k=0; k< eff_nR_Perp_size; k++)
					{		gamma_sum  += Rho_t[s-1][nR_Perp[k]];	}
				// RECALL : rho_inverse[s] = 1/Rho_avg[s]
				gamma[s][i] = (gamma_sum /(eff_nR_Perp_size))*rho_inverse[s-1];	
			}

			if(gamma[s][i] > 1)
			{	gamma[s][i] = 1.0;	} //Ensuring gamma is not greater than 1.
			else if(gamma[s][i] < 0)
			{	gamma[s][i] = 0.0;	} //Ensuring gamma is not less than 0.
			else if(isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
			{
				errorMutex.lock();
				std::cerr << "Gamma is NaN with value: " << gamma[s][i] << " For [s, thr, i]\t" << s << " , " << omp_get_thread_num() 
						<< " , " << i << " ] with Rho_avg[0]: " << Rho_avg[0] << " and Rho_avg[1]: " << Rho_avg[1] << " and Rho_avg[2]: " << Rho_avg[2] 
						<< " \n and Rho_t[0][i]: " << Rho_t[0][i] << " and Rho_t[1][i]: " << Rho_t[1][i] << " and Rho_t[2][i]: " << Rho_t[2][i] 
						<< " \n and Nveg_frac: " << nVeg_frac << " and r_max: " << r_max << " and rfrac.size(): " << rfrac.size() << " and nR_Perp.size(): " << nR_Perp.size()
						<< std::endl;
				errorMutex.unlock();
			}
		}
		//vector <int>().swap(nR_Perp);
	}
}


//Calculates gamma for 3 Sp Rietkerk model (details in PDF) [ 0 < ]
void calc_gamma_3Sp(const vector<pair<int, int>>& centralNeighboringSites, D2Vec_Double &Rho_t, D2Vec_Double &gamma,  
			double (&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, vector <std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L)
{	
	if(nVeg_frac < 0.35)
		nVeg_frac = 0.35; //Minimum fraction of vegetation cover.
	
	int new_length = int(nVeg_frac*nVeg_frac*centralNeighboringSites.size());
	new_length < 1 ? new_length = 1 : new_length = new_length; // At least one nearest neighbor (the site itself) is considered.

	auto range_CentralNeighboringSites = std::ranges::subrange(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);
	
	// Make a copy of the centralNeighboringSites, resized to the fraction nVeg_frac*nVeg_frac.
	//vector<pair<int, int>> centralNeighboringSitesCopy(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);

	std::vector<int> nR_Perp(new_length);


	int L2 = L*L;
	double eps = 1.0e-12; //Small number to avoid division by zero.

	double rho_inverse[SpB];// ={1/eps, 1/eps, 1/eps}; //Inverse of average density.
	for(int s=0; s< SpB; s++)
		Rho_avg[s] >= eps ? rho_inverse[s] = 1/Rho_avg[s] : rho_inverse[s] = 1/eps; //Avoid division by zero.
	// Avoid division by zero.
	for( int i=0; i < L2; i++)
	{
		mutex errorMutex; // Mutex to make error messages thread-safe

		int c_i = int(i/L); int c_j = i%L; //Current x and y coordinates of site i.
		generateNeighboringSitesFromCentral(range_CentralNeighboringSites, nR_Perp, c_i, c_j, L);

		//if(nVeg_frac < 1)
		//	nR_Perp.resize(int(nVeg_frac*nVeg_frac*nR_Perp.size()));
		// In this case, r_effective is reduced to nVeg_frac*r_max, hence the reduction in the number of nearest neighbours.
		// rfrac is sorted in descending order.
		double fr_prev = 1.0; //Fraction of perception radius to max radius.
		for(const auto& frac: rfrac)
		{
			int s = frac.second;
			double fr = frac.first; //Fraction of perception radius to max radius.

			int eff_nR_Perp_size = int(fr*fr*nR_Perp.size());
			// If fr < 1, reduce the number of nearest neighbors to the fraction of the original number of nearest neighbors.

			if( fr == 0.0)
			{	continue;	} // No perception radius for this species (indicates vegetation species)
			if( eff_nR_Perp_size < 1 && fr > 0)
			{	eff_nR_Perp_size = 1;	} // At least one nearest neighbor (the site itself) is considered.
			if(Rho_avg[s] < eps)
				{	gamma[s][i] = 0.0; continue;	} // No species left, value assigned doesn't matter.
			else if(Rho_avg[s-1] < eps)
				{	gamma[s][i] = 0.0;	continue; } // No resource left, consumers will advect to extinction.

			/** 
			if(dtV_counter[s].first != 1 && dtV_counter[s].second != 0)
			{	continue;	} 
			//Skip this species if it is not updated at this time step (depending on advection timescale)
			//*/
			double gamma_sum = 0.0; //Sum of gamma values for this species at site i (needed for SIMD reduction).
			// Gamma for grazers (indexed 1)
			if (s == 1)
			{	
				if( Rho_avg[2] < eps)
				{	
					#pragma omp simd reduction(+:gamma_sum)
					for(int k=0; k< eff_nR_Perp_size; k++)
					{		gamma_sum += Rho_t[s-1][nR_Perp[k]];	}
				} // No predators left
				else
				{
					//if(Rho_avg[0] == 0)
					//{	gamma[s][i] = 0.0;	continue; } // No vegetation 
					//double rho_pred1 = 1/Rho_avg[2];
					#pragma omp simd reduction(+:gamma_sum)
					for(int k=0; k< eff_nR_Perp_size; k++)
					{
								
						gamma_sum += Rho_t[s-1][nR_Perp[k]]*(1 -Rho_t[2][nR_Perp[k]]*rho_inverse[2]);	
					}
				// RECALL : rho_inverse[s] = 1/Rho_avg[s]
				}
				gamma[s][i] = (gamma_sum*rho_inverse[s-1])/(eff_nR_Perp_size);	

			}
			if(s == 2)
			{
				
				if( Rho_avg[0] < eps)
				{	
					#pragma omp simd reduction(+:gamma_sum)	
					for(int k=0; k< eff_nR_Perp_size; k++)
					{		gamma_sum += Rho_t[s-1][nR_Perp[k]];	}
				} // No vegetation left
				else
				{
					//if(Rho_avg[1] < eps)
					//{	gamma[s][i] = 0.0;	continue; } // No grazers left.
					//double rho_veg1 = 1/Rho_avg[0];
					#pragma omp simd reduction(+:gamma_sum)
					for(int k=0; k< eff_nR_Perp_size; k++)
					{
						gamma_sum += Rho_t[s-1][nR_Perp[k]]*(1 -Rho_t[0][nR_Perp[k]]*rho_inverse[0]);	
					}
				// RECALL : rho_inverse[s] = 1/Rho_avg[s]
				}
				gamma[s][i] = (gamma_sum*rho_inverse[s-1])/(eff_nR_Perp_size);		
			}

			if(gamma[s][i] > 1)
			{	gamma[s][i] = 1.0;	} //Ensuring gamma is not greater than 1.
			else if(gamma[s][i] < 0)
			{	gamma[s][i] = 0.0;	} //Ensuring gamma is not less than 0.
			else if(isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
			{
				errorMutex.lock();
				std::cerr << "Gamma is NaN with value: " << gamma[s][i] << " For [s, thr, i]\t" << s << " , " << omp_get_thread_num() 
						<< " , " << i << " ] with Rho_avg[0]: " << Rho_avg[0] << " and Rho_avg[1]: " << Rho_avg[1] << " and Rho_avg[2]: " << Rho_avg[2] 
						<< " \n and Rho_t[0][i]: " << Rho_t[0][i] << " and Rho_t[1][i]: " << Rho_t[1][i] << " and Rho_t[2][i]: " << Rho_t[2][i] 
						<< " \n and Nveg_frac: " << nVeg_frac << " and r_max: " << r_max << " and rfrac.size(): " << rfrac.size() << " and nR_Perp.size(): " << nR_Perp.size()
						<< std::endl;
				errorMutex.unlock();
			}
		}
	}
	vector <int>().swap(nR_Perp); //vector <pair<int, int>>().swap(centralNeighboringSitesCopy); //Free up memory.
}

void rietkerk_Dornic_2D_MultiSp(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0, 
	double (&D)[Sp], double (&v)[SpB], double (&K)[3], double sigma[], double a_st, double a_end, double a_c, double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double pR[], 
	int (&dtV)[SpB], double clow[], double dt, double dx, double dP, int r, int g, FFTW3_CentralPlanner& fft_central_plan, double Gstar /* =-1.*/, double Vstar /* = -1.*/)
{
	int thrID; //Unique ID for each thread, used for saving files.
	double epsilon = 1.0e-12; //Small number to avoid division by zero.
	double perc =0.015; //Percentage of high density patches.
	//Store max value in pR[] in a variable r_max.
	//NOTE: SpB = Sp - 2; //Number of species excluding water terms.
	double r_max_effective = (*max_element(pR, pR + SpB))/dx;
	// Store pR[]/dx as a fraction of r_max_effective in r_frac[].
	//Initialise vector<std::pair<int, int>> r_frac with initial size SpB and set all values to 0.0.
	vector <std::pair<double, int>> r_frac(SpB, {0.0, 0});
 
	for(int s=0; s< SpB; s++)
	{	r_frac[s] = {pR[s]/(dx*r_max_effective), s};	}

	// Sort the r_frac[] array in descending order and store the indices in r_sort[].
	sort(r_frac.begin(), r_frac.end(), [](std::pair<double, int> &a, std::pair<double, int> &b) 
				{ return a.first > b.first; });
	std::vector<std::pair<int, int>> origin_Neighbourhood = computeNeighboringSitesCentral(int(r_max_effective));

	int nR_Perp_length = origin_Neighbourhood.size();

	#if defined(BARRACUDA) && SPB > 1
	
	// CUDA code, set up memory pointers and allocate memory on the GPU (parralel processing for advdiff_FKE_MultiSp() function)
	double* d_Rho_dt, * d_Rho_tsar, * d_gamma; std::pair<double, double>* d_v_eff;
	//double *d_sigD_ScaleBounds, *d_sigvD_ScaleBounds; // Pointers to the sigma bounds for the Gaussian stencils.
	//double *mu_vel_prefac; // Defined as constant arrays now.
	// RECALL: SpB = Sp - 2; Sp_NV = Sp - 3; // Number of consumers.
	// Allocate memory on the GPU for the various arrays.
	cudaMalloc(&d_Rho_dt, SpB * g * g * sizeof(double)); // Stores the curr population densities of consumers.
	cudaMalloc(&d_Rho_tsar, Sp_NV * g * g * sizeof(double)); // Stores the updated population densities of consumers.
	cudaMalloc(&d_gamma, Sp_NV * g * g * sizeof(double)); // Stores the gamma values for the consumers.
	cudaMalloc(&d_v_eff, Sp_NV * g * g * sizeof(std::pair<double, double>)); // Stores the effective velocities for the consumers.
	//cudaMalloc(&d_sigD_ScaleBounds, Sp_NV*sizeof(double)); // Stores the sigma bounds for the Gaussian stencils.
	//cudaMalloc(&d_sigvD_ScaleBounds, Sp_NV*sizeof(double)); // Stores the sigma bounds for the Gaussian stencils.
	//cudaMalloc(&mu_vel_prefac, Sp_NV*sizeof(double)); // Stores the mu_vel_prefactor[] array.

	// Check if memory allocation was successful.
	if (d_Rho_dt == NULL || d_Rho_tsar == NULL || d_gamma == NULL || d_v_eff == NULL)
	{
		std::cerr << "Memory allocation failed on the GPU for the various arrays in rietkerk_DorFPE_2D_MultiSp() function. Exiting .... \n";
		exit(1);
	}
	// Copy origin_Neighbourhood to global memory on the GPU.
	int2* d_origin_Neighbourhood;
	std::vector<int2> central_neighbours(origin_Neighbourhood.size());
    for (int k = 0; k < origin_Neighbourhood.size(); k++) 
	{
        central_neighbours[k].x = origin_Neighbourhood[k].first;
        central_neighbours[k].y = origin_Neighbourhood[k].second;
    }
	cudaMalloc(&d_origin_Neighbourhood, central_neighbours.size() * sizeof(int2));
	cudaMemcpy(d_origin_Neighbourhood, central_neighbours.data(), 
				central_neighbours.size() * sizeof(int2), cudaMemcpyHostToDevice);
	
	double *d_rho_inverse; int *d_eff_sizes;
	cudaMalloc(&d_rho_inverse, SpB * sizeof(double));
	cudaMalloc(&d_eff_sizes, SpB * sizeof(int));

	vector <DoubleIntPair> dr_frac_copy(SpB, { 0.0, 0 });
	for(int s = 0; s < SpB; s++)
	{
		dr_frac_copy[s].first = r_frac[s].first; dr_frac_copy[s].second = r_frac[s].second;
	}
	//cudaMemcpyToSymbol(d_r_frac, dr_frac_copy.data(), 
	//			dr_frac_copy.size() * sizeof(DoubleIntPair), 0, cudaMemcpyHostToDevice);

	cudaError_t err;

	err = copyToDeviceConstantMemory_GammaSweepTerms(dr_frac_copy.data());

	if (err != cudaSuccess) {
		cerr << "CUDA error in copyToDeviceConstantMemory_AdvTerms: " << cudaGetErrorString(err) << std::endl;
		exit(EXIT_FAILURE);
	}

	if (d_origin_Neighbourhood == NULL || d_rho_inverse == NULL || d_eff_sizes == NULL)
	{
		std::cerr << "Memory allocation failed on the GPU for the movement arrays in dP_Dornic_2D_MultiSp() function. Exiting .... \n";
		exit(1);
	}
	//Finally report Gaussian stencil values on device using reportGaussianStencilValues() function.
	if (omp_get_thread_num() == 1)
	{
		reportPerceptionArrays_Interface(d_origin_Neighbourhood, d_rho_inverse, d_eff_sizes, nR_Perp_length, r_max_effective, g);
	}
	#endif


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

	int frame_tot_iter = frame_tmeas.size(); // Stores total number of iterations (i.e. time steps logged) for frame saving.

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
	//vector <double> gamma_i(SpB, 0.0); //Stores gamma for each species.
	

	beta[0] = -M[0] - 4*D[0]/(dx2);
	lambda_exp[0]= exp( (beta[0])*dt); lambda[0] = 2*(beta[0])/(sigma[0]*sigma[0]*(lambda_exp[0] -1.0));
	diff_coefficient[0] = D[0]/(dx2); sigma2_1[0] = 1/(sigma[0]*sigma[0]);
	for(int s=1; s < SpB; s++)
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
	for(int s= SpB; s < Sp; s++)
		diff_coefficient[s] = D[s]/(dx2);
	//double Gstar = M[2]/((E[2] -M[2]*H[1][2])*A[1][2]); // MFT Coexistance Steady state for grazer.
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

		D2Vec_Double gamma(SpB, vector<double> (g*g, 0.0)); //Stores gamma for each species at each site.

		// 2D vector of pairs of dim: SpBx(L*l), called v. Stores the velocity of each species at each site.
		D2Vec_Pair_Double v_eff(SpB, std::vector<std::pair<double, double>>(g*g, {0.0, 0.0}));

		// 2D vectors of dim: Sx(L*l). Used for RK4 integration (of linear terms) of Rho_dt , after Dornic update.
		D2Vec_Double K1(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K2(Sp, vector<double> (g*g, 0.0));
    	D2Vec_Double K3(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K4(Sp, vector<double> (g*g, 0.0));

		D2Vec_Double Rho_M(tot_iter, vector<double> (Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		//double Rho_M[tot_iter][Sp2] ={0.0}; //Stores <rho>x values per time step for a single replicate.
		D2Vec_Double Rho_Mov(tot_iter, vector<double> (SpB4, 0.0)); //Stores averaged gamma and velocity values per time step for a single replicate.

		double Rhox_avg[Sp] = {0.0}; double Rhox_num[Sp] = {0.0}; //Stores spatial average and number of occupied sites at given t.

		//init_randbistableframe(Rho_dt, g*g, a, a_c,  perc, chigh, clow); // Returns a frame with random speckles of high and low density.

		///** EXPRT-TK INITIALISATIONS 
		// MFT PERTURBATION BASED FRAME INITIALISATION
		#if defined(INIT) && INIT == 0
			// HOMOGENEOUS FRAME INITIALISATION
			init_exprtk_homogenousMFTframe(Rho_dt, g*g, a, a_c, clow); // Returns a frame with random speckles of high and low density.
		#elif defined(INIT) && INIT == 2
			// BURN-IN FRAME INITIALISATION (EXPRTK)

			init_burnin_wrapper(Rho_dt, a, a_c, dP, perc,  g, j, clow);

			vector<double> temp_vec= {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species 's'
			double init_veg_num = occupied_sites_of_vector(temp_vec, g*g); //Finds number of occupied at given t.
			vector<double>().swap(temp_vec); //Flush temp_vec out of memory.
	
			if(init_veg_num == 0.0)
			{
				stringstream m1_1;
				m1_1 << "RUN-TIME WARNING: ZERO Active BURN-IN veg sites for Thread Rank:\t " << omp_get_thread_num() 
				<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << " CONSIDER RE-RUNNING OVER A DIFFERENT RANGE.\n"; 
				cout << m1_1.str(); cerr << m1_1.str();
			}
		#else
			// DEFAULT FRAME INITIALISATION (RANDOM MFT BASED SPECKLES)
			init_exprtk_randbiMFTframe(Rho_dt, g*g, a, a_c, dP, perc, clow); // Returns a frame with random speckles of high and low density.
		#endif
		//*/
		/** OLD FRAME INITIALISATION FOR 3 SPB.
		#if SPB == 3
		// NOTE: clow[] is used to store the fractional change from MFT steady state for some species, corresponds to c_spread[] in init_exprtk_randbiMFTframe(...).
		init_exprtk_randbiMFTframe_OLD(Rho_dt, g*g, a, a_c, dP, perc, clow); // Returns a frame with random speckles of high and low density.
		#endif
		// */

		/**  GAUSSIAN FRAME INITIALISATION  
		// Initialise vector <double> amp to elements of clow[].
		vector <double> amp(Sp, 500.0);  //Setting amplitude of gaussian distributions of vegetation to 500.
		for (int s=0; s< Sp; s++)
			amp[s] = clow[s+ Sp];
		
		vector <double> sd{g/8.0, g/16.0}; //Setting standard deviation of gaussian distributions of vegetation to g/8 and g/16.

		init_gaussframe(Rho_dt, g*g, sd, amp); 
		//*/
	
		// The first SpB species stored in Rho_dt are to be initialised on frame as Gaussian distributions.
		// Species 0 is the vegetation, Species 1 is the grazer and Species 2 is the predator.
		// Species 0 should be centered near the top right corner of the grid, Species 1 near the bottom right corner and Species 2 near the bottom left corner. */
		
		/** // GRADIENT FRAME INITIALISATION
		#if defined(LOCUST_WIND)
		// Define map for gradient initialisation of grazer (indexed 1) and predator (indexed 2) species.
		
		std::map<int, std::pair<double, double>> init_grad_map;
		init_grad_map[1] = {0.0, 1.0}; //init_grad_map[2] = {0.0, 0.0}; //Initialise to 0.0.

		// Initialise Rho_dt with gradient for grazer and predator species.
		init_gradientframe(Rho_dt, g*g, init_grad_map, "R"); //Initialise Rho_dt with gradient for Grazer from left to right.

		#endif
		//*/

		/** // BURN-IN FRAME INITIALISATION (OLD)
		
		stringstream  rain, rep_j, grid; 
		rain  << a; rep_j << int(j%2); grid << g;
    	string initcsv_filename= "Burning_In_Frames/FRAME_RAND_Basic_P_c_DP_G_" + grid.str() +"_T_91201_dt_0.2_alph_0.00833333_cgmax_0.0208333_a_" 
		+ rain.str() +"_dx_0.1_R_"  + rep_j.str() + ".csv";
		const vector<int> initcsv_columns = {2, 3, 4}; //Columns to be read from csv file.
		D2Vec_Double const_ind_val = {{1, dP/1000.0 }, {2, dP/5000.0}}; 
		init_csvconstframe(Rho_dt, const_ind_val, initcsv_filename, initcsv_columns, g*g); //Initialise Rho_dt with burn-in frame.

		// **/
		
		//*/

		for(int i=0; i < g*g; i++)
		{
			for(int s=0; s < Sp; s++)
			{
				DRho[s][i] = Rho_dt[s][i]; //Assigning initial conditions
				Rho_tsar[s][i] =Rho_dt[s][i];
			}
			//for(int s=1; s < SpB; s++)
			//   gamma[s][i] = 1.0; // Set gamma to 0.0 for debugging.
		}
		
		stringstream m1_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_2 << "Initial Conditions:\t Per:\t" << perc << "\t C_Lo[0], C_Lo[1], C_Lo[2]:\t " 
		<< clow[0] << "," << clow[1] << "," << clow[2] << ",\t R: " << a << "\n"; //cout << m1.str();
		if(omp_get_thread_num()== 1)
			cout << m1_2.str();  
		errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close(); 
		

		poisson_distribution<int> poisson; gamma_distribution <double> gamma_distr; 
		uniform_real_distribution<double> unif(0.0, 1.0); normal_distribution<double> norm(0.0, 1.0);

		//Define thrID as a random number between 0 and 10000, iff j =0.
		thrID = (j == 0) ? int(round(unif(rng)*10000)) : thrID; //Random number between 0 and 1000.
		double t=0; int index = 0;  //Initialise t
		int frame_index=0; int po=1; int so =1; int lo=1; int counter =0;
		vector <std::pair<int, int>> dtV_counter(SpB, {0, 0}); 
		//Keeps track of dtV[s] and counter%dtV[s] for each species s. Used to determine if species s is updated at time t.
		// Note counter := t/dt, so counter%dtV[s] = 0 implies t%dtV[s] = 0.

		for (int s=0; s< SpB; s++){	
			if(dtV[s] == 1)
				dtV_counter[s] = {1, 0}; //If dtV[s] = 1, then update species s at every time step.
			else
			{	dtV_counter[s] = {dtV[s], 0};	} //else update species s at every dtV[s] time step.
		}

		while( t < t_max + dt)
		{
			//Update dtV_counter for each species s.
			for (int s=1; s< SpB; s++)
			{
				if(dtV[s] != 1)
					dtV_counter[s].second = (counter)%dtV_counter[s].first;
			}

			// Basic Book-keeping below. Prints prelims out to file at given time points. 
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

				//Next store spatial average of gamma and velocity for each species at given t.
				#if SPB > 1
				double mov_avg, mov_var;
				for(int s=1; s< SpB; s++)
				{
					vector <double> temp= {gamma[s].begin(),gamma[s].end()}; //gamma for species 's'
					mov_avg = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
					mov_var = variance_of_vector(temp, g*g); //Finds spatial average of densities at given t.
					//rhox_num = occupied_sites_of_vector(temp, g*g); //Finds number of occupied at given t.
					Rho_Mov[index][4*s] = mov_avg; Rho_Mov[index][4*s +1] = mov_var; //Stores spatial average of gamma for species s at given t.
					vector<double>().swap(temp); //Flush temp out of memory.0

					// v_eff is a 2D vector of pairs of dim: SpBx(L*l), called v. Stores the velocity of each species at each site.
					// Compute the spatial average of the velocity in the x and y directions for each species.
					vector <double> temp_vx; //v_eff for species 's'
					vector <double> temp_vy; //v_eff for species 's'
					for (const auto& vel : v_eff[s]) 
					{
						temp_vx.push_back(vel.first); temp_vy.push_back(vel.second);
					}
					mov_avg = mean_of_vector(temp_vx, g*g); //Finds spatial average of densities at given t.
					Rho_Mov[index][4*s +2] = mov_avg; //Stores spatial average of gamma for species s at given t.
					
					mov_avg = mean_of_vector(temp_vy, g*g); //Finds spatial average of densities at given t.
					Rho_Mov[index][4*s +3] = mov_avg; //Stores spatial average of gamma for species s at given t.

					vector<double>().swap(temp_vx); vector<double>().swap(temp_vy); //Flush temp out of memory.
				}
				#endif

				vector <double> temp_alt= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
				rhox_num = occupied_sites_of_vector(temp_alt, g*g); //Finds number of occupied at given t.
				vector<double>().swap(temp_alt); //Flush temp out of memory.0

				// Saving frames to file at given time points.
				if( index == int(tot_iter*0.6) || index == int(tot_iter*0.85) || index == int(tot_iter*0.9) || index == int(tot_iter*0.95) || index == int(tot_iter-1))
				{
					int break_flag = save_prelimfileswrapper(index, tot_iter, j, thrID, t, dt, dx, Rho_dt, DRho, Rho_M, Rho_Mov, rho_rep_avg_var,
						t_meas, gamma, v_eff, t_max, a, c, gmax, alpha, rW, W0, D, v, K, sigma, a_st, a_end, a_c, A,H,E,M, 
						dP, r, g, Gstar /* =-1.*/, Vstar /* = -1.*/);
					if( break_flag == 1)
						break; //Finish negotiations as all basal species are dead.
				}

				//vector <double> temp_1= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
				//double rhox_DR = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.
				//vector<double>().swap(temp_1); //Flush temp out of memory.0

				//temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
				//double  rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
				//vector<double>().swap(temp_1); //Flush temp out of memory.
				
				index+=1;
				/**
				stringstream m4_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m4_1 << "AFTER SAVING FRAME, STATUS UPDATE AT TIME [t,i]\t" << t_meas[index] << " with # non-zero sites in VEG Rho_dt_0:  " << rhox_DR 
				<< " AND # Non-Zero in VEG DRho_0:  " << rhox_dt <<  "\n"; cout << m4_1.str();
				errout.open(thr, std::ios_base::app); errout << m4_1.str(); errout.close(); **/
				//int counter_nzero_alpha=0;
			}
			// BASIC BOOK-KEEPING. Saves frame to file at given time points.
			if(t == 0 || t >= frame_tmeas[frame_index] -dt/2.0 && t < frame_tmeas[frame_index] +dt/2.0)
			{
				// Saving frames to file at given time points.
				if(frame_index >= frame_tot_iter -10 &&  frame_index <= frame_tot_iter-1  ||  t >= 60000 && t <= 150000 
					|| t >= 50 && t <= 15000 || t== 0)
				{
					save_framefileswrapper(frame_index, frame_tot_iter, j, thrID, t, dt, dx, Rho_dt, DRho, rho_rep_avg_var,
					frame_tmeas, gamma, v_eff, t_max, a, c, gmax, alpha, rW, W0, D, v, K, sigma, a_st, a_end, a_c, A,H,E,M, 
					dP, r, g, Gstar /* =-1.*/, Vstar /* = -1.*/);
				}
				frame_index+=1;
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
					/**
					stringstream m6_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
					m6_1 << "DORNIC VEG RHO* " << Rho_dt[s][i] << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
					<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i  << "\t and POISSON ENTRY:\t" << ziggy << "\t and GAMMA ENTRY:\t" << gru << "\n"; 
					cout << m6_1.str(); errout.open(thr, std::ios_base::app); errout << m6_1.str(); errout.close(); **/
					if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "VEG HIT ROCK BOTTOM WITH: Rho*"<< s <<"[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i << 
						"\t and GAMMA ENTRY:\t" << gru << "\n"; cout << m6.str(); cerr << m6.str();
					}
					//Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
				}
			} // End of Vegetation Integration

			
			// Book-keeping for determining gamma_i for higher order species. Calculating Rho averages per species at each time step.
			// Compute ONLY IF SPB > 1
			#if SPB > 1
			for(int s=0; s < SpB; s++)
			{
				vector <double> temp= {DRho[s].begin(),DRho[s].end()}; //Rho_dt for species '0'
				Rhox_avg[s] = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
				vector<double>().swap(temp); //Flush temp out of memory.0
			}
			vector <double> temp_veg= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
			double rhox_num_veg = occupied_sites_of_vector(temp_veg, g*g); //Finds number of occupied sites at given t.
			vector<double>().swap(temp_veg); //Flush temp out of memory.
			//double nR_fac = 1;
			double nR_fac = 1 - rhox_num_veg/(g*g); //Factor to reduce the number of neighbours for gamma_i estimation
			if (nR_fac < 0.35)
			{	nR_fac = 0.35; }
			int nR_Perp_length = int(nR_fac*nR_fac*origin_Neighbourhood.size());
			nR_Perp_length < 1 ? nR_Perp_length = 1 : nR_Perp_length = nR_Perp_length; 
			// At least one nearest neighbor (the site itself) is considered.
			#endif

			
			// NOTING THAT THE DORNIC INTEGRATION IS DONE FOR VEGETATION ONLY, WITH MEAN AND OCCUPANCY FOR VEGETATION SHOWED
			#if defined(DEBUG)
			stringstream m6_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
			vector <double> temp_veg_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
			double rhox_mean_veg = mean_of_vector(temp_veg_1, g*g); //Finds spatial average of densities at given t.
			double rhox_num_veg1 = occupied_sites_of_vector(temp_veg_1, g*g); //Finds number of occupied sites at given t.
			vector<double>().swap(temp_veg_1); //Flush temp out of memory.
			m6_2 << "DORNIC INTEGRATION FOR VEG DONE,  AT TIME [t, thr, a]\t" << t << " , " << omp_get_thread_num() << " , " << a 
			<< " with # non-zero sites in VEG Rho_dt_0:  " << rhox_num_veg1 << " AND MEAN VEG DENSITY:  " << rhox_mean_veg << "\n"; cout << m6_2.str();
			cerr << m6_2.str();
			#endif
			//*/
			
			#if defined(BARRACUDA) && SPB > 1
			//Transfer Rho_dt, gamma and v_eff to device memory.
			
			for (int s = 0; s < 1; s++)
			{	cudaMemcpy(d_Rho_dt, Rho_dt[s].data(), g * g * sizeof(double), cudaMemcpyHostToDevice);		}
			
			int rho_offset = 0;
			for (int s = 1; s < SpB; s++)
			{
				cudaMemcpy(d_Rho_dt + rho_offset + g * g, Rho_dt[s].data(), g * g * sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(d_gamma + rho_offset, gamma[s].data(), g * g * sizeof(double), cudaMemcpyHostToDevice);
				//transfer_castedVec_to_device(gamma[s].data(), d_gamma + rho_offset, g * g);
				cudaMemcpy(&d_v_eff[(s - 1) * g * g], v_eff[s].data(), g * g * sizeof(std::pair<double, double>), cudaMemcpyHostToDevice);
				rho_offset += g * g;
			}
			// Reset Rho_tsar on device memory to 0.0
			//cudaMemset(d_Rho_tsar, 0.0, Sp_NV * g * g * sizeof(double));

			if (counter % 50000 == 1)
			{
				stringstream m6_2;     //To make cout thread-safe as well as non-garbled due to race conditions
				m6_2 << "CUDA TRANSFER COMPLETE FOR HIGHER ORDER SPECIES AT TIME [t, thr]\t" << t << " , " << omp_get_thread_num() << "\n";
				cout << m6_2.str(); cerr << m6_2.str();
				//if (omp_get_thread_num() == 2)
				//	reportPerceptionArrays_Interface(d_origin_Neighbourhood, d_rho_inverse, d_eff_sizes, nR_Perp_length, r_max_effective, g);
			}

			// Launch CUDA kernel to compute gamma and velocity.
			#if SPB == 3
				calc_gamma_vel_NonRefugia_CUDA(d_origin_Neighbourhood, d_Rho_dt, d_gamma, d_v_eff,
				d_rho_inverse, d_eff_sizes, Rhox_avg, r_frac, dtV_counter, nR_Perp_length, r_max_effective, g);
			#endif
			// Copy results back to host memory.
			rho_offset = 0;
			for (int s = 1; s < SpB; s++)
			{
				cudaMemcpy(gamma[s].data(), d_gamma + rho_offset, g * g * sizeof(double), cudaMemcpyDeviceToHost);
				//transfer_castedVec_to_host(d_gamma + rho_offset, gamma[s].data(), g*g);
				cudaMemcpy(v_eff[s].data(), &d_v_eff[(s - 1) * g * g], g * g * sizeof(std::pair<double, double>), cudaMemcpyDeviceToHost);
				rho_offset += g * g;
			}
			#else
			//Compute gamma and velocity for higher order species using CPU implementation.
				#if SPB == 3
					//calc_gamma_3Sp_NonRefugia(origin_Neighbourhood, DRho, gamma, Rhox_avg, r_frac, dtV_counter, nR_fac, r_max_effective, g);
					calc_gamma_3Sp_NonRefugia_FFT(DRho, gamma, Rhox_avg, r_frac, dtV_counter, nR_fac, r_max_effective, g, fft_central_plan);
				#elif SPB == 2
					//calc_gamma_2Sp_NonRefugia(origin_Neighbourhood, DRho, gamma, Rhox_avg, r_frac, dtV_counter, nR_fac, r_max_effective, g);
					calc_gamma_2Sp_NonRefugia_FFT(DRho, gamma, Rhox_avg, r_frac, dtV_counter, nR_fac, r_max_effective, g, fft_central_plan);
				#endif 
				//Calculates gamma for each species at each site.
			#endif


			#if defined(DEBUG)
			if(counter%50000 == 0)
			{
				stringstream m6_3;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m6_3 << "GAMMA CALCULATION DONE AT TIME [t, thr, a, j]\t" << t << " , " << omp_get_thread_num() << " , " << a << " , " << j << "\n"; 
				cout << m6_3.str(); cerr << m6_3.str();
			}
			#endif

			/** // BLOCK FOR SAVING GAMMA FRAMES
			if(t == 0 || t >= t_meas[index-1] -dt/2.0 && t < t_meas[index-1] +dt/2.0)// && index >= tot_iter -9 &&  index <= tot_iter)
			{
				// Saving gamma values to file at given time points.
				stringstream L, tm ,d3, p1, rini, a1, a2, Dm1, Dm2, dix, dimitri, geq, veq;

  				L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  				rini << j; Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2];
				a1 << a_st; a2  << a_end; dimitri  << dP; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;

				string filenamePattern = gamma_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str()  + "_D1_"+ Dm1.str() + "_D2_"+ Dm2.str() + "_dx_"+ dix.str() + "_R_";
				
				string parendir = "";
				// Creating a file instance called output to store output data as CSV.
				if(Gstar != -1)
					parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str();
				else if(Vstar != -1)
					parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str();
				else
					parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();

				string gammaheader = "a_c,  x,  Gam0(x; t), Gam1(x; t), Gam2(x; t) \n"; //Header for output frame.
				save_frame(gamma, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
				
			}
			// */

			#if defined(DEBUG) // STATUS UPDATE ABOUT VEGETATION
			///** STATUS UPDATE
			if(counter%50000 == 0)
			{

				vector <double>  temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
				double rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, thr]\t" << t << " , " << omp_get_thread_num()  << "\t VEG DORNIC DONE with non-zero elements in DRho_0:  " 
				<< rhox_num_veg1 << "\t and non-zero elements in Rho_dt_0:  " << rhox_dt << "\t and INDEX VAL, T_MEAS:\t" << index <<  " , " << t_meas[index] << "\n"; cout << m5_1.str();
				errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			//*/
			#endif
			// Now for higher order species.
			for(int s=1; s< SpB; s++)
            {
				if(Rhox_avg[s] == 0.0)
					continue;
				//If the average density of a species is 0, then skip the Dornic integration for that species.
				
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					//Finally check if gamma_i > 1, if so set it to 1.
					if(gamma[s][i] > 1)
					{	gamma[s][i] = 1; }

					/** COUNTER GAMMA MESSAGE
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j 
						<< "\t GAMMA FOUND WITH GAMMA_I:" << gamma[s][i] << "\n"; cout << m5_2.str();
					} */
					
					if( gamma[s][i] < 0.0)
					{
						stringstream m5_3;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_3 << "GAMMA_I FALLS BELOW 0,  AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j  
						<< "\t WITH GAMMA_I:" << gamma[s][i] << " for species: " << s << "\n"; cout << m5_3.str(); cerr << m5_3.str(); 
					}
					#if defined(DEBUG) // Check if gamma_i is NaN or Inf, if so, SAVE GAMMA & exit the program.
					if(isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
					{
						stringstream m5_3;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_3 << "GAMMA_I BLOWS THE LID,  AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j  
						<< "\t WITH GAMMA_I:" << gamma[s][i] << " for species: " << s << " with Rho_avg[0]: " << Rhox_avg[0] << " and Rho_avg[1]: " << Rhox_avg[1] << " and Rho_avg[2]: " << Rhox_avg[2] 
					  	<< " \n and Rho_t[0][i]: " << Rho_dt[0][i] << " and Rho_t[1][i]: " << Rho_dt[1][i] << " and Rho_t[2][i]: " << Rho_dt[2][i]<< "\n"; cout << m5_3.str();
						cerr << m5_3.str();

						stringstream L, tm ,d3, p1, rini, a1, a2, Dm1, Dm2, dix, dimitri, geq, veq;

  						L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  						rini << j; Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2];
						a1 << a_st; a2  << a_end; dimitri  << dP; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;

						string ErrGamfilenamePattern = "/ERROR_GAMMA_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
						+ "_a_" + p1.str()  + "_D1_"+ Dm1.str() + "_D2_"+ Dm2.str() + "_dx_"+ dix.str() + "_R_";
						string parendir = "";

						// Creating a file instance called output to store output data as CSV.
						if(Gstar != -1)
							parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str();
						else if(Vstar != -1)
							parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str();
						else
							parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();

						string gammaheader = "a_c,  x,  Gam0(x; t), Gam1(x; t), Gam2(x; t) \n"; //Header for output frame.
						save_frame(Rho_dt, parendir, ErrGamfilenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
						exit(4);
						
					}
					#endif

					//Next sampling the advection vector from v[s].
					//unif = uniform_real_distribution<double>(0.0, 1.0);

					#if defined(LOCUST_WIND) 
					// Set grazer advection vector to max of absolute value of wind speed and v[s].
					// Sampling the advection vector from v[s], every dtV[s] time steps.
					if(dtV[s] == 1 || dtV_counter[s].second == 0)
					{
						
						if (s == 1 && (wind_speed.first != 0.0 || wind_speed.second != 0.0)) 
						{
    						// Use wind speed if it has at least one non-zero component
    						v_eff[s][i].first = wind_speed.first;	v_eff[s][i].second = wind_speed.second;
						}
						// Regular sampling of advection vector for other species
						else
						{
							double ran = unif(rng);
							v_eff[s][i].first = v[s]*cos(2*PI*ran); v_eff[s][i].second = v[s]*sin(2*PI*ran); //Advection vector (vx, vy).
						}
						
					}
					#else
						double ran = unif(rng);
						//double vx = v[s]*cos(2*PI*ran); double vy = v[s]*sin(2*PI*ran); //Advection vector.
						v_eff[s][i].first = v[s]*cos(2*PI*ran); v_eff[s][i].second = v[s]*sin(2*PI*ran); //Advection vector (vx, vy).
					#endif


					//double vx = v_eff[s][i].first; double vy = v_eff[s][i].second; //Advection vector.
					double vx_abs = v_eff[s][i].first*sgn(v_eff[s][i].first); double vy_abs = v_eff[s][i].second*sgn(v_eff[s][i].second); //Absolute value of advection vector.

					//RECALL: dx1_2 = 1/(dx2); dxb2 = dx/2.0;

					//diff_eff[s].first = (D[s] -(dxb2)*(vx_abs*(1 - dt*dx1*vx_abs)))*dx1_2;
					//diff_eff[s].second = (D[s] -(dxb2)*(vy_abs*(1 - dt*dx1*vy_abs)))*dx1_2;
					double gamma_prime_si = (1 - gamma[s][i]); // Probability of advection.
					///** CORRECT DsC VERSION 
					diff_eff[s].first = (D[s] -(dxb2)*(gamma_prime_si*vx_abs*(1 - dt*dx1*gamma_prime_si*vx_abs)))*dx1_2;
					diff_eff[s].second = (D[s] -(dxb2)*(gamma_prime_si*vy_abs*(1 - dt*dx1*gamma_prime_si*vy_abs)))*dx1_2;
					//*/
					//diff_coefficient[s] = (D[s] -(dx/2.0)*(v[s]*(1 - dt*vdx[s])))/(dx2);
					// Advection leads to excess diffusion, hence the correction term.

					if(diff_eff[s].first < diff_coefficient[0])
						diff_eff[s].first = diff_coefficient[0];
					
					if(diff_eff[s].second < diff_coefficient[0])
						diff_eff[s].second = diff_coefficient[0];
					// If the effective diffusion coefficient is less than D[0]/dx2, set it to D[0]/dx2.

					#if defined(DEBUG)
					if(counter%50000 == 0 && counter == i% (16*g))
					{ 
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j 
						<< "\t DORNIC VALUE OF SP " << s << " UPDATED WITH RHO*(I,t): " << DRho[s][i] 
						<< " and (vx, vy):  [" <<  v_eff[s][i].first << " , " << v_eff[s][i].second << "] \n"; cout << m5_2.str();
						cerr << m5_2.str();
					}
					#endif
					

					// RECALL:	nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 
					// NOTE: dx1 is the inverse of dx. sgn_index(vx) returns 0 if vx is positive, else 1.
					// Recall nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1;
					// ALSO NOTE: diff_coefficient[s] = D[s]/(dx*dx); diff_eff[s].first = D_X_eff[s]/(dx*dx); diff_eff[s].second = D_Y_eff[s]/(dx*dx);
					
					// Standard diffusion term.
					//double alpha_i = gamma[s][i]*diff_coefficient[s]*(1*(DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
					//				DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]])) 
					//				+dx1*(1-gamma[s][i])*(vx_abs*DRho[s][nR2[i][1][sgn_index(vx)]]+ vy_abs*DRho[s][nR2[i][0][sgn_index(vy)]]);

					// Advection corrected diffusion term.
					double alpha_i = gamma[s][i]*(diff_eff[s].second*(DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]]) 
												+ diff_eff[s].first*(DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]])) 
									+dx1*(gamma_prime_si)*(vx_abs*DRho[s][nR2[i][1][sgn_index(v_eff[s][i].first)]]+ vy_abs*DRho[s][nR2[i][0][sgn_index(v_eff[s][i].second)]]);
					
					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					//Beta depends on the gamma_i value at each lattice site for each species, hence needs to be updated at each iteration.
					// Standard diffusion term.
					//beta[s] = -M[s] - 4*gamma[s][i]*diff_coefficient[s] - (1-gamma[s][i])*(vx_abs + vy_abs)*dx1; // NOTE: diff_coefficient[s] = D[s]/(dx*dx)
					// Advection corrected diffusion term.
					beta[s] = -M[s] - 2*gamma[s][i]*(diff_eff[s].first + diff_eff[s].second) - (1-gamma[s][i])*(vx_abs + vy_abs)*dx1;

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
						"\t ,GAMMA ENTRY:\t" << gru << "\t (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]"
						"\t and (VX, VY):  [" << v_eff[s][i].first << " , " << v_eff[s][i].second <<  " ]" <<"\n";
						cout << m6.str(); cerr << m6.str();
						exit(3);
					}
					//Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
					///** COUNTER STATUS UPDATE
					#if defined(DEBUG)
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j 
						<< "\t DORNIC VALUE OF SP " << s << " UPDATED WITH RHO*(I,t): " << Rho_dt[s][i] 
						<< " and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" << "\n"; cout << m5_2.str();
						errout.open(thr, std::ios_base::app); errout << m5_2.str(); errout.close();
					}
					#endif
					//*/
				}
			}
			/** // Finally, update Rho_tsar with the new values of Rho_dt for the abiotic species.
			for(int s=SpB; s< Sp; s++)
			{
				for(int i=0;i<Rho_dt[0].size();i++)
					Rho_tsar[s][i] = Rho_dt[s][i];		
					// This is the rho* value (refer to Dornic et 2005)
			} //**/

			
			
			if(counter%50000 == 0)
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j << "\t GRAZER + PRED DORNIC DONE, FOR " << prefix + "\n"; 
				cout << m5_1.str(); cerr << m5_1.str();
			}
			
			// NOTE: diff_coefficient[s] = D[s]/(dx*dx)
			// Finally RK4 integration of the remaining terms.
			#if SPB == 3
				RK4_Integrate_Stochastic_MultiSp(Rho_dt, Rho_tsar, K1, K2, K3, K4, nR2, a, c, gmax, alpha, rW, W0, diff_coefficient, K, A,H,E, t, dt, dx, g);
			#elif SPB == 2
				RK4_Integrate_Stochastic_2Sp(Rho_dt, Rho_tsar, K1, K2, K3, K4, nR2, a, c, gmax, alpha, rW, W0, diff_coefficient, K, A,H,E, t, dt, dx, g);
			#elif SPB == 1
				RK4_Integrate_Stochastic_1Sp(Rho_dt, Rho_tsar, K1, K2, K3, K4, nR2, a, c, gmax, alpha, rW, W0, diff_coefficient, K, t, dt, dx, g);
			#endif
			//RK4_Integrate_Stochastic_MultiSp(Rho_dt, Rho_tsar, K1, K2, K3, K4, nR2, a, c, gmax, alpha, rW, W0, diff_coefficient, K, A,H,E, t, dt, dx, g);
			
			if(counter%50000 == 0)
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num()  << " , " << j << "\t RK4 INTEGRATION OF REMAINDER TERMS DONE." 
				<< "\n"; cout << m5_1.str(); cerr << m5_1.str();
			}
			
			for(int s=0; s< Sp; s++)
			{
					
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					// CHECK FOR NAN AND INF
					if( (Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true) and so==1)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO "<<s<<"  FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i] << "\t and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" <<  "\n"; 
						cout << m6.str(); cerr << m6.str(); lo == -1;
					}
					// */ 

					if( isnan(Rho_dt[s][i] == true))
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO "<<s<<"  NAAN PARATHA:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<< "\t and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" <<  "\n"; 
						cout << m6.str(); cerr << m6.str();
						lo == -1;
					}
					// Set all lattice sites that are less than 10^-12 to 0.0.
					if(Rho_dt[s][i] < epsilon && Rho_dt[s][i] >= -1.0e-80)
						Rho_dt[s][i] = 0.0;	
					DRho[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
					//Rho_tsar[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
				}
			}

			t+=dt; //Update timestep
			counter+=1; //Update counter

			if(counter%40000 ==1)
			{
				//size_t currentSize = getCurrentRSS( ); //Check Memory usage.
				//size_t peakSize    = getPeakRSS( ); //Find memory usage and time taken for every 50000 iterations.
				auto end_t = high_resolution_clock::now();
				auto elapsed_min = duration_cast<minutes>(end_t-start_t);
				// Reset the clock
				start_t = high_resolution_clock::now();
				stringstream m7;
				m7 << "STATUS UPDATE AT [t, a, thr, j]\t " << t << " , "  << a << " , " << omp_get_thread_num() << " , " << j 
				//<< " ,with current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) 
				<< " and time taken for " << counter << " iterations: " << elapsed_min.count() << " min." <<  "\n"; cout << m7.str();
				errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();
				//Report syst every 50000 iterations.
			}

			if(lo == -1)
			{
				stringstream m6_1;     //To make cout thread-safe as well as non-garbled
				m6_1 << "Saving error frames and exiting program .... \n"; cout << m6_1.str(); cerr << m6_1.str();

				// SAVING ERROR Rho_dt, DRho, and gamma vectors TO FILE
				stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, Dm2, alph, w0t, aij, hij, dix, dimitri, sig0, sig1, geq, veq;

				L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
				rini << j; Dm0 << D[0]*pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2]; 
				a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1]; dimitri  << dP; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;
				//gm << setprecision(3) << gmax; w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; 
				// Three replicates are over.

				string frame_parendir = "";

				// Creating a file instance called output to store output data as CSV.
				if(Gstar != -1)
					frame_parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str();
				else if(Vstar != -1)
					frame_parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str();
				else
					frame_parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();

				// SAVING RHO_DT FRAME
				string filenamePattern_Rhodt = "/ERROR_RHO_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
				+ "_a_" + p1.str()  + "_dx_"+ dix.str() + "_R_";
				//Next, designate the naive filename
				string filename = filenamePattern_Rhodt + rini.str() + ".csv";
				save_frame(Rho_dt, frame_parendir, filenamePattern_Rhodt, a, a_st, a_end, t, dt, dx, dP, j, g);

				// SAVING DRHO FRAME
				string filenamePattern_DRho = "/ERROR_DRHO_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
				+ "_a_" + p1.str()  + "_dx_"+ dix.str() + "_R_";
				//Next, designate the naive filename
				filename = filenamePattern_DRho + rini.str() + ".csv";
				save_frame(DRho, frame_parendir, filenamePattern_DRho, a, a_st, a_end, t, dt, dx, dP, j, g);

				///**  SAVING GAMMA FRAME
				string filenamePattern_Gamma = "/ERROR_GAMMA_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
				+ "_a_" + p1.str()  + "_dx_"+ dix.str() + "_R_";
				//Next, designate the naive filename
				filename = filenamePattern_Gamma + rini.str() + ".csv";
					
				//string gammaheader = "a_c,  x,  Gam0(x; t), Gam1(x; t) \n"; //Header for output frame.
				save_frame(gamma, frame_parendir, filenamePattern_Gamma, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
				//*/
				
				// SAVING PRELIMINARY FRAME
				// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
				// and save to file.
				D2Vec_Double rho_rep_avg_var_temp(index, vector<double> (Sp4_1, 0.0)); 
				
				//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
				//std::copy(rho_rep_avg_var.begin(), rho_rep_avg_var.begin() + index, rho_rep_avg_var_temp.begin());
				var_mean_incremental_surv_runs(rho_rep_avg_var_temp, Rho_M, index, 0);

				string filenamePattern_Prelim = "/ERROR_SERIES_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
				"_D1_"+ Dm1.str() + "_R_";
				//Next, designate the naive filename

				string prelim_parendir= "";
				if(Gstar != -1)
					prelim_parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str() + "/TimeSeries";
				else if(Vstar != -1)
					prelim_parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str() + "/TimeSeries";
				else
					prelim_parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/TimeSeries";

				save_prelimframe(rho_rep_avg_var_temp, prelim_parendir, filenamePattern_Prelim, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);
				exit(3);
			}
		} // End of time loop.

		vector<vector<double>>().swap(gamma); vector<vector<double>>().swap(K1); 
		vector<vector<double>>().swap(K2); vector<vector<double>>().swap(K3); vector<vector<double>>().swap(K4);
		vector<vector <double>>().swap(Rho_dt); vector<vector <double>>().swap(DRho); vector<vector <double>>().swap(Rho_tsar);
		// Freeing up memory.
		if(int(omp_get_thread_num()) == 5)
		{
			size_t currentSize = getCurrentRSS( ); //Check Memory usage.
			size_t peakSize    = getPeakRSS( );
			stringstream m7;
			m7 << "In replicate " << j << ", current size in MB: " << currentSize/(1024.0*1024.0) 
			<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << "\n"; cout << m7.str();
			errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();
		}

		if(j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "INITIAL UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << tot_iter << "\n"; cout << m7.str();
			for(int i =0; i< tot_iter; i++)
			{	
				for(int s =0; s< SpB; s++)
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
			stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, rini_prev, Dm, cgm,  sig0, geq, veq, jID;
			
			a1 << a_st; a2 << a_end;
  			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(5) << a; dimitri << dP;
  			rini << j; Dm << setprecision(4) << D[2]; cgm << c*gmax; sig0 << sigma[0]; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;

			string parendir ="";
			if(Gstar != -1)
				parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str();
			else if(Vstar != -1)
				parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str() + "/TimeSeries";
			else
				parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();

			double ran_jid = (unif(rng)*1000.0)/1000.0; jID << ran_jid; // Random number between 0 and 1.
			
			string filenamePattern = prelim_prefix + jID.str() +"_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
			"_D2_"+ Dm.str() + "_R_";

			save_prelimframe(rho_rep_avg_var, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, true);
		}

	} // End of r loop.

	// Free up all device memory.
	#if defined(BARRACUDA)
	cudaFree(d_Rho_dt); cudaFree(d_Rho_tsar); cudaFree(d_gamma); cudaFree(d_v_eff); 
	cudaFree(d_origin_Neighbourhood); cudaFree(d_rho_inverse); cudaFree(d_eff_sizes);
	#endif

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


void first_order_critical_exp_delta_stochastic_MultiSp(int div, double t_max, double a_start, double a_end, double a_c,  double c, double gmax, double alpha,
	double rW, double W0,  double (&D)[Sp], double (&v)[SpB], double (&K)[3], double sigma[], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double pR[], int (&dtV)[SpB], double clo[],
	double dt, double dx, double dP, int r,  int g, double Gstar /*= -1.*/, double Vstar /*= -1.*/)
{
	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//double mean[Sp] = {1.0, 0.05}; double sd[Sp] = {0.25, 0.0125};
	//init_randframe(Rho_0, Sp,  g*g); //Returns Rho_0 with a full initial frame filled with 0.2.

	//set_density_dependentmortality(M[2]); //Sets the density dependent mortality rate.

	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "NOSTRA" <<endl;
  	// The pspace to iterate over.
	// Store linspace values in frame_tmeas.
	

    vector <double> t_measure = logarithm10_time_bins(t_max, 1*dt);
	//vector <double> t_measure = linspace(40, 1200, 30);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  	// Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

	//vector <double> t_measure = linspace(40.0, 1000.0, 25); // Computes and returns linearly distributed points from t= 20 to 240 (the latter rounded down to 1 decimal place
	t_measure.insert(t_measure.begin(), 0.0);

	

	//Remove last element of t_measure
	//t_measure.pop_back();
	
	// Remove all elements of t_measure that lie between 750 and 1750.
	//t_measure.erase(std::remove_if(t_measure.begin(), t_measure.end(), [](double x){return (x > 600.0 && x < 1500.0);}), t_measure.end());


	//vector <double> t_linearwindow = linspace(600, 1500, 19); // Computes and returns linearly distributed points from t= 20 to 240 (the latter rounded down to 1 decimal place

	// Insert the linear window to t_measure vector after the element in t_measure that is just less than the first element in t_linearwindow.
	//auto it = std::upper_bound(t_measure.begin(), t_measure.end(), t_linearwindow[0]);
	//t_measure.insert(it, t_linearwindow.begin(), t_linearwindow.end());

	// Sort the t_measure vector in ascending order and remove duplicates (if any).
	sort( t_measure.begin(), t_measure.end() );

	t_measure.erase( unique( t_measure.begin(), t_measure.end() ), t_measure.end() );

	// Copy t_measure values to frame_tmeas [PRELIM AND FRAME TMEAS ARE THE SAME]
	frame_tmeas = t_measure;

	// Use switchsort_and_bait() to swap values in frame_tmeas b/w [50, 400] with linspace(50, 400, 8) 
	// and [400, 1200] with linspace(400, 1200, 9).
	frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 50, 400, 8, "linspace", true);
	frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 400, 1200, 9, "linspace", true);
	frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 2000, 10000, 4, "linspace", true);

	// For correlation measurement, add linear window from 20 to 6000, with 300 points (spaced 20 apart).
	//frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 20, 6000, 300, "linspace", true);

	// For spreading experiments, add linear window from 40 to 1200, with 30 points (spaced 40 apart).
	//frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 40, 1200, 30, "linspace", true);


	if(t_measure[t_measure.size()-1] > 90000)
	{

		// Remove all elements of frame_tmeas that lie between 50000 and 100000.
		//frame_tmeas.erase(std::remove_if(frame_tmeas.begin(), frame_tmeas.end(), 
		//[](double x){return (x > 50000.0 && x < 100000.0);}), frame_tmeas.end());
		// Create linear window from 50000 to 100000, with 26 points (spaced 2500 apart).
		vector <double> t_frame_linearwindow = linspace(80000, 120000, 11); 
		// Rounded down to 1 decimal place.
		// Insert the linear window to frame_tmeas vector after the element in frame_tmeas that is just less than the first element in t_frame_linearwindow.
		auto it = std::upper_bound(frame_tmeas.begin(), frame_tmeas.end(), t_frame_linearwindow[0]);
		frame_tmeas.insert(it, t_frame_linearwindow.begin(), t_frame_linearwindow.end());
		// Sort the frame_tmeas vector in ascending order and remove duplicates (if any).
		sort( frame_tmeas.begin(), frame_tmeas.end() );

		// Create a linear window from 80000 to 150000, spaced 500 apart (141 points).
		vector <double> t_frame_linearwindow2 = linspace(80000, 120000, 81);
		// Insert into t_measure vector after the element in t_measure that is just less than the first element in t_frame_linearwindow2.
		auto it2 = std::upper_bound(t_measure.begin(), t_measure.end(), t_frame_linearwindow2[0]);
		t_measure.insert(it2, t_frame_linearwindow2.begin(), t_frame_linearwindow2.end());
		// Sort the t_measure vector in ascending order and remove duplicates (if any).
		sort( t_measure.begin(), t_measure.end() );	
	}
	// Remove all duplicate elements in frame_tmeas and t_measure.
	frame_tmeas.erase( unique( frame_tmeas.begin(), frame_tmeas.end() ), frame_tmeas.end() );
	t_measure.erase( unique( t_measure.begin(), t_measure.end() ), t_measure.end() );

  	cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  	for (int i=0; i< t_measure.size(); i++)
  	{
		cout << t_measure[i] << " ";
  	} 	cout << endl;

	cout << "Values in frame_tmeas (which is of total size " << frame_tmeas.size() << ") are :" <<endl;
  	for (int i=0; i< frame_tmeas.size(); i++)
  	{
		cout << frame_tmeas[i] << " ";
  	} 	cout << endl;

	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd); // Seed the generator
	std::uniform_int_distribution<int> unif_int(0, 30); //Use this to generate random numbers between 0 and 30.

	int r_max = (*max_element(pR, pR + SpB)) / dx; // Maxmimum perception radius across all species in units of lattice sites.
	cout << "Maximum Perception Radius in Lattice Sites is: " << r_max << endl;
	
  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = unif_int(rng)* microToSeconds;     //5 seconds
	// Random delay to stagger the start of the parallel threads.
	usleep(static_cast<int>(delay1)); //Sleep for delay1 microseconds
	//std::this_thread::sleep_for(std::chrono::microseconds(static_cast<int>(delay1))); //Sleep for delay1 microseconds
    cout<<"Delay 1 in progress... ("<<delay1/microToSeconds<<"s)"<<endl;
	std::vector<vector <double>> vec;
  	// Stores collated output from parallel method calls in proper ascending order of p values.

  	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );

	cout << "On initialisation, current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl;

	auto start = high_resolution_clock::now();
	//fftw_init_threads();
	
	// Create temporary workspace for kernel precomputation
	if(SpB > 1)
	{
		// Create temporary planner and then use it to precompute all FFT kernels for KHats
		FFTW3_CentralPlanner fiveyrplan(g); precompute_FFT_KHat_kernels(fiveyrplan, g, r_max);
	} //Auto destructor will clean up planner after this block.

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
	  FFTW3_CentralPlanner* fftw_central_planner = nullptr;
	
	  #pragma omp critical
	  {
		#if SPB > 1
		fftw_central_planner = new FFTW3_CentralPlanner(g); // Central planner for multi-species simulations.
		#else
		fftw_central_planner = new FFTW3_CentralPlanner(1); // Central planner is not necessary for single species, but kept for code uniformity.
		#endif
	  }

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
		rietkerk_Dornic_2D_MultiSp(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, rW, W0, D, v, K , sigma, a_start, a_end, a_c, A,H,E,M, pR, dtV, clo, dt, dx, dP, r, g, *fftw_central_planner, Gstar, Vstar);
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
	  // Finally clean up the FFTW3 Central Planner
	  delete fftw_central_planner;
  	} //End of Pragma.
  	cout << "Finished all parallel computations." << endl;

	// Finally clean up the shared precomputed FFTW3 kernel cache (kHat_kernel_FFT_cache inline std::map<int, fftw_complex*>)
	// Do this by iterating over all the keys in the map and freeing up the allocated fftw_complex* arrays.
	for (auto const& fft_cache_pair_R : kHat_kernel_FFT_cache)
		fftw_free(fft_cache_pair_R.second);
	
    //fftw_cleanup_threads();

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	//Store a uniformly distributed random number between 0 and 1 (rounded to 3 decimal places).
	double id = round(unif(rng)* 1000.0) / 1000.0; //Random number between 0 and 1.
	//Round to 3 decimal places.
	stringstream L, coco, tm ,d3, p1, p2, rini, Dm0, Dm1, aij, hij, cgm, alph, dix, dimitri, Sig0, geq, veq, ID;

	//double Gstar = M[2]/((E[2] -M[2]*H[1][2])*A[1][2]); // MFT estimate of Grazer density at coexistence.

  	L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << a_start; p2  << a_end; rini << r; 
	alph << alpha; cgm << c*gmax; coco << setprecision(4) << c;
  	Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; Sig0 << sigma[0]; dix << setprecision(2) << dx; 
	aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;
	ID << id;

	string savefile = stat_prefix + L.str()+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D1_"+ Dm1.str() + "_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str();
	if(Gstar != -1)
		savefile += "_Geq_" + geq.str() /* + "_ID_" + ID.str() */ + "_R_"+ rini.str() + ".csv";
	else if(Vstar != -1)
		savefile += "_Veq_" + veq.str() /* + "_ID_" + ID.str() */ + "_R_"+ rini.str() + ".csv";
	else
		savefile += /* + "_ID_" + ID.str() */ + "_R_"+ rini.str() + ".csv";

	ofstream output_1stdp;
  	// Creating a file instance called output to store output data as CSV.
	output_1stdp.open(savefile);

	cout << "Save file name: \n";
	cout << savefile <<endl;

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

	string header = prelimheader;
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