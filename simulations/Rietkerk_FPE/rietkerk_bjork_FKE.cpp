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
	if (T(0) < val)
		return 0;
	else
		return 1;
	//Template designed specifically for use in determining the index for the UFDM Scheme.
	//Specifically, the notation follows entirely from determine_neighbours_Sq4().
	//Left (or up) neighbour (index 0) if val > 0, right (or down) neighbour (index 1) if val < 0.
}

// Template function to return smallest positive value in a vector
template <typename T> T min_positive_vec(const std::vector<T>& vec)
{
	T min = std::numeric_limits<T>::max();
	for (const T& val : vec)
	{
		if (val > 0 && val < min)
		{
			min = val;
		}
	}
	return min;
}

//Template function to return smallest positive value in a dynamic array
template <typename T> T min_positive_arr(T* arr, int size)
{
	T min = std::numeric_limits<T>::max();
	for (int i = 0; i < size; i++)
	{
		if (arr[i] > 0 && arr[i] < min)
		{
			min = arr[i];
		}
	}
	return min;
}


//--------------------Defining Primitives------------------------------------//

void increase_stack_limit(long long stack_size) {

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
	int d = a + b + c; cout << d << endl;
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
		}
		else {
			throw std::runtime_error("Key not found: " + key);
		}
	}

	return result;
}



bool comparePairs(const pair<int, int>& a, const pair<int, int>& b)
{
	return (a.first * a.first + a.second * a.second) < (b.first * b.first + b.second * b.second);
}

void set_Prefix(string& user_prefix, double mG /** = 0.0 */, double mP /** = 0.0 */)
{
	//Sets the global prefix for the output files.
	prefix = user_prefix;
	// If mG and mP are integers, convert them to integers, else round them to 2 decimal places.
	mG = (mG == int(mG)) ? int(mG) : round(mG * 100) / 100;
	mP = (mP == int(mP)) ? int(mP) : round(mP * 100) / 100;
	// This is done to ensure that the string representation of mG and mP are consistent across different functions.
	ostringstream mG_str, mP_str;
	// Format the string representation of mG and mP to 2 decimal places at most and remove trailing zeros.
	mG_str << std::defaultfloat << mG; mP_str << std::defaultfloat << mP;
	string SpB_str = "";
	// If non-zero values of mG and mP are provided, include them in the directory structure.
	if (mG > 0.0 && mP > 0.0)
		SpB_str = std::to_string(SpB) + "Sp_mG_" + mG_str.str() + "_mP_" + mP_str.str();
	else if (mG > 0.0 && mP == 0.0)
		SpB_str = std::to_string(SpB) + "Sp_mG_" + mG_str.str();
	else
		SpB_str = std::to_string(SpB) + "Sp";
	frame_folder = "../Data/Rietkerk/Frames/Stochastic/" + SpB_str + "/" + prefix + "_";  //Folder to store frames.
	prelim_folder = "../Data/Rietkerk/Prelims/Stochastic/" + SpB_str + "/" + prefix + "_"; //Folder to store preliminary data.
	stat_prefix = "../Data/Rietkerk/Stochastic/1stOrderCC_Rietkerk_" + prefix + "_P_c_G_"; //Header for frame files.
}

void set_input_Prefix(string& user_inputframepath, string& user_prefix, double user_a_c, double user_dP, double user_Geq /* = -1*/, double user_input_T /* = -1*/)
{
	//Sets the global prefix for the INPUT files (for Burn-In Frames)
	input_frame_subfolder = user_inputframepath; //Folder to store input data.
	input_frame_parenfolder = "../Input/Rietkerk/" + std::to_string(SpB) + "Sp"; //Folder to store input data.

	// Update map input_keys with the new values
	input_keys["outPRE"] = user_prefix;
	input_keys["a_c"] = std::to_string(user_a_c);
	input_keys["dP"] = std::to_string(user_dP);
	if (user_Geq >= 0)
		input_keys["Geq"] = std::to_string(user_Geq);

	if (user_input_T >= 0)
		input_keys["T"] = std::to_string(user_input_T);

}

void set_global_system_params(double dt, double dx)
{
	//Sets the global parameters for the simulation.
	//These are the parameters that are user-defined and are used in the simulation.
	//These are set by the user in the main function

	dt2 = dt / 2.0; dt6 = dt / 6.0;
	dx2 = dx * dx; dx1 = 1.0 / dx; dx1_2 = 1.0 / (dx * dx); double dxb2 = dx / 2.0;
}

void set_global_FKE_params(double(&D)[Sp], double(&v)[SpB], double(&dtAdv)[SpB], double Dt_Analytic, double dx, int L, int CI /** = 4 */)
{
	//Verify Dt_Analytic is equal to smallest non-zero dtAdv, if not, NOTIFY USER and set Dt_Analytic to smallest non-zero dtAdv.
	// Do this using the min_positive function.
	if (Dt_Analytic != min_positive_arr(dtAdv, SpB))
	{
		cout << "WARNING: Dt_Analytic is not equal to the smallest non-zero dtAdv. Setting Dt_Analytic to the smallest non-zero dtAdv.\n";
		Dt_Analytic = min_positive_arr(dtAdv, SpB);
	}

	//Sets the global parameters for the FKE model.
	Dt_ay = Dt_Analytic;


	// Set other global system-wide parameters.
	dt2 = Dt_ay / 2.0; dt6 = Dt_ay / 6.0; L2 = L * L;
	dx2 = dx * dx; dx1 = 1.0 / dx; dx1_2 = 1.0 / (dx * dx); double dxb2 = dx / 2.0;

	// Creating CUDA arrays for the global parameters.
	#if defined(BARRACUDA)
	int2 sigma_vD_ScBounds[Sp_NV];
	vector<int> sig_D_GaussSize(Sp_NV, 0); vector<int> sig_vDXY_GaussSize(Sp_NV, 0);
	// The size of the Gaussian terms for the distance FKE.
	//These are set in a cumulative manner, 
	//s.t sig_D_GaussSize[s] = sig_D_GaussSize[s-1] + 2*sig_D_ScaledBounds[s] + 1
	// and sig_D_GaussSize[0] = 0;
	#endif


	for (int s = 1; s < SpB; s++)
	{
		#if defined(MV_INVARIANCE) && MV_INVARIANCE == 1
		//beta_vel[s] = 3.0 / dtAdv[s]; // Velocity parameter, ~ 3.0/Dt_analytical
		beta_vel[s] = 3.0 / Dt_Analytic; // Velocity parameter, ~ 3.0/Dt_analytical
		#else // Default case with displacement invariance.
		beta_vel[s] = 3.0 / Dt_Analytic; // Velocity parameter, ~ 3.0/Dt_analytical
		#endif
		exp_beta_Dt[s] = exp(-beta_vel[s] * Dt_Analytic); exp_2beta_Dt[s] = exp(-2 * beta_vel[s] * Dt_Analytic);
		// ASSUMING SYMMETRIC VARIANCES for velocity OU process.
		sigma_v[s] = { std::sqrt(std::sqrt(2) * beta_vel[s] * 1 * std::abs(v[s])), std::sqrt(std::sqrt(2) * beta_vel[s] * 1 * std::abs(v[s])) };
		sigma_D[s] = std::sqrt(2 * D[s] * Dt_Analytic);
		double sigma_vx = (sigma_v[s].first / beta_vel[s]) * std::sqrt((Dt_ay + (1 / (2 * beta_vel[s])) * (1 - exp_2beta_Dt[s])
			- (2 / beta_vel[s]) * (1 - exp_beta_Dt[s])));
		double sigma_vy = (sigma_v[s].second / beta_vel[s]) * std::sqrt((Dt_ay + (1 / (2 * beta_vel[s])) * (1 - exp_2beta_Dt[s])
			- (2 / beta_vel[s]) * (1 - exp_beta_Dt[s])));
		// Note: sigma_vx and sigma_vy are the SDs for the Gaussian contribution to the distance FKE by the velocity OU process.
		// ASSUMPTION: The velocity OU process is isotropic, so sigma_vx = sigma_vy.
		sigma_vD[s] = { sigma_vx, sigma_vy };

		sig_D_ScaledBounds[s] = int(CI * sigma_D[s] / dx); // The bounds for the scaled distance OU terms, generally SET to 4*sigma_D/dx.
		// The bounds for the scaled velocity OU terms for the distance FKE, scaled by dx, BY DEFAULT SET to 4*sigma_vD.
		sig_vD_ScaledBounds[s] = { int(CI * sigma_vx / dx), int(CI * sigma_vy / dx) };

		// Finally save prefactors to save computation time.
		mu_vel_prefactor[s] = (1 / dx) * (Dt_ay);// - ((1 - exp_beta_Dt[s]) / beta_vel[s]));

		#if defined(BARRACUDA)
		sigma_vD_ScBounds[(s - 1)] = make_int2(int(CI * sigma_vx / dx), int(CI * sigma_vy / dx));
		if (s >= 2)
		{
			sig_D_GaussSize[(s - 1)] = sig_D_GaussSize[(s - 2)] + sig_D_ScaledBounds[(s - 1)] + 1;
			sig_vDXY_GaussSize[(s - 1)] = sig_vDXY_GaussSize[(s - 2)] + sig_vD_ScaledBounds[(s - 1)].first + 1;
		}
		#endif
	}
	

	size_t size_gauss_stencil_D = 0; size_t size_gauss_stencil_VXY = 0;
	//Stores the size of the Gaussian stencil for the diffusion and advection terms FOR CUDA MEMORY ALLOCATION.
	for (int s = 1; s < SpB; s++)
	{
		//gaussian_stencil_D[s].resize(sig_D_ScaledBounds[s] + 1, 0.0); //Resize the Gaussian stencil for the diffusion term.
		//gaussian_stencil_VXY[s].resize(sig_vD_ScaledBounds[s].first + 1, 0.0); //Resize the Gaussian stencil for the advection term.
		//First calculate the Gaussian stencil for the diffusion term using std::vector<double> precompute_gaussian(double sigma, int Rbound)
		gaussian_stencil_D[s] = precompute_gaussian(sigma_D[s] * dx1, sig_D_ScaledBounds[s]);
		// And for the advection term [RECALL: sig_vD_ScaledBounds[s].first = 2*sig_vD_ScaledBounds[s].second]
		// IMPORTANT: We assume the advection contribution is symmetric in x and y directions
		gaussian_stencil_VXY[s] = precompute_gaussian(sigma_vD[s].first * dx1, sig_vD_ScaledBounds[s].first);
		size_gauss_stencil_D += gaussian_stencil_D[s].size(); size_gauss_stencil_VXY += gaussian_stencil_VXY[s].size();
	}

	#if defined(BARRACUDA)
	// Set the global parameters for the CUDA kernels.
	// Copy sig_D_ScaledBounds, sig_vD_ScaledBounds, mu_vel_prefactor from 1st index onwards to the device.
	cudaError_t err;

	err = copyToDeviceConstantMemory_AdvTerms(sig_D_ScaledBounds.data() + 1, sigma_vD_ScBounds,
		mu_vel_prefactor.data() + 1, sig_D_GaussSize.data(), sig_vDXY_GaussSize.data());

	if (err != cudaSuccess) {
		cerr << "CUDA error in copyToDeviceConstantMemory_AdvTerms: " << cudaGetErrorString(err) << std::endl;
		exit(EXIT_FAILURE);
	}

	// If both size_gauss_stencil_D and size_gauss_stencil_VXY are lesser than
	// MAX_GAUSS_D_STENCIL_SIZE and MAX_GAUSS_VXY_STENCIL_SIZE respectively,
	// copy the Gaussian stencils to the device constant memory.
	if (size_gauss_stencil_D <= MAX_GAUSS_D_STENCIL_SIZE && size_gauss_stencil_VXY <= MAX_GAUSS_VXY_STENCIL_SIZE)
	{
		err = copyToDeviceConstantMemory_GaussianStencils(size_gauss_stencil_D , size_gauss_stencil_VXY);

		if (err != cudaSuccess) {
			cerr << "CUDA error in copyToDeviceConstantMemory_GaussianStencils: " << cudaGetErrorString(err) << std::endl;
			exit(EXIT_FAILURE);
		}
		// Set the flag to indicate that constant memory is being used for Gaussian stencils.
		use_constmem_gaussstencils = true;
	}
	else
	{
		// Set the flag to indicate that constant memory is NOT being used for Gaussian stencils.
		use_constmem_gaussstencils = false;
		stringstream cuda_m0;
		cuda_m0 << "WARNING!!: Gaussian stencil sizes exceed maximum allowed sizes for constant memory in CUDA kernels. \n"
		<< "Size of Gaussian stencil for D: " << size_gauss_stencil_D << " (Max allowed: " << MAX_GAUSS_D_STENCIL_SIZE << ") \n"
		<< "Size of Gaussian stencil for VXY: " << size_gauss_stencil_VXY << " (Max allowed: " << MAX_GAUSS_VXY_STENCIL_SIZE << ") \n"
		<< "Falling back to using global memory for Gaussian stencils in CUDA kernels." << std::endl;
		cout << cuda_m0.str(); cerr << cuda_m0.str();

		//exit(EXIT_FAILURE);
		usleep(static_cast<int>(5000000)); //Sleep for 5 seconds to allow user to read the warning.
	}

#endif


}
void set_global_predator_params(double Km)
{
	/** Sets the global parameters for the predator.
	// Km is the carrying capacity of the predator.
	//These are set by the user in the main function **/
	K_P = Km; // Carrying Capacity of the predator
	K_P1 = 1.0 / Km; // Inverse of the carrying capacity of the predator
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

void init_fullframe(D2Vec_Double& array, int Sp, int size)
{
	//Returns array with all elements initialised to 1
	for (int s = 0; s < Sp; s++)
	{
		for (int j = 0; j < size; j++)
		{
			array[s][j] = 1;
		}
	}
}

void init_constframe(D2Vec_Double& array, int Sp, int size, double constant[])
{
	//Returns array with all elements initialised to a constant value: "const"
	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < size; i++)
		{
			array[s][i] = constant[s];
		}
	}
}


//Initialises a frame (D2Vec_Double vector) with variable percent of the elements randomly set to const c_high, others set to const val c_low.
void init_randconstframe(D2Vec_Double& array, int Sp, int size, double perc, double c_high[], double c_low[])
{

	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);
	//Returns array with all elements initialised to a constant value: "const"

	for (int i = 0; i < size; i++)
	{
		uniform_real_distribution<double> unif;
		unif = uniform_real_distribution<double>(0.0, 1.0);
		double r = unif(rng);
		for (int s = 0; s < Sp; s++)
		{
			if (r <= perc)
				array[s][i] = c_high[s];
			else
				array[s][i] = c_low[s];
		}
	}

}

//Initialises a frame (D2Vec_Double vector) with variable percent of the elements randomly set to const c_high, others set to const val c_low.
// The first Sp elements of c_high and c_low are used if R < R_c, otherwise the last Sp elements are used.
void init_randbistableframe(D2Vec_Double& array, int size, double R, double R_c, double perc, double c_high[], double c_low[])
{
	mutex errorMutex; // Mutex to make error messages thread-safe
	if (Sp != 5)
	{
		std::lock_guard<std::mutex> lock(errorMutex);
		std::cerr << "Error: Number of species must be 3 for this function." << std::endl;
	}
	//c_high and c_low are arrays of size 2*Sp, with the values of the constants for each species.
	// If R < R_c, then the first Sp elements of c_high and c_low are passed to init_randconstframe() to initialise the frame.
	// If R >= R_c, then the last Sp elements of c_high and c_low are passed to init_randconstframe() to initialise the frame.
	if (R < R_c)
		init_randconstframe(array, Sp, size, perc, c_high, c_low);

	else
		init_randconstframe(array, Sp, size, perc, c_high + Sp, c_low + Sp);

}

void init_exprtk_homogenousMFTframe(D2Vec_Double& array, int size, double R, double R_c, double c_spread[])
{
	/** SUMMARY: Same as init_exprtk_randbiMFTframe(), but with the following differences:
	 * dP and perc are not considered and all elements of the frame are set to the same value.
	 * These values are given by cspread[s]*MFT_Vec_CoexExpr[s] for each species s.
	 */

	mutex errorMutex; // Mutex to make error messages thread-safe
	// Create expression objects for each species
	std::vector<exprtk::expression<double>> expressions(2 * Sp);
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
	for (int s = 0; s < 2 * Sp; s++)
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
		for (size_t i = 0; i < size; i++)
		{
			for (size_t s = 0; s < Sp; s++)
			{
				// If the MFT condition before the critical point is 0, set the species to c_spread[s]
				// Otherwise, set the species to the MFT value multiplied by the constant values in c_spread
				if (expressions[s].value() == 0)
					array[s][i] = c_spread[s];
				else
					array[s][i] = c_spread[s] * expressions[s].value();
			}
		}
	}
	else if (R >= R_c)
	{
		for (size_t i = 0; i < size; i++)
		{
			// Set the other species to the constant values
			for (size_t s = 0; s < Sp; s++)
			{
				if (expressions[s + Sp].value() != 0)
					array[s][i] = c_spread[s + Sp] * expressions[s + Sp].value();
				else
					array[s][i] = c_spread[s + Sp];
			}
			// Set the other species to their MFT values multiplied by the constant values in c_spread
		}
	}
}


void init_exprtk_randbiMFTframe(D2Vec_Double& array, int size, double R, double R_c, double dP, double perc, double c_spread[])
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
	std::vector<exprtk::expression<double>> expressions(2 * Sp);
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
	for (int s = 0; s < 2 * Sp; s++)
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
		double highval_veg = std::max(dP, c_spread[0] * (expressions[0].value()));
		/** Use the max function to calculate the maximum of the two values, dP and c_spread[Sp]*expressions[0].value()
		// and assign it to the first species. */
		for (size_t i = 0; i < size; i++)
		{
			double r = unif(rng);
			if (r <= perc)
				array[0][i] = highval_veg; // perc fraction of the vegetation species are set to dP
			else
				array[0][i] = c_spread[0] * expressions[0].value(); // Set the other sites of the first species to the MFT value.

			for (size_t s = 1; s < Sp; s++)
			{
				// If the MFT condition before the critical point is 0, set the species to c_spread[s]
				// Otherwise, set the species to the MFT value multiplied by the constant values in c_spread
				if (expressions[s].value() == 0)
					array[s][i] = c_spread[s];
				else
					array[s][i] = c_spread[s] * expressions[s].value();
			}
		}
	}
	else if (R >= R_c)
	{
		double highval_veg = std::max(dP, c_spread[Sp] * (expressions[Sp].value()));
		/** Use the max function to calculate the maximum of the two values, dP and c_spread[Sp]*expressions[0].value()
		// and assign it to the first species. */
		for (size_t i = 0; i < size; i++)
		{
			double r = unif(rng);
			if (r <= perc)
				array[0][i] = highval_veg;
			else
				array[0][i] = expressions[Sp].value(); // Set the other sites of the first species to the MFT value.

			// Set the other species to the constant values
			for (size_t s = 1; s < Sp; s++)
			{
				if (expressions[s + Sp].value() != 0)
					array[s][i] = c_spread[s + Sp] * expressions[s + Sp].value();
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

void init_exprtk_readCSVcolumns_frame(D2Vec_Double& array, vector<int>& const_index, const string& parendir, const string& filenamePattern,
	const std::vector<std::string>& read_cols, double a, double a_c, int L, int r, double c_spread[])
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
	unif_int = uniform_int_distribution<int>(0, maxRepNo);

	r = r + unif_int(rng); // Randomly select a replicate number
	r = r % (maxRepNo + 1); // Ensure that r is within bounds
	// File name template: "/FRAME_T_{Tval}_a_{a_c}_R_{Rval}.csv"
	string basefilename = filenamePattern + std::to_string(r) + ".csv";

	// Loop through all files in the parent directory matching the pattern.
	for (const auto& entry : fs::directory_iterator(parendir))
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

		// Verify that all columns in read_cols are present in the CSV header
		for (const auto& col : read_cols)
		{
			if (header_map.find(col) == header_map.end())
			{
				stringstream m0;
				m0 << "Column missing in CSV file: " << col << " in file " << filename
					<< "with [a, thr, r]: " << a << " , " << omp_get_thread_num() << " , " << r << std::endl;
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
								array[Sp - 2][row] = std::stod(cell); // Store value in array
							else if (read_cols[i] == "O(x;t)")
								array[Sp - 1][row] = std::stod(cell); // Store value in array
							else
								array[csv_col - 2][row] = std::stod(cell); // Store value in array
						}
						catch (const std::invalid_argument& e)
						{
							std::lock_guard<std::mutex> lock(errorMutex);
							std::cerr << "Error: Unable to convert cell value to double. Column index: " << i
								<< " , for file: " << filename << " on thread #" << omp_get_thread_num()
								<< "\n with the offending line: \n" << line_ss.str() << std::endl;
							//Check if read_cols[i] == "W(x;t)" or "O(x;t)"
							if (read_cols[i] == "W(x;t)")
								array[Sp - 2][row] = 0; // Store value in array
							else if (read_cols[i] == "O(x;t)")
								array[Sp - 1][row] = 0; // Store value in array
							else
								array[csv_col - 2][row] = 0; // Store value in array
						}
						catch (const std::out_of_range& e)
						{
							std::lock_guard<std::mutex> lock(errorMutex);
							std::cerr << "Error: Cell value out of range for double. Column index: " << i
								<< " , for file: " << filename << " on thread #" << omp_get_thread_num()
								<< "\n with the offending line: \n" << line_ss.str() << std::endl;
							if (read_cols[i] == "W(x;t)")
								array[Sp - 2][row] = 0; // Store value in array
							else if (read_cols[i] == "O(x;t)")
								array[Sp - 1][row] = 0; // Store value in array
							else
								array[csv_col - 2][row] = 0; // Store value in array
						}
						catch (const std::exception& e)
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
		if (const_index.size() != 0)
		{

			// Create expression objects for each species
			std::vector<exprtk::expression<double>> expressions(2 * Sp);
			exprtk::parser<double> parser;

			// Create a local symbol table, copying the global symbol table
			exprtk::symbol_table<double> local_symbol_table;// = global_symbol_table;

			uniform_real_distribution<double> unif;
			unif = uniform_real_distribution<double>(0.0, 1.0);
			int thr = omp_get_thread_num();

			local_symbol_table.add_variable("a", a);
			local_symbol_table.add_variable("a_c", a_c);

			// Compile the expressions
			for (int s = 0; s < 2 * Sp; s++)
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
			for (auto& s : const_index)
			{
				if (a < a_c)
				{
					for (size_t i = 0; i < L * L; i++)
					{
						// If the MFT condition before the critical point is 0, set the species to c_spread[s]
						// Otherwise, set the species to the MFT value multiplied by the constant values in c_spread
						if (expressions[s].value() == 0)
							array[s][i] = c_spread[s];
						else
							array[s][i] = c_spread[s] * expressions[s].value();
					}
				}
				else if (a >= a_c)
				{
					// Set the other species to the constant values
					for (size_t i = 0; i < L * L; i++)
					{
						if (expressions[s + Sp].value() != 0)
							array[s][i] = c_spread[s + Sp] * expressions[s + Sp].value();
						else
							array[s][i] = c_spread[s + Sp];
					}
				}
			} // End of for loop for const_index

		} // End of if(const_index.size() != 0)


		inputFile.close();
	} // End of for loop for directory_iterator
}

void init_burnin_wrapper(D2Vec_Double& Rho_dt, double a, double a_c, double dP, double perc, int L, int r, double c_spread[])
{
	stringstream  rain, dPO, Lgrid, t_val; rain << a; Lgrid << L; dPO << dP; t_val << 76000;
	map<string, string> local_input_keys = input_keys; // Copy input_keys to local_input_keys.
	local_input_keys["T"] = t_val.str(); local_input_keys["a"] = rain.str(); local_input_keys["dP"] = dPO.str();
	local_input_keys["L"] = Lgrid.str();
	//Update the values of T, a, dP and L in local_input_keys.
	string local_input_frame_subfolder = input_frame_subfolder; //Copy input_frame_subfolder to local_input_frame_subfolder.
	// Formatted string for input_frame_subfolder.
	string local_input_folder = format_str(local_input_frame_subfolder, local_input_keys); //Format input_folder_pattern with local_input_keys.

	string initcsv_parendir = input_frame_parenfolder + local_input_folder; //Parent directory for csv files.
	vector<int> const_species; //Vector to store species to be initialised with constant values.
	string csv_filename_pattern; //Filename pattern for csv files.
	vector<string> initcsv_columns_temp; //Columns to be read from csv file.
	// Check if local_input_folder contains "HEXBLADE"
	if (local_input_folder.find("HEXBLADE") != std::string::npos || local_input_folder.find("TRUEHEX") != std::string::npos)
	{
		///** INPUT WHEN READING ARTIFICIAL HEX FILES (Filename patterns: FRAME_T_0_a_0_R_{}.csv) 
		csv_filename_pattern = input_prefix + "0_a_0_R_";
		initcsv_columns_temp = { "P(x;t)" };

		for (int s = 1; s < Sp; s++)
			const_species.push_back(s); //Species to be initialised with constant values (ALL EXCEPT VEGETATION).
		//*/
	}
	else
	{
		///**INPUT Filename Patterns are of the form: FRAME_T_{}_a_{}_R_{}.csv (WHEN READING RIETKERK FRAMES)
		csv_filename_pattern = input_prefix + t_val.str() + "_a_" + rain.str() + "_R_";
		initcsv_columns_temp = { "P(x;t)", "W(x;t)", "O(x;t)" }; //Columns to be read from csv file.
		for (int s = 1; s < SpB; s++)
			const_species.push_back(s); //Species to be initialised with constant values.
		//*/
	}

	//Copy the contents of initcsv_columns_temp to initcsv_columns.
	const vector<string> initcsv_columns = initcsv_columns_temp;

	// Remove initcsv_columns_temp from memory.
	vector<string>().swap(initcsv_columns_temp);

	//Check first if initcsv_parendir exists, if not throw an error and exit.
	if (!std::filesystem::exists(initcsv_parendir))
	{
		stringstream m1_1;
		m1_1 << "RUN-TIME ERROR: Parent Directory " << initcsv_parendir <<
			"\nfor Initialisation CSV Files does not exist for Thread Rank:\t " << omp_get_thread_num()
			<< "  with a_value:\t" << a << "\t and Replicate:\t" << r << " EXITING." << endl;
		cout << m1_1.str(); cerr << m1_1.str();
		exit(5);
	}

	init_exprtk_readCSVcolumns_frame(Rho_dt, const_species, initcsv_parendir, csv_filename_pattern, initcsv_columns, a, a_c, L, r, c_spread);

	vector<double> temp_vec = { Rho_dt[0].begin(),Rho_dt[0].end() }; //Rho_dt for species 's'
	double init_veg_num = occupied_sites_of_vector(temp_vec, L * L); //Finds number of occupied at given t.
	vector<double>().swap(temp_vec); //Flush temp_vec out of memory.

	if (init_veg_num == 0.0)
	{
		stringstream m1_1;
		m1_1 << "RUN-TIME WARNING: ZERO Active BURN-IN veg sites for Thread Rank:\t " << omp_get_thread_num()
			<< "  with a_value:\t" << a << "\t and Replicate:\t" << r << " CONSIDER RE-RUNNING OVER A DIFFERENT RANGE.\n";
		cout << m1_1.str();
	}

	stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
	m1_1 << "Initial Conditions: BURN-IN # Active Veg sites:\t" << init_veg_num << " for Thread Rank:\t " << omp_get_thread_num()
		<< "  with a_value:\t" << a << "\t and Replicate:\t" << r << "\n";
	cout << m1_1.str();
}

void init_csvconstframe(D2Vec_Double& array, D2Vec_Double& const_ind_val, const string& parendir, const std::string& filename, const vector<int>& columns, int size)
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
	int columnIndex = columns[columns.size() - 1];
	while (std::getline(inputFile, line))
	{
		std::stringstream ss(line);
		std::string cell;

		//Check if ss starts with "a_c", if so, skip the line.
		if (ss.str().find("a_c") != std::string::npos)
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
					result[i - 2].push_back(value);
				}
				catch (const std::invalid_argument& e)
				{
					std::lock_guard<std::mutex> lock(errorMutex);
					std::cerr << "Error: Unable to convert cell value to double. Column index: " << i
						<< " , for file: " << filename << " on thread #" << omp_get_thread_num()
						<< "\n with the offending line: \n" << ss.str() << std::endl;
					result[i - 2].push_back(0);
				}
				catch (const std::out_of_range& e)
				{
					std::lock_guard<std::mutex> lock(errorMutex);
					std::cerr << "Error: Cell value out of range for double. Column index: " << i
						<< " , for file: " << filename << " on thread #" << omp_get_thread_num()
						<< "\n with the offending line: \n" << ss.str() << std::endl;
					result[i - 2].push_back(0);
				}
			}
		}

	}

	// First store result in array.

	for (int s = 0; s < 1; s++)
	{
		for (int i = 0; i < size; i++)
		{
			array[s][i] = result[s][i + 1];
		}
	}

	for (int i = 0; i < size; i++)
	{
		array[Sp - 2][i] = result[1][i + 1];    // Storing W in the second last row of array
		array[Sp - 1][i] = result[2][i + 1];    // Storing O in the last row of array
	}

	// const_ind_val is a 2D vector of size (Sp - columns.size()) x 2, which stores the constant index 
	//for each column not being extracted (and thus not in result) and the constant value to be substituted in the array at that index value.
	// These values are pre-assigned in the main function.

	for (int i = 0; i < const_ind_val.size(); i++)
	{
		int index = const_ind_val[i][0];
		double value = const_ind_val[i][1];

		for (int j = 0; j < size; j++)
		{
			array[index][j] = value;
		}
	}
	// Close the file
	inputFile.close();

	//Free up memory allocated to result
	vector <vector<double>>().swap(result);
}


void init_randframe(D2Vec_Double& array, int Sp, int size, double mean[], double sd[])
{
	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);
	//Returns array with all elements initialised to random values drawn from a Gaussian distribution with mean mu and SD sigma,
	// with all values 
	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < size; i++)
		{
			normal_distribution<double> norm;
			norm = normal_distribution<double>(mean[s], sd[s]);
			array[s][i] = norm(rng);
			if (array[s][i] < 0)
				array[s][i] = 0.0; //Values less than 0 are unrealistic and not permissible.
		}
	}
}


void init_quarterframe(D2Vec_Double& array, int Sp, int g, double c1[], double c2[], double c3[], double c4[])
{
	//Returns array with elements initialised to a constant values: "const1", "const2", "const3", "const4" in clockwise quadrants
	// starting from the top left.
	int l2 = static_cast<int>(g / 2);
	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < l2; i++)
		{
			int il2 = i + l2;
			for (int j = 0; j < l2; j++)
			{
				array[s][i * g + j] = c1[s]; //Filling up top-left quadrant
				array[s][i * g + l2 + j] = c2[s]; //Filling up top-right quadrant
				array[s][il2 * g + l2 + j] = c3[s]; //Filling up bottom-right quadrant
				array[s][il2 * g + j] = c4[s]; //Filling up bottom-right quadrant
			}
		}
	}
}

void init_gaussframe(D2Vec_Double& array, int size, vector<double>& sd, vector<double>& amp, int s_gauss /* = SpB*/)
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
	D2Vec_Double gauss_frame(s_gauss, vector<double>(g * g, 0.0));
	// Make cx, cy as dynamic arrays of size s_gauss.
	double* cx = new double[s_gauss]; double* cy = new double[s_gauss];
	// Fill cx and cy arrays with Cartesian coordinates of the centre of the Gaussian distribution for each species, starting from the top-left corner.
	for (int i = 0; i < s_gauss; i++)
	{
		// cx and cy reprsent Cartesian coordinates in a counter-clockwise manner starting from the top-left corner of the array.
		// i.e cx =[ 0, 1, 1, 0, 0, 1 ....] and cy = [0, 0, 1, 1, 0, 0 ....]
		cx[i] = (i % 4) % 3 == 0 ? int(g / 8) : int(7 * g / 8);
		cy[i] = (i % 4) < 2 ? int(g / 8) : int(7 * g / 8);
	}
	//If size of sd[] is NOT s_gauss, fill the remaining elements with the last value of sd[], else ignore the extra elements.

	if (sd.size() < s_gauss)
	{
		double last_sd = sd[sd.size() - 1];
		for (int i = sd.size(); i < s_gauss; i++)
			sd.push_back(last_sd);
	}
	else if (sd.size() > s_gauss)
	{
		sd.resize(s_gauss);
	}

	init_gaussiantear(gauss_frame, s_gauss, g, cx, cy, sd, amp);

	//Initialise the constant frame.
	D2Vec_Double const_frame(Sp - s_gauss, vector<double>(g * g, 0.0));
	for (int i = 0; i < Sp - s_gauss; i++)
	{
		for (int j = 0; j < g * g; j++)
		{
			const_frame[i][j] = amp[s_gauss + i];
		}
	}

	//Combine the two frames to get the final frame.
	for (int s = 0; s < s_gauss; s++) {
		for (int i = 0; i < g * g; i++)
			array[s][i] = gauss_frame[s][i];

	}
	for (int s = 0; s < Sp - s_gauss; s++) {
		for (int i = 0; i < g * g; i++)
			array[s + s_gauss][i] = const_frame[s][i];
	}

	std::lock_guard<std::mutex> lock(errorMutex);
	ostringstream oss;
	oss << "Gaussian frame initialised with " << s_gauss << " species.\n";
	oss << "Amp values:\n";
	for (int i = 0; i < Sp; i++)
		oss << amp[i] << " ";
	oss << "\n";
	oss << "SD values:\n";
	for (int i = 0; i < s_gauss; i++)
		oss << " ( " << i << " , " << sd[i] << " ) \t";
	oss << "\n";
	oss << "cx values & cy values:\n";
	for (int i = 0; i < s_gauss; i++)
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

void init_gaussiantear(D2Vec_Double& array, int Sp_lim, int g, double cx[], double cy[], const vector<double>& sd, const vector<double>& amp)
{
	/**
	 * Returns array with elements initialised to a 2D Gaussian tear drop with mean (cx, cy) and standard deviation sd.
	 * array is a 2D vector of size Sp x g^2. cx, cy, amp & sd are arrays of size Sp, with the values of the constants for each species.
	 */
	for (int s = 0; s < Sp_lim; s++)
	{
		double A = amp[s] / (2 * PI * sd[s] * sd[s]);
		for (int i = 0; i < g; i++)
		{
			for (int j = 0; j < g; j++)
				array[s][i * g + j] = amp[s] * exp(-((i - cx[s]) * (i - cx[s]) + (j - cy[s]) * (j - cy[s])) / (2 * sd[s] * sd[s]));
		}
	}
}

void init_solitarytear(D2Vec_Double& array, int Sp, int length)
{
	//Returns array with only central element initalised to 1, all others 0.
	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < length * length; i++)
			array[s][i] = 0;
		int lhalf = static_cast<int>(length / 2.0);
		array[s][length * lhalf + lhalf] = 1;
	}
}

double mean_of_array(double array[], int size) {

	double sum = 0.0;

	for (int i = 0; i < size; ++i) {
		sum += array[i];
	}
	return sum / (double)size;
}

double standard_deviation_of_array(double array[], int size) {

	double mean = mean_of_array(array, size);
	double sum = 0.0;

	for (int i = 0; i < size; ++i) {
		sum += (array[i] - mean) * (array[i] - mean);
	}

	double variance = sum / (double)size;

	return sqrt(variance);
}

double mean_of_vector(vector<double> array, int size) {

	double sum = 0.0;

	for (int i = 0; i < size; ++i) {
		sum += array[i];
	}
	return sum / (double)size;
}

double variance_of_vector(vector<double> array, int size) {
	double mean = mean_of_vector(array, size);
	float sum = 0.0;

	for (int i = 0; i < size; ++i) {
		sum += (array[i] - mean) * (array[i] - mean);
	}
	double variance = sum / (double)size;
	return variance;
}

double standard_deviation_of_vector(vector<double> array, int size) {

	double mean = mean_of_vector(array, size);
	float sum = 0.0;

	for (int i = 0; i < size; ++i) {
		sum += (array[i] - mean) * (array[i] - mean);
	}
	double variance = sum / (double)size;
	return sqrt(variance);
}

double occupied_sites_of_vector(vector<double> array, int size) {
	//Calculates number of occupied sites in a vector.
	double sum = 0.0;

	for (int i = 0; i < size; ++i) {
		if (array[i] > 0)
			sum += 1;
	}
	return sum;
}

auto meansq_spread_of_vector(vector<double> array, int g, int c_x, int c_y) {
	//Calculates number of occupied sites in a vector. c_x and c_y represent the Cartesian coordinates of the central site.
	struct coord {        // Declare a local structure 
		double i1; int i2;
	};
	double sq_dist = 0.0; int sum = 0;
	for (int i = 0; i < g * g; ++i) {
		if (array[i] > 0)
		{
			int p = int(i / g); int q = i % g; //Indices of occupied site.
			sum += 1;
			sq_dist += (p - c_x) * (p - c_x) + (q - c_y) * (q - c_y); //Calcuates norm-2 distance.
		}
	}
	return coord{ sq_dist, sum };
}

void generic_SPBmean_surv_runs(D2Vec_Double& t_avg_var_rep_N, const D2Vec_Double& X_curr, int size, int spcolnum, int j)
{
	/** X_curr contains time-series data (indexed as (N(t), <rho(t)>x)) from the current surviving replicate (given by r)
	 * with a length (number of rows) "size", and spcolnum is the columns per species in X_curr (so X_curr will have spcolnum*SpB columns).
	 *  rep_avg_var stores the time(t), running density average and variance of X_curr, as well as N(t) and number of survivng replicates
	 * for each time-step from previous surviving replicates ( <= r-1) respectively.
	 * In this function, X_curr is used to update the running avg and variance of rep_avg_var to account for current surviving replicate r(t).
	 * Note rep_avg_Var[t][3] stores number of surviving runs at t.
	*/
	for (int s = 0; s < SpB; s++)
	{
		for (int i = 0; i < size; i++)
		{
			// For each species s, save the values of X_curr in the corresponding columns of t_avg_var_rep_N, iterating using spcolnum.
			for (int k = 0; k < spcolnum; k++)
			{
				t_avg_var_rep_N[i][spcolnum * s + k + 1] = X_curr[i][spcolnum * s + k];
			}
		}
	}
}

void var_mean_incremental_surv_runs(D2Vec_Double& t_avg_var_rep_N, const D2Vec_Double& X_curr, int size, int j)
{
	/** X_curr contains time-series data (indexed as (N(t), <rho(t)>x)) from the current surviving replicate (given by r) with a length "size".
	 *  rep_avg_var stores the time(t), running density average and variance of X_curr, as well as N(t) and number of survivng replicates
	 * for each time-step from previous surviving replicates ( <= r-1) respectively.
	 * In this function, X_curr is used to update the running avg and variance of rep_avg_var to account for current surviving replicate r(t).
	 * Note rep_avg_Var[t][3] stores number of surviving runs at t.
	*/
	vector<double> mean_prev(size);
	for (int s = 0; s < SpB; s++)
	{
		for (int i = 0; i < size; i++)
		{
			int r = t_avg_var_rep_N[i][4 * s + 5] + 1; //Number of surviving replicates
			if (X_curr[i][2 * s + 0] == 0.0)
			{	//Non-surviving run. No-updates required. Also <rho(t)>x == 0, it remains 0 for t' > t.
				break;
			}
			else
			{	//Number of surviving run at time t.
				mean_prev[i] = t_avg_var_rep_N[i][4 * s + 3]; //Stores X_n-1 mean data as temp variable.
				if (t_avg_var_rep_N[i][4 * s + 5] == 0.0)
				{	//First surviving run encountered at time t, encountered here.
					t_avg_var_rep_N[i][4 * s + 3] = X_curr[i][2 * s + 1]; t_avg_var_rep_N[i][4 * s + 4] = 0.0;
					t_avg_var_rep_N[i][4 * s + 6] = X_curr[i][2 * s + 0];
				}
				else
				{	//Not the first surviving run, which means r = rep_avg_var[i][3]+1 >=2
					t_avg_var_rep_N[i][4 * s + 3] = (X_curr[i][2 * s + 1] + (r - 1) * t_avg_var_rep_N[i][4 * s + 3]) / r; //Note r>= 1.
					t_avg_var_rep_N[i][4 * s + 6] = (X_curr[i][2 * s + 0] + (r - 1) * t_avg_var_rep_N[i][4 * s + 6]) / r; //Note r>= 1.
					//Calculating and storing incremental variance (V^2_n) [For formula derivation, refer to notes.]
					t_avg_var_rep_N[i][4 * s + 4] = (t_avg_var_rep_N[i][4 * s + 3] - mean_prev[i]) * (t_avg_var_rep_N[i][4 * s + 3] - mean_prev[i]) +
						(1 / (r - 1)) * ((r - 2) * t_avg_var_rep_N[i][4 * s + 4] +
							(X_curr[i][2 * s + 1] - t_avg_var_rep_N[i][4 * s + 3]) * (X_curr[i][2 * s + 1] - t_avg_var_rep_N[i][4 * s + 3]));
				}
				t_avg_var_rep_N[i][4 * s + 5] = r; //At the end, update number of succesful runs.
			}
		}
	}
	for (int i = 0; i < size; i++)
	{	//For soil water (W) and surface water (0) respectively
		t_avg_var_rep_N[i][1] = (X_curr[i][2 * (Sp - 2) + 1] + (j)*t_avg_var_rep_N[i][1]) / (j + 1);
		t_avg_var_rep_N[i][2] = (X_curr[i][2 * (Sp - 1) + 1] + (j)*t_avg_var_rep_N[i][2]) / (j + 1);
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

	for (int i = 0; i < num - 1; ++i)
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
	double delta = (logMax - logMin) / (log_points - 1);
	if (log_points == 0) { return logspaced; }
	if (log_points == 1)
	{
		logspaced.push_back(static_cast<double>(start_in));
		return logspaced;
	}
	for (int i = 0; i < log_points - 1; ++i)
	{
		double v = pow(logarithmicBase, logMin + delta * i);
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

template<typename T>std::vector<T> switchsort_and_bait(std::vector<T> vals, T low /**= 0 */, T high /**= 0 */,
	int num_points /**= 50 */, string insert_type /** = "None" */, bool erase /**= false*/)
{
	// Function to switch and bait elements in a vector.
	if (erase)
		vals.erase(std::remove_if(vals.begin(), vals.end(), [low, high](T x) {return x >= low && x <= high;}), vals.end());

	// Next check if low < high, if not, swap the values.
	if (low > high) {
		T temp = low; low = high; high = temp;
	}
	if (insert_type == "linspace")
	{
		std::vector<double> lin_spaced = linspace(low, high, num_points);
		vals.insert(vals.end(), lin_spaced.begin(), lin_spaced.end());
	}
	else if (insert_type == "lnspace")
	{
		std::vector<double> log_spaced = lnspace(low, high, num_points);
		vals.insert(vals.end(), log_spaced.begin(), log_spaced.end());
	}
	else
		std::cout << "Invalid insert_type. Sorting and returning original vector. \n";

	std::sort(vals.begin(), vals.end());
	vals.erase(unique(vals.begin(), vals.end()), vals.end()); //Remove duplicates if any.
	return vals;
}

std::vector<double> logarithmic_time_bins(double t_max, double dt)
{	// Computes and returns ln-distributed points from t= e^2 to ln(t_max) (the latter rounded down to 1 decimal place)
	std::vector<double> logspaced;
	if (t_max <= exp(2.0))
	{
		cout << "GOUT!!" << endl;
		return logspaced; // No points recorded if t_max < e^{2}
	}

	double tend_one = static_cast<double>(floor(log(t_max) * 10.)) / 10.0;
	// Finds the closest power of e (rounded down to one decimal place) from t_max.
	double pow = 2.0;
	while (pow <= tend_one)
	{
		double t_raw = exp(pow); // Getting raw time-stamp.
		//Rounding down raw time-stamp to dt.
		double t_fin = static_cast<int>(t_raw / dt) * dt;

		if (t_fin < exp(2))
			t_fin = t_fin + dt; //Ensure t_fin is strictly above e^{2}

		if (logspaced.size() > 0 && t_fin != logspaced[logspaced.size() - 1])
		{
			logspaced.push_back(t_fin);
		} // This is to prevent duplicate values from arising in logspaced.
		else if (logspaced.size() == 0)
			logspaced.push_back(t_fin);
		pow += 0.04; //25 time-stamps per decade.
	}
	cout << "TROUT! " << logspaced.size() << endl;
	return logspaced;


}

std::vector<double> logarithm10_time_bins(double t_max, double dt)
{	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place)
	std::vector<double> logspaced;
	if (t_max <= 1)
	{
		cout << "GOUT!!" << endl;
		return logspaced; // No points recorded if t_max < e^{2}
	}

	double tend_one = static_cast<double>(floor(log10(t_max) * 10.)) / 10.0;
	// Finds the closest power of e (rounded down to one decimal place) from t_max.
	double power = 0.0;
	while (power <= tend_one)
	{
		double t_raw = pow(10, power); // Getting raw time-stamp.
		//Rounding down raw time-stamp to dt.
		double t_fin = static_cast<int>(t_raw / dt) * dt;

		if (logspaced.size() == 1)
			cout << "Feq " << logspaced[0] << " HODOR " << t_fin << " modor " << t_fin - logspaced[0] << endl;

		if (t_fin <= 1.0)
			t_fin = t_fin + dt; //Ensure t_fin is strictly above 1.0

		if (logspaced.size() > 0 && t_fin != logspaced[logspaced.size() - 1])
		{
			logspaced.push_back(t_fin);
		} // This is to prevent duplicate values from arising in logspaced.
		else if (logspaced.size() == 0)
			logspaced.push_back(t_fin);
		power += 0.04; //25 time-stamps per decade.
	}
	int n = logspaced.size();
	for (int i = n - 1; i > 0; i--)
	{	//Rooting out duplicates.
		if (logspaced[i] < logspaced[i - 1] + dt / 2.0)
			logspaced.erase(logspaced.begin() + i);
	}
	cout << "TROUT! " << logspaced.size() << endl;
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
		if (regex_match(matchStr, match, pattern))
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

int theborderwall(D2Vec_Double& Rho_t, int g)
{
	//Checks if the border sites are occupied.
	for (int s = 0; s < Sp; s++)
	{
		int g2 = g * (g - 1);
		for (int i = 0; i < g; i++)
		{
			if (Rho_t[s][i] > 0 || Rho_t[s][g * i] > 0 || Rho_t[s][g2 + i] > 0 || Rho_t[s][g * i + g - 1] > 0)
				return 1; //Top, left, bottom and right borders are occupied by active sites.
		}

	}
	return 0; //No active sites detected at the edges.
}

void determine_neighbours_R2(int g, D3Vec_Int& neighbours_R2)
{
	//Stores neighbouring sites (Van-Neumann Radius of 2) of each cell site in a 
	// S*(g^2)*2*5 matrix.
	int a = 0; int b = 0; //Useful for accounting purposes
	for (int i = 0; i < g * g; i++)
	{
		int p = int(i / g); int q = i % g;
		//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
		for (int j = 0; j < 5; j++)
		{
			int k = j - 2;
			a = ((p + k) % g + g) % g; // Eq to (p+k)mod g. -2 <= k <= 2
			b = ((q + k) % g + g) % g; // Eq to (q+k)mod g. -2 <= k <= 2

			//Recall, first row stores vertical neighbours, second row stores horizontal neighbours.
			neighbours_R2[i][0][j] = a * g + q; neighbours_R2[i][1][j] = p * g + b;
		}
	}

}

void determine_neighbours_Sq4(int g, D3Vec_Int& neighbours_Sq4)
{

	//Stores neighbouring sites of species s and at site i in a square (Norm 1) radius of 1

	// Neighbours_Sq8 has dimensions S*(L*L)*2*2 where S is the number of species (=2) in this case
	// First row stores vertical neighbours, second row stores horizontal neighbours.
	int a = 0; int b = 0; //Useful for accounting purposes

	for (int i = 0; i < g * g; i++)
	{
		int p = int(i / g); int q = i % g;
		//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.

		int y_up = ((p - 1) % g + g) % g; int y_down = ((p + 1) % g + g) % g; // Eq to (p+ky)mod g. -1 <= ky <= 1
		int x_left = ((q - 1) % g + g) % g; int x_right = ((q + 1) % g + g) % g; // Eq to (q+kx)mod g. -1 <= kx <= 1
		neighbours_Sq4[i][0][0] = y_up * g + q; neighbours_Sq4[i][0][1] = y_down * g + q;
		neighbours_Sq4[i][1][0] = p * g + x_left; neighbours_Sq4[i][1][1] = p * g + x_right;
	}
}


void determine_neighbours_Sq8(int g, D3Vec_Int& neighbours_Sq8)
{

	//Stores neighbouring sites of species s and at site i in a square (Norm 1) radius of 1

	// Neighbours_Sq8 has dimensions S*(L*L)*3*3 where S is the number of species (=2) in this case
	int a = 0; int b = 0; //Useful for accounting purposes

	for (int i = 0; i < g * g; i++)
	{
		int p = int(i / g); int q = i % g;
		//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
		for (int y = 0; y < 3; y++)
		{
			int ky = y - 1; a = ((p + ky) % g + g) % g; // Eq to (p+ky)mod g. -1 <= ky <= 1
			for (int x = 0; x < 3; x++)
			{
				int kx = x - 1;
				b = ((q + kx) % g + g) % g; // Eq to (q+kx)mod g. -1 <= kx <= 1
				//Recall, for given s and i first (third) row stores vertical neighbours, second (fourth) row stores horizontal neighbours.
				neighbours_Sq8[i][y][x] = a * g + b;
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
	int maxDistanceIndex = int(f * f * numNeighbors);
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

// Precompute Gaussian stencil, of size Rbound+1, for given sigma. 
std::vector<double> precompute_gaussian(double sigma, int Rbound)
{
	std::vector<double> stencil(Rbound + 1);
	double sum = 0.0;
	double factor = 1.0 / (std::sqrt(2 * PI) * sigma);
	for (int x = 0; x <= Rbound; ++x) {
		stencil[x] = factor * std::exp(-0.5 * (x * x) / (sigma * sigma));
		sum += (x > 0) ? 2 * stencil[x] : stencil[0];
	}
	// Renormalise stencil to ensure total probability = 1.
	for (int x = 0; x <= Rbound; ++x)
		stencil[x] /= sum;

	return stencil;
}



//----------------------------- Regular Rietkerk Integration Machinery --------------------------------------------//


// ALT INTEGRATION MACHINERY WHERE MIXING OF TIMESTEPS DUE TO LAPLACIAN MAY BE AN ISSUE.



void f_2D(D2Vec_Double& f, D2Vec_Double& Rho_M, D3Vec_Int& nR2, double a, double c, double gmax, double alpha, double d, double rW, double W0,
	double D[], double K[], double t, double dt, double dx2, int g)
{
	//Vector function that updates an array containing (dP/dt, dW/dt, DO/dt) for each site in the lattice.
	//Based on the Rietkerk model for plant vegetation dynamics.

	for (int i = 0; i < g * g; i++)
	{
		//Equations for the density of plants, soil water and surface water at each site.
		//Note that the Laplacian is calculated using reflective boundary conditions.
		//The Laplacian is calculated using the 5-point stencil method.
		//Equivalent to dP/dt = c*g_max*P*W/(W+K1) - d*P + D*(Laplacian of P)
		f[0][i] = c * gmax * Rho_M[1][i] * Rho_M[0][i] / (Rho_M[1][i] + K[1]) - d * Rho_M[0][i] + (D[0] / (dx2)) * (Rho_M[0][nR2[i][0][1]] + Rho_M[0][nR2[i][2][1]]
			+ Rho_M[0][nR2[i][1][0]] + Rho_M[0][nR2[i][1][2]] - 4 * Rho_M[0][i]);
		//Equivalent to dW/dt = alpha*(P+K2*W)/(P+K2) - rW*W + D*(Laplacian of W)
		f[1][i] = alpha * (Rho_M[0][i] + K[2] * W0) / (Rho_M[0][i] + K[2]) * Rho_M[2][i] - gmax * Rho_M[1][i] * Rho_M[0][i] / (Rho_M[1][i] + K[1]) - rW * Rho_M[1][i]
			+ (D[1] / (dx2)) * (Rho_M[1][nR2[i][0][1]] + Rho_M[1][nR2[i][2][1]] + Rho_M[1][nR2[i][1][0]] + Rho_M[1][nR2[i][1][2]] - 4 * Rho_M[1][i]);
		//Equivalent to dO/dt = a - alpha*(P+K*W)/(P+K)*O + D*(Laplacian of O)
		f[2][i] = a - alpha * (Rho_M[0][i] + K[2] * W0) / (Rho_M[0][i] + K[2]) * Rho_M[2][i] + (D[2] / (dx2)) * (Rho_M[2][nR2[i][0][1]] + Rho_M[2][nR2[i][2][1]]
			+ Rho_M[2][nR2[i][1][0]] + Rho_M[2][nR2[i][1][2]] - 4 * Rho_M[2][i]);

	}



}


//Modify the function below to perform RK4 integration of the 2D Vector Rho_t, with a laplacian and time step of dt and spatial step dx.
//The function should also take in the parameters of the model, and the time t at which the integration is being performed.
//The function should also take in the number of grid points g, and the number of species Sp.
//The function should also take in the number of threads being used, and the thread number.
//The function should also take in the number of iterations performed so far, and the number of iterations to be performed.


void RK4_Integrate(D2Vec_Double& Rho_t, D2Vec_Double& Rho_tsar, D3Vec_Int& nR2, double a, double c, double gmax, double alpha, double d, double rW,
	double W0, double D[], double K[], double t, double dt, double dx, int g)
{


	// Rho_t stores value of rho at time t.

	// Rho_tstar stores value of rho*, provided by sampling Gaussian.
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	double dt6 = dt / 6.0; double dt2 = dt / 2.0; double dx2 = dx * dx;

	D2Vec_Double K1(Sp, vector<double>(g * g, 0.0)); D2Vec_Double K2(Sp, vector<double>(g * g, 0.0));
	D2Vec_Double K3(Sp, vector<double>(g * g, 0.0)); D2Vec_Double K4(Sp, vector<double>(g * g, 0.0));
	D2Vec_Double Rho_M(Sp, vector<double>(g * g, 0.0));

	//double K1[Sp][g*g] ={0.0}; double K2[Sp][g*g] ={0.0}; double K3[Sp][g*g] ={0.0}; double K4[Sp][g*g] ={0.0};
	//double Rho_M[Sp] ={0.0, 0.0};


	f_2D(K1, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t, dt, dx2, g); //K1 updated.

	for (int i = 0; i < g * g; i++)
	{
		for (int s = 0; s < Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K1[s][i];
	}

	f_2D(K2, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt2, dt, dx2, g); //K2 updated.

	for (int i = 0; i < g * g; i++)
	{
		for (int s = 0; s < Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K2[s][i];
	}
	f_2D(K3, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt2, dt, dx2, g); //K3 updated.
	for (int i = 0; i < g * g; i++)
	{
		for (int s = 0; s < Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt)*K3[s][i];
	}
	f_2D(K4, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt, dt, dx2, g); //K4 updated.

	for (int i = 0; i < g * g; i++)
	{
		for (int s = 0; s < Sp; s++)
		{
			Rho_t[s][i] += (dt6) * (K1[s][i] + 2.0 * K2[s][i] + 2.0 * K3[s][i] + K4[s][i]);


			if (Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || isnan(Rho_t[s][i]) == true)
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

void RK4_Wrapper_2D(D2Vec_Double& Rho, vector <double>& t_meas, double t_max, double a, double c, double gmax,
	double alpha, double d, double rW, double W0, double D[], double K[], double a_st, double a_end, double dt, double dx, double dP, int r, int g)
{
	// Rho_0 stores value of rho at time t=0.
	// Rho_t stores value of rho at time t.
	// Rho_tstar stores value of rho*, provided by sampling Gaussian.
	// Rho_M stores value of rho at intermediate time steps.

	//D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0));
	//double Rho_tsar[Sp][g*g] ={0.0};
	//double Rho_M[Sp][g*g] ={0.0, 0.0};

	D3Vec_Int nR2(g * g, D2Vec_Int(3, vector<int>(3, 0)));
	//nR2 is 3D [(L*L)x3x3] vector initialised to 0, which stores neighbours of all sites i (second dim) in a ball of NORM 1 (SQUARE) radius 1 (2ND & 3RD dim).
	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_Sq8(g, nR2); //Assigning neigbours.
	int tot_iter = static_cast<int>(1.0 / dt) + 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	t_max = t_meas[t_meas.size() - 1]; //Only compute until end of sampling events

	//t_meas[0] = 100; t_meas[1] = 200; t_meas[2] = 300; 

	//double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0};
	D2Vec_Double rho_rep_avg_var(tot_iter, vector<double>(Sp4_1, 0.0));

	int ind = 0;
	for (double t = 0; t <= 1; t += dt)
	{
		rho_rep_avg_var[ind][0] = t; ind += 1;
	}
	for (int i = 0; i < t_meas.size(); i++)
	{
		rho_rep_avg_var[ind + i][0] = t_meas[i];
	}  //Setting up time index of rho_t_N_R2_rep

	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
	m0 << "Parameters:\t a:\t" << a << "\t c:\t " << c << "\t dt,dx:\t " << dt << "," << dx << "\t gmax:\t " << gmax
		<< "\t and  Alpha:\t " << alpha << endl; //cout << m1.str();
	cout << m0.str();
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	double beta = c * gmax - d;
	double p0istar, p0jstar, p0mstar; // Analytic steady state values.

	if (a < 1)
	{
		p0istar = 0.0; p0jstar = a / rW; p0mstar = a / (alpha * W0);
	}
	else
	{
		p0jstar = d * K[1] / beta;
		p0istar = (c / d) * (a - rW * p0jstar);
		p0mstar = (a / alpha) * (p0istar + K[2]) / (p0istar + K[2] * W0);
	}

	for (int j = 0; j < r; j++)
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


		D2Vec_Double Rho_dt(Sp, vector<double>(g * g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double Rho_tsar(Sp, vector<double>(g * g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		//double Rho_M[tot_iter][Sp2] ={0.0}; //Stores <rho>x values per time step for a single replicate.
		D2Vec_Double Rho_M(tot_iter, vector<double>(Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		double perc = 0.015; double c_high[Sp] = { p0istar + dP, p0jstar, p0mstar }; double c_low[Sp] = { p0istar, p0jstar, p0mstar };

		init_randconstframe(Rho_dt, Sp, g * g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.


		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
		m1 << "Initial Conditions:\t Per:\t" << perc << "\t C_High[0], C_High[1], C_High[2]:\t " << c_high[0] << "," << c_high[1] << "," << c_high[2]
			<< "\t C_Lo[0], C_Lo[1], C_Lo[2]:\t " << c_low[0] << "," << c_low[1] << "," << c_low[2] << ",\t R: " << a << endl; //cout << m1.str();
		cout << m1.str();
		//std::ofstream errout; //Save Error Logs
		//std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
		errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();

		for (int i = 0; i < g * g; i++)
		{
			for (int s = 0; s < Sp; s++)
			{
				Rho_tsar[s][i] = Rho_dt[s][i];
			}
		}

		double t = 0; int index = 0;  //Initialise t
		int iter_index = 0; int po = 1; int so = 1; int lo = 1;

		while (t < t_max + dt)
		{
			if (t <= 1.0 || t >= t_meas[index] - dt / 2.0 && t < t_meas[index] + dt / 2.0)
			{

				double rhox_avg, rhox_num; //Stores spatial average and number of occupied sites at given t.
				for (int s = 0; s < Sp; s++)
				{
					vector <double> temp = { Rho_dt[s].begin(),Rho_dt[s].end() }; //Rho_dt for species 's'
					rhox_avg = mean_of_vector(temp, g * g); //Finds spatial average of densities at given t.
					rhox_num = occupied_sites_of_vector(temp, g * g); //Finds number of occupied at given t.
					Rho_M[iter_index][2 * s] = rhox_num; Rho_M[iter_index][2 * s + 1] = rhox_avg;

					vector<double>().swap(temp); //Flush temp out of memory.0

					if (rhox_num == 0.0 and s == 0)
						break; //Finish negotiations as all basal species are dead.
				}
				iter_index += 1;//Rho_M.push_back({t, rhox_avg});

				if (t >= 100 && t <= 2500 || t == 0 || iter_index >= tot_iter - 3 && iter_index <= tot_iter)
				{
					//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.

					stringstream L, tm, d3, p1, rini, gm, a1, a2, Dm0, Dm1, alph, w0t, dix, dimitri;

					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
					rini << j; Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; gm << setprecision(3) << gmax;
					w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; dimitri << setprecision(3) << dP;
					a1 << a_st; a2 << a_end;
					// Three replicates are over.

					ofstream frame_dp;
					// Creating a file instance called output to store output data as CSV
					frame_dp.open("../Data/Rietkerk/Frames/Non-Stochastic/" + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "/FRAME_RAND_Basic_P_c_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str()
						+ "_D0_" + Dm0.str() + "_D1_" + Dm1.str() + "_dx_" + dix.str() + "_alph_" + alph.str() + "_gmax_" + gm.str() + "_R_" + rini.str() + ".csv");

					// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|
					frame_dp << "a_c,  x,  P(x; tmax), W(x; tmax), O(x; tmax) \n";
					for (int i = 0; i < g * g; i++)
					{
						frame_dp << a << "," << i << "," << Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] << endl;

					}

					frame_dp.close();


				}

				if (t > 1.0)
					index += 1;

			}

			for (int s = 0; s < Sp; s++)
				Rho_tsar[s] = Rho_dt[s]; //Assigning initial conditions

			RK4_Integrate(Rho_dt, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t, dt, dx, g);

			for (int s = 0; s < Sp; s++)
			{
				for (int i = 0; i < g * g; i++)
				{
					if (Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true and so == 1)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num()
							<< " at time:\t:" << t << " with Rho(t-dt,i):   " << Rho_tsar[s][i] << "\t and Analytic Integration Term:\t"
							<< Rho_tsar[s][i] << endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}

				}

			}


			t += dt;



		} //End of Time Loop


		if (int(omp_get_thread_num()) == 5)
		{
			size_t currentSize = getCurrentRSS(); //Check Memory usage.
			size_t peakSize = getPeakRSS();
			stringstream m7;
			m7 << "In replicate " << j << ", current size in MB: " << currentSize / (1024.0 * 1024.0)
				<< " and Peak Size (in MB): " << peakSize / (1024.0 * 1024.0) << endl; cout << m7.str();
			cerr << m7.str();
		}

		if (j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
			m7 << "INITIAL UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num()
				<< " THE SIZE OF RHO_M: " << tot_iter << endl; cout << m7.str();
			for (int i = 0; i < tot_iter; i++)
			{
				for (int s = 0; s < 1; s++)
				{
					rho_rep_avg_var[i][4 * s + 3] = Rho_M[i][2 * s + 1]; // Avg density of frame at given time point i.
					rho_rep_avg_var[i][4 * s + 4] = 0.0; // Variance in avg frame density.
					rho_rep_avg_var[i][4 * s + 6] = Rho_M[i][2 * s]; // Number of active sites in frame at given time point i.

					if (Rho_M[i][2 * s] > 0) //Ensuring only surviving runs are considered
					{
						rho_rep_avg_var[i][4 * s + 5] = 1;
					} //One surviving run as of now
					else if (Rho_M[i][2 * s] == 0) //Ensuring only surviving runs are considered
					{
						rho_rep_avg_var[i][4 * s + 5] = 0;
					} //No surviving run as of now
				}
				rho_rep_avg_var[i][1] = Rho_M[i][2 * 1 + 1]; // Avg density of frame at given time point i.
				rho_rep_avg_var[i][2] = Rho_M[i][2 * 2 + 1]; // Avg density of frame at given time point i.



			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			// Second or higher replicate, use incremental advances.
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter, j);
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		if ((j + 1) % 3 == 1)
		{	//Start logging at every multiple of 3
			stringstream L, tm, d3, p1, a1, a2, dimitri, rini, Dm, gm;

			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dimitri << dP;
			rini << j; Dm << setprecision(3) << D[0]; gm << gmax;
			// Three replicates are over.

			ofstream output_dp;
			// Creating a file instance called output to store output data as CSV.
			output_dp.open("../Data/Rietkerk/Prelims/Non-Stochastic/" + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "/PRELIM_AGGRAND_Basic_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
				"_D0_" + Dm.str() + "_gmax_" + gm.str() + "_R_" + rini.str() + ".csv");

			// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
			output_dp << " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs, # Active Sites, \n";
			for (int i = 0; i < tot_iter; i++)
			{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |

				output_dp << a << "," << j << "," << g << "," << rho_rep_avg_var[i][0] << "," << rho_rep_avg_var[i][1] << "," << rho_rep_avg_var[i][2] << ",";
				for (int s = 0; s < 1; s++)
				{
					output_dp << rho_rep_avg_var[i][4 * s + 3] << "," <<
						rho_rep_avg_var[i][4 * s + 4] << "," << rho_rep_avg_var[i][4 * s + 5] << "," << rho_rep_avg_var[i][4 * s + 6] << ",";
				}
				output_dp << endl;

			}

			output_dp.close();
		}




	} //End of Replicate Loop


	size_t currentSize = getCurrentRSS(); //Check Memory usage.
	size_t peakSize = getPeakRSS();
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize / (1024.0 * 1024.0)
		<< " and Peak Size (in MB): " << peakSize / (1024.0 * 1024.0) << endl; cout << m9.str(); cerr << m9.str();


	for (int i = 0; i < tot_iter; i++)
	{	// Recall Rho is: | 	a		|    L 		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |

		Rho.push_back({ a, static_cast<double>(r), static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1],
					rho_rep_avg_var[i][2] });
		for (int s = 0;s < Sp - 2; s++)
		{


			Rho[i].insert(Rho[i].end(), { rho_rep_avg_var[i][4 * s + 3], rho_rep_avg_var[i][4 * s + 4],
				rho_rep_avg_var[i][4 * s + 5], rho_rep_avg_var[i][4 * s + 6] });
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
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize / (1024.0 * 1024.0)
		<< " and Peak Size (in MB): " << peakSize / (1024.0 * 1024.0) << endl; cout << m10.str(); cerr << m10.str();
}



void first_order_critical_exp_delta(int div, double t_max, double a_start, double a_end, double c, double gmax, double alpha,
	double d, double rW, double W0, double D[], double K[], double dt, double dx, double dP, int r, int g)
{


	//double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2
	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//double mean[Sp] = {1.0, 0.05}; double sd[Sp] = {0.25, 0.0125};
	//init_randframe(Rho_0, Sp,  g*g); //Returns Rho_0 with a full initial frame filled with 0.2.
	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "NOSTRA" << endl;
	// The pspace to iterate over.
	vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

	cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" << endl;
	for (int i = 0; i < t_measure.size(); i++)
	{
		cout << t_measure[i] << " ";
	} 	cout << endl;

	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
	const int microToSeconds = 1000000;
	const double delay1 = 5 * microToSeconds;     //5 seconds

	cout << "Delay 1 in progress... (" << delay1 / microToSeconds << "s)" << endl;
	std::vector<vector <double>> vec;
	// Stores collated output from parallel method calls in proper ascending order of p values.

	size_t currentSize = getCurrentRSS(); //Check Memory usage.
	size_t peakSize = getPeakRSS();

	cout << "On initialisation, current size in MB: " << currentSize / (1024.0 * 1024.0) << " and Peak Size (in MB): " << peakSize / (1024.0 * 1024.0) << endl;

	auto start = high_resolution_clock::now();

	int nProcessors = omp_get_max_threads();
	if (nProcessors > 25)
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
		for (int i = 0; i < a_space.size(); i++)
		{

			//init_randconstframe(Rh0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
			//type="Gam";
			stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
			message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() << endl;
			cout << message.str();
			std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
			/**
			 * Namely CExpRho_a is structured as:
			 * | 	a		|    t 		|     <<Rho1(t)>>x,r			|    Var[<Rho1(t)>x],r    |
			**/

			RK4_Wrapper_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K, a_start, a_end, dt, dx, dP, r, g);
			//expanded_percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_space[i], b, c, D, sigma, dt, dx, r, g);
			//crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

			vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
			vector<vector <double>>().swap(CExpRho_a);
			// Remove dynamic vectors from memory and free up allocated space.
		}

		#pragma omp for schedule(static) ordered
		for (int i = 0; i < omp_get_num_threads(); i++)
		{
		#pragma omp ordered
			vec.insert(vec.end(), vec_private.begin(), vec_private.end());
			// Inserting critical exponent data for each grid size in order.
			stringstream message3;
			size_t currSize = getCurrentRSS(); //Check Memory usage.
			size_t peaksSize = getPeakRSS();

			message3 << "Is this happening?\n" << "In replicate " << i << ", current size in MB: " << currSize / (1024.0 * 1024.0)
				<< " and Peak Size (in MB): " << peaksSize / (1024.0 * 1024.0) << endl;
			cout << message3.str();

		}
		vector<vector <double>>().swap(vec_private);
		// Remove dynamic vectors from memory and free up allocated space.
	} //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, coco, tm, d3, p1, p2, rini, Dm0, Dm1, gm, alph, dix, dimitri;

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



	L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << a_start; p2 << a_end; rini << r; alph << alpha;
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
	output_1stdp.open("../Data/Rietkerk/Non-Stochastic/1stOrderCC_Rietkerk_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_D0_" + Dm0.str() + "_D1_" + Dm1.str() +
		"_a1_" + p1.str() + "_a2_" + p2.str() + "_dx_" + dix.str() + "_gmax_" + gm.str() + "_alpha_" + alph.str() + "_R_" + rini.str() + ".csv");

	cout << "Save file name: " << endl;
	cout << "../Data/Rietkerk/Non-Stochastic/1stOrderCC_Rietkerk_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_D0_" + Dm0.str() + "_D1_" + Dm1.str() +
		"_a1_" + p1.str() + "_a2_" + p2.str() + "_dx_" + dix.str() + "_gmax_" + gm.str() + "_alpha_" + alph.str() + rini.str() + ".csv";

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

	string header = " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs, # Active Sites";;
	output_1stdp << header << "\n";
	cout << "The vector elements are: " << endl;
	cout << header << "\n";

	for (int i = 0; i < vec.size(); i++)
	{
		for (int j = 0; j < vec[i].size() - 1; j++)
			output_1stdp << setprecision(16) << vec[i][j] << ",";
		output_1stdp << setprecision(16) << vec[i][vec[i].size() - 1] << endl;
		if (i % (1000) == 1)
		{
			for (int j = 0; j < vec[i].size() - 1; j++)
				cout << setprecision(16) << vec[i][j] << ",";
			cout << setprecision(16) << vec[i][vec[i].size() - 1] << endl;
		}
	}
	output_1stdp.close();
}


// ------------------------------------- Stochastic Integration Machinery ------------------------------------- //


// ALT INTEGRATION MACHINERY WHERE MIXING OF TIMESTEPS DUE TO LAPLACIAN MAY BE AN ISSUE.


void f_2Dor(D2Vec_Double& f, D2Vec_Double& Rho_M, D3Vec_Int& nR2, double a, double c, double gmax,
	double alpha, double rW, double W0, double(&Dxd2)[Sp], double(&K)[3], double t, double dt, double g)
{
	//Vector function that updates an array containing (dP/dt, dW/dt, DO/dt) for each site in the lattice.
	//Based on the Rietkerk model for plant vegetation dynamics with the Dornic twist where linear and stoch term for vegetation are already taken care of.
	// The following values are precomputed to reduce time complexity.
	double cgmax = c * gmax; double K2W0 = K[2] * W0;
	for (int i = 0; i < g * g; i++)
	{
		//Equations for the density of plants, soil water and surface water at each site.
		//Note that the Laplacian is calculated using reflective boundary conditions.
		//The Laplacian is calculated using the 5-point stencil method.
		/**
		 * NOTE!!!!!:  Equivalent to dP/dt = c*g_max*P*W/(W+K1)
		 *  Linear and stochastic terms taken care of by Dornic integration routine previously.
		**/

		f[0][i] = cgmax * Rho_M[1][i] * Rho_M[0][i] / (Rho_M[1][i] + K[1]);
		//Equivalent to dW/dt = alpha*(P+K2*W)/(P+K2) - rW*W + D*(Laplacian of W)
		f[1][i] = alpha * (Rho_M[0][i] + K[2] * W0) / (Rho_M[0][i] + K[2]) * Rho_M[2][i] - gmax * Rho_M[1][i] * Rho_M[0][i] / (Rho_M[1][i] + K[1]) - rW * Rho_M[1][i]
			+ (Dxd2[1]) * (Rho_M[1][nR2[i][0][0]] + Rho_M[1][nR2[i][0][1]] + Rho_M[1][nR2[i][1][0]] + Rho_M[1][nR2[i][1][1]] - 4 * Rho_M[1][i]);
		//Equivalent to dO/dt = a - alpha*(P+K*W)/(P+K)*O + D*(Laplacian of O)
		f[2][i] = a - alpha * (Rho_M[0][i] + K[2] * W0) / (Rho_M[0][i] + K[2]) * Rho_M[2][i] + (Dxd2[2]) * (Rho_M[2][nR2[i][0][0]] + Rho_M[2][nR2[i][0][1]]
			+ Rho_M[2][nR2[i][1][0]] + Rho_M[2][nR2[i][1][1]] - 4 * Rho_M[2][i]);

	}

}


void RK4_Integrate_Stochastic_1Sp(D2Vec_Double& Rho_t, D2Vec_Double& Rho_tsar, D2Vec_Double& K1, D2Vec_Double& K2, D2Vec_Double& K3, D2Vec_Double& K4,
	D3Vec_Int& nR2, double a, double c, double gmax, double alpha, double rW, double W0, double(&Dxd2)[Sp], double(&K)[3], double t, double dt, double dx, int g)
{

	/** // OLD DECLARATIONS of K1, K2, K3, K4 and Rho_M. Skipped as declaration and assignment is done in rietkerk_Dornic_2D_MultiSp
	 * and doing it here takes O(n^2) time each time this function is called.
	*/

	// RECALL: Dxd2[s] = D[s]/dx2
	f_2Dor(K1, Rho_t, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, t, dt, g); //K1 updated.

	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < g * g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K1[s][i];
	}

	f_2Dor(K2, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, t + dt2, dt, g); //K2 updated.

	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < g * g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K2[s][i];
	}
	f_2Dor(K3, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, t + dt2, dt, g); //K3 updated.
	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < g * g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt)*K3[s][i];
	}
	f_2Dor(K4, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, t + dt, dt, g); //K4 updated.

	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < g * g; i++)
		{
			Rho_t[s][i] += (dt6) * (K1[s][i] + 2.0 * K2[s][i] + 2.0 * K3[s][i] + K4[s][i]);

			if (Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || isnan(Rho_t[s][i]) == true)
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

void f_DorFKE_2Sp(D2Vec_Double& f, D2Vec_Double& Rho_M, D3Vec_Int& nR2, double a, double c, double gmax, double alpha, double rW, double W0,
	double(&Dxd2)[Sp], double(&K)[3], double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double t, double dt, double dx1_2, double g)
{
	//cout << "Placeholder" << endl;
	//Vector function that updates an array containing ( dP/dt, dW/dt, DO/dt,  dG/dt) for each site in the lattice.
	//Based on the Rietkerk model for plant vegetation dynamics with the Dornic twist where linear and stoch term for vegetation (and grazer) are already taken care of.
	// The following values are precomputed to reduce time complexity.
	double cgmax = c * gmax; double K2W0 = K[2] * W0; double A01H01 = A[0][1] * H[0][1]; double EA1 = E[1] * A[0][1];

	for (int i = 0; i < g * g; i++)
	{
		//Equations for the density of plants, soil water, surface water and grazers at each site.
		/** NOTE!!!!!:  Equivalent to dP/dt = c*g_max*P*W/(W+K1) - aij*hij*P*G/(1+aij*hij*V)
		 *  Linear and stochastic terms taken care of by Dornic integration routine previously.
		**/
		f[0][i] = cgmax * Rho_M[2][i] * Rho_M[0][i] / (Rho_M[2][i] + K[1]) - (A[0][1] * Rho_M[0][i] / (1 + A01H01 * Rho_M[0][i])) * Rho_M[1][i];
		//Equivalent to dG/dt = e*(aij*V*G)/(1+aij*hij*V) -mj*G  [NO PREDATION].
		f[1][i] = (EA1 * Rho_M[0][i] / (1 + A01H01 * Rho_M[0][i])) * Rho_M[1][i] - M[1] * Rho_M[1][i];
		//Equivalent to dW/dt = alpha*(P+K2*W0)/(P+K2)*O - rW*W + D*(Laplacian of W)
		f[2][i] = alpha * (Rho_M[0][i] + K[2] * W0) / (Rho_M[0][i] + K[2]) * Rho_M[3][i] - gmax * Rho_M[2][i] * Rho_M[0][i] / (Rho_M[2][i] + K[1]) - rW * Rho_M[2][i]
			+ (Dxd2[2]) * (Rho_M[2][nR2[i][0][0]] + Rho_M[2][nR2[i][0][1]] + Rho_M[2][nR2[i][1][0]] + Rho_M[2][nR2[i][1][1]] - 4 * Rho_M[2][i]);
		//Equivalent to dO/dt = a - alpha*(P+K*W)/(P+K)*O + D*(Laplacian of O)
		f[3][i] = a - alpha * (Rho_M[0][i] + K[2] * W0) / (Rho_M[0][i] + K[2]) * Rho_M[3][i] + (Dxd2[3]) * (Rho_M[3][nR2[i][0][0]] + Rho_M[3][nR2[i][0][1]]
			+ Rho_M[3][nR2[i][1][0]] + Rho_M[3][nR2[i][1][1]] - 4 * Rho_M[3][i]);
	}
}

void RK4_Integrate_Stochastic_2Sp(D2Vec_Double& Rho_t, D2Vec_Double& Rho_tsar, D2Vec_Double& K1, D2Vec_Double& K2, D2Vec_Double& K3, D2Vec_Double& K4, D3Vec_Int& nR2,
	double a, double c, double gmax, double alpha, double rW, double W0, double(&Dxd2)[Sp], double(&K)[3], double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double t, double dt, double dx, int g)
{
	/** // OLD DECLARATIONS of K1, K2, K3, K4 and Rho_M. Skipped as declaration and assignment is done in rietkerk_Dornic_2D_MultiSp
	 * and doing it here takes O(n^2) time each time this function is called.
	*/

	// RECALL: Dxd2[s] = D[s]/dx2
	f_DorFKE_2Sp(K1, Rho_t, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t, dt, dx1_2, g); //K1 updated.

	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < g * g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K1[s][i];
	}

	f_DorFKE_2Sp(K2, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t + dt2, dt, dx1_2, g); //K2 updated.

	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < g * g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K2[s][i];
	}
	f_DorFKE_2Sp(K3, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t + dt2, dt, dx1_2, g); //K3 updated.
	for (int s = 0; s < Sp; s++)
	{
		for (int i = 0; i < g * g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt)*K3[s][i];
	}
	f_DorFKE_2Sp(K4, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t + dt, dt, dx1_2, g); //K4 updated.

	for (int s = 0; s < Sp; s++)
	{
		bool isKoed = false;
		for (int i = 0; i < g * g; i++)
		{
			Rho_t[s][i] += (dt6) * (K1[s][i] + 2.0 * K2[s][i] + 2.0 * K3[s][i] + K4[s][i]);

			if (Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || isnan(Rho_t[s][i]) == true)
			{
				if (isKoed == false)
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
	}
}

// Calculate the gamma values for 2Species model with Vegetation and Grazers.
void calc_gamma_vel_2Sp_NonRefugia(const vector<pair<int, int>>& centralNeighboringSites, std::mt19937_64& rng, D2Vec_Double& Rho_t, D2Vec_Double& gamma,
	D2Vec_Pair_Double& v_eff, double(&v)[SpB], double(&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, double nVeg_frac, int r_max, int L)
{
	if (nVeg_frac < 0.35)
		nVeg_frac = 0.35; //Minimum fraction of vegetation cover.

	//int L2 = L*L;
	double eps = 1.0e-12; //Small number to avoid division by zero.

	int new_length = int(nVeg_frac * nVeg_frac * centralNeighboringSites.size());
	// Make a copy of the centralNeighboringSites, resized to the fraction nVeg_frac*nVeg_frac.
	auto range_CentralNeighboringSites = std::ranges::subrange(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);
	//vector<pair<int, int>> centralNeighboringSitesCopy(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);

	std::vector<int> nR_Perp(new_length);

	double rho_inverse[SpB]; //Inverse of average density.
	for (int s = 0; s < SpB; s++)
		Rho_avg[s] >= eps ? rho_inverse[s] = 1 / Rho_avg[s] : rho_inverse[s] = 1 / eps; //Avoid division by zero.
	// Avoid division by zero.

	for (int i = 0; i < L2; i++)
	{
		mutex errorMutex; // Mutex to make error messages thread-safe

		int c_i = int(i / L); int c_j = i % L; //Current x and y coordinates of site i.

		generateNeighboringSitesFromCentral(range_CentralNeighboringSites, nR_Perp, c_i, c_j, L);
		//generateNeighboringSitesFromCentral(centralNeighboringSites, nR_Perp, c_i, c_j, L);

		//if(nVeg_frac < 1)
		//	nR_Perp.resize(int(nVeg_frac*nVeg_frac*nR_Perp.size()));
		// In this case, r_effective is reduced to nVeg_frac*r_max, hence the reduction in the number of nearest neighbours.
		// rfrac is sorted in descending order.
		double fr_prev = 1.0; //Fraction of perception radius to max radius.
		for (const auto& frac : rfrac)
		{
			int s = frac.second;
			double fr = frac.first; //Fraction of perception radius to max radius.

			int eff_nR_Perp_size = int(fr * fr * nR_Perp.size());

			if (fr == 0.0)
			{
				continue;
			} // No perception radius for this species (indicates vegetation species)
			if (Rho_avg[s] < eps)
			{
				gamma[s][i] = 1.0; continue;
			} // No species left, value assigned doesn't matter.
			else if (Rho_avg[s - 1] < eps)
			{
				gamma[s][i] = 0.0;	continue;
			} // No resource left, consumers will advect to extinction.

// Gamma for grazers (indexed 1)
			for (int k = 0; k < eff_nR_Perp_size; k++)
			{
				gamma[s][i] += Rho_t[s - 1][nR_Perp[k]];
			}

			gamma[s][i] = (gamma[s][i] * rho_inverse[s - 1]) / (eff_nR_Perp_size);

			if (gamma[s][i] > 1)
			{
				gamma[s][i] = 1.0;
			} //Ensuring gamma is not greater than 1.
			else if (gamma[s][i] < 0)
			{
				gamma[s][i] = 0.0;
			} //Ensuring gamma is not less than 0.
			else if (isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
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

void f_DorFKE_3Sp(D2Vec_Double& f, D2Vec_Double& Rho_M, D3Vec_Int& nR2, double a, double c, double gmax, double alpha, double rW, double W0,
	double(&Dxd2)[Sp], double(&K)[3], double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double t, double dt, double dx1_2, double g)
{
	//Vector function that updates an array containing ( dP/dt, dW/dt, DO/dt,  dG/dt, dPr/dt) for each site in the lattice.
	//Based on the Rietkerk model for plant vegetation dynamics with the Dornic twist where linear and stoch term for vegetation (and grazer) are already taken care of.

	// The following values are precomputed to reduce time complexity.
	double cgmax = c * gmax; double K2W0 = K[2] * W0; double A01H01 = A[0][1] * H[0][1]; double A12H12 = A[1][2] * H[1][2];
	// The following are the values of E*A and A*H for each species, which reduces time complexity by avoiding repeated calculations.
	double EA1 = E[1] * A[0][1]; double EA2 = E[2] * A[1][2]; int L2 = g * g; //L2 is the total number of sites in the lattice.
	//vector<double> EA(SpB);
	//for(int s=1; s< SpB; s++)
	//	EA[s] = E[s]*A[s-1][s];

#pragma omp simd
	for (int i = 0; i < L2; i++)
	{
		//Equations for the density of plants, soil water, surface water and grazers at each site.
		//Note that the Laplacian is calculated using reflective boundary conditions.
		//The Laplacian is calculated using the 5-point stencil method.
		/**
		 * NOTE!!!!!:  Equivalent to dP/dt = c*g_max*P*W/(W+K1) - aij*hij*P*G/(1+aij*hij*V)
		 *  Linear and stochastic terms taken care of by Dornic integration routine previously.
		**/

		// RECALL: Dxd2[s] = D[s]/dx2 AND K_P1 = Inverse of Carrying Capacity of Predator (defined as Km in main())
		f[0][i] = cgmax * Rho_M[3][i] * Rho_M[0][i] / (Rho_M[3][i] + K[1]) - (A[0][1] * Rho_M[0][i] / (1 + A01H01 * Rho_M[0][i])) * Rho_M[1][i];
		//Equivalent to dG/dt = ej*(aij*V*P)/(1+aij*hij*V)  -ajm*G*P/(1+ajm*hjm*G) -Mj*G
		f[1][i] = (EA1 * Rho_M[0][i] / (1 + A01H01 * Rho_M[0][i])) * Rho_M[1][i] - M[1] * Rho_M[1][i]
			- A[1][2] * Rho_M[1][i] / (1 + A12H12 * Rho_M[1][i]) * Rho_M[2][i];
		//Equivalent to dPr/dt = em*(ajm*G*Pr)/(1+ajm*hjm*G) - Mm*Pr
		f[2][i] = (EA2 * Rho_M[1][i] / (1 + A12H12 * Rho_M[1][i])) * Rho_M[2][i] - M[2] * Rho_M[2][i];
		//Equivalent to dW/dt = alpha*(P+K2*W0)/(P+K2)*O - g_max*P*W/(W+K1) - rW*W + D*(Laplacian of W)
		f[3][i] = alpha * (Rho_M[0][i] + K2W0) / (Rho_M[0][i] + K[2]) * Rho_M[4][i] - gmax * Rho_M[3][i] * Rho_M[0][i] / (Rho_M[3][i] + K[1]) - rW * Rho_M[3][i]
			+ (Dxd2[3]) * (Rho_M[3][nR2[i][0][0]] + Rho_M[3][nR2[i][0][1]] + Rho_M[3][nR2[i][1][0]] + Rho_M[3][nR2[i][1][1]] - 4 * Rho_M[3][i]);
		//Equivalent to dO/dt = a - alpha*(P+K2*W)/(P+K2)*O + D*(Laplacian of O)
		f[4][i] = a - alpha * (Rho_M[0][i] + K2W0) / (Rho_M[0][i] + K[2]) * Rho_M[4][i] + (Dxd2[4]) * (Rho_M[4][nR2[i][0][0]] + Rho_M[4][nR2[i][0][1]]
			+ Rho_M[4][nR2[i][1][0]] + Rho_M[4][nR2[i][1][1]] - 4 * Rho_M[4][i]);

	}
}

void RK4_Integrate_Stochastic_MultiSp(D2Vec_Double& Rho_t, D2Vec_Double& Rho_tsar, D2Vec_Double& K1, D2Vec_Double& K2, D2Vec_Double& K3, D2Vec_Double& K4, D3Vec_Int& nR2,
	double a, double c, double gmax, double alpha, double rW, double W0, double(&Dxd2)[Sp], double(&K)[3], double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double t, double dt, double dx, int g)
{
	//double dt6 = dt/6.0; double dt2 = dt/2.0; double dx1_2 = 1/(dx*dx);

	/** // OLD DECLARATIONS of K1, K2, K3, K4 and Rho_M. Skipped as declaration and assignment is done in rietkerk_Dornic_2D_MultiSp
	 * and doing it here takes O(n^2) time each time this function is called.
	D2Vec_Double K1(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K2(Sp, vector<double> (g*g, 0.0));
	D2Vec_Double K3(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K4(Sp, vector<double> (g*g, 0.0));
	D2Vec_Double Rho_M(Sp, vector<double> (g*g, 0.0));
	*/

	// RECALL: Dxd2[s] = D[s]/dx2
	#if SPB == 3
	f_DorFKE_3Sp(K1, Rho_t, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t, dt, dx1_2, g); //K1 updated.
	#elif SPB == 2
	f_DorFKE_2Sp(K1, Rho_t, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t, dt, dx1_2, g); //K1 updated.
	#endif

	for (int s = 0; s < Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for (int i = 0; i < g * g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K1[s][i];
	}

	#if SPB == 3
	f_DorFKE_3Sp(K2, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t + dt2, dt, dx1_2, g); //K2 updated.
	#elif SPB == 2
	f_DorFKE_2Sp(K2, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t + dt2, dt, dx1_2, g); //K2 updated.
	#endif

	for (int s = 0; s < Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for (int i = 0; i < g * g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt2)*K2[s][i];
	}
	#if SPB == 3
	f_DorFKE_3Sp(K3, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t + dt2, dt, dx1_2, g); //K3 updated.
	#elif SPB == 2
	f_DorFKE_2Sp(K3, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t + dt2, dt, dx1_2, g); //K3 updated.
	#endif
	for (int s = 0; s < Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for (int i = 0; i < g * g; i++)
			Rho_tsar[s][i] = Rho_t[s][i] + (dt)*K3[s][i];
	}
	#if SPB == 3
	f_DorFKE_3Sp(K4, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t + dt, dt, dx1_2, g); //K4 updated.
	#elif SPB == 2
	f_DorFKE_2Sp(K4, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, Dxd2, K, A, H, E, M, t + dt, dt, dx1_2, g); //K4 updated.
	#endif

	for (int s = 0; s < Sp; s++)
	{
		#pragma omp simd // SIMD optimization for better performance.
		for (int i = 0; i < g * g; i++)
		{
			Rho_t[s][i] += (dt6) * (K1[s][i] + 2.0 * K2[s][i] + 2.0 * K3[s][i] + K4[s][i]);
		}
		#ifdef DEBUG
		bool isKoed = false; // Flag to indicate if the RK4 integration was "KO'ed" due to invalid values.
		for (int i = 0; i < g * g; i++)
		{
			if (Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || isnan(Rho_t[s][i]) == true)
			{
				if (isKoed == false)
				{
					//std::ofstream errout; //Save Error Logs
					//std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
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

void save_prelimframe(D2Vec_Double& Rho_t, const string& parendir, const string& filenamePattern, double a, double a_st, double a_end, double t,
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

		m1 << "File exists. New filename: " << filename << " \n and maxRepNo: " << maxRepNo << endl;
		cout << m1.str();
		basefilename = filenamePattern + std::to_string(maxRepNo + 1) + ".csv";
	}



	std::ostringstream oss;
	if (header != "")
		oss << header;
	else
		oss << frame_header;
	// Using default header (defined in header file) if no header is provided.

//tot_iter is the number of rows in the Rho_t matrix.
	double tot_rows = Rho_t.size();
	for (int i = 0; i < tot_rows; i++)
	{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |
		oss << a << "," << r + 1 << "," << g << "," << Rho_t[i][0] << "," << Rho_t[i][1] << "," << Rho_t[i][2] << ",";
		for (int s = 0; s < SpB; s++)
		{
			oss << Rho_t[i][4 * s + 3] << "," << Rho_t[i][4 * s + 4] << "," << Rho_t[i][4 * s + 5] << "," << Rho_t[i][4 * s + 6] << ",";
		}
		oss << "\n";
	}

	FILE* fp = fopen(filename.c_str(), "w");
	if (fp)
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

	if (r > 2 && delete_previous)
	{
		// Delete the previous file (if it exists) to save space.
		// Reset rini value.
		stringstream rini_prev; rini_prev << r - 2;
		// Get substring from the string filenamePattern starting after the sequence "_DP_G_" to the end of the string.
		std::smatch m;
		std::regex_search(filenamePattern, m, std::regex("_DP_G_"));
		string filenamePattern_prev = prelim_prefix + ".*_DP_G_" + m.suffix().str() + rini_prev.str() + ".csv";
		//stringstream m7; m7 << "For (r,a,t)" << r << " , " << a << " , " << t 
		//<< " Previous filename pattern: " << filenamePattern_prev << "with m: " << m.str() << "\n"; cout << m7.str();

		// Using regex to delete the file with the previous replicate number.
		std::regex filePattern(filenamePattern_prev);

		// Iterate over the files in the parent directory (path_to_dir) and delete the file if it matches the pattern.
		for (const auto& entry : fs::directory_iterator(parendir))
		{
			string matchStr = "/" + entry.path().filename().string();
			if (std::regex_search(matchStr, filePattern))
			{
				if (fs::exists(entry.path()) && fs::is_regular_file(entry.path()))
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

void save_frame(D2Vec_Double& Rho_t, const string& parendir, const string& filenamePattern, double a, double a_st, double a_end, double t, double dt, double dx, double dP, int r, int g, string header /* =""*/, bool overwrite /* = false*/)
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
	if (header != "")
		oss << header;
	else
		oss << frame_header;
	// Using default header (defined in header file) if no header is provided.

	double Sp = Rho_t.size();
	for (int i = 0; i < g * g; i++)
	{
		oss << a << "," << i;
		for (int s = 0; s < Sp; s++)
			oss << "," << Rho_t[s][i];
		oss << "\n";
	}

	FILE* fp = fopen(filename.c_str(), "w");
	if (fp)
	{
		fprintf(fp, "%s", oss.str().c_str());
		fclose(fp);
	}
	else
	{
		std::cerr << "Error: Could not open file " << filename << std::endl;
	}

}


void save_framefileswrapper(int index, int tot_iter, int j, int thrID, double t, double dt, double dx, D2Vec_Double& Rho_dt, D2Vec_Double& DRho,
	D2Vec_Double& rho_rep_avg_var, vector<double>& t_meas, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, double t_max,
	double a, double c, double gmax, double alpha, double rW, double W0, double(&D)[Sp], double(&v)[SpB],
	double(&K)[3], double sigma[], double a_st, double a_end, double a_c, double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB],
	double dP, int r, int g, double Gstar /* =-1.*/, double Vstar /* = -1.*/)
{

		//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.

		stringstream L, tm, d3, p1, rini, gm, a1, a2, Dm0, Dm1, Dm2, alph, w0t, aij, hij, dix, dimitri, sig0, sig1, geq, veq, thr_ID;

		L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
		rini << j; Dm0 << D[0] * pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2]; thr_ID << thrID;
		a1 << a_st; a2 << a_end; sig0 << sigma[0]; sig1 << sigma[1]; dimitri << dP; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;
		//gm << setprecision(3) << gmax; w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; 
		// Three replicates are over.

		string parendir = "";

		// Creating a file instance called output to store output data as CSV.
		if (Gstar != -1)
			parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Geq_" + geq.str();
		else if (Vstar != -1)
			parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Veq_" + veq.str();
		else
			parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str();

		string filenamePattern = frame_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
			+ "_a_" + p1.str() + "_D1_" + Dm1.str() + "_jID_" + thr_ID.str() + "_dx_" + dix.str() + "_R_";
		//Next, designate the naive filename
		string filename = filenamePattern + rini.str() + ".csv";

		#if SPB > 1
		string gammafilenamePattern = gamma_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
			+ "_a_" + p1.str() + "_D1_" + Dm1.str() + "_jID_" + thr_ID.str() + "_dx_" + dix.str() + "_R_";
		#endif

		///** // SAVE ALL FRAMES
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
		if (t > 1500 && t < 6000)
		{	//Only save one in two frames here.
			if (index % 2 == 0)
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
int save_prelimfileswrapper(int index, int tot_iter, int j, int thrID, double t, double dt, double dx, D2Vec_Double& Rho_dt, D2Vec_Double& DRho, D2Vec_Double& Rho_M, D2Vec_Double& Rho_Mov,
	D2Vec_Double& rho_rep_avg_var, vector<double>& t_meas, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, double t_max,
	double a, double c, double gmax, double alpha, double rW, double W0, double(&D)[Sp], double(&v)[SpB],
	double(&K)[3], double sigma[], double a_st, double a_end, double a_c, double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB],
	double dP, int r, int g, double Gstar /* =-1.*/, double Vstar /* = -1.*/)
{
	// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
	// and save to file.
	D2Vec_Double rho_rep_avg_var_temp(index + 1, vector<double>(Sp4_1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
	//std::copy(rho_rep_avg_var.begin(), rho_rep_avg_var.begin() + index, rho_rep_avg_var_temp.begin());
	var_mean_incremental_surv_runs(rho_rep_avg_var_temp, Rho_M, index + 1, 0);

	#if SPB > 1 // Save movement time-series for each species in Rho_Mov.
	D2Vec_Double gamma_mov_avg_temp(index + 1, vector<double>(int(Rho_Mov[0].size()) + 1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
	// Get cols_per_sp by dividing the number of columns in Rho_Mov by the number of species (SpB).
	int cols_per_sp = int(Rho_Mov[0].size() / SpB);
	generic_SPBmean_surv_runs(gamma_mov_avg_temp, Rho_Mov, index + 1, cols_per_sp, 0);
	//Time is stored in the first column.
	for (int i = 0; i < index + 1; i++)
		gamma_mov_avg_temp[i][0] = t_meas[i];

	#endif


	//Time is stored in the first column.
	for (int i = 0; i < index + 1; i++)
		rho_rep_avg_var_temp[i][0] = t_meas[i];

	//Finally save to file.
	stringstream L, tm, d3, p1, a1, a2, dimitri, rini, Dm, geq, veq, jID; // cgm, sig0;
	a1 << a_st; a2 << a_end;
	L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(5) << a; dimitri << dP; jID << thrID;
	rini << j; Dm << setprecision(4) << D[2]; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;

	string parendir = "";
	if (Gstar != -1)
		parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Geq_" + geq.str() + "/TimeSeries";
	else if (Vstar != -1)
		parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Veq_" + veq.str() + "/TimeSeries";
	else
		parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "/TimeSeries";

	//double ran_jid = (unif(rng)*1000.0)/1000.0; jID << ran_jid; // Random number between 0 and 1.

	string filenamePattern = replicate_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
		"_jID_" + jID.str() + "_R_";

	save_prelimframe(rho_rep_avg_var_temp, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);

	#if SPB > 1 // Save movement time-series for each species in Rho_Mov.
	string filenamePattern_mov = movement_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
		"_jID_" + jID.str() + "_R_";
	save_prelimframe(gamma_mov_avg_temp, parendir, filenamePattern_mov, a, a_st, a_end, t, dt, dx, dP, j, g, movmheader, false, false);
	#endif

	// Check if no vegetation is left at this current index. If so, break out of the while loop.
	if (Rho_M[index][0] == 0.0)
	{
		stringstream m3_1;
		m3_1 << "RUN-TIME WARNING: ZERO Active VEG sites for TIME [t, a, thr, j]\t"
			<< t << " , " << a << " , " << omp_get_thread_num() << " , " << j << " Skipping to next replicate .... \n";
		cout << m3_1.str(); cerr << m3_1.str();

		//Next save future preliminary frames with 0.0 values for population densities after this time point.
		int index_prelim_range[4] = { int(tot_iter * 0.85), int(tot_iter * 0.9), int(tot_iter * 0.95), int(tot_iter - 1) };
		int iter_index = 0;
		for (int iter_index = 0; iter_index < 4; iter_index++)
		{
			int new_index = index_prelim_range[iter_index];
			if (index > new_index)
				continue; //Skip to next iteration if index is greater than the current range.
			// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
			// and save to file.
			D2Vec_Double rho_rep_avg_var_future_temp(new_index + 1, vector<double>(Sp4_1, 0.0));
			//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
			var_mean_incremental_surv_runs(rho_rep_avg_var_future_temp, Rho_M, new_index + 1, 0);
			stringstream tm_future; tm_future << t_meas[new_index]; // Time at which the frame is saved.

			//Time is stored in the first column.
			for (int i = 0; i < new_index + 1; i++)
				rho_rep_avg_var_future_temp[i][0] = t_meas[i];

			string filenamePattern_future = replicate_prefix + L.str() + "_T_" + tm_future.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
				"_jID_" + jID.str() + "_R_";

			//Finally save to file.
			save_prelimframe(rho_rep_avg_var_future_temp, parendir, filenamePattern_future, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);

			#if SPB > 1 // Save movement time-series for each species in Rho_Mov.
			D2Vec_Double gamma_mov_avg_future_temp(new_index + 1, vector<double>(int(Rho_Mov[0].size()) + 1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
			// Get cols_per_sp by dividing the number of columns in Rho_Mov by the number of species (SpB).
			int cols_per_sp = int(Rho_Mov[0].size() / SpB);
			generic_SPBmean_surv_runs(gamma_mov_avg_future_temp, Rho_Mov, new_index + 1, cols_per_sp, 0);
			//Time is stored in the first column.
			for (int i = 0; i < new_index + 1; i++)
				gamma_mov_avg_future_temp[i][0] = t_meas[i];

			string filenamePattern_mov_future = movement_prefix + L.str() + "_T_" + tm_future.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
				"_jID_" + jID.str() + "_R_";
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


// A wrapper function to save the files for each thread, called from the main simulation loop in rietkerk_Dornic_2D_MultiSp()
int save_fileswrapper(int index, int tot_iter, int j, int thrID, double t, double dt, double dx, D2Vec_Double& Rho_dt, D2Vec_Double& DRho, D2Vec_Double& Rho_M, D2Vec_Double& Rho_Mov,
	D2Vec_Double& rho_rep_avg_var, vector<double>& t_meas, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, double t_max, double a, double c, double gmax, double alpha, double rW, double W0,
	double(&D)[Sp], double(&v)[SpB], double(&K)[3], double sigma[], double a_st, double a_end, double a_c, double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB],
	double dP, int r, int g, double Gstar /* =-1.*/, double Vstar /* = -1.*/)
{

	// FRAME SAVING
	if (index >= tot_iter - 10 && index <= tot_iter - 1 || t >= 60000 && t <= 150000 || t >= 50 && t <= 12000
		|| t == 0)
	{
		//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.

		stringstream L, tm, d3, p1, rini, gm, a1, a2, Dm0, Dm1, Dm2, alph, w0t, aij, hij, dix, dimitri, sig0, sig1, geq, veq, thr_ID;

		L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
		rini << j; Dm0 << D[0] * pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2]; thr_ID << thrID;
		a1 << a_st; a2 << a_end; sig0 << sigma[0]; sig1 << sigma[1]; dimitri << dP; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;
		//gm << setprecision(3) << gmax; w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; 
		// Three replicates are over.

		string filenamePattern = frame_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
			+ "_a_" + p1.str() + "_D1_" + Dm1.str() + "_jID_" + thr_ID.str() + "_dx_" + dix.str() + "_R_";
		//Next, designate the naive filename
		string filename = filenamePattern + rini.str() + ".csv";
		string parendir = "";

		// Creating a file instance called output to store output data as CSV.
		if (Gstar != -1)
			parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Geq_" + geq.str();
		else if (Vstar != -1)
			parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Veq_" + veq.str();
		else
			parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str();

		/** // SAVE ALL FRAMES
		save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);

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
				stringstream m3;
				m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
				cout << m3.str();
			}
		} //*/

		///** PREVIOUS FRAME SAVING STRATEGY (DT ~ 0.1)
		if (t < 760)
		{
			//Only save one in three frames here.
			if (index % 3 == 2 || t == 0)
			{
				save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);

				stringstream m3;
				m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
				cout << m3.str();
			}
		}
		else if (t >= 760 && t < 6000)
		{	//Only save one in two frames here.
			if (index % 2 == 0)
			{
				save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
				stringstream m3;
				m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
				cout << m3.str();
			}
		}//*/
		else
		{
			//Save all frames here.
			save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
			stringstream m3;
			m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
			cout << m3.str();
		}
		//*/



	}
	//

	// BLOCK FOR CALCULATING AND TEMP PRELIMINARY FRAMES
	if (index == int(tot_iter * 0.85) || index == int(tot_iter * 0.9) || index == int(tot_iter * 0.95) || index == int(tot_iter - 1))
	{
		// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
		// and save to file.
		D2Vec_Double rho_rep_avg_var_temp(index + 1, vector<double>(Sp4_1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
		//std::copy(rho_rep_avg_var.begin(), rho_rep_avg_var.begin() + index, rho_rep_avg_var_temp.begin());
		var_mean_incremental_surv_runs(rho_rep_avg_var_temp, Rho_M, index + 1, 0);

		#if SPB > 1 // Save movement time-series for each species in Rho_Mov.
		D2Vec_Double gamma_mov_avg_temp(index + 1, vector<double>(int(Rho_Mov[0].size()) + 1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
		// Get cols_per_sp by dividing the number of columns in Rho_Mov by the number of species (SpB).
		int cols_per_sp = int(Rho_Mov[0].size() / SpB);
		generic_SPBmean_surv_runs(gamma_mov_avg_temp, Rho_Mov, index + 1, cols_per_sp, 0);

		//Time is stored in the first column.
		for (int i = 0; i < index + 1; i++)
			gamma_mov_avg_temp[i][0] = t_meas[i];

		#endif
		//Time is stored in the first column.
		for (int i = 0; i < index + 1; i++)
			rho_rep_avg_var_temp[i][0] = t_meas[i];

		//Finally save to file.
		stringstream L, tm, d3, p1, a1, a2, dimitri, rini, Dm, geq, veq, jID; // cgm, sig0;
		a1 << a_st; a2 << a_end;
		L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(5) << a; dimitri << dP;
		rini << j; Dm << setprecision(4) << D[2]; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;

		string parendir = "";
		if (Gstar != -1)
			parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Geq_" + geq.str() + "/TimeSeries";
		else if (Vstar != -1)
			parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Veq_" + veq.str() + "/TimeSeries";
		else
			parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "/TimeSeries";

		//double ran_jid = (unif(rng)*1000.0)/1000.0; jID << ran_jid; // Random number between 0 and 1.

		string filenamePattern = replicate_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
			"_jID_" + jID.str() + "_R_";

		save_prelimframe(rho_rep_avg_var_temp, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);

		#if SPB > 1 // Save movement time-series for each species in Rho_Mov.
		string filenamePattern_mov = movement_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
			"_jID_" + jID.str() + "_R_";
		save_prelimframe(gamma_mov_avg_temp, parendir, filenamePattern_mov, a, a_st, a_end, t, dt, dx, dP, j, g, movmheader, false, false);
		#endif


		// Check if no vegetation is left at this current index. If so, break out of the while loop.
		if (Rho_M[index][0] == 0.0)
		{
			stringstream m3_1;
			m3_1 << "RUN-TIME WARNING: ZERO Active VEG sites for TIME [t, a, thr, j]\t"
				<< t << " , " << a << " , " << omp_get_thread_num() << " , " << j << " Skipping to next replicate .... \n";
			cout << m3_1.str(); cerr << m3_1.str();

			//Next save future preliminary frames with 0.0 values for population densities after this time point.
			int index_prelim_range[4] = { int(tot_iter * 0.85), int(tot_iter * 0.9), int(tot_iter * 0.95), int(tot_iter - 1) };
			int iter_index = 0;
			for (int iter_index = 0; iter_index < 4; iter_index++)
			{
				int new_index = index_prelim_range[iter_index];
				if (index > new_index)
					continue; //Skip to next iteration if index is greater than the current range.
				// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
				// and save to file.
				D2Vec_Double rho_rep_avg_var_future_temp(new_index + 1, vector<double>(Sp4_1, 0.0));
				//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
				var_mean_incremental_surv_runs(rho_rep_avg_var_future_temp, Rho_M, new_index + 1, 0);
				stringstream tm_future; tm_future << t_meas[new_index]; // Time at which the frame is saved.

				//Time is stored in the first column.
				for (int i = 0; i < new_index + 1; i++)
					rho_rep_avg_var_future_temp[i][0] = t_meas[i];

				string filenamePattern_future = replicate_prefix + L.str() + "_T_" + tm_future.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
					"_D2_" + Dm.str() + "_R_";

				//Finally save to file.
				save_prelimframe(rho_rep_avg_var_future_temp, parendir, filenamePattern_future, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);
				vector<vector<double>>().swap(rho_rep_avg_var_future_temp); //Flush temp out of memory.

				#if SPB > 1 // Save movement time-series for each species in Rho_Mov.
				D2Vec_Double gamma_mov_avg_future_temp(new_index + 1, vector<double>(int(Rho_Mov[0].size()) + 1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
				// Get cols_per_sp by dividing the number of columns in Rho_Mov by the number of species (SpB).
				int cols_per_sp = int(Rho_Mov[0].size() / SpB);
				generic_SPBmean_surv_runs(gamma_mov_avg_future_temp, Rho_Mov, new_index + 1, cols_per_sp, 0);
				//Time is stored in the first column.
				for (int i = 0; i < new_index + 1; i++)
					gamma_mov_avg_future_temp[i][0] = t_meas[i];

				string filenamePattern_mov_future = movement_prefix + L.str() + "_T_" + tm_future.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
					"_jID_" + jID.str() + "_R_";
				save_prelimframe(gamma_mov_avg_future_temp, parendir, filenamePattern_mov_future, a, a_st, a_end, t, dt, dx, dP, j, g, movmheader, false, false);
				vector<vector<double>>().swap(gamma_mov_avg_future_temp); //Flush temp out of memory.
				#endif
			}
			return 1; //Indicate that the current replicate has no vegetation left at this time point (and break out of the while loop).
		}

		vector<vector<double>>().swap(rho_rep_avg_var_temp); //Flush temp out of memory.
	}
	// */

	return 0;


}


//Calculates gamma and velocities for 3 Sp Rietkerk model (details in PDF, 
// ASSUMING PREDATORS AREN'T DETERRED BY HIGH LOCAL VEGETATION DENSITY)
void calc_gamma_vel_3Sp_NonRefugia(const vector<pair<int, int>>& centralNeighboringSites, std::mt19937_64& rng, D2Vec_Double& Rho_t, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, double(&v)[SpB],
	double(&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, vector <std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L)
{
	std::uniform_real_distribution<double> unif(0.0, 1.0);

	if (nVeg_frac < 0.35)
		nVeg_frac = 0.35; //Minimum fraction of vegetation cover.

	//int L2 = L*L;
	double eps = epsilon; //Small number to avoid division by zero.

	int new_length = int(nVeg_frac * nVeg_frac * centralNeighboringSites.size());
	new_length < 1 ? new_length = 1 : new_length = new_length; // At least one nearest neighbor (the site itself) is considered.
	// Make a copy of the centralNeighboringSites, resized to the fraction nVeg_frac*nVeg_frac.
	auto range_CentralNeighboringSites = std::ranges::subrange(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);
	//vector<pair<int, int>> centralNeighboringSitesCopy(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);

	std::vector<int> nR_Perp(new_length);

	double rho_inverse[SpB];// ={1/eps, 1/eps, 1/eps}; //Inverse of average density.
	for (int s = 0; s < SpB; s++)
		Rho_avg[s] >= eps ? rho_inverse[s] = 1 / Rho_avg[s] : rho_inverse[s] = 1 / eps; //Avoid division by zero.
	// Avoid division by zero.

	for (int i = 0; i < L2; i++)
	{
		mutex errorMutex; // Mutex to make error messages thread-safe

		int c_i = int(i / L); int c_j = i % L; //Current x and y coordinates of site i.

		generateNeighboringSitesFromCentral(range_CentralNeighboringSites, nR_Perp, c_i, c_j, L);

		//if(nVeg_frac < 1)
		//	nR_Perp.resize(int(nVeg_frac*nVeg_frac*nR_Perp.size()));
		// In this case, r_effective is reduced to nVeg_frac*r_max, hence the reduction in the number of nearest neighbours.
		// rfrac is sorted in descending order.
		double fr_prev = 1.0; //Fraction of perception radius to max radius.
		for (const auto& frac : rfrac)
		{
			int s = frac.second;
			double fr = frac.first; //Fraction of perception radius to max radius.

			int eff_nR_Perp_size = int(fr * fr * nR_Perp.size());

			/**
			if(dtV_counter[s].first != 1 && dtV_counter[s].second != 0)
			{	continue;	}
			//*/
			//Skip this species if it is not updated at this time step (depending on advection timescale)

			//Next, update the effective velocity for the species.
			double ran = unif(rng);
			//double vx = v[s]*cos(2*PI*ran); double vy = v[s]*sin(2*PI*ran); //Advection vector.
			v_eff[s][i].first = v[s] * cos(2 * PI * ran); v_eff[s][i].second = v[s] * sin(2 * PI * ran); //Advection vector (vx, vy).

			if (fr == 0.0)
			{
				continue;
			} // No perception radius for this species (indicates vegetation species)
			if (eff_nR_Perp_size < 1 && fr > 0)
			{
				eff_nR_Perp_size = 1;
			} // At least one nearest neighbor (the site itself) is considered.
			if (Rho_avg[s] < eps)
			{
				gamma[s][i] = 1.0; continue;
			} // No species left, value assigned doesn't matter.
			else if (Rho_avg[s - 1] < eps)
			{
				gamma[s][i] = 0.0;	continue;
			} // No resource left, consumers will advect to extinction.

			// Gamma for grazers (indexed 1)
			double gamma_sum = 0.0; //Sum of gamma values for this species at site i (needed for SIMD reduction).
			if (s == 1)
			{

				if (Rho_avg[2] < eps)
				{
					#pragma omp simd reduction(+:gamma_sum )
					for (int k = 0; k < eff_nR_Perp_size; k++)
					{
						gamma_sum += Rho_t[s - 1][nR_Perp[k]];
					}
				} // No predators left
				else
				{
					//double rho_pred1 = 1/Rho_avg[2];
					#pragma omp simd reduction(+:gamma_sum )
					for (int k = 0; k < eff_nR_Perp_size; k++)
					{
						gamma_sum += Rho_t[s - 1][nR_Perp[k]] * (1 - Rho_t[2][nR_Perp[k]] * rho_inverse[2]);
					}
				}
				// RECALL : rho_inverse[s] = 1/Rho_avg[s]
				gamma[s][i] = (gamma_sum * rho_inverse[s - 1]) / (eff_nR_Perp_size);

			}
			if (s == 2)
			{
				#pragma omp simd reduction(+:gamma_sum )
				for (int k = 0; k < eff_nR_Perp_size; k++)
				{
					gamma_sum += Rho_t[s - 1][nR_Perp[k]];
				}
				// RECALL : rho_inverse[s] = 1/Rho_avg[s]
				gamma[s][i] = (gamma_sum / (eff_nR_Perp_size)) * rho_inverse[s - 1];
			}

			if (gamma[s][i] > 1)
			{
				gamma[s][i] = 1.0;
			} //Ensuring gamma is not greater than 1.
			else if (gamma[s][i] < 0)
			{
				gamma[s][i] = 0.0;
			} //Ensuring gamma is not less than 0.
			else if (isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
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
void calc_gamma_vel_3Sp(const vector<pair<int, int>>& centralNeighboringSites, std::mt19937_64& rng, D2Vec_Double& Rho_t, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, double(&v)[SpB],
	double(&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, vector <std::pair<int, int>>& dtV_counter, double nVeg_frac, int r_max, int L)
{
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	if (nVeg_frac < 0.35)
		nVeg_frac = 0.35; //Minimum fraction of vegetation cover.

	int new_length = int(nVeg_frac * nVeg_frac * centralNeighboringSites.size());
	new_length < 1 ? new_length = 1 : new_length = new_length; // At least one nearest neighbor (the site itself) is considered.

	auto range_CentralNeighboringSites = std::ranges::subrange(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);

	// Make a copy of the centralNeighboringSites, resized to the fraction nVeg_frac*nVeg_frac.
	//vector<pair<int, int>> centralNeighboringSitesCopy(centralNeighboringSites.begin(), centralNeighboringSites.begin() + new_length);

	std::vector<int> nR_Perp(new_length);


	//int L2 = L*L;
	double eps = 1.0e-12; //Small number to avoid division by zero.

	double rho_inverse[SpB];// ={1/eps, 1/eps, 1/eps}; //Inverse of average density.
	for (int s = 0; s < SpB; s++)
		Rho_avg[s] >= eps ? rho_inverse[s] = 1 / Rho_avg[s] : rho_inverse[s] = 1 / eps; //Avoid division by zero.
	// Avoid division by zero.
	for (int i = 0; i < L2; i++)
	{
		mutex errorMutex; // Mutex to make error messages thread-safe

		int c_i = int(i / L); int c_j = i % L; //Current x and y coordinates of site i.
		generateNeighboringSitesFromCentral(range_CentralNeighboringSites, nR_Perp, c_i, c_j, L);

		//if(nVeg_frac < 1)
		//	nR_Perp.resize(int(nVeg_frac*nVeg_frac*nR_Perp.size()));
		// In this case, r_effective is reduced to nVeg_frac*r_max, hence the reduction in the number of nearest neighbours.
		// rfrac is sorted in descending order.
		double fr_prev = 1.0; //Fraction of perception radius to max radius.
		for (const auto& frac : rfrac)
		{
			int s = frac.second;
			double fr = frac.first; //Fraction of perception radius to max radius.

			int eff_nR_Perp_size = int(fr * fr * nR_Perp.size());
			// If fr < 1, reduce the number of nearest neighbors to the fraction of the original number of nearest neighbors.

			/**
			if(dtV_counter[s].first != 1 && dtV_counter[s].second != 0)
			{	continue;	}
			//*/
			//Skip this species if it is not updated at this time step (depending on advection timescale)

			//Next, update the effective velocity for the species.
			double ran = unif(rng);
			//double vx = v[s]*cos(2*PI*ran); double vy = v[s]*sin(2*PI*ran); //Advection vector.
			v_eff[s][i].first = v[s] * cos(2 * PI * ran); v_eff[s][i].second = v[s] * sin(2 * PI * ran); //Advection vector (vx, vy).

			if (fr == 0.0)
			{
				continue;
			} // No perception radius for this species (indicates vegetation species)
			if (eff_nR_Perp_size < 1 && fr > 0)
			{
				eff_nR_Perp_size = 1;
			} // At least one nearest neighbor (the site itself) is considered.
			if (Rho_avg[s] < eps)
			{
				gamma[s][i] = 1.0; continue;
			} // No species left, value assigned doesn't matter.
			else if (Rho_avg[s - 1] < eps)
			{
				gamma[s][i] = 0.0;	continue;
			} // No resource left, consumers will advect to extinction.

// Gamma for grazers (indexed 1)
			if (s == 1)
			{
				if (Rho_avg[2] < eps)
				{
					for (int k = 0; k < eff_nR_Perp_size; k++)
					{
						gamma[s][i] += Rho_t[s - 1][nR_Perp[k]];
					}
				} // No predators left
				else
				{
					//if(Rho_avg[0] == 0)
					//{	gamma[s][i] = 0.0;	continue; } // No vegetation 
					//double rho_pred1 = 1/Rho_avg[2];
					for (int k = 0; k < eff_nR_Perp_size; k++)
					{
						gamma[s][i] += Rho_t[s - 1][nR_Perp[k]] * (1 - Rho_t[2][nR_Perp[k]] * rho_inverse[2]);
					}
					// RECALL : rho_inverse[s] = 1/Rho_avg[s]
				}
				gamma[s][i] = (gamma[s][i] * rho_inverse[s - 1]) / (eff_nR_Perp_size);

			}
			if (s == 2)
			{

				if (Rho_avg[0] < eps)
				{
					for (int k = 0; k < eff_nR_Perp_size; k++)
					{
						gamma[s][i] += Rho_t[s - 1][nR_Perp[k]];
					}
				} // No vegetation left
				else
				{
					//if(Rho_avg[1] < eps)
					//{	gamma[s][i] = 0.0;	continue; } // No grazers left.
					//double rho_veg1 = 1/Rho_avg[0];
					for (int k = 0; k < eff_nR_Perp_size; k++)
					{
						gamma[s][i] += Rho_t[s - 1][nR_Perp[k]] * (1 - Rho_t[0][nR_Perp[k]] * rho_inverse[0]);
					}
					// RECALL : rho_inverse[s] = 1/Rho_avg[s]
				}
				gamma[s][i] = (gamma[s][i] * rho_inverse[s - 1]) / (eff_nR_Perp_size);
			}

			if (gamma[s][i] > 1)
			{
				gamma[s][i] = 1.0;
			} //Ensuring gamma is not greater than 1.
			else if (gamma[s][i] < 0)
			{
				gamma[s][i] = 0.0;
			} //Ensuring gamma is not less than 0.
			else if (isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
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

// Analytical solution of FPE for movement of grazer and predator species with periodic boundary conditions.
void advdiff_FKE_MultiSp(D2Vec_Double& Rho_t, D2Vec_Double& Rho_tsar, D2Vec_Double& gamma, D2Vec_Pair_Double& v_eff, int g)
{
	int L2 = g*g;
	for (int s = 1; s < SpB; s++)
	{
		for (size_t i = 0; i < L2; i++)
		{
			if (Rho_t[s][i] < epsilon)
			{
				continue;
			} //No species left.

			int c_i = int(i / g); int c_j = i % g; //Current x and y coordinates of site i.

			double gamma2_si = gamma[s][i] * gamma[s][i]; // gamma^2, used to reduce computation.

			// General case where both advection and diffusion terms are present. O(L^2*5) ---> o(L^2*l^2)
			// Define bounds for the stencil
			int min_p = c_i - sig_D_ScaledBounds[s]; int max_p = c_i + sig_D_ScaledBounds[s];
			int min_q = c_j - sig_D_ScaledBounds[s]; int max_q = c_j + sig_D_ScaledBounds[s];

			//Next define the mean coordinates for the advection term.
			int mu_x = static_cast<int>(round(c_i + v_eff[s][i].first * (mu_vel_prefactor[s])) + g) % g; // Assuming initial condition v_x(t) = v_y(t) =0
			//RECALL  mu_vel_prefactor[s]*v_eff[s][i].first=  v_eff[s][i].first*(1/dx)*(Dt_ay - ((1 - exp_beta_Dt[s])/beta_vel[s]); 
			int mu_y = static_cast<int>(round(c_j + v_eff[s][i].second * (mu_vel_prefactor[s])) + g) % g; // Assuming initial condition v_x(t) = v_y(t) =0
			//int mu_y = static_cast<int>(round(c_j + (1/dx)*v_eff[s][i].second*(Dt_ay - ((1 - exp_beta_Dt[s])/beta_vel[s]) )) + g)  % g;

			int min_p_vx = mu_x - sig_vD_ScaledBounds[s].first; int max_p_vx = mu_x + sig_vD_ScaledBounds[s].second;
			int min_q_vy = mu_y - sig_vD_ScaledBounds[s].first; int max_q_vy = mu_y + sig_vD_ScaledBounds[s].second;


			// Update C_t[][] for the diffusion term using the Gaussian stencil
			for (int p = min_p; p <= max_p; ++p)
			{
				int px_mod = ((p + g) % g) * g;
				// Apply periodic boundary conditions
				int dx = std::abs(p - c_i); //int x = (p + g) % g; 
				for (int q = min_q; q <= max_q; ++q)
				{
					//int y = (q + g) % g;
					int dy = std::abs(q - c_j);
					int i_new = px_mod + (q + g) % g;
					Rho_tsar[s][i_new] += (gamma2_si)*Rho_t[s][i] * (gaussian_stencil_D[s][dx] * gaussian_stencil_D[s][dy]);
				}

				// Next update C_t[][] using joint probability of advection and diffusion terms
				for (int q = min_q_vy; q <= max_q_vy; ++q)
				{
					//int y = (q + g) % g;
					int dy = std::abs(q - mu_y);
					int i_new = px_mod + (q + g) % g;
					Rho_tsar[s][i_new] += (gamma[s][i] - gamma2_si) * Rho_t[s][i] * (gaussian_stencil_D[s][dx] * gaussian_stencil_VXY[s][dy]);
				}
			}
			// Update C_t[][] for the advection term using the Gaussian stencil (plus other mixture term)
			for (int p = min_p_vx; p <= max_p_vx; ++p)
			{
				int px_mod = ((p + g) % g) * g;
				// Apply periodic boundary conditions
				int dx = std::abs(p - mu_x); //int x = (p + g) % g; 
				for (int q = min_q; q <= max_q; ++q)
				{
					//int y = (q + g) % g;
					int dy = std::abs(q - c_j);
					int i_new = px_mod + (q + g) % g;
					Rho_tsar[s][i_new] += (gamma[s][i] - gamma2_si) * Rho_t[s][i] * (gaussian_stencil_VXY[s][dx] * gaussian_stencil_D[s][dy]);
				}

				for (int q = min_q_vy; q <= max_q_vy; ++q)
				{
					//int y = (q + g) % g;
					int dy = std::abs(q - mu_y);
					int i_new = px_mod + (q + g) % g;
					Rho_tsar[s][i_new] += (1 - gamma[s][i]) * (1 - gamma[s][i]) * Rho_t[s][i] * (gaussian_stencil_VXY[s][dx] * gaussian_stencil_VXY[s][dy]);
				}
			}

		}

		// Update Rho_t[][] with the new values from Rho_tsar[][]
		#pragma omp simd
		for (int i = 0; i < L2; i++)
		{
			Rho_t[s][i] = Rho_tsar[s][i];
		}

	}
}

void rietkerk_DorFPE_2D_MultiSp(D2Vec_Double& Rho, vector <double>& t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0,
	double(&D)[Sp], double(&v)[SpB], double(&K)[3], double sigma[], double a_st, double a_end, double a_c, double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double pR[],
	int(&dtV)[SpB], double clow[], double dt, double dx, double dP, int r, int g, double Gstar /* =-1.*/, double Vstar /* = -1.*/)
{

	#if defined(BARRACUDA) && SPB > 1
	// CUDA code, set up memory pointers and allocate memory on the GPU (parralel processing for advdiff_FKE_MultiSp() function)
	double* d_Rho_dt, * d_Rho_tsar, * d_gamma; std::pair<double, double>* d_v_eff;
	//double *d_sigD_ScaleBounds, *d_sigvD_ScaleBounds; // Pointers to the sigma bounds for the Gaussian stencils.
	//double *mu_vel_prefac; // Defined as constant arrays now.
	// RECALL: SpB = Sp - 2; Sp_NV = Sp - 3; // Number of consumers.
	// Allocate memory on the GPU for the various arrays.
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(&d_Rho_dt, Sp_NV * g * g * sizeof(double)); // Stores the curr population densities of consumers.
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaMalloc failed for d_Rho_dt in rietkerk_DorFPE_2D_MultiSp() function, with error: " 
		<< cudaGetErrorString(cudaStatus) << endl; //exit(1);
	}
	cudaStatus = cudaMalloc(&d_Rho_tsar, Sp_NV * g * g * sizeof(double)); // Stores the updated population densities of consumers.
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaMalloc failed for d_Rho_tsar in rietkerk_DorFPE_2D_MultiSp() function, with error: " 
		<< cudaGetErrorString(cudaStatus) << endl; //exit(1);
	}
	cudaStatus =cudaMalloc(&d_gamma, Sp_NV * g * g * sizeof(double)); // Stores the gamma values for the consumers.
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaMalloc failed for d_gamma in rietkerk_DorFPE_2D_MultiSp() function, with error: " 
		<< cudaGetErrorString(cudaStatus) << endl; //exit(1);
	}
	cudaStatus = cudaMalloc(&d_v_eff, Sp_NV * g * g * sizeof(std::pair<double, double>)); // Stores the effective velocities for the consumers.
	if (cudaStatus != cudaSuccess) {
		cerr << "cudaMalloc failed for d_v_eff in rietkerk_DorFPE_2D_MultiSp() function, with error: " 
		<< cudaGetErrorString(cudaStatus) << endl; //exit(1);
	}
	//cudaMalloc(&d_sigD_ScaleBounds, Sp_NV*sizeof(double)); // Stores the sigma bounds for the Gaussian stencils.
	//cudaMalloc(&d_sigvD_ScaleBounds, Sp_NV*sizeof(double)); // Stores the sigma bounds for the Gaussian stencils.
	//cudaMalloc(&mu_vel_prefac, Sp_NV*sizeof(double)); // Stores the mu_vel_prefactor[] array.

	// Check if memory allocation was successful.
	if (d_Rho_dt == NULL || d_Rho_tsar == NULL || d_gamma == NULL || d_v_eff == NULL)
	{
		cudaError_t err = cudaGetLastError();
		std::cerr << "Memory allocation failed on the GPU for the various arrays in rietkerk_DorFPE_2D_MultiSp() function. Exiting .... \n";
		std::cerr << "With CUDA Error: " << cudaGetErrorString(err) << std::endl;
		exit(1);
	}

	//Next create flattened vectors for Rho and Rho_tsar.
	//vector<double> Rho_dtNVflat(Sp_NV*g*g, 0.0); vector<double> Rho_dtNVflat(Sp_NV*g*g, 0.0); vector<double> gamma_NVflat(Sp_NV*g*g, 0.0); vector<double> v_eff_NVflat(Sp_NV*g*g, 0.0);

	#endif


	int thrID; //Unique ID for each thread, used for saving files.
	//RECALL: double epsilon = 1.0e-10; //Small number to avoid division by zero.
	double perc = 0.015; //Percentage of high density patches.
	//Store max value in pR[] in a variable r_max.
	//NOTE: SpB = Sp - 2; //Number of species excluding water terms.
	double r_max_effective = (*max_element(pR, pR + SpB)) / dx;
	// Store pR[]/dx as a fraction of r_max_effective in r_frac[].
	//Initialise vector<std::pair<int, int>> r_frac with initial size SpB and set all values to 0.0.
	vector <std::pair<double, int>> r_frac(SpB, { 0.0, 0 });

	for (int s = 0; s < SpB; s++)
	{
		r_frac[s] = { pR[s] / (dx * r_max_effective), s };
	}

	// Sort the r_frac[] array in descending order and store the indices in r_sort[].
	sort(r_frac.begin(), r_frac.end(), [](std::pair<double, int>& a, std::pair<double, int>& b)
		{ return a.first > b.first; });
	std::vector<std::pair<int, int>> origin_Neighbourhood = computeNeighboringSitesCentral(r_max_effective);

	if (omp_get_thread_num() == 1)
	{
		for (int i = 0; i < r_frac.size(); i++)
		{
			stringstream m0;	//To make cout thread-safe as well as non-garbled
			m0 << "r_frac[" << i << "]:\t" << r_frac[i].first << "\t and r_sort[" << i << "]:\t" << r_frac[i].second << endl;
			cout << m0.str();
		}
	}

	// Next, pre-calculate the various Gaussian stencils for the diffusion and advection terms.
	//D2Vec_Double gaussian_stencil_D(SpB, vector<double>(1, 0.0)); //Stores the Gaussian stencil for the diffusion term.
	//D2Vec_Double gaussian_stencil_VXY(SpB, vector<double>(1, 0.0)); //Stores the Gaussian stencil for the advection term.
	//Initialise the Gaussian stencil for the diffusion and advection terms.
	size_t size_gauss_stencil_D = 0; size_t size_gauss_stencil_VXY = 0;
	//Stores the size of the Gaussian stencil for the diffusion and advection terms FOR CUDA MEMORY ALLOCATION.
	for (int s = 1; s < SpB; s++)
	{
		size_gauss_stencil_D += gaussian_stencil_D[s].size(); size_gauss_stencil_VXY += gaussian_stencil_VXY[s].size();
		//Next, calculate the Gaussian stencil for the diffusion term using std::vector<double> precompute_gaussian(double sigma, int Rbound)
		//gaussian_stencil_D[s] = precompute_gaussian(sigma_D[s] * dx1, sig_D_ScaledBounds[s]);
		// And for the advection term [RECALL: sig_vD_ScaledBounds[s].first = 2*sig_vD_ScaledBounds[s].second]
		// IMPORTANT: We assume the advection contribution is symmetric in x and y directions
		//gaussian_stencil_VXY[s] = precompute_gaussian(sigma_vD[s].first * dx1, sig_vD_ScaledBounds[s].first);	
	}

	#if defined(BARRACUDA) && SPB >1
	
	double* devloc_gaussian_stencil_D, * devloc_gaussian_stencil_VXY;
	// Allocate memory on the GPU for the Gaussian stencils if either size_gauss_stencil_D or size_gauss_stencil_VXY
	// is greater than MAX_GAUSS_D_STENCIL_SIZE or MAX_GAUSS_VXY_STENCIL_SIZE respectively.
	if (size_gauss_stencil_D > MAX_GAUSS_D_STENCIL_SIZE || size_gauss_stencil_VXY > MAX_GAUSS_VXY_STENCIL_SIZE)
	{
		use_constmem_gaussstencils = false; // Switch to local device memory for Gaussian stencils.
		stringstream m0;
		m0 << "WARNING!!! Gaussian stencil size exceeds maximum allowed sizes (" << MAX_GAUSS_D_STENCIL_SIZE << " , " 
		<<  MAX_GAUSS_VXY_STENCIL_SIZE <<  ")for CUDA memory allocation with found sizes: [" << size_gauss_stencil_D << " , " 
		<< size_gauss_stencil_VXY << "]  in rietkerk_DorFPE_2D_MultiSp() function. Switching to local device memory! .... \n";

		cout << m0.str(); cerr << m0.str();
		cudaError_t allocerr, cudaStatus;
		allocerr = cudaMalloc(&devloc_gaussian_stencil_D, size_gauss_stencil_D * sizeof(double)); // Stores the Gaussian stencil for the diffusion term.
		allocerr = cudaMalloc(&devloc_gaussian_stencil_VXY, size_gauss_stencil_VXY * sizeof(double)); // Stores the Gaussian stencil for the advection term.

		size_t stencilD_offset = 0; size_t stencilVXY_offset = 0;
		for (int s = 1; s < SpB; s++)
		{
			// Copy the Gaussian stencils for the diffusion and advection terms to the GPU device.
			cudaStatus = cudaMemcpy(devloc_gaussian_stencil_D + stencilD_offset, gaussian_stencil_D[s].data(), gaussian_stencil_D[s].size() * sizeof(double), cudaMemcpyHostToDevice);
			stencilD_offset += gaussian_stencil_D[s].size();
			if (cudaStatus != cudaSuccess) {
				cerr << "cudaMemcpy failed for devloc_gaussian_stencil_D in rietkerk_DorFPE_2D_MultiSp() function, with error: "
					<< cudaGetErrorString(cudaStatus) << endl; //exit(1);
			}
			cudaStatus = cudaMemcpy(devloc_gaussian_stencil_VXY + stencilVXY_offset, gaussian_stencil_VXY[s].data(), gaussian_stencil_VXY[s].size() * sizeof(double), cudaMemcpyHostToDevice);
			stencilVXY_offset += gaussian_stencil_VXY[s].size();
			if (cudaStatus != cudaSuccess) {
				cerr << "cudaMemcpy failed for devloc_gaussian_stencil_VXY in rietkerk_DorFPE_2D_MultiSp() function, with error: "
					<< cudaGetErrorString(cudaStatus) << endl; //exit(1);
			}
		}
		//Finally report Gaussian stencil values on device using reportGaussianStencilValues() function.
		if (omp_get_thread_num() == 1)
		{
			stringstream m0;	//To make cout thread-safe as well as non-garbled
			for (int s = 1; s < SpB; s++) // Noting values Transferred to device!
				m0 << "gaussian_stencil_D[" << s << "]:\t { " << gaussian_stencil_D[s].data() << " }\t and gaussian_stencil_VXY[" << s << "]:\t {" << gaussian_stencil_VXY[s].size() << " }\n";
			cout << m0.str(); cerr << m0.str();
			reportGaussianStencilValues<double>(devloc_gaussian_stencil_D, devloc_gaussian_stencil_VXY);
		}
		// Check if memory allocation was successful.
		if (devloc_gaussian_stencil_D == NULL || devloc_gaussian_stencil_VXY == NULL)
		{
			std::cerr << "Memory allocation failed on the GPU for the Gaussian stencils in rietkerk_DorFPE_2D_MultiSp() function. Exiting .... \n";
			//
			cudaError_t err = cudaGetLastError();
			std::cerr << "With CUDA Error: " << cudaGetErrorString(err) << std::endl;
			exit(1);
		}
	}
	#endif

	if (omp_get_thread_num() == 1) // NOTING GAUSSIAN STENCIL SIZES!
	{
		for (int i = 0; i < SpB; i++)
		{
			stringstream m0;	//To make cout thread-safe as well as non-garbled
			m0 << "gaussian_stencil_D[" << i << "]:\t" << gaussian_stencil_D[i].size() << "\t and gaussian_stencil_VXY[" << i << "]:\t" << gaussian_stencil_VXY[i].size() << "\n";
			cout << m0.str(); cerr << m0.str();
		}
	}

	//Initialise variables for the RK4 integration.
	D3Vec_Int nR2(g * g, D2Vec_Int(2, vector<int>(2, 0)));
	//nR2 is 3D [(L*L)x2x2] vector initialised to 0, which stores neighbours of all sites i (second dim) in a ball of NORM 1 (SQUARE) radius 1 (2ND & 3RD dim).
	// nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 

	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_Sq4(g, nR2); //Assigning neigbours.
	//Append the value 0.0 to the start of t_meas.
	//t_meas.insert(t_meas.begin(), 0.0); //Inserting 0.0 at the start of t_meas.
	int tot_iter = t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	t_max = t_meas[t_meas.size() - 1]; //Only compute until end of sampling events

	int frame_tot_iter = frame_tmeas.size(); // Stores total number of iterations (i.e. time steps logged) for frame saving.

	//double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0}; 
	D2Vec_Double rho_rep_avg_var(tot_iter, vector<double>(Sp4_1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively. rho_rep_avg_var[0][0] = 0;
	//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively. rho_rep_avg_var[0][0] = 0;
	for (int i = 0; i < t_meas.size(); i++)
	{
		rho_rep_avg_var[i][0] = t_meas[i];
	}  //Setting up time index of rho_t_N_R2_rep


	//Defining Dornic variables for integration schema.
	//double bdt = b*dt; double cdt = c*dt; 
	//double dx2 = dx*dx; double dx1 = 1/dx; double dx1_2 = 1/(dx2); double dxb2 = dx/2.0; 
	double diff_coefficient[Sp] = { 0.0 }; double beta[Sp] = { 0.0 }; double lambda_exp[Sp] = { 0.0 }; double lambda[Sp] = { 0.0 }; double sigma2_1[Sp] = { 0.0 };
	vector<std::pair<double, double>> diff_eff(SpB, { 0.0, 0.0 }); //Stores effective diffusion coefficient for each species after correcting for advection.
	//vector <double> gamma_i(SpB, 0.0); //Stores gamma for each species.


	beta[0] = -M[0] - 4 * D[0] / (dx2);
	lambda_exp[0] = exp((beta[0]) * dt); lambda[0] = 2 * (beta[0]) / (sigma[0] * sigma[0] * (lambda_exp[0] - 1.0));
	diff_coefficient[0] = D[0] / (dx2); sigma2_1[0] = 1 / (sigma[0] * sigma[0]);
	for (int s = 1; s < SpB; s++)
	{
		diff_coefficient[s] = D[s] / (dx2);
		sigma2_1[s] = 1 / (sigma[s] * sigma[s]);
		//diff_coefficient[s] = (D[s] -(dx/2.0)*(v[s]*(1 - dt*vdx[s])))/(dx2);
		// Advection leads to excess diffusion, hence the correction term.
		if (diff_coefficient[s] <= 0)
		{
			stringstream m0;	//To make cout thread-safe as well as non-garbled due to race conditions.
			m0 << "You Fucked Up With Species:\t" << s << "\t with D_eff[s]:\t " << diff_coefficient[s]
				<< "\t dt,dx:\t " << dt << "," << dx << "\t v:\t " << v[s] << endl; //cout << m1.str();
			cout << m0.str();
		}
	}
	for (int s = SpB; s < Sp; s++)
		diff_coefficient[s] = D[s] / (dx2);
	//double Gstar = M[2]/((E[2] -M[2]*H[1][2])*A[1][2]); // MFT Coexistance Steady state for grazer.
	double alpha_prime = diff_coefficient[0] * 4; //Assume each Rho in neighbourhood for alpha estimation averages 1.
	double poiss_ru = lambda[0] * lambda_exp[0] * 2.5; //Mean
	double mu_nought = 2.0 * alpha_prime / (sigma[0] * sigma[0]) + poiss_ru;  // For Gamma, Beta =1. So Mean = 2*alpha/(sigma^2) + lambda = alpha/beta
	poiss_ru += 5 * sqrt(poiss_ru); //Mean of Poisson Sampler + 5 SD.
	mu_nought += 6 * sqrt(mu_nought); //Mean of  Mean of Gamma Sampler + 6 SD.

	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
	m0 << "Parameters:\t a:\t" << a << "\t c:\t " << c << "\t dt,dx:\t " << dt << "," << dx << "\t gmax:\t " << gmax
		<< "\n Lambda:\t" << lambda[0] << "\t Beta:\t " << beta[0] << "\t dt,dx:\t " << dt << "," << dx << "\t Poisson Cut-Off  For Veg:\t " << poiss_ru
		<< "\t and  Gamma Cut-Off For Veg:\t " << mu_nought << "\t and  Alpha  For Veg:\t " << alpha
		<< "\n with MFT biomass density of Grazer at Coexistance = " << Gstar << " kg/km^2\n" << endl; //cout << m1.str();
	omp_get_thread_num() == 1 ? cout << m0.str() : cout << "";
	omp_get_thread_num() == 1 ? cerr << m0.str() : cerr << "";



	//const vector<int> initcsv_columns = {2, 3, 4}; //Columns to be read from csv file.

	// Const value of dp/50000 for Rho indexed 1 (i.e. GRAZER). Other initial values determined by burn-in csv files through init_csvconstframe(...)

	auto start_t = high_resolution_clock::now(); //Start the clock
	for (int j = 0; j < r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
		m1 << "Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num()
			<< "\t with Total Iterations:\t " << tot_iter << endl; cout << m1.str();


		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd + j);
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number

		D2Vec_Double Rho_dt(Sp, vector<double>(g * g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double DRho(Sp, vector<double>(g * g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		D2Vec_Double Rho_tsar(Sp, vector<double>(g * g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.

		D2Vec_Double gamma(SpB, vector<double>(g * g, 0.0)); //Stores gamma for each species at each site.

		// 2D vector of pairs of dim: SpBx(L*l), called v. Stores the velocity of each species at each site.
		D2Vec_Pair_Double v_eff(SpB, std::vector<std::pair<double, double>>(g * g, { 0.0, 0.0 }));

		// 2D vectors of dim: Sx(L*l). Used for RK4 integration (of linear terms) of Rho_dt , after Dornic update.
		D2Vec_Double K1(Sp, vector<double>(g * g, 0.0)); D2Vec_Double K2(Sp, vector<double>(g * g, 0.0));
		D2Vec_Double K3(Sp, vector<double>(g * g, 0.0)); D2Vec_Double K4(Sp, vector<double>(g * g, 0.0));

		D2Vec_Double Rho_M(tot_iter, vector<double>(Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		D2Vec_Double Rho_Mov(tot_iter, vector<double>(SpB4, 0.0)); //Stores averaged gamma and velocity values per time step for a single replicate.

		double Rhox_avg[Sp] = { 0.0 }; double Rhox_num[Sp] = { 0.0 }; //Stores spatial average and number of occupied sites at given t.

		//init_randbistableframe(Rho_dt, g*g, a, a_c,  perc, chigh, clow); // Returns a frame with random speckles of high and low density.
		///** GENERIC EXPRTK INITIALISATIONS
		// MFT PERTURBATION BASED FRAME INITIALISATION
		#if defined(INIT) && INIT == 0
		// HOMOGENEOUS FRAME INITIALISATION
		init_exprtk_homogenousMFTframe(Rho_dt, g * g, a, a_c, clow); // Returns a frame with random speckles of high and low density.
		#elif defined(INIT) && INIT == 2
		// BURN-IN FRAME INITIALISATION (EXPRTK)

		init_burnin_wrapper(Rho_dt, a, a_c, dP, perc, g, j, clow);

		vector<double> temp_vec = { Rho_dt[0].begin(),Rho_dt[0].end() }; //Rho_dt for species 's'
		double init_veg_num = occupied_sites_of_vector(temp_vec, g * g); //Finds number of occupied at given t.
		vector<double>().swap(temp_vec); //Flush temp_vec out of memory.

		if (init_veg_num == 0.0)
		{
			stringstream m1_1;
			m1_1 << "RUN-TIME WARNING: ZERO Active BURN-IN veg sites for Thread Rank:\t " << omp_get_thread_num()
				<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << " CONSIDER RE-RUNNING OVER A DIFFERENT RANGE.\n";
			cout << m1_1.str(); cerr << m1_1.str();
		}
		#else
		// DEFAULT FRAME INITIALISATION (RANDOM MFT BASED SPECKLES)
		init_exprtk_randbiMFTframe(Rho_dt, g * g, a, a_c, dP, perc, clow); // Returns a frame with random speckles of high and low density.
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
		for (int s = 0; s < Sp; s++)
			amp[s] = clow[s + Sp];

		vector <double> sd{ g / 8.0, g / 16.0 }; //Setting standard deviation of gaussian distributions of vegetation to g/8 and g/16.

		init_gaussframe(Rho_dt, g * g, sd, amp); //**/

		// The first SpB species stored in Rho_dt are to be initialised on frame as Gaussian distributions.
		// Species 0 is the vegetation, Species 1 is the grazer and Species 2 is the predator.
		// Species 0 should be centered near the top right corner of the grid, Species 1 near the bottom right corner and Species 2 near the bottom left corner. */


		//*/

		for (int i = 0; i < g * g; i++)
		{
			for (int s = 0; s < Sp; s++)
			{
				DRho[s][i] = Rho_dt[s][i]; //Assigning initial conditions
				Rho_tsar[s][i] = Rho_dt[s][i];
			}
		}

		stringstream m1_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
		m1_2 << "Initial Conditions:\t Per:\t" << perc << "\t C_Lo[0], C_Lo[1], C_Lo[2]:\t "
			<< clow[0] << "," << clow[1] << "," << clow[2] << ",\t R: " << a << "\n"; //cout << m1.str();
		if (omp_get_thread_num() == 1)
			cout << m1_2.str();  //Print initial conditions for first thread.

		poisson_distribution<int> poisson; gamma_distribution <double> gamma_distr;
		uniform_real_distribution<double> unif(0.0, 1.0); normal_distribution<double> norm(0.0, 1.0);

		//Define thrID as a random number between 0 and 10000, iff j =0.
		thrID = (j == 0) ? int(round(unif(rng) * 10000)) : thrID; //Random number between 0 and 1000.

		double t = 0; int index = 0;  //Initialise t
		int frame_index = 0; int po = 1; int so = 1; int lo = 1; int counter = 0;
		vector <std::pair<int, int>> dtV_counter(SpB, { 0, 0 });
		//Keeps track of dtV[s] and counter%dtV[s] for each species s. Used to determine if species s is updated at time t.
		// Note counter := t/dt, so counter%dtV[s] = 0 implies t%dtV[s] = 0.

		for (int s = 0; s < SpB; s++) {
			if (dtV[s] == 1)
				dtV_counter[s] = { 1, 0 }; //If dtV[s] = 1, then update species s at every time step.
			else
			{
				dtV_counter[s] = { dtV[s], 0 };
			} //else update species s at every dtV[s] time step.
		}

		while (t < t_max + dt)
		{
			/** Update dtV_counter for each species s.
			for (int s=1; s< SpB; s++)
			{
				if(dtV[s] != 1)
					dtV_counter[s].second = (counter)%dtV_counter[s].first;
			}
			//*/

			// Basic Book-keeping below. Prints frames out to file at given time points. 
			//Also updates number and avg density of active sites at given time points.
			if (t == 0 || t >= t_meas[index] - dt / 2.0 && t < t_meas[index] + dt / 2.0)
			{
				double rhox_avg, rhox_num; //Stores spatial average and number of occupied sites at given t.
				for (int s = 0; s < Sp; s++)
				{
					vector <double> temp = { Rho_dt[s].begin(),Rho_dt[s].end() }; //Rho_dt for species 's'
					rhox_avg = mean_of_vector(temp, g * g); //Finds spatial average of densities at given t.
					rhox_num = occupied_sites_of_vector(temp, g * g); //Finds number of occupied at given t.
					Rho_M[index][2 * s] = rhox_num; Rho_M[index][2 * s + 1] = rhox_avg;

					vector<double>().swap(temp); //Flush temp out of memory.0

					//if(rhox_num == 0.0 and s== 0)
					//    break; //Finish negotiations as all basal species are dead.
				}

				//Next store spatial average of gamma and velocity for each species at given t.
				#if SPB > 1
				double mov_avg, mov_var;
				for (int s = 1; s < SpB; s++)
				{
					vector <double> temp = { gamma[s].begin(),gamma[s].end() }; //gamma for species 's'
					mov_avg = mean_of_vector(temp, g * g); //Finds spatial average of densities at given t.
					mov_var = variance_of_vector(temp, g * g); //Finds spatial average of densities at given t.
					//rhox_num = occupied_sites_of_vector(temp, g*g); //Finds number of occupied at given t.
					Rho_Mov[index][4 * s] = mov_avg; Rho_Mov[index][4 * s + 1] = mov_var; //Stores spatial average of gamma for species s at given t.
					vector<double>().swap(temp); //Flush temp out of memory.0

					// v_eff is a 2D vector of pairs of dim: SpBx(L*l), called v. Stores the velocity of each species at each site.
					// Compute the spatial average of the velocity in the x and y directions for each species.
					vector <double> temp_vx; //v_eff for species 's'
					vector <double> temp_vy; //v_eff for species 's'
					for (const auto& vel : v_eff[s])
					{
						temp_vx.push_back(vel.first); temp_vy.push_back(vel.second);
					}
					mov_avg = mean_of_vector(temp_vx, g * g); //Finds spatial average of densities at given t.
					Rho_Mov[index][4 * s + 2] = mov_avg; //Stores spatial average of gamma for species s at given t.

					mov_avg = mean_of_vector(temp_vy, g * g); //Finds spatial average of densities at given t.
					Rho_Mov[index][4 * s + 3] = mov_avg; //Stores spatial average of gamma for species s at given t.

					vector<double>().swap(temp_vx); vector<double>().swap(temp_vy); //Flush temp out of memory.
				}
				#endif

				vector <double> temp_alt = { DRho[0].begin(),DRho[0].end() }; //Rho_dt for species '0'
				rhox_num = occupied_sites_of_vector(temp_alt, g * g); //Finds number of occupied at given t.
				vector<double>().swap(temp_alt); //Flush temp out of memory.0

				// Saving frames to file at given time points.
				// BLOCK FOR CALCULATING AND TEMP PRELIMINARY FRAMES
				if (index == int(tot_iter * 0.85) || index == int(tot_iter * 0.9) || index == int(tot_iter * 0.95) || index == int(tot_iter - 1))
				{
					int break_flag = save_prelimfileswrapper(index, tot_iter, j, thrID, t, dt, dx, Rho_dt, DRho, Rho_M, Rho_Mov, rho_rep_avg_var,
						t_meas, gamma, v_eff, t_max, a, c, gmax, alpha, rW, W0, D, v, K, sigma, a_st, a_end, a_c, A, H, E, M,
						dP, r, g, Gstar /* =-1.*/, Vstar /* = -1.*/);
					if (break_flag == 1)
						break; //Finish negotiations as all basal species are dead.
				}
				index += 1;
				/** // STATUS UPDATE AFTER SAVING FRAME
				vector <double> temp_1= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
				double rhox_DR = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
				double  rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
				vector<double>().swap(temp_1); //Flush temp out of memory.0
				stringstream m4_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m4_1 << "AFTER SAVING FRAME, STATUS UPDATE AT TIME [t,i]\t" << t_meas[index] << " with # non-zero sites in VEG Rho_dt_0:  " << rhox_DR
				<< " AND # Non-Zero in VEG DRho_0:  " << rhox_dt <<  "\n"; cout << m4_1.str();
				 **/
				 //int counter_nzero_alpha=0;
			}
			// BASIC BOOK-KEEPING. Saves frame to file at given time points.
			if (t == 0 || t >= frame_tmeas[frame_index] - dt / 2.0 && t < frame_tmeas[frame_index] + dt / 2.0)
			{
				// Saving frames to file at given time points.
				if (frame_index >= frame_tot_iter - 10 && frame_index <= frame_tot_iter - 1 || t >= 50000 && t <= 150000 || t >= 2 && t <= 15000
					|| t == 0)
				{
					save_framefileswrapper(frame_index, frame_tot_iter, j, thrID, t, dt, dx, Rho_dt, DRho, rho_rep_avg_var,
					frame_tmeas, gamma, v_eff, t_max, a, c, gmax, alpha, rW, W0, D, v, K, sigma, a_st, a_end, a_c, A, H, E, M,
					dP, r, g, Gstar /* =-1.*/, Vstar /* = -1.*/);	
				}
				frame_index += 1;
			}

			// INTEGRATE ALL REACTION TERMS FIRST (DORNIC FOR VEG LINEAR + STOCHASTIC, RK4 FOR OTHERS)
			
			// First, Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC. ONLY LIVING MATTER EXPERIENCES DORNIC INTEGRATION.
			//First up vegetation.

			for (int s = 0; s < 1; s++)
			{
				for (int i = 0; i < Rho_dt[0].size();i++)
				{
					double alpha_i = diff_coefficient[s] * ((DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
						DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]]));

					if (alpha_i == 0 && DRho[s][i] == 0)
					{
						continue;
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					//counter_nzero_alpha +=1;

					double mu = -1.0 + 2.0 * alpha_i * sigma2_1[s]; // NOTE: sigma2_1 = 1/(sigma^2)
					double ziggy = lambda[s] * lambda_exp[s] * Rho_dt[s][i];
					double gru; // Stores the Gamma random variable.
					if (ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if (ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng) * sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if (gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng) * sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i] = gauss / lambda[s];
					}
					else
					{
						gamma_distr = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i] = gamma_distr(rng) / lambda[s];
					}

					if (isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "VEG HIT ROCK BOTTOM WITH: Rho*" << s << "[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num()
							<< " at time:\t:" << t << " with Rho" << s << "(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i <<
							"\t and GAMMA ENTRY:\t" << gru << "\n"; cout << m6.str(); cerr << m6.str();
					}
					//Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
				}
			} // End of Vegetation Integration

			#if defined(DEBUG) ///** STATUS UPDATE FOR VEG DORNIC INTEGRATION
			if (counter % 200 == 0)
			{
				vector <double> temp_1 = { Rho_dt[0].begin(),Rho_dt[0].end() }; //Rho_dt for species '0'
				double rhox_dt_veg = occupied_sites_of_vector(temp_1, L2); //Finds number of occupied at given t.q
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr]\t" << t << " , " << a << " , " << omp_get_thread_num() << "\t VEG DORNIC DONE with non-zero elements in DRho_0:  "
					<< rhox_num_veg << "\t and non-zero elements in Rho_dt_0:  " << rhox_dt_veg << "\t and INDEX VAL, T_MEAS:\t" << index << " , " << t_meas[index] << "\n"; cout << m5_1.str();
				cerr << m5_1.str();
			}
			//*/
			#endif

			// Next, do the RK4 integration of all species without movement terms of grazer and predator.

			#if SPB == 3
				RK4_Integrate_Stochastic_MultiSp(Rho_dt, Rho_tsar, K1, K2, K3, K4, nR2, a, c, gmax, alpha, rW, W0, diff_coefficient, K, A, H, E, M, t, dt, dx, g);
			#elif SPB == 2
				RK4_Integrate_Stochastic_2Sp(Rho_dt, Rho_tsar, K1, K2, K3, K4, nR2, a, c, gmax, alpha, rW, W0, diff_coefficient, K, A, H, E, M, t, dt, dx, g);
			#elif SPB == 1
				RK4_Integrate_Stochastic_1Sp(Rho_dt, Rho_tsar, K1, K2, K3, K4, nR2, a, c, gmax, alpha, rW, W0, diff_coefficient, K, t, dt, dx, g);
			#endif

			if (counter % 50000 == 0) // STATUS UPDATE FOR RK4 INTEGRATION OF REMAINDER TERMS
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j << "\t RK4 INTEGRATION OF REMAINDER TERMS DONE."
					<< "\n";
				omp_get_thread_num() == 1 ? cout << m5_1.str() : cout << "";
				cerr << m5_1.str();
			}

			#if defined(DEBUG)
			if (counter % 200 == 0) // STATUS UPDATE FOR GRAZER + PRED MOVEMENT INTEGRATION
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j << "\t GRAZER + PRED MOVEMENT INTEGRATION DONE, FOR " << prefix + "\n";
				cout << m5_1.str(); cerr << m5_1.str();
			}
			#endif
			// ALL REACTION TERMS INTEGRATED, NOW WORK ON MOVEMENT TERMS OF GRAZER + PREDATOR.
			
			// NEXT, calculate gamma and v_eff for each species at each site.
			// Compute ONLY IF SPB > 1
			#if SPB > 1
			for (int s = 0; s < SpB; s++)
			{
				vector <double> temp = { Rho_dt[s].begin(),Rho_dt[s].end() }; //Rho_dt for species '0'
				Rhox_avg[s] = mean_of_vector(temp, L2); //Finds spatial average of densities at given t.
				vector<double>().swap(temp); //Flush temp out of memory.0
			}
			vector <double> temp_veg = { Rho_dt[0].begin(),Rho_dt[0].end() }; //Rho_dt for species '0'
			double rhox_num_veg = occupied_sites_of_vector(temp_veg, L2); //Finds number of occupied sites at given t.
			vector<double>().swap(temp_veg); //Flush temp out of memory.
			double nR_fac = 1 - rhox_num_veg / (L2); //Factor to reduce the number of neighbours for gamma_i estimation
			if (nR_fac < 0.35)
			{
				nR_fac = 0.35;
			}
			#endif

			#if SPB == 3
				calc_gamma_vel_3Sp_NonRefugia(origin_Neighbourhood, rng, Rho_dt, gamma, v_eff, v, Rhox_avg, r_frac, dtV_counter, nR_fac, r_max_effective, g);
			#elif SPB == 2
				calc_gamma_vel_2Sp_NonRefugia(origin_Neighbourhood, rng, Rho_dt, gamma, v_eff, v, Rhox_avg, r_frac, nR_fac, r_max_effective, g);
			#endif 
			//Calculates gamma for each species at each site.

			// Finally, update the movement terms of grazer and predator using advdiff_FKE_MultiSp()
			
			#if SPB > 1
			if (counter % 40000 == 120)
			{
				// Report average densities of all species at given time.
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "DEBUG UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j << "\t AVERAGE DENSITIES BEFORE ADV/DIFF ARE: \n";
				for (int s = 0; s < SpB; s++)
				{
					vector <double> temp = { Rho_dt[s].begin(),Rho_dt[s].end() }; //Rho_dt for species 's'
					Rhox_avg[s] = mean_of_vector(temp, g * g); //Finds spatial average of densities at given t.
					double sd = standard_deviation_of_vector(temp, g * g); //Finds spatial average of densities at given t.
					vector<double>().swap(temp); //Flush temp out of memory.
					m5_1 << "Species " << s << ":\t" << Rhox_avg[s] << " +/- " << sd << "\n";
				}
				omp_get_thread_num() == 1 ? cout << m5_1.str() : cout << ""; //
				cerr << m5_1.str();

			}
			#if defined(BARRACUDA)
			// Temporary code for CUDA implementation.
			//Transfer Rho_dt, gamma and v_eff to device memory.
			int rho_offset = 0;
			for (int s = 1; s < SpB; s++)
			{
				cudaMemcpy(d_Rho_dt + rho_offset, Rho_dt[s].data(), g * g * sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(d_gamma + rho_offset, gamma[s].data(), g * g * sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(&d_v_eff[(s - 1) * g * g], v_eff[s].data(), g * g * sizeof(std::pair<double, double>), cudaMemcpyHostToDevice);
				rho_offset += g * g;
			}
			// Reset Rho_tsar on device memory to 0.0
			cudaMemset(d_Rho_tsar, 0.0, Sp_NV * g * g * sizeof(double));

			// CALL CUDA KERNEL
			launchAdvDiffKernel_MultiSp(d_Rho_dt, d_Rho_tsar, d_gamma, d_v_eff, g, devloc_gaussian_stencil_D, devloc_gaussian_stencil_VXY);

			// Swap Rho_tsar and Rho_dt on device memory.
			//double *temp = d_Rho_tsar; d_Rho_tsar = d_Rho_dt; d_Rho_dt = temp;
			// Transfer Rho_dt back to host memory (Rho_dt[s] for s > 0)
			rho_offset = 0;
			for (int s = 1; s < SpB; s++)
			{
				cudaMemcpy(Rho_dt[s].data(), d_Rho_tsar + rho_offset, g * g * sizeof(double), cudaMemcpyDeviceToHost);
				rho_offset += g * g;
			}
			#else
			// Standard CPU implementation.
			// First set Rho_tsar to 0.0
			for (int s = 1; s < SpB; s++)
			{
				for (int i = 0; i < Rho_dt[0].size(); i++)
					Rho_tsar[s][i] = 0.0;
			}
			advdiff_FKE_MultiSp(Rho_dt, Rho_tsar, gamma, v_eff, g);
			#endif

			if (counter % 40000 == 120)
			{
				// Report average densities of all species at given time.
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "DEBUG UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j << "\t AVERAGE DENSITIES AFTER ADV/DIFF ARE: \n";
				for (int s = 0; s < SpB; s++)
				{
					vector <double> temp = { Rho_dt[s].begin(),Rho_dt[s].end() }; //Rho_dt for species 's'
					Rhox_avg[s] = mean_of_vector(temp, g * g); //Finds spatial average of densities at given t.
					double sd = standard_deviation_of_vector(temp, g * g); //Finds spatial average of densities at given t.
					vector<double>().swap(temp); //Flush temp out of memory.
					m5_1 << "Species " << s << ":\t" << Rhox_avg[s] << " +/- " << sd << "\n";
				}
				omp_get_thread_num() == 1 ? cout << m5_1.str() : cout << "";
				cerr << m5_1.str();
			}
			#endif


			if (counter % 50000 == 0) // STATUS UPDATE FOR GRAZER + PRED MOVEMENT INTEGRATION
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j
					<< "\t GRAZER + PRED MOVEMENT INTEGRATION DONE, FOR " << prefix + "\n";
				cout << m5_1.str(); cerr << m5_1.str();
			}

			#if defined(DEBUG)
			if (counter % 200 == 0) // STATUS UPDATE FOR GRAZER + PRED MOVEMENT INTEGRATION
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j
					<< "\t GRAZER + PRED MOVEMENT INTEGRATION DONE, FOR " << prefix + "\n";
				cout << m5_1.str(); cerr << m5_1.str();
			}
			#endif


			for (int s = 0; s < Sp; s++)
			{


				for (int i = 0;i < Rho_dt[0].size();i++)
				{
					// CHECK FOR NAN AND INF
					if (Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true and so == 1)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO " << s << "  FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index [i, t, a, thr]:  " << i << " , " << t << " , "
							<< a << " , " << omp_get_thread_num() << "\t with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Movement Integration Term:\t"
							<< Rho_tsar[s][i];
						if (s < SpB)
							m6 << "\t and (VX, VY, GAMMA):  [" << v_eff[s][i].first << " , " << v_eff[s][i].second << " , " << gamma[s][i] << " ]" << "\n";
						else
							m6 << "\n";
						cout << m6.str(); cerr << m6.str(); lo = -1;
					}
					// */ 

					#if defined(DEBUG) // FLAG FOR CHECKING RHO_DT ERRORS
					if (counter % 20000 == 0 && counter == i)
					{
						stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_1 << "SUCCESFUL CHECKING FOR RHO_DT ERRORS AT TIME [t, a, s, i]\t" << t << " , " << a << " , " << s << " , " << i << "\n";
						cout << m5_1.str(); cerr << m5_1.str();
					}
					#endif

					if (isnan(Rho_dt[s][i] == true))
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO " << s << "  NAAN PARATHA:\t" << Rho_dt[s][i] << "\t at index [i, t, a, thr]:  " << i << " , " << t << " , "
							<< a << " , " << omp_get_thread_num() << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Movement Integration Term:\t"
							<< Rho_tsar[s][i];
						if (s < SpB)
							m6 << "\t and (VX, VY, GAMMA):  [" << v_eff[s][i].first << " , " << v_eff[s][i].second << " , " << gamma[s][i] << " ]" << "\n";
						else
							m6 << "\n";
						cout << m6.str(); cerr << m6.str(); lo = -1;
					}

					if (s < SpB)
					{
						if (gamma[s][i] < 0.0)
						{
							stringstream m6;      //To make cout thread-safe as well as non-garbled due to race conditions.
							m6 << "GAMMA_I FALLS BELOW 0,  AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j
								<< "\t WITH GAMMA_I:" << gamma[s][i] << " for species: " << s << "\n"; cout << m6.str(); cerr << m6.str();
						}
					}

					#if defined(DEBUG) // FLAG FOR CHECKING GAMMA < 0 ERRORS
					if (counter % 20000 == 0 && counter == i)
					{
						stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_1 << "SUCCESFUL CHECKING FOR GAMMA < 0 ERRORS AT TIME [t, a, s, i]\t" << t << " , " << a << " , " << s << " , " << i << "\n";
						cout << m5_1.str(); cerr << m5_1.str();
					}
					#endif

					#if defined(DEBUG) // Check if gamma_i is NaN or Inf, if so, SAVE GAMMA & exit the program.
					if (s < SpB)
					{
						if (isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
						{
							stringstream m5_3;     //To make cout thread-safe as well as non-garbled due to race conditions.
							m5_3 << "GAMMA_I BLOWS THE LID,  AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j
								<< "\t WITH GAMMA_I:" << gamma[s][i] << " for species: " << s << " with Rho_avg[0]: " << Rhox_avg[0] << " and Rho_avg[1]: " << Rhox_avg[1] << " and Rho_avg[2]: " << Rhox_avg[2]
								<< " \n and Rho_t[0][i]: " << Rho_dt[0][i] << " and Rho_t[1][i]: " << Rho_dt[1][i] << " and Rho_t[2][i]: " << Rho_dt[2][i] << "\n"; cout << m5_3.str();
							cerr << m5_3.str();

							stringstream L, tm, d3, p1, rini, a1, a2, Dm1, Dm2, dix, dimitri, geq, veq;

							L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
							rini << j; Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2];
							a1 << a_st; a2 << a_end; dimitri << dP; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;

							string ErrGamfilenamePattern = "/ERROR_GAMMA_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
								+ "_a_" + p1.str() + "_D1_" + Dm1.str() + "_D2_" + Dm2.str() + "_dx_" + dix.str() + "_R_";
							string parendir = "";

							// Creating a file instance called output to store output data as CSV.
							if (Gstar != -1)
								parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Geq_" + geq.str();
							else if (Vstar != -1)
								parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Veq_" + veq.str();
							else
								parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str();

							string gammaheader = "a_c,  x,  Gam0(x; t), Gam1(x; t), Gam2(x; t) \n"; //Header for output frame.
							save_frame(Rho_dt, parendir, ErrGamfilenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
							exit(4);

						}
					}
					if (counter % 20000 == 0 && counter == i)
					{
						stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_1 << "SUCCESFUL CHECKING FOR GAMMA INF ERRORS AT TIME [t, a, s, i]\t" << t << " , " << a << " , " << s << " , " << i << "\n";
						cout << m5_1.str(); cerr << m5_1.str();
					}
					// */
					#endif


					// Set all lattice sites that are less than 10^-12 to 0.0.
					if (Rho_dt[s][i] < epsilon && Rho_dt[s][i] >= 0.0)
						Rho_dt[s][i] = 0.0;
					DRho[s][i] = Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
					//Rho_tsar[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
				}
			}

			t += dt; //Update timestep
			counter += 1; //Update counter

			if (counter % 40000 == 1)
			{
				//Find memory usage and time taken for every 50000 iterations.
				//size_t currentSize = getCurrentRSS( ); //Check Memory usage.
				//size_t peakSize    = getPeakRSS( );
				auto end_t = high_resolution_clock::now();
				auto elapsed_min = duration_cast<minutes>(end_t - start_t);
				// Reset the clock
				start_t = high_resolution_clock::now();
				stringstream m7;
				m7 << "STATUS UPDATE AT [t, a, thr, j]\t " << t << " , " << a << " , " << omp_get_thread_num() << " , " << j
					//<< " ,with current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) 
					<< " and time taken for " << counter << " iterations: " << elapsed_min.count() << " min." << "\n"; cout << m7.str(); cerr << m7.str();
				//Report syst every 50000 iterations.
			}

			if (lo == -1)
			{
				stringstream m6_1;     //To make cout thread-safe as well as non-garbled
				m6_1 << "Saving error frames and exiting program .... \n"; cout << m6_1.str(); cerr << m6_1.str();

				// SAVING ERROR Rho_dt, DRho, and gamma vectors TO FILE
				stringstream L, tm, d3, p1, rini, gm, a1, a2, Dm0, Dm1, Dm2, alph, w0t, aij, hij, dix, dimitri, sig0, sig1, geq, veq;

				L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
				rini << j; Dm0 << D[0] * pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2];
				a1 << a_st; a2 << a_end; sig0 << sigma[0]; sig1 << sigma[1]; dimitri << dP; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;
				//gm << setprecision(3) << gmax; w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; 
				// Three replicates are over.

				string frame_parendir = "";

				// Creating a file instance called output to store output data as CSV.
				if (Gstar != -1)
					frame_parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Geq_" + geq.str();
				else if (Vstar != -1)
					frame_parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Veq_" + veq.str();
				else
					frame_parendir = frame_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str();

				// SAVING RHO_DT FRAME
				string filenamePattern_Rhodt = "/ERROR_RHO_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
					+ "_a_" + p1.str() + "_dx_" + dix.str() + "_R_";
				//Next, designate the naive filename
				string filename = filenamePattern_Rhodt + rini.str() + ".csv";
				save_frame(Rho_dt, frame_parendir, filenamePattern_Rhodt, a, a_st, a_end, t, dt, dx, dP, j, g);

				// SAVING DRHO FRAME
				string filenamePattern_DRho = "/ERROR_DRHO_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
					+ "_a_" + p1.str() + "_dx_" + dix.str() + "_R_";
				//Next, designate the naive filename
				filename = filenamePattern_DRho + rini.str() + ".csv";
				save_frame(DRho, frame_parendir, filenamePattern_DRho, a, a_st, a_end, t, dt, dx, dP, j, g);

				/** // SAVING GAMMA FRAME
				string filenamePattern_Gamma = "/ERROR_GAMMA_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
				+ "_a_" + p1.str()  + "_dx_"+ dix.str() + "_R_";
				//Next, designate the naive filename
				filename = filenamePattern_Gamma + rini.str() + ".csv";

				string gammaheader = "a_c,  x,  Gam0(x; t), Gam1(x; t) \n"; //Header for output frame.
				save_frame(gamma, frame_parendir, filenamePattern_Gamma, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
				*/

				// SAVING PRELIMINARY FRAME
				// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
				// and save to file.
				D2Vec_Double rho_rep_avg_var_temp(index, vector<double>(Sp4_1, 0.0));

				//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
				//std::copy(rho_rep_avg_var.begin(), rho_rep_avg_var.begin() + index, rho_rep_avg_var_temp.begin());
				var_mean_incremental_surv_runs(rho_rep_avg_var_temp, Rho_M, index, 0);

				string filenamePattern_Prelim = "/ERROR_SERIES_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
					"_D1_" + Dm1.str() + "_R_";
				//Next, designate the naive filename

				string prelim_parendir = "";
				if (Gstar != -1)
					prelim_parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Geq_" + geq.str() + "/TimeSeries";
				else if (Vstar != -1)
					prelim_parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Veq_" + veq.str() + "/TimeSeries";
				else
					prelim_parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "/TimeSeries";

				save_prelimframe(rho_rep_avg_var_temp, prelim_parendir, filenamePattern_Prelim, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);
				exit(3);
			}
		} // End of time loop.

		vector<vector<double>>().swap(gamma); vector<vector<double>>().swap(K1);
		vector<vector<double>>().swap(K2); vector<vector<double>>().swap(K3); vector<vector<double>>().swap(K4);
		vector<vector <double>>().swap(Rho_dt); vector<vector <double>>().swap(DRho); vector<vector <double>>().swap(Rho_tsar);
		// Freeing up memory.
		if (int(omp_get_thread_num()) == 5)
		{
			size_t currentSize = getCurrentRSS(); //Check Memory usage.
			size_t peakSize = getPeakRSS();
			stringstream m7;
			m7 << "In replicate " << j << ", current size in MB: " << currentSize / (1024.0 * 1024.0)
				<< " and Peak Size (in MB): " << peakSize / (1024.0 * 1024.0) << "\n"; cout << m7.str();
		}

		if (j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
			m7 << "INITIAL UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num()
				<< " THE SIZE OF RHO_M: " << tot_iter << "\n"; cout << m7.str();
			for (int i = 0; i < tot_iter; i++)
			{
				for (int s = 0; s < SpB; s++)
				{
					rho_rep_avg_var[i][4 * s + 3] = Rho_M[i][2 * s + 1]; // Avg density of frame at given time point i.
					rho_rep_avg_var[i][4 * s + 4] = 0.0; // Variance in avg frame density.
					rho_rep_avg_var[i][4 * s + 6] = Rho_M[i][2 * s]; // Number of active sites in frame at given time point i.

					if (Rho_M[i][2 * s] > 0) //Ensuring only surviving runs are considered
					{
						rho_rep_avg_var[i][4 * s + 5] = 1;
					} //One surviving run as of now
					else if (Rho_M[i][2 * s] == 0) //Ensuring only surviving runs are considered
					{
						rho_rep_avg_var[i][4 * s + 5] = 0;
					} //No surviving run as of now
				}
				rho_rep_avg_var[i][1] = Rho_M[i][2 * (Sp - 2) + 1]; // Avg density of frame at given time point i.
				rho_rep_avg_var[i][2] = Rho_M[i][2 * (Sp - 1) + 1]; // Avg density of frame at given time point i.

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

		if ((j + 1) % 2 == 1 || (j + 1) % 2 == 0)
		{	//Start logging at every multiple of 2
			stringstream L, tm, d3, p1, a1, a2, dimitri, rini, rini_prev, Dm, cgm, sig0, geq, veq, jID;

			a1 << a_st; a2 << a_end;
			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(5) << a; dimitri << dP;
			rini << j; Dm << setprecision(4) << D[2]; cgm << c * gmax; sig0 << sigma[0]; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;

			string parendir = "";
			if (Gstar != -1)
				parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Geq_" + geq.str();
			else if (Vstar != -1)
				parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str() + "_Veq_" + veq.str() + "/TimeSeries";
			else
				parendir = prelim_folder + a1.str() + "-" + a2.str() + "_dP_" + dimitri.str();

			double ran_jid = (unif(rng) * 1000.0) / 1000.0; jID << ran_jid; // Random number between 0 and 1.

			string filenamePattern = prelim_prefix + jID.str() + "_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str() +
				"_D2_" + Dm.str() + "_R_";

			save_prelimframe(rho_rep_avg_var, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, true);
		}

	} // End of r loop.

	// Free up all device memory.
	#if defined(BARRACUDA)
	cudaFree(d_Rho_dt); cudaFree(d_Rho_tsar); cudaFree(d_gamma);
	cudaFree(d_v_eff);
	// Free Gaussian stencils only if not using constant memory.
	if (!use_constmem_gaussstencils) {	cudaFree(devloc_gaussian_stencil_D); cudaFree(devloc_gaussian_stencil_VXY); }
	
	#endif

	size_t currentSize = getCurrentRSS(); //Check Memory usage.
	size_t peakSize = getPeakRSS();
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize / (1024.0 * 1024.0)
		<< " and Peak Size (in MB): " << peakSize / (1024.0 * 1024.0) << endl; cout << m9.str();

	for (int i = 0; i < tot_iter; i++)
	{	// Recall Rho is: | 	a		|    r 		|     L 		|    t 		|     <<W(t)>>x,r			|      <<O(t)>>x,r		|
		//..     <<Rho0(t)>>x,r			|    Var0[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |....

		Rho.push_back({ a, static_cast<double>(r), static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1],
					rho_rep_avg_var[i][2] });
		for (int s = 0;s < Sp - 2; s++)
		{
			Rho[i].insert(Rho[i].end(), { rho_rep_avg_var[i][4 * s + 3], rho_rep_avg_var[i][4 * s + 4],
			rho_rep_avg_var[i][4 * s + 5], rho_rep_avg_var[i][4 * s + 6] });
		}
	}
	stringstream m10;
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize / (1024.0 * 1024.0)
		<< " and Peak Size (in MB): " << peakSize / (1024.0 * 1024.0) << endl; cout << m10.str(); cerr << m10.str();
}


void first_order_critical_exp_delta_stochastic_MultiSp(int div, double t_max, double a_start, double a_end, double a_c, double c, double gmax, double alpha,
	double rW, double W0, double(&D)[Sp], double(&v)[SpB], double(&K)[3], double sigma[], double(&A)[SpB][SpB], double(&H)[SpB][SpB], double(&E)[SpB], double(&M)[SpB], double pR[], int(&dtV)[SpB], double clo[],
	double dt, double dx, double dP, int r, int g, double Gstar /*= -1.*/, double Vstar /*= -1.*/)
{
	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//double mean[Sp] = {1.0, 0.05}; double sd[Sp] = {0.25, 0.0125};
	//init_randframe(Rho_0, Sp,  g*g); //Returns Rho_0 with a full initial frame filled with 0.2.

	//set_density_dependentmortality(M[2]); //Sets the density dependent mortality rate.

	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "NOSTRA" << endl;
	// The pspace to iterate over.
	//vector <double> t_measure = logarithm10_time_bins(t_max, 1 * dt);
	//vector <double> t_measure = linspace(40, 1240, 21);
	vector <double> t_measure = linspace(50, 2000, 41);
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
	sort(t_measure.begin(), t_measure.end());

	t_measure.erase(unique(t_measure.begin(), t_measure.end()), t_measure.end());

	// Copy t_measure values to frame_tmeas [PRELIM AND FRAME TMEAS ARE THE SAME]
	frame_tmeas = t_measure;

	// Use switchsort_and_bait() to swap values in frame_tmeas b/w [50, 400] with linspace(50, 400, 8) 
	// and [400, 1200] with linspace(400, 1200, 9).
	//frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 50, 400, 8, "linspace", true);
	//frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 400, 1200, 9, "linspace", true);
	//Used for spreading analysis.
	frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 0, 40, 11, "linspace", true);
	frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 40, 240, 11, "linspace", true);
	frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 250, 1000, 16, "linspace", true);
	frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 1100, 2000, 10, "linspace", true);

	//if (g > 250)
	//{
	//	frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 0, 5000, 6, "linspace", true);
	//	frame_tmeas = switchsort_and_bait<double>(frame_tmeas, 6000, 10000, 3, "linspace", true);
	//}


	if (t_measure[t_measure.size() - 1] > 90000)
	{

		// Remove all elements of frame_tmeas that lie between 50000 and 100000.
		//frame_tmeas.erase(std::remove_if(frame_tmeas.begin(), frame_tmeas.end(), 
		//[](double x){return (x > 50000.0 && x < 100000.0);}), frame_tmeas.end());
		// Create linear window from 50000 to 100000, with 26 points (spaced 2500 apart).
		vector <double> t_frame_linearwindow = linspace(50000, 100000, 26);
		// Rounded down to 1 decimal place.
		// Insert the linear window to frame_tmeas vector after the element in frame_tmeas that is just less than the first element in t_frame_linearwindow.
		auto it = std::upper_bound(frame_tmeas.begin(), frame_tmeas.end(), t_frame_linearwindow[0]);
		frame_tmeas.insert(it, t_frame_linearwindow.begin(), t_frame_linearwindow.end());
		// Sort the frame_tmeas vector in ascending order and remove duplicates (if any).
		sort(frame_tmeas.begin(), frame_tmeas.end());
	}
	// Remove all duplicate elements in frame_tmeas.
	frame_tmeas.erase(unique(frame_tmeas.begin(), frame_tmeas.end()), frame_tmeas.end());

	cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" << endl;
	for (int i = 0; i < t_measure.size(); i++)
	{
		cout << t_measure[i] << " ";
	} 	cout << endl;

	cout << "Values in frame_tmeas (which is of total size " << frame_tmeas.size() << ") are :" << endl;
	for (int i = 0; i < frame_tmeas.size(); i++)
	{
		cout << frame_tmeas[i] << " ";
	} 	cout << endl;

	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd); // Seed the generator
	std::uniform_int_distribution<int> unif_int(0, 30); //Use this to generate random numbers between 0 and 30.


	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
	const int microToSeconds = 1000000;
	const double delay1 = unif_int(rng) * microToSeconds;     //5 seconds

	// Random delay to stagger the start of the parallel threads.
	usleep(static_cast<int>(delay1)); //Sleep for delay1 microseconds
	//std::this_thread::sleep_for(std::chrono::microseconds(static_cast<int>(delay1))); //Sleep for delay1 microseconds


	cout << "Delay 1 in progress... (" << delay1 / microToSeconds << "s)" << endl;
	std::vector<vector <double>> vec;
	// Stores collated output from parallel method calls in proper ascending order of p values.

	size_t currentSize = getCurrentRSS(); //Check Memory usage.
	size_t peakSize = getPeakRSS();

	cout << "On initialisation, current size in MB: " << currentSize / (1024.0 * 1024.0) << " and Peak Size (in MB): " << peakSize / (1024.0 * 1024.0) << endl;

	auto start = high_resolution_clock::now();

	int nProcessors = omp_get_max_threads();
	if (nProcessors > 32)
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
		for (int i = 0; i < a_space.size(); i++)
		{

			//init_randconstframe(Rh0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
			//type="Gam";
			stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
			message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() << endl;
			cout << message.str();
			std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
			/**
			 * Namely CExpRho_a is structured as:
			 * | 	a		|    t 		|     <<Rho1(t)>>x,r			|    Var[<Rho1(t)>x],r    |
			**/
			rietkerk_DorFPE_2D_MultiSp(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, rW, W0, D, v, K, sigma, a_start, a_end, a_c, A, H, E, M, pR, dtV, clo, dt, dx, dP, r, g, Gstar, Vstar);
			//RK4_Wrapper_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , a_start, a_end, dt, dx, dP, r, g);

			vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
			vector<vector <double>>().swap(CExpRho_a);
			// Remove dynamic vectors from memory and free up allocated space.
		}

		#pragma omp for schedule(static) ordered
		for (int i = 0; i < omp_get_num_threads(); i++)
		{
			#pragma omp ordered
			vec.insert(vec.end(), vec_private.begin(), vec_private.end());
			// Inserting critical exponent data for each grid size in order.
			stringstream message3;
			size_t currSize = getCurrentRSS(); //Check Memory usage.
			size_t peaksSize = getPeakRSS();

			message3 << "Is this happening?\n" << "In thread: " << i << ", current size in MB: " << currSize / (1024.0 * 1024.0)
				<< " and Peak Size (in MB): " << peaksSize / (1024.0 * 1024.0) << endl;
			cout << message3.str();

		}
		vector<vector <double>>().swap(vec_private);
		// Remove dynamic vectors from memory and free up allocated space.
	} //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;


	std::uniform_real_distribution<double> unif(0.0, 1.0);
	//Store a uniformly distributed random number between 0 and 1 (rounded to 3 decimal places).
	double id = round(unif(rng) * 1000.0) / 1000.0; //Random number between 0 and 1.
	//Round to 3 decimal places.
	stringstream L, coco, tm, d3, p1, p2, rini, Dm0, Dm1, aij, hij, cgm, alph, dix, dimitri, Sig0, geq, veq, ID;

	//double Gstar = M[2]/((E[2] -M[2]*H[1][2])*A[1][2]); // MFT estimate of Grazer density at coexistence.

	L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << a_start; p2 << a_end; rini << r;
	alph << alpha; cgm << c * gmax; coco << setprecision(4) << c;
	Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; Sig0 << sigma[0]; dix << setprecision(2) << dx;
	aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; geq << setprecision(5) << Gstar; veq << setprecision(5) << Vstar;
	ID << id;

	string savefile = stat_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_D1_" + Dm1.str() + "_a1_" + p1.str() + "_a2_" + p2.str() + "_dx_" + dix.str();
	if (Gstar != -1)
		savefile += "_Geq_" + geq.str() /* + "_ID_" + ID.str() */ + "_R_" + rini.str() + ".csv";
	else if (Vstar != -1)
		savefile += "_Veq_" + veq.str() /* + "_ID_" + ID.str() */ + "_R_" + rini.str() + ".csv";
	else
		savefile += /* + "_ID_" + ID.str() */ +"_R_" + rini.str() + ".csv";

	ofstream output_1stdp;
	// Creating a file instance called output to store output data as CSV.
	output_1stdp.open(savefile);

	cout << "Save file name: \n";
	cout << savefile << endl;

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

	string header = prelimheader;
	output_1stdp << header << "\n";
	cout << "The vector elements are: " << endl;
	cout << header << "\n";

	for (int i = 0; i < vec.size(); i++)
	{
		for (int j = 0; j < vec[i].size() - 1; j++)
			output_1stdp << setprecision(16) << vec[i][j] << ",";
		output_1stdp << setprecision(16) << vec[i][vec[i].size() - 1] << endl;
		if (i % (1000) == 1)
		{
			for (int j = 0; j < vec[i].size() - 1; j++)
				cout << setprecision(16) << vec[i][j] << ",";
			cout << setprecision(16) << vec[i][vec[i].size() - 1] << endl;
		}
	}
	output_1stdp.close();
}
