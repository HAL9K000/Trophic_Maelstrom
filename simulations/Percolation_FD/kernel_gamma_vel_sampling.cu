// CUDA kernel for calculating gamma values with optimized memory access
// Designed for L=128, r_max~60 (neighborhood size ~11,000)

#define DEFINE_CUDA_CONSTANTS // Define the CUDA constant arrays only once through this .cu file in the header file
#include "MultiSPDP.h"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
//#include <utility>
#include <cstdio>

// Constant memory for neighborhood offsets (faster than global memory)
// Max 64KB, so we can store ~8000 int2 pairs

//#ifdef DEFINE_CUDA_CONSTANTS
//__constant__ std::pair<double, int> d_r_frac[CuSpB]; // Species perception radius fractions
//#endif

cudaError_t copyToDeviceConstantMemory_GammaSweepTerms(const DoubleIntPair* dr_frac)
{
    cudaError_t err;

    err = cudaMemcpyToSymbol(d_r_frac, dr_frac, CuSpB * sizeof(DoubleIntPair), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpyToSymbol failed for d_sigD_ScaleBounds: %s\n", cudaGetErrorString(err)); return err;
    }
    else {
        printf("d_sigD_ScaleBounds copied to constant memory with values:\n");
        for (int i = 0; i < CuSpB; i++)
            printf("(%f, %d)\t", dr_frac[i].first, dr_frac[i].second);
    }
    return err;
}


// build A-HAT on device: A = Rho_prev  OR  Rho_prev * (1 - Rho_pred * rho_inv_pred) for s==1
__global__ void build_A_kernel_3Sp(const float* __restrict__ d_Rho_t, float* __restrict__ d_A,
                               int s, int L2, float rho_inv_pred, bool pred_absent)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= L2) return;

    // species layout: d_Rho_t[s*L2 + i]
    const float rho_prev = d_Rho_t[(s-1) * L2 + i];
    if (s == 1 && !pred_absent) {
        float rho_pred = d_Rho_t[2 * L2 + i];
        d_A[i] = rho_prev * (1.0f - rho_pred * rho_inv_pred);
    } else {
        d_A[i] = rho_prev;
    }
}

// complex multiply a_hat * k_hat -> prod_hat (Nc complex points)
__global__ void complex_pointwise_mul_kernel(cufftComplex* __restrict__ a_hat,
                                             const cufftComplex* __restrict__ k_hat,
                                             cufftComplex* __restrict__ prod_hat,
                                             int Nc)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= Nc) return;
    cufftComplex A = a_hat[idx];
    cufftComplex K = k_hat[idx];
    cufftComplex P;
    P.x = A.x * K.x - A.y * K.y;
    P.y = A.x * K.y + A.y * K.x;
    prod_hat[idx] = P;
}

// normalize inverse-FFT result and write gamma: conv_real is unnormalized c2r result.
// gamma_out = clamp( (conv_real / N / kernel_sum) * rho_inv_resource )
__global__ void normalize_and_write_gamma_kernel(const float* __restrict__ conv_real,
                                                 float* __restrict__ gamma_out,
                                                 int N, float kernel_sum, float rho_inv_resource)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    float val = conv_real[i] / (float) N;
    float gamma = (kernel_sum > 0.0f) ? (val / kernel_sum) * rho_inv_resource : 0.0f;
    if (!isfinite(gamma)) gamma = 0.0f;
    if (gamma > 1.0f) gamma = 1.0f;
    else if (gamma < 0.0f) gamma = 0.0f;
    gamma_out[i] = gamma;
}



// Device kernel: Each thread processes one spatial site. Naive O[L2*K2] kernel sweep algorithm implementation.
__global__ void calc_gamma_3Sp_NonRefugia_kernel(
    const double* __restrict__ d_Rho_dt,      // [SpB][L*L] flattened
    double* __restrict__ d_gamma,            // [SpB][L*L] flattened
    std::pair<double, double>* d_v_eff,      // [SpB][L*L][2] effective velocities
    int2* d_origin_Neighbourhood, // [nR_Perp_length] neighborhood offsets
    const double* __restrict__ d_rho_inverse, // [SpB]
    //std::pair<double, int>* d_r_frac,  // [SpB][2] species indices, radius fractions sorted by radius
    const int* __restrict__ d_eff_sizes,      // [SpB] effective neighborhood sizes
    int L,
    int L2,
    int max_neighbors,
    double eps, double inv_eps)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= L2) return;
    
    int c_i = i / L;
    int c_j = i % L;

    double rho_inv_veg = d_rho_inverse[0];
    double rho_inv_grazer = d_rho_inverse[1];
    double rho_inv_pred = d_rho_inverse[2];
    //double inv_eps = 1/(2.0*eps);
    
    // Process each species (RECALL d_gamma stores ONLY Consumer species so of size [SpB-1]*[L2])
    for (int sp_idx = 0; sp_idx < CuSpNV; sp_idx++) 
    {
        int s = d_r_frac[sp_idx].second; // Should be 1 or 2 (Grazer or Predator)
        int gamma_s = s-1; // Should be 0 or 1 
        double fr = d_r_frac[sp_idx].first;
        int eff_size = d_eff_sizes[sp_idx];
        
        if (fr == 0.0 || s == 0) 
        {   continue;   } // Indicates vegetation, gamma not needed.
        if( eff_size < 1)
        {
            eff_size = 1; // At least one site nearest neighbour (the site itself is considered).
        }
        if(d_rho_inverse[s] > inv_eps || d_rho_inverse[s-1] > inv_eps)    
        {
            //Indicates no species or no resources left respectively, assume consumers advect to extinction.
            d_gamma[gamma_s * L2 + i] = 0.0;
            continue;
        }
        // Check for zero densities
        //const double* rho_s = d_Rho_dt + s * L2;
        const double* rho_prev = d_Rho_dt + (s - 1) * L2;
        const double* rho_pred = d_Rho_dt + 2 * L2;
        double gamma_sum = 0.0;

        #pragma unroll 8
        for (int k = 0; k < eff_size; ++k) // Iterating over all neighbours of a site
        {
            int2 offset = d_origin_Neighbourhood[k];
            int nx = c_i + offset.x; nx += (nx < 0) * L; nx -= (nx >= L) * L; //Branchless wrapping.
            int ny = c_j + offset.y; ny += (ny < 0) * L; ny -= (ny >= L) * L;   //Branchless wrapping.
            //int ny = (c_j + offset.y + L) % L;
            const int k_idx = nx*L +ny;
            // Species 1 (Grazers)
            if (s == 1) 
            {
                double rho_veg = __ldg(&rho_prev[k_idx]);
                if(d_rho_inverse[2] > inv_eps) // Predators no longer exist, gamma sum is simpler.
                {   gamma_sum += rho_veg;    }
                else //Full calculation otherwise
                {
                    // Vectorized reduction with coalesced memory access
                    double rho_pred_val = __ldg(&rho_pred[k_idx]);
                    gamma_sum += rho_veg * (1.0 - rho_pred_val * rho_inv_pred);    
                }  
            }
            // Species 2 (Predators)
            else if (s == 2) 
            {   gamma_sum += __ldg(&rho_prev[k_idx]);   }
        }
        //double rho_inv_veg = d_rho_inverse[0];
        if( s == 1)
        {   d_gamma[gamma_s * L2 + i] = fmin(1.0, fmax(0.0, (gamma_sum * rho_inv_veg) / eff_size));   }
        if( s == 2)
        {   d_gamma[gamma_s * L2 + i] = fmin(1.0, fmax(0.0, (gamma_sum * rho_inv_grazer/ eff_size)));   }

        // d_gammma is a probability and must be capped between 0 and 1.
        //if(d_gamma[gamma_s * L2 + i] > 1)
        //{  d_gamma[gamma_s * L2 + i] =1; }
        //else if(d_gamma[gamma_s * L2 + i] < 0)
        //{  d_gamma[gamma_s * L2 + i] =0; }    
    }
}

//Calculates gamma and velocities for 3 Sp Rietkerk model (details in PDF, 
// ASSUMING PREDATORS AREN'T DETERRED BY HIGH LOCAL VEGETATION DENSITY)
void calc_gamma_vel_NonRefugia_CUDA(int2* d_origin_Neighbourhood, double* d_Rho_dt, double* d_gamma, std::pair<double, double>* d_v_eff, double *d_rho_inverse,
	int* d_eff_sizes, double(&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, vector <std::pair<int, int>>& dtV_counter, int nR_Perp_length, int r_max, int L)
{
	//if(nVeg_frac < 0.35)
	//	nVeg_frac = 0.35; //Minimum fraction of vegetation cover.
	int L2 = L*L;
	double eps = 1.0e-12; //Small number to avoid division by zero.
	double inv_eps = 1/(2.0*eps);

    cudaError_t cudaStatus;
	// Calculate and copy rho_inverse to device.
    vector<double> rho_inverse(SpB);
	//double rho_inverse[SpB];// ={1/eps, 1/eps, 1/eps}; //Inverse of average density.
	for(int s=0; s< SpB; s++)
    {   rho_inverse[s] = (Rho_avg[s] >= eps) ?  (1/Rho_avg[s]) : (1/eps); }

    #if defined(DEBUG_CUDA)
    printf("Launching calc_gamma_vel_NonRefugia_CUDA kernel...\n");
    printf("Reporting parameters L2, nR_Perp_length, r_max, L, eps, inv_eps: %d , %d , %d , %d, %f, %f \n", L2, nR_Perp_length, r_max, L, eps, inv_eps);
    printf("Rho_avg values:\n");
    for (int s = 0; s < SpB; s++) 
        printf("[%d, %f]\t", s, Rho_avg[s]);
    printf("Rho_inverse values:\n");
    for (int s = 0; s < SpB; s++)
        printf("[%d, %f]\t", s, rho_inverse[s]);
    printf("\n");
    #endif
    
	// Avoid division by zero.
	cudaStatus = cudaMemcpy(d_rho_inverse, rho_inverse.data(), SpB * sizeof(double), cudaMemcpyHostToDevice);

    #if defined(DEBUG_CUDA)
    if (cudaStatus != cudaSuccess) {
        printf("cudaMemcpy to d_rho_inverse failedwith error: %s\n", cudaGetErrorString(cudaStatus));
    }
    #endif

	// Next calculate and transfer effective nR_Perp sizes to d_eff_sizes.
	//int eff_nR_sizes[SpB];
    vector<int> eff_nR_sizes(SpB); 
	//for(const auto& frac: d_rfrac)
    for(int fr_idx=0; fr_idx < SpB; fr_idx++)
	{
		int s = rfrac[fr_idx].second;
		double fr = rfrac[fr_idx].first; //Fraction of perception radius to max radius.
		eff_nR_sizes[s]= int(fr*fr*nR_Perp_length);
	}
	cudaStatus = cudaMemcpy(d_eff_sizes, eff_nR_sizes.data(), SpB * sizeof(int), cudaMemcpyHostToDevice);

    #if defined(DEBUG_CUDA)
    if (cudaStatus != cudaSuccess) {
        printf("cudaMemcpy to d_eff_sizes failed, with error: %s\n", cudaGetErrorString(cudaStatus));
    }
    reportPerceptionArrays << <1, 1 >> > (d_origin_Neighbourhood, d_rho_inverse, d_eff_sizes, nR_Perp_length, r_max, L);
    #endif

	//for(int s=0; s< SpB; s++)
	//	double s = int(d_rfrac[s].first*d_rfrac[s].first*nR_Perp.size());

    

	// Launch kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (L2 + threadsPerBlock - 1) / threadsPerBlock;
    //int sharedMemSize = nR_Perp_length * sizeof(int);

	//cudaError_t err = cudaFuncSetAttribute(calc_gamma_3Sp_NonRefugia_kernel<<<blocksPerGrid, threadsPerBlock>>>,
    //    cudaFuncAttributeMaxDynamicSharedMemorySize, sharedMemSize);
	//if (err != cudaSuccess) 
	//{
	//	cerr << "CUDA error in setting cudaFuncSetAttribute: " << cudaGetErrorString(err) << std::endl;
	//	exit(EXIT_FAILURE);
	//}
    #if SPB == 3
	    calc_gamma_3Sp_NonRefugia_kernel<<<blocksPerGrid, threadsPerBlock>>>(
            d_Rho_dt, d_gamma, d_v_eff, d_origin_Neighbourhood, d_rho_inverse,
            d_eff_sizes, L, L2, nR_Perp_length, eps, inv_eps);
    #endif

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) 
    {   printf("CUDA kernel error: %s\n", cudaGetErrorString(err));    }
    
    //cudaDeviceSynchronize();
}



//Function that reports values in Movement Arrays.
__global__ void reportPerceptionArrays(int2* d_origin_Neighbourhood, double *d_rho_inverse, 
    int* d_eff_sizes, int nR_Perp_length, int r_max, int L)
{
    printf("Rho Inverse Values: \n");
    for (int s = 0; s < CuSpB; s++)
    {
        printf("(%d, %f)\t", s, d_rho_inverse[s]);
    }
    printf("\n");

    printf("Effective neighbourhood sizes (d_eff_sizes): \n");
    for (int s = 0; s < CuSpB; s++)
    {
        printf("(%d, %d)\t", s, d_eff_sizes[s]);  
    }
    printf("\n");

    printf("Contrast with max perception neighborhood size: %d \n", nR_Perp_length);

    printf("d_r_frac values: \n");
    for (int s = 0; s < CuSpB; s++)
    {
        printf("[%d, %f, %d]\t", s, d_r_frac[s].first, d_r_frac[s].second);    
    }
    printf("\n");
    int k_step = (int)(nR_Perp_length / 20);
    printf("Some central neighboring sites offsets: \n");
    for (int k = 0; k < nR_Perp_length; k += k_step)
    {
        int2 offset = d_origin_Neighbourhood[k];
        printf("k=%d: (%d, %d)\n", k, offset.x, offset.y);
    }
    printf("\n");
}

void reportPerceptionArrays_Interface(int2* d_origin_Neighbourhood, double *d_rho_inverse, 
    int* d_eff_sizes, int nR_Perp_length, int r_max, int L)
{
    reportPerceptionArrays << <1, 1 >> > (d_origin_Neighbourhood, d_rho_inverse, d_eff_sizes, nR_Perp_length, r_max, L);
}