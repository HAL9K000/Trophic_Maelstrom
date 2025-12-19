#define DEFINE_CUDA_CONSTANTS // Define the CUDA constant arrays only once through this .cu file in the header file

#include "rietkerk_bjork_basic.h"
#include <cstdio>
#include <vector>



// AtomicAdd function analogue for double precision floating point numbers for devices with compute capability < 6.0
#if __CUDA_ARCH__ < 600
__device__ double atomicAdd_Double(double* address, double val)
{
    unsigned long long int* address_as_ull =
        (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(val +
                __longlong_as_double(assumed)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

cudaError_t copyToDeviceConstantMemory_AdvTerms(const int* sig_D_ScaledBounds, const int2* sigma_vD_ScBounds,
    const double* mu_vel_prefactor, const int* size_gauss_D, const int* size_gauss_VXY)
{
    cudaError_t err;

    err = cudaMemcpyToSymbol(d_sigD_ScaleBounds, sig_D_ScaledBounds, CuSpNV * sizeof(int), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpyToSymbol failed for d_sigD_ScaleBounds: %s\n", cudaGetErrorString(err)); return err;
    }
    else {
        printf("d_sigD_ScaleBounds copied to constant memory with values:\n");
        for (int i = 0; i < CuSpNV; i++)
            printf("%d\t", sig_D_ScaledBounds[i]);
    }

    err = cudaMemcpyToSymbol(d_sigvD_ScaleBounds, sigma_vD_ScBounds, CuSpNV * sizeof(int2), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpyToSymbol failed for d_sigvD_ScaleBounds: %s\n", cudaGetErrorString(err)); return err;
    }
    else {
        printf("\nd_sigvD_ScaleBounds copied to constant memory with values:\n");
        for (int i = 0; i < CuSpNV; i++)
            printf("[ %d , %d , %d ]\t", i, sigma_vD_ScBounds[i].x, sigma_vD_ScBounds[i].y);

    }

    err = cudaMemcpyToSymbol(d_mu_vel_prefactor, mu_vel_prefactor, CuSpNV * sizeof(double), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpyToSymbol failed for d_mu_vel_prefactor: %s\n", cudaGetErrorString(err));
        return err;
    }
    else {
        printf("\nd_mu_vel_prefactor copied to constant memory with values:\n");
        for (int i = 0; i < CuSpNV; i++)
            printf("%f\t", mu_vel_prefactor[i]);
    }

    err = cudaMemcpyToSymbol(d_size_gauss_D, size_gauss_D, CuSpNV * sizeof(int), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpyToSymbol failed for d_size_gauss_D: %s\n", cudaGetErrorString(err));
        return err;
    }
    else {
        printf("\nd_size_gauss_D copied to constant memory with values:\n");
        for (int i = 0; i < CuSpNV; i++)
            printf("%d\t", size_gauss_D[i]);
    }

    err = cudaMemcpyToSymbol(d_size_gauss_VXY, size_gauss_VXY, CuSpNV * sizeof(int), 0, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpyToSymbol failed for d_size_gauss_VXY: %s\n", cudaGetErrorString(err)); return err;
    }
    else {
        printf("\nd_size_gauss_VXY copied to constant memory with values:\n");
        for (int i = 0; i < CuSpNV; i++)
            printf("%d\t", size_gauss_VXY[i]);
    }

    return cudaSuccess;
}


cudaError_t copyToDeviceConstantMemory_GaussianStencils(size_t total_size_gauss_D, size_t total_size_gauss_VXY)
{
    cudaError_t err;

    size_t stencilD_offset = 0; size_t stencilVXY_offset = 0;
	for (int s = CuSpB - CuSpNV; s < SpB; s++)
	{
        printf("Copying Gaussian stencils for species %d to constant memory:\n", s);
        std::vector<float> temp_D(gaussian_stencil_D[s].begin(), gaussian_stencil_D[s].end());
        //printf("Temp_D:\t"); for (float val : temp_D) { printf("%f\t", val);}
		// Copy the Gaussian stencils for the diffusion and advection terms to the GPU device.
		err = cudaMemcpyToSymbol(d_gaussian_stencil_D, temp_D.data(), temp_D.size() * sizeof(float), stencilD_offset*sizeof(float), cudaMemcpyHostToDevice);
		stencilD_offset += gaussian_stencil_D[s].size();
        std::vector<float> temp_VXY(gaussian_stencil_VXY[s].begin(), gaussian_stencil_VXY[s].end());
        //printf("\nTemp_VXY:\t"); for (float val : temp_VXY) { printf("%f\t", val);}
		err = cudaMemcpyToSymbol(d_gaussian_stencil_VXY, temp_VXY.data(), temp_VXY.size()* sizeof(float), stencilVXY_offset*sizeof(float), cudaMemcpyHostToDevice);
		stencilVXY_offset += gaussian_stencil_VXY[s].size();
	}
    if (err != cudaSuccess) {
        printf("cudaMemcpyToSymbol failed for Gaussian Stencils: %s\n", cudaGetErrorString(err)); return err;
    }
    else 
    {
        printf("d_gaussian_stencil_D copied to constant memory with values:\n");
        for (int i = 0; i < CuSpNV; i++)
        {
            printf("For species %d,  sigD_ScaleBounds = %d\n", i, sig_D_ScaledBounds[i+(CuSpB - CuSpNV)]);
            for (int j = 0; j <= sig_D_ScaledBounds[i+(CuSpB - CuSpNV)]; j++) {
                printf("%f\t", gaussian_stencil_D[i+(CuSpB - CuSpNV)][j]);
            }
            printf("\n");
        }

        printf("d_gaussian_stencil_VXY copied to constant memory with values:\n");
        for (int i = 0; i < CuSpNV; i++)
        {
            printf("For species %d, d_sigvD_ScaleBounds = ( %d , %d )\n", i,
                sig_vD_ScaledBounds[i+(CuSpB - CuSpNV)].first, sig_vD_ScaledBounds[i+(CuSpB - CuSpNV)].second);
            for (int j = 0; j <= sig_vD_ScaledBounds[i+(CuSpB - CuSpNV)].first; j++) {
                printf("%f\t", gaussian_stencil_VXY[i+(CuSpB - CuSpNV)][j]);
            }
            printf("\n");
        }
    }
    //Assert if sums of sigma bounds match total sizes, if not return error.
    if (stencilD_offset != total_size_gauss_D) {
        printf("Error: Sum of sigD_ScaleBounds %zu does not match total_size_gauss_D: %zu\n", stencilD_offset, total_size_gauss_D);
        return cudaErrorInvalidValue;
    }
    if (stencilVXY_offset != total_size_gauss_VXY) {
        printf("Error: Sum of sigvD_ScaleBounds %zu does not match total_size_gauss_VXY: %zu \n", stencilVXY_offset, total_size_gauss_VXY);
        return cudaErrorInvalidValue;
    }
    
    return cudaSuccess;
}

// CUDA kernel to update the concentration field
template<bool useConstantGaussMemory> // Template parameter to choose between constant memory or passed pointers for Gaussian stencils.
__global__ void updateMovement(double* Rho_t, double* Rho_tsar, double* gamma, std::pair<double, double>* v_eff, int L,
    double* gaussian_stencil_D /**  = nullptr */, double* gaussian_stencil_VXY /**  = nullptr */)
{
    int L2 = L * L;
    // Thread indices
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;

    for (int s = 0; s < CuSpNV; s++)
    {
        // IMPORTANT: s = 0 is GRAZER, s = 1 is PREDATOR
        // s has been shifted by -1 compared to CPU code to save memory.

        // Iterate over lattice sites in a grid-stride loop
        for (int i = idx; i < L * L; i += stride)
        {
            int c_i = int(i / L);  // Row index
            int c_j = i % L;  // Column index

            

            if (Rho_t[s * L2 + i] < epsilon)
                continue;

            int sL2 = s * L2; int sL2_i = sL2 + i; // Precomputing species indices to save time.
            double gamma_si = gamma[sL2_i]; double gamma_2_si = gamma_si * gamma_si;
            double gammap_2_si = (1 - gamma_si) * (1 - gamma_si);

            // Diffusion stencil bounds
            int min_p = c_i - d_sigD_ScaleBounds[s]; int max_p = c_i + d_sigD_ScaleBounds[s];
            int min_q = c_j - d_sigD_ScaleBounds[s]; int max_q = c_j + d_sigD_ScaleBounds[s];

            // Advection mean coordinates
            int mu_x = (int)(round(c_i + v_eff[s * L2 + i].first * (d_mu_vel_prefactor[s])) + L) % L; // Assuming initial condition v_x(t) = v_y(t) =0

            int mu_y = (int)(round(c_j + v_eff[s * L2 + i].second * (d_mu_vel_prefactor[s])) + L) % L;

            // Advection stencil bounds
            /**  RECALL: d_sigvD_Scalebound is __constant__ int2 array of size CuSpNV
             and corresponds to sig_vD_ScaledBounds[s] from CPU code **/
            int min_p_vx = mu_x - d_sigvD_ScaleBounds[s].x; int max_p_vx = mu_x + d_sigvD_ScaleBounds[s].x;
            int min_q_vy = mu_y - d_sigvD_ScaleBounds[s].y; int max_q_vy = mu_y + d_sigvD_ScaleBounds[s].y;

            // Update diffusion term
            for (int p = min_p; p <= max_p; ++p)
            {
                // Apply periodic boundary conditions
                int dx = std::abs(p - c_i); // Calculate distance from mean //int x = (p + L) % L; 
                //int pLmodL = ((p + L) % L)*L; // Precompute p index for efficiency 
                int pLmodL = p; pLmodL += (pLmodL < 0)*L; pLmodL -= (pLmodL >= L)*L; pLmodL *= L; // Precompute p index for efficiency
                for (int q = min_q; q <= max_q; ++q)
                {
                    int qmodL = q; qmodL += (qmodL < 0)*L; qmodL -= (qmodL >= L)*L; // Apply periodic BCs [Equivalent to (q + L) % L]
                    int dy = abs(q - c_j); int i_new = pLmodL + qmodL;
                    //C_t[x * Nx + y] += (1 - a1) * (1 - a1) * C[i * Nx + j] * gaussian_stencil_D[dx] * gaussian_stencil_D[dy];
                    if constexpr (useConstantGaussMemory) { // Using __constant__ memory (d_gaussian_stencil_*) for Gaussian stencils
                        atomicAdd(&Rho_tsar[sL2 + i_new], gamma_2_si * Rho_t[sL2_i]
                            * (d_gaussian_stencil_D[d_size_gauss_D[s] + dx] * d_gaussian_stencil_D[d_size_gauss_D[s] + dy]));
                    } 
                    else { // Using passed pointers for Gaussian stencils
                    atomicAdd(&Rho_tsar[s * L2 + i_new], (gamma[s * L2 + i]) * (gamma[s * L2 + i]) * Rho_t[s * L2 + i]
                        * (gaussian_stencil_D[d_size_gauss_D[s] + dx] * gaussian_stencil_D[d_size_gauss_D[s] + dy]));
                    }
                }
                // Next update d_C_t[] using joint probability of advection and diffusion terms
                for (int q = min_q_vy; q <= max_q_vy; ++q)
                {
                    int qmodL = q; qmodL += (qmodL < 0)*L; qmodL -= (qmodL >= L)*L; //[Equivalent to (q + L) % L] but faster
                    int dy = std::abs(q - mu_y); int i_new = pLmodL + qmodL;
                    //C_t[x][y] += (1-a1)*(a1)*C[i][j]*(gaussian_stencil_D[dx]*gaussian_stencil_vy[dy]);
                    if constexpr (useConstantGaussMemory) { // Using __constant__ memory (d_gaussian_stencil_*) for Gaussian stencils
                        atomicAdd(&Rho_tsar[sL2 + i_new], (gamma_si - gamma_2_si) * Rho_t[sL2_i] *
                            (d_gaussian_stencil_D[d_size_gauss_D[s] + dx] * d_gaussian_stencil_VXY[d_size_gauss_VXY[s] + dy]));
                    } 
                    else { // Using passed pointers for Gaussian stencils
                    atomicAdd(&Rho_tsar[s * L2 + i_new], (gamma[s * L2 + i]) * (1 - gamma[s * L2 + i]) * Rho_t[s * L2 + i] *
                        (gaussian_stencil_D[d_size_gauss_D[s] + dx] * gaussian_stencil_VXY[d_size_gauss_VXY[s] + dy]));
                    }
                }
            }

            // Update advection term
            for (int p = min_p_vx; p <= max_p_vx; ++p)
            {
                //int x = (p + L) % L;
                int dx = abs(p - mu_x); //int pLmodL = ((p + L) % L)*L; // Precompute p index for efficiency 
                int pLmodL = p; pLmodL += (pLmodL < 0)*L; pLmodL -= (pLmodL >= L)*L; pLmodL *= L; // Precompute p index for efficiency
                // Probability of advection in x-axis and diffusion in y-direction
                for (int q = min_q; q <= max_q; ++q)
                {
                    int qmodL = q; qmodL += (qmodL < 0)*L; qmodL -= (qmodL >= L)*L; //[Equivalent to (q + L) % L]
                    int dy = std::abs(q - c_j); int i_new = pLmodL + qmodL;
                    //C_t[x][y] += (a1)*(1-a1)*C[i][j]*(gaussian_stencil_vx[dx]*gaussian_stencil_D[dy]);
                    if constexpr (useConstantGaussMemory) { // Using __constant__ memory (d_gaussian_stencil_*) for Gaussian stencils
                        atomicAdd(&Rho_tsar[sL2 + i_new], (gamma_si - gamma_2_si) * Rho_t[sL2_i]
                            * (d_gaussian_stencil_VXY[d_size_gauss_VXY[s] + dx] * d_gaussian_stencil_D[d_size_gauss_D[s] + dy]));
                    } 
                    else { // Using passed pointers for Gaussian stencils
                    atomicAdd(&Rho_tsar[s * L2 + i_new], (1 - gamma[s * L2 + i]) * (gamma[s * L2 + i]) * Rho_t[s * L2 + i]
                        * (gaussian_stencil_VXY[d_size_gauss_VXY[s] + dx] * gaussian_stencil_D[d_size_gauss_D[s] + dy]));
                    }
                }

                // Update pure advection term joint probability
                for (int q = min_q_vy; q <= max_q_vy; ++q)
                {
                    int qmodL = q; qmodL += (qmodL < 0)*L; qmodL -= (qmodL >= L)*L; //[Equivalent to (q + L) % L] but faster
                    int dy = abs(q - mu_y); int i_new = pLmodL + qmodL;
                    //C_t[x * Nx + y] += a1 * a1 * C[i * Nx + j] * gaussian_stencil_vx[dx] * gaussian_stencil_vy[dy];
                    if constexpr (useConstantGaussMemory) { // Using __constant__ memory (d_gaussian_stencil_*) for Gaussian stencils
                        atomicAdd(&Rho_tsar[sL2 + i_new], gammap_2_si * Rho_t[sL2_i]
                            * (d_gaussian_stencil_VXY[d_size_gauss_VXY[s] + dx] * d_gaussian_stencil_VXY[d_size_gauss_VXY[s] + dy]));
                    } 
                    else { // Using passed pointers for Gaussian stencils
                    atomicAdd(&Rho_tsar[s * L2 + i_new], (1 - gamma[s * L2 + i]) * (1 - gamma[s * L2 + i]) * Rho_t[s * L2 + i]
                        * (gaussian_stencil_VXY[d_size_gauss_VXY[s] + dx] * gaussian_stencil_VXY[d_size_gauss_VXY[s] + dy]));
                    }
                }
            }
        }
    }
}

//Function that reports values in Gaussian stencils.
template <typename T, bool useConstantGaussMemory>
__global__ void reportGaussianStencils(T* gaussian_stencil_D /**  = nullptr */, T* gaussian_stencil_VXY /**  = nullptr */)
{
    printf("Gaussian Stencil D: \n");
    for (int s = 0; s < CuSpNV; s++)
    {
        //printf("Species %d: ", s);
        //int size_D = d_size_gauss_D[s];
        printf("For species %d, size_gauss_D = %d, d_sigD_ScaleBounds = %d\n", s, d_size_gauss_D[s], d_sigD_ScaleBounds[s]);
        for (int i = 0; i <= d_sigD_ScaleBounds[s]; i++) 
        {
            if constexpr (useConstantGaussMemory) { // Using __constant__ memory (d_gaussian_stencil_*) for Gaussian stencils
                printf("%f\t", d_gaussian_stencil_D[d_size_gauss_D[s] + i]);
            } 
            else { // Using passed pointers for Gaussian stencils
            printf("%f\t", gaussian_stencil_D[d_size_gauss_D[s] + i]);
            }
        }
        printf("\n");
    }
    // Next as a sanity check, print all Gaussian stencil D values using the size of gaussian_stencil_D
    printf("All Gaussian Stencil D values: \n");
    for (int i = 0; i < d_size_gauss_D[CuSpNV - 1] + d_sigD_ScaleBounds[CuSpNV - 1] + 1; i++) {
         if constexpr (useConstantGaussMemory) printf("%f\t", d_gaussian_stencil_D[i]); 
         else printf("%f\t", gaussian_stencil_D[i]);
         // Reporting all values from either constant memory or passed pointer.
    }

    printf("Gaussian Stencil VXY: \n");
    for (int s = 0; s < CuSpNV; s++) {
        //printf("Species %d: ", s);
        //int size_VXY = d_size_gauss_VXY[s];
        printf("For species %d, size_gauss_VXY = %d, d_sigvD_ScaleBounds = ( %d , %d )\n", s, d_size_gauss_VXY[s],
            d_sigvD_ScaleBounds[s].x, d_sigvD_ScaleBounds[s].y);
        for (int i = 0; i <= d_sigvD_ScaleBounds[s].x; i++) {
            if constexpr (useConstantGaussMemory) printf("%f\t", d_gaussian_stencil_VXY[d_size_gauss_VXY[s] + i]); 
            else printf("%f\t", gaussian_stencil_VXY[d_size_gauss_VXY[s] + i]);
        }
        printf("\n");
    }

    // Next as a sanity check, print all Gaussian stencil VXY values using the size of gaussian_stencil_VXY
    printf("All Gaussian Stencil VXY values: \n");
    for (int i = 0; i < d_size_gauss_VXY[CuSpNV - 1] + d_sigvD_ScaleBounds[CuSpNV - 1].y + 1; i++) {
        if constexpr (useConstantGaussMemory) printf("%f\t", d_gaussian_stencil_VXY[i]);
        else  printf("%f\t", gaussian_stencil_VXY[i]);
        // Reporting all values from either constant memory or passed pointer.
    }
    printf("\n");
}

template <typename T>
void reportGaussianStencilValues(T* d_gaussian_stencil_D, T* d_gaussian_stencil_VXY)
{
    reportGaussianStencils <T, false> << <1, 1 >> > (d_gaussian_stencil_D, d_gaussian_stencil_VXY);
}

// Explicit instantiation for types so that linker can find the definitions in host code.
template void reportGaussianStencilValues<float>(float*, float*);
template void reportGaussianStencilValues<double>(double*, double*);

// Function that reports all the values stored in __constant__ memory
__global__ void reportConstantMemory_AdvTerms()
{
    printf("Printing values stored in __constant__ memory \n");
    printf("(Sp, d_size_gauss_D): \t");
    for (int i = 0; i < CuSpNV; i++) {
        printf("( %d , %d )\t", i, d_size_gauss_D[i]);
    }
    printf("\n");

    printf("d_size_gauss_VXY: \t");
    for (int i = 0; i < CuSpNV; i++) {
        printf("( %d , %d )\t", i, d_size_gauss_VXY[i]);
    }
    printf("\n");

    printf("d_mu_vel_prefactor: ");
    for (int i = 0; i < CuSpNV; i++) {
        printf("%f\t", d_mu_vel_prefactor[i]);
    }
    printf("\n");

    printf("d_sigD_ScaleBounds: ");
    for (int i = 0; i < CuSpNV; i++) {
        printf("%d\t", d_sigD_ScaleBounds[i]);
    }
    printf("\n");

    printf("(Sp, d_sigvD_ScaleBounds.X, d_sigvD_ScaleBounds.Y): ");
    for (int i = 0; i < CuSpNV; i++) {
        printf("[ %d , %d , %d ]\t", i, d_sigvD_ScaleBounds[i].x, d_sigvD_ScaleBounds[i].y);
    }
    printf("\n");


}

void reportConstantMemory()
{
    reportConstantMemory_AdvTerms << <1, 1 >> > ();
    if (use_constmem_gaussstencils) reportGaussianStencils <float, true> << <1, 1 >> > (); // Assuming float type for Gaussian stencils in constant memory.
}


// Function to launch the FKE- Advection Diffusion CUDA kernel.
// NOTE: DEFAULT ARGUMENTS ASSUME USE OF CONSTANT MEMORY FOR GAUSSIAN STENCILS!
void launchAdvDiffKernel_MultiSp(double* Rho_t, double* Rho_tsar, double* gamma, std::pair<double, double>* v_eff, int L,
    double* gaussian_stencil_D /**  = nullptr */, double* gaussian_stencil_VXY /**  = nullptr */)
{
    // Calculate grid and block sizes
    int blockSize = 256;
    int numBlocks = (L * L + blockSize - 1) / blockSize;

    // Launch the CUDA kernel with appropriate template instantiation
    if (use_constmem_gaussstencils) {
        updateMovement<true> << <numBlocks, blockSize >> > (Rho_t, Rho_tsar,  gamma, v_eff, L);
    }
    else {
        updateMovement<false> << <numBlocks, blockSize >> > (Rho_t, Rho_tsar, gamma, v_eff, L, gaussian_stencil_D, gaussian_stencil_VXY);
    }

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) 
    {   fprintf(stderr, "CUDA kernel error: %s\n", cudaGetErrorString(err));    }
    
    //cudaDeviceSynchronize();
}