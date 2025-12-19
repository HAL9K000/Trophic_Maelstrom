// CUDA kernel for calculating gamma values with optimized memory access
// Designed for L=128, r_max~60 (neighborhood size ~11,000)

#define DEFINE_CUDA_CONSTANTS // Define the CUDA constant arrays only once through this .cu file in the header file
#include "rietkerk_bjork_basic.h"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cstdio>

// Constant memory for neighborhood offsets (faster than global memory)
// Max 64KB, so we can store ~8000 int2 pairs
__constant__ int2 d_neighborOffsets[8192];

// Device kernel: Each thread processes one spatial site
__global__ void calc_gamma_kernel(
    const double* __restrict__ d_Rho_t,      // [SpB][L*L] flattened
    double* __restrict__ d_gamma,            // [SpB][L*L] flattened
    const double* __restrict__ d_rho_inverse, // [SpB]
    const int* __restrict__ d_rfrac_indices,  // [SpB] species indices sorted by radius
    const double* __restrict__ d_rfrac_values, // [SpB] radius fractions
    const int* __restrict__ d_eff_sizes,      // [SpB] effective neighborhood sizes
    int L,
    int L2,
    int max_neighbors,
    double eps)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= L2) return;
    
    int c_i = i / L;
    int c_j = i % L;
    
    // Shared memory for neighborhood indices (reduce global memory access)
    extern __shared__ int s_neighbors[];
    
    // Compute neighborhood indices for this site (reuse across species)
    for (int k = threadIdx.x; k < max_neighbors; k += blockDim.x) {
        int2 offset = d_neighborOffsets[k];
        int nx = (c_i + offset.x + L) % L;
        int ny = (c_j + offset.y + L) % L;
        s_neighbors[k] = nx * L + ny;
    }
    __syncthreads();
    
    // Process each species
    for (int sp_idx = 0; sp_idx < 3; sp_idx++) 
    {
        int s = d_rfrac_indices[sp_idx];
        double fr = d_rfrac_values[sp_idx];
        int eff_size = d_eff_sizes[sp_idx];
        
        if (fr == 0.0 || eff_size < 1) 
        {
            if (i == 0 && s > 0) 
            {   d_gamma[s * L2 + i] = 0.0;  }
            continue;
        }
        
        // Check for zero densities
        const double* rho_s = d_Rho_t + s * L2;
        const double* rho_prev = d_Rho_t + (s - 1) * L2;
        
        if (s == 0) 
        {   continue;  }     // Vegetation doesn't need gamma
        
        double gamma_sum = 0.0;
        
        // Species 1 (Grazers)
        if (s == 1) 
        {
            const double* rho_pred = d_Rho_t + 2 * L2;
            double rho_inv_pred = d_rho_inverse[2];
            
            // Vectorized reduction with coalesced memory access
            #pragma unroll 8
            for (int k = 0; k < eff_size; k++) 
            {
                int neighbor_idx = s_neighbors[k];
                double rho_veg = rho_prev[neighbor_idx];
                double rho_pred_val = rho_pred[neighbor_idx];
                gamma_sum += rho_veg * (1.0 - rho_pred_val * rho_inv_pred);
            }
            
            double rho_inv_veg = d_rho_inverse[0];
            d_gamma[s * L2 + i] = fmin(1.0, fmax(0.0, (gamma_sum * rho_inv_veg) / eff_size));
        }
        // Species 2 (Predators)
        else if (s == 2) {
            #pragma unroll 8
            for (int k = 0; k < eff_size; k++) 
            {
                int neighbor_idx = s_neighbors[k];
                gamma_sum += rho_prev[neighbor_idx];
            }
            
            double rho_inv_grazer = d_rho_inverse[1];
            d_gamma[s * L2 + i] = fmin(1.0, fmax(0.0, (gamma_sum * rho_inv_grazer/ eff_size)));
        }
    }
}

// Optimized kernel when neighborhood fits entirely in constant memory
__global__ void calc_gamma_kernel_small(
    const double* __restrict__ d_Rho_t,
    double* __restrict__ d_gamma,
    const double* __restrict__ d_rho_inverse,
    int L,
    int L2,
    int eff_size_s1,
    int eff_size_s2,
    double eps)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= L2) return;
    
    int c_i = i / L;
    int c_j = i % L;
    
    double rho_inv_veg = d_rho_inverse[0];
    double rho_inv_grazer = d_rho_inverse[1];
    double rho_inv_pred = d_rho_inverse[2];
    
    // Species 1 (Grazers)
    double gamma_sum_s1 = 0.0;
    #pragma unroll 4
    for (int k = 0; k < eff_size_s1; k++) {
        int2 offset = d_neighborOffsets[k];
        int nx = (c_i + offset.x + L) % L;
        int ny = (c_j + offset.y + L) % L;
        int neighbor_idx = nx * L + ny;
        
        double rho_veg = d_Rho_t[neighbor_idx];
        double rho_pred = d_Rho_t[2 * L2 + neighbor_idx];
        gamma_sum_s1 += rho_veg * (1.0 - rho_pred * rho_inv_pred);
    }
    d_gamma[L2 + i] = fmin(1.0, fmax(0.0, 
        (gamma_sum_s1 * rho_inv_veg) / eff_size_s1));
    
    // Species 2 (Predators)
    double gamma_sum_s2 = 0.0;
    #pragma unroll 4
    for (int k = 0; k < eff_size_s2; k++) {
        int2 offset = d_neighborOffsets[k];
        int nx = (c_i + offset.x + L) % L;
        int ny = (c_j + offset.y + L) % L;
        int neighbor_idx = nx * L + ny;
        
        gamma_sum_s2 += d_Rho_t[L2 + neighbor_idx];
    }
    d_gamma[2 * L2 + i] = fmin(1.0, fmax(0.0, 
        (gamma_sum_s2 / eff_size_s2) * rho_inv_grazer));
}

// Host wrapper function
extern "C" void calc_gamma_3Sp_NonRefugia_CUDA(
    const std::vector<std::pair<int, int>>& centralNeighboringSites,
    const std::vector<std::vector<double>>& Rho_t,
    std::vector<std::vector<double>>& gamma,
    const double Rho_avg[3],
    const std::vector<std::pair<double, int>>& rfrac,
    const std::vector<std::pair<int, int>>& dtV_counter,
    double nVeg_frac,
    int r_max,
    int L)
{
    const int SpB = 3;
    const int L2 = L * L;
    const double eps = 1.0e-12;
    
    // Prepare neighborhood data
    int new_length = int(nVeg_frac * nVeg_frac * centralNeighboringSites.size());
    
    // Copy neighborhood offsets to constant memory
    std::vector<int2> h_offsets(new_length);
    for (int k = 0; k < new_length; k++) {
        h_offsets[k].x = centralNeighboringSites[k].first;
        h_offsets[k].y = centralNeighboringSites[k].second;
    }
    cudaMemcpyToSymbol(d_neighborOffsets, h_offsets.data(), 
                       new_length * sizeof(int2));
    
    // Flatten and copy Rho_t to device
    std::vector<double> h_Rho_t_flat(SpB * L2);
    for (int s = 0; s < SpB; s++) {
        std::copy(Rho_t[s].begin(), Rho_t[s].end(), 
                  h_Rho_t_flat.begin() + s * L2);
    }
    
    double *d_Rho_t, *d_gamma, *d_rho_inverse;
    cudaMalloc(&d_Rho_t, SpB * L2 * sizeof(double));
    cudaMalloc(&d_gamma, SpB * L2 * sizeof(double));
    cudaMalloc(&d_rho_inverse, SpB * sizeof(double));
    
    cudaMemcpy(d_Rho_t, h_Rho_t_flat.data(), SpB * L2 * sizeof(double), 
               cudaMemcpyHostToDevice);
    
    // Prepare rho_inverse
    double h_rho_inverse[SpB];
    for (int s = 0; s < SpB; s++) {
        h_rho_inverse[s] = (Rho_avg[s] >= eps) ? (1.0 / Rho_avg[s]) : (1.0 / eps);
    }
    cudaMemcpy(d_rho_inverse, h_rho_inverse, SpB * sizeof(double), 
               cudaMemcpyHostToDevice);
    
    // Calculate effective sizes
    int h_eff_sizes[SpB];
    double h_rfrac_values[SpB];
    int h_rfrac_indices[SpB];
    
    for (int sp_idx = 0; sp_idx < SpB; sp_idx++) 
    {
        int s = rfrac[sp_idx].second;
        double fr = rfrac[sp_idx].first;
        h_rfrac_values[sp_idx] = fr;
        h_rfrac_indices[sp_idx] = s;
        
        int eff_size = int(fr * fr * new_length);
        h_eff_sizes[sp_idx] = (eff_size < 1 && fr > 0) ? 1 : eff_size;
    }
    
    int *d_eff_sizes, *d_rfrac_indices;
    double *d_rfrac_values;
    cudaMalloc(&d_eff_sizes, SpB * sizeof(int));
    cudaMalloc(&d_rfrac_indices, SpB * sizeof(int));
    cudaMalloc(&d_rfrac_values, SpB * sizeof(double));
    
    cudaMemcpy(d_eff_sizes, h_eff_sizes, SpB * sizeof(int), 
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_rfrac_indices, h_rfrac_indices, SpB * sizeof(int), 
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_rfrac_values, h_rfrac_values, SpB * sizeof(double), 
               cudaMemcpyHostToDevice);
    
    // Launch kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (L2 + threadsPerBlock - 1) / threadsPerBlock;
    int sharedMemSize = new_length * sizeof(int);
    
    // Choose kernel based on neighborhood size
    if (new_length <= 4096) {
        // Use smaller, faster kernel
        calc_gamma_kernel_small<<<blocksPerGrid, threadsPerBlock>>>(
            d_Rho_t, d_gamma, d_rho_inverse, L, L2,
            h_eff_sizes[1], h_eff_sizes[2], eps);
    } else {
        calc_gamma_kernel<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(
            d_Rho_t, d_gamma, d_rho_inverse, d_rfrac_indices, d_rfrac_values,
            d_eff_sizes, L, L2, new_length, eps);
    }
    
    // Copy results back
    std::vector<double> h_gamma_flat(SpB * L2);
    cudaMemcpy(h_gamma_flat.data(), d_gamma, SpB * L2 * sizeof(double), 
               cudaMemcpyDeviceToHost);
    
    for (int s = 0; s < SpB; s++) 
    {
        std::copy(h_gamma_flat.begin() + s * L2, 
                  h_gamma_flat.begin() + (s + 1) * L2, 
                  gamma[s].begin());
    }
    
    // Cleanup
    cudaFree(d_Rho_t);
    cudaFree(d_gamma);
    cudaFree(d_rho_inverse);
    cudaFree(d_eff_sizes);
    cudaFree(d_rfrac_indices);
    cudaFree(d_rfrac_values);
}