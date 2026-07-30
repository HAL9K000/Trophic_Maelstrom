#include <string>
#include <stdlib.h>
#include <stdio.h>

inline const int Sp = 2; //Total number of species in the system.
inline const std::string frame_header = "a_c,  x,  P(x; t), G(x; t) \n"; 
//Header for frame files.
inline const std::string gammaheader = "a_c,  x,  GAM[P(x; t)], GAM[G(x; t)] \n"; 
//Header for gamma frame.
inline const std::string gmmclusterfreq_header = "CLUSTER_SIZE[P(x; t)],CLUSTER_COUNTS[P(x; t)],CLUSTER_FREQ[P(x; t)]," 
			"CLUSTER_SIZE[G(x; t)],CLUSTER_COUNTS[G(x; t)],CLUSTER_FREQ[G(x; t)],\n"; 
//Header for GMM cluster frequency files.
inline std::string prelimheader = " a , r, L, t , <<P(x; t)>_x>_r, Var[<P(x; t)>_x]_r, # Surviving Runs P(x; t),"
					" # Active Sites P(x; t), <<G(x; t)>_x>_r, Var[<G(x; t)>_x]_r, # Surviving Runs G(x; t), # Active Sites G(x; t), \n";
//Header for preliminary files.

inline std::string movmheader = " a , r, L, t , <GAM[P(x; t)]>_x, Var[<GAM[P(x; t)]>_x], <vx[P(x; t)]>_x,  <vy[P(x; t)]>_x,"
					" <GAM[G(x; t)]>_x, Var[<GAM[G(x; t)]>_x], <vx[G(x; t)]>_x,  <vy[G(x; t)]>_x, \n";
//Header for movement files.

inline std::string critgmmheader = " a , r, L, t , CLUS_DENSITY{P(x; t)}, OCCUPIED_SITE_FRAC{P(x; t)},"
					" MEAN_SQ_DIST_LCENTER{P(x; t)}, MEAN_SQ_DIST_UCOM{P(x; t)}, MEAN_CLUSTER_SIZE{P(x; t)}, VAR_CLUSTER_SIZE{P(x; t)}, NUM_CLUSTERS{P(x; t)},"
	" CLUS_DENSITY{G(x; t)}, OCCUPIED_SITE_FRAC{G(x; t)}, MEAN_SQ_DIST_LCENTER{G(x; t)}, MEAN_SQ_DIST_UCOM{G(x; t)}, MEAN_CLUSTER_SIZE{G(x; t)}, VAR_CLUSTER_SIZE{G(x; t)}, NUM_CLUSTERS{G(x; t)}, \n";
//Header for GMM classified prelim files.

#if defined(__CUDACC__) || defined(BARRACUDA)
// Creating constexpr variables to be used in the CUDA kernels.
inline const int CuSpB = 2; // Sp, biotic species in the system.
inline const int CuSpNV = 1; // Sp - 1, biotic species in the system.
#endif