#include <string>
#include <stdlib.h>
#include <stdio.h>

inline const int Sp = 1; //Total number of species in the system.
inline const std::string frame_header = "a_c,  x,  P(x; t) \n";
inline const std::string gammaheader ="\n"; 
//Header for frame files.
inline const std::string gmmclusterfreq_header = "CLUSTER_SIZE[P(x; t)],CLUSTER_COUNTS[P(x; t)],CLUSTER_FREQ[P(x; t)],\n";
//Header for GMM cluster frequency files.
inline std::string prelimheader = " a , r, L, t , <<P(x; t)>_x>_r, Var[<P(x; t)>_x]_r, # Surviving Runs P(x; t),"
					" # Active Sites P(x; t), \n";
//Header for preliminary files.

inline std::string movmheader ="\n";
inline std::string critgmmheader = " a , r, L, t , CLUS_DENSITY{P(x; t)}, OCCUPIED_SITE_FRAC{P(x; t)},"
	" MEAN_SQ_DIST_LCENTER{P(x; t)}, MEAN_SQ_DIST_UCOM{P(x; t)}, MEAN_CLUSTER_SIZE{P(x; t)}, VAR_CLUSTER_SIZE{P(x; t)}, NUM_CLUSTERS{P(x; t)}, \n";
//Header for GMM classified prelim files.

#if defined(__CUDACC__) || defined(BARRACUDA)
// Creating constexpr variables to be used in the CUDA kernels.
inline const int CuSpB = 1; // Sp, biotic species in the system.
inline const int CuSpNV = 0; // Sp - 1, biotic species in the system.
#endif