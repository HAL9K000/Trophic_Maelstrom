#include <string>
#include <stdlib.h>
#include <stdio.h>

inline const int Sp = 3; //Total number of species in the system.
inline const std::string frame_header = "a_c,  x,  P(x; t), W(x; t), O(x; t) \n"; 
//Header for frame files.
inline std::string prelimheader = " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<P(x; t)>_x>_r, Var[<P(x; t)>_x]_r, # Surviving Runs P(x; t),"
					" # Active Sites P(x; t), \n";
//Header for preliminary files.
inline const std::string gammaheader="\n";
inline std::string movmheader ="\n";

#if defined(__CUDACC__) || defined(BARRACUDA)
// Creating constexpr variables to be used in the CUDA kernels.
inline const int CuSpB = 1; // Sp - 2, biotic species in the system.
inline const int CuSpNV = 0; // Sp - 3, biotic species in the system.
#endif