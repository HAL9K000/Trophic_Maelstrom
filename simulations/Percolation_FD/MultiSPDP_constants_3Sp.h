#include <string>
#include <stdlib.h>
#include <stdio.h>

inline const int Sp = 3; //Total number of species in the system.
inline const std::string frame_header = "a_c,  x,  P(x; t), G(x; t), Pr(x; t)\n"; 
//Header for frame files.
inline const std::string gammaheader = "a_c,  x,  GAM[P(x; t)], GAM[G(x; t)], GAM[Pr(x; t)] \n"; 
//Header for gamma frame.
inline std::string prelimheader = " a , r, L, t , <<P(x; t)>_x>_r, Var[<P(x; t)>_x]_r, # Surviving Runs P(x; t),"
					" # Active Sites P(x; t), <<G(x; t)>_x>_r, Var[<G(x; t)>_x]_r, # Surviving Runs G(x; t), # Active Sites G(x; t),"
					"<<Pr(x; t)>_x>_r, Var[<Pr(x; t)>_x]_r, # Surviving Runs Pr(x; t), # Active Sites Pr(x; t), \n";
//Header for preliminary files.

inline std::string movmheader = " a , r, L, t , <GAM[P(x; t)]>_x, Var[<GAM[P(x; t)]>_x], <vx[P(x; t)]>_x,  <vy[P(x; t)]>_x,"
					" <GAM[G(x; t)]>_x, Var[<GAM[G(x; t)]>_x], <vx[G(x; t)]>_x,  <vy[G(x; t)]>_x,"
					" <GAM[Pr(x; t)]>_x, Var[<GAM[Pr(x; t)]>_x], <vx[Pr(x; t)]>_x,  <vy[Pr(x; t)]>_x, \n";
//Header for movement files.

#if defined(__CUDACC__) || defined(BARRACUDA)
// Creating constexpr variables to be used in the CUDA kernels.
inline const int CuSpB = 3; // Sp, biotic species in the system.
inline const int CuSpNV = 2; // Sp - 1, biotic species in the system.
#endif