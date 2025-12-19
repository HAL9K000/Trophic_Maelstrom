import numpy as np
import matplotlib.pyplot as plt
import pandas as pan
import seaborn as sea
import cmath
import math
import numpy.linalg as la
import scipy.optimize as opt


# Define global variables.
epsilon = 1e-4; # Tolerance for checking if a number is close to zero.

c = 10000; gmax = 0.05*pow(10, -3.0)/24.0; d = 0.25/24.0; alpha =0.2/24.0; W0 = 0.2; rW = 0.2/24.0; # From Bonachela et al 2015
dx= 0.1 ; dt = 0.1; #From Bonachela et al 2015 (in km)
d0 = 0.00025/24.0; d1=0.0298; d2 = 0.00025/24.0; d3= 0.025/24.0; #From Bonachela et al 2015 (in km^2/hr)
k0= 0; k1 = 5; k2 =5000;

# Scaling factors for the grazer parameters.
mj_scaling = 1; aij_scaling = 25; ej_scaling = 1;
init_g_scaling = 1e-7*mj_scaling; # Scaling factor for initial grazer density in kg/m^2

# Scaling factors for the predator parameters.
mm_scaling = 1; ajm_scaling = 150; em_scaling = 1;
init_p_scaling = 1e-7*mm_scaling; # Scaling factor for initial predator density in kg/m^2

# Parmeters for grazer (Pawar & Co)
mG = 20.0; # Mass of producer in kg
mP = 100.0; # Mass of predator in kg
aij = aij_scaling*3.6*pow(10.0, -6.08)*pow(mG, -0.37); # in km^2/(hr kg)
hij = 1; #Handling time in hrs
ej =ej_scaling*0.45; mj = mj_scaling*0.061609*pow(mG, -0.25)/8760.0; # Mortality rate in hr^{-1}
# Parmeters for predator (Pawar & Co)
ajm = ajm_scaling*3.6*pow(10.0, -6.08)*pow(mP, -0.37); # in km^2/(hr kg)
hjm = 1; #Handling time in hrs
em =em_scaling*0.85; mm = mm_scaling*0.061609*pow(mP, -0.25)/8760.0; # Mortality rate in hr^{-1}


''' # Parameters for grazer (Kefi and Brose)
#Assuming mass of producer (mi) = 1 kg, mass of grazer (mj) = 10 kg, Mass of predator (mm) =100 kg
ej = 0.45; mj = 0.314*pow(10, -0.25); y =8; B0 =0.5
aij = (mj*y)/B0; hij = 1/(mj*y);
# Parmeters for predator (Kefi and Brose)
em = 0.85; mm = 0.314*pow(100, -0.25); y =8; B0 =0.5
ajm = (mj*y)/B0; 
hjm = 1/(mj*y); #Handling time
# '''

print("Parameters:\t")
print("c = %g, gmax = %g, d = %g, alpha = %g, W0 = %g, rW = %g" %(c, gmax, d, alpha, W0, rW))
print("cgmax = %g, d = %g" %(c*gmax, d))
print("K1 = %g, K2 = %g" %(k1, k2))
print("Grazer Parameters:\t")
print("aij = %g, hij = %g, ej = %g, mj = %g" %(aij, hij, ej, mj))
print("Predator Parameters:\t")
print("ajm = %g, hjm = %g, em = %g, mm = %g" %(ajm, hjm, em, mm))


# Define derived parameters/terms.

kappa = c*gmax - d;
Vstar_graz = mj/((ej- mj*hij)*aij);
Wstar_veg = d*k1/kappa;
Gstar_coex=  mm/((em- mm*hjm)*ajm);


print("Derived parameters:\t")
print("Kappa = %g, Gstar_coex= %g,  Vstar_graz = %g, Wstar_veg = %g" %(kappa, Gstar_coex, Vstar_graz, Wstar_veg))

R_trans = d*rW*k1/kappa;
omega= rW*k1 + gmax*Vstar_graz;
R_plus_min = 2*rW*k1 + omega + 2*rW*k1*math.sqrt(1 +omega/(rW*k1));
R_minus_min = 2*rW*k1 + omega - 2*rW*k1*math.sqrt(1 +omega/(rW*k1));


def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        exponent = int(math.floor(math.log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)

''' 
This script will plot certain functions of a variable "R" with defined functional form, across a range of values of R (specified by user).
Each of these functions will be plotted on the same graph, with the x-axis being the range of R values, 
and the y-axis being the value of the function.
Pretty colour palette will be used to distinguish between the functions (using Seaborn).
Some of the functions can return imaginary values for certain values of R, in which only the real part will be plotted (as a dashed line).
'''

# Define functional forms of the equilibrium values of the variables as functions of R for scipy.optimize.root to find the roots.
def saturating_exponential(x, A, b,c):
    return A*(1 - np.exp(-b*(x-c)))

def decaying_exponential(x, A, b, c, D):
    return A*np.exp(-b*(x-c)) + D

def linear(x, m, c):
    return m*x + c

def quadratic(x, A, B, C):
    return A*x**2 + B*x + C

def cubic(x, A, B, C, D):
    return A*x**3 + B*x**2 + C*x + D

# Define the functions to be plotted.
# First functions for equilbirum values with no vegetation and grazer.

def Vstar_noveg(R):
    return np.zeros_like(R)

def Wstar_noveg(R):
    W = R/rW
    return W

def Ostar_noveg(R):
    O = R/(alpha*W0)
    return O

def Gstar_noveg(R): 
    return np.zeros_like(R)

def Pstar_noveg(R): 
    return np.zeros_like(R)

# Next functions for equilibria with vegetation but no grazer.

def Vstar_veg(R):
    V = R - rW*Wstar_veg
    V *= (c/d)
    return V

def Wastar_veg(R):
    return Wstar_veg*np.ones_like(R)

def Ostar_veg(R):
    V = Vstar_veg(R)
    #num = V + k2; denom = V  + k2*W0;
    O = R/(alpha)*(V + k2)/(V+ k2*W0);# O *= R/(alpha)
    return O

def Gstar_veg(R): 
    return np.zeros_like(R)

def Pstar_veg(R): 
    return np.zeros_like(R)

# Next functions for equilibria with grazer and vegetation (coexistence) but no predator.

def Vegstar_graz(R):
    return Vstar_graz*np.ones_like(R)

def Ostar_graz(R):
    V = Vegstar_graz(R)
    #num = V + k2; denom = V  + k2*W0;
    #O = num/denom; O *= R/(alpha)
    O = R/(alpha)*(V + k2)/(V+ k2*W0); #O *= R/(alpha)
    return O

# Wstar_coex and Gstar_coex can be calculated from the other equilibria.
# They can also return imaginary values for certain R, so we will not plot them for those values.

def Wstar_graz(R):
    V = Vegstar_graz(R)
    b = R  - gmax*V -rW*k1; #b.astype(complex);
    riyal = b/(2*rW); #riyal.astype(complex);
    Wstar_plus = riyal + np.emath.sqrt(b**2 + 4*rW*k1*R)/(2*rW)
    Wstar_minus = riyal - np.emath.sqrt(b**2 + 4*rW*k1*R)/(2*rW)
    return Wstar_plus, Wstar_minus

def Gstar_graz(R, Wstar_plus, Wstar_minus):
    V = Vegstar_graz(R)
    common= (ej*c/mj)*R - (ej*d/mj)*V; common.astype(complex);
    Gstar_plus = common- (ej*rW*c/mj)*Wstar_plus
    Gstar_minus = common - (ej*rW*c/mj)*Wstar_minus
    return Gstar_plus, Gstar_minus

def Pstar_graz(R): 
    return np.zeros_like(R)

# Next functions for equilibria with grazer, vegetation and predator (coexistence).

def Grazstar_coex(R):
    return Gstar_coex*np.ones_like(R)

def Wstar_coex(R, Vstar_1, Vstar_2, Vstar_3):
   G= Grazstar_coex(R)
   Wstar_1 = R/rW - (d/(c*rW))*Vstar_1 -  (aij*Vstar_1*G)/(1 +aij*hij*Vstar_1)
   Wstar_2 = R/rW - (d/(c*rW))*Vstar_2 -  (aij*Vstar_2*G)/(1 +aij*hij*Vstar_2)
   Wstar_3 = R/rW - (d/(c*rW))*Vstar_3 -  (aij*Vstar_3*G)/(1 +aij*hij*Vstar_3)
   return Wstar_1, Wstar_2, Wstar_3

def Ostar_coex(R, Vstar_1, Vstar_2, Vstar_3):
    G= Grazstar_coex(R)
    Ostar_1 = R/(alpha)*(Vstar_1 + k2)/(Vstar_1 + k2*W0)
    Ostar_2 = R/(alpha)*(Vstar_2 + k2)/(Vstar_2 + k2*W0)
    Ostar_3 = R/(alpha)*(Vstar_3 + k2)/(Vstar_3 + k2*W0)
    return Ostar_1, Ostar_2, Ostar_3

def Predstar_coex(R, Vstar_1, Vstar_2, Vstar_3):
    G= Grazstar_coex(R)
    Pstar_1 = (em/mm)*G*((ej*aij*Vstar_1)/(1 + aij*hij*Vstar_1) - mj)
    Pstar_2 = (em/mm)*G*((ej*aij*Vstar_2)/(1 + aij*hij*Vstar_2) - mj)
    Pstar_3 = (em/mm)*G*((ej*aij*Vstar_3)/(1 + aij*hij*Vstar_3) - mj)
    return Pstar_1, Pstar_2, Pstar_3

# Define both the range of R values and the values of Vstar_coex for which the functions will be plotted.
# This is done by reading a csv file where the first column is the R values 
# and the second, third and fourth columns are the Vstar_coex values (Vstar_1, Vstar_2, Vstar_3) for which the functions will be plotted.
# The size of R and Vstar_coex should be the same.


def eigenval_check(R, Vstar, Wstar, Ostar, Gstar, Pstar):
    # Check the eigenvalues of the Jacobian matrix at the equilibria.
    # If the real part of the eigenvalues are negative, the equilibria are stable.
    #Input: R, Vstar, Wstar, Ostar, Gstar are arrays of equilibrium values of the variables.
    #Output: Returns a boolean array of the same length as R, where True indicates stability and False indicates instability.
    #Output: Also returns the eigenvalues of the Jacobian matrix at the equilibria as a 2D array of size (len(R), 4).

    # Define the Jacobian matrix.
    J = np.zeros((len(R), 5, 5), dtype=complex)
    J[:, 0, 0] = -d + c*gmax*Wstar/(Wstar + k1) -aij*Gstar/np.emath.power(1 + aij*hij*Vstar, 2)
    J[:, 0, 1] = c*gmax*k1*Vstar/np.emath.power(Wstar + k1, 2)
    J[:, 0, 2] = 0
    J[:, 0, 3] = -aij*Vstar/(1 + aij*hij*Vstar)
    J[:, 0, 4] = 0

    J[:, 1, 0] = alpha*k2*(1-W0)*Ostar/np.emath.power(Vstar + k2, 2) -gmax*Wstar/(Wstar + k1)
    J[:, 1, 1] = - gmax*k1*Vstar/np.emath.power(Wstar + k1, 2) - rW
    J[:, 1, 2] = alpha*(Vstar + k2*W0)/(Vstar + k2)
    J[:, 1, 3] = 0
    J[:, 1, 4] = 0

    J[:, 2, 0] = -alpha*k2*(1-W0)*Ostar/np.emath.power(Vstar + k2, 2)
    J[:, 2, 1] = 0
    J[:, 2, 2] = -alpha*(Vstar + k2*W0)/(Vstar + k2)
    J[:, 2, 3] = 0
    J[:, 2, 4] = 0

    J[:, 3, 0] = ej*aij*Gstar/np.emath.power(1 + aij*hij*Vstar, 2)
    J[:, 3, 1] = 0
    J[:, 3, 2] = 0
    J[:, 3, 3] = ej*aij*Vstar/(1 + aij*hij*Vstar) - mj -ajm*Pstar/np.emath.power(1 + ajm*hjm*Gstar, 2)
    J[:, 3, 4] = -ajm*Gstar/(1 + ajm*hjm*Gstar)

    J[:, 4, 0] = 0
    J[:, 4, 1] = 0
    J[:, 4, 2] = 0
    J[:, 4, 3] = em*ajm*Pstar/np.emath.power(1 + ajm*hjm*Gstar, 2)
    J[:, 4, 4] = em*ajm*Gstar/(1 + ajm*hjm*Gstar) - mm

    print("At R = ", R[0], "Vstar, Wstar, Ostar, Gstar: ", 
          Vstar[0], "  ", Wstar[0], "  ", Ostar[0], "  ", Gstar[0])
    print("Jacobian matrix at R =  " +str(R[0]))
    print(J[0, 0, :]); print(J[0, 1, :]); print(J[0, 2, :]); print(J[0, 3, :]); print("\n")
    print("At R = ", R[len(R)//2], "Vstar, Wstar, Ostar, Gstar: ",
            Vstar[len(R)//2], "  ", Wstar[len(R)//2], "  ", Ostar[len(R)//2], "  ", Gstar[len(R)//2])
    print("Jacobian matrix at R = " +str(R[len(R)//2]))
    print(J[len(R)//2, 0, :]); print(J[len(R)//2, 1, :]); print(J[len(R)//2, 2, :]); 
    print(J[len(R)//2, 3, :]); print("\n")

    # Calculate the eigenvalues of the Jacobian matrix.
    eigvals = np.zeros((len(R), 5), dtype=complex)
    for i in range(len(R)):
        try:
            eigvals[i] = la.eigvals(J[i])
        except la.LinAlgError:
            print("LinAlgError at R = ", R[i])
            print("Jacobian matrix: ", J[i, :, :])
            print("Vstar, Wstar, Ostar, Gstar, Pstar: ", Vstar[i], Wstar[i], Ostar[i], Gstar[i], Pstar[i])
            eigvals[i] = np.array([np.nan, np.nan, np.nan, np.nan], dtype=complex)

        #eigvals[i] = la.eigvals(J[i])
    
    # Check the stability of the equilibria.
    # If the real part of the eigenvalues are negative, the equilibria are stable.
    # Return a boolean array of the same length as R, where True indicates stability and False indicates instability.
    return np.all(eigvals.real < 0, axis=1), eigvals

try:
    R_Vcoex = pan.read_csv(f'Allo_a_{aij_scaling:g}_{ajm_scaling:g}_m_{mj_scaling:g}_{mm_scaling:g}_e_{ej_scaling:g}_{em_scaling:g}_{mG:g}_{mP:g}_0_0.3_3000_Trial_Analytical_soln.csv', 
                           converters={'Soln2': lambda s: complex(s.replace('*I', 'j')),
                                       'Soln3': lambda s: complex(s.replace('*I', 'j'))})
    
    '''R_Vcoex = pan.read_csv(f'Allo_{mG:g}_{mP:g}_0_0.3_1000_Trial_Analytical_solutions.csv', 
                           converters={'Soln2': lambda s: complex(s.replace('*I', 'j')),
                                       'Soln3': lambda s: complex(s.replace('*I', 'j'))})'''
    R = np.array(R_Vcoex['R']).astype(float)
    print(R)
    # Convert columns to numpy arrays, converting to float if string.
    Vstarcoex_1 = np.array(R_Vcoex['Soln1']); Vstarcoex_2 = np.array(R_Vcoex['Soln2']); Vstarcoex_3 = np.array(R_Vcoex['Soln3'])
    # Convert columns to numpy arrays, converting to float if string.
    Vstarcoex_1 = Vstarcoex_1.astype(complex); Vstarcoex_2 = Vstarcoex_2.astype(complex); Vstarcoex_3 = Vstarcoex_3.astype(complex)
    
    print("Read R and Vstar_coex from file.")
except Exception as e:
    print("Error reading file: %s" % e)
    print("Exiting.")
    exit(0)
finally:
    # Verify that the size of R and Vstar_coex is the same.
    if len(R) != len(Vstarcoex_1) or len(R) != len(Vstarcoex_2) or len(R) != len(Vstarcoex_3):
        print("Error: Size of R and Vstar_coex should be the same.")
        print("Exiting.")
        exit(1)
    else:
        print("Size of R and Vstar_coex is the same.")
        print("Continuing.")

def ss_plot():
    # Plot the various functions defined above over the range of R values.
    # Use Seaborn to make the plot look pretty.
    # Use a different colour for each function.
    # Use a legend and axis labels to make the plot easy to understand.
    # Use ax.text to label the functions on the plot.
    # Wstar_graz, Gstar_graz, Wstar_coex, Ostar_coex, Predstar_coex can return 
    #imaginary values for certain R, in which case only the real part will be plotted with dashed lines.

    print("R_min_+, R_min_-, R_c = (%g, %g, %g)" %(R_plus_min, R_minus_min, R_trans))

    # Create 5 subplots, 2 rows and 3 columns, with the plots on the lower row aligned to the center.
    fig, ax = plt.subplots(2, 3, figsize=(16, 8), sharex=True)
    sea.set(style='whitegrid')
    sea.set_palette('husl')
    fig.suptitle(r'SS Soln for R, Allo $m_j, m_m = %g ,  %g $  $a_{ij}, a_{jm} = %g ,  %g $ $ e_j, e_m = %g ,  %g$' 
                 %(mj_scaling, mm_scaling, aij_scaling, ajm_scaling, ej_scaling, em_scaling))

    #Extinction Equilibria
    ax[0, 0].plot(R, Vstar_noveg(R), label='$V^*_\mathrm{Ext}$', color='cornflowerblue')
    ax[0, 1].plot(R, Wstar_noveg(R), label='$W^*_\mathrm{Ext}$', color='cornflowerblue')
    ax[0, 2].plot(R, Ostar_noveg(R), label='$O^*_\mathrm{Ext}$', color='cornflowerblue')
    ax[1, 0].plot(R, Gstar_noveg(R), label='$G^*_\mathrm{Ext}$', color='cornflowerblue')
    ax[1, 1].plot(R, Pstar_noveg(R), label='$P^*_\mathrm{Ext}$', color='cornflowerblue')

    #Vegetation Equilibria
    ax[0, 0].plot(R, Vstar_veg(R), label='$V^*_\mathrm{Veg}$', color='lightseagreen')
    ax[0, 1].plot(R, Wastar_veg(R), label='$W^*_\mathrm{Veg}$', color='lightseagreen')
    ax[0, 2].plot(R, Ostar_veg(R), label='$O^*_\mathrm{Veg}$', color='lightseagreen')
    ax[1, 0].plot(R, Gstar_veg(R), label='$G^*_\mathrm{Veg}$', color='lightseagreen')
    ax[1, 1].plot(R, Pstar_veg(R), label='$P^*_\mathrm{Veg}$', color='lightseagreen')

    #Grazer Equilibria
    ax[0, 0].plot(R, Vegstar_graz(R), label='$V^*_\mathrm{Grazer}$', color='orchid')
    ax[0, 2].plot(R, Ostar_graz(R), label='$O^*_\mathrm{Grazer}$', color='orchid')
    ax[1, 1].plot(R, Pstar_graz(R), label='$P^*_\mathrm{Grazer}$', color='orchid')

    Wstar_plus, Wstar_minus = Wstar_graz(R)
    Gstar_plus, Gstar_minus = Gstar_graz(R, Wstar_plus, Wstar_minus)

    # Plot only the real part of Wstar_graz and Gstar_graz for R < R_minus_min where they are imaginary (AS DASHED LINES).

    mask1 = np.logical_and(np.isreal(Wstar_plus), R <= R_minus_min)
    ax[0, 1].plot(R[mask1], Wstar_plus[mask1], label='$W^*_\mathrm{Grazer+}$', color='violet')
    ax[0, 1].plot(R[~mask1], np.real(Wstar_plus[~mask1]), color='violet', linestyle='dashed')

    mask2 = np.logical_and(np.isreal(Wstar_minus), R <= R_minus_min)
    ax[0, 1].plot(R[mask2], Wstar_minus[mask2], label='$W^*_\mathrm{Grazer-}$', color='plum')
    ax[0, 1].plot(R[~mask2], np.real(Wstar_minus[~mask2]), color='plum', linestyle='dashed')

    #mask1 = np.isreal(Gstar_plus); mask2 = np.isreal(Gstar_minus);
    ax[1, 0].plot(R[mask1], Gstar_plus[mask1], label='$G^*_\mathrm{Grazer+}$', color='violet')
    ax[1, 0].plot(R[~mask1], Gstar_plus[~mask1],  color='violet', linestyle='dashed', alpha=0.5)
    ax[1, 0].plot(R[mask2], Gstar_minus[mask2], label='$G^*_\mathrm{Grazer-}$', color='plum')
    ax[1, 0].plot(R[~mask2], Gstar_minus[~mask2], color='plum', linestyle='dashed', alpha=0.5)

    # Plot only the real part of Wstar_graz and Gstar_graz for R > R_+_min where they are imaginary (AS DASHED LINES).

    mask1 = np.logical_and(np.isreal(Wstar_plus), R >= R_plus_min)
    mask2 = np.logical_and(np.isreal(Wstar_minus), R >= R_plus_min)
    ax[0, 1].plot(R[mask1], Wstar_plus[mask1], color='violet')
    ax[0, 1].plot(R[~mask1], np.real(Wstar_plus[~mask1]), color='violet', linestyle='dashed')
    ax[0, 1].plot(R[mask2], Wstar_minus[mask2], color='plum')
    ax[0, 1].plot(R[~mask2], np.real(Wstar_minus[~mask2]), color='plum', linestyle='dashed')

    #mask1 = np.isreal(Gstar_plus); mask2 = np.isreal(Gstar_minus);
    ax[1, 0].plot(R[mask1], Gstar_plus[mask1], color='violet')
    ax[1, 0].plot(R[~mask1], Gstar_plus[~mask1],  color='violet', linestyle='dashed', alpha=0.5)
    ax[1, 0].plot(R[mask2], Gstar_minus[mask2], color='plum')
    ax[1, 0].plot(R[~mask2], Gstar_minus[~mask2], color='plum', linestyle='dashed', alpha=0.5)

    #Coexistence Equilibria
    ax[1,0].plot(R, Grazstar_coex(R), label='$G^*_\mathrm{Coexist}$', color='hotpink')

    print(R[: 10])
    print("Vstar_coex: ")
    print(Vstarcoex_1[: 10])
    print(Vstarcoex_2[: 10])
    print(Vstarcoex_3[: 10])
    Wstar_1, Wstar_2, Wstar_3 = Wstar_coex(R, Vstarcoex_1, Vstarcoex_2, Vstarcoex_3)
    Ostar_1, Ostar_2, Ostar_3 = Ostar_coex(R, Vstarcoex_1, Vstarcoex_2, Vstarcoex_3)
    Pstar_1, Pstar_2, Pstar_3 = Predstar_coex(R, Vstarcoex_1, Vstarcoex_2, Vstarcoex_3)

    # Plot only the real part of Vstar_coex, Wstar_coex, Ostar_coex  and Predstar_coex for R where they are imaginary (AS DASHED LINES).

    mask1 = np.isreal(Vstarcoex_1); mask2 = np.isreal(Vstarcoex_2); mask3 = np.isreal(Vstarcoex_3);
    #If Vstar_coex is imaginary, so are Wstar_coex, Ostar_coex and Predstar_coex.

    ax[0, 0].plot(R[mask1], Vstarcoex_1[mask1], label='$V^*_\mathrm{Coexist1}$', color='indigo')
    ax[0, 0].plot(R[~mask1], np.real(Vstarcoex_1[~mask1]), color='indigo', linestyle='dashed')
    ax[0, 0].plot(R[mask2], Vstarcoex_2[mask2], label='$V^*_\mathrm{Coexist2}$', color='navy')
    ax[0, 0].plot(R[~mask2], np.real(Vstarcoex_2[~mask2]), color='navy', linestyle='dashed')
    ax[0, 0].plot(R[mask3], Vstarcoex_3[mask3], label='$V^*_\mathrm{Coexist3}$', color='hotpink')
    ax[0, 0].plot(R[~mask3], np.real(Vstarcoex_3[~mask3]), color='hotpink', linestyle='dashed')

    ax[0, 1].plot(R[mask1], Wstar_1[mask1], label='$W^*_\mathrm{Coexist1}$', color='indigo')
    ax[0, 1].plot(R[~mask1], np.real(Wstar_1[~mask1]), color='indigo', linestyle='dashed')
    ax[0, 1].plot(R[mask2], Wstar_2[mask2], label='$W^*_\mathrm{Coexist2}$', color='navy')
    ax[0, 1].plot(R[~mask2], np.real(Wstar_2[~mask2]), color='navy', linestyle='dashed')
    ax[0, 1].plot(R[mask3], Wstar_3[mask3], label='$W^*_\mathrm{Coexist3}$', color='hotpink')
    ax[0, 1].plot(R[~mask3], np.real(Wstar_3[~mask3]), color='hotpink', linestyle='dashed')

    ax[0, 2].plot(R[mask1], Ostar_1[mask1], label='$O^*_\mathrm{Coexist1}$', color='indigo')
    ax[0, 2].plot(R[~mask1], np.real(Ostar_1[~mask1]), color='indigo', linestyle='dashed')
    ax[0, 2].plot(R[mask2], Ostar_2[mask2], label='$O^*_\mathrm{Coexist2}$', color='navy')
    ax[0, 2].plot(R[~mask2], np.real(Ostar_2[~mask2]), color='navy', linestyle='dashed')
    ax[0, 2].plot(R[mask3], Ostar_3[mask3], label='$O^*_\mathrm{Coexist3}$', color='hotpink')
    ax[0, 2].plot(R[~mask3], np.real(Ostar_3[~mask3]), color='hotpink', linestyle='dashed')

    ax[1, 1].plot(R[mask1], Pstar_1[mask1], label='$P^*_\mathrm{Coexist1}$', color='indigo')
    ax[1, 1].plot(R[~mask1], np.real(Pstar_1[~mask1]), color='indigo', linestyle='dashed')
    ax[1, 1].plot(R[mask2], Pstar_2[mask2], label='$P^*_\mathrm{Coexist2}$', color='navy')
    ax[1, 1].plot(R[~mask2], np.real(Pstar_2[~mask2]), color='navy', linestyle='dashed')
    ax[1, 1].plot(R[mask3], Pstar_3[mask3], label='$P^*_\mathrm{Coexist3}$', color='hotpink')
    ax[1, 1].plot(R[~mask3], np.real(Pstar_3[~mask3]), color='hotpink', linestyle='dashed')

    # Convert the y-axis of ax[1,0] and ax[1,1] to log scale.
    ax[1,0].set_yscale('log'); ax[1,1].set_yscale('log'); ax[0, 1].set_yscale('log');
    # Plot the vertical lines for R_min_+, R_min_-, R_c
    for i in range(2):
      for j in range(3):
        axes = ax[i,j]
        axes.axvline(x=R_plus_min, color='tomato', linestyle='-.' )#, label='$R^+_{min}$')
        axes.axvline(x=R_minus_min, color='tomato', linestyle='-.')# , label='$R^-_{min}$')
        axes.axvline(x = R_trans, color='gray', linestyle=':', label= '$R_c$')
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axes.set_xlabel('R $(mm/hr)$')
        axes.set_ylabel('Equilibrium values $kg m^{-2}$')
        axes.tick_params(axis='both')

        # Find max and min values of y for each subplot.
        ymax = max(axes.get_ylim()); xmax = max(axes.get_xlim())
        axes.text(R_plus_min - xmax/50.0, ymax/5.0, '$R^{+}_{min}$', rotation=90, color='tomato')
        axes.text( R_minus_min - xmax/50.0, ymax/5.0, '$R^{-}_{min}$', rotation=90, color='tomato')
        axes.text(R_trans - xmax/50.0, ymax/10.0, '$R_{c}$', rotation=90, color='grey')

    # Delete the last subplot.
    for L in ax[0,2].get_xaxis().get_majorticklabels():
        L.set_visible(True)
    #ax[0,2].get_xaxis().get_majorticklabels().set_visible(True)
    fig.delaxes(ax[1,2])

    #Finally, name the subplots.
    ax[0,0].set_ylim(-1000, 200000); ax[0,0].set_title('$V^*$ vs R')
    ax[0,1].set_title('$W^*$ vs R')
    ax[0,2].set_ylim(bottom= -100); ax[0,2].set_title('$O^*$ vs R')
    ax[1,0].set_title('$G^*$ vs R')
    ax[1,1].set_title('$P^*$ vs R')

    plt.tight_layout()
    figure = plt.gcf() # get current figure
    figure.set_size_inches(16, 10)
    plt.savefig('aij_m_scalings/SS PLOT ALLO A --- %g %g M --- %g %g E --- %g %g.png' 
                %(aij_scaling, ajm_scaling, mj_scaling, mm_scaling, ej_scaling, em_scaling), dpi=1000)
    plt.show()
    plt.close()




def stable_ss_plot():
    # Plot the various functions defined above over the range of R values.
    # Use Seaborn to make the plot look pretty.
    # Use a different colour for each function.
    # Use a legend and axis labels to make the plot easy to understand.
    # Use ax.text to label the functions on the plot.
    # Wstar_coex and Gstar_coex can return imaginary values for certain R, so we will not plot them for those values.

    print("R_min_+, R_min_-, R_c = (%g, %g, %g)" %(R_plus_min, R_minus_min, R_trans))

    #R= np.array([0.04, 0.25])
    # Extinction equilibria
    Veq_noveg = Vstar_noveg(R); Weq_noveg = Wstar_noveg(R); Oeq_noveg = Ostar_noveg(R); Geq_noveg = Gstar_noveg(R); Peq_noveg = Pstar_noveg(R);
    # Vegetation only equilibria
    Veq_veg = Vstar_veg(R); Weq_veg = Wastar_veg(R); Oeq_veg = Ostar_veg(R); Geq_veg = Gstar_veg(R); Peq_veg = Pstar_veg(R);
    # Grazer equilibria
    Veq_graz = Vegstar_graz(R); Oeq_graz = Ostar_graz(R); Peq_graz = Pstar_graz(R);
    Wstar_plus, Wstar_minus = Wstar_graz(R)
    Gstar_plus, Gstar_minus = Gstar_graz(R, Wstar_plus, Wstar_minus)
    # True Coexistence equilibria.
    
    Geq_coex = Grazstar_coex(R)
    Weq_coex_1, Weq_coex_2, Weq_coex_3 = Wstar_coex(R, Vstarcoex_1, Vstarcoex_2, Vstarcoex_3)
    Oeq_coex_1, Oeq_coex_2, Oeq_coex_3 = Ostar_coex(R, Vstarcoex_1, Vstarcoex_2, Vstarcoex_3)
    Peq_coex_1, Peq_coex_2, Peq_coex_3 = Predstar_coex(R, Vstarcoex_1, Vstarcoex_2, Vstarcoex_3)


    # Check the stability of the equilibria using the eigenvalue check function.
    # If the real part of the eigenvalues are negative, the equilibria are stable.
    # We will plot the stable equilibria only.
    checkstable_noveg, eigenval_noveg = eigenval_check(R, Veq_noveg, Weq_noveg, Oeq_noveg, Geq_noveg, Peq_noveg)
    print("Stability profile of extinction equilibria: ")
    for i in range(0, len(R), len(R)//25):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_noveg[i], "with stability: ", checkstable_noveg[i], "\n\n")
    
    checkstable_veg, eigenval_veg = eigenval_check(R, Veq_veg, Weq_veg, Oeq_veg, Geq_veg, Peq_veg)
    print("Stability profile of vegetation equilibria: ")
    for i in range(0, len(R), len(R)//25):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_veg[i], "with stability: ", checkstable_veg[i], "\n\n")

    checkstable_graz_plus, eigenval_graz_plus = eigenval_check(R, Veq_graz, Wstar_plus, Oeq_graz, Gstar_plus, Peq_graz)
    print("Stability profile of Grazer equilibria (Wstar_plus): ")
    for i in range(0, len(R), len(R)//25):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_graz_plus[i], "with stability: ", checkstable_graz_plus[i], "\n\n")
    
    checkstable_graz_minus, eigenval_graz_minus = eigenval_check(R, Veq_graz, Wstar_minus, Oeq_graz, Gstar_minus, Peq_graz)
    print("Stability profile of Grazer equilibria (Wstar_minus): ")
    for i in range(0, len(R), len(R)//10):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_graz_minus[i], "with stability: ", checkstable_graz_minus[i], "\n\n")
    
    checkstable_graz = np.logical_or(checkstable_graz_plus, checkstable_graz_minus) 
    # Coexistence equilibria for V and O are stable if either Wstar_plus or Wstar_minus are stable.
    print("Number of true values in checkstable_graz_plus: ", np.sum(checkstable_graz_plus))
    print("Number of true values in checkstable_graz_minus: ", np.sum(checkstable_graz_minus))
    print("Number of true values in checkstable_graz: ", np.sum(checkstable_graz))

    checkstable_coex1, eigenval_coex1 = eigenval_check(R, Vstarcoex_1, Weq_coex_1, Oeq_coex_1, Geq_coex, Peq_coex_1)
    print("Stability profile of Coexistence equilibria 1: ")
    for i in range(0, len(R), len(R)//25):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_coex1[i], "with stability: ", checkstable_coex1[i], "\n\n")
    
    checkstable_coex2, eigenval_coex2 = eigenval_check(R, Vstarcoex_2, Weq_coex_2, Oeq_coex_2, Geq_coex, Peq_coex_2)
    print("Stability profile of Coexistence equilibria 2: ")
    for i in range(0, len(R), len(R)//25):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_coex2[i], "with stability: ", checkstable_coex2[i], "\n\n")
    
    checkstable_coex3, eigenval_coex3 = eigenval_check(R, Vstarcoex_3, Weq_coex_3, Oeq_coex_3, Geq_coex, Peq_coex_3)
    print("Stability profile of Coexistence equilibria 3: ")
    for i in range(0, len(R), len(R)//25):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_coex3[i], "with stability: ", checkstable_coex3[i], "\n\n")
    
    checkstable_coex = np.logical_or(checkstable_coex1, checkstable_coex2, checkstable_coex3)
    print("Number of true values in checkstable_coex1: ", np.sum(checkstable_coex1))
    print("Number of true values in checkstable_coex2: ", np.sum(checkstable_coex2))
    print("Number of true values in checkstable_coex3: ", np.sum(checkstable_coex3))
    print("Number of true values in checkstable_coex: ", np.sum(checkstable_coex))
    

    
    fig, ax = plt.subplots(2, 3, figsize=(16, 8), sharex=True)
    sea.set(style='whitegrid')
    sea.set_palette('husl')
    ax[0,0].plot(R[checkstable_noveg], Veq_noveg[checkstable_noveg], label='$V^*_\mathrm{Ext}$', alpha=0.75)
    ax[0,1].plot(R[checkstable_noveg], Weq_noveg[checkstable_noveg], label='$W^*_\mathrm{Ext}$', alpha=0.75)
    ax[0,2].plot(R[checkstable_noveg], Oeq_noveg[checkstable_noveg], label='$O^*_\mathrm{Ext}$', alpha=0.75)
    ax[1,0].plot(R[checkstable_noveg], Geq_noveg[checkstable_noveg], label='$G^*_\mathrm{Ext}$', alpha=0.75)
    ax[1,1].plot(R[checkstable_noveg], Peq_noveg[checkstable_noveg], label='$P^*_\mathrm{Ext}$', alpha=0.75)

    ax[0,0].plot(R[checkstable_veg], Veq_veg[checkstable_veg], label='$V^*_\mathrm{Veg}$', alpha=0.75)
    ax[0,1].plot(R[checkstable_veg], Weq_veg[checkstable_veg], label='$W^*_\mathrm{Veg}$', alpha=0.75)
    ax[0,2].plot(R[checkstable_veg], Oeq_veg[checkstable_veg], label='$O^*_\mathrm{Veg}$', alpha=0.75)
    ax[1,0].plot(R[checkstable_veg], Geq_veg[checkstable_veg], label='$G^*_\mathrm{Veg}$', alpha=0.75)
    ax[1,1].plot(R[checkstable_veg], Peq_veg[checkstable_veg], label='$P^*_\mathrm{Veg}$', alpha=0.75)
    

    ax[0,0].plot(R[checkstable_graz], Veq_graz[checkstable_graz], label='$V^*_\mathrm{Grazer}$', alpha=0.75)
    ax[0,2].plot(R[checkstable_graz], Oeq_graz[checkstable_graz], label='$O^*_\mathrm{Grazer}$', alpha=0.75)
    ax[1,1].plot(R[checkstable_graz], Peq_graz[checkstable_graz], label='$P^*_\mathrm{Grazer}$', alpha=0.75)

    ax[0,1].plot(R[checkstable_graz_plus], Wstar_plus[checkstable_graz_plus], label='$W^*_\mathrm{Grazer+}$', alpha=0.75)
    ax[0,1].plot(R[checkstable_graz_minus], Wstar_minus[checkstable_graz_minus], label='$W^*_\mathrm{Grazer-}$', alpha=0.75)

    ax[1,0].plot(R[checkstable_graz_plus], Gstar_plus[checkstable_graz_plus], label='$G^*_\mathrm{Grazer+}$', alpha=0.75)
    ax[1,0].plot(R[checkstable_graz_minus], Gstar_minus[checkstable_graz_minus], label='$G^*_\mathrm{Grazer-}$', alpha=0.75)

    ax[1,0].plot(R[checkstable_coex], Geq_coex[checkstable_coex], label='$G^*_\mathrm{Coexist}$', alpha=0.75)

    ax[0,0].plot(R[checkstable_coex1], Vstarcoex_1[checkstable_coex1], label='$V^*_\mathrm{Coexist1}$', alpha=0.75)
    ax[0,0].plot(R[checkstable_coex2], Vstarcoex_2[checkstable_coex2], label='$V^*_\mathrm{Coexist2}$', alpha=0.75)
    ax[0,0].plot(R[checkstable_coex3], Vstarcoex_3[checkstable_coex3], label='$V^*_\mathrm{Coexist3}$', alpha=0.75)

    ax[0,1].plot(R[checkstable_coex1], Weq_coex_1[checkstable_coex1], label='$W^*_\mathrm{Coexist1}$', alpha=0.75)
    ax[0,1].plot(R[checkstable_coex2], Weq_coex_2[checkstable_coex2], label='$W^*_\mathrm{Coexist2}$', alpha=0.75)
    ax[0,1].plot(R[checkstable_coex3], Weq_coex_3[checkstable_coex3], label='$W^*_\mathrm{Coexist3}$', alpha=0.75)

    ax[0,2].plot(R[checkstable_coex1], Oeq_coex_1[checkstable_coex1], label='$O^*_\mathrm{Coexist1}$', alpha=0.75)
    ax[0,2].plot(R[checkstable_coex2], Oeq_coex_2[checkstable_coex2], label='$O^*_\mathrm{Coexist2}$', alpha=0.75)
    ax[0,2].plot(R[checkstable_coex3], Oeq_coex_3[checkstable_coex3], label='$O^*_\mathrm{Coexist3}$', alpha=0.75)

    ax[1,1].plot(R[checkstable_coex1], Peq_coex_1[checkstable_coex1], label='$P^*_\mathrm{Coexist1}$', alpha=0.75)
    ax[1,1].plot(R[checkstable_coex2], Peq_coex_2[checkstable_coex2], label='$P^*_\mathrm{Coexist2}$', alpha=0.75)
    ax[1,1].plot(R[checkstable_coex3], Peq_coex_3[checkstable_coex3], label='$P^*_\mathrm{Coexist3}$', alpha=0.75)

    ax[0,1].set_yscale('log'); ax[1,1].set_yscale('log');
    ax[0,0].set_ylim(-1000, 200000); ax[0,0].set_title('Locally Stable $V^*$ vs R'); ax[0,1].set_title('Locally Stable $W^*$ vs R')
    ax[0,2].set_ylim(-20, 200); ax[0,2].set_title('Locally Stable $O^*$ vs R'); ax[1,0].set_title('Locally Stable $G^*$ vs R')
    ax[1,1].set_title('Locally Stable $P^*$ vs R')

    # Set title for overall plot
    fig.suptitle(r'Stable Equilibria for Rietkerk 3Sp Model, Allometric, $m_G = %g kg, m_P = %g kg$; $m_j, m_m = %g ,  %g $;  $a_{ij}, a_{jm} = %g ,  %g $; $ e_j, e_m = %g ,  %g$' 
                 %(mG, mP, mj_scaling, mm_scaling, aij_scaling, ajm_scaling, ej_scaling, em_scaling))
    
    # Finally add three dashed horizontal lines to show values of R_plus_min, R_minus_min and R_c ON ALL SUBPLOTS with appropriate axes.text labels (rotated 90)

    # Convert the y-axis of ax[1,0] and ax[1,1] to log scale.
    #ax[1,0].set_yscale('log'); ax[1,1].set_yscale('log'); ax[0, 1].set_yscale('log');
    # Plot the vertical lines for R_min_+, R_min_-, R_c
    for i in range(2):
      for j in range(3):
        axes = ax[i,j]
        axes.axvline(x=R_plus_min, color='tomato', linestyle='-.' )#, label='$R^+_{min}$')
        axes.axvline(x=R_minus_min, color='tomato', linestyle='-.')# , label='$R^-_{min}$')
        axes.axvline(x = R_trans, color='gray', linestyle=':', label= '$R_c$')
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axes.set_xlabel('R $(mm/hr)$')
        axes.set_ylabel('Equilibrium values $kg m^{-2}$')
        axes.tick_params(axis='both')

        # Find max and min values of y for each subplot.
        ymax = max(axes.get_ylim()); xmax = max(axes.get_xlim())
        axes.text(R_plus_min - xmax/50.0, ymax/5.0, '$R^{+}_{min}$', rotation=90, color='tomato')
        axes.text( R_minus_min - xmax/50.0, ymax/5.0, '$R^{-}_{min}$', rotation=90, color='tomato')
        axes.text(R_trans - xmax/50.0, ymax/10.0, '$R_{c}$', rotation=90, color='grey')

    # Delete the last subplot.
    for L in ax[0,2].get_xaxis().get_majorticklabels():
        L.set_visible(True)
    #ax[0,2].get_xaxis().get_majorticklabels().set_visible(True)
    fig.delaxes(ax[1,2])

    plt.tight_layout()
    figure = plt.gcf() # get current figure
    figure.set_size_inches(15, 10)
    plt.savefig('aij_m_scalings/STABLE ALLO_3Sp_ss_plot mG -- %g mP -- %g A --- %g %g M --- %g %g E --- %g %g.png' 
                %(mG, mP, aij_scaling, ajm_scaling, mj_scaling, mm_scaling, ej_scaling, em_scaling), dpi=1000)
    plt.show()
    plt.close()


def find_functional_forms_ss():
    '''Find the functional forms of the stable steady state equilibria for the Rietkerk model with 3 species as functions of R.
    This is done by first determining the equilibrium values of the model for the given R values, and then checking the stability of these equilibria.
    Thereafter, these values are fitted to a functional form using curve_fit from scipy.optimize.'''

    print("R_min_+, R_min_-, R_c = (%g, %g, %g)" %(R_plus_min, R_minus_min, R_trans))

    #R= np.array([0.04, 0.25])
    '''
    # Extinction equilibria
    Veq_noveg = Vstar_noveg(R); Weq_noveg = Wstar_noveg(R); Oeq_noveg = Ostar_noveg(R); Geq_noveg = Gstar_noveg(R); Peq_noveg = Pstar_noveg(R);
    # Vegetation only equilibria
    Veq_veg = Vstar_veg(R); Weq_veg = Wastar_veg(R); Oeq_veg = Ostar_veg(R); Geq_veg = Gstar_veg(R); Peq_veg = Pstar_veg(R);
    # Grazer equilibria
    Veq_graz = Vegstar_graz(R); Oeq_graz = Ostar_graz(R); Peq_graz = Pstar_graz(R);
    Wstar_plus, Wstar_minus = Wstar_graz(R)
    Gstar_plus, Gstar_minus = Gstar_graz(R, Wstar_plus, Wstar_minus)
    '''
    # True Coexistence equilibria.

    
    Geq_coex = Grazstar_coex(R)
    Weq_coex_1, Weq_coex_2, Weq_coex_3 = Wstar_coex(R, Vstarcoex_1, Vstarcoex_2, Vstarcoex_3)
    Oeq_coex_1, Oeq_coex_2, Oeq_coex_3 = Ostar_coex(R, Vstarcoex_1, Vstarcoex_2, Vstarcoex_3)
    Peq_coex_1, Peq_coex_2, Peq_coex_3 = Predstar_coex(R, Vstarcoex_1, Vstarcoex_2, Vstarcoex_3)

    # Check the stability of the equilibria using the eigenvalue check function.
    # If the real part of the eigenvalues are negative, the equilibria are stable.
    '''
    checkstable_noveg, eigenval_noveg = eigenval_check(R, Veq_noveg, Weq_noveg, Oeq_noveg, Geq_noveg, Peq_noveg)
    checkstable_veg, eigenval_veg = eigenval_check(R, Veq_veg, Weq_veg, Oeq_veg, Geq_veg, Peq_veg)
    checkstable_graz_plus, eigenval_graz_plus = eigenval_check(R, Veq_graz, Wstar_plus, Oeq_graz, Gstar_plus, Peq_graz)
    checkstable_graz_minus, eigenval_graz_minus = eigenval_check(R, Veq_graz, Wstar_minus, Oeq_graz, Gstar_minus, Peq_graz)
    checkstable_graz = np.logical_or(checkstable_graz_plus, checkstable_graz_minus) 
    '''
    # Coexistence equilibria for V and O are stable if either Wstar_plus or Wstar_minus are stable.
    checkstable_coex1, eigenval_coex1 = eigenval_check(R, Vstarcoex_1, Weq_coex_1, Oeq_coex_1, Geq_coex, Peq_coex_1)
    checkstable_coex2, eigenval_coex2 = eigenval_check(R, Vstarcoex_2, Weq_coex_2, Oeq_coex_2, Geq_coex, Peq_coex_2)
    checkstable_coex3, eigenval_coex3 = eigenval_check(R, Vstarcoex_3, Weq_coex_3, Oeq_coex_3, Geq_coex, Peq_coex_3)

    print("Stability profile of Coexistence equilibria 3: ")
    for i in range(0, len(R), len(R)//25):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_coex3[i], "with stability: ", checkstable_coex3[i], "\n\n")
    
    checkstable_coex = np.logical_or(checkstable_coex1, checkstable_coex2, checkstable_coex3)
    print("Number of true values in checkstable_coex1: ", np.sum(checkstable_coex1))
    print("Number of true values in checkstable_coex2: ", np.sum(checkstable_coex2))
    print("Number of true values in checkstable_coex3: ", np.sum(checkstable_coex3))
    print("Number of true values in checkstable_coex: ", np.sum(checkstable_coex))

    # Splice checkstable_coex3 for indices where corresponding R values are > R_trans.
    upper_R_trans_checkstable_coex3 = np.logical_and(checkstable_coex3, Peq_coex_3 > 0)
    # Splice checkstable_coex3 for indices where corresponding R values are > R_minus_min
    upper_R_minus_min_checkstable_coex3 = np.logical_and(checkstable_coex3, R > R_minus_min)


    checkstable_coex = np.logical_or(checkstable_coex1, checkstable_coex2, checkstable_coex3)
    # Coexistence 3 Eq for Vstar has a linear form. 
    popt_coex3_Vst, pcov_coex3_Vst = opt.curve_fit(linear, R[checkstable_coex3], Vstarcoex_3[checkstable_coex3])
    # Coexistence 3 Eq for Wstar has a linear form.
    #popt_coex3_Wst, pcov_coex3_Wst = opt.curve_fit(linear, R[checkstable_coex3], Weq_coex_3[checkstable_coex3], p0= [-3.36, 5.18])#, bounds=([-np.inf, 5.18], [-3.3, 5.2]))
    # Coexistence 3 Eq for Wstar has a decaying exponential form.
    popt_coex3_Wst, pcov_coex3_Wst = opt.curve_fit(decaying_exponential, R[upper_R_minus_min_checkstable_coex3], Weq_coex_3[upper_R_minus_min_checkstable_coex3], 
                                                   p0=[5, 0.1, R_minus_min, 1])#, bounds=( [ 0, 0, 0, 0], [np.inf, np.inf, np.inf, 1]))
    # Coexistence 1 for Ostar has a linear form.
    popt_coex3_Ost, pcov_coex3_Ost = opt.curve_fit(linear, R[checkstable_coex1], Oeq_coex_3[checkstable_coex1])
    # Coexistence 3 for Predator has a saturating exponential form, with guess values for A ~ 10^5.5, 
    popt_coex3_Pst, pcov_coex3_Pst = opt.curve_fit(saturating_exponential, R[upper_R_trans_checkstable_coex3], Peq_coex_3[upper_R_trans_checkstable_coex3], 
                                                   p0=[pow(10, 5.5), 0.1, R_trans], bounds=( [ 1000, 0, R_trans - epsilon], [pow(10, 6.5), 50, R_trans + epsilon]))

    print("Coexistence 3 Eq for Vstar: ", popt_coex3_Vst)
    print("Coexistence 3 Eq for Wstar: ", popt_coex3_Wst)
    print("Coexistence 1 Eq for Ostar: ", popt_coex3_Ost)
    print("Coexistence 3 Eq for Predator: ", popt_coex3_Pst)


    # Plot the fitted functions for the stable equilibria.
    fig, ax = plt.subplots(2, 2, figsize=(16, 12), sharex=True)
    sea.set(style='whitegrid')
    sea.set_palette('husl')
    ax[0,0].plot(R[checkstable_coex3], Vstarcoex_3[checkstable_coex3], label='$V^*_\mathrm{Coexist3}$', alpha=0.75)
    ax[0,0].plot(R[checkstable_coex3], linear(R[checkstable_coex3], *popt_coex3_Vst), label=r"$V^* = %g\times R + %g$" %(popt_coex3_Vst[0], popt_coex3_Vst[1]), alpha=0.75)
    ax[0,0].set_title(' $V^*$ vs R')
    ax[0,0].set_ylabel('$V^*$')

    ax[0,1].plot(R[checkstable_coex3], Weq_coex_3[checkstable_coex3], label='$W^*_\mathrm{Coexist3}$', alpha=0.75)
    #ax[0,1].plot(R[checkstable_coex3], linear(R[checkstable_coex3], *popt_coex3_Wst), label=r"$W^* = %g\times R + %g$" %(popt_coex3_Wst[0], popt_coex3_Wst[1]), alpha=0.75)
    ax[0,1].plot(R[checkstable_coex3], decaying_exponential(R[checkstable_coex3], *popt_coex3_Wst), label=r"$W^* = %g + %g e^{-%g(R - %g)}$" %(popt_coex3_Wst[3], popt_coex3_Wst[0], popt_coex3_Wst[1], popt_coex3_Wst[2]), alpha=0.75)
    ax[0,1].set_title('$W^*$ vs R')
    ax[0,1].set_ylabel('$W^*$')

    ax[1,0].plot(R[checkstable_coex1], Oeq_coex_3[checkstable_coex1], label='$O^*_\mathrm{Coexist3}$', alpha=0.75)
    ax[1,0].plot(R[checkstable_coex1], linear(R[checkstable_coex1], *popt_coex3_Ost), label=r"$O^* = %g\times R + %g$" %(popt_coex3_Ost[0], popt_coex3_Ost[1]), alpha=0.75)
    ax[1,0].set_title(' $O^*$ vs R')

    ax[1,1].plot(R[checkstable_coex3], Peq_coex_3[checkstable_coex3], label='$P^*_\mathrm{Coexist3}$', alpha=0.75)
    ax[1,1].plot(R[checkstable_coex3], saturating_exponential(R[checkstable_coex3], *popt_coex3_Pst), label=r"$P^* = $" +sci_notation(popt_coex3_Pst[0], 5, 5, 5)+ r"$(1 - e^{-%g(R - %g)})$" %(popt_coex3_Pst[1], popt_coex3_Pst[2]), alpha=0.75)
    ax[1,1].set_title('$P^*$ vs R')

    # Convert the y-axis of ax[0,1] and ax[1,1] to log scale.
    ax[1,1].set_yscale('log'); ax[0,1].set_yscale('log')
    # Plotting the vertical lines for R_min_+, R_min_-, R_c.
    for i in range(2):
      for j in range(2):
        axes = ax[i,j]
        axes.axvline(x=R_plus_min, color='tomato', linestyle='-.' )#, label='$R^+_{min}$')
        axes.axvline(x=R_minus_min, color='tomato', linestyle='-.')# , label='$R^-_{min}$')
        axes.axvline(x = R_trans, color='gray', linestyle=':', label= '$R_c$')
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axes.set_xlabel('R $(mm/hr)$')
        axes.set_ylabel('Equilibrium values $kg m^{-2}$')
        axes.tick_params(axis='both')

        # Find max and min values of y for each subplot.
        ymax = max(axes.get_ylim()); xmax = max(axes.get_xlim())
        axes.text(R_plus_min - xmax/50.0, ymax/5.0, '$R^{+}_{min}$', rotation=90, color='tomato')
        axes.text( R_minus_min - xmax/50.0, ymax/5.0, '$R^{-}_{min}$', rotation=90, color='tomato')
        axes.text(R_trans - xmax/50.0, ymax/10.0, '$R_{c}$', rotation=90, color='grey')

        # Set legends on lower right corner of the plot.
        axes.legend(loc='center left', bbox_to_anchor=(0.25, 0.2))

    
    figure = plt.gcf() # get current figure

    # Set title for overall plot
    fig.suptitle(r'Functional Forms of Stable Eq for Rietkerk 3Sp Model, Allometric, $m_G = %g kg, m_P = %g kg$; $m_j, m_m = %g ,  %g $;  $a_{ij}, a_{jm} = %g ,  %g $; $ e_j, e_m = %g ,  %g$' 
                 %(mG, mP, mj_scaling, mm_scaling, aij_scaling, ajm_scaling, ej_scaling, em_scaling))
    #plt.tight_layout()

    plt.savefig(f'aij_m_scalings/STANDARD Functional_Forms_ALLO_3Sp_ss_plot mG --- {mG:g} mP --- {mP:g} A --- {aij_scaling:g} {ajm_scaling:g} M --- {mj_scaling:g} {mm_scaling:g} E --- {ej_scaling:g} {em_scaling:g}.png', dpi=600)
    plt.show()
    plt.close()















ss_plot()
stable_ss_plot()
find_functional_forms_ss()













