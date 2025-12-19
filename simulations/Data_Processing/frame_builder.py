'''
This script visualizes the data in the form of a video. It reads the data from the csv file and creates a video of the data.

The script is divided into the following sections:
1. Importing the required libraries
2. Reading the data from the csv file
   The data is stored in a directory structure as follows:
   {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a_val}_R_{R}.csv
   where terms in {} are variables that are provided as input to the script.
   In these csvs, the columns are as follows:
   a_c	  x	  P(x; t)	 G(x; t)	 Pr(x; t)	 W(x; t)	 O(x; t) 
   where the 2nd column onwards (indexed from 0) are the data columns and represent the species concentration (for each species) at each grid point.
   The columns are of length {g*g} (excluding the header row) where g is the grid size.

3. To visualise the data, the script reads the data from the csv files and stores them in a Pandas DataFrame.
   The data for each frame is stored in a separate DataFrame, and is visualised by reshaping the data into a 2D grid.
   The data is then plotted using Seaborn's heatmap function.
   Two types of plots are created:
    a. A plot of the concentration of each species.
    b. Subplots of the concentration of each species.
    These plots are saved as images in the {out_dir} and {outdir/Combined} directory as png files.

4. The images of the individual plots and the subplots are then seperately combined to create two videos using the cv2.videoWriter. 
   The videos are saved in the {out_dir/Videos} directory.
5. Unless otherwise specified, the images are then deleted from the {out_dir} and {outdir/Combined} directories.

'''

# Importing the required libraries
import os
import sys
import time
import math
import pandas as pan
import numpy as np
import seaborn as sea
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import cv2
import shutil
import argparse
import glob
import powerlaw as pwl
import regex as re
import fnmatch
from pathlib import Path
import adjustText as adjT
from collections import defaultdict

os.environ["OMP_NUM_THREADS"] = '7' # To avoid KMeans memory leak

import matplotlib.ticker as mtick
import matplotlib.colors as mcolours
from matplotlib.collections import LineCollection
from sklearn.cluster import KMeans
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


# User inputs.
SPB = 3; # Number of species in the model
#in_dir = f"../Data/DP/Reorganised_Frames/Stoc/{SPB}Sp/DPCORR-20_100_WIN-4000-DsB6OTH2/" # FOR DP-CORR
#in_dir = f"../Data/DP/Reorganised_Frames/Stoc/{SPB}Sp/DDParam_MFT/"
#in_dir = f"../Data/DP/Reorganised_Frames/Stoc/{SPB}Sp/DPKVScParam_GAMTST/"
in_dir = f"../Data/DP/Reorganised_Frames/Stoc/{SPB}Sp/DPParam_20_100_DSB6UNITY/"
#out_dir = f"../../Images/{SPB}Sp/DPParam_MFT/" # FOR DP
out_dir = f"../../Images/{SPB}Sp/DPParam_20_100_DSB6UNITY/" # FOR DP
#out_dir = f"../../Images/{SPB}Sp/DPKVScParam_SCALEMFT/" # FOR DP
#out_dir = f"../../Images/{SPB}Sp/DPKVScParam_20_100_GAMTST/" # FOR DP

#in_dir = f"../Data/LocRiet/Reorganised_Frames/Stoc/{SPB}Sp/ASCALE_1g_150g_LOC_DsC_TST/"
#out_dir = f"../../Images/{SPB}Sp/ASCALE_1g_150g_LOC_DsC_TST/"

#in_dir = f"../Data/Rietkerk/Reorganised_Frames/Stoc/{SPB}Sp/VSCALE_GAU_cU_20_100_FPE/"
#out_dir = f"../../Images/{SPB}Sp/VSCALE_GAU_cU_20_100/"

#in_dir = f"../Data/Rietkerk/Reorganised_Frames/Stoc/{SPB}Sp/GAM_GAU_20_100_TSTFD/"
#out_dir = f"../../Images/{SPB}Sp/FD_GAU_20_100_GAMTST/"

#in_dir = f"../Data/Rietkerk/Reorganised_Frames/Stoc/{SPB}Sp/ASCALE_20_100_HsX/"
#out_dir = f"../../Images/{SPB}Sp/ASCALE_20_100_HsX/"

#in_dir = f"../Data/Rietkerk/Reorganised_Frames/Stoc/{SPB}Sp/CORR-20_100_WIN-4000-HsXA0/"
#out_dir = f"../../Images/{SPB}Sp/CORR-20_100_WIN-4000/"
Path(out_dir).mkdir(parents=True, exist_ok=True)
#prefixes = ["DIC-NREF-1.1HI", "DIC-NREF-0.5LI", "DIC-NREF-0.1LI"]

#g = 128;  dP = 10000; Geq = 4.802  ; R_max= -1; #4.802    0.19208   0.096039         #0.2991     7.4774
#g = 128;  dP = 5; Geq = 5  ; R_max= -1;
g = 512;  dP = 1; Geq = 0 ; R_max= -1; #0.0032013  
# Geq (for LOC case) 2.2007 0.088027 [2 SP] 2.2786 0.091143
# Geq (for small ratio case) 3.0029 0.12012 [2 SP] 4.1569 0.16627
T_vals =[]; TS_vals =[];
a_vals = []    

#Making the output directories if they don't exist
#Path(out_dir + 'Videos/').mkdir(parents=True, exist_ok=True)
#Path(out_dir + 'Combined/').mkdir(parents=True, exist_ok=True)
#Path(out_dir + 'Singular/').mkdir(parents=True, exist_ok=True)
hex_list= [#['#D8F3DC', '#B7E4C7', '#95D5B2', '#74C69D', '#57CC99', '#38A3A5', '#006466', '#0B525B', '#235842', '#1B4332'],
            ['#DDBEA9', '#FFE8D6', '#D8F3DC', '#95D5B2', '#57CC99', '#38A3A5', '#006466', '#0B525B', '#235842', '#1B4332'], 
           ['#B7094C', '#A01A58', '#892B64', '#723C70', '#5C4D7D', '#455E89', '#2E6F95', '#1780A1', '#1987a9', '#0091AD'],
           ["#e6f2ff", "#ccdcff", "#b3beff", "#9a99f2", "#8b79d9", "#805ebf", "#6f46a6", "#60308c", "#511f73", "#431259"],            
           ["#eebba0", "#eaac8b", "#e88c7d", "#e77c76", "#e56b6f", "#b56576", "#915f78", "#6d597a", "#355070", "#283c55"],
           ["#ffff3f" ,"#eeef20", "#dddf00", "#d4d700", "#bfd200", "#aacc00", "#80b918", "#55a630", "#2b9348", "#007f5f"],
           ["#80ffdb", "#72efdd", "#64dfdf", "#56cfe1", "#48bfe3", "#4ea8de", "#5390d9", "#5e60ce", "#6930c3", "#7400b8"],
           ["#e3f2fd", "#bbdefb", "#90caf9", "#64b5f6", "#42a5f5", "#2196f3", "#1e88e5", "#1976d2", "#1565c0", "#0d47a1"],
           ["#e8b2c3", "#ffc4d6", "#ffa6c1", "#ff87ab", "#ff5d8f" ,"#ff4f84", "#ff447d", "#ff3d78", "#ff3471", "#ff2a6a"],           
           ['#B80068', "#BB0588", "#CC00AA", "#CC00CC", "#BC00DD", "#B100E8", "#A100F2", "#8900F2", '#802fe2', '#6A00F4']]
           #["#c4fff9", "#9ceaef", "#68d8d6", "#3dccc7", "#07beb8"]]

paired_hex_list =[ "#669900", "#99CC33", "#006699", "#3399CC", "#990066", "#CC3399",  "#A83535", "#D75959", "#FF9900", "#FFCC00"]
pastelpair_hex_list = ["#80b918", '#0091AD', "#2196f3", "#ff4f84", '#2E6F95', "#A100F2"]
gammapair_hex_list =['#2E6F95', '#ff4f84']
# Increasingly stormier blues
trueblues_hex_list = ["#64B5F6", "#4EA8DE", "#5390D9", "#4180B7", '#2E6F95',  '#176A7e', '#1B5A6B', '#224B60', '#283C55', '#2E2F4A']

float_list = [0, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 1]
vmax =[500, 35, 40]
cbar_width = 0.02  # Width of each colorbar
cbar_gap = 0.0025    # Gap between successive colorbars

def redirect_todevnull():
    # Redirects the stdout and stderr to /dev/null
    old_stdout = sys.stdout; old_stderr = sys.stderr
    sys.stdout = open(os.devnull, "w"); sys.stderr = open(os.devnull, "w");
    return old_stdout, old_stderr

def set_figsize(SPB):
    #Sets the figure size based on the number of species in the model.
    if SPB == 1:
        return (8, 6)
    else:
        return (8 + 6*(SPB-1), 6)
    
def set_figsizeGrid(p, q):
    # Sets the figure size to accomodate square subplots (i.e. 1x1, 2x2, 3x3 etc), based on the number of species in the model.
    if p == 1 and q == 1:
        return (8, 6)
    else:
        return (8 + 6*(q-1), 6 + 6*(p-1))
    
def set_figsizeSquare(SPB):
    # Sets the figure size to accomodate square subplots (i.e. 1x1, 2x2, 3x3 etc), based on the number of species in the model.
    if SPB == 1:
        return (8, 6)
    else:
        square_SPB_Size = int(math.ceil(math.sqrt(SPB))) # Calculate the square size based on the number of species
        return (8 + 6*(square_SPB_Size-1), 6 + 6*(square_SPB_Size-1))

def power_lawfit(x, A, exp):
        return A * (x ** (-exp))

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

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]


def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolours.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp


# Create a generic function that accepts a numeric value (float, integer, complex) and returns if it is an integer or not.
def is_integer(num):
    if isinstance(num, (int)):
        return True
    elif isinstance(num, float):
        return num.is_integer()
    elif isinstance(num, complex):
        return num.real.is_integer() and num.imag.is_integer()
    else:
        return False

''' # Summary of auto_guess_avals_tvals(...)
    NOTE: This only works when the indir structure is of the form: 
    indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/
    This function attempts to recreate a_vals (and DEPENDENT t_vals nested within) from crawling relevant txt files in the above directory structure.
    Default behaviour is to guess a common set of  a_vals and t_vals for all sub-directories in indir + PREFIX.
    When a_vals is provided, it returns a_vals (unmodified) and a guessed set of t_vals 
    that exists in indir + PREFIX + f"/L_{g}_a_{a}/ for all such a vals
    When t_vals is passed as None, it only returns a_vals.
'''
def auto_guess_avals_tvals(indir, PREFIX, a_vals = [], t_vals= []):

    if len(a_vals) == 0:
        try:
            with open(indir + PREFIX + "/a_vals.txt", 'r') as f:
                a_vals = [float(line.strip()) for line in f]
                print(f"List of a values: {a_vals}")
                a_vals = [int(a) if a.is_integer() else a for a in a_vals]
        except FileNotFoundError:
            print(f"a_vals.txt not found in the input directory {indir + PREFIX} and no a_vals specified. Terminating the program.")

            #a_vals = None; t_vals = None;
            return (None, None)
    #print(t_vals)        
    if len(t_vals) == 0:
        for a in a_vals:
            print(f"List of a values: {a_vals}")
            #a = int(a) if a.is_integer() else a
            a = int(a) if is_integer(a) else a
            try:
                with open(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_vals.txt", 'r') as f:
                    t_vals.extend([float(line.strip()) for line in f])
            except FileNotFoundError:
                print(f"T_vals.txt not found in the input directory {indir + PREFIX} and no T_vals specified for a = {a}.\n Skipping this a-value.")
                continue
        t_vals= [int(T) if T.is_integer() else T for T in t_vals]
        # Sort and remove duplicates
        t_vals = sorted(list(set(t_vals)))
        # End of T loop
    return (a_vals, t_vals)
    
''' # Summary of home_video(...)
This function creates a video of the png data in the directory dir which has the structure:
{filedir}/{prefixes}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/{pngformat}
The data is stored in the png files with the format specified in the pngformat parameter.
The video is saved in the out_dir directory with the name specified in the videoname parameter.
The function reads the data from the png files and stores them in a list of images.
Optionally, one can provide the relative path (from parendir) to the text files delineating a_vals, T_vals and maxR (if these values are not specified/ set to 0)
by the relative paths given by avals_txtdir, Tvals_txtdir and maxR_txtdir. 
The same applies for pngdir (relative path to png files for given a, T, dP, Geq etc values)
'''
def home_video(file_parendir, out_dir, prefixes, a_vals, T_vals, maxR, minR = 0, pngformat= "CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png", pngdir =  "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/",
               videoname="guess", video_relpath="guess", avals_txtdir = "{Pre}/a_vals.txt", Tvals_txtdir = "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_vals.txt",
                maxR_txtdir = "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt"):

    min_T = 0; max_T = 0; inputmax_R = maxR;
    if len(T_vals) > 0:
        min_T= min(T_vals); max_T = max(T_vals);# These are used to determine savepath and video name.

    for Pre in prefixes:
        img_array = []

        if len(a_vals) == 0:
            try:
                with open(file_parendir + avals_txtdir.format(Pre=Pre, g=g, dP= dP, Geq= Geq), 'r') as f:
                    a_vals = [float(line.strip()) for line in f]
            except FileNotFoundError:
                print(f"a_vals.txt not found in the input directory and no a_vals specified for {Pre}. skipping this prefix.")
                continue
        for a in a_vals:
            a = int(a) if float(a).is_integer() else a
            if len(T_vals) == 0:
                try:
                    with open(file_parendir + Tvals_txtdir.format(Pre=Pre, g=g, dP= dP, Geq= Geq, a=a), 'r') as f:
                        T_vals = [float(line.strip()) for line in f]
                    if(max(T_vals) > max_T):
                        max_T = max(T_vals)
                    if(min(T_vals) < min_T):
                        min_T = min(T_vals)
                except FileNotFoundError:
                    print(f"T_vals.txt not found in the input directory and no T_vals specified for a = {a}, {Pre}. Skipping this a-value.")
                    continue
            print(f"List of T values: {T_vals}")
            
            for T in T_vals:
                T = int(T) if float(T).is_integer() else T
                print(f"Processing a = {a}, T = {T} at dP = {dP}, Geq = {Geq}")
                if inputmax_R <= 0:
                    # This indicates a case where the number of frames is not known apriori and must be read from the data.
                    # This data is stored in the first line of the following text file:
                    # in_dir/PREFIX/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt
                    try:
                        with open(file_parendir + maxR_txtdir.format(Pre=Pre, g=g, dP= dP, Geq= Geq, a=a, T=T), 'r') as f:
                            maxR = int(f.readline().strip())

                    except FileNotFoundError:
                        print(f"maxR.txt not found in the input directory for a = {a}, T = {T}. Skipping this a-value.")
                        continue
                
                for R in range(minR, maxR):
                    try:
                        img = cv2.imread(file_parendir +pngdir.format(s=0, Pre=Pre, g=g, dP=dP, Geq=Geq, a=a, T=T, R=R, 
                                        Tmin = min_T, Tmax = max_T, amin = a_vals[0], amax = a_vals[-1])
                                         + pngformat.format(s=0, g=g, dP=dP, Geq=Geq,  a=a, T=T, R=R, Tmin = min_T, Tmax = max_T, amin = a_vals[0], amax = a_vals[-1]))
                        img_array.append(img)
                    except FileNotFoundError:
                        print(f"Image not found for a = {a}, T = {T}, R = {R}, Pref = {Pre}, in Path: \n"  + pngdir.format(s=0, Pre=Pre, L=g, dP=dP, Geq=Geq, a=a, T=T, R=R, amin = a_vals[0], amax = a_vals[-1]) 
                                         + pngformat.format(s=0, g=g, dP=dP, Geq=Geq,  a=a, T=T, R=R, amin = a_vals[0], amax = a_vals[-1]) + " \n Skipping....")
                        continue
                # End of R loop
            # End of T loop
        # End of a loop
        print(f"Number of images: {len(img_array)}")
        #Remove images from the list if they are None.
        img_array = [img for img in img_array if img is not None]
        print(f"Number of valid images detected: {len(img_array)}")
        if(len(img_array) > 0 and img_array[0] is not None):
            height, width, layers = img_array[0].shape; size = (width, height)
            if videoname == "guess":
                #Extract string preceding the first underscore in the pngformat string.
                Prefix = re.sub("/", "", Pre)
                if inputmax_R <= 0:
                    videoname = f"{Prefix}_Tmax_{max_T}.mp4"
                else:
                    videoname = f"{Prefix}_Tmax_{max_T}_Rmin_{minR}_Rmax_{maxR}.mp4"
            
            if video_relpath == "guess":
                # Modify video save path depending on number of a_vals and T_vals encountered.
                if(len(a_vals) > 1):
                    video_relpath = out_dir + f'{Pre}/Videos/'+ pngformat.split("_")[0] + f'/{a_vals[0]}_{a_vals[-1]}/' 
                else:
                    video_relpath = out_dir + f'{Pre}/Videos/'+ pngformat.split("_")[0] + f'/{a_vals[0]}/'
                if(len(T_vals) > 1):
                    video_relpath += f'Tmin_{min_T}_Tmax_{max_T}/'
                else:
                    video_relpath += f'Tmax_{max_T}/'
            else:
                video_relpath = out_dir + video_relpath.format(Pre=Pre, g=g, dP=dP, Geq=Geq, a=a_vals[-1], T= max_T, R= maxR,  
                                                               amin= a_vals[0], amax = a_vals[-1], Tmin = min_T, Tmax = max_T)
            #Path(video_path).mkdir(parents=True, exist_ok=True) # Create the directory if it doesn't exist.
            video_name = video_relpath + videoname
            os.makedirs(video_relpath, exist_ok=True)
            out = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'mp4v'), 1, size)

            print(f"Saving {videoname} to Relative Path :\n {video_relpath}")

            for i in range(len(img_array)):
                out.write(img_array[i])
            out.release()
        else:
            print(f"No images found for {Pre}. Skipping....")
        
    # End of Prefix loop


def frame_visualiser(in_dir, out_dir, PREFIX, a_vals, T_vals, maxR, minR =0, plt_gamma= False, delpng = False):
    # If a_vals and T_vals are not provided, the script reads the relevant text files and extracts the a_vals and T_vals from them.
    # List of a_vals provided {in_dir}/{PREFIX}/a_vals.txt
    # List of T_vals provided {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/T_vals.txt

    input_avals = a_vals; input_tvals = T_vals; input_Rmax = maxR
    print(f"Processing data for {PREFIX}...")
    #global a_vals, T_vals
    
    if input_avals == []:
        try:
            with open(in_dir + PREFIX + "/a_vals.txt", 'r') as f:
                a_vals = [float(line.strip()) for line in f]
                a_vals = [int(a) if a.is_integer() else a for a in a_vals]
        except FileNotFoundError:
            print("a_vals.txt not found in the input directory and no a_vals specified. Terminating the program.")
            return
    print(f"List of a values: {a_vals}")
    for a in a_vals:
        # Convert a and T to int format if they are integers.
        a = int(a) if float(a).is_integer() else a
        if input_tvals == []:
            try:
                print(f"Reading T_vals.txt for {PREFIX}, L= {g}, dP= {dP}, a = {a}...")
                with open(in_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_vals.txt", 'r') as f:
                    T_vals = [float(line.strip()) for line in f]
            except FileNotFoundError:
                print(f"T_vals.txt not found in the input directory {in_dir + PREFIX} and no T_vals specified for a = {a}. Skipping this a-value.")
                continue
        print(f"List of T values: {T_vals}")
        for T in T_vals:
            # Convert a and T to int format if they are integers. 
            T = int(T) if float(T).is_integer() else T
            print(f"Processing a = {a}, T = {T} at dP = {dP}, Geq = {Geq}")
            if input_Rmax <= 0:
                # This indicates a case where the number of frames is not known apriori and must be read from the data.
                # This data is stored in the first line of the following text file:
                # in_dir/PREFIX/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt
                try:
                    #print(f"Reading maxR.txt for {PREFIX}, L= {g}, dP= {dP}, a = {a}, T = {T}...")
                    with open(in_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt", 'r') as f:
                        maxR = int(f.readline().strip())
                except FileNotFoundError:
                    print(f"maxR.txt not found in the input directory for a = {a}, T = {T}. Skipping this a-value.")
                    continue
            
            savepath = out_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/"
            Path(savepath + "CombinedImg/").mkdir(parents=True, exist_ok=True)

            print(f"Creating images for a = {a}, T = {T} at dP = {dP}, Geq = {Geq}, maxR = {maxR}...")

            for R in range(minR, maxR):
            # Reading the data from the csv file
                # Convert a and T to int format if they are integers.
                try:
                    data = pan.read_csv(in_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/FRAME_T_{T}_a_{a}_R_{R}.csv", header= 0)
                except FileNotFoundError:
                    print(f"Data not found for a = {a}, T = {T}, R = {R}. Skipping....")
                    continue

                # Extracting the data columns
                data = data.iloc[:, 2:]
                # Check if second species is present in the data (i.e. the second column is not 0 for all rows)
                if SPB ==2 and data.iloc[:, 1].sum() == 0:
                    print(f"Skipping a = {a}, T = {T}, R = {R} as second species is absent.")
                    continue
                # If SPB = 3, check if second and third species are present in the data.
                '''
                if SPB == 3 and (data.iloc[:, 1].sum() == 0 or data.iloc[:, 2].sum() == 0):
                    print(f"Skipping a = {a}, T = {T}, R = {R} as second/third species is absent.")
                    continue
                #'''
                # If SPB = 3, check if second species are present in the data.
                '''
                if SPB == 3 and (data.iloc[:, 1].sum() == 0):
                    print(f"Skipping a = {a}, T = {T}, R = {R} as second species is absent.")
                    continue
                #'''

                fig, axs = plt.subplots(1, SPB, figsize = set_figsize(SPB))
                ''' # Individual plots as opposed to combined ones.
                for s in range(SPB):
                    data_Sp = np.array(data.iloc[:, s]).reshape(g, g)
                    # Plot data first individually and then add to subplot.
                    fig, ax =plt.subplots(figsize=(9, 6))
                    ax = sea.heatmap(data_Sp, ax=ax, cmap=get_continuous_cmap(hex_list[s]), cbar=True)
                    ax.set_title(r'$\rho_%g(t = %g)$,' % (s, T))
                    ax.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    ax.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    plt.savefig(out_dir + f"Singular/Rho_s_{s}_a_{a}_T_{T}_n_{R}.png")
                    plt.close()
                '''
                # Creating combinbed subplots
                #fig, axs = plt.subplots(1, SPB, figsize = set_figsize(SPB))
                for s in range(SPB):
                    ax = axs[s] if SPB > 1 else axs
                    data_Sp = np.array(data.iloc[:, s]).reshape(g, g)
                    ax = sea.heatmap(data_Sp , ax=ax, vmin=0, cmap=get_continuous_cmap(hex_list[s], float_list),  cbar=True)
                    ax.set_title(r'$\rho_%g(t = %g)$,' % (s, T))
                    ax.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    ax.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                #Set super title
                fig.suptitle(r'$R = %g$, $dP = %g$, $n =%d$' % (a, dP, R))
                
                plt.savefig(savepath + f"CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png")
                plt.close()

                '''
                # Creating subplots
                fig, axs = plt.subplots(1, 2, figsize = (12, 6))
                for s in range(2):
                    ax = axs[s]
                    data_Sp = np.array(data.iloc[:, s+3]).reshape(g, g)
                    ax = sea.heatmap(data_Sp , ax=ax, cmap=get_continuous_cmap(hex_list[s+3]),  cbar=True)
                    ax.set_title(r'$\rho_%g(t = %g)$,' % (s+3, T))
                    ax.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                    ax.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                #Set super title
                fig.suptitle(r'$R = %g$, $dP = %g$, $n =%d$' % (a, dP, R))
                
                plt.savefig(savepath + f"CombinedImg/WatConc_a_{a}_T_{T}_n_{R}.png")
                plt.close()
                '''
                
                if plt_gamma:
                    data_gamma = pan.read_csv(in_dir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/GAMMA_T_{T}_a_{a}_R_{R}.csv", header= 0)
                    data_gamma = data_gamma.iloc[:, 2:]
                    fig, axs = plt.subplots(2, SPB, figsize = (20, 12))

                    for s in range(SPB):
                        ax1 = axs[0, s]
                        ax2 = axs[1, s]
                        data_Sp = np.array(data.iloc[:, s]).reshape(g, g)
                        data_GAM = np.array(data_gamma.iloc[:, s]).reshape(g, g)
                        ax1 = sea.heatmap(data_Sp, ax=ax1, cmap=get_continuous_cmap(hex_list[s], float_list), cbar=True)#, vmax = vmax[s])
                        ax1.set_title(r'$\rho_%g(t = %g)$,' % (s, T))
                        ax1.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                        ax1.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                        ax2 = sea.heatmap(data_GAM, ax=ax2, cmap=get_continuous_cmap(hex_list[4]), cbar=True, vmin=0, vmax = 1)
                        ax2.set_title(r'$\gamma_%g(t = %g)$,' % (s, T))
                        ax2.xaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))
                        ax2.yaxis.set(major_locator=mtick.MultipleLocator(g / 16), major_formatter = mtick.FormatStrFormatter('%g'))

                    #Set super title
                    fig.suptitle(r'$R = %g$, $dP = %g$, $n =%d$' % (a, dP, R))
                    plt.savefig(savepath + f"CombinedImg/BioGammaConc_a_{a}_T_{T}_n_{R}.png")
                    plt.close()
                
            # End of R loop
            plt.close()
        # End of T loop

        # Combining images to create video for each a value.
        # Images are arranged to create video in increasing order of T, followed by increasing order of R.
        img_array = []
        
        """
        home_video(out_dir, out_dir, [PREFIX], input_avals, input_tvals, maxR=  input_Rmax, 
               pngformat= "CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png", pngdir= "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/",
               videoname = "guess")
        """
    


'''# Summary Description of get_FRAME_EQdata(...) (called by analyse_EQdata(...))
This function creates a multi-Indexed DataFrame with the equilibrium data for each species.
The indices are as follows:
Prefix, a, T
The columns are determined by reading the header of tab-delineated txt file MEAN_STD_Surviving_runs.txt
which is located in the directory {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/
These column names should be the same for all such files in all such directories.
The second line of the file contains the values for the columns.
All the data is stored in the DataFrame and returned.
# If a_vals and T_vals are not provided, the script reads the relevant text files and extracts the a_vals and T_vals from them.
# List of a_vals provided {in_dir}/{PREFIX}/a_vals.txt
# List of T_vals provided {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/T_vals.txt
'''
def get_FRAME_EQdata(indir, PREFIX, a_vals, T_vals, filename = "MEAN_STD_Surviving_Runs.txt"):

    a_vals, T_vals = auto_guess_avals_tvals(indir, PREFIX, a_vals, T_vals)
    if a_vals == None:
        print("Exiting..."); return;
    print(f"List of a values: {a_vals}")
    print(f"List of T values: {T_vals}")
    
    # Create a multi-indexed DataFrame
    #data = pan.MultiIndex.from_product([[PREFIX], a_vals, T_vals], names = ['Prefix', 'a', 'T'])
    data = pan.DataFrame()

    # Now populate the DataFrame with the data from the files
    for a in a_vals:
        for T in T_vals:
            try:
                # Read tab-delineated txt file
                df = pan.read_table(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/" + filename, header= 0)
                #data = pan.read_csv(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/" + filename, header= 0)
                # Use the header row as the column names and the second row as the data

                df.columns = df.columns.astype(str)
                #print(f"Data for a = {a}, T = {T}:")

                df.index = pan.MultiIndex.from_tuples([(PREFIX, a, T)], names = ['Prefix', 'a', 'T'])

                '''#print(df.info())
                #print(df.head())
                #print(df.columns)'''

                #Assign the data df to the appropriate index in the data DataFrame
                data = pan.concat([data, df], axis = 0).sort_index(axis = 0)
                #print(data.columns)
                '''
                df = df.iloc[0, :]
                print(df.info()); print(data.shape);
                print(df.columns)
                df = df.to_frame().transpose()
                # Check and set column names
                if not hasattr(df, 'columns'):
                    df.columns = df.columns.astype(str)
                # Set column names to the first row of the data
                df.columns = df.iloc[0, :]
                
                df.index = pan.MultiIndex.from_tuples([(PREFIX, a, T)], names = ['Prefix', 'a', 'T'])
                '''
            except FileNotFoundError:
                print(f"Data not found for a = {a}, T = {T}. Skipping....")
                continue
    return data




def get_FRAME_FFTPOWdata(indir, PREFIX, a_vals, T_vals, a_scale=1, 
                         include_col_labels= ["P(x; t)" , "G(x; t)", "Pr(x; t)"], filename = "FFT_POWERspectra.csv"):

    a_vals, T_vals = auto_guess_avals_tvals(indir, PREFIX, a_vals, T_vals)
    if a_vals == None:
        print(f"No Avals found for {PREFIX} Exiting..."); return;
    a_scaled_vals = [a*a_scale for a in a_vals]
    print(f"List of a values: {a_vals}")
    if (a_scale != 1):
        print(f"List of scaled a values: {a_scaled_vals}")
    print(f"List of T values: {T_vals}")
    
    # Create a multi-indexed DataFrame
    #data = pan.MultiIndex.from_product([[PREFIX], a_vals, T_vals], names = ['Prefix', 'a', 'T'])
    data = pan.DataFrame()

    # Now populate the DataFrame with the data from the files
    for a in a_vals:
        for T in T_vals:
            try:
                # Read tab-delineated txt file
                print(f"Reading FFT PS File at" + indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/FFT_PowerSpect/" + filename)
                # Report files in the directory
                print(os.listdir(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/FFT_PowerSpect/"))
                #"..\Data\Remote\Rietkerk\Reorg_Frames\1Sp\StdParam_MFT\DiC-STD\L_128_a_0.048\dP_10000\Geq_5\T_91201\FFT_PowerSpect\FFT_POWERspectra.csv"
                #"../Data/Remote/Rietkerk/Reorg_Frames/1Sp/StdParam_MFT/DiC-STD/L_128_a_0.048/dP_10000/Geq_5/T_91201/FFT_PowerSpect/FFT_POWERspectra.csv"
                df = pan.read_csv(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/FFT_PowerSpect/" + filename, header= 0, encoding='utf-8')
                print(f"Data for a = {a}, T = {T}:")
                #data = pan.read_csv(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/" + filename, header= 0)
                # Use the header row as the column names and the second row as the data

                print(df.head())
                print(df.columns)

                df.columns = df.columns.astype(str); df.columns = df.columns.str.strip()
                #print(f"Data for a = {a}, T = {T}:")
                col_labels = include_col_labels.copy(); #col_labels.remove("t");
                df = df[[col for col in df.columns if any([label in col for label in col_labels])] ]


                df.index = pan.MultiIndex.from_tuples([(PREFIX, a*a_scale, T)]*len(df), names = ['Prefix', 'a', 'T'])

                #'''#print(df.info())
                print(df.head())
                print(df.columns)#'''

                #Assign the data df to the appropriate index in the data DataFrame
                data = pan.concat([data, df], axis = 0).sort_index(axis = 0)
                #print(data.columns)

            except FileNotFoundError:
                print(f"Data not found for a = {a}, T = {T}. Skipping....")
                continue
    return data


''' # Summary Description of get_FRAME_CORRdata(...) (called by analyse_CORRdata(...))
This function creates a multi-Indexed DataFrame with the correlation data of species across different time steps.
The indices are as follows:
Prefix, a, T0, T1, maxR (Default corrfilename without HarmonicPeaks) or 
Prefix, a, T0, Rmax (if corrfilename contains HarmonicPeaks)
The columns are determined by reading the header of comma-delineated csv file {corrType}_{analType}_TD_{t-delay}.csv
If more than one matching file is used, use the file with the largest TD_{t-delay} value, where t-delay is an integer
representing the time delay between the two time steps, i.e. t0 - t1 = t-delay.
The file is located in the directory {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T0}/2DCorr/
Each such input file will have the following columns:
t0, t1, t-delay, maxR (in fixed order), followed by the correlation data for each species:
{analytype}[{var1}]_R_0, ... {analytype}[{var1}]_R_{R1}, {analytype}[{var2}]_R_0, ... {analytype}[{varN}]_R_{R1}, AVG[{analytype}[{var1}]] ... AVG[{analytype}[{varN}]]
where {var} is one of the variable names. Note that the columns are not in fixed order, i.e. var1, var2, ... varN can be in any order.
and R1 is the actual max replicate number and is generally less than maxR.

This script will read all the files in the directory and create a multi-indexed DataFrame with the data.
The data is stored in the DataFrame with the following structure:
Indices: Prefix, a, T0, T1, maxR (WHERE maxR = R1!!!)
Columns: t-delay, AVG[{var1}], ... AVG[{varN}], {analytype}[{var1}]_R_0, ... {analytype}[{var1}]_R_{R1}, {analytype}[{var2}]_R_0, ... {analytype}[{varN}]_R_{R1}

The data is stored in the DataFrame and returned, alongside var_names, which is a list of the variable names, 
as determined by reading the columns of the input files.

'''
def get_FRAME_CORRdata(indir, PREFIX, a_vals, focus_T_vals, corrsubdir="2DCorr/", corrfilename = "{corrType}_{analType}_TD_*.csv", 
                       corrType = "Auto", analType = "NCC", a_scale=1 ):
    
    a_vals, _ = auto_guess_avals_tvals(indir, PREFIX, a_vals, focus_T_vals)
    if a_vals == None:
        print("Exiting..."); return;
    a_scaled_vals = [a*a_scale for a in a_vals]
    print(f"List of a values: {a_vals}")
    if (a_scale != 1):
        print(f"List of scaled a values: {a_scaled_vals}")

    # Create a multi-indexed DataFrame
    corrdata = pan.DataFrame()
    var_names = []
    # Now populate the DataFrame with the data from the files
    for a in a_vals:
        for T0 in focus_T_vals:
            try:
                # Read the csv file corresponding to the T0 value
                # In case of multiple files, use the one with the largest t-delay value.
                #print(f"Reading {corrfilename} for a = {a}, T0 = {T0}...")

                matched_files = glob.glob(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T0}/{corrsubdir}" + 
                                          corrfilename.format(corrType=corrType, analType=analType))
                #print(f"Matched files: {matched_files}")
                if len(matched_files) == 0:
                    print(f"No files found for a = {a}, T0 = {T0}. Skipping....")
                    continue
                if len(matched_files) > 1:
                    # Sort the files by t-delay value (the last part of the filename)
                    matched_files.sort(key=lambda x: int(x.split("_")[-1].split(".")[0]), reverse=True)
                    print(f"Multiple files found for {corrType}_{analType}_TD a = {a}, T0 = {T0}. Using the file with the largest t-delay value: {os.path.basename(matched_files[0])}")
                    
                # Read the file
                df = pan.read_csv(matched_files[0], header=0)
                df.columns = df.columns.astype(str)
                df.columns = df.columns.str.strip()
                # Get the variable names from the columns of the input file.


            
                # Create the multi-index for the DataFrame
                # Comprising of Prefix, a, T0, T1, Rmax (if default corrfilename without HarmonicPeaks is used)
                # Else comprising of Prefix, a, T0, Rmax (if corrfilename contains HarmonicPeaks)
                # Note that T1 is not used in the case of HarmonicPeaks, but is used in the case of other correlation types.
                # a is scaled by a_scale, T0 and T1 are assigned right from the file.
                # Rmax is determined by reading the maximum value corresponding to _R_* in the df.columns
                if (re.search(r'HarmonicPeaks', corrfilename)):
                    Rmax = max([int(re.search(r'_R_(\d+)', col).group(1)) for col in df.columns if re.search(r'_R_', col)]) + 1
                else:
                    Rmax = max([int(col.split("_R_")[-1]) for col in df.columns if "_R_" in col]) + 1

                # Check if corrfilename contains the string "HarmonicPeaks"
                if (re.search(r'HarmonicPeaks', corrfilename)):
                    # If so, assign the [Prefix, a, T0, Rmax] to the DataFrame MultiIndex
                    df.index = pan.MultiIndex.from_tuples([(PREFIX, a*a_scale, T0, Rmax)]*len(df), names = ['Prefix', 'a', 'T0', 'maxR'])
                else:
                    # Assign the index to the DataFrame
                    df.index = pan.MultiIndex.from_tuples( [(PREFIX, a*a_scale, df.loc[i, "t0"], df.loc[i, "t1"], Rmax) 
                                                        for i in range(len(df))], names = ['Prefix', 'a', 'T0', 'T1', 'maxR'])
                
                    # Remove the columns t0, t1, maxR from the DataFrame
                    df = df.drop(columns=["t0", "t1"])

                # Concatenate the DataFrame with the existing data
                corrdata = pan.concat([corrdata, df], axis=0).sort_index(axis=0)

                # The list of unique variable names can be determined from the columns AVG[{analytype}[{var1}]] ... AVG[{analytype}[{varN}]]
                match_vars_nongamma = r'AVG\[' + f'{analType}{{' +  r'(.*?)\}\]'
                #match_vars_gamma = r'AVG\[' + f'{analType}\[GAM{{' +  r'(.*?)\]\]\]'
                file_var_names = list(set([re.search(match_vars_nongamma, col).group(1) for col in df.columns if re.search(match_vars_nongamma, col)]))
                #file_var_names.extend(list(set(["GAM[" + re.search(match_vars_gamma, col).group(1) +"]" for col in df.columns if re.search(match_vars_gamma, col)])))
                # Add unique variable names to the list of variable names
                var_names.extend(file_var_names)
                # Remove duplicates from the list of variable names
                var_names = list(set(var_names))


            
            except FileNotFoundError as e:
                print(f"{e} for a = {a}, T0 = {T0}. Skipping....")
                continue
    
    #print(df.columns)
    if (re.search(r'HarmonicPeaks', corrfilename) and not corrdata.empty):
        # If the corrfilename contains HarmonicPeaks, then do the following:
        # Remove columns (BUT NOT INDEXES) that contain the string "Rmax" "t0" or "t1"
        corrdata = corrdata.loc[:, ~corrdata.columns.str.contains("Rmax|t0|t1")]
        # Next, for each column that starts with "FFT-PEAK[*" or "FFT-PEAK-REPAVG{*" insert a new column right
        # after which is called "PEAK-TIMEPERIOD[*" or "PEAK-TIMEPERIOD-REPAVG{*" which is just the reciprocal
        # of the preceding FFT-PEAK col values (REPRESENT DIVISION BY 0 as np.nan)
        for col in corrdata.columns:
            if re.search(r'FFT-PEAK\[', col):
                # Get the peak time period by taking the reciprocal of the peak frequency
                peak_time_period = 1 / corrdata[col]
                # Create a new column name
                new_col_name = col.replace("FFT-PEAK[", "PEAK-TIMEPERIOD[")
                # Insert the new column right after the FFT-PEAK column
                corrdata.insert(corrdata.columns.get_loc(col) + 1, new_col_name, peak_time_period)
            elif re.search(r'FFT-PEAK-REPAVG\{', col):
                # Get the peak time period by taking the reciprocal of the peak frequency
                peak_time_period = 1 / corrdata[col]
                # Create a new column name
                new_col_name = col.replace("FFT-PEAK-REPAVG{", "PEAK-TIMEPERIOD-REPAVG{")
                # Insert the new column right after the FFT-PEAK column
                corrdata.insert(corrdata.columns.get_loc(col) + 1, new_col_name, peak_time_period)
            elif re.search(r'FFT-PEAK-REPVAR\{', col):
                # Get the peak time period by taking the reciprocal of the peak frequency
                peak_time_period = 1 / corrdata[col]
                # Create a new column name
                new_col_name = col.replace("FFT-PEAK-REPVAR{", "PEAK-TIMEPERIOD-REPVAR{")
                # Insert the new column right after the FFT-PEAK column
                corrdata.insert(corrdata.columns.get_loc(col) + 1, new_col_name, peak_time_period)

        # Finally replace all instances of np.inf with np.nan in the DataFrame
        corrdata.replace([np.inf, -np.inf], np.nan, inplace=True)
    return (corrdata, var_names)


""" # Summary Description of get_FRAME_CLUSTERdata(...) (called by analyse_CLUSTERdata(...))
This function creates a multi-Indexed DataFrame with the clustered frequency data for each species.
The indices are as follows:
Prefix, a, T
The columns are determined by reading the header of comma-delineated csv file CLUSTERED_FREQUENCIES.csv
which is located in the directory {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/{clustsubdir}/
The data is stored in the DataFrame and returned.
"""
def get_FRAME_CLUSTERdata(indir, PREFIX, a_vals, T_vals, a_scale =1, freqfilename = "CLUSTERED_FREQUENCIES.csv", 
                                     clustsubdir= "CLUST/{binType}_{nbins}/", binType="KMeans", nbins=2):
    
    clustsubdir = clustsubdir.format(binType=binType, nbins=nbins)
    a_vals, T_vals = auto_guess_avals_tvals(indir, PREFIX, a_vals, T_vals)
    if a_vals == None:
        print("Exiting..."); return;
    a_scaled_vals = [a*a_scale for a in a_vals]
    print(f"List of a values: {a_vals}")
    if (a_scale != 1):
        print(f"List of scaled a values: {a_scaled_vals}")
    # Create a multi-indexed DataFrame
    clustdata = pan.DataFrame()
    # Now populate the DataFrame with the data from the files
    for a in a_vals:
        for T in T_vals:
            # Read the file
            try:
                df = pan.read_csv(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/{clustsubdir.format(binType=binType, nbins=nbins)}" 
                        + freqfilename, header=0)
                df.columns = df.columns.astype(str)

                # If CDF is True, then for each column that starts with "CLUSTER_FREQ[",
                # replace the column name with "CLUSTER_FREQ_CDF[{var}]" where {var} is the variable name.
                # This replacement is done as follows:
                # CDF[X=x] = P(X <= x) = 1 - P(X > x) = 1 - CDF[X>x]

                for col in df.columns:
                    if col.startswith("CLUSTER_FREQ["):
                        # Compute the CDF: CDF[X=x] = cumulative sum from the left
                        df["CLUSTER_CDF[" + col[13:-1] + "]"] = df[col].cumsum()/ df[col].sum()
                    

                #print(f"Data for a = {a}, T = {T}:")
                df.index = pan.MultiIndex.from_tuples([(PREFIX, a*a_scale, T)]*len(df), names=['Prefix', 'a', 'T'])
                # Assign the data df to the appropriate index in the clustdata DataFrame
                clustdata = pan.concat([clustdata, df], axis=0).sort_index(axis=0)
            except FileNotFoundError:
                print(f"Cluster Frequency data not found for a = {a}, T = {T}. Skipping....")
                continue

    # Finally, get unique var_names by extracting them from the columns of the DataFrame as follows:
    # "CLUSTER_FREQ[{var}]" where {var} is the variable name

    #var_names = list(set([re.search(r'CLUSTER_FREQ\[(.*?)\]', col).group(1) 
    #        for col in clustdata.columns if re.search(r'CLUSTER_FREQ\[', col)]))
    # Next replace all 0 values in CLUSTER_FREQ[{var}] columns with np.nan
    for col in clustdata.columns:
        if re.search(r'CLUSTER_FREQ\[(.*?)\]', col):
            clustdata[col].replace(0, np.nan, inplace=True)

    return clustdata


''' # Summary Description of get_FRAME_POTdata(...) (called by analyse_POTdata(...))
This function creates a multi-Indexed DataFrame with the potential data for each species.
The indices are as follows:
Prefix, a, T, maxR
The columns are determined by reading the header of comma-delineated csv file Pot_Well.csv
which is located in the directory {in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/Pot_Well/
The column names should be the same for all such files in all such directories.
The data is stored in the DataFrame and returned.
maxR is determined by reading the first line of the text file maxR.txt located in the path:
{in_dir}/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/T_{T}/

Additionally, the local minima is provided in the csv file LOCAL_MINIMA.csv located in the same directory as Pot_Well.csv.
If read_local_minima is set to True, the local minima data is read and stored in another DataFrame and returned.
The local minima data is stored in the Multi-Indexed DataFrame with the same indices as the potential data.

Returns: (data, local_minima) if read_local_minima is set to True, else (data, None)
'''

def get_FRAME_POTdata(indir, PREFIX, a_vals, T_vals, a_scale =1, return_localminima= True, potfilename = "Pot_Well.csv", minimafilename = "LOCAL_MINIMA.csv", 
                      potsubdir = "Pot_Well/"):

    a_vals, T_vals = auto_guess_avals_tvals(indir, PREFIX, a_vals, T_vals)
    if a_vals == None:
        print("Exiting..."); return;
    print(f"List of a values: {a_vals}")
    print(f"List of T values: {T_vals}")
    a_scaled_vals = [a*a_scale for a in a_vals]
    if (a_scale != 1):
        print(f"List of scaled a values: {a_scaled_vals}")
    # Create a multi-indexed DataFrame
    data = pan.DataFrame()
    local_minima = pan.DataFrame()

    # Now populate the DataFrame with the data from the files
    for a in a_vals:
        for T in T_vals:
            # First Read maxR from maxR.txt

            try:
                with open(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt", 'r') as f:
                    maxR = int(f.readline ().strip())
            except FileNotFoundError:
                print(f"maxR.txt not found in the input directory for a = {a}, T = {T}. Skipping this a-value.")
                continue
            try:
                # Read comma-delineated csv file
                df = pan.read_csv(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/{potsubdir}" + potfilename, header= 0)
                # Use the header row as the column names and the second row as the data
                df.columns = df.columns.astype(str)
                #print(f"Data for a = {a}, T = {T}:")
                df.index = pan.MultiIndex.from_tuples([(PREFIX, a*a_scale, T , maxR)]*len(df), names = ['Prefix', 'a', 'T', 'MaxR'])

                #Assign the data df to the appropriate index in the data DataFrame
                data = pan.concat([data, df], axis = 0).sort_index(axis = 0)

            except FileNotFoundError:
                #print(f" At location: {indir + PREFIX + f'/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/{potsubdir}' + potfilename}")
                print(f"Data not found for a = {a}, T = {T}. Skipping....")
                continue

            if return_localminima:
                # Read comma-delineated csv file
                try:
                    df = pan.read_csv(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/{potsubdir}" + minimafilename, header= 0)
                    # Use the header row as the column names and the second row as the data
                    df.columns = df.columns.astype(str)
                    #print(f"Data for a = {a}, T = {T}:")

                    df.index = pan.MultiIndex.from_tuples([(PREFIX, a*a_scale, T , maxR)]*len(df), names = ['Prefix', 'a', 'T', 'MaxR'])

                    #Assign the data df to the appropriate index in the data DataFrame
                    local_minima = pan.concat([local_minima, df], axis = 0).sort_index(axis = 0)

                except FileNotFoundError:
                    print(f"Local Minima data not found for a = {a}, T = {T}. Skipping....")
                    continue

    return (data, local_minima)


def get_PRELIM_EQdata(indir, PREFIX, a_vals, tseries_vals, a_scale =1,
                      include_col_labels= ["t",  "<<P(x; t)>_x>_r" , "<<G(x; t)>_x>_r", "<<Pr(x; t)>_x>_r"], 
                      meanfilename = "MEAN_{serType}_T_{TS}.csv", serType = "TSERIES"):

    a_vals, _ = auto_guess_avals_tvals(indir, PREFIX, a_vals, [0])
    if a_vals == None:
        print("Exiting..."); return;
    #print(f"List of TS values: {tseries_vals}")
    a_scaled_vals = [a*a_scale for a in a_vals]
    if len(tseries_vals) == 0 or tseries_vals == None:
        for a in a_vals:
            a = int(a) if a.is_integer() else a
            # Generate tseries_vals if not provided.
            try:
                with open(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/TimeSeries/TimeSeries_vals.txt", 'r') as f:
                    tseries_vals.extend([float(line.strip()) for line in f])
            except FileNotFoundError:
                print(f"TimeSeries_vals.txt not found in the input directory and no TimeSeries_vals specified for a = {a}. Skipping this a-value.")
                continue
        # Sort and remove duplicates
        tseries_vals = [int(T) if T.is_integer() else T for T in tseries_vals]; tseries_vals = sorted(list(set(tseries_vals)))

                                                                                                    
    print(f"List of a values: {a_vals}")
    if (a_scale != 1):
        print(f"List of scaled a values: {a_scaled_vals}")
    print(f"List of TS values: {tseries_vals}")


    
    
    # Create a multi-indexed DataFrame
    #data = pan.MultiIndex.from_product([[PREFIX], a_vals, T_vals], names = ['Prefix', 'a', 'T'])
    data = pan.DataFrame()

    # Now populate the DataFrame with the data from the files
    for a in a_vals:
        for TS in tseries_vals:

            # Try reading Mean values from meanfilename, only include the columns that fully contain the string in include_col_labels.
            try:
                df = pan.read_csv(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/TimeSeries/" + meanfilename.format(Pre=PREFIX, g=g, dP=dP, Geq=Geq, a=a, serType= serType, TS=TS), header= 0)
                df.columns = df.columns.astype(str)
                df.columns = [col.strip() for col in df.columns]
                
                df.index = pan.MultiIndex.from_tuples([(PREFIX, a*a_scale, -1, TS)]*len(df), names = ['Prefix', 'a', 'R', 'Tmax'])
                # Only include the columns that fully contain any of the strings in include_col_labels barring "t"
                mean_col_labels = include_col_labels.copy(); mean_col_labels.remove("t")
                df = df[[col for col in df.columns if any([label in col for label in mean_col_labels])] + ["t"]]

                data = pan.concat([data, df], axis = 0).sort_index(axis = 0)
            except FileNotFoundError:
                print(f"Mean data not found for a = {a}, TS = {TS}. Skipping this a-value....")
                print(f" With missing filename: {indir + PREFIX + f'/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/TimeSeries/' + meanfilename.format(Pre=PREFIX, g=g, dP=dP, Geq=Geq, a=a, serType= serType, TS=TS)}")
                continue
            selected_files = glob.glob(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/TimeSeries/{serType}_T_{TS}_*.csv")
            # First get maxR for this a and TS value
            try:
                if serType == "MOVSERIES":
                    with open(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/TimeSeries/movmaxR_T_{TS}.txt", 'r') as f:
                        maxR = int(f.readline().strip())
                else:
                    #Default is TSERIES
                    with open(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/TimeSeries/maxR_T_{TS}.txt", 'r') as f:
                        maxR = int(f.readline().strip())
            except FileNotFoundError:
                maxR = len(selected_files)
                if serType == "MOVSERIES":
                    print(f"movmaxR_T_{TS}.txt not found in the input directory for a = {a}, TS = {TS}. Attempting to read from file selection as ... {maxR}.")
                else:
                    print(f"maxR_T_{TS}.txt not found in the input directory for a = {a}, TS = {TS}. Attempting to read from file selection as ... {maxR}.")
            for f in selected_files:
                try:
                    # Read time-series csv files.
                    df = pan.read_csv(f, header= 0)
                    #data = pan.read_csv(indir + PREFIX + f"/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/" + filename, header= 0)
                    # Use the header row as the column names and the second row as the data

                    df.columns = df.columns.astype(str)
                    # Remove white spaces from column names
                    df.columns = [col.strip() for col in df.columns]
                    #print(f"Data for a = {a}, TS = {TS}:")
                    #print(f"df columns: {df.columns}")
                    #print(f"include_col_labels: {include_col_labels}")

                    # Only include the columns specified in include_col_labels.
                    try:
                        df = df[include_col_labels]
                    except KeyError:
                        print(f"KeyError: Columns not found in the data for a = {a}, TS = {TS}. Skipping this file....")
                        continue

                    # Check if "t" is all zeros, if so, replace with "t" column values from data slice indexed by [(PREFIX, a, -1, TS)]
                    if df["t"].sum() == 0:
                        df["t"] = data.loc[(PREFIX, a*a_scale, -1, TS), "t"].values
                    Rval = int(re.search(r"R_(\d+)", f).group(1))

                    # Assign the data to the appropriate index in the DataFrame (for the length of the data)
                    df.index = pan.MultiIndex.from_tuples([(PREFIX, a*a_scale, Rval, TS)]*len(df), names = ['Prefix', 'a', 'R', 'Tmax'])

                    #df.index = pan.MultiIndex.from_tuples([(PREFIX, a*a_scale, Rval, TS)], names = ['Prefix', 'a', 'R', 'Tmax'])


                    '''#print(df.info())
                    #print(df.head())
                    #print(df.columns)'''

                    #Assign the data df to the appropriate index in the data DataFrame
                    data = pan.concat([data, df], axis = 0).sort_index(axis = 0)
                    #print(data.columns)
                    '''
                    df = df.iloc[0, :]
                    print(df.info()); print(data.shape);
                    print(df.columns)
                    df = df.to_frame().transpose()
                    # Check and set column names
                    if not hasattr(df, 'columns'):
                        df.columns = df.columns.astype(str)
                    # Set column names to the first row of the data
                    df.columns = df.iloc[0, :]
                    
                    df.index = pan.MultiIndex.from_tuples([(PREFIX, a, T)], names = ['Prefix', 'a', 'T'])
                    '''
                except FileNotFoundError:
                    print(f"Data not found for a = {a}, TS = {TS}. Skipping....")
                    continue
            # End of R loop



        # End of T loop
    # End of a loop
    return data



'''# Summary Description of analyse_timeseriesData(...)
This function analyses the timeseries data for each species and creates plots for each species for each Prefix and a value, with T as the x-axis.
NOTE: This function assumes that the data is stored in the FILENAME the following format:
R_max   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]_R_0   ...   AVG[{var}]_R_{R_max}
where {var} is one of the species names.
'''
def analyse_FRAME_timeseriesData(indir, out_dir, prefixes=[], a_vals=[], T_vals =[], filename = "MEAN_REPLICATES.txt"):

    savefilename = re.sub(r'\.txt$', '', filename)
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data = pan.DataFrame()
    sea.set_palette("husl")
    colours = [hex_list[i][-1] for i in range(len(hex_list))]

    init_a_vals = a_vals; init_T_vals = T_vals
    
    for Pre in prefixes:
        
        data = get_FRAME_EQdata(indir, Pre, init_a_vals, [], filename)
        # We are interested in the AVG columns for each species.
        # Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        # Data has columns for each species given by {var}, of the form:
        # R_max   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]_R_0   ...   AVG[{var}]_R_{R_max}
        if data.empty:
            print(f"No data found for {Pre}. Skipping....")
            continue
        print("====================================================================")
        print(f"Data for {Pre}:")
        print("====================================================================\n")
        #print(data.info())
        print(data.head())
        print(data.tail())
        print(data.columns)
        print(data.index)
        #print(data.describe())

        # Find Tmax as the maximum value of T in the data.
        Tmax = data.index.get_level_values('T').max()
        Tmax = int(Tmax) if float(Tmax).is_integer() else Tmax
        savecsvdir = out_dir + f"{Pre}/TimeSeries/"
        Path(savecsvdir).mkdir(parents=True, exist_ok=True)

        # Save the data to a csv file.
        data.to_csv(savecsvdir + f"{savefilename}_dP_{dP}_Geq_{Geq}.csv")
        #Concatenate the data to the combined data
        combined_data = pan.concat([combined_data, data], axis = 0)
        # Use the data to create plots.
        # The plots are created for each species for each Prefix and a value, with T as the x-axis.
        # For each species, create a plot of T vs all the data for that species as lines, with the mean across surviving and all replicates as darker lines.
        # The means for individual replicates are lighter lines.
        
        # Try determining number of species by counting all columns with _SURV in the name.
        try:
            num_species = len([col for col in data.columns if "_SURV" in col])
        except IndexError:
            print(f"No columns with _SURV found in the data. Not plotting data for Prefix {Pre}. Skipping....")
            num_species = None
            continue

        # Find column index of first column with _R_0 in the name. If not found, print a warning and set to 0.
        try:
            R0_index = data.columns.get_loc([col for col in data.columns if "_R_0" in col][0])
        except IndexError:
            print(f"No column with _R_0 found in the data. Not plotting individual replicate data for Prefix {Pre}.")
            R0_index = None
        print(f"R0_index = {R0_index}, num_species = {num_species}")
        # Get the list of a values and T values from the data.
        a_vals = data.index.get_level_values('a').unique().to_list()
        print(a_vals, T_vals)
        for a in a_vals:
            a = int(a) if float(a).is_integer() else a
            savepngdir = out_dir + f"{Pre}/TimeSeries/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T_{Tmax}/"
            Path(savepngdir).mkdir(parents=True, exist_ok=True)
            fig, axs = plt.subplots(1, SPB, figsize = set_figsize(SPB))
            for s in range(SPB):
                data_A = data.loc[(slice(None), a, slice(None)), :]
                # Plotting for a fixed a value, species and Prefix for all T values.
                ax = axs[s] if SPB > 1 else axs
                # RECALL, First column is R_max, then the AVG columns for surviving replicates, then AVG columns for all replicates, then AVG columns for each individual replicate.
                # Plotting mean of surviving replicates
                ax.plot(data_A.index.get_level_values('T'), data_A.iloc[:, 1+ s], label = r"$\mu_{surv}(\rho_{%g})$" %(s),
                        color = colours[s], linestyle = 'solid', linewidth = 2, alpha = 0.9)
                # Plotting mean of all replicates
                ax.plot(data_A.index.get_level_values('T'), data_A.iloc[:, 1 + num_species +  s], label = r"$\mu_{all}(\rho_{%g})$" %(s),
                        color = colours[s], linestyle = 'dashed', linewidth = 1.35, alpha = 0.75)
                # Plotting mean of individual replicates (suppress labels that represent columns with solely NaN values).
                if R0_index is not None and num_species is not None:
                    for r in range(R0_index, len(data_A.columns), num_species):
                        if not data_A.iloc[:, r].isnull().all():
                            ax.plot(data_A.index.get_level_values('T'), data_A.iloc[:, r + s],
                                    color ="grey", linestyle = 'solid', alpha = 0.45)
                ax.set_title(r"$\langle \rho_{%g}(t) \rangle_{x}$ vs $t$" %(s)  + f" at a = {a}")
                ax.set_xlabel(r'Time $(d)$')
                ax.set_ylabel(r'$\langle \rho_{%g}(t) \rangle_{x}$' % (s))
                ax.legend()
            #Set ymax for axs[3] to be 500.
            if(SPB == 3):
                axs[2].set_ylim(-5, 500)
            fig.suptitle(r"$\mu(\rho(x, t))$" + f" For {Pre} at a = {a}, dP = {dP}, Geq = {Geq}")
            plt.savefig(savepngdir + f"TimeSeries_a_{a}_dP_{dP}_Geq_{Geq}.png")
            #plt.show()
            plt.close()
        # End of a loop
        # Use home_video function to create a video of the plots.
        home_video(out_dir, out_dir, prefixes=[f"{Pre}/TimeSeries"], a_vals= a_vals, T_vals=[Tmax], maxR= 1, 
                   pngformat= "TimeSeries_a_{a}_dP_{dP}_Geq_{Geq}.png", pngdir= "{Pre}/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/", 
                   videoname = f"TimeSeries_{Pre}_Tmax_{Tmax}.mp4", video_relpath= "{Pre}/Videos/{amin}-{amax}/")
        # Note that there is only one png per a value, so the video will be a simple slideshow (and R_max = 1).
        input("Press F to Continue...")
    # End of Prefix loop
    '''
    # Finally if prefixes > 1, use the combined data to create plots of T vs. the mean of the data for each species, with Prefix as hue.
    if len(prefixes) > 1:
        # Estimate number of species by counting all columns with _SURV in the name.
        try:
            num_species = len([col for col in combined_data.columns if "_SURV" in col])
        except IndexError:
            print(f"No columns with _SURV found in the data. Skipping....")
            num_species = None
            return
        # Estimate R0_index by finding the first column with _R_0 in the name. If not found, print a warning and set to 0.  
        try:
            R0_index = combined_data.columns.get_loc([col for col in combined_data.columns if "_R_0" in col][0])
        except IndexError:
            print(f"No column with _R_0 found in the data. Not plotting individual replicate data for multiple prefixes.")
            R0_index = None
        for a in a_vals:
            fig, axs = plt.subplots(1, SPB, figsize = set_figsize(SPB))
            for s in range(SPB):
                data_A = combined_data.loc[(slice(None), a, slice(None)), :]
                ax = axs[s]
                sea.lineplot(data = data_A, x = 'T', y = data_A.iloc[:, 1+ s], hue = 'Prefix', ax = ax, alpha = 1)
                sea.lineplot(data = data_A, x = 'T', y = data_A.iloc[:,1 + num_species + s], hue = 'Prefix', ax = ax, alpha = 0.75)
                if R0_index is not None and num_species is not None:
                    for r in range(R0_index, len(data_A.columns), num_species):
                        if not data_A.iloc[:, r].isnull().all():
                            sea.lineplot(data = data_A, x = 'T', y = data_A.iloc[:, r + s], hue = 'Prefix', ax = ax, alpha = 0.35)
                ax.set_title(r"$\langle \rho_{%g}(t) \rangle_{x}$ vs $t$" %(s) + f" at a = {a}")
                ax.set_xlabel(r'Time $(d)$')
                ax.set_ylabel(r'$\langle \rho_{%g}(t) \rangle_{x}$' % (s))
                ax.legend()
            fig.suptitle(r"$\mu$ of" + f" Species For All Prefixes at a = {a}, dP = {dP}, Geq = {Geq}")
            Path(out_dir + "TimeSeries/Combined/").mkdir(parents=True, exist_ok=True)
            plt.savefig(out_dir + f"TimeSeries/Combined/TimeSeries_AllPrefixes_a_{a}_dP_{dP}_Geq_{Geq}.png")
            plt.show()
            plt.close()
        # End of a loop
        # Use home_video function to create a video of the plots.
        home_video(out_dir, out_dir, ["TimeSeries/Combined"], a_vals, [Tmax], 1, pngformat= "TimeSeries_AllPrefixes_a_{a}_dP_{dP}_Geq_{Geq}.png", videoname = f"TimeSeries_All_T_{Tmax}.mp4")
    '''

'''# Summary Description of analyse_FRAME_FFTPOWERdata(...)
This function analyses the FFT POWER data for each species and creates plots for each species for each Prefix and a value, averaged over provided T values.
It creates heatmaps of the power spectra for each species for each Prefix and a value, with frequency as the y-axis and a as the x-axis.
NOTE: This function assumes that the data is stored in the FILENAME the following format:
FREQ[{var1}] ... FREQ[{varN}]   POWER[{var1}]_MEAN_ALL ... POWER[{varN}]_MEAN_ALL POWER[{var1}]_MEAN_SURV ... POWER[{varN}]_MEAN_SURV k_R_0 ... POWER[{varN}]_R_{R_max}
where {var} is one of the species names.
It creates a Multi-Index DataFrame with Indices: Prefix, a, T, R and columns: FREQ[{var}], POWER[{var}]_MEAN_ALL, POWER[{var}]_MEAN_SURV, POWER[{var}]_R_{R_max}


'''

def analyse_FRAME_FFTPOWERdata(indir, out_dir, prefixes=[], a_vals=[], T_vals =[], Tavg_window_index = [-100000, 0], filename = "FFT_POWERspectra.csv", show_fundamental_freq=False, a_scaling=1,
                                 var_labels= [ "P(x; t)" , "G(x; t)", "Pr(x; t)"]):
    
    savefilename = re.sub(r'\.csv$', '', filename)
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data = pan.DataFrame()
    sea.set_palette("husl")
    colours = [hex_list[i][-1] for i in range(len(hex_list))]

    

    init_a_vals = a_vals; init_T_vals = T_vals
    
    if (Tavg_window_index[1] < Tavg_window_index[0]):
        print(f"Invalid Tavg_window_index = {Tavg_window_index}. Swapping values...")
        Tavg_window_index = [Tavg_window_index[1], Tavg_window_index[0]]
    

    for Pre in prefixes:

        savedir = out_dir + f"{Pre}/FFT_PS/"
        Path(savedir).mkdir(parents=True, exist_ok=True)
        
        data = get_FRAME_FFTPOWdata(indir, Pre, init_a_vals, init_T_vals, include_col_labels= var_labels, a_scale= a_scaling, filename = filename)
        # We are interested in the AVG columns for each species.
        # Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        # Data has columns for each species given by {var}, of the form:
        # FREQ[{var1}] ... FREQ[{varN}]   POWER[{var1}]_MEAN_ALL ... POWER[{varN}]_MEAN_ALL POWER[{var1}]_MEAN_SURV ... POWER[{varN}]_MEAN_SURV 
        # k_R_0 ... POWER[{varN}]_R_{R_max}
        if data.empty:
            print(f"No data found for {Pre}. Skipping....")
            continue
        print("====================================================================")
        print(f"Data for {Pre}:")
        print("====================================================================\n")
        #print(data.info())
        print(data.head())
        print(data.tail())
        print(data.columns)
        print(data.index)
        #print(data.describe())

        # Find Tmax as the maximum value of T in the data.
        Tmax = data.index.get_level_values('T').max()
        Tmax = int(Tmax) if float(Tmax).is_integer() else Tmax
        savecsvdir = out_dir + f"{Pre}/FFT_POWER/"
        Path(savecsvdir).mkdir(parents=True, exist_ok=True)

        # Save the data to a csv file.
        data.to_csv(savedir + f"{savefilename}_dP_{dP}_Geq_{Geq}.csv")
        #Concatenate the data to the combined data
        combined_data = pan.concat([combined_data, data], axis = 0)

        # Try determining number of species by counting all columns with _SURV in the name.
        try:
            num_species = len([col for col in data.columns if "_SURV" in col])
        except IndexError:
            print(f"No columns with _SURV found in the data. Not plotting data for Prefix {Pre}. Skipping....")
            num_species = None
            continue
        
        print(f"num_species = {num_species}")
        Avals_scaled_list = sorted(data.index.get_level_values('a').unique().to_list())
        print(Avals_scaled_list)
        Tvals_list = sorted(data.index.get_level_values('T').unique().to_list())
        Tmax = Tvals_list[-1]
        if(Tavg_window_index[0] < 0 and Tavg_window_index[1] <= 0):
            Twin_min = Tmax + Tavg_window_index[0]; Twin_max = Tmax + Tavg_window_index[1]
        else:
            Twin_min = Tavg_window_index[0]; Twin_max = Tavg_window_index[1]
        Tavg = data

        # Get log10 of the all the POWER[{var}]_MEAN_ALL and POWER[{var}]_MEAN_SURV columns.
        for s in range(num_species):
            Tavg["POWER[" + var_labels[s] + "]_MEAN_ALL"] = np.log10(Tavg["POWER[" + var_labels[s] + "]_MEAN_ALL"])
            Tavg["POWER[" + var_labels[s] + "]_MEAN_SURV"] = np.log10(Tavg["POWER[" + var_labels[s] + "]_MEAN_SURV"])

        # Remove rows with FREQ[{var}] > 80
        #Tavg = Tavg[Tavg["FREQ[" + var_labels[0] + "]"] <= 80]

        # Iterating over the Pre, a indices, average over the last Tavg_window_index values of T and assign to new DataFrame.
        ''' FIX THIS LATER
        Tavg = pan.DataFrame()
        for a in Avals_scaled_list:
                
                print(f"Averaging data for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}:")
                # Try averaging only if Twin_max is less than or equal to Tmax and Twin_min is greater than or equal to 0.
                if  Twin_min < 0:
                    print(f"Skipping averaging for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}.")
                    continue
                # Average over the last Tavg_window_index values of t (for average surviving/all replicates) (R = -1)
                Twin_data = data[(data.index.get_level_values('a') == a) & (data.index.get_level_values('T') >= Twin_min) 
                                 & (data.index.get_level_values('T') <= Twin_max)]
                df = Twin_data.groupby(['Prefix', 'a', 'T'])
                if df.empty:
                    print(f"No data found for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}. Skipping....")
                    continue
                # Create a new Multi-index as {Prefix, a, Twin_min, Twin_max}
                df.index = pan.MultiIndex.from_tuples([(Pre, a, Twin_min, Twin_max)]*len(df), names = ['Prefix', 'a', 'Tmin', 'Tmax'])

                Tavg = pan.concat([Tavg, df], axis = 0).sort_index(axis = 0)
            # End of a loop
        
        print(f"Data for {Pre} after averaging over {Tavg_window_index[0]} --- {Tavg_window_index[1]} values of T:")
        print(Tavg.head())
        print(Tavg.tail())
        print(Tavg.columns)
        print(Tavg.index)
        '''
        # Save the averaged data to a csv file.
        Tavg.to_csv(savedir + f"{savefilename}_dP_{dP}_Geq_{Geq}_Tavg_{Twin_min}_{Twin_max}.csv")
        # Use the data to create plots.

        # Create a heatmap of the power spectra for each species for each Prefix and a value, with frequency as the y-axis and a as the x-axis.
        # Do this using pcolor_mesh on {num_species} subplots, with a as the x-axis and frequency as the y-axis.
        # For each species, plot the mean power spectra for all replicates and surviving replicates.
        # The mean power spectra for individual replicates are not plotted.



        fig_SURV, axs_SURV = plt.subplots(1, num_species, figsize = set_figsize(num_species))
        for s in range(num_species):

            ax_SURV = axs_SURV[s] if num_species > 1 else axs_SURV
            # Plotting the mean power spectra for all replicates
            data_A = Tavg.loc[(slice(None), slice(None), T_vals[-1]), :]
            #data_A = Tavg.loc[(slice(None), slice(None), slice(None), slice(None)), :]
            # Create a mesh grid of a and FREQ[{var}] values, as X and Y respectively, with POWER[{var}]_MEAN_SURV OR  POWER[{var}]_MEAN_ALL as Z.
            X =data_A.index.get_level_values('a').unique()
            Y=data_A["FREQ[" + var_labels[s] + "]"].unique()
            Z_SURV = data_A["POWER[" + var_labels[s] + "]_MEAN_SURV"].values.reshape(len(Y), -1, order='F')
            Z_ALL = data_A["POWER[" + var_labels[s] + "]_MEAN_ALL"].values.reshape(len(data_A.index.get_level_values('a').unique()), len(data_A["FREQ[" + var_labels[s] + "]"].unique()))
            #Z_SURV = data_A["POWER[" + var_labels[s] + "]_MEAN_SURV"].values.reshape(len(data_A["FREQ[" + var_labels[s] + "]"].unique()), len(data_A.index.get_level_values('a').unique()))
            # Plot the heatmap for the mean power spectra for surviving replicates
            print(Z_SURV.shape, X.shape, Y.shape)
            print(Z_SURV)
            plot = ax_SURV.pcolormesh(X, Y, Z_SURV, cmap=sea.color_palette("viridis", as_cmap=True), shading='nearest')
            # If s = 0, set vmin and vmax for it.
            if s == 0:
                plot.set_clim(vmin = 0, vmax = 5)
            #ax_SURV.set_title(f"Power Spectra of {var_labels[s]}" )
            ax_SURV.set_xlabel(r'Rainfall R $(mm/hr)$', fontsize=14)
            ax_SURV.set_ylabel(r'Spatial Frequency ($px^{-1}$ (px $ \approx 0.1 km$))', fontsize=14)
            #axs_SURV.set_yscale('log')
            ax_SURV.set_ylim(ax_SURV.get_ylim()[::-1])
            ax_SURV.yaxis.set_major_locator(mtick.MultipleLocator(5))
            ax_SURV.tick_params(axis='x', labelsize=12) ; ax_SURV.tick_params(axis='y', labelsize=13)
            ax_SURV.set_ylim(25,0)

            if(show_fundamental_freq):
                peaks_df_cols  = [ str(k) + str(type) for k in X for type in ['_peak_freqs', '_peak_mags'] ]
                # Create a new dataframe to store the peak frequencies and magnitudes
                peaks_df = pan.DataFrame(columns=peaks_df_cols)
                # Next replace the NaN values in Z_SURV with 0
                Z_SURV = np.nan_to_num(Z_SURV, nan=0.0)  # Replace NaN values with 0
                # Use signal.find_peaks to find the first harmonic of the frequency spectrum (i.e the peak with highest magnitude)
                # Do so for each column in the data
                
                for i in range(Z_SURV.shape[1]):
                    max_signal = np.nanmax(Z_SURV[:, i]) # Get the maximum signal in the column
                    print(f"Maximum signal in column {i} for {var_labels[s]} at a = {X[i]} is {max_signal:.2f}")
                    threshold = max_signal * 0.5 # Set the threshold to 50% of the maximum signal
                    peaks, _ = find_peaks(Z_SURV[:, i], height=threshold) # Find the peaks in the column
                    print(f"Found {len(peaks)} peaks in column {i} for {var_labels[s]} at a = {X[i]} with threshold = {threshold:.2f}")
                    # Get the peak frequencies from index of Y
                    peak_freqs = Y[peaks]  # Get the frequencies corresponding to the peaks
                    peak_mags = Z_SURV[peaks, i] # Get the magnitudes of the peaks

                    maxima_data = sorted(zip(peak_freqs, peak_mags), key=lambda x: x[1], reverse=True)
                    # Store the highest peak frequency and magnitude in the peaks_df dataframe
                    if len(peaks) > 0:
                        peaks_df[str(X[i]) + '_peak_freqs'] = [maxima_data[0][0]]
                        peaks_df[str(X[i]) + '_peak_mags'] = [maxima_data[0][1]]

                print(peaks_df)
                for i in range(peaks_df.shape[0]):
                    # Get the peak frequency and magnitude for the current row
                    peak_freqs = peaks_df.iloc[i, ::2].to_numpy()
                    peak_mags = peaks_df.iloc[i, 1::2].to_numpy()
                    # Plot the peak frequencies as a dashed line
                    #ax_SURV.plot(X, peak_freqs, linestyle='--', color='tomato', linewidth=1.5, label='First Harmonic Peak' if i == 0 else "")
                    # Plot the second harmonic peaks as a dashed line
                    #ax_SURV.plot(X, peak_freqs * 2, linestyle='--', color='lightsalmon', linewidth=1.5, label='Second Harmonic Peak' if i == 0 else "")

            # Only show y-axis values between 0 and 30
            #ax_SURV.set_ylim(0, 30)

            # Set x axis ticks to display a-values at the center of the bins.
            #aval = data_A.index.get_level_values('a').unique()
            #ax_SURV.set_xticks( [aval[i] + (aval[i+1] - aval[i])/2 for i in range(len(aval) - 1)] + [aval[-1] + (aval[-1] - aval[-2])/2])

            # Add a colorbar to the plot
            cbar = fig_SURV.colorbar(plot, ax = ax_SURV)
            cbar.ax.tick_params(labelsize=13); cbar.set_label(r'$log_{10}(Power)$', fontsize=14);
        
        #fig_SURV.suptitle(r" Power Spectra vs Rainfall ($R$)" + f" For {Pre} at dP = {dP}, Geq = {Geq} averaged over T = {Twin_min} -- {Twin_max}")
        fig_SURV.suptitle(r" Power Spectra vs Rainfall ($R$)" + f" For {Pre} at dP = {dP}, Geq = {Geq} at Equilibrium")
        plt.tight_layout()
        plt.savefig(savedir + f"FFT_POWER_NOPEAK_{Pre}_dP_{dP}_Geq_{Geq}_Tavg_{T_vals[-1]}_BWIN_0.5_DEF.png", dpi=400, transparent=True)
        plt.show(); plt.close()



""" Summary Description of analyse_FRAME_CORRdata(...)
Analyse and visualize 2D auto- and cross-correlation data from simulation outputs.
This function processes correlation data files for multiple simulation prefixes, computes summary statistics,
generates various plots (line plots and violin plots) for both auto- and cross-correlation data, saves the results,
and optionally creates videos from the generated plots.
Parameters
----------
indir : str
    Input directory containing simulation data subdirectories.
out_dir : str
    Output directory where processed data and plots will be saved.
prefixes : list of str, optional
    List of simulation prefixes (subdirectory names) to process. If empty, all subdirectories in `indir` are used.
a_vals : list, optional
    List of 'a' parameter values to process. Used for filtering and plotting.
focus_T_vals : list, optional
    List of T0 (time) values to focus on for analysis and plotting.
analType : str, default "NCC"
    Type of analysis to perform (e.g., "NCC", "Morans", "AMI").
filename : str, default "{corrType}_{analType}_TD_*.csv"
    Filename pattern for correlation data files to load.
var_labels : str or list, default "deduce"
    Variable labels for plotting. If "deduce", labels are inferred from the data.
a_scaling : float, default 1
    Scaling factor applied to 'a' values when processing data.
Notes
-----
- The function expects correlation data files to be in a specific format, with MultiIndex DataFrames and columns
    for time delay, average correlation, and replicate values.
- Generates and saves line plots and violin plots for each variable, a value, and T0 value.
- Saves processed data as CSV files in the output directory.
- Optionally creates violin plots of mean correlations across replicates for each variable and a value.
Returns
-------
None
    All results are saved to disk; nothing is returned.
"""

def analyse_FRAME_CORRdata(indir, out_dir, prefixes=[], a_vals=[], focus_T_vals =[], analType = "NCC", 
                           filename = "{corrType}_{analType}_TD_*.csv", var_labels= "deduce", a_scaling = 1):
    autosavefilename = f"AUTOCORR-2D_{analType}_fTmax_{focus_T_vals[-1]}"
    crosssavefilename = f"CROSSCORR-2D_{analType}_fTmax_{focus_T_vals[-1]}"
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data_auto = pan.DataFrame(); combined_data_cross = pan.DataFrame()
    sea.set_palette("husl")

    colours = [hex_list[i][-1] for i in range(len(hex_list))]
    colours_med = [hex_list[i][-4] for i in range(len(hex_list))]
    colours_light = [hex_list[i][1] for i in range(len(hex_list))]
    init_a_vals = a_vals; init_focusT_vals = focus_T_vals

    for Pre in prefixes:

        auto_data, autovariables = get_FRAME_CORRdata(indir, Pre, init_a_vals, init_focusT_vals, 
                    corrfilename = filename, a_scale= a_scaling, analType = analType, corrType="Auto")
        '''' NOTE: autovariables is a list of the variables in the data, which is used for plotting if var_labels = "deduce".
        Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        The data is a MultiIndex DataFrame with Indices: Prefix, a, T0, T1, Rmax and columns:
        "t-delay", "AVG[{analytype}[{var}]], "{analytype}[{var}]_R_0, .... "{analytype}[{var}]_R_{R1}"
        where {var} is one of the species names, t-delay is the time delay (T0 - T1) and {analytype} is the type of analysis (e.g. NCC, Morans, AMI etc.)
        and 0 < R1 <= Rmax (i.e. the number of actual replicates reported in the data may be less than Rmax)
        '''

        cross_data, crossvariables = get_FRAME_CORRdata(indir, Pre, init_a_vals, init_focusT_vals, 
                    corrfilename = filename, a_scale= a_scaling, analType = analType, corrType="Cross")
        ''' NOTE: crossvariables is a list of the variables in the data, which is used for plotting if var_labels = "deduce".
        Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        The data is a MultiIndex DataFrame with Indices: Prefix, a, T0, T1, Rmax and columns:
        "t-delay", "AVG[{analytype}[{var}]], "{analytype}[{var}_R_0, .... "{analytype}[{var}]_R_{R1}"
        where {var} is one of the species names, t-delay is the time delay (T0 - T1) and {analytype} is the type of analysis (e.g. NCC, Morans, AMI etc.)'''

        if auto_data.empty and cross_data.empty:
            print(f"No correlation data found for {Pre}. Skipping....")
            continue
        if auto_data.empty:
            print(f"WARNING: No auto correlation data found for {Pre}.")
        if cross_data.empty:
            print(f"WARNING: No cross correlation data found for {Pre}.")
        

        print("====================================================================")
        print(f"Auto-data for {Pre}:")
        print("====================================================================\n")

        #print(auto_data.info())
        print(auto_data.head())
        print(auto_data.tail())
        print(auto_data.columns)
        print(auto_data.index)

        #print(auto_data.describe())
        print("\n====================================================================")
        print(f"Cross-data for {Pre}:")
        print("====================================================================\n")
        #print(cross_data.info())
        print(cross_data.head())
        print(cross_data.tail())
        print(cross_data.columns)
        print(cross_data.index)
        #print(cross_data.describe())

        # Also report found auto/cross variables.
        print("-------------------------------------------------------------------------\n")
        print(f"Auto variables deduced: {autovariables}");
        print(f"Cross variables deduced: {crossvariables}")
        print("-------------------------------------------------------------------------\n")

        # Save the data to a csv file.
        savecsvdir = out_dir + f"{Pre}/2DCorrelation/" 
        Path(savecsvdir).mkdir(parents=True, exist_ok=True)
        auto_data.to_csv(savecsvdir + f"{autosavefilename}_dP_{dP}_Geq_{Geq}.csv")
        cross_data.to_csv(savecsvdir + f"{crosssavefilename}_dP_{dP}_Geq_{Geq}.csv")

        # Get max t-delay value from the data.
        max_tdelay_auto = auto_data["t-delay"].max(); max_tdelay_auto = round(max_tdelay_auto) #if float(max_tdelay_auto).is_integer() else max_tdelay_auto
        max_tdelay_cross = cross_data["t-delay"].max(); max_tdelay_cross = round(max_tdelay_cross) #if float(max_tdelay_cross).is_integer() else max_tdelay_cross

        # For custom plots
        if SPB == 3 and len(autovariables) >= 5:
            autovariables = ['P(x; t)', 'G(x; t)', 'Pr(x; t)', 'GAM[G(x; t)]', 'GAM[Pr(x; t)]']
        if SPB == 3 and len(crossvariables) >= 20:
            crossvariables = ['P(x; t);G(x; t)', 'G(x; t);P(x; t)', 'G(x; t);Pr(x; t)', 'Pr(x; t);G(x; t)',
                            'GAM[G(x; t)];G(x; t)', 'G(x; t);GAM[G(x; t)]', 'GAM[G(x; t)];Pr(x; t)', 'Pr(x; t);GAM[G(x; t)]',
                            'GAM[Pr(x; t)];Pr(x; t)', 'Pr(x; t);GAM[Pr(x; t)]', 'GAM[Pr(x; t)];G(x; t)', 'G(x; t);GAM[Pr(x; t)]',
                            'GAM[G(x; t)];P(x; t)', 'P(x; t);GAM[G(x; t)]', 'GAM[G(x; t)];GAM[Pr(x; t)]', 'GAM[Pr(x; t)];GAM[G(x; t)]',
                            'P(x; t);Pr(x; t)', 'Pr(x; t);P(x; t)', 'GAM[Pr(x; t)];P(x; t)', 'P(x; t);GAM[Pr(x; t)]']
    
        # Making auto-correlation plots for each a value and T0 value first.
        if not auto_data.empty:
            # Working on auto_data 
            a_scaled_vals = auto_data.index.get_level_values('a').unique().to_list()
            focus_T_vals = auto_data.index.get_level_values('T0').unique().to_list()


            for T0 in focus_T_vals:
                T0 = int(T0) if float(T0).is_integer() else T0

                violin_data_all = defaultdict(list)  # Keys: (variable, a), (variable, '$\langle a\rangle$')
                # Stores violin plot data for each autovariable and T0 value.
                avg_viol_data = defaultdict(lambda: defaultdict(dict))  # Structure: {variable: {a: DataFrame}}

                saveviolpngdir = out_dir + f"{Pre}/2DCorrelation/BoxViolin/L_{g}/dP_{dP}/Geq_{Geq}/T0_{T0}/"
                # Create the directory if it does not exist.
                Path(saveviolpngdir).mkdir(parents=True, exist_ok=True);

                for a in a_scaled_vals:
                    a = int(a) if float(a).is_integer() else a
                    savepngdir_auto = out_dir + f"{Pre}/2DCorrelation/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T0_{T0}/"
                    
                    # Create the directory if it does not exist.
                    Path(savepngdir_auto).mkdir(parents=True, exist_ok=True);
                    try:
                        auto_data_A = auto_data.loc[(slice(None), a, T0, slice(None), slice(None)), :]
                        #cross_data_A = cross_data.loc[(slice(None), a, T0, slice(None), slice(None)), :]
                    except KeyError:
                        print(f"No data found for {Pre} at a = {a}, T0 = {T0}. Skipping....")
                        continue

                    # Plotting for a fixed a value, species and Prefix for all R values.
                    # Create subplots of shape (p, p) where p is round(sqrt(len(autovariables)))
                    p_auto = len(autovariables)// SPB + 1; q_auto = SPB
                    fig_auto, axs_auto = plt.subplots(p_auto, q_auto, figsize = set_figsizeGrid(p_auto, q_auto), sharex = True, sharey = True)

                    # Flatten the axs array for easier indexing.
                    axs_auto = axs_auto.flatten(); #axs_cross = axs_cross.flatten()
                    # Loop over the variables and plot the data.
                    for s in range(len(autovariables)):
                        # Plotting the auto data
                        
                        # First plotting the AVG data vs t-delay
                        ax_auto = axs_auto[s] if len(autovariables) > 1 else axs_auto
                        # For each species, plot the first and third columns
                        print(f"Plotting MEAN {analType}[{autovariables[s]}] for {Pre} at a = {a}, T0 = {T0}.")
                        ax_auto.plot(auto_data_A["t-delay"], auto_data_A["AVG[" + analType + "{" + autovariables[s] + "}]"], label = r"$\mu_{all}(" + autovariables[s] + ")$",
                                color = colours[s % len(colours)], linestyle = 'solid', linewidth = 2, alpha = 0.85)
                        # Plotting the mean of all replicates

                        # Get Rmax from data_A
                        Rmax = int(auto_data_A.index.get_level_values('maxR').max())
                        # Plotting the mean of surviving replicates (if the column exists)
                        for R in range(0, Rmax):
                            # Check if the column exists and is not empty before plotting.
                            if f"{analType}{{{autovariables[s]}}}_R_{R}" in auto_data_A.columns and not auto_data_A[f"{analType}{{{autovariables[s]}}}_R_{R}"].isnull().all():
                                #print(f"Plotting R {analType}[{autovariables[s]}]_R_{R} for {Pre} at a = {a}, T0 = {T0}.")
                                ax_auto.plot(auto_data_A["t-delay"], auto_data_A[f"{analType}{{{autovariables[s]}}}_R_{R}"],
                                        color = "grey", linestyle = 'solid', alpha = 0.45)
                        # Set the title and labels for the plot.
                        ax_auto.set_title(r"$\langle$ " + analType + "[" + autovariables[s] + r"] $\rangle_{x}$ vs $\Delta(T)$" + f" at a = {a}, T0 = {T0}")
                        ax_auto.set_xlabel(r'Time Delay $\Delta(T)$ $(hr)$')
                        ax_auto.set_ylabel(r'$\langle$' + analType + "[" + autovariables[s] + r'] $\rangle_{x}$' )
                        ax_auto.legend()

                        # Next, generate violin plot data for the current variable and a value.
                        avg_viol_data[autovariables[s]][str(a)] = auto_data_A[["AVG[" + analType + "{" + autovariables[s] + "}]", "t-delay"]].dropna()
                        violin_data_all[(autovariables[s], str(a))].extend(avg_viol_data[autovariables[s]][str(a)]["AVG[" + analType + "{" + autovariables[s] + "}]"].values.tolist())

                        # Remove existing Multindex index for avg_viol_data[autovariables[s]][str(a)] and set it to a single index
                        # given by the t-delay values, and then drop the t-delay column.
                        avg_viol_data[autovariables[s]][str(a)].reset_index(drop=True, inplace=True)
                        # Set the index to the t-delay values.
                        avg_viol_data[autovariables[s]][str(a)].set_index("t-delay", inplace=True)
                        # Drop the t-delay column.
                        #avg_viol_data[autovariables[s]][str(a)].drop(columns = ["t-delay"], inplace=True)

                        #print(avg_viol_data[autovariables[s]][str(a)].head())

                    # End of s loop
                    fig_auto.suptitle(r"$\langle$ " + analType + r" $\rangle_{x}$  vs $\Delta(T)$" + f" For {Pre} at a = {a}, T0 = {T0}, dP = {dP}, Geq = {Geq}")
                    plt.tight_layout()
                    plt.savefig(savepngdir_auto + f"{autosavefilename}_TDelay_{max_tdelay_auto}.png")
                    plt.close()
                    #plt.show(); plt.close()
                # End of a loop

                # Finally get the (variable, 'avg') data for the violin plot, storing it in  violin_data_all.
                if avg_viol_data:
                    for var in autovariables:
                        # Get the average of the data for all a values.
                        avg_viol_data[var]['avg'] = pan.concat([avg_viol_data[var][str(a)] for a in a_scaled_vals], axis=1).mean(axis=1)
                        # Store the data in violin_data_all.
                        violin_data_all[(var, r'$\langle a\rangle$')].extend(avg_viol_data[var]['avg'].values.tolist())

                    # Generate violin_df and plot the violin plot for each variable.
                    violin_df_all = pan.DataFrame([ {"variable": var, "a": a, "value": value} 
                                                for (var, a), valslist in violin_data_all.items() for value in valslist])
                    
                    # Create a violin plot for each variable.
                    p_violin = len(autovariables)// SPB + 1; q_violin = SPB
                    fig_violin, axs_violin = plt.subplots(p_violin, q_violin, figsize = set_figsizeGrid(p_violin, q_violin), sharex = True, sharey = True)
                    # Ensures that <a> is always the first variable in the subplot followed by the rest of the a vals.
                    cat_var_plot_order = [r'$\langle a\rangle$'] + sorted(x for x in violin_df_all['a'].unique() if x != r'$\langle a\rangle$')
                    # Flatten the axs array for easier indexing.
                    axs_violin = axs_violin.flatten()
                    # Loop over the variables and plot the data.
                    for s in range(len(autovariables)):
                        # Plotting the auto data
                        ax_violin = axs_violin[s] if len(autovariables) > 1 else axs_violin
                        
                        print(f"Plotting AUTO VIOLIN {analType}[{autovariables[s]}] for {Pre} at T0 = {T0}.")
                        # First plot mean.
                        if not violin_df_all[(violin_df_all['variable'] == autovariables[s]) & (violin_df_all['a'] == r'$\langle a\rangle$')].empty:
                            sea.violinplot(data=violin_df_all[(violin_df_all['variable'] == autovariables[s]) & 
                                    (violin_df_all['a'] == r'$\langle a\rangle$')], x='a', y='value', ax=ax_violin, cut=0, bw_adjust=0.5, bw_method='scott',
                                    color=colours[s % len(colours)], inner='box', linewidth=1.5, order=cat_var_plot_order, alpha=0.85, fill=True, density_norm="area")
                        # Next plot violin plots for each a value.
                        sea.violinplot(data=violin_df_all[(violin_df_all['variable'] == autovariables[s]) 
                                    & (violin_df_all['a'] != r'$\langle a\rangle$')], x='a', y='value', ax=ax_violin,  order=cat_var_plot_order, 
                                    color=colours_med[s % len(colours_med)], inner='box', linewidth=1.1, alpha=0.85, cut=0, bw_adjust=0.5, bw_method='scott', density_norm="area")
                        # Add a vertical grey line seperating the mean violin plot (i.e. at a = r'$\langle a\rangle$') from the rest.
                        ax_violin.axvline(x=0.5, color='grey', linestyle='--', linewidth=1.5, alpha=0.75)
                        # Set the title and labels for the plot.
                        ax_violin.set_title(f" AVG {analType}[{autovariables[s]}] vs a")
                        ax_violin.set_xlabel(r'R $(mm/d)$')
                        ax_violin.set_ylabel(r'$\langle$' + analType + "[" + autovariables[s] + r'] $\rangle_{R}$' )
                        
                    # End of s loop
                    fig_violin.suptitle(r"$\langle$ " + analType + r" $\rangle_{R}$" + f" For {Pre} at T0 = {T0}, dP = {dP}, Geq = {Geq}")
                    plt.tight_layout()
                    plt.savefig(saveviolpngdir + f"{autosavefilename}_TDelay_{max_tdelay_auto}_Violin.png")
                    plt.show()
                    plt.close()
            # End of T0 loop          
        # Now plotting the cross data
            
        if not cross_data.empty:
            a_scaled_vals = auto_data.index.get_level_values('a').unique().to_list()
            focus_T_vals = auto_data.index.get_level_values('T0').unique().to_list()
            for T0 in focus_T_vals:
                T0 = int(T0) if float(T0).is_integer() else T0

                violin_data_all = defaultdict(list)  # Keys: (variable, a), (variable, '$\langle a\rangle$')
                # Stores violin plot data for each crossvariable and T0 value.
                avg_viol_data_cross = defaultdict(lambda: defaultdict(dict))  # Structure: {variable: {a: DataFrame}}

                saveviolpngdir_cross = out_dir + f"{Pre}/2DCorrelation/BoxViolin/L_{g}/dP_{dP}/Geq_{Geq}/T0_{T0}/"
                # Create the directory if it does not exist.
                Path(saveviolpngdir_cross).mkdir(parents=True, exist_ok=True);
                for a in a_scaled_vals:
                    a = int(a) if float(a).is_integer() else a
                    savepngdir_cross = out_dir + f"{Pre}/2DCorrelation/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T0_{T0}/"
                    Path(savepngdir_cross).mkdir(parents=True, exist_ok=True)
                    try:
                        cross_data_A = cross_data.loc[(slice(None), a, T0, slice(None), slice(None)), :]
                    except KeyError:
                        print(f"No CROSS data found for {Pre} at a = {a}, T0 = {T0}. Skipping....")
                        continue
                    # Plotting for a fixed a value, species and Prefix for all R values.
                    # Create subplots of shape (p, p) where p is round(sqrt(len(crossvariables)))
                    p_cross = len(crossvariables)// SPB + 1; q_cross = SPB
                    #p_cross = len(crossvariables)//4 + 1; q_cross = 4
                    fig_cross, axs_cross = plt.subplots(p_cross, q_cross, figsize = set_figsizeGrid(p_cross, q_cross), sharex = True, sharey = True)
                    #fig_cross, axs_cross = plt.subplots(p_cross, p_cross, figsize = set_figsizeSquare(len(crossvariables)), sharex = True, sharey = True)
                    # Flatten the axs array for easier indexing.
                    axs_cross = axs_cross.flatten()
                    # Loop over the variables and plot the data.
                    for s in range(len(crossvariables)):
                        # Plotting the cross data
                        
                        # First plotting the AVG data vs t-delay
                        ax_cross = axs_cross[s] if len(crossvariables) > 1 else axs_cross
                        # For each species, plot the first and third columns
                        print(f"Plotting MEAN {analType}{{{crossvariables[s]}}}] for {Pre} at a = {a}, T0 = {T0}.")
                        ax_cross.plot(cross_data_A["t-delay"], cross_data_A["AVG[" + analType + "{" + crossvariables[s] + "}]"], 
                                label = r"$\mu_{all}(" + crossvariables[s] + ")$", color = colours[s % len(colours)],
                                linestyle = 'solid', linewidth = 2, alpha = 0.85)
                        # Plotting the mean of all replicates

                        # Get Rmax from data_A
                        Rmax = int(cross_data_A.index.get_level_values('maxR').max())
                        # Plotting the mean of surviving replicates (if the column exists)
                        for R in range(0, Rmax):
                            # Check if the column exists and is not empty before plotting.
                            if f"{analType}{{{crossvariables[s]}}}_R_{R}" in cross_data_A.columns and not cross_data_A[f"{analType}{{{crossvariables[s]}}}_R_{R}"].isnull().all():
                                #print(f"Plotting R {analType}[{crossvariables[s]}]_R_{R} for {Pre} at a = {a}, T0 = {T0}.")
                                ax_cross.plot(cross_data_A["t-delay"], cross_data_A[f"{analType}{{{crossvariables[s]}}}_R_{R}"],
                                        color = "grey", linestyle = 'solid', alpha = 0.45)
                        # Set the title and labels for the plot.
                        ax_cross.set_title(r"$\langle$ " + analType + "[" + crossvariables[s] + r"] $\rangle_{x}$ vs $\Delta(T)$" + f" at a = {a}, T0 = {T0}")
                        ax_cross.set_xlabel(r'Time Delay $\Delta(T)$ $(hr)$')
                        ax_cross.set_ylabel(r'$\langle$' + analType + "[" + crossvariables[s] + r'] $\rangle_{x}$' )
                        ax_cross.legend()

                        avg_viol_data_cross[crossvariables[s]][str(a)] = cross_data_A[["AVG[" + analType + "{" + crossvariables[s] + "}]", "t-delay"]].dropna()
                        violin_data_all[(crossvariables[s], str(a))].extend(avg_viol_data_cross[crossvariables[s]][str(a)]["AVG[" 
                                                            + analType + "{" + crossvariables[s] + "}]"].values.tolist())

                        # Remove existing Multindex index for avg_viol_data[crossvariables[s]][str(a)] and set it to a single index
                        # given by the t-delay values, and then drop the t-delay column.
                        avg_viol_data_cross[crossvariables[s]][str(a)].reset_index(drop=True, inplace=True)
                        # Set the index to the t-delay values.
                        avg_viol_data_cross[crossvariables[s]][str(a)].set_index("t-delay", inplace=True)
                    # End of s loop
                    
                    fig_cross.suptitle(r"$\langle$ " + analType + r" $\rangle_{x}$" + f" For {Pre} at a = {a}, T0 = {T0}, dP = {dP}, Geq = {Geq}")
                    plt.tight_layout()
                    plt.savefig(savepngdir_cross + f"{crosssavefilename}_TDelay_{max_tdelay_cross}.png")
                    plt.close()
                    #plt.show(); plt.close()
                # End of a loop
            
                # Finally get the (variable, 'avg') data for the violin plot, storing it in  violin_data_all.
                if avg_viol_data_cross:
                    for var in crossvariables:
                        # Get the average of the data for all a values.
                        avg_viol_data_cross[var]['avg'] = pan.concat([avg_viol_data_cross[var][str(a)] for a in a_scaled_vals], axis=1).mean(axis=1)
                        # Store the data in violin_data_all.
                        violin_data_all[(var, r'$\langle a\rangle$')].extend(avg_viol_data_cross[var]['avg'].values.tolist())

                    # Generate violin_df and plot the violin plot for each variable.
                    violin_df_all = pan.DataFrame([ {"variable": var, "a": a, "value": value} 
                                                for (var, a), valslist in violin_data_all.items() for value in valslist])
                    
                    # Create a violin plot for each variable.
                    #p_violin = len(crossvariables)// SPB + 1; q_violin = SPB
                    p_violin = len(crossvariables)// 4; q_violin = 4
                    fig_violin, axs_violin = plt.subplots(p_violin, q_violin, figsize = set_figsizeGrid(p_violin, q_violin), sharex = True, sharey = True)
                    # Ensures that <a> is always the first variable in the subplot followed by the rest of the a vals.
                    cat_var_plot_order = [r'$\langle a\rangle$'] + sorted(x for x in violin_df_all['a'].unique() if x != r'$\langle a\rangle$')
                    # Flatten the axs array for easier indexing.
                    axs_violin = axs_violin.flatten()
                    # Loop over the variables and plot the data.
                    for s in range(len(crossvariables)):
                        # Plotting the cross data
                        ax_violin = axs_violin[s] if len(crossvariables) > 1 else axs_violin
                        
                        print(f"Plotting CROSS VIOLIN {analType}[{crossvariables[s]}] for {Pre} at T0 = {T0}.")
                        # First plot mean.
                        if not violin_df_all[(violin_df_all['variable'] == crossvariables[s]) & (violin_df_all['a'] == r'$\langle a\rangle$')].empty:
                            sea.violinplot(data=violin_df_all[(violin_df_all['variable'] == crossvariables[s]) & 
                                    (violin_df_all['a'] == r'$\langle a\rangle$')], x='a', y='value', ax=ax_violin, cut=0, bw_adjust=0.5, bw_method='scott',
                                    color=colours[s % len(colours)], inner='box', linewidth=1.5, order=cat_var_plot_order, alpha=0.85, fill=True, density_norm="area")
                        # Next plot violin plots for each a value.
                        sea.violinplot(data=violin_df_all[(violin_df_all['variable'] == crossvariables[s]) 
                                    & (violin_df_all['a'] != r'$\langle a\rangle$')], x='a', y='value', ax=ax_violin,  order=cat_var_plot_order, 
                                    color=colours_med[s % len(colours_med)], inner='box', linewidth=1.1, alpha=0.85, cut=0, bw_adjust=0.5, bw_method='scott', density_norm="area")
                        # Add a vertical grey line seperating the mean violin plot (i.e. at a = r'$\langle a\rangle$') from the rest.
                        ax_violin.axvline(x=0.5, color='grey', linestyle='--', linewidth=1.5, alpha=0.75)
                        # Set the title and labels for the plot.
                        ax_violin.set_title(f" AVG {analType}[{crossvariables[s]}] vs a")
                        ax_violin.set_xlabel(r'R $(mm/d)$')
                        ax_violin.set_ylabel(r'$\langle$' + analType + "[" + crossvariables[s] + r'] $\rangle_{R}$' )
                        
                    # End of s loop
                    fig_violin.suptitle(r"$\langle$ " + analType + r" $\rangle_{R}$" + f" For {Pre} at T0 = {T0}, dP = {dP}, Geq = {Geq}")
                    plt.tight_layout()
                    plt.savefig(saveviolpngdir_cross + f"{crosssavefilename}_TDelay_{max_tdelay_cross}_Violin.png")
                    plt.show()
                    plt.close()
            # End of T0 loop


        # Use home_video function to create a video of the plots.
        home_video(out_dir, out_dir, prefixes=[f"{Pre}/2DCorrelation"], a_vals= a_scaled_vals, T_vals=focus_T_vals, maxR= 1, 
                   pngformat= "%s_TDelay_%g.png" %(autosavefilename, max_tdelay_auto), 
                   pngdir= "{Pre}/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T0_{T}/", 
                   videoname = f"{autosavefilename}_{Pre}_TDmax_{max_tdelay_auto}.mp4", video_relpath= "{Pre}/Videos/{amin}-{amax}/T0_{Tmax}/")
        # Note that there is only one png per a value, so the video will be a simple slideshow (and R_max = 1).

        # Use home_video function to create a video of the cross plots.
        home_video(out_dir, out_dir, prefixes=[f"{Pre}/2DCorrelation"], a_vals= a_scaled_vals, T_vals=focus_T_vals, maxR= 1,
                   pngformat= "%s_TDelay_%g.png" %(crosssavefilename, max_tdelay_cross),
                   pngdir= "{Pre}/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T0_{T}/",
                   videoname = f"{crosssavefilename}_{Pre}_TDmax_{max_tdelay_cross}.mp4", video_relpath= "{Pre}/Videos/{amin}-{amax}/T0_{Tmax}/")
        input("Press F to Continue...")
        

def analyse_FRAME_FFTCORRdata(indir, out_dir, prefixes=[], a_vals=[], focus_T_vals =[], analType = "NCC", corrsubdir="2DCorr/FFT/",
                           filename = "HarmonicPeaks_{corrType}_{analType}_TD_*.csv", var_labels= "deduce", harm_peaks_index=[0,1,2], a_scaling = 1):
    
    autosavefilename = f"FFTPIQUES_AUTOCORR-2D_{analType}_fTmax_{focus_T_vals[-1]}"
    crosssavefilename = f"FFTPIQUES_CROSSCORR-2D_{analType}_fTmax_{focus_T_vals[-1]}"

    # Check if harm_peaks_index is a list or array with integer values.
    if not ( isinstance(harm_peaks_index, (list, np.ndarray)) 
            and all(isinstance(i, int) for i in harm_peaks_index)):
        print(f"harm_peaks_index is not a valid list or array of integers. Using default harm_peaks_index = [0, 1, 2].")
        harm_peaks_index = [0, 1, 2]  # Default values if not provided or invalid.


    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data_auto = pan.DataFrame(); combined_data_cross = pan.DataFrame()
    sea.set_palette("husl")

    colours = [hex_list[i][-1] for i in range(len(hex_list))]
    colours_med = [hex_list[i][-4] for i in range(len(hex_list))]
    colours_light = [hex_list[i][1] for i in range(len(hex_list))]
    init_a_vals = a_vals; init_focusT_vals = focus_T_vals

    for Pre in prefixes:

        auto_FFTPeakdata, auto_FFTPeakvariables = get_FRAME_CORRdata(indir, Pre, init_a_vals, init_focusT_vals, 
                    corrsubdir=corrsubdir, corrfilename = filename, a_scale= a_scaling, analType = analType, corrType="Auto")
        '''' NOTE: autovariables is a list of the variables in the data, which is used for plotting if var_labels = "deduce".
        Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        The data is a MultiIndex DataFrame with Indices: Prefix, a, T0, Rmax and columns:
        "FFT-PEAK[AVG[{analytype}[{var}]]", "PEAK-TIMEPERIOD[AVG[{analytype}[{var}]]", "FFT-PEAK-VALS[AVG[{analytype}[{var}]]",
        "FFT-PEAK[{analytype}[{var}]_R_0], .... "FFT-PEAK-VALS[{analytype}[{var}]_R_{R1}]"
        where {var} is one of the species names, and {analytype} is the type of analysis (e.g. NCC, Morans, AMI etc.)
        and 0 < R1 <= Rmax (i.e. the number of actual replicates reported in the data may be less than Rmax)
        '''

        cross_FFTPeakdata, cross_FFTPeakvariables = get_FRAME_CORRdata(indir, Pre, init_a_vals, init_focusT_vals, 
                    corrsubdir=corrsubdir, corrfilename = filename, a_scale= a_scaling, analType = analType, corrType="Cross")
        ''' NOTE: crossvariables is a list of the variables in the data, which is used for plotting if var_labels = "deduce".
        Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        The data is a MultiIndex DataFrame with Indices: Prefix, a, T0, Rmax and columns:
        "FFT-PEAK[AVG[{analytype}[{var}]]", "PEAK-TIMEPERIOD[AVG[{analytype}[{var}]]", "FFT-PEAK-VALS[AVG[{analytype}[{var}]]",
        "FFT-PEAK[{analytype}[{var}]_R_0], .... "FFT-PEAK-VALS[{analytype}[{var}]_R_{R1}]"
        where {var} is one of the species names, t-delay is the time delay (T0 - T1) and {analytype} is the type of analysis (e.g. NCC, Morans, AMI etc.)'''

        if auto_FFTPeakdata.empty and cross_FFTPeakdata.empty:
            print(f"No data found for {Pre}. Skipping....")
            continue
        print("====================================================================")
        print(f"Auto-data for {Pre}:")
        print("====================================================================\n")
        print(auto_FFTPeakdata.head())
        print(auto_FFTPeakdata.tail())
        print(auto_FFTPeakdata.columns)
        print(auto_FFTPeakdata.index)
        print("====================================================================")
        print(f"Cross-data for {Pre}:")
        print("====================================================================\n")
        print(cross_FFTPeakdata.head())
        print(cross_FFTPeakdata.tail())
        print(cross_FFTPeakdata.columns)
        print(cross_FFTPeakdata.index)

        # Report found variables in the data.
        print("-----------------------------------------------------------------------")
        print(f"Deduced {len(auto_FFTPeakvariables)} auto variables : {auto_FFTPeakvariables}")
        print(f"Deduced {len(cross_FFTPeakvariables)} cross variables : {cross_FFTPeakvariables}")
        print("-----------------------------------------------------------------------\n")

        # Save the data to a csv file.
        savecsvdir = out_dir + f"{Pre}/2DCorrelation/" 
        Path(savecsvdir).mkdir(parents=True, exist_ok=True)
        auto_FFTPeakdata.to_csv(savecsvdir + f"{autosavefilename}_dP_{dP}_Geq_{Geq}.csv")
        cross_FFTPeakdata.to_csv(savecsvdir + f"{crosssavefilename}_dP_{dP}_Geq_{Geq}.csv")

        # For custom plots
        if SPB == 3 and len(auto_FFTPeakvariables) >= 5:
            auto_FFTPeakvariables = ['P(x; t)', 'G(x; t)', 'Pr(x; t)', 'GAM[G(x; t)]', 'GAM[Pr(x; t)]']
        if SPB == 3 and len(cross_FFTPeakvariables) >= 20:
            cross_FFTPeakvariables = ['P(x; t);G(x; t)', 'G(x; t);P(x; t)', 'G(x; t);Pr(x; t)', 'Pr(x; t);G(x; t)',
                            'GAM[G(x; t)];G(x; t)', 'G(x; t);GAM[G(x; t)]', 'GAM[G(x; t)];Pr(x; t)', 'Pr(x; t);GAM[G(x; t)]',
                            'GAM[Pr(x; t)];Pr(x; t)', 'Pr(x; t);GAM[Pr(x; t)]', 'GAM[Pr(x; t)];G(x; t)', 'G(x; t);GAM[Pr(x; t)]',
                            'GAM[G(x; t)];P(x; t)', 'P(x; t);GAM[G(x; t)]', 'GAM[G(x; t)];GAM[Pr(x; t)]', 'GAM[Pr(x; t)];GAM[G(x; t)]',
                            'P(x; t);Pr(x; t)', 'Pr(x; t);P(x; t)', 'GAM[Pr(x; t)];P(x; t)', 'P(x; t);GAM[Pr(x; t)]'] 
            

        # Making scatter plots for the data, in the form of a grid of subplots (reprsenting the variables).
        # This scatter plot represents FFT Time Periods vs a value for each variable (with individual replicates in a light grey shade).
        if not auto_FFTPeakdata.empty:
            # Working on auto_data 
            a_scaled_vals = auto_FFTPeakdata.index.get_level_values('a').unique().to_list()
            focus_T_vals = auto_FFTPeakdata.index.get_level_values('T0').unique().to_list()


            for T0 in focus_T_vals:
                T0 = int(T0) if float(T0).is_integer() else T0

                # If harm_peaks_index is not None or an empty list/array, then iterate over the indices "i"
                # and extract the "ith" ( i % 10) rows from the auto_FFTPeakdata DataFrame.
                    
                for harm_indx in harm_peaks_index:
                    # Extract the ith (i % 10) rows from the auto_FFTPeakdata DataFrame.
                    try:
                        auto_data_T0 = auto_FFTPeakdata.loc[(slice(None), slice(None), T0, slice(None)), :]
                    except KeyError:
                        print(f"No AUTO data found for {Pre} at T0 = {T0}. Skipping....")
                        continue
                    # Extract the ith (i % 10) rows from the auto_FFTPeakdata DataFrame.
                    # If harm_indx is -1, then we simply plot all the data for that T0 value.
                    if harm_indx == -1:
                        print(f"Plotting ALL AUTO PEAK data (Superimposed) for {Pre} at T0 = {T0}....")
                    if harm_indx != -1:
                        try:
                            auto_data_T0 = auto_data_T0.iloc[[(i % 10) == harm_indx for i in range(len(auto_data_T0))], :]  # Note the double brackets
                        except IndexError:
                            print(f"No AUTO data found for {Pre} at T0 = {T0} for harm_peaks_index = {harm_indx}. Skipping....")
                            continue
                
                    
                    # Create save directory for the plots.
                    savepngdir_auto = out_dir + f"{Pre}/2DCorrelation/BoxViolin/L_{g}/dP_{dP}/Geq_{Geq}/T0_{T0}/FFTTime/"
                    Path(savepngdir_auto).mkdir(parents=True, exist_ok=True)
                    # Create subplots of shape (p, p) where p is round(sqrt(len(auto_FFTPeakvariables)))
                    #p_auto = len(auto_FFTPeakvariables)// SPB + 1; q_auto = SPB
                    p_auto = len(auto_FFTPeakvariables)//4 + 1; q_auto = 4
                    fig_auto, axs_auto = plt.subplots(p_auto, q_auto, figsize = set_figsizeGrid(p_auto, q_auto), sharex = True, sharey = True)
                    print(f"Auto_data_T0 index: {auto_data_T0.index}")
                    # Flatten the axs array for easier indexing.
                    axs_auto = axs_auto.flatten()
                    # Loop over the variables and plot the data.
                    for s in range(len(auto_FFTPeakvariables)):
                        # Plotting the auto data
                        ax_auto = axs_auto[s] if len(auto_FFTPeakvariables) > 1 else axs_auto

                        if(s< 0):
                            # Print the contents of the columns that correspond to PEAK-TIMEPERIOD-REPAVG{...} or PEAK-TIMEPERIOD[AVG[...]]
                            col_repavg = f"PEAK-TIMEPERIOD-REPAVG{{{analType}{{{auto_FFTPeakvariables[s]}}}}}"
                            col_avg = f"PEAK-TIMEPERIOD[AVG[{analType}{{{auto_FFTPeakvariables[s]}}}]]"
                            if col_repavg in auto_data_T0.columns and col_avg in auto_data_T0.columns:
                                print(f"Index | {col_repavg} | {col_avg}")
                                for idx, (repavg_val, avg_val) in zip(auto_data_T0.index, zip(auto_data_T0[col_repavg], auto_data_T0[col_avg])):
                                    print(f"{idx} | {repavg_val} | {avg_val}")
                            elif col_repavg in auto_data_T0.columns:
                                print(f"Index | {col_repavg}")
                                for idx, repavg_val in zip(auto_data_T0.index, auto_data_T0[col_repavg]):
                                    print(f"{idx} | {repavg_val}")
                            elif col_avg in auto_data_T0.columns:
                                print(f"Index | {col_avg}")
                                for idx, avg_val in zip(auto_data_T0.index, auto_data_T0[col_avg]):
                                    print(f"{idx} | {avg_val}")
                            #Next print the entire df
                            #print("\nFull DataFrame for auto_data_T0:")
                            #print(auto_data_T0.head())

                            # Show the columns that correspond to the variables.
                            input("Press Enter to continue...")
                        
                        

                        # For each species, plot the second and third columns
                        print(f"Plotting AUTO {analType}{{{auto_FFTPeakvariables[s]}}}] for {Pre} at T0 = {T0}.")
                        # First plotting a seaborn violinplot 
                        # First plotting the average of all time periods corresponding to individual replicates vs a value.
                        ax_auto.scatter(auto_data_T0.index.get_level_values('a'),
                                        auto_data_T0["PEAK-TIMEPERIOD-REPAVG{" + analType + "{" + auto_FFTPeakvariables[s] + r"}}"],
                                        label = r"$\delta T[\mu_{R}(FFT_{i}($" + auto_FFTPeakvariables[s] + r"$))]$", color = colours[s % len(colours)],
                                        marker = 'o', s = 20, alpha = 0.85)

                        #Next plotting the time period of AVG (average of all replicates timeseries) vs a value.
                        ax_auto.scatter(auto_data_T0.index.get_level_values('a'), 
                                        auto_data_T0["PEAK-TIMEPERIOD[AVG[" + analType + "{" + auto_FFTPeakvariables[s] + "}]]"],
                                        label = r"$\delta T[FFT(\mu_{R}($" + auto_FFTPeakvariables[s] + r"$))]$", color = colours_med[s % len(colours_med)],
                                        marker = 'o', s = 20, alpha = 0.75)
                        
                        violin_df_T0 = pan.DataFrame()  # DataFrame to store the violin plot data for this T0 value.
                        # Plotting the individual replicates in a light grey shade.
                        Rmax = int(auto_data_T0.index.get_level_values('maxR').max())
                        for R in range(0, Rmax):
                            # Check if the column exists and is not empty before plotting.
                            if ( f"PEAK-TIMEPERIOD[{analType}{{{auto_FFTPeakvariables[s]}}}_R_{R}]" in auto_data_T0.columns 
                                and not auto_data_T0[f"PEAK-TIMEPERIOD[{analType}{{{auto_FFTPeakvariables[s]}}}_R_{R}]"].isnull().all() ):
                                #print(f"Plotting R {analType}[{auto_FFTPeakvariables[s]}]_R_{R} for {Pre} at a = {a}, T0 = {T0}.")
                                ax_auto.scatter(auto_data_T0.index.get_level_values('a'),
                                                auto_data_T0[f"PEAK-TIMEPERIOD[{analType}{{{auto_FFTPeakvariables[s]}}}_R_{R}]"],
                                                color = "grey", marker = 'o', s = 15, alpha = 0.35)
                                # Join the data for violin plot to violin_df_T0 DataFrame.
                                tmp_df = auto_data_T0[[f"PEAK-TIMEPERIOD[{analType}{{{auto_FFTPeakvariables[s]}}}_R_{R}]"]].rename(columns={f"PEAK-TIMEPERIOD[{analType}{{{auto_FFTPeakvariables[s]}}}_R_{R}]": "value"})
                                # Next set index of tmp_df to the a value.
                                tmp_df['a'] = auto_data_T0.index.get_level_values('a')
                                # Drop the index and reset it.
                                tmp_df.reset_index(drop=True, inplace=True)
                                violin_df_T0 = pan.concat([violin_df_T0, tmp_df], ignore_index=True)
                                # Set the index to the a value and variable name.
                                #violin_df_T0['a'] = auto_data_T0.index.get_level_values('a')
                        # Report the replicate wide mean, median and std for given "a" value using the violin_df_T0 DataFrame.
                        if not violin_df_T0.empty:
                            print(f"Violin plot statistics for {Pre}, T0 = {T0}:")
                            print(violin_df_T0.groupby('a')['value'].agg(['mean', 'median', 'std']))
                        sea.violinplot(data=violin_df_T0, x='a', y='value', ax=ax_auto,  native_scale= True, fill=False, color=colours_med[s % len(colours_med)], inner='box',
                                    linewidth=1.1, alpha=0.85, cut=0, bw_adjust=0.5, bw_method='scott', density_norm="area")
                        
                        # TODO: Implement 95% CI as infill area for the scatter plot.
                        err = 2*np.sqrt(auto_data_T0["PEAK-TIMEPERIOD-REPVAR{" + analType + "{" + auto_FFTPeakvariables[s] + r"}}"])
                        #ax_auto.fill_between(auto_data_T0.index.get_level_values('a'), 
                        #    auto_data_T0["PEAK-TIMEPERIOD-REPAVG{" + analType + "{" + auto_FFTPeakvariables[s] + r"}}"] - err,
                        #    auto_data_T0["PEAK-TIMEPERIOD-REPAVG{" + analType + "{" + auto_FFTPeakvariables[s] + r"}}"] + err, 
                        #    color = colours_light[s % len(colours_light)], alpha= 0.25)
                        # Set the title and labels for the plot.
                        ax_auto.set_title(r"$\delta T($ " + analType + "[" + auto_FFTPeakvariables[s] + r"]) vs a")
                        ax_auto.set_xlabel(r'R $(mm/hr)$')
                        ax_auto.set_ylabel(r'Recovery Time Period $\delta T$')
                        # Set lower ylim as -5
                        ax_auto.set_ylim(bottom=-1)
                        ax_auto.legend()
                    # End of s loop
                    fig_auto.suptitle(r"Recovery Time Periods for " + analType + f" data at T0 = {T0} for {Pre}, dP = {dP}, Geq = {Geq}")
                    plt.tight_layout()
                    plt.savefig(savepngdir_auto + f"{autosavefilename}_Peak_{harm_indx}.png")
                    plt.show(); plt.close()
                # End of harm_indx loop
            # End of T0 loop
        
        # Now plotting the cross data
        if not cross_FFTPeakdata.empty:
            a_scaled_vals = cross_FFTPeakdata.index.get_level_values('a').unique().to_list()
            focus_T_vals = cross_FFTPeakdata.index.get_level_values('T0').unique().to_list()

            for T0 in focus_T_vals:
                T0 = int(T0) if float(T0).is_integer() else T0

                # Extracting the cross data for the current T0 value.
                try:
                    cross_data_T0 = cross_FFTPeakdata.loc[(slice(None), slice(None), T0, slice(None), slice(None)), :]
                except KeyError:
                    print(f"No CROSS data found for {Pre} at T0 = {T0}. Skipping....")
                    continue
                    
                for harm_indx in harm_peaks_index:
                    # Extract the ith (i % 10) rows from the auto_FFTPeakdata DataFrame.
                    try:
                        cross_data_T0 = cross_FFTPeakdata.loc[(slice(None), slice(None), T0, slice(None)), :]
                    except KeyError:
                        print(f"No CROSS data found for {Pre} at T0 = {T0}. Skipping....")
                        continue
                    # Extract the ith (i % 10) rows from the auto_FFTPeakdata DataFrame.
                    # If harm_indx is -1, then we simply plot all the data for that T0 value.
                    if harm_indx == -1:
                        print(f"Plotting ALL CROSS PEAK data (Superimposed) for {Pre} at T0 = {T0}....")
                    if harm_indx != -1:
                        try:
                            cross_data_T0 = cross_data_T0.iloc[[(i % 10) == harm_indx for i in range(len(cross_data_T0))], :]
                        except IndexError:
                            print(f"No CROSS data found for {Pre} at T0 = {T0} for harm_peaks_index = {harm_indx}. Skipping....")
                            continue
                
                    # Create save directory for the plots.
                    savepngdir_cross = out_dir + f"{Pre}/2DCorrelation/BoxViolin/L_{g}/dP_{dP}/Geq_{Geq}/T0_{T0}/FFTTime/"
                    Path(savepngdir_cross).mkdir(parents=True, exist_ok=True)
                    # Create subplots of shape (p, p) where p is round(sqrt(len(cross_FFTPeakvariables)))
                    #p_cross = len(cross_FFTPeakvariables)// SPB + 1; q_cross = SPB
                    p_cross = len(cross_FFTPeakvariables)//4 + 1; q_cross = 4
                    fig_cross, axs_cross = plt.subplots(p_cross, q_cross, figsize = set_figsizeGrid(p_cross, q_cross), sharex = True, sharey = True)

                    # Flatten the axs array for easier indexing.
                    axs_cross = axs_cross.flatten()
                    # Loop over the variables and plot the data.
                    for s in range(len(cross_FFTPeakvariables)):
                        # Plotting the cross data
                        ax_cross = axs_cross[s] if len(cross_FFTPeakvariables) > 1 else axs_cross
                        # For each species, plot the second and third columns
                        print(f"Plotting CROSS {analType}{{{cross_FFTPeakvariables[s]}}}] for {Pre} at T0 = {T0}.")
                        # First plotting the average of all time periods corresponding to individual replicates vs a value.
                        ax_cross.scatter(cross_data_T0.index.get_level_values('a'),
                                        cross_data_T0["PEAK-TIMEPERIOD-REPAVG{" + analType + "{" + cross_FFTPeakvariables[s] + r"}}"],
                                        label = r"$\delta T[\mu_{R}(FFT_{i}($" + cross_FFTPeakvariables[s] + r"$))]$", color = colours[s % len(colours)],
                                        marker = 'o', s = 20, alpha = 0.85)
                        # Next plotting the time period of AVG (average of all replicates timeseries) vs a value.
                        ax_cross.scatter(cross_data_T0.index.get_level_values('a'),
                                        cross_data_T0["PEAK-TIMEPERIOD[AVG[" + analType + "{" + cross_FFTPeakvariables[s] + "}]]"],
                                        label = r"$\delta T[FFT(\mu_{R}($" + cross_FFTPeakvariables[s] + r"$))]$", color = colours_med[s % len(colours_med)],
                                        marker = 'o', s = 20, alpha = 0.75)
                        # Plotting the individual replicates in a light grey shade.
                        Rmax = int(cross_data_T0.index.get_level_values('maxR').max())
                        for R in range(0, Rmax):
                            # Check if the column exists and is not empty before plotting.
                            if ( f"PEAK-TIMEPERIOD[{analType}{{{cross_FFTPeakvariables[s]}}}_R_{R}]" in cross_data_T0.columns 
                                and not cross_data_T0[f"PEAK-TIMEPERIOD[{analType}{{{cross_FFTPeakvariables[s]}}}_R_{R}]"].isnull().all() ):
                                #print(f"Plotting R {analType}[{cross_FFTPeakvariables[s]}]_R_{R} for {Pre} at a = {a}, T0 = {T0}.")
                                ax_cross.scatter(cross_data_T0.index.get_level_values('a'),
                                                cross_data_T0[f"PEAK-TIMEPERIOD[{analType}{{{cross_FFTPeakvariables[s]}}}_R_{R}]"],
                                                color = "grey", marker = 'o', s = 15, alpha = 0.35)
                                
                        
                        # TODO: Implement 95% CI as infill area for the scatter plot.
                        err = 2*np.sqrt(cross_data_T0["PEAK-TIMEPERIOD-REPVAR{" + analType + "{" + cross_FFTPeakvariables[s] + r"}}"])
                        ax_auto.fill_between(cross_data_T0.index.get_level_values('a'), 
                            cross_data_T0["PEAK-TIMEPERIOD-REPAVG{" + analType + "{" + cross_FFTPeakvariables[s] + r"}}"] - err,
                            cross_data_T0["PEAK-TIMEPERIOD-REPAVG{" + analType + "{" + cross_FFTPeakvariables[s] + r"}}"] + err, 
                            color = colours_light[s % len(colours_light)], alpha= 0.25)
                        # Set the title and labels for the plot.
                        ax_cross.set_title(r"$\delta T($ " + analType + "[" + cross_FFTPeakvariables[s] + r"]) vs a")
                        ax_cross.set_xlabel(r'R $(mm/hr)$')
                        ax_cross.set_ylim(bottom=-1)  # Set lower ylim as -5
                        ax_cross.set_ylabel(r'Interaction Time Scales $\delta T$')
                        ax_cross.legend()
                    # End of s loop
                    fig_cross.suptitle(r"Interaction Time Scales for " + analType + f" data at T0 = {T0} for {Pre}, dP = {dP}, Geq = {Geq}")
                    plt.tight_layout()
                    plt.savefig(savepngdir_cross + f"{crosssavefilename}.png")
                    plt.show(); plt.close()
                # End of T0 loop
        # End of cross data check




def analyse_FRAME_POTdata(indir, out_dir, prefixes=[], a_vals=[], T_vals =[], find_minima= True, filename = "Pot_Well.csv", 
                          minimafilename = "LOCAL_MINIMA.csv", var_labels= [ "P(x; t)" , "G(x; t)", "Pr(x; t)", "O(x; t)"], a_scaling = 1):

    savefilename = re.sub(r'\.csv$', '', filename)
    saveminimafilename = re.sub(r'\.csv$', '', minimafilename)
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data = pan.DataFrame()
    sea.set_palette("husl")
    colours = [hex_list[i][-1] for i in range(len(hex_list))]

    init_a_vals = a_vals; init_T_vals = T_vals
    

    for Pre in prefixes:
        
        data, local_minima = get_FRAME_POTdata(indir, Pre, init_a_vals, init_T_vals, a_scaling, find_minima, filename, minimafilename, "Pot_Well/")
        # We are interested in the AVG columns for each species.
        # Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        # Data has columns for each species given by {var}, of the form:
        # R_max   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]_R_0   ...   AVG[{var}]_R_{R_max}
        if data.empty:
            print(f"No data found for {Pre}. Skipping....")
            continue
        if local_minima.empty and find_minima:
            print(f"WARNING: No local minima data found for {Pre}. Skipping....")
            continue
        
        print("====================================================================")
        print(f"Data for {Pre}:")
        print("====================================================================\n")
        #print(data.info())
        print(data.head())
        print(data.tail())
        print(data.columns)
        print(data.index)
        #print(data.describe())

        # Find Tmax as the maximum value of T in the data.
        Tmax = data.index.get_level_values('T').max()
        Tmax = int(Tmax) if float(Tmax).is_integer() else Tmax
        savecsvdir = out_dir + f"{Pre}/Pot_Well/"
        Path(savecsvdir).mkdir(parents=True, exist_ok=True)

        # Saving minima data to a csv file if it exists.
        if not local_minima.empty:
            local_minima.to_csv(savecsvdir + f"{saveminimafilename}_dP_{dP}_Geq_{Geq}.csv")
            
        # Save the data to a csv file.
        data.to_csv(savecsvdir + f"{savefilename}_dP_{dP}_Geq_{Geq}.csv")
        # Use the data to create plots.

        # The plots are created for each species for each Prefix, a value and T value, for each species.

        a_scaled_vals = data.index.get_level_values('a').unique().to_list()
        t_vals_list = data.index.get_level_values('T').unique().to_list()
        for a in a_scaled_vals:
            a = int(a) if float(a).is_integer() else a
            for T in t_vals_list:
                T = int(T) if float(T).is_integer() else T
                savepngdir = out_dir + f"{Pre}/Pot_Well/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/"
                Path(savepngdir).mkdir(parents=True, exist_ok=True)

                try:
                    data_A = data.loc[(slice(None), a, T, slice(None)), :]
                except KeyError:
                    print(f"No data found for {Pre} at a = {a}, T = {T}. Skipping....")
                    continue
                fig, axs = plt.subplots(1, len(var_labels), figsize = set_figsize(len(var_labels)))
                for s in range(len(var_labels)):
                    
                    # Plotting for a fixed a value, species and Prefix for all R values.
                    ax = axs[s]
                    # For each species, plot the first and third columns
                    ax.plot(data_A["BINS[" + var_labels[s] + "]"], data_A["POTENTIAL_WELL[" + var_labels[s] + "]"], label = r"$V(" + var_labels[s] + ")$",
                            color = colours[s], linestyle = 'solid', linewidth = 2, alpha = 0.85)
                    if find_minima:
                        # If Find_minima is set to True, plot the local minima data, as diamond points.
                        ax.plot(local_minima.loc[(slice(None), a, T, slice(None)), "LOCAL_MINIMA[" + var_labels[s] + "]"], 
                                local_minima.loc[(slice(None), a, T, slice(None)), "MIN_POT[" + var_labels[s] + "]"], 
                                marker = 'D', color = colours[s], linestyle = 'None', alpha = 0.9, label = "Local Minima")
                    ax.set_title(r"$V(" + var_labels[s] + ")$ vs " + r"$\langle $" + var_labels[s] + r"$(x, t) \rangle_{x}$" + f" at a = {a}, T = {T}")
                    ax.set_xlabel(r'$\rho( ' + var_labels[s] + ')$')
                    ax.set_ylabel(r'$V(' + var_labels[s] + ')$')
                    ax.legend()
                # Noting maxR for the plot title.
                maxR = data_A.index.get_level_values('MaxR').max()
                fig.suptitle(r"$V(\rho(x, t))$" + f" For {Pre} at a = {a}, dP = {dP}, Geq = {Geq}, T = {T}, n = {maxR}")
                plt.savefig(savepngdir + f"{savefilename}_a_{a}_T_{T}_dP_{dP}_Geq_{Geq}.png")
                #plt.show()
                plt.close()
                #print(f"Plot saved to {savepngdir + f'{savefilename}_a_{a}_T_{T}_dP_{dP}_Geq_{Geq}.png'}")
            # End of T loop
        # End of a loop
        # Use home_video function to create a video of the plots.
        home_video(out_dir, out_dir, prefixes=[f"{Pre}/Pot_Well"], a_vals= a_scaled_vals, T_vals=T_vals, maxR= 1, 
                   pngformat= "%s_a_{a}_T_{T}_dP_{dP}_Geq_{Geq}.png" %(savefilename), pngdir= "{Pre}/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/", 
                   videoname = f"{savefilename}_{Pre}_Tmax_{Tmax}.mp4", video_relpath= "{Pre}/Videos/{amin}-{amax}/")
        # Note that there is only one png per a value, so the video will be a simple slideshow (and R_max = 1).
        input("Press F to Continue...")
 

""" # Summary description of the analyse_FRAME_CLUSTEREDISTdata(...)
Reads frequency distribution of clustered data from csv file given by filename, for various prefixes, a and T values.
Then plots the data as lineplots with size of cluster on the x-axis and frequency on the y-axis (a value  serving as hue for different lines).
This plot is made for each prefix and T value.
"""
def analyse_FRAME_CLUSTEREDISTdata(indir, out_dir, prefixes=[], a_vals=[], T_vals =[], filename = "CLUSTERED_FREQUENCIES.csv", 
            clustsubdir ="CLUST/{binType}_{nbins}/", binType="KMeans", nbins=2, var_labels= [ "P(x; t)" , "G(x; t)", "Pr(x; t)", "O(x; t)"],
             plot_binframes = False, CDF=True, a_c=None,  a_scaling = 1):
    
    savefilename = re.sub(r'\.csv$', '', filename) + f"a_c_{a_c[0]}_FINAL" if a_c is not None else re.sub(r'\.csv$', '', filename) + "_FINAL"
    clustsubdir = clustsubdir.format(binType=binType, nbins=nbins)
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data_auto = pan.DataFrame(); combined_data_cross = pan.DataFrame()
    sea.set_palette("husl")

    colours = [hex_list[i][-1] for i in range(len(hex_list))]
    colours_med = [hex_list[i][-4] for i in range(len(hex_list))]
    colours_light = [hex_list[i][1] for i in range(len(hex_list))]

    # Replace "#ff447d" in colours_med with "#FFCC00" if it exists.
    colours_med = [c if c != "#ff447d" else "#FFCC00" for c in colours_med]

    #colours_med = trueblues_hex_list 

    init_a_vals = a_vals; init_T_vals = T_vals

    for Pre in prefixes:

        data = get_FRAME_CLUSTERdata(indir, Pre, init_a_vals, init_T_vals, a_scaling, freqfilename=filename, 
                                     clustsubdir= clustsubdir, binType=binType, nbins=nbins)
        # Returns a MultiIndex Dataframe with Indices: Prefix, a, T and columns for each species given by {var}, of the form:
        # CLUSTER_SIZE[{var1}], CLUSTER_COUNTS[{var1}], CLUSTER_FREQ[[{var1}] ... CLUSTER_FREQ[[{varN}]
        # where {var} is one of the species names.
        if data.empty:
            print(f"No data found for {Pre}. Skipping....")
            continue

        print("====================================================================")
        print(f"Data for {Pre}:")
        print("====================================================================\n")
        #print(data.info())
        print(data.head())
        print(data.tail())
        print(data.columns)
        print(data.index)

        savecsvdir = out_dir + f"{Pre}/{clustsubdir}/"
        Path(savecsvdir).mkdir(parents=True, exist_ok=True)
        # Save the data to a csv file.
        data.to_csv(savecsvdir + f"{savefilename}_L_{g}_dP_{dP}_Geq_{Geq}.csv")

        a_scaled_vals = data.index.get_level_values('a').unique().to_list()
        print(f"Unique a values found in data: {a_scaled_vals}")
        t_vals_list = data.index.get_level_values('T').unique().to_list()
        amin = min(a_scaled_vals); amax = max(a_scaled_vals)
        Tmax = max(t_vals_list)
        Tmax = int(Tmax) if float(Tmax).is_integer() else Tmax

        # If a_c is provided, and is iterable, then scale each value by a_scaling.
        if a_c is not None:
            if isinstance(a_c, (list, tuple, np.ndarray)):
                a_c = [ac * a_scaling for ac in a_c]
            else:
                # Make a_c into a singleton list.
                a_c = [a_c * a_scaling]
            print(f"Critical a value(s) provided: {a_c}")

        for T in t_vals_list:
            T = int(T) if float(T).is_integer() else T

            # We will plot the data CLUSTER_SIZE[{var}] vs CLUSTER_FREQ[{var}] for each species {var} 
            # and a value as hue, with a separate plot for each Prefix and T value.
            savepngdir = out_dir + f"{Pre}/{clustsubdir}/L_{g}/dP_{dP}/Geq_{Geq}/T_{T}/"
            Path(savepngdir).mkdir(parents=True, exist_ok=True)

            #if a_c is not None:
            #    all_txts=[] # List to hold all text annotations for adjustText

            data_T = data.loc[(slice(None), slice(None), T), :] # Extract data for the current T value.
            fig, axs = plt.subplots(1, len(var_labels), figsize= set_figsize(len(var_labels)))#=(7.5,6))
            for s in range(len(var_labels)):
                # First extract the columns for the current species {var}, and drop any NaN values.
                col_size = f"CLUSTER_SIZE[{var_labels[s]}]";
                col_freq = f"CLUSTER_CDF[{var_labels[s]}]" if CDF else f"CLUSTER_FREQ[{var_labels[s]}]"
                col_counts= f"CLUSTER_COUNTS[{var_labels[s]}]"
                if col_size not in data_T.columns or col_freq not in data_T.columns or col_counts not in data_T.columns:
                    print(f"Columns {col_size} or {col_freq} or {col_counts} not found in data for {Pre} at T = {T}. Skipping....")
                    continue
                data_T_var = data_T[[col_size, col_counts, col_freq]].dropna()

                # Plotting the data for the current species {var} and T value, across all a values.
                ax = axs[s] if len(var_labels) > 1 else axs
                for a in a_scaled_vals:
                    a = int(a) if float(a).is_integer() else a
                    try:
                        data_T_var_a = data_T_var.loc[(slice(None), a), :]
                    except KeyError:
                        print(f"Data for a = {a} not found in data for {Pre} at T = {T}. Skipping....")
                        continue
                    # Plotting the CLUSTER_SIZE vs CLUSTER_FREQ for the current species {var} and a value.
                    if CDF:
                        
                        # If a_c is provided and matches the current a value, estimate and plot the power-law fit.
                        if a_c is not None and a in a_c:

                            # Plot the CCDF normally.
                            ax.plot(data_T_var_a[col_size], 1 - data_T_var_a[col_freq], #label = f"a = {a}",
                            #color = "#915f78", linestyle = 'solid', linewidth = 1.3, alpha = 0.6)
                            color = colours_med[a_scaled_vals.index(a) % len(colours_med)], linestyle = 'solid', linewidth = 1.3, alpha = 0.6)
                            # Superimpose scatter points for the data.
                            ax.scatter(data_T_var_a[col_size], 1 - data_T_var_a[col_freq], label = sci_notation(a,3, None, -3), #f"{a:.4g}", #color = "#915f78",
                                #color = "#915f78",
                                color = colours_med[a_scaled_vals.index(a) % len(colours_med)],
                                marker = 'o', s = 15, alpha = 0.65)
                            
                            counts = data_T_var_a[col_counts].astype(int)
                            if counts.sum() == 0:
                                print(f"WARNING: 0 counts detected for {Pre} at a = {a}, T = {T}. Skipping power-law fit....")
                                continue
                            samples = pan.Series(np.repeat(data_T_var_a[col_size].values, counts.values.astype(int)))
                            fit = pwl.Fit(samples, discrete=True, verbose=True) #,xmin=3.8e3, xmin=7e3) #xmin=5, xmax=0.6e4
                            ##, xmax=5e4) #xmin=1.9e2, xmax=1.535e3) xmin=1e2, xmax=7e3)
                            print(f"Power-law fit for {Pre} at a = {a}, T = {T}: alpha = {fit.alpha:.4g}, xmin = {fit.xmin}, D = {fit.D:.4g}")

                            # Plot the power-law fit as a dashed line.
                            #fit.power_law.plot_ccdf(ax=ax, color=colours[a_scaled_vals.index(a) % len(colours)], linestyle='dashed', linewidth=1.5, 
                            #                  label=r'$N_{S}(T_{eq})[s \geq S] \propto S^{-(\alpha - 1)}, \alpha = %.4g$' % (fit.alpha))
                            #fit.alpha = 2.08
                            # Manually plot the power-law fit for the CCDF.
                            # Note for a power-law distribution P(x) = Cx^(-alpha), the CCDF is given by P(X > x) ~ (C/(alpha - 1))*x^(-(alpha - 1))
                            # We can estimate C from our data as: C = {(1 - data_T_var_a[col_freq]).iloc[data_T_var_a[col_size].searchsorted(fit.xmin)] * (fit.xmin**(fit.alpha - 1))}*(fit.alpha - 1)
                            C = ((1 - data_T_var_a[col_freq]).iloc[data_T_var_a[col_size].searchsorted(fit.xmin)] * (fit.xmin**(fit.alpha - 1)))*(fit.alpha - 1)
                            #x_fit = np.linspace(fit.xmin, data_T_var_a[col_size].max()*0.1, 100) - 50
                            x_fit = np.linspace(fit.xmin, 5e4, 100)  #-35
                            y_fit = (C/(fit.alpha - 1))*x_fit**(-(fit.alpha - 1))*1.02  # Slightly scale up y_fit to avoid overlap with data points.
                            ax.plot(x_fit, y_fit, color="indianred", linestyle='dashed', linewidth=3)
                            # Add axis text with automatically placed text box using adjustText
                            txts = []
                            txts.append(ax.text(0.05, 0.05, r'$P_{S}(T_{eq})[s \geq S] \propto S^{(1 -\tau)}, \tau = %.4g \pm %.3f$' % (fit.alpha, fit.sigma),
                                    transform=ax.transAxes, fontsize=14.5, color="indianred",
                                    bbox=dict(boxstyle='round,pad=0.25', edgecolor='indianred', facecolor='white', alpha=0.6)))
                            adjT.adjust_text(txts, ax=ax)
                            #ax.text(0.025, 0.025, r'$N_{S}(T_{eq})[s \geq S] \propto S^{-(\alpha - 1)}, \alpha = %.4g \pm %.3f$' % (fit.alpha, fit.sigma),
                            #        fontsize=10, color="indianred", 
                            #        verticalalignment='bottom', horizontalalignment='left',
                            #        bbox=dict(boxstyle='round,pad=0.25', edgecolor='indianred', facecolor='white', alpha=0.6))
                        else:
                            # Plot the CCDF normally.
                            ax.plot(data_T_var_a[col_size], 1 - data_T_var_a[col_freq], #label = f"a = {a}",
                            color = colours_med[a_scaled_vals.index(a) % len(colours_med)], linestyle = 'solid', linewidth = 1.25, alpha = 0.6)
                            # Superimpose scatter points for the data.
                            ax.scatter(data_T_var_a[col_size], 1 - data_T_var_a[col_freq], label = sci_notation(a,3, None, -3), #f"{a:.4g}",
                                       color = colours_med[a_scaled_vals.index(a) % len(colours_med)], marker = 'o', s = 15, alpha = 0.65)

                    else:
                        ax.plot(data_T_var_a[col_size], data_T_var_a[col_freq], #label = f"a = {a}",
                            color = colours[a_scaled_vals.index(a) % len(colours)], linestyle = 'solid', linewidth = 1.25, alpha = 0.6)
                        # Superimpose scatter points for the data.
                        ax.scatter(data_T_var_a[col_size], data_T_var_a[col_freq], label = f"a = {a:.4g}", color = colours_med[a_scaled_vals.index(a) % len(colours_med)],
                               marker = 'o', s = 10, alpha = 0.5)
                # Set the title and labels for the plot.
                #ax.set_title(r"$P_{S}(T_{eq}) $" + var_labels[s] + f" at T = {T}")
                ax.set_xlabel(r'Cluster Size ($S$)', fontsize=16)
                if CDF:
                    ax.set_ylabel(r'$P_{S}(T_{eq})[s \geq S]$',  fontsize=16)
                else:
                    ax.set_ylabel(r'Normalised Frequency $P_{S}(T_{eq})$',  fontsize=16)
                # Set both x and y axes to log scale.
                ax.set_xscale('log'); ax.set_yscale('log')
                #Set x and y ticks fontsize
                ax.tick_params(axis='x', labelsize=14) ;ax.tick_params(axis='y', labelsize=14)
                # Set min y lim to 1e-4
                ax.set_ylim(bottom=1e-4, top=1)
                ax.legend(fontsize=11.5, title=r"Vegetation Growth Rate $a$ (hr$^{-1}$)", title_fontsize=12)
            # End of s loop
            # If a_c is provided, use adjustText to avoid overlapping text annotations.
            
            # Generate fig suptitle with Prefix, dP, Geq and T values.
            fig.suptitle(r"Cluster Size Distribution for " + f"{Pre} at T = {T}, dP = {dP}, Geq = {Geq}")
            plt.tight_layout()
            if(CDF):
                plt.savefig(savepngdir + f"{savefilename}_CCDF_a_{amin}-{amax}.png", transparent=True, dpi = 800)
            else:
                plt.savefig(savepngdir + f"{savefilename}_a_{amin}-{amax}.png", transparent=True, dpi = 800)
            plt.show(); plt.close()

        # End of T loop
        # Use home_video function to create a video of the plots.
        if CDF:
            pngformat = "%s_CCDF_a_{amin}-{amax}.png" %(savefilename)
        else:
            pngformat = "%s_a_{amin}-{amax}.png" %(savefilename)
        home_video(out_dir, out_dir, prefixes=[f"{Pre}/{clustsubdir}"], a_vals= [amin, amax], T_vals=T_vals, maxR= 1,
                     pngformat= pngformat, pngdir= "{Pre}/L_{g}/dP_{dP}/Geq_{Geq}/T_{T}/", 
                     videoname = f"{savefilename}_{Pre}_Tmax_{Tmax}.mp4", video_relpath= "{Pre}/Videos/{amin}-{amax}/")






def analyse_PRELIMS_TIMESERIESdata(indir, out_dir, prefixes=[], a_vals=[], tseries_vals =[], serType = "TSERIES",
                                   meanfilename = "Mean_{serType}_T_{TS}.csv", var_labels= [ "<<P(x; t)>_x>_r" , "<<G(x; t)>_x>_r", "<<Pr(x; t)>_x>_r"], a_c=None, a_scaling = 1):
    savefilename = re.sub(r'\.csv$', '', meanfilename.format(Pre="PREFIX", g="g", dP="dP", Geq="Geq", a="a", serType = serType, TS=tseries_vals[0]))
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data = pan.DataFrame()
    sea.set_palette("husl")
    colours = [hex_list[i][-1] for i in range(len(hex_list))]

    init_a_vals = a_vals; init_TS_vals = tseries_vals;
    include_col_labels = ["t"] + var_labels
    # Replace "(x; t)" with "(x,t)" in include_col_labels and append "$" at the start and end of each label for proper LaTeX formatting.
    dict_varlabels = {}
    for label in include_col_labels:
        new_label = label.replace("(x; t)", r"(x,t)"); new_label = new_label.replace("<", r"\langle ").replace(">", r" \rangle")
        dict_varlabels[label] = r"$" + new_label + r"$"

    for Pre in prefixes:
        

        data = get_PRELIM_EQdata(indir, Pre, init_a_vals, init_TS_vals, a_scaling, meanfilename = meanfilename, include_col_labels= include_col_labels, serType = serType)
        # We are interested in the AVG columns for each species.
        # Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        # Data has columns for each species given by {var}, of the form:
        # AVG[{var}]_ALL ... VAR[{var}]_ALL ...   AVG[{var}]_SURV   ...  VAR[{var}]_SURV   ... t ... {var}
        # where {var} is one of the species names. The dataframe is multi-indexed by Prefix, a, R, Tmax,
        # where R = -1 for the mean values.
        if data.empty:
            print(f"No data found for {Pre}. Skipping....")
            continue
        print("====================================================================")
        print(f"Data for {Pre}:")
        print("====================================================================\n")
        #print(data.info())
        print(data.head())
        print(data.tail())
        print(data.columns)
        print(data.index)
        #print(data.describe())

        # Find Tmax as the maximum value of T in the data.
        Tmax = data.index.get_level_values('Tmax').max()
        Tmax = int(Tmax) if float(Tmax).is_integer() else Tmax
        savecsvdir = out_dir + f"{Pre}/TimeSeries/"
        Path(savecsvdir).mkdir(parents=True, exist_ok=True)

        # Save the data to a csv file.
        data.to_csv(savecsvdir + f"{savefilename}_dP_{dP}_Geq_{Geq}.csv")
        #Concatenate the data to the combined data
        combined_data = pan.concat([combined_data, data], axis = 0)
        # Use the data to create plots.
        # The plots are created for each species for each Prefix and a value, with T as the x-axis.
        # For each species, create a plot of T vs all the data for that species as lines, with the mean across surviving and all replicates as darker lines.
        # The means for individual replicates are lighter lines.
        
        # Try determining number of species by counting all columns with _SURV in the name.
        try:
            num_species = len([col for col in data.columns if "_SURV" in col])/3
        except IndexError:
            print(f"No columns with _SURV found in the data. Not plotting data for Prefix {Pre}. Skipping....")
            num_species = None
            continue
        '''
        # Find column index of first column with _R_0 in the name. If not found, print a warning and set to 0.
        try:
            R0_index = data.columns.get_loc([col for col in data.columns if "_R_0" in col][0])
        except IndexError:
            print(f"No column with _R_0 found in the data. Not plotting individual replicate data for Prefix {Pre}.")
            R0_index = None
        '''
        print(f"num_species = {num_species}")
        # Get the list of a values and T values from the data.
        a_scaled_vals = data.index.get_level_values('a').unique().to_list()
        #Force a to be represented in decimal form rather than scientific notation (so 0.00001 is 0.00001 rather than 1e-5).
        tseries_vals = data.index.get_level_values('Tmax').unique().to_list()
        print(a_scaled_vals, tseries_vals)

        # If a_c is provided, and is iterable, then scale each value by a_scaling.
        if a_c is not None:
            if isinstance(a_c, (list, tuple, np.ndarray)):
                a_c = [ac * a_scaling for ac in a_c]
            else:
                # Make a_c into a singleton list.
                a_c = [a_c * a_scaling]
            print(f"Critical a value(s) provided: {a_c}")
        
        for a in a_scaled_vals:
            #a = int(a) if float(a).is_integer() else a
            #Force a to be represented in decimal form rather than scientific notation.
            
            for TS in tseries_vals:
                savepngdir = out_dir + f"{Pre}/TimeSeries/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T_{TS}/"
                Path(savepngdir).mkdir(parents=True, exist_ok=True)
                fig, axs = plt.subplots(1, len(var_labels), figsize = set_figsize(len(var_labels)))
                for s in range(len(var_labels)):
                    data_A_TS = data.loc[(slice(None), a, slice(None), TS), :]
                    # Plotting for a fixed a value, Tmax value species and Prefix for all R and t values.
                    ax = axs[s] if SPB > 1 else axs
                    # RECALL, _SURV and _ALL columns are given by multi_index for R = -1.
                    data_all = data_A_TS.loc[(slice(None), a, -1, TS), :]
                    # Plotting mean of surviving replicates ( with _SURV in the name)
                    if serType == "MOVSERIES":
                        label_surv = r"$\mu_{surv}($" +f"{var_labels[s]}" +r")$"
                        label_all = r"$\mu_{all}($" +f"{var_labels[s]}" +r")$"
                    else:
                        label_surv = r"$\mu_{surv}(\rho_{%g})$" %(s)
                        label_all = r"$\mu_{all}(\rho_{%g})$" %(s) #Default label for TSERIES
                    ax.plot(data_all["t"], data_all["AVG[" + var_labels[s] + "]_SURV"], label = r"$\mu_{surv}($" + f"{var_labels[s]}" +r"$)$", color = colours[s], linestyle = 'solid', linewidth = 2, alpha = 0.9)
                    # Next, if a in a_c, use power-law fitting to plot the expected scaling behaviour.
                    if a_c is not None and a in a_c and s== 0:
                        # Estimate power-law fit between surviving mean data over time
                        pwl_data = data_all[["t", "AVG[" + var_labels[s] + "]_SURV"]].dropna()
                        pwl_data= pwl_data[pwl_data["AVG[" + var_labels[s] + "]_SURV"] > 0]
                        # Remove all negative entries from pwl_data
                        tmin=200; tmax =10000
                        pwl_data = pwl_data[pwl_data["t"] >= tmin]  if tmin is not None else pwl_data
                        pwl_data = pwl_data[pwl_data["t"] < tmax]  if tmax is not None else pwl_data
                        # Create fit:
                        try:
                            # Initial guesses.
                            A0 = pwl_data["AVG[" + var_labels[s] + "]_SURV"].iloc[0]; exp = 0.451
                            popt, pcov = curve_fit(power_lawfit, pwl_data["t"].values, pwl_data["AVG[" + var_labels[s] + "]_SURV"].values, p0=[A0, exp], maxfev=10000)

                            C_fit, delta_fit = popt
                            C_err, delta_err = np.sqrt(np.diag(pcov))
                            print(f"SciPy fit for {Pre} at a={a}, TS={TS}, species={var_labels[s]}:")
                            print(f"    C = {C_fit:.4g} ± {C_err:.4g}")
                            print(f"    δ = {delta_fit:.4g} ± {delta_err:.4g}")
                            x_fit= np.linspace(pwl_data["t"].min(), pwl_data["t"].max(), 100 )
                            y_fit = power_lawfit(x_fit, C_fit, delta_fit)
                            print(f"Plotting...")
                            ax.plot(x_fit, y_fit, color="indianred", linestyle='dashed', linewidth=2.25)

                            # Add axis text with automatically placed text box using adjustText
                            txts = []
                            txts.append(ax.text(0.05, 0.05, r'$P_{S}(\rho_{surv}($' + f"{dict_varlabels[s]}" +r"$[t]) \propto S^{(-\delta)}, \delta = %.4g \pm %.3f$" % (delta_fit, delta_err),
                                transform=ax.transAxes, fontsize=13, color="indianred",
                                bbox=dict(boxstyle='round,pad=0.25', edgecolor='indianred', facecolor='white', alpha=0.6)))
                            adjT.adjust_text(txts, ax=ax)
                        except Exception as e:
                            print(f"SciPy curve_fit failed for a={a}, TS={TS}, var={var_labels[s]} with Error: {e}")

                    
                    # Plotting mean of all replicates ( with _ALL in the name)
                    ax.plot(data_all["t"], data_all["AVG[" + var_labels[s] + "]_ALL"], label = r"$\mu_{all}($" + f"{var_labels[s]}" +r"$)$", color = colours[s], linestyle = 'dashed', linewidth = 1.35, alpha = 0.75)

                    # Infill variance for surviving replicates
                    err = 2*np.sqrt(data_all["VAR[" + var_labels[s] + "]_SURV"])
                    ax.fill_between(data_all["t"], data_all["AVG[" + var_labels[s] + "]_SURV"] - err, 
                                    data_all["AVG[" + var_labels[s] + "]_SURV"] + err, color = colours[s], alpha = 0.25)
                    
                    # Now plot the mean of individual replicates (suppress labels that represent columns with solely NaN values).
                    for R in data_A_TS.index.get_level_values('R').unique():
                        if R == -1: continue
                        data_R = data_A_TS.loc[(slice(None), a, R, TS), :]
                        ax.plot(data_R["t"], data_R[var_labels[s]], color ="grey", linestyle = 'solid', alpha = 0.2)
                    # End of R loop
                    #Set y-axis as log scale.
                    #ax.set_yscale('log')
                    #ax.set_xscale('log')
                    #ax.set_xlim(1e2, 1.25e5)
                    # Display x-axis labels to be 
                    
                    if serType == "MOVSERIES":
                        #ax.set_xlim(7.5e4, 1.5e5)
                        ax.set_title(f"{var_labels[s]}" + r" vs $t$"  + f" at a = {a}")
                        ax.set_ylabel(f"{var_labels[s]}")
                        #ax.set_ylim(0-1e-2, 1.01)
                    else:
                        ax.set_yscale('log')
                        ax.set_xscale('log')
                        ax.set_title(r"$\langle \rho_{%g}(t) \rangle_{x}$ vs $t$" %(s)  + f" at a = {a}")
                        ax.set_ylabel(r'$\langle \rho_{%g}(t) \rangle_{x}$' % (s))
                    ax.set_xlabel(r'Time $(d)$')
                    ax.legend()
                # End of s loop
                #Set ymax for axs[3] to be 500.
                #if(SPB == 3):
                #    axs[2].set_ylim(-5, 500)
                #axs[2].set_ylim(-5, 500)
                fig.suptitle(r"$\mu(\rho(x, t))$" + f" For {Pre} at a = {a}, dP = {dP}, Geq = {Geq}")
                plt.savefig(savepngdir + f"TimeSeries_{serType}_a_{a}_dP_{dP}_Geq_{Geq}.png")
                #plt.show()
                plt.close()
            # End of TS loop
        
        # End of a loop
        # Use home_video function to create a video of the plots.
        home_video(out_dir, out_dir, prefixes=[f"{Pre}/TimeSeries"], a_vals= a_scaled_vals, T_vals=[Tmax], maxR= 1, 
                   pngformat= f"TimeSeries_{serType}" + "_a_{a}_dP_{dP}_Geq_{Geq}.png", pngdir= "{Pre}/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/", 
                   videoname = f"TimeSeries_{serType}_{Pre}_Tmax_{Tmax}.mp4", video_relpath= "{Pre}/Videos/{amin}-{amax}/")
        # Note that there is only one png per a value, so the video will be a simple slideshow (and R_max = 1).
        
        input("Press F to Continue...")



'''# Summary Description of analyse_PRELIMS_TRAJECTORYdata(...)
This function analyses the trajectory data for each species and creates plots for each species for each Prefix, a value and TS (Time Series) value.
These plots are trajectory plots (for each individual replicate) of the various species (designated by the x_label, y_label) over time with the hue and the size of the points given by
the hue_label and size_label respectively.

NOTE: This function assumes that the data is stored in the meanfilename in the following format (which is the format output by the get_PRELIM_EQdata(...) function):
AVG[{var}]_ALL ... VAR[{var}]_ALL ...   AVG[{var}]_SURV   ...  VAR[{var}]_SURV   ... t ... {var}
where {var} is one of the species names. The dataframe is multi-indexed by Prefix, a, R, Tmax, which when read in correspond to the first four columns of the data.
The multi-index R = -1 corresponds to the mean values. If the data is not found, the function will call get_PRELIM_EQdata(...) to generate the data.

By plotting the {var} columns given by the x_label, y_label, hue_label and size_label, we can create a trajectory plot for each species for each Prefix, a value, TS and R value.
The hue_label and size_label are optional, and if not provided, the hue and size of the points will be set to the same value for all points.
Additionally, trajectories for indivudual replicates are plotted both as distinct plots and as a combined plot for all replicates in [minR, maxR] (with R >= 0) in all cases.

NOTE: reduce_points (by default None) can be set to an integer value to reduce the number of points plotted for each trajectory by that factor (e.g. reduce_points = 10 will plot every 10th point).
      A special case is reduce_points = 0 or reduce_points = "auto", which will automatically calculate a suitable reduction factor based on the binned size (if size_label is provided) or the total number of points (if size_label is not provided).

ALSO NOTE: maxR= -1 will plot all encountered R values in data, while minR= -1 will additionally plot the mean values across all R values.
'''

def analyse_PRELIMS_TRAJECTORYdata(indir, out_dir, prefixes=[], a_vals=[], tseries_vals =[],  x_label = ["<<G(x; t)>_x>_r"], y_label = ["<<Pr(x; t)>_x>_r"],
        hue_label = ["<<P(x; t)>_x>_r"], size_label = [], maxR=-1, minR=0, T_window = [200, 5000], meanfilename = "Mean_TSERIES_T_{TS}_dP_{dP}_Geq_{Geq}.csv", 
        reduce_points=None, plot_mean_trajectory= False, a_scaling = 1):
    
    # Set a white grid background for the plots.
    sea.set_style("darkgrid")
    
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data = pan.DataFrame()
    sea.set_palette("husl")
    colours = [hex_list[i][-1] for i in range(len(hex_list))]

    init_a_vals = a_vals; init_TS_vals = tseries_vals;
    include_col_labels = x_label + y_label + hue_label + size_label
    #Remove duplicates from include_col_labels if any.
    include_col_labels = list(dict.fromkeys(include_col_labels))
    print(f"include_col_labels = {include_col_labels}")
    # Create dictionary to map species names to more readable labels for plotting.
    dict_varlabels = {
        "<<P(x; t)>_x>_r": r"Mean Vegetation Density $\langle V(t) \rangle_X$, ($kg\: km^{-2}$)", "<<G(x; t)>_x>_r": r"Mean Grazer Density $\langle G \rangle_X$, ($kg\: km^{-2}$)",
        "<<Pr(x; t)>_x>_r": r"Mean Predator Density $\langle Pr(t) \rangle_X$, ($kg\: km^{-2}$)", "<GAM[G(x; t)]>_x": r"$\langle \gamma_G(t) \rangle_x$",
        "<GAM[Pr(x; t)]>_x": r"$\langle \gamma_{Pr}(t) \rangle_x$"}

    for Pre in prefixes:
        readfilename =  meanfilename.format(Pre=Pre, g=g, dP=dP, Geq=Geq, TS = tseries_vals[0])#, TS =tseries_vals[0])
        readcsvdir = out_dir + f"{Pre}/TimeSeries/"
        # We are interested in the AVG columns for each species.
        # Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        # Data has columns for each species given by {var}, of the form:
        # t   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]   ...
        # where {var} is one of the species names. The dataframe is multi-indexed by Prefix, a, R, TS_val,
        # where R = -1 for the mean values.
        #Check if the file exists
        if not os.path.exists(readcsvdir + readfilename):
            print(f"File {readcsvdir + readfilename} not found. Generating data for {Pre}....")
            data = get_PRELIM_EQdata(indir, Pre, init_a_vals, init_TS_vals, a_scaling, include_col_labels= include_col_labels)
            if data.empty:
                print(f"No data found for {Pre}. Skipping....")
                continue
        else:
            # Data has columns for each species given by {var}, of the form:
            # t   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]   ...
            # where {var} is one of the species names. The dataframe is multi-indexed by Prefix, a, R, TS_val,
            # where R = -1 for the mean values.
            data = pan.read_csv(readcsvdir + readfilename, header = 0, index_col = [0, 1, 2, 3])
        '''
        print("====================================================================")
        print(f"Data for {Pre}:")
        print("====================================================================\n")
        #print(data.info())
        print(data.head())
        print(data.tail())
        print(f"Shape = {data.shape}")
        print(data.columns)
        print(data.index)
        print("====================================================================\n")
        '''

        # Next we need to filter the data to only include the T values (given by the "t" column) in the T_window.
        # If T_window is not provided, we will use the entire range of T values.
        # If no "t" column values fit the T_window, we will skip the Prefix.
        if T_window is not None:
            
            # Filter the MultiIndex dataframe to include only rows where the "t" column values are in the T_window.
            # The index is multi-indexed by Prefix, a, R, TS_val, so we need to filter by the "t" column values.
            #data.groupby(level=0).filter(lambda x: x['t'].between(T_window[0], T_window[-1]).any())
            #data.groupby(level=0).filter(lambda x: (x["t"] >= T_window[0]) & (x["t"] <= T_window[-1])).any()
            #filtered_data = data.groupby(level=0).filter(lambda x: ((x["t"] >= T_window[0]) & (x["t"] <= T_window[-1])).any())
            #filtered_data =  data.loc[(slice(None), slice(None), slice(None), slice(None)), :].loc[(data["t"] >= T_window[0]) & (data["t"] <= T_window[1])]

            data_reset = data.reset_index()
            print(data_reset.index)
            print(data_reset.columns)
            # Filter rows based on the condition for column "t"
            filtered_data = data_reset[(data_reset["t"] >= T_window[0]) & (data_reset["t"] <= T_window[1])]
            # Next, restore the MultiIndex for the filtered data.
            filtered_data.set_index(["Prefix", "a", "R", "Tmax"], inplace = True)
            if filtered_data.empty:
                print(f"No data found for {Pre} in the T range {T_window}. Skipping....")
                continue
        
        # Next create a new column called "Hue" which will be the hue for the scatterplot.
        # If hue_label is not provided, we will set the hue to be the same for all points.
        # If hue_label is provided, we wil bin the values in the hue_label column into 10 bins and assign a hue to each bin.
        # If size_label is not provided, we will set the size to be the same for all points.

        if len(hue_label) > 0:
            # If hue_label is provided, we will bin the values in the hue_label column into 10 bins and assign a hue to each bin.
            # If the hue_label column is not found, we will set the hue to be the same for all points.
            if hue_label[0] in filtered_data.columns:
                # First find the maximum and minimum values in the hue_label column.
                #max_hue = filtered_data[hue_label[0]].max(); min_hue = filtered_data[hue_label[0]].min()
                # Drop all rows with NaN values in the hue_label column.
                kmeans = KMeans(n_clusters=10).fit(filtered_data[hue_label[0]].dropna().values.reshape(-1, 1))
                bin_centers = np.sort(kmeans.cluster_centers_.flatten())

                print(f"Clustered Bin Centers For {hue_label[0]} = {bin_centers}")
                bin_edges = np.sort([bin_centers[i] - (bin_centers[i+1] - bin_centers[i])/2 for i in range(len(bin_centers)-1)] + [bin_centers[-1] - (bin_centers[-1] - bin_centers[-2])/2] + [bin_centers[-1] + (bin_centers[-1] - bin_centers[-2])/2])
                # Next set all negative bin_edges to 0+ index(bin_edges).
                #print(f"Bin Edges For {hue_label[0]} = {bin_edges}")
                # Find minimum positive bin_edge (this is done to avoid negative bin_edges which are more likely to occur when the replicate is nearly flat.
                min_bin_edge = min([b for b in bin_edges if b > 0])
                bin_edges[np.where(bin_edges < 0)] = np.where(bin_edges < 0)/(1 + np.sum(np.where(bin_edges < 0))) * min_bin_edge
                print(f"Bin Edges For {hue_label[0]} = {bin_edges}")
                # If the hue_label column is found, we will bin the values into 10 bins and assign a hue to each bin.
                filtered_data["Hue"] = pan.cut(filtered_data[hue_label[0]], bins = bin_edges, labels = [bin_centers[i] for i in range(len(bin_centers))])
                #filtered_data["Hue"] = pan.cut(filtered_data[hue_label[0]], bins = 10, labels = [ min_hue + i*(max_hue - min_hue)/10 for i in range(10)])
            else:
                # If the hue_label column is not found, we will set the hue to be the same for all points.
                filtered_data["Hue"] = 1
        else:
            # If hue_label is not provided, we will set the hue to be the same for all points.
            filtered_data["Hue"] = 1

        # Next create a new column called "Size(size_label)" which will be the size for the individual points of the scatterplot.
        # If size_label is not provided, we will set the size to be the same for all points.
        # If size_label is provided, we will bin the values in the size_label column into 6 bins and assign a size to each bin.
        if len(size_label) > 0:
            if size_label[0] in filtered_data.columns:
                # First find the maximum and minimum values in the size_label column.
                kmeans = KMeans(n_clusters=6).fit(filtered_data[size_label[0]].values.reshape(-1, 1))
                bin_centers = np.sort(kmeans.cluster_centers_.flatten())
                bin_centers = np.array([50*round(b/50) for b in bin_centers])
                print(f"Clustered Bin Centers For {size_label[0]} = {bin_centers}")
                # Get the extent of the bin_centers.
                # Get the bins from the bin_centers.
                bin_edges = np.sort([bin_centers[i] - (bin_centers[i+1] - bin_centers[i])/2 for i in range(len(bin_centers)-1)] + [bin_centers[-1] - (bin_centers[-1] - bin_centers[-2])/2] + [bin_centers[-1] + (bin_centers[-1] - bin_centers[-2])/2])
                # Use the above bins to bin the values in the size_label column into 6 bins and assign a size to each bin.
                bin_edges[0] = 0 if bin_edges[0] < 0 else bin_edges[0]
                # Set non-unique bin_edges to be unique.
                bin_edges = np.unique(bin_edges)
                filtered_data[f"Size({size_label[0]})"] = pan.cut(filtered_data[size_label[0]], bins = bin_edges, labels = [bin_centers[i] for i in range(len(bin_edges)-1)])
            else:
                # If the size_label column is not found, we will set the size to be the same for all points.
                filtered_data[f"Size({size_label[0]})"] = 1
        else:
            # If size_label is not provided, we will set the size to be the same for all points.
            filtered_data[f"Size({size_label[0]})"] = 1

        # Next drop all columns that do not match the x_label, y_label, hue_label and size_label.
        filtered_data = filtered_data[[*x_label, *y_label, *hue_label, *size_label, "Hue", f"Size({size_label[0]})"]]
        # Next if two or more columns are identical, drop all but one of them.
        if filtered_data.columns.duplicated().any():
            print(f"WARNING: Found duplicated column named {filtered_data.columns[filtered_data.columns.duplicated()]} in filtered_data for {Pre}." 
                  "Dropping duplicates....")
            filtered_data = filtered_data.loc[:,~filtered_data.columns.duplicated()]

        # Finally, if reduce_points is set, reduce the number of points plotted for each trajectory by that factor.
        
        if reduce_points is not None:
            # First check if reduce_points is an integer (!= 0), if so, throw out everything but every reduce_points-th point.
            if isinstance(reduce_points, int) and reduce_points > 1:
                print(f"Reducing points by factor of {reduce_points} for {Pre}....")
                filtered_data = filtered_data.groupby(level=['Prefix', 'a', 'R', 'Tmax']).apply(lambda x: x.iloc[::reduce_points, :]).reset_index(level=[0,1,2,3], drop=True)
            elif reduce_points == 0 or reduce_points == "auto":
                # Automatically calculate a suitable reduction factor to reduce points to ~50 based 
                # on Kmeans binning of the "t" column.
                try:
                    print(f"Automatically reducing points for {Pre}....")
                    # First check if "t" column exists in data.
                    if "t" in filtered_data.columns:
                        # Check if number of unique points in "t" column is less than or equal to 50, if so, no reduction is needed.
                        if filtered_data["t"].nunique() <= 50:
                            print(f"Number of unique points in 't' column is less than or equal to 50 for {Pre}. No reduction needed.")
                        else:
                            # Then use Kmeans to cluster the "t" column into 50 clusters.
                            kmeans = KMeans(n_clusters=50).fit(filtered_data["t"].values.reshape(-1, 1))
                            bin_centers = np.sort(kmeans.cluster_centers_.flatten())
                            print(f"Clustered Auto t Centers = {bin_centers}")
                            # Next for each bin center, find the closest point in the "t" column and keep that point.
                            tvals_to_keep = []
                            for center in bin_centers:
                                tval = filtered_data["t"].iloc[(filtered_data["t"] - center).abs().argsort()[:1]].values[0]
                                tvals_to_keep.append(tval)                      # Add the first and last points to the indices to keep and then drop duplicates.
                            tvals_to_keep.append(filtered_data["t"].min()); tvals_to_keep.append(filtered_data["t"].max())
                            tvals_to_keep = list(set(tvals_to_keep))
                            print(f"Trying to keep t values: {tvals_to_keep} for {Pre}....")
                            filtered_data = filtered_data[filtered_data["t"].isin(tvals_to_keep)]
                            # Reset the index to remove gaps in the index after filtering.
                            #filtered_data = filtered_data.reset_index(level=[0,1,2,3], drop=True)
                            print(f"Automatically Reduced number of points to {filtered_data.shape[0]} for {Pre} with retained time values\n: {filtered_data['t'].unique()}")
                except Exception as e:
                    print(f"ERROR: {e} encountered while trying to automatically reduce points for {Pre}. Skipping point reduction....")
            else:
                print(f"WARNING: 't' column not found in data for {Pre}. Cannot automatically reduce points. Skipping point reduction....")

        # Finally, as a sanity check, report the number of rows for each MultiIndex group after point reduction.
        group_sizes = filtered_data.groupby(level=['Prefix', 'a', 'R', 'Tmax']).size()
        print(f"Number of points per trajectory after point reduction for {Pre}:\n{group_sizes}")



        print("====================================================================")
        print(f"Data for {Pre}: AFTER FILTERING")
        print("====================================================================\n")
        #print(data.info())
        print(filtered_data.head())
        print("____________________________________________________________________\n")
        print(filtered_data.tail())
        print("____________________________________________________________________\n")
        #Get data shape
        print(f"Shape = {filtered_data.shape}")
        print(filtered_data.columns)
        print("____________________________________________________________________\n")
        print(filtered_data.index)

        # Save the data to a csv file.
        savecsvdir = out_dir + f"{Pre}/Trajectories/"
        Path(savecsvdir).mkdir(parents=True, exist_ok=True)
        filtered_data.to_csv(savecsvdir + f"{Pre}_Trajectories_T_{T_window[0]}_{T_window[1]}.csv")
        # Otherwise no filtering is needed.
        # Get the list of a values and T values from the data.
        a_scaled_vals = filtered_data.index.get_level_values('a').unique().to_list()
        #Force a to be represented in decimal form rather than scientific notation (so 0.00001 is 0.00001 rather than 1e-5).
        tseries_vals = filtered_data.index.get_level_values('Tmax').unique().to_list()
        Rvals = filtered_data.index.get_level_values('R').unique().to_list()
        if minR not in Rvals:
            print(f"WARNING: minR = {minR} not found in R values of encountered data. Setting minR = -1.")
            minR = -1
        prelim_data_maxR = max(Rvals)
        Rmin = -1 if plot_mean_trajectory else min(min(Rvals)+1, minR)
        print(a_scaled_vals, tseries_vals, Rmin, prelim_data_maxR)
        for a_scaled in a_scaled_vals:
            #a = int(a) if float(a).is_integer() else a
            #Force a to be represented in decimal form rather than scientific notation.
            
            for TS in tseries_vals:
                if T_window is None:
                    savepngdir = out_dir + f"{Pre}/Trajectories/L_{g}/a_{a_scaled}/dP_{dP}/Geq_{Geq}/T_0_T_{TS}/"
                else:
                    savepngdir = out_dir + f"{Pre}/Trajectories/L_{g}/a_{a_scaled}/dP_{dP}/Geq_{Geq}/T_{T_window[0]}_T_{T_window[1]}/"
                Path(savepngdir).mkdir(parents=True, exist_ok=True)
                # Get maximum R value for given a and TS values, and iterate over all R values.
                Rmax = data.loc[(slice(None), a_scaled, slice(None), TS), :].index.get_level_values('R').max()
                if maxR != -1:
                    Rmax = min(Rmax, maxR)
                if Rmax < maxR:
                    print(f"WARNING: Encountered maxR = {Rmax} < Provided maxR = {maxR}. Setting maxR  = {Rmax}.")
                    Rmax = min(Rmax, maxR)
                fig_combined, axs_combined = plt.subplots(1, len(x_label), figsize = set_figsize(len(x_label)))
                print(f"Working on trajectories for {Pre} at a = {a_scaled}, T = {TS}, Rmax = {Rmax}....")
                for R in range(Rmin, Rmax +1):
                    # Get the data for the given a, TS and R values.
                    data_R = filtered_data.loc[(slice(None), a_scaled, R, TS), :]
                    if data_R.empty:
                        print(f"No data found for {Pre} at a = {a_scaled}, R = {R}, T = {TS}. Skipping....")
                        continue
                    # Plotting for a fixed a value, Tmax value species and Prefix for all R and t values.
                    fig, axs = plt.subplots(1, len(x_label), figsize = set_figsize(len(x_label)))
                    for s in range(len(x_label)):
                        # Plotting the trajectory for the given species.
                        ax = axs[s] if len(x_label) > 1 else axs
                        # Plotting the trajectory for the given species.
                        sea.scatterplot(data = data_R, x = x_label[s], y = y_label[s], hue = "Hue", size = f"Size({size_label[0]})",  
                                        sizes=(150,25),
                                        palette= hex_list[R % len(hex_list)], ax = ax, alpha = 0.6)
                        # Annotate each point with the time value.
                        # Set the y axis to log scale.
                        ax.set_yscale('log') #;ax.set_xscale('log')
                        print(f"Working on smart annotations at R= {R}...")
                        txts = []; sysout, syserr = redirect_todevnull();
                        for i in range(data_R.shape[0]):
                            # Get the time value for each point to 1 decimal place.
                            txt = f"{data_R['t'].iloc[i]:.1f}"
                            x = data_R[x_label[s]].iloc[i]; y = data_R[y_label[s]].iloc[i];
                            txts.append(ax.text(x, y, txt, fontsize = 'x-small', color = 'darkslategrey', alpha = 0.4))
                            #ax.annotate(txt, (data_R[x_label[s]].iloc[i], data_R[y_label[s]].iloc[i]), fontsize = 'x-small', color = 'grey', alpha = 0.25)
                        adjT.adjust_text(txts, arrowprops=dict(arrowstyle="->", color='darkslategrey', lw=0.5, alpha = 0.4))
                        sys.stdout = sysout; sys.stderr = syserr    
                        #ax.set_title(f"Individual Trajectory for {Pre} at a = {a_scaled}, R = {R}, T = {TS}")
                        #ax.set_xlabel(x_label[s]); ax.set_ylabel(y_label[s])
                        # Check if the x_label[s] and y_label[s] are in the dict_varlabels, and use the corresponding value if found.
                        ax.set_xlabel(dict_varlabels[x_label[s]]) if x_label[s] in dict_varlabels else x_label[s]
                        ax.set_ylabel(dict_varlabels[y_label[s]]) if y_label[s] in dict_varlabels else y_label[s]
                        
                        # Only show 10% of the legend items.
                        handles, labels = ax.get_legend_handles_labels()
                        
                        formatted_labels = [f"{float(lbl):.2f}" if lbl.replace('.', '', 1).isdigit() else lbl for lbl in labels]
                        
                        # Replace f"Size({size_label[0]})" with f"{size_label[0]}" in the formatted labels.
                        formatted_labels = [lbl.replace(f"Size({size_label[0]})", f"{size_label[0]}") for lbl in formatted_labels]
                        print(f"Formatted Labels = {formatted_labels}")
                        '''
                        #print(f"Handles = {handles}")
                        hue_labels = formatted_labels[: 12] # Ten bins for the hue values.
                        hue_handles = handles[: 12]
                        size_labels = formatted_labels[12:] # Remaining labels for the size values.
                        # Convert the size labels to float values.
                        size_labels = np.array([float(lbl) for lbl in size_labels])
                        size_handles = handles[12:]
                        #print(f"Hue Labels = {hue_labels},\n Hue Handles = {hue_handles}")
                        # Define custom size labels
                        desired_size_labels = np.linspace(data_R[size_label[s]].min(), data_R[size_label[s]].max(), 6, dtype=int).tolist()
                        print(f"Desired Size Labels = {desired_size_labels}")
                        desired_size_handles = [size_handles[np.argmin(np.abs(size_labels - d))] for d in desired_size_labels]
                        # Create size handles based on the size labels (these are matplotlib.lines.Line2D objects)
                        #size_handles = [mlines.Line2D([], [], marker='o', linestyle='',  markersize= int(3*size/size_labels[0]),  label=f"{s}", color='gray') for size in size_labels]
                        #print(f"Size Handles = {desired_size_handles}")

                        
                        handles = hue_handles + desired_size_handles; labels = hue_labels + desired_size_labels

                        print(f"Handles = {handles}, \n Labels = {labels}")
                        '''
                        #ax.legend(handles[::int(len(handles)/10)], formatted_labels[::int(len(formatted_labels)/10)], fontsize = 'x-small',
                        ax.legend(handles, formatted_labels, fontsize = 'x-small',
                                  loc = 'upper right', bbox_to_anchor=(1.1, 1.0), borderaxespad=0.)

                       
                        # Plotting the trajectory for the given species for all R values in the same plot.
                        ax_combined = axs_combined[s] if len(x_label) > 1 else axs_combined
                        sea.scatterplot(data = data_R, x = x_label[s], y = y_label[s], hue = "Hue", size = size_label[s], palette= hex_list[R % len(hex_list)], 
                                        ax = ax_combined, alpha = 0.75, legend= False)
                        #ax_combined.set_xlabel(x_label[s]); ax_combined.set_ylabel(y_label[s])
                        ax_combined.set_xlabel(dict_varlabels[x_label[s]]) if x_label[s] in dict_varlabels else x_label[s]
                        ax_combined.set_ylabel(dict_varlabels[y_label[s]]) if y_label[s] in dict_varlabels else y_label[s]
                        # Set the y axis to log scale.
                        #ax_combined.set_yscale('log') #ax_combined.set_xscale('log')
                        
                        # Get the current position of the axis
                        pos = ax_combined.get_position(); 
                        # Designing the colourbars for the hue values.
                        cbar_width = 0.015; cbar_gap = 0.003;
                        # To ensure cbars fit within the figure, set the left position of the cbar as follows
                        fixed_left_pos =  1 - (cbar_width + cbar_gap) * max(Rmax+1, len(hex_list))
                        # Creating colourbars for the hue values, only for the total number of distinct palettes available.
                        if R < len(hex_list):
                            sm = plt.cm.ScalarMappable(cmap=get_continuous_cmap(hex_list[R], float_list))
                            sm.set_array([])
                            # Manually define the colorbar position (normalize to figure size)
                            
                            left_pos_R = fixed_left_pos + R * (cbar_width + cbar_gap)  # Adjust left position for each colorbar
                            cbar_ax = fig_combined.add_axes([left_pos_R, pos.y0, cbar_width, pos.height])  # [left, bottom, width, height]
                            cbar = fig_combined.colorbar(sm, ax=ax_combined, orientation='vertical', cax=cbar_ax)

                            # Don't set labels for the colourbar.
                            cbar.ax.tick_params(size=0);cbar.set_ticks([]);
                            if R == Rmax or R == len(hex_list) - 1:
                                # For the last cbar, set the ticks and labels.
                                cbar.set_ticks([0, 1]); cbar.set_ticklabels([f"{float_list[0]:d}", f"{float_list[-1]:d}"])
                                # Finally, adjust the position of the axis to ensure all cbars within the figure.
                                axcomb_width = fixed_left_pos + cbar_width;
                                # Update the axis position to constrain it within the desired width
                                ax_combined.set_position([pos.x0, pos.y0, axcomb_width - pos.x0, pos.height])

                    
                    # End of s loop
                    if T_window is not None:
                        video_tvals = [T_window[0], T_window[1]]
                    else:
                        video_tvals = [0 , TS]
                    
                    if R == -1:
                        fig.suptitle(f"Mean Trajectory for {Pre} at a = {a_scaled}, R = {R}, T = [{video_tvals[0]}, {video_tvals[1]}]")
                    else:
                        fig.suptitle(f"Individual Trajectory for {Pre} at a = {a_scaled}, R = {R}, T = [{video_tvals[0]}, {video_tvals[1]}]")
                    
                    plt.savefig(savepngdir + f"IndivTrajectory_a_{a_scaled}_R_{R}_T_{video_tvals[-1]}.png")
                    #plt.show()
                    plt.close()
                # End of R loop
                if T_window is not None:
                    video_tvals = [T_window[0], T_window[1]]
                    fig_combined.suptitle(f"Combined Trajectory for {Pre} at a = {a_scaled}, R = {R}, T = [{T_window[0]}, {T_window[1]}]")
                else:
                    video_tvals = [0 , TS]
                    fig_combined.suptitle(f"Combined Trajectory for {Pre} at a = {a_scaled}, R = {R}, T = [0, {TS}]")
                
                # Plotting the combined trajectory for the given species for all R values in the same plot.
                plt.savefig(savepngdir + f"CombinedTrajectory_a_{a_scaled}_T_{video_tvals[-1]}.png")
                #plt.show()
                plt.close()

                # Use home_video function to create a video of the plots for each a value.
                video_tvals = [0] + tseries_vals if T_window is None else [T_window[0], T_window[1]]
                home_video(out_dir, out_dir, prefixes=[f"{Pre}/Trajectories"], a_vals= [a_scaled], T_vals=video_tvals, maxR= Rmax, minR= Rmin,
                          pngformat= "IndivTrajectory_a_{a}_R_{R}_T_{Tmax}.png", pngdir= "{Pre}/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T_{Tmin}_T_{Tmax}/",
                            videoname = f"IndivTrajectory_{Pre}_Tmin_{video_tvals[0]}_Tmax_{video_tvals[-1]}.mp4", video_relpath= "{Pre}/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/Videos/")
                
            # End of TS loop
            

        # End of a loop
        # Use home_video function to create a video of the individual trajectory plots for all a.
        #video_tvals = [0] + tseries_vals if T_window is None else [T_window[0], T_window[1]]
        #home_video(out_dir, out_dir, prefixes=[f"{Pre}/Trajectories"], a_vals= a_scaled_vals, T_vals=video_tvals, maxR= prelim_data_maxR, minR= Rmin,
        #           pngformat= "IndivTrajectory_a_{a}_R_{R}_T_{T}.png", pngdir= "{Pre}/L_{g}/a_{a}/dP_{dP}/Geq_{Geq}/T_{Tmin}_T_{Tmax}/", 
        #           videoname = f"IndivTrajectory_{Pre}_Tmax_{video_tvals[-1]}.mp4", video_relpath= "{Pre}/Videos/{amin}-{amax}/")


def analyse_FRAME_EQdata(indir, out_dir, prefixes=[], a_vals=[], T_vals =[], Tavg_window_index = [-100000, 0], filename = "MEAN_STD_Surviving_Runs.txt"):

    savefilename = re.sub(r'\.txt$', '', filename)
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data = pan.DataFrame()
    sea.set_palette("husl")
    colours = [hex_list[i][-1] for i in range(len(hex_list))]
    for Pre in prefixes:

        savedir = out_dir + f"{Pre}/PhaseDiagrams/"
        Path(savedir).mkdir(parents=True, exist_ok=True)
        
        data = get_FRAME_EQdata(indir, Pre, a_vals, [], filename)
        # Data has columns for each species, of the form:
        # AVG[P(x; t)]	VAR[P(x; t)]	COUNT[P(x; t)]	AVG[G(x; t)]	VAR[G(x; t)]	COUNT[G(x; t)] ...
        # We are interested in the AVG and VAR columns for each species.
        if data.empty:
            print(f"No data found for {Pre}. Skipping....")
            continue
        print(f"Data for {Pre}:")
        #print(data.info())
        print(data.head())
        print(data.tail())
        print(data.columns)
        print(data.index)
        #print(data.describe())
        # Save the data to a csv file
        data.to_csv(savedir + f"{savefilename}_dP_{dP}_Geq_{Geq}.csv")
        #Concatenate the data to the combined data
        combined_data = pan.concat([combined_data, data], axis = 0)
        Avals_list = sorted(data.index.get_level_values('a').unique().to_list())
        #TSvals_list = sorted(data.index.get_level_values('Tmax').unique().to_list())
        Tavg = pan.DataFrame()

        for a in Avals_list:
            # Get Tmax for this a value
            Tmax = data.loc[(slice(None), a, slice(None)), :].index.get_level_values('T').max()
            if(Tavg_window_index[0] < 0 and Tavg_window_index[1] <= 0):
                Twin_min = Tmax + Tavg_window_index[0]; Twin_max = Tmax + Tavg_window_index[1]
            else:
                Twin_min = Tavg_window_index[0]; Twin_max = Tavg_window_index[1]
            print(f"Averaging data for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}:")
            # Try averaging only if Twin_max is less than or equal to Tmax and Twin_min is greater than or equal to 0.
            if  Twin_min < 0:
                print(f"Skipping averaging for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}.")
                continue
            df = data.loc[(slice(None), a, slice(Twin_min, Twin_max)), :].groupby(['Prefix', 'a']).mean()
            if df.empty:
                print(f"No data found for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}. Skipping....")
                continue
            df.index = pan.MultiIndex.from_tuples([(Pre, a, Twin_max)], names = ['Prefix', 'a', 'Twindow'])
            Tavg = pan.concat([Tavg, df], axis = 0).sort_index(axis = 0)
        # End of a loop
        # Iterating over the Pre, a indices, average over the last Tavg_window_index values of T and assign to new DataFrame.


        print(f"Data for {Pre} after averaging over last {Tavg_window_index[1]} values of T:")
        print(Tavg.head())
        print(Tavg.tail())
        print(Tavg.columns)
        print(Tavg.index)

        #Save the data to a csv file
        Tavg.to_csv(savedir + f"DEBUG_Twin_{Twin_min}_{Twin_max}_dP_{dP}_Geq_{Geq}.csv")

        # Use Tavg to create plots of a vs. the mean of the data for each species, with STD as infill.
        fig, axs = plt.subplots(1, SPB, figsize = set_figsize(SPB))
        for s in range(SPB):
            ax = axs[s] if SPB > 1 else axs
            ax.scatter(Tavg.index.get_level_values('a'), Tavg.iloc[:, 3*s], label = r"$\mu_{%g}$" %(s), color = colours[s], s = 15, alpha = 0.75)
            err = 2 * np.sqrt(Tavg.iloc[:, 3*s + 1]) # Get 95% confidence interval
            ax.fill_between(Tavg.index.get_level_values('a'), Tavg.iloc[:, 3*s] - err, Tavg.iloc[:, 3*s] + err, alpha = 0.35, color = colours[s])
            # Annote the points with the values in data_T.iloc[:, 3*s+2] (which is the COUNT column of the data)
            for i, txt in enumerate(Tavg.iloc[:, 3*s + 2]):
                ax.annotate(int(txt), (Tavg.index.get_level_values('a')[i], Tavg.iloc[i, 3*s]))
            ax.set_title(r" $ \langle \langle \rho_{%g}(x, t = %g) \rangle_{x} \rangle_{r} $" % (s, Twin_max))
            ax.set_xlabel(r'R $(mm/hr)$')
            ax.set_ylabel(str(data.columns[3*s]))
            ax.legend()
        fig.suptitle(r"$\mu$ and $\sigma$ of" + f" Species For {Pre}, dP = {dP}, Geq = {Geq} B/W {Twin_min} -- {Twin_max}")
        plt.savefig(savedir + f"MeanSTD_{Pre}_amin_{Avals_list[0]}_amax_{Avals_list[-1]}_Twin_{Twin_min}_{Twin_max}_dP_{dP}_Geq_{Geq}.png")
        plt.show()
        plt.close()
        

        ''' # Previous version of the code, where we iterate over T values and create plots for each T value.
        # Use the data to create plots.
        # For each T-value, create a plot of a vs. the mean of the data for each species, with STD as infill.
        for T in T_vals:
            T = int(T) if float(T).is_integer() else T
            data_T = data.loc[(slice(None), slice(None), T), :]
            fig, axs = plt.subplots(1, SPB, figsize = set_figsize(SPB))
            for s in range(SPB):
                ax = axs[s] if SPB > 1 else axs
                ax.scatter(data_T.index.get_level_values('a'), data_T.iloc[:, 3*s], label = r"$\rho_{%g}$" %(s), color = colours[s], s = 15, alpha = 0.75)
                err = 2 * np.sqrt(data_T.iloc[:, 3*s + 1]) # Get 95% confidence interval
                ax.fill_between(data_T.index.get_level_values('a'), data_T.iloc[:, 3*s] - err, data_T.iloc[:, 3*s] + err, alpha = 0.35, color = colours[s])
                # Annote the points with the values in data_T.iloc[:, 3*s+2] (which is the COUNT column of the data)
                for i, txt in enumerate(data_T.iloc[:, 3*s + 2]):
                    ax.annotate(txt, (data_T.index.get_level_values('a')[i], data_T.iloc[i, 3*s]))
                ax.set_title(r" $ \langle \langle \rho_{%g}(x, t = %g) \rangle_{x} \rangle_{r} $" % (s, T))
                ax.set_xlabel(r'R $(mm/hr)$')
                ax.set_ylabel(str(data.columns[3*s]))
                ax.legend()
            fig.suptitle(r"$\mu$ and $\sigma$ of" + f" Species For {Pre} at T = {T}, dP = {dP}, Geq = {Geq}")
            plt.savefig(savedir + f"MeanSTD_{Pre}_T_{T}_dP_{dP}_Geq_{Geq}.png")
            plt.show()
            plt.close()
        '''

    savedir = out_dir + "PhaseDiagrams/"
    Path(savedir).mkdir(parents=True, exist_ok=True)
    # Finally if prefixes > 1, use the combined data to create plots of a vs. the mean of the data for each species, with STD as infill, and Prefix as hue.
    if len(prefixes) > 1:
        for T in T_vals:
            T = int(T) if float(T).is_integer() else T
            combined_data_T = combined_data.loc[(slice(None), slice(None), T), :]
            fig, axs = plt.subplots(1, SPB, figsize = set_figsize(SPB))
            for s in range(SPB):
                ax = axs[s] if SPB > 1 else axs
                sea.lineplot(data = combined_data_T, x = 'a', y = combined_data_T.iloc[:, 3*s], hue = 'Prefix', ax = ax)
                err = 2 * np.sqrt(combined_data_T.iloc[:, 3*s + 1])
                ax.fill_between(combined_data_T.index.get_level_values('a'), combined_data_T.iloc[:, 3*s] - err, combined_data_T.iloc[:, 3*s] + err, alpha = 0.35)
                # Annote the points with the values in data_T.iloc[:, 3*s+2] (which is the COUNT column of the data)
                #for i, txt in enumerate(combined_data.iloc[:, 3*s + 2]):
                #    ax.annotate(txt, (combined_data.index.get_level_values('a')[i], combined_data.iloc[i, 3*s]))
                ax.set_title(r" $ \langle \langle \rho_{%g}(x, t = %g) \rangle_{x} \rangle_{r}$" % (s, T))
                ax.set_xlabel('R $(mm/hr)$')
                ax.set_ylabel(str(data.columns[3*s]))
                ax.legend()
            fig.suptitle(r"$\mu$ and $\sigma$ of" + f" For All Prefixes at T = {T}, dP = {dP}, Geq = {Geq}")
            plt.savefig(savedir + f"MeanSTD_T_{T}_AllPrefixes.png")
            plt.show()
            plt.close()


def analyse_PRELIMS_EQdata(indir, out_dir, prefixes=[], a_vals=[], TS_vals =[], Tavg_window_index = [-100000, 0], serType = "TSERIES",
                        meanfilename = "Mean_{serType}_T_{TS}.csv",
                        var_labels= [ "<<P(x; t)>_x>_r" , "<<G(x; t)>_x>_r", "<<Pr(x; t)>_x>_r"], plt_violin=True, a_scaling =1):

    
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(indir, '*/'))]
    combined_data = pan.DataFrame()
    sea.set_palette("husl")
    colours = [hex_list[i][-1] for i in range(len(hex_list))]
    colours_med = [hex_list[i][-4] for i in range(len(hex_list))]
    colours_light = [hex_list[i][1] for i in range(len(hex_list))]
    #colours_med = gammapair_hex_list
    include_col_labels = ["t"] + var_labels

    if serType == "MOVSERIES":
        colours_med = gammapair_hex_list
        dict_var_labels = {"<GAM[G(x; t)]>_x": r'$\langle\gamma[G(x; t)]\rangle_x$' , "<GAM[Pr(x; t)]>_x": r'$\langle\gamma[Pr(x; t)]\rangle_x$'}
    elif serType == "TSERIES":
        dict_var_labels = {"<<P(x; t)>_x>_r": r'$\langle V(x, t) \rangle_{x}$' , "<<G(x; t)>_x>_r": r'$\langle G(x, t) \rangle_{x}$',
                            "<<Pr(x; t)>_x>_r": r'$\langle Pr(x, t) \rangle_{x}$'}
    for Pre in prefixes:

        savefilename = re.sub(r'\.csv$', '', meanfilename.format(Pre= Pre, g= g, dP= dP, Geq= Geq, TS= TS_vals[0], serType= serType))

        savedir = out_dir + f"{Pre}/PhaseDiagrams/"
        Path(savedir).mkdir(parents=True, exist_ok=True)
        
        data = get_PRELIM_EQdata(indir, Pre, a_vals, TS_vals, a_scaling, meanfilename= meanfilename, include_col_labels= include_col_labels, serType = serType)
        # Note that the number of columns can vary for each a and T value, as R_max can vary, so we need to handle this.
        # Data has columns for each species given by {var}, of the form:
        # t   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]   ...
        # where {var} is one of the species names. The dataframe is multi-indexed by Prefix, a, R, Tmax,
        # where R = -1 for the mean values.
        if data.empty:
            print(f"No data found for {Pre}. Skipping....")
            continue
        print(f"Data for {Pre}:")
        #print(data.info())
        print(data.head())
        print(data.tail())
        print(data.columns)
        print(data.index)
        #print(data.describe())
        # Save the data to a csv file
        data.to_csv(savedir + f"{savefilename}_dP_{dP}_Geq_{Geq}.csv")
        #Concatenate the data to the combined data
        combined_data = pan.concat([combined_data, data], axis = 0)
        Avals_scaled_list = sorted(data.index.get_level_values('a').unique().to_list())
        print(Avals_scaled_list)
        TSvals_list = sorted(data.index.get_level_values('Tmax').unique().to_list())
        Tavg = pan.DataFrame(); Tbox = pan.DataFrame()
        # Iterating over the Pre, a indices, average over the last Tavg_window_index values of T and assign to new DataFrame.
        for a in Avals_scaled_list:
            for TS in TSvals_list:
                # Get Tmax for this a value
                Tmax = TS
                if(Tavg_window_index[0] < 0 and Tavg_window_index[1] <= 0):
                    Twin_min = Tmax + Tavg_window_index[0]; Twin_max = Tmax + Tavg_window_index[1]
                else:
                    Twin_min = Tavg_window_index[0]; Twin_max = Tavg_window_index[1]
                print(f"Averaging data for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}:")
                # Try averaging only if Twin_max is less than or equal to Tmax and Twin_min is greater than or equal to 0.
                if  Twin_min < 0:
                    print(f"Skipping averaging for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}.")
                    continue
                # Average over the last Tavg_window_index values of t (for average surviving/all replicates) (R = -1)
                Twin_data = data[(data.index.get_level_values('a') == a) & (data.index.get_level_values('Tmax') == TS) 
                                 & ( data["t"] >= Twin_min) & (data["t"] <= Twin_max)]
                df = Twin_data.groupby(['Prefix', 'a', 'R', 'Tmax']).mean()
                if df.empty:
                    print(f"No data found for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}. Skipping....")
                    continue
                # Create a new Multi-index list of tuples with R = -1 for the mean values.
                index = []; index.extend([(Pre, a, i, Twin_min, Twin_max, TS) for i in sorted(df.index.get_level_values('R').unique())])
                df.index = pan.MultiIndex.from_tuples(index, names = ['Prefix', 'a', 'R', 'Twin_min', 'Twin_max', 'Tmax'])

                Tavg = pan.concat([Tavg, df], axis = 0).sort_index(axis = 0)
            # End of a loop

        print(f"Data for {Pre} after averaging over {Tavg_window_index[0]} --- {Tavg_window_index[1]} values of T:")
        print(Tavg.head())
        print(Tavg.tail())
        print(Tavg.columns)
        print(Tavg.index)

        #Save the data to a csv file
        Tavg.to_csv(savedir + f"DEBUG_{savefilename}_Twin_{Twin_min}_{Twin_max}_dP_{dP}_Geq_{Geq}.csv")

        '''If plt_violin is True, then create violin plots of the data for each species.
        This is by iterating over the Pre, a indices, grabbing col data for each species lying between Twin_min and Twin_max,
        and concatenating the data to a single Df, which is then used to create the violin plots.
        The data is then saved to a csv file.
        '''

        if( plt_violin):
            for a in Avals_scaled_list:
                for TS in TSvals_list:
                    Tmax = TS
                    if(Tavg_window_index[0] < 0 and Tavg_window_index[1] <= 0):
                        Twin_min = Tmax + Tavg_window_index[0]; Twin_max = Tmax + Tavg_window_index[1]
                    else:
                        Twin_min = Tavg_window_index[0]; Twin_max = Tavg_window_index[1]
                    print(f"Grabbing box/violin plot data for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}:")
                    if  Twin_min < 0:
                        print(f"Skipping averaging for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}.")
                        continue
                    Twin_data = data[(data.index.get_level_values('a') == a) & (data.index.get_level_values('Tmax') == TS)
                                     & ( data["t"] >= Twin_min) & (data["t"] <= Twin_max)]
                    if Twin_data.empty:
                        print(f"No data found for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max}. Skipping....")
                        continue
                    # Extend the Multi-index list of Twin_data to include Twin_min and Twin_max.
                    index =[]; index.extend([(Pre, str(a), i, Twin_min, Twin_max, TS) for Pre, a, i, _ in Twin_data.index])
                    Twin_data.index = pan.MultiIndex.from_tuples(index, names = ['Prefix', 'a', 'R', 'Twin_min', 'Twin_max', 'Tmax'])
                    # Concatenate the data to the Tbox DataFrame.
                    Tbox = pan.concat([Tbox, Twin_data], axis = 0).sort_index(axis = 0)
            # End of a loop

            #Concatenate all the entries in Tbox to itself, and ensure that the Multiindex for these entries have "$\sum(a)$"
            # as the entry for the a value.
            df_tmp = Tbox.copy()
            df_tmp.index = pan.MultiIndex.from_tuples([(Pre, r"$\sum(a)$", i, Twin_min, Twin_max, TS) 
                for Pre, a, i, Twin_min, Twin_max, TS in df_tmp.index], names = ['Prefix', 'a', 'R', 'Twin_min', 'Twin_max', 'Tmax'])
            
            Tbox = pan.concat([Tbox, df_tmp], axis = 0).sort_index(axis = 0)
            

            print(f"TBox data for {Pre} b/w {Twin_min} and {Twin_max}:")
            print(Tbox.head())
            print(Tbox.tail())
            print(Tbox.columns)
            print(Tbox.index)
            # Save the data to a csv file
            Tbox.to_csv(savedir + f"BOX-DEB_{savefilename}_Twin_{Twin_min}_{Twin_max}_dP_{dP}_Geq_{Geq}.csv")


        

        # Use Tavg to create plots of a vs. the mean of the data for each species, with STD as infill.
        fig, axs = plt.subplots(1, len(var_labels), figsize = set_figsize(len(var_labels)))
        for TS in TSvals_list:
            Tavg_TS = Tavg.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), TS), :]
            # For all entries where "AVG[<<P(x; t)>_x>_r]_SURV" is 0, set corresponding 
            # "AVG[<<G(x; t)>_x>_r]_ALL", "AVG[<<G(x; t)>_x>_r]_SURV",  "AVG[<<G(x; t)>_x>_r]", VAR[<<G(x; t)>_x>_r]_SURV, VAR[<<G(x; t)>_x>_r]_ALL, VAR[<<G(x; t)>_x>_r] to 0.
            #Tavg_TS.loc[Tavg_TS["AVG[" + var_labels[0] + "]_SURV"] == 0, ["AVG[" + var_labels[1] + "]_ALL", "AVG[" + var_labels[1] + "]_SURV", "AVG[" + var_labels[1] + "]", "VAR[" + var_labels[1] + "]_SURV", "VAR[" + var_labels[1] + "]_ALL", "VAR[" + var_labels[1] + "]"]] = 0

            # Set ALL ENTRIES of Tavg_TS with a < 0.0018 to 0.
            #Tavg_TS.loc[Tavg_TS.index.get_level_values('a') < 0.00178, :] = 0
            for s in range(len(var_labels)):

                # First plot scatter plots of Surviving and All replicates for each species.
                ax = axs[s] if SPB > 1 else axs
                # RECALL, _SURV and _ALL columns are given by multi_index for R = -1.
                Tavg_TS_all = Tavg_TS.loc[(slice(None), slice(None), -1, slice(None), slice(None), TS), :]
                # Plotting mean of surviving replicates ( with _SURV in the name)
                ax.scatter(Tavg_TS_all.index.get_level_values('a'), Tavg_TS_all["AVG[" + var_labels[s] + "]_SURV"], label = r"$\mu_{surv}($" +f"{dict_var_labels[var_labels[s]]}" + r"$)$", color = colours[s], s = 25, alpha = 0.9, marker = 's')
                # Annotate the points with percentage of surviving replicates.
                # This is given by the ratio of the "R_SURV[{var}]" column to the maximum value of R ( which is max value of R for the Prefix, a, Twindow, TS).
                for i, Rsurv in enumerate(Tavg_TS_all["R_SURV[" + var_labels[s] + "]"]):
                    # Get maximum value of R for the Prefix, a, Twindow, TS.
                    Rmax = Tavg_TS.loc[(slice(None), Tavg_TS_all.index.get_level_values('a')[i], slice(None), slice(None), slice(None), TS), :].index.get_level_values('R').max() + 1;
                    txt = f"{100*Rsurv/Rmax:.1f} %"

                    ax.annotate(txt, (Tavg_TS_all.index.get_level_values('a')[i], Tavg_TS_all["AVG[" + var_labels[s] + "]_SURV"][i]), fontsize = 'x-small', color = 'grey', alpha = 0.35)
                # Infill variance for surviving replicates
                err = 2*np.sqrt(Tavg_TS_all["VAR[" + var_labels[s] + "]_SURV"])
                ax.fill_between(Tavg_TS_all.index.get_level_values('a'), Tavg_TS_all["AVG[" + var_labels[s] + "]_SURV"] - err,
                                Tavg_TS_all["AVG[" + var_labels[s] + "]_SURV"] + err, color = colours[s], alpha = 0.25)
            
                # Plotting mean of all replicates ( with _ALL in the name)
                ax.scatter(Tavg_TS_all.index.get_level_values('a'), Tavg_TS_all["AVG[" + var_labels[s] + "]_ALL"], label = r"$\mu_{all}($" +f"{dict_var_labels[var_labels[s]]}" +r"$)$", color = colours_med[s], s = 30, alpha = 0.8, marker = 'D', facecolor = 'none')

                # Plotting mean of individual replicates (ignoring R = -1 entries)
                for R in Tavg_TS.index.get_level_values('R').unique():
                    if R == -1: 
                        continue # Skip R = -1 entries (which represent the mean values)
                    print(f"Plotting data for {Pre} at a = {a} for T in range {Twin_min} -- {Twin_max} for R = {R}")
                    Tavg_TS_R = Tavg_TS.loc[(slice(None), slice(None), R, slice(None), slice(None), TS), :]
                    ax.scatter(Tavg_TS_R.index.get_level_values('a'), Tavg_TS_R[var_labels[s]], color ="grey", s = 15, alpha = 0.25)
                # End of R loop
                ax.set_ylim(bottom = -0.05); #ax.set_xlim(right= 0.014)
                if serType == "MOVSERIES"  and var_labels[s] in dict_var_labels.keys():
                    ax.set_title(r" $ \langle$" +f"{dict_var_labels[var_labels[s]]}" + r" $\rangle_{t} $")
                    ax.set_ylabel(r" " +f"{dict_var_labels[var_labels[s]]}" + r" ", fontsize= 16)
                elif serType == "TSERIES" and var_labels[s] in dict_var_labels.keys():
                    ax.set_title(r" $ \langle$" +f"{dict_var_labels[var_labels[s]]}" + r" $\rangle_{t} $")
                    ax.set_ylabel(r"Avg Density, " +f"{dict_var_labels[var_labels[s]]}" + r"$\quad (kg\:km^{-2})$", fontsize=15)
                else:
                    ax.set_title(r" $ \langle \langle \rho_{%g}(x, t) \rangle_{x} \rangle_{t} $" % (s))
                    ax.set_ylabel(r" $ \langle \langle \rho_{%g}(x, t) \rangle_{x} \rangle_{t} $" % (s))
                #ax.set_title(r" $ \langle \langle \rho_{%g}(x, t) \rangle_{x} \rangle_{t} $" % (s))
                #ax.set_ylabel(r" $ \langle \langle \rho_{%g}(x, t) \rangle_{x} \rangle_{t} $" % (s))
                ax.set_xlabel(r'Rainfall $R \: (mm/hr)$', fontsize=15) # For Rietkerk
                #ax.set_xlabel(r'a $(hr^{-1})$') # For DP
                # Set x-axis to log scale (for DP)
                #ax.set_xscale('log')
                # Set x-axis ticks to display four significant figures (for DP)
                #ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.4g}'))
                
                ax.legend(fontsize=13)

                # Add vertical line at a = 0.0018
                #ax.axvline(x = 0.0018, color = 'red', linestyle = '--', alpha = 0.5)

                # Add a red dot for the mean of all replicates (with _ALL in the name) at a = 0.000201 and a = 0.0025.
                #ax.scatter([0.000201, 0.0025], [Tavg_TS_all.loc[(slice(None), 0.000201, -1, slice(None), slice(None), TS), "AVG[" + var_labels[s] + "]_ALL"], 
                #            Tavg_TS_all.loc[(slice(None), 0.0025, -1, slice(None), slice(None), TS), "AVG[" + var_labels[s] + "]_ALL"]], color = 'gold', s = 15, alpha = 0.9) 
            # End of s loop
            Twin_min = Tavg_TS_all.index.get_level_values("Twin_min").min(); Twin_max = Tavg_TS_all.index.get_level_values("Twin_max").max()
            fig.suptitle(r"$\mu$ and $\sigma$ of" + f" Species For {Pre}, dP = {dP}, Geq = {Geq} B/W {Twin_min} -- {Twin_max}")
            plt.savefig(savedir + f"MeanSTD_{serType}_amin_{Avals_scaled_list[0]}_amax_{Avals_scaled_list[-1]}_Twin_{Twin_min}_{Twin_max}_dP_{dP}_Geq_{Geq}.png", dpi =600, transparent=True)
            plt.show()
            plt.close()
        # End of TS loop

        if not Tbox.empty:

            # Drop rows where col entries for cols "var_labels[0]" and "var_labels[1]" are BOTH EITHER 1 or 0.
            if serType == "MOVSERIES":
                #pass
                #Tbox = Tbox[~( ((Tbox[var_labels[0]] == 0) | (Tbox[var_labels[0]] > 0.9)) & ((Tbox[var_labels[1]] == 0) | (Tbox[var_labels[1]] == 1)) )]
                #Tbox = Tbox[~((Tbox[var_labels[0]] > 0.8))]; Tbox = Tbox[~((Tbox[var_labels[1]] == 0) | (Tbox[var_labels[1]] == 1))]
                #Tbox = Tbox[~((Tbox[var_labels[1]] == 0))]
                Tbox = Tbox[~((Tbox[var_labels[0]] == 0) | (Tbox[var_labels[1]] == 0))] #For Rietkerk
                # For DP, for all rows corresponding to a = 0.0017, set the col entries for col "var_labels[1]" to 0.
                #Tbox.loc[(Tbox.index.get_level_values('a') == "0.0017"), var_labels[1]] = 0
            elif serType == "TSERIES":
                # Drop rows where col entries for col "var_labels[0]" (vegetation) is 0.
                Tbox = Tbox[~(Tbox[var_labels[0]] == 0)]
                # Also drop rows where col entries for either "var_labels[1]" (grazer) or "var_labels[2]" (predator) is 0.
                Tbox = Tbox[~((Tbox[var_labels[1]] == 0) | (Tbox[var_labels[2]] == 0))]
            # Plotting the box/violin plots for each species.
            # NOTE: Seaborn's violinplot function does NOT support NON-UNIQUE multi-indexes, so we need to reset the index.
            Tbox_flat = Tbox.reset_index()
            fig, axs = plt.subplots(1, len(var_labels), figsize = set_figsize(len(var_labels)))
            # Set the category order to ensure $\sum(a)$ is the first category.
            
            # Replace all Tbox_flat['a'] numeric strings with their numeric 5 significant figure representation
            #Tbox_flat['a'] = Tbox_flat['a'].replace({r"\$\\sum\(a\)\$": "$\sum(a)$"}, regex=True)
            Tbox_flat['a'] = Tbox_flat['a'].apply(lambda x: f"{float(x):.5g}" if x != "$\sum(a)$" else x)
            # Print unique values of Tbox_flat['a]
            print(f"Unique values of Tbox_flat['a']: {Tbox_flat['a'].unique()}")

            #cat_var_plot_order = ["$\sum(a)$"] + [f"{float(x):.5g}" for x in Tbox.index.get_level_values('a').unique() if x != "$\sum(a)$"]
            cat_var_plot_order = sorted([f"{float(x):.5g}" for x in Tbox.index.get_level_values('a').unique() if x != "$\sum(a)$"])
            cat_dtype = pan.api.types.CategoricalDtype(categories=cat_var_plot_order, ordered=True)
            Tbox_flat['a'] = Tbox_flat['a'].astype(cat_dtype)
            print(f"Category order for Tbox_flat['a']: {cat_var_plot_order}")
            for s in range(len(var_labels)):
                ax = axs[s] if SPB > 1 else axs
                # First plotting violin plot for '$\sum(a)$' and then violin plot for the rest of the data.
                #sea.violinplot(data = Tbox_flat[(Tbox_flat['a'] == "$\sum(a)$")], x = 'a', y = var_labels[s], 
                #               ax = ax, order = cat_var_plot_order, color = colours[s% len(colours)], cut=0, bw_adjust=0.5, bw_method='scott',
                #               inner='box', linewidth = 1.5, alpha = 0.85, density_norm='area')
                
                # Next plotting the violin plot for the rest of the data, excluding '$\sum(a)$', and all zero entries.
                #if var_labels[s] == "<GAM[G(x; t)]>_x":
                #    sea.violinplot(data = Tbox_flat[(Tbox_flat['a'] != "$\sum(a)$")], x = 'a', y = var_labels[s], 
                #               ax = ax, order = cat_var_plot_order, color=colours_med[s % len(colours_med)], cut=0, bw_adjust=0.5, bw_method='scott',
                #               inner='box', linewidth = 1.5, alpha = 0.85, density_norm='area')
                #               #, inner_kws=dict(alpha = 0.85, linewidth=1.1, box_width= 2.75))
                #elif var_labels[s] == "<GAM[Pr(x; t)]>_x":
                #    sea.violinplot(data = Tbox_flat[(Tbox_flat['a'] != "$\sum(a)$")], x = 'a', y = var_labels[s], 
                #               ax = ax, order = cat_var_plot_order, color=colours_med[s % len(colours_med)], cut=0, bw_adjust=0.5, bw_method='scott',
                #               inner='box', linewidth = 1.5, alpha = 0.85, density_norm='area') #, inner_kws=dict(marker="o", color="#6f46a6"))
                #               #, inner_kws=dict(alpha = 0.85, linewidth=1.1, box_width= 2.75, marker="o", color="#6f46a6"))
                #else:
                sea.violinplot(data = Tbox_flat[(Tbox_flat['a'] != "$\sum(a)$")], x = 'a', y = var_labels[s], 
                               ax = ax, order = cat_var_plot_order, color=colours_med[s % len(colours_med)], cut=0, bw_adjust=0.5, bw_method='scott',
                               inner='box', linewidth = 1.5, alpha = 0.85, density_norm='area')
                               # , inner_kws=dict(alpha = 0.85, linewidth=1.1, box_width= 2.75))
                # Add vertical line seperating sum data from individual a data.
                #ax.axvline(x=0.5, color='grey', linestyle='--', linewidth=1.5, alpha=0.75)
               
                # Set x-axis ticks to display 5 significant figures if numeric for DP

                '''
                if serType == "MOVSERIES":
                    #ax.set_xticks([0.5] + [i for i in Tbox.index.get_level_values('a').unique() if i != "$\sum(a)$"])
                    #ax.set_xticklabels([r"$\sum(a)$"] + [f"{i:.5g}" for i in Tbox.index.get_level_values('a').unique() if i != "$\sum(a)$"])
                else:
                    #ax.set_xticks(np.arange(len(Tbox.index.get_level_values('a').unique())))
                    #ax.set_xticklabels([r"$\sum(a)$"] + [f"{i:.5g}" for float(i) in Tbox.index.get_level_values('a').unique() if i != "$\sum(a)$" else pass])
                #'''              
                #ax.set_ylim(bottom = -0.05); ax.set_xlim(right= 0.014)
                if serType == "MOVSERIES" and var_labels[s] in dict_var_labels.keys():
                    #ax.set_title(r" $ P[$ " +f"{dict_var_labels[var_labels[s]]}" + r"$ ]$")
                    #ax.set_ylabel(r" $ P[$ " +f"{dict_var_labels[var_labels[s]]}" + r"$ ]$", fontsize=16)
                    if var_labels[s] == "<GAM[G(x; t)]>_x":
                        # In this case have the y-axis ticks appear on the left side of the plot using twin axis, with 
                        # colours_med[s % len(colours_med)] as the tick colour.
                        #ax.set_ylim(0.0, 0.6) #FOR DP
                        ax.set_ylim(0.0, 0.2) # For Rietkerk
                        ax_left = ax.twinx()
                        ax_left.set_ylim(ax.get_ylim())      # Mirror limits (optional)
                        ax_left.yaxis.set_label_position("left"); ax_left.yaxis.tick_left();
                        ax_left.tick_params(axis='y', labelsize=14, colors=colours_med[s % len(colours_med)])
                        #ax_left.tick_params(axis='x', labelsize=14); 
                        ax_left.set_ylabel(r" $ P[$ " +f"{dict_var_labels[var_labels[s]]}" + r"$ ]$", fontsize=16, color=colours_med[s % len(colours_med)])
                        #ax_left.set_ylim(bottom=0.0)
                        #ax_left.set_ylim(0.0, 0.6)

                    if var_labels[s] == "<GAM[Pr(x; t)]>_x":
                        # In this case have the y-axis ticks appear on the right side of the plot using twin axis, with 
                        # colours_med[s % len(colours_med)] as the tick colour.
                        #ax.set_ylim(0.5, 1.0) #FOR DP
                        ax.set_ylim(0.2, 1.0) # For Rietkerk
                        ax_right = ax.twinx()
                        ax_right.set_ylim(ax.get_ylim())      # Mirror limits (optional)
                        ax_right.yaxis.set_label_position("right"); ax_right.yaxis.tick_right();
                        ax_right.tick_params(axis='y', labelsize=14, colors=colours_med[s % len(colours_med)])
                        #ax_right.tick_params(axis='x', labelsize=14); 
                        ax_right.set_ylabel(r" $ P[$ " +f"{dict_var_labels[var_labels[s]]}" + r"$ ]$", fontsize=16, color=colours_med[s % len(colours_med)])
                        # For DP, set y-axis limits between 0.5 and 1.0
                        #ax_right.set_ylim(top=1.0)
                        #ax_right.set_ylim(0.6, 1.0); 
                    # Finally remove all y-axis ticks, labels for clarity (For Dual plots)
                    ax.yaxis.set_ticks_position('none'); ax.yaxis.set_ticklabels([]); ax.set_ylabel('');
                elif serType == "TSERIES" and var_labels[s] in dict_var_labels.keys():
                    ax.set_title(r" $ P[$ " +f"{dict_var_labels[var_labels[s]]}" + r"$ ]$")
                    ax.set_ylabel(r"Avg Density, " +f"{dict_var_labels[var_labels[s]]}" + r"$\quad (kg\:km^{-2})$", fontsize=15)
                else:
                    ax.set_title(r" $ P[ \rho_{%g}(x, t) ] $" % (s))
                    ax.set_ylabel(r" $ P[ \rho_{%g}(x, t) ] $" % (s))
                #ax.set_title(r" $ \langle \langle \rho_{%g}(x, t) \rangle_{x} \rangle_{t} $" % (s))
                #ax.set_ylabel(r" $ \langle \langle \rho_{%g}(x, t) \rangle_{x} \rangle_{t} $" % (s))
                ax.set_xlabel(r'Rainfall $R \: (mm/hr)$', fontsize=16) # For Rietkerk
                #ax.set_xlabel(r'Growth Rate $a \: (hr^{-1})$', fontsize=16) # For DP
                ax.tick_params(axis='x', labelsize=14); #ax.tick_params(axis='y', labelsize=14)
                ax.set_xticklabels(cat_var_plot_order, rotation=30, fontsize=14, ha='right')

                # If cat_var_plot_order has > 10 entries, only display every second tick label and rotate by 45 degrees.
                if len(cat_var_plot_order) > 8:
                    ax.set_xticks(range(len(cat_var_plot_order)))
                    ax.set_xticklabels(cat_var_plot_order, rotation=45, fontsize=14, ha='right')
                    ax.tick_params(axis='y', labelsize=14)
                    #for i, label in enumerate(ax.get_xticklabels()):
                    #    if i % 2 == 0:
                    #        label.set_visible(False)
                #ax.set_xticks(range(len(cat_var_plot_order)))
                #ax.set_xticklabels(cat_var_plot_order)

                # Display current x-axis tick labels
                tick_labels = ax.get_xticklabels(); print(f"Current x-axis tick labels: {tick_labels}")
            
            fig.suptitle(f"Eq Distribution of Species For {Pre}, dP = {dP}, Geq = {Geq} B/W {Twin_min} -- {Twin_max}")
            plt.tight_layout()
            plt.savefig(savedir + f"ALT3_VIOLIN_{serType}_amin_{Avals_scaled_list[0]}_amax_{Avals_scaled_list[-1]}_Twin_{Twin_min}_{Twin_max}_dP_{dP}_Geq_{Geq}.png", dpi =600, transparent=True)
            plt.show(); plt.close()

    # End of Pre loop

''' # Summary Description of multiplot_PRELIMS_CompareEQ(...)

This function essentially creates combined plots of multiple dataframes generated and saved to csv files by the 
analyse_PRELIMS_EQdata(...) function. NOTE: The analyse_PRELIMS_EQdata(...) function saves the data in the following format:
Savedir = out_dir + f"{Pre}/PhaseDiagrams/" where {Pre} is the Prefix of the data.
Filename = f"DEBUG_{savefilename}_L_{g}_dP_{dP}_Geq_{Geq}.csv"
where {savefilename} is the filename of the mean data file, and {g}, {dP}, {Geq} are the values of the parameters L, dP, Geq respectively.
The data is a Multi-Index Pandas file saved in the format:
Prefix, a, R, Tmax, Twin_min, Twin_max, AVG[{var}]_ALL, ..., VAR[{var}]_ALL, ...,  
    AVG[{var}]_SURV, ..., VAR[{var}]_SURV, ..., R_SURV[{var}], ..., t, {var} ...
where {var} are the species names provided in the var_labels list in the analyse_PRELIMS_EQdata(...) function.

Prefix, a, R, Tmax, Twin_min, Twin_max are the indices of the data, where R = -1 for the mean values.

Given the data saved in this format, the multiplot_PRELIMS_CompareEQ(...) function creates combined plots of the data for each species, 
doing so in a similar fashion to the analyse_PRELIMS_EQdata(...) function, but for multiple {Pre}, {g}, {dP} OR {Geq} values.

multiplot_PRELIMS_CompareEQ(...) takes the following arguments:
indir: The input directory where the data is saved.
out_dir: The output directory where the plots will be saved.
compare_list: A dictionary of the form {"Prefix": [Pre_string_vals], "g": [g_vals], "Prefix": [dP_vals], "Prefix": [Geq_vals]} 
where the values are lists of the parameters Pre, g, dP, Geq respectively.
meanfilename: The filename of the mean data file, which is used to load the data. 
If provided as "guess", the function will attempt to guess the filename based on comparisons to matching files in the directory as 
DEBUG_*_{compare_list.keys()}_{compare_list.values()}.csv. for each key, value pair in compare_list (other than "Prefix").

var_labels: A list of the species names in the data, which are used to create the plots.

It will then create distinct combined plots of the data for each compare_list key, iterating over the values in the compare_list values.
Plots will be generated in an identical manner to the analyse_PRELIMS_EQdata(...) function.

It will do so by reading the DEBUG_*_{compare_list.keys()}_*.csv files in the indir directory, for each key in compare_list.keys().
It will read these as Pandas Multi-Index DataFrames, and then concatenate them to create a combined DataFrame for each key,
wuth a new Multi-Index level for the key. It will then create plots of the data in the same manner as the analyse_PRELIMS_EQdata(...) function,
iterating over the values in the compare_list values for each key.

For each key, the combined plot will be saved to the out_dir directory in the format:
f"COMPARE_{compare_list.keys()}_MeanSTD_{key}_amin_{Avals_list[0]}_amax_{Avals_list[-1]}_Twin_{Twin_min}_{Twin_max}.png"

'''

def get_combined_df(combined_data, indir, compare_key, compare_vals, prefixes, init_meanfilename = "guess",
                                 var_labels= [ "<<P(x; t)>_x>_r" , "<<G(x; t)>_x>_r", "<<Pr(x; t)>_x>_r"]):

    # Try to read the data from the matching file in the indir directory.
    # If the file does not exist, skip to the next value in the compare_vals list.
    # If there are multiple files that match the pattern, provide an error message and skip to the next value.
    
    meanfilename = init_meanfilename # Dummy variable to store the meanfilename for each iteration.
    if compare_key == "Prefixes":
        for Pre in compare_vals:
            if  init_meanfilename == "guess":
                meanfilename = f"DEBUG_MEAN_TSERIES*.csv"
            meanfilename = meanfilename.format(Pre=Pre)#, val=val)
            try:
                files = glob.glob(indir.format(Pre=Pre) + meanfilename)
            except FileNotFoundError:
                print(f"File {indir.format(Pre=Pre) + meanfilename} not found. Skipping...")
                continue

            if len(files) > 1:
                    print(f"Multiple files found for {indir.format(Pre=Pre) + meanfilename}. Skipping...")
                    continue
            try:
                print(f"Reading file {files[0]}... with meanfilename {meanfilename}")
                data = pan.read_csv(files[0], index_col = [0, 1, 2, 3, 4, 5])
                # The data is a Multi-Index with the following index levels: Prefix, a, R, Twin_min, Twin_max, Tmax
            except FileNotFoundError:
                    print(f"File {indir.format(Pre=Pre)+ meanfilename} not found. Skipping...")
                    continue
            except pan.errors.ParserError:
                    print(f"Error reading file {indir.format(Pre=Pre) + meanfilename}. Skipping...")
                    continue
            data[compare_key] = Pre
            data.set_index(compare_key, append = True, inplace = True) # Add the compare_key as a new index level.

            combined_data = pan.concat([combined_data, data], axis = 0)
    else:
        for val in compare_vals:
            for Pre in prefixes:
                if init_meanfilename == "guess":
                    meanfilename = f"DEBUG_*{compare_key}_{val}_*csv"
                meanfilename = meanfilename.format(Pre=Pre, val=val)

                print(f"Working on {compare_key} comparison with value {val} for {Pre}...")
                print(f"Looking for file {indir.format(Pre=Pre) + meanfilename}...")
                try:
                    files = glob.glob(indir.format(Pre=Pre) + meanfilename)
                    #print(files)
                except FileNotFoundError:
                    print(f"File {indir.format(Pre=Pre) + meanfilename} not found. Skipping...")
                    continue
                if len(files) > 1:
                    print(f"Multiple files found for {indir.format(Pre=Pre) + meanfilename}. Skipping...")
                    continue
                try:
                    print(f"Reading file {files[0]}... with meanfilename {meanfilename}")
                    data = pan.read_csv(files[0], index_col = [0, 1, 2, 3, 4, 5])
                    # The data is a Multi-Index with the following index levels: Prefix, a, R, Twin_min, Twin_max, Tmax
                except FileNotFoundError:
                    print(f"File {indir.format(Pre=Pre)+ meanfilename} not found. Skipping...")
                    continue
                except pan.errors.ParserError:
                    print(f"Error reading file {indir.format(Pre=Pre) + meanfilename}. Skipping...")
                    continue
                # The data is a Multi-Index with the following index levels: Prefix, a, R, Twin_min, Twin_max, Tmax
                # First, add the compare_key as a new Multi-Index level.
                data[compare_key] = val
                data.set_index(compare_key, append = True, inplace = True) # Add the compare_key as a new index level.

                print(f"Data for {compare_key} comparison with value {val} for {Pre}:")
                print(data.head())
                print(data.tail())
                print(data.columns)

                combined_data = pan.concat([combined_data, data], axis = 0)
                print(f'Length of combined_data: {len(combined_data)}')


    return combined_data

def multiplot_PRELIMS_CompareEQ(indir, out_dir, compare_list, prefixes =[], meanfilename = "guess", plt_violin = True,
                                 var_labels= [ "<<P(x; t)>_x>_r" , "<<G(x; t)>_x>_r", "<<Pr(x; t)>_x>_r"]):
    
    colours = [hex_list[i][abs(6-i)] for i in range(len(hex_list))]# for j in range(0, len(hex_list[i]), 2)]
    #colours = [hex_list[i][-1] for i in range(len(hex_list))]
    colours_med = [hex_list[i][2] for i in range(len(hex_list))]
    colours_light = [hex_list[i][1] for i in range(len(hex_list))]
    
    # Get the list of keys and values from the compare_list dictionary.
    keys = list(compare_list.keys())
    # Iterate over each key in the compare_list dictionary, creating combined plots of the data for each key.
    # If "Prefix" is a key and prefixes is empty, set the prefixes list to the values for the key.
    meanfilename_init = meanfilename

    access_dir = indir + "{Pre}/PhaseDiagrams/"

    if "Prefixes" in keys and len(prefixes) == 0:
        prefixes = compare_list["Prefixes"]

    if "Prefixes" not in keys and len(prefixes) ==0:
        print("Error: No Prefix values provided in compare_list or prefixes list. Please provide a list of Prefix values.")
        return

    for compare_key in keys:
        # Get the values for the key from the compare_list dictionary.
        compare_vals = compare_list[compare_key]
        # If the key is "Prefixes", set the prefixes list to the compare_vals list.
        if len(compare_vals) == 0:
            print(f"No values provided for {compare_key}. Skipping...")
            continue
        # Create a combined DataFrame for the key, iterating over the values in the compare_vals list.
        combined_data = pan.DataFrame()
        # Get the combined DataFrame for the key.
        combined_data = get_combined_df(combined_data, access_dir, compare_key, compare_vals, prefixes, meanfilename, var_labels)
        # Verify that the combined data is not empty.
        if combined_data.empty:
            print(f"No data found for {compare_key} comparison with values {compare_vals}. Skipping...")
            continue

        # If there are multiple unique values for the "Twin_min" and "Twin_max" index levels, throw an error and skip to the next key.
        if len(combined_data.index.get_level_values('Twin_min').unique()) > 1 or len(combined_data.index.get_level_values('Twin_max').unique()) > 1:
            print(f"Multiple inconsistent values detected for Twin_min or Twin_max {compare_key} comparison with values {compare_vals}. Skipping...")
            continue

        # Show the combined data.
        print(f"Combined data for {compare_key} comparison with values {compare_vals}:")
        print(combined_data.head())
        print(combined_data.tail())
        print(combined_data.columns)
        print(combined_data.index)  

        ''' NOTE: The combined data is a Multi-Index DataFrame with the following index levels:
        Prefix, a, R, Twin_min, Twin_max, Tmax, {compare_key}  
        where R = -1 for the mean values.''' 

        
        
        # Next, create plots of the data for each species, iterating over the values in the compare_vals list as hues.
        # The plots will be saved to the out_dir directory.
        # compare_vals are stored in the {compare_key} index level of the combined data.
        # If "Prefix" column has multiple unique values, create a combined plot of the data for each species, with Prefix and ({compare_key} if applicable) as distinct hues.
        # If "Prefix" column has only one unique value, create a combined plot of the data for each species, with {compare_key} as the hue.

        savedir = access_dir.format(Pre = prefixes[0]) + f"CombinedPhaseDiagrams/"
        if compare_key == "Prefixes":
            savedir = out_dir + f"CombinedPhaseDiagrams/"
        Path(savedir).mkdir(parents=True, exist_ok=True)

        # Save the combined data to a csv file.
        combined_data.to_csv(savedir + f"Combined_{compare_key}_MeanSTD_{compare_vals}.csv")
        
        # Get the unique values for the Prefix index level.
        Prefix_vals = combined_data.index.get_level_values('Prefix').unique().to_list()
        # Get the unique values for the {compare_key} index level.
        compare_vals = combined_data.index.get_level_values(compare_key).unique().to_list()
        # Get the unique values for the a index level.
        Avals_list = combined_data.index.get_level_values('a').unique().to_list()

        TSvals_list = sorted(combined_data.index.get_level_values('Tmax').unique().to_list())

        for TS in TSvals_list:
            fig, axs = plt.subplots(1, len(var_labels), figsize = set_figsize(len(var_labels)))
            
            for val in compare_vals:
                Tavg_TSval = combined_data.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), TS, val), :]
                for s in range(len(var_labels)):
                    # First plot scatter plots of Surviving and All replicates for each species.
                    ax = axs[s] if len(var_labels) > 1 else axs
                    # RECALL, _SURV and _ALL columns are given by multi_index for R = -1.
                    Tavg_TSval_all = Tavg_TSval.loc[(slice(None), slice(None), -1, slice(None), slice(None), TS, slice(None)), :]
                    # Plotting mean of surviving replicates ( with _SURV in the name)
                    for Pre in prefixes:
                        print(f"Plotting data for {Pre} for {compare_key} = {val} at T = {TS}")
                        #if val in Prefix_vals:
                            # Represents case where compare_key is "Prefix", h
                        #    Pre = val
                        if compare_vals == Prefix_vals:
                            # Represents case where compare_key is "Prefix", hence val is the Prefix.
                            # Here, we only plot data in the case where Pre == val.
                            if Pre != val:
                                continue
                        Tavg_TSPreval_all = Tavg_TSval_all.loc[( Pre, slice(None), slice(None), slice(None), slice(None), TS, val), :]
                        ax.scatter(Tavg_TSPreval_all.index.get_level_values('a'), Tavg_TSPreval_all["AVG[" + var_labels[s] + "]_SURV"], label = r"$\mu_{{surv}}(\rho_{%g})$" %(s) + f" {compare_key} = {val}", 
                                   color = colours[1*len(compare_vals) + compare_vals.index(val)], s = 15, alpha = 0.9, marker = 's')
                        
                        # This is given by the ratio of the "R_SURV[{var}]" column to the maximum value of R ( which is max value of R for the Prefix, a, Twindow, TS).
                        # Infill variance for surviving replicates
                        err = 1*np.sqrt(Tavg_TSPreval_all["VAR[" + var_labels[s] + "]_SURV"])
                        ax.fill_between(Tavg_TSPreval_all.index.get_level_values('a'), Tavg_TSPreval_all["AVG[" + var_labels[s] + "]_SURV"] - err,
                                Tavg_TSPreval_all["AVG[" + var_labels[s] + "]_SURV"] + err, color = colours_med[1*len(compare_vals) + compare_vals.index(val)], alpha = 0.25)
            
                        # Plotting mean of all replicates ( with _ALL in the name)
                        #ax.scatter(Tavg_TSPreval_all.index.get_level_values('a'), Tavg_TSPreval_all["AVG[" + var_labels[s] + "]_ALL"], label = r"$\mu_{all}(\rho_{%g})$" %(s) + f" {compare_key} = {val}", 
                        #           color = colours[s*len(compare_vals) + compare_vals.index(val)], s = 15, alpha = 0.75, marker = 'D', facecolor = 'none')

                        # Plotting mean of individual replicates (ignoring R = -1 entries)
                        for R in Tavg_TSval.index.get_level_values('R').unique():
                            if R == -1: 
                                continue # Skip R = -1 entries (which represent the mean values)
                            Tavg_TSPreval_R = Tavg_TSval.loc[( Pre, slice(None), R, slice(None), slice(None), TS, val), :]
                            ax.scatter(Tavg_TSPreval_R.index.get_level_values('a'), Tavg_TSPreval_R[var_labels[s]], color ="grey", s = 15, alpha = 0.1)
                        # End of R loop
                    # End of Pre loop
                    ax.set_title(r" $ \langle \langle \rho_{%g}(x, t) \rangle_{x} \rangle_{t} $" % (s))
                    ax.set_xlabel(r'R $(mm/hr)$')
                    ax.set_ylabel(r" $ \langle \langle \rho_{%g}(x, t) \rangle_{x} \rangle_{t} $" % (s))
                    # Set lower limit of y-axis to 0.
                    ax.set_ylim(bottom = -50)
                    if(s == 2):
                        ax.set_ylim(top = 100)
                    #Set the legend in the lower right corner.
                    ax.legend(loc = 'lower right')
                # End of s loop 
            # End of val loop
            Twin_min = Tavg_TSval_all.index.get_level_values("Twin_min").min(); Twin_max = Tavg_TSval_all.index.get_level_values("Twin_max").max()
            fig.suptitle(r"$\mu$ and $\sigma$ of" + f" Species For Various {compare_key} B/W {Twin_min} -- {Twin_max}")
            plt.savefig(savedir + f"HsX2005_COMPARE_{compare_key}_MeanSTD_amin_{Avals_list[0]}_amax_{Avals_list[-1]}_Twin_{Twin_min}_{Twin_max}.png")

            plt.show()
            plt.close()
        # End of TS loop

        # Next, create box/violin plots of the data for each species, iterating over the values in the compare_vals list as hues.
        if(plt_violin):
            pass
            # Extract data for box/violin plots.

    # End of compare_key loop




            
''' # Summary Description of recursive_copydir(...)
This function recursively copies files from the source directory to the destination directory.
src: source directory
dst: destination directory
include_filetypes: list of filetypes to be copied
symlinks: If True, symlinks are copied as symlinks. If False, the files are copied as hardlinks.
'''    
def recursive_copydir(src, dst, include_filetypes = ["*.txt"], 
                      exclude_filetypes =["*.csv", "*.png", "*.jpg", "*.jpeg", "*.mp4"], symlinks=False):
    for root, dirs, files in os.walk(src):
        for file in files:
            if any(fnmatch.fnmatch(file, pattern) for pattern in include_filetypes):
                if not any(fnmatch.fnmatch(file, pattern) for pattern in exclude_filetypes):
                    src_file = os.path.join(root, file)
                    dst_file = os.path.join(dst, os.path.relpath(src_file, src))
                    os.makedirs(os.path.dirname(dst_file), exist_ok=True)
                    if symlinks:
                        if os.path.exists(dst_file):
                            os.remove(dst_file)
                        os.symlink(src_file, dst_file)
                    else:
                        shutil.copy2(src_file, dst_file)

    print("\n\n__________________________________________________________________________________________\n\n")
    print(f"Done copying files from {src} to {dst}.")
    print("\n\n__________________________________________________________________________________________\n\n")


# ============================= PRELIM DATA ANALYSIS FUNCTIONS =============================



#
recursive_copydir(in_dir, out_dir, include_filetypes = ["*.txt"], exclude_filetypes =["*.png", "*.jpg", "*.jpeg", "*.mp4"], symlinks=False)
a_vals = [1.7, 1.72, 1.7225, 1.725, 1.7275, 1.73, 1.735, 1.74, 1.76, 1.78, 1.8]
#[0.02, 0.022, 0.023, 0.0235, 0.024, 0.0245, 0.025, 0.0255, 0.026, 0.028, 0.03, 0.034, 0.036, 0.038, 0.0415, 0.042, 0.044, 0.046, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.06, 0.065, 0.07]
# #1.7, 1.72, 1.7225, 1.725, 1.7275, 1.73, 1.7325, 1.735, 1.74, 1.76, 1.78, 1.8] #1.75, 1.8, 20] #0.026, 0.034, 0.0415, 0.042, 0.05, 0.052] #, 0.057 , 0.06] #0.051, 0.053, 0.055]; 
#0.022, 0.023, 0.0235, 0.024, 0.245, 0.25, 0.255, 0.028, 0.038, 0.044, 0.046, 0.051, 0.053, 0.054, 0.055 ]  
#1.7, 1.78, 1.8, 1.84, 1.9, 2, 5, 20 ; 
#1.66, 1.67, 1.68, 1.686, 1.688, 1.69, 1.692, 1.694, 1.696, 1.698, 1.7, 1.72, 1.74, 1.76, 1.78, 1.8 ;
a_scaling = 0.001 #0.001
T_vals=[120226]
#T_vals= [0, 76000, 80000, 82000, 84000, 86000, 88000, 90000, 91201.1]#0, 80000, 84000,  88000, 92000, 96000, 98000, 100000, 104000, 108000, 112000, 116000, 120000, 120226] # [0, 80000.1, 82000.1, 84000.1, 86000.1, 88000.1, 90000.1, 92000.1, 94000.1, 96000.1, 100000.1]
#80000, 82000.1, 84000, 86000, 88000, 90000, 92000, 94000, 96000, 100000, 120226, 144544]#[0, 75857.7, 83176.3, 91201, 100000, 109648, 120226, 131826]
#T_vals= [0, 63.03, 109.56, 144.54, 190.52, 251.13, 331.1, 436.48, 575.41, 758.56, 831.71, 999.9, 1202.19, 1445.4, 1737.78, 2089.23, 2511.85, 3019.94, 3630.77, 4365.13, 5247.99, 6309.49, 6918.23, 7585.71, 8317.54, 9120.1, 9999.99]
#T_vals=[0, 63095.7, 69182.9, 75857.7, 83176.3, 91201, 99999.9, 109647, 120226, 131826, 144544, 158489, 173780, 190546, 208930, 229087];
# TVALS WHEN DT= 0.1
#T_vals=[0, 63095.7, 69183, 75857.7, 83176.3, 91201, 100000, 109647, 120226, 131826, 144544, 158489, 173780, 190546];
# TVALS WHEN DT= 0.12
#T_vals=[ 0, 20.88, 30.12, 39.72, 50.04, 99.96, 150, 200.04, 249.96, 300, 350.04, 399.96, 500.04, 600, 699.96, 800.04, 900, 999.96, 1100.04, 1200, 1318.2, 1445.4, 1737.72, 2089.2, 2511.84, 3019.92, 3630.72, 4365.12, 5274.96, 6309.48, 6918.24, 7585.68, 9120.1, 9999.96, 10964.8, 12022.6, 13182.5, 14454.4 ];#0, 63095.6, 69183, 75857.6, 83176.3, 91201, 100000, 109647, 120226, 131826, 144544, 158489, 173780, 190546];
#[0, 20.88, 30.12, 39.72, 50.04 99.96, 150, 200.04, 249.96, 300, 350.04, 399.96, 500.04, 600, 699.96, 800.04, 900, 999.96, 1100.04, 1200, 1318.2, 1445.4, 1737.72, 2089.2, 2511.84, 3019.92, 3630.72, 4365.12, 5274.96, 6309.48, 6918.24, 7585.68, 9120.1, 9999.96, 10964.8, 12022.6, 13182.5, 14454.4] 
#[0, 33, 47.85, 68.75, 99.55, 144.1, 208.45, 250.8, 301.95, 363, 436.15, 524.7, 600.05, 649.99, 700.04, 749.98, 800.03, 849.97, 900.02, 949.96, 1000.01, 1049.95, 1100, 1150.05, 1199.99, 1250.04, 1299.98, 1350.03, 1399.97, 1450.02, 1499.96, 1584.55, 1737.45, 1905.2, 7585.6] 
#[0, 144.54, 190.52, 251.13, 331.1, 436.48, 575.41, 758.56, 831.71, 999.9, 1202.19, 1445.4, 1737.78, 2089.23, 2511.85, 3019.94, 3630.77, 4365.13, 5247.99, 6309.49, 6918.23, 7585.71, 9120.1, 9999.99, 10964.7, 91201, 63095.7, 69183.1, 75857.7, 91201,  99999.9, 10964.3, 120226, 131826, 144544, 158489]
# 3SP INIT [0, 575.3, 758.45, 911.9, 1096.15, 1317.8, 1584.55, 1905.2, 2290.75, 2753.85, 3311, 3980.7, 4786.1, 5754.1, 6309.05, 6917.9, 7585.6, 8317.1, 9120.1, 9999.55, 91201, 99999.9, 10964.3, 120226, 131826, 144544, 158489]
#T_vals= [ 109648, 120226, 131826, 144544, 158489, 173780, 190546, 208930]
#T_vals=[0, 63095.7, 69182.9, 75857.7, 83176.3, 91201, 99999.9, 109647, 120226, 131826, 144544, 158489, 173780, 190546, 208930]
#T_vals= [158489, 173780, 190546, 208930, 229087, 251189, 275423, 301995, 331131, 363078]
# TVALS FOR GAMMA CHECK:
#T_vals= [0, 82000.1, 84000, 86000, 88000, 90000, 92000, 94000, 96000, 98000, 100000, 144544, 158489, 173780, 190546]
#T_vals= [0, 76000, 80000, 82000.1, 84000, 86000, 88000, 90000, 92000, 94000, 96000, 98000, 100000, 144544, 190546]
#T_vals =[4000.04, 3980.02, 3960, 3900.05, 3799.95, 3699.96, 3499.98, 3399.99, 3200.01, 3000.03, 
#        2800.05, 2599.96, 2399.98, 2200, 2000.02, 1800.04, 1599.95, 1399.97, 1199.99, 1000.01,
#        800.03, 600.05, 499.95, 399.96, 340.01, 330, 319.99, 299.97, 230.01, 220, 209.99, 199.98, 120.01, 110, 99.99, 0]
#        #800.03, 600.05, 499.95, 399.96, 299.97, 220, 199.98, 110, 99.99, 0]

#../Data/DP/Reorganised_Frames/Stoc/3Sp/DPParam_20_100_DsB6HX/DsB6-HX2001-UA0A0-1UNI
#../Data/DP/Reorganised_Frames/Stoc/3Sp/DPParam_20_100_DSB6HX/DsB6-HX2001-UA0A0-1UNI
#prefixes = ["K-S0B6X1-UNITY", "K-S-1B6X1-UNITY" ]#, "HsX2001-UA0A0-5E2UNI"] 
#"DsB6L9-UA0A0-1UNI", "DiC-L9B6-UNITY", "DsB6-UA0A0-1UNI"
#"HsX2001-UA125A0-5E2UNI", "HsX2005-UA125A0-5E2UNI", "HsX2001-UA0A0-5E2UNI"]
#prefixes = ["DsC-HXL2010-UA125A0-5E2UNI", "DsCGm-HXL2010-UA125A0-5E2UNI", "DsCGm-HXR03015-UA125A0-5E2UNI",  "DsC-HXR03015-UA125A0-5E2UNI"]
#"DsC-HXR03015-UA125A0-5E2UNI", "DsC-HXR2010-UA125A0-5E2UNI", "DsC-HXL03015-UA125A0-5E2UNI", "DsC-HXL2010-UA125A0-5E2UNI", "DsC-HXC03015-UA125A0-5E2UNI", "DsC-HXC2010-UA125A0-5E2UNI", "DsCGm-HXR03015-UA125A0-5E2UNI", "DsCGm-HXL03015-UA125A0-5E2UNI", "DsCGm-HXL2010-UA125A0-5E2UNI", "HX03002-UA125A0-1UNI",  "HX03002-UA125A0-5E2UNI"]

#"HsX2001-UA125A0-5E2UNI", "HsX2005-UA125A0-5E2UNI", "HsX2001-UA125A125-5E2UNI", "HsX2005-UA125A125-5E2UNI"] 
#"HX03002-UA125A125-1UNI",  "HX03002-UA125A125-5E2UNI", "DsCX-TST-HX03002-UA125A125-5E2UNI" ]
#"DisC-UA125A125-1UNI", "DsC-UA125A125-1UNI", "DsCX-TST-HX2001-UA125A125-5E2UNI" ]
#"DsC-UA0A0-S1DTUNI", "DsCX-TST-UA0A0-S1DTUNI",  "DsC-UA125A0-S1DTUNI", "DsC-UA125A0-S1DT5E2UNI", "DsCX-TST-UA125A0-S1DTUNI", "DsC-UA125A125-S1DTUNI", "DsC-UA125A125-S1DT5E2UNI", "DsCX-TST-UA125A125-S1DTUNI" ]
#"B6-UA0A0-1UNI", "B6-UA125A0-1UNI", "B6-UA125A125-1UNI"]
#"DPB6-HX2001-UA0A0-1UNI", "DPB6-UA0A0-REGPRG-1UNI",  "DPB6-HX2001-UA125A0-1UNI", "DPB6-UA125A0-REGPRG-1UNI", "DPB6-UA125A125-REGPRG-1UNI"]
#"B6-HX2001-UA0", "B6-HX2005-UA0", "B6-HX2001-UA125", "B6-HX2005-UA125"]
#"RTK-UA125-LOPRG-1UNI", "RTK-UA125-REGPRG-1UNI" ]
#"RTK-UA125A0-REGPRG-1UNI", "RTK-UA125A0-REGPRG-5E2UNI",  "RTK-UA125A0-LOPRG-1UNI", "RTK-UA125A0-LOPRG-5E2UNI"]
#"HX1001-UA125A0-5E2UNI", "HX1005-UA125A0-5E2UNI", "HX2001-UA125A0-5E2UNI", "HX2005-UA125A0-5E2UNI", "HX3001-UA125A0-5E2UNI"]           
#"DiC-UA125A0-1UNI", "DiC-UA125A0-5E2UNI"]
#"DiC-UA125A0-S1DTUNI", "DiC-UA125A0-S1DT5E2UNI"]   #"DiC-UA125A0-1UNI", "DiC-UA125A0-5E2UNI"]#"DiC-UA125A2-1UNI", "DiC-UA125A2-5E2UNI"]   #"DiC-UA125A0-1UNI", "DiC-UA125A0-5E2UNI"]    #"DiC-G0A1A0-1LI"]  
#"DiC-B6-UNITY"] #"DiC-B6-UNTY" 
#"DiC-B6-MFTEQ"]#"DiC-STD"]#,"DiC-S7LI", "DiC-0.1LI"]
prefixes = ["DsB6L9-UA0A0-1UNI"]
#prefixes = ["cuFKD-GAU-SLV1V1", "cuFKD-GAU-SLV2V1", "cuFKD-GAU-SLV5V4", "cuFKV-GAU-SLV1V1", "cuFKV-GAU-SLV2V1", "cuFKV-GAU-SLV5V4"]
#prefixes = ["HsX2001-UA0A0-5E2UNI", "HsX2001-UA125A0-5E2UNI", "HsX2005-UA0A0-5E2UNI", "HsX2005-UA125A0-5E2UNI"]
#"HsX2001-UA0A0-5E2UNI", "HsX2001-UA125A0-5E2UNI", "HsX2005-UA0A0-5E2UNI", "HsX2005-UA125A0-5E2UNI"]#, "HsX2001-UA125A125-5E2UNI"] #] "HsX2001-UA0A0-5E2UNI", "HsX2001-UA125A0-5E2UNI", "HsX2005-UA125A0-5E2UNI", "HsX2005-UA0A0-5E2UNI"
#prefixes = [ "CORR-DsB6-UA0A0-1UNI"] #"CORR-DsB6-UA0A0-1UNI"] #, "CORR-HsX2001-UA125A0"]#, "CORR-HsX2001-UA125A125", "CORR-HsX2005-UA125A125"]
# PREFIXES FOR SMALL BODY SIZE METAPOPLN

#LOCUSTS
#prefixes = ["DsCW10-HXR2010-UA0A0", "DsCW10-HXR2010-UA125A0", "DsC-HXR2010-UA0A0",  "DsC-HXR2010-UA125A0",
#            "DsCW10-HXR03015-UA0A0", "DsCW10-HXR03015-UA125A0", "DsC-HXR03015-UA0A0",  "DsC-HXR03015-UA125A0" ]
#prefixes = ["DsCW10-HXR2010-UA125", "DsCW10-HXR03015-UA125", "DsC-HXR2010-UA125"]
#prefixes = ["DsCW10-HXR2010-UA0A0-DT1", "DsCW10-HXR2010-UA125A0-DT1", "DsCW5-HXR2010-UA0A0-DT1", "DsCW5-HXR2010-UA125A0-DT1",
#            "DsC-HXR2010-UA0A0-DT1", "DsC-HXR2010-UA125A0-DT1" ]
#["DsCW10-HXR2010-UA125A125", "DsC-HXR2010-UA125A125"]#  "DsCW10-HXR03015-UA125A125", "DsC-HXR03015-UA125A125"]


# GENERIC SMALL MAMMAL
#prefixes = ["DsC-HXR03015-UA0-1UNI", "DsC-HXR2010-UA0-1UNI"] #"DsC-HXR03015-UA0-1UNI",
#prefixes = ["DsC-HXR03015-UA0A0-1UNI", "DsC-HXR2010-UA0A0-1UNI", "DsC-HXR03015-UA125A0-1UNI", "DsC-HXR2010-UA125A0-1UNI"] #"DsC-HXR03015-UA0-1UNI",

TS_vals = [120226]  #33113.1] 36307.7]#, 131826] 145750] #190546] #57544] #69183.1] # #[229087] #[208930] #91201]  #[190546]; #[109648];
#a_vals = []#0.034, 0.048, 0.054]; 
#T_vals = []#0, 91201, 190546, 208930, 229087]
TCorr_vals =[4000.04]

# col1 = "GAM[G(x; t)]"; col2 = "G(x; t)"; NCC_index_key = f"NCC-Index{{{col1};{col2}}}_{g}"

print("Note the following set values:")
print(f"Prefixes: {prefixes}")
print(f"TS_vals: {TS_vals}")
print(f"T_vals: {T_vals}")
print(f"a_vals: {a_vals}")
print(f"SPB: {SPB}, g: {g}, dP: {dP}, Geq: {Geq}, maxR: {R_max}")
print(f"in_dir: {in_dir} \nout_dir: {out_dir}")
#print(f"NCC_index_key: {NCC_index_key}")
input("Press F to pay respects...")
print("\n\n__________________________________________________________________________________________\n\n")

''' FOR VIDEOS OF ALL REPLICATES
for Pre in prefixes:
    frame_visualiser(in_dir, out_dir, Pre, a_vals, T_vals, maxR= -1, plt_gamma= False, delpng = False)
    print(f"Done with making frames for {Pre}")
    for a in a_vals:
        print(f"Making video for {Pre} at a = {a} \n\n")
        home_video(out_dir, out_dir, [Pre], [a], T_vals, maxR= -1, 
                   pngformat= "CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png", pngdir= "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/",
                     maxR_txtdir = "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt", videoname = "guess",
                     video_relpath= "{Pre}/Videos/Conc/{a}/{Tmin}-{Tmax}/")
#'''

''' FOR VIDEOS OF INDIVIDUAL REPLICATES
for Pre in prefixes:
    for R in range(0, 3):
        frame_visualiser(in_dir, out_dir, Pre, a_vals, T_vals, maxR= R+1, minR= R, plt_gamma= True, delpng = False)
        print(f"Done with making frames for {Pre}")
        for a in a_vals:
            print(f"Making video for {Pre} at a = {a} \n\n")
            home_video(out_dir, out_dir, [Pre], [a], T_vals= T_vals, maxR= R+1, minR= R,
                    pngformat= "CombinedImg/BioGammaConc_a_{a}_T_{T}_n_{R}.png", pngdir= "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/",
                        maxR_txtdir = "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt", videoname = "guess",
                        video_relpath= "{Pre}/Videos/Conc/{a}/{Tmin}-{Tmax}/")
#'''

''' FOR VIDEOS OF ALL REPLICATES ACROSS A AND FOR GIVEN T
for Pre in prefixes:
    frame_visualiser(in_dir, out_dir, Pre, a_vals, T_vals, maxR= -1, plt_gamma= False, delpng = False)
    print(f"Done with making frames for {Pre}")
    #for a in a_vals:
    #print(f"Making video for {Pre} at a = {a} \n\n")
    home_video(out_dir, out_dir, [Pre], a_vals, T_vals, maxR= -1, 
            pngformat= "CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png", pngdir= "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/",
            maxR_txtdir = "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt", videoname = "guess",
            video_relpath= "{Pre}/Videos/Conc/{amin}-{amax}/{Tmin}-{Tmax}/")
#'''

''' FOR SPREADING TESTS, WITH PREFIX = NREF-GAU/ REF-GAU
for Pre in prefixes:
    for R in range(0, 1):
        frame_visualiser(in_dir, out_dir, Pre, a_vals, T_vals, maxR= R+1, minR= R, plt_gamma= True, delpng = False)
        #frame_visualiser(in_dir, out_dir, Pre, a_vals, T_vals, maxR= -1, plt_gamma= True, delpng = False)
    print(f"Done with making frames for {Pre}")
    for a in a_vals:
        print(f"Making video for {Pre} at a = {a} \n\n")
        home_video(out_dir, out_dir, [Pre], [a], T_vals= T_vals, maxR= R+1, minR= R, 
                   pngformat= "CombinedImg/BioGammaConc_a_{a}_T_{T}_n_{R}.png", pngdir= "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/",
                     maxR_txtdir = "{Pre}/L_{g}_a_{a}/dP_{dP}/Geq_{Geq}/T_{T}/maxR.txt", videoname = "guess",
                     video_relpath= "{Pre}/Videos/Conc/{a}/{Tmin}-{Tmax}/")
#'''
#frame_visualiser(in_dir, out_dir, prefixes[0], a_vals, T_vals, maxR= -1, plt_gamma= False, delpng = False)
#home_video(out_dir, out_dir, prefixes, a_vals, T_vals, R_max, pngformat= "CombinedImg/BioConc_a_{a}_T_{T}_n_{R}.png", videoname = "guess")

#prefixes = ["DIC-NREF-1.1HI", "DIC-NREF-0.5LI", "DIC-NREF-0.1LI"]
#analyse_FRAME_EQdata(in_dir, out_dir, prefixes, a_vals, T_vals, Tavg_window_index = [150000, 240000], filename = "MEAN_STD_Surviving_Runs.txt")
#analyse_timeseriesData(in_dir, out_dir, prefixes, a_vals, T_vals, filename = "MEAN_REPLICATES.txt")
if SPB == 3:
    variable_labels = [ "<<P(x; t)>_x>_r" , "<<G(x; t)>_x>_r", "<<Pr(x; t)>_x>_r"]
    var_frame_labels = ["P(x; t)", "G(x; t)", "Pr(x; t)"]
    move_var_labels = [ "<GAM[G(x; t)]>_x" , "<GAM[Pr(x; t)]>_x"]
    Tavg_win_index = [190000, 200000]; Traj_win_index =[100, 3000]; Tavg_win_index_move= [190000, 200000];#[80000, 150000]; #[100000, 240000]
    # Window for averaging over last Tavg_win_index values of T (if negative, then average over last |Tavg_win_index| values of T)
    # If Tavg_win_index[0] < 0 and Tavg_win_index[1] <= 0, then average over last Tavg_win_index[0] + Tmax to Tavg_win_index[1] + Tmax values of T.
    # If Tavg_win_index[1] >= Tavg_win_index[0] >= 0, then average over last Tavg_win_index[0] to Tavg_win_index[1] values of T.
elif SPB == 2:
    variable_labels = [ "<<P(x; t)>_x>_r" , "<<G(x; t)>_x>_r"]; 
    var_frame_labels = ["P(x; t)" , "G(x; t)"]
    move_var_labels = [ "<GAM[G(x; t)]>_x"]
    Tavg_win_index = [160000, 200000]; Traj_win_index =[1000, 10000]; Tavg_win_index_move= [80000, 150000];#[150000, 240000] #[50000, 100000] #
elif SPB == 1:
    variable_labels = ["<<P(x; t)>_x>_r"]; var_frame_labels = ["P(x; t)"];  Tavg_win_index = [120000, 125000]; Traj_win_index =[200, 5000]; #[65000, 100000]
    move_var_labels = [ "<GAM[G(x; t)]>_x"]; Tavg_win_index_move= [120000, 125000];

analyse_PRELIMS_TIMESERIESdata(in_dir, out_dir, prefixes, a_vals, TS_vals, serType = "TSERIES", meanfilename = "MEAN_TSERIES_T_{TS}.csv", a_c=[1.72, 1.725, 1.7275, 1.73], var_labels= variable_labels, a_scaling= a_scaling)
#analyse_PRELIMS_EQdata(in_dir, out_dir, prefixes, a_vals, TS_vals, Tavg_window_index = Tavg_win_index, meanfilename = "MEAN_TSERIES_T_{TS}.csv", var_labels= variable_labels, a_scaling= a_scaling)
#analyse_PRELIMS_EQdata(in_dir, out_dir, prefixes, a_vals, TS_vals, Tavg_window_index = Tavg_win_index, serType = "TSERIES", meanfilename = "MEAN_TSERIES_T_{TS}.csv", var_labels= variable_labels, plt_violin= True, a_scaling= a_scaling)
'''
for a in a_vals:
    analyse_PRELIMS_TRAJECTORYdata(in_dir, out_dir, prefixes, [a], TS_vals, size_label = ["t"], maxR= 5, minR=0, T_window = Traj_win_index, meanfilename = "Mean_TSERIES_T_{TS}_dP_{dP}_Geq_{Geq}.csv", a_scaling= a_scaling, 
        x_label= ["<<G(x; t)>_x>_r"], y_label= ["<<Pr(x; t)>_x>_r"], hue_label=["<<P(x; t)>_x>_r"])
#'''

# For FRAME DATA
#analyse_FRAME_FFTPOWERdata(in_dir, out_dir, prefixes, a_vals, T_vals, filename = "FFT_POWERspectra.csv", Tavg_window_index = Tavg_win_index, show_fundamental_freq=True, var_labels= var_frame_labels, a_scaling= a_scaling)
#analyse_FRAME_POTdata(in_dir, out_dir, prefixes, a_vals, T_vals, find_minima= True, filename = "Pot_Well.csv", 
#                          minimafilename = "LOCAL_MINIMA.csv", var_labels= [ "P(x; t)" , "G(x; t)", "Pr(x; t)"], a_scaling= a_scaling)
#analyse_FRAME_CLUSTEREDISTdata(in_dir, out_dir, prefixes, a_vals, T_vals, filename = "CLUSTERED_FREQUENCIES.csv", clustsubdir ="CLUST/{binType}_{nbins}/",
#        binType="GMM", nbins=2,  var_labels= var_frame_labels, plot_binframes=False, CDF=True, a_c=[1.725], a_scaling= a_scaling)

# For MOVEMENT DATA
#analyse_PRELIMS_TIMESERIESdata(in_dir, out_dir, prefixes, a_vals, TS_vals, serType = "MOVSERIES", meanfilename = "MEAN_{serType}_T_{TS}.csv", var_labels= move_var_labels, a_scaling= a_scaling)
#analyse_PRELIMS_EQdata(in_dir, out_dir, prefixes, a_vals, TS_vals, Tavg_window_index = Tavg_win_index_move, serType = "MOVSERIES",  
#                       meanfilename = "MEAN_{serType}_T_{TS}.csv", var_labels= move_var_labels, plt_violin= True, a_scaling= a_scaling)

#analyse_FRAME_POTdata(in_dir, out_dir, prefixes, a_vals, T_vals, find_minima= True, filename = "Gamma_Pot_Well.csv", 
#                          minimafilename = "Gamma_LOCAL_MINIMA.csv", var_labels= [ "GAM[G(x; t)]" , "GAM[Pr(x; t)]"])
''' # For MEAN and Timeseries of CORR DATA.
analType= ["NCC", "ZNCC", "AMI", "MI", "BVMoransI"]
for anal in analType:
    analyse_FRAME_CORRdata(in_dir, out_dir, prefixes, a_vals, TCorr_vals, analType = anal, 
                           filename = "{corrType}_{analType}_TD_*.csv", var_labels= "deduce", a_scaling= a_scaling)
#'''
''' # For HARMONIC PEAK (FFT ANALYSIS) of CORR DATA.
analType= ["NCC", "ZNCC", "AMI", "MI", "BVMoransI"]
for anal in analType:
    print(f"Analysing {anal} FFT data..."); time.sleep(1)
    analyse_FRAME_FFTCORRdata(in_dir, out_dir, prefixes, a_vals, TCorr_vals, analType = anal, corrsubdir="2DCorr/FFT/",
            harm_peaks_index = [0], filename = "HarmonicPeaks_{corrType}_{analType}_TD_*.csv", 
            var_labels= "deduce", a_scaling= a_scaling)
#'''

#compare_list = {"Prefixes": [], "g": [], "dP": [100, 10000], "Geq": []}
#compare_list = {"Prefixes": ["HsX2005-UA125A0-5E2UNI", "HsX2005-UA125A125-5E2UNI"], "g": [], "dP": [], "Geq": []}
#multiplot_PRELIMS_CompareEQ(out_dir, out_dir, compare_list, prefixes = prefixes, meanfilename = "guess", var_labels= variable_labels)