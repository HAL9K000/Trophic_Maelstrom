import numpy as np
import os
import pandas as pan
import regex as re
from pathlib import Path
import glob
import shutil
import sys
import argparse
import copy
import warnings
import time

import scipy.stats as stats

from legacy_glow_up import *

# Show the first FutureWarning that occurs in the script, then ignore all future FutureWarnings.
warnings.simplefilter(action='once', category=FutureWarning)
warnings.simplefilter(action='once', category=Warning)
warnings.simplefilter(action='once', category=pan.errors.PerformanceWarning)
'''
This script reorganises the directory structure of the data files in the Rietkerk model.
The original directory structure (filepaths) are as follows:
root_dir/{PREFIX}*_dP_{dP}_Geq_{Geq}/TimeSeries/PRELIM_TSERIES_*_G_{g}_T_{T}_*_a_{a_val}_*_R_{R}.csv

The new directory structure (filepaths) will be as follows:
out_dir/{PREFIX}/L_{g}_a_{a_val}/dP_{dP}/Geq_{Geq}/TimeSeries/PRELIM_T_{T}_a_{a_val}_R_{R}.csv

The script will:
1. Find all immediate subdirectories in root_dir.
2. For each subdirectory, find all unique values of a_val in the files in the subdirectory.
3. For each value of a_val, find all unique values of T in the files in the subdirectory.
4. For each value of T, copy all files in the subdirectory to the new directory structure in out_dir.
5. If a file with the same name already exists in the new directory, rename the file to avoid conflicts.

This script ASSUMES that the files in the original directory structure have the following format,
and WILL NOT WORK if the files are named differently:
I. Contains "_G_{L}_T_{Time}_" in the filename, which is followed by
II. "_a_{a_val}_" in the filename, which is followed by
III. "_R_{R}.csv" in the filename.

The script accepts the following arguments:
1. prefixes: A list of prefixes to be used in the subdirectory names in out_dir.
2. root_dir: The root directory containing the original data files.
3. out_dir_noprefix: The output directory where the reorganised data files will be stored.
4. dP: The value of dP in the subdirectory names in root_dir.
5. Geq: The value of Geq in the subdirectory names in root_dir (if "NA", Geq is not used in the subdirectory name).
6. L: A list of values of L to be used in the subdirectory names in out_dir.
7. indx_vals_t: Extract n largest values of T in each sub-directory if indx_vals_t = -n, n smallest values of T if indx_vals_t = n.

'''

prefixes =["DiC-NREF-1.1HI", "DiC-NREF-0.5LI", "DiC-NREF-0.1LI"]#, "DiC-NEW"]
#prefixes =["DiC-NREF-HI", "DiC-NREF-LI"]
#prefixes =["", "DiC", "BURNIN", "DiC-BURNIN", "DDM-DiC", "DDM-DiC-BURNIN"]
root_dir = "../Data/Amarel/Rietkerk/Prelims/Stochastic/3Sp/"
#out_dir_noprefix = "../Data_Processing/StdParam_20_100_Test/"
#root_dir = out_dir_noprefix

out_dir_noprefix = "../Data/Remote/Rietkerk/Reorg_Frames/3Sp/StdParam_20_100_MFTEQ/"

dP = 10000
Geq = 4.802 # Optional. If Geq is not used in the subdirectory name, set Geq = "NA".
Veq = 7.4774  # Optional. If Veq is not used in the subdirectory name, set Veq = "NA".
L= [128]
tmin = 100000; tmax = None; 
# Optional. If tmin and tmax are not provided, set tmin = None and tmax = None.
indx_vals_t = -5; dt =0.11
#Extract n largest values of T (within the subset [tmin, tmax] if provided) if indx_vals_t = -n, 
# n smallest values of T if indx_vals_t = n.
CPU_Ncores = 1; # Number of CPU cores to be used by the script (monothreaded by default).
dynamic_inspect = False;    # Set to True to stop script periodically and inspect values of parameters and outputs.


def set_prelims_inputs():
    parser = argparse.ArgumentParser(description="Set Prelim Inputs for the script.")

    # Boolean flag for default input values.
    parser.add_argument("--default", action="store_true", help="Use default input values for the script.")
    # Boolean flag for dynamic input values.
    parser.add_argument("--dynamic", action="store_true", help="Use dynamic input values for the script.")
    #Number of CPU cores to be used by the script (monothreaded by default).
    parser.add_argument("--CPUCores", type=int, help="Number of CPU cores to be used by the script (monothreaded by default).")
    # List of prefixes to be used in the subdirectory names in out_dir.
    parser.add_argument("--prefixes", nargs="+", help="List of prefixes to be used in the subdirectory names in out_dir.")
    # Root directory containing the original data files.
    parser.add_argument("--indir", help="Root directory containing the original data files.")
    # Output directory where the reorganised data files will be stored.
    parser.add_argument("--outdir", help="Output directory where the reorganised data files will be stored.")
    # Value of dP in the subdirectory names in root_dir.
    parser.add_argument("--dP", type=int, help="Value of dP in the subdirectory names in root_dir.")
    # Value of Geq in the subdirectory names in root_dir (if "NA", Geq is not used in the subdirectory name).
    parser.add_argument("--Geq", help="Value of Geq in the subdirectory names in root_dir, use 'NA' if not used.")
    # List of values of L to be used in the subdirectory names in out_dir.
    parser.add_argument("--Veq", help="Value of Veq in the subdirectory names in root_dir, use 'NA' if not used.")
    # List of values of L to be used in the subdirectory names in out_dir.
    parser.add_argument("--L", nargs="+", help="List of values of L to be used in the subdirectory names in out_dir.")
    # Extract n largest values of T in each sub-directory if indx_vals_t = -n, n smallest values of T if indx_vals_t = n.
    parser.add_argument("--indx_vals_t", type=int, help="Extract n largest values of T in each sub-directory if indx_vals_t = -n, n smallest values of T if indx_vals_t = n.")
    parser.add_argument("--dt", type=float, help="Value of dt in the subdirectory names in root_dir.")

    # Parse the arguments.
    args = parser.parse_args()

    global prefixes, root_dir, out_dir_noprefix, dP, Geq, Veq, L, indx_vals_t, dt, dynamic_inspect, CPU_Ncores

    # Set the input values.
    if args.default:
        print("Using default input values for the script.")
        dynamic_inspect = True
        return
    if args.dynamic:
        print("Using dynamic input values for the script.")
        print("\n=====================================================================================")
        print("Current values are:")
        print("prefixes: " + str(prefixes))
        print("root_dir: " + root_dir)
        print("out_dir_noprefix: " + out_dir_noprefix)
        print("dP: " + str(dP))
        print("Geq: " + str(Geq))
        print("L: " + str(L))
        print("indx_vals_t: " + str(indx_vals_t))
        print(f"dt: {dt}")
        print(f"Also note, using {CPU_Ncores}/{os.cpu_count()} CPU cores for the script.")
        print("=====================================================================================\n")
        prefixes = input("Enter a list of prefixes to be used in the subdirectory names in out_dir: ").split()
        root_dir = input("Enter the root directory containing the original data files: ")
        out_dir_noprefix = input("Enter the output directory where the reorganised data files will be stored: ")
        dP = int(input("Enter the value of dP in the subdirectory names in root_dir: "))
        try:
            Geq = float(input("Enter the value of Geq in the subdirectory names in root_dir, use 'NA' if not used: "))
        except ValueError:
            print("Setting to NA..."); Geq = "NA"
        L = input("Enter a list of values of L to be used in the subdirectory names in out_dir: ").split()
        L = [int(l) for l in L]
        indx_vals_t = int(input("Enter the value of indx_vals_t: "))
        dt = float(input("Enter the value of dt in the subdirectory names in root_dir: "))
        CPU_Ncores = int(input("Enter the number of CPU cores to be used by the script (monothreaded by default): "))
        dynamic_inspect = True
        return
    
    else:
        if args.prefixes:
            prefixes = args.prefixes
            #Remove newline, tab and trailing or leading whitespaces from prefixes.
            prefixes = [pre.strip("\n\t ").lstrip().rstrip() for pre in prefixes]
        if args.indir:
            root_dir = args.indir
        if args.outdir:
            out_dir_noprefix = args.outdir
        if args.dP:
            try:
                dP = int(args.dP)
            except ValueError:
                print("Need to provide an integer value for dP. Terminating script..."); sys.exit(1)
        if args.Geq:
            try:
                Geq = float(args.Geq)
                Geq = int(Geq) if Geq.is_integer() else Geq
            except ValueError:
                print("Setting Geq to NA..."); Geq = "NA"
            #Geq = float(args.Geq) if isinstance(args.Geq, (float, int)) else "NA"
        if args.Veq:
            try:
                Veq = float(args.Veq)
                Veq = int(Veq) if Veq.is_integer() else Veq
            except ValueError:
                print("Setting Veq to NA..."); Veq = "NA"
            #Veq = float(args.Veq) if isinstance(args.Veq, (float, int)) else "NA"
        if args.L:
            L = args.L; L = [int(l) for l in L]
        if args.indx_vals_t:
            try:
                indx_vals_t = int(args.indx_vals_t)
            except ValueError:
                print("Setting indx vals to -10..."); indx_vals_t = -10;
        if args.dt:
            try:
                dt = float(args.dt)
            except ValueError:
                print("No valid value for dt provided. Terminating script..."); sys.exit(1)
        if args.CPUCores:
            try:
                CPU_Ncores = int(args.CPUCores)
                CPU_Ncores = CPU_Ncores if CPU_Ncores <= os.cpu_count() else os.cpu_count()
            except ValueError:
                print("Invalid core selection. Setting CPU cores to 1..."); CPU_Ncores = 1

        print("prefixes: " + str(prefixes))
        print("root_dir: " + root_dir)
        print("out_dir_noprefix: " + out_dir_noprefix)
        print("dP: " + str(dP))
        print("Geq: " + str(Geq))
        print("Veq: " + str(Veq))
        print("L: " + str(L))
        print("indx_vals_t: " + str(indx_vals_t))
        print(f"dt: {dt}")
        print(f"Also note, using {CPU_Ncores}/{os.cpu_count()} CPU cores for the script.")

        return



def post_process(prefixes= []):
    #Assigning prefixes to all subdirectories in out_dir_noprefix if prefixes is empty.
    if len(prefixes) == 0:
        prefixes = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(out_dir_noprefix, '*/'))]
    
    for pre in prefixes:

        workdir = os.path.join(out_dir_noprefix, pre)
        # Recursovely navigate to each sub-directory within out_dir_noprefix/{pre} that contains files of the form "MOVSERIES*.csv".
        # These subdirectories may be immediate or nested. For each such non-empty nested sub-directory, use glob to find all files of the form "MOVSERIES*.csv".
        # Pass this list of files to other functions for further processing.
        all_valid_subdirpaths = [dir for dir, _, files in os.walk(workdir) if len(glob.glob(dir + "/MOVSERIES*.csv")) > 0]
        for subdirpath in all_valid_subdirpaths:
            files = glob.glob(subdirpath + "/MOVSERIES*.csv")

            # FIRST, getting maximum replicate number (maxR) for each value of T in subdirpath.
            # Read list of T values from TimeSeries_vals.txt in subdirpath (if it exists), else create it using the files.
            T_vals = []
            if os.path.exists(subdirpath + "/TimeSeries_vals.txt"):
                with open(subdirpath + "/TimeSeries_vals.txt", "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        T_vals.append(float(line.strip()))
                f.close()
            else:
                #Using find_vals to extract T values from files.
                find_vals(files, [r'T_[\d]*[.][\d]+', r'T_[\d]+'], T_vals)
                # Remove T_ from T_vals and sort them in ascending order.
                T_vals = sorted([float(re.sub( "T_", "", v)) for v in T_vals])
                # Only consider unique values of T.
                T_vals =np.unique(T_vals)
                # For each value of T, find Rmax of corresponding files in files.
            for t in T_vals:
                t = int(t) if float(t).is_integer() else t
                selectedfiles = glob.glob(subdirpath + "/MOVSERIES*_T_" + str(t) + "*.csv")

                # Using find_vals to extract R values from files.
                movmaxR = max([int(re.findall(r'[\d]+' , f)[-1]) for f in selectedfiles])

                # Save movmaxR to a file in subdirpath.
                try:
                    with open(subdirpath + "/movmaxR_T_" + str(t) + ".txt", "w") as f:
                        f.write(str(movmaxR+1))
                    f.close()
                except Exception as e:
                    print("Error: Could not write movmaxR to */" + os.path.basename(subdirpath) 
                            + "/movmaxR_T_" + str(t) + ".txt with error message: \n" + str(e))
                
                # Elementwise average, variance and surviving runs for each value of T in subdirpath using gen_MEAN_INDVL_Prelimsfiledata()
                df_avgtimeseries = gen_MEAN_INDVL_Prelimsfiledata(selectedfiles, subdirpath, "csv", tmax =t, dt= dt, include_col_labels= ["t", "<GAM[P(x; t)]>_x", "<vx[P(x; t)]>_x", "<vy[P(x; t)]>_x",
                                "<GAM[G(x; t)]>_x", "<vx[G(x; t)]>_x", "<vy[G(x; t)]>_x", "<GAM[Pr(x; t)]>_x", "<vx[Pr(x; t)]>_x", "<vy[Pr(x; t)]>_x"] , 
                                handle_nonstandardtime_fileconflicts= "interpolate")
                try:
                    if df_avgtimeseries is not None:
                        df_avgtimeseries.to_csv(subdirpath + f"/MEAN_MOVSERIES_T_{t}.csv", index=False, header=True)
                except Exception as e:
                    print("Error: Could not write MEAN_MOVSERIES to " + subdirpath + f"/MEAN_MOVSERIES_T_{t}.csv with error message: \n" + str(e))
        print(f"Done post-processing MOVSERIES for prefix: {pre}...")
        # NOW DEALING WITH TSERIES FILES.
        # Recursively navigate to each sub-directory within out_dir_noprefix/{pre} that contains files of the form "TSERIES*.csv".
        # These subdirectories may be immediate or nested.
        # For each such non-empty nested sub-directory, use glob to find all files of the form "TSERIES*.csv".
        # Pass this list of files to other functions for further processing.
        
        all_valid_subdirpaths = [dir for dir, _, files in os.walk(workdir) if len(glob.glob(dir + "/TSERIES*.csv")) > 0]
        for subdirpath in all_valid_subdirpaths:
            files = glob.glob(subdirpath + "/TSERIES*.csv")

            # FIRST, getting maximum replicate number (maxR) for each value of T in subdirpath.
            # Read list of T values from TimeSeries_vals.txt in subdirpath (if it exists), else create it using the files.
            T_vals = []
            if os.path.exists(subdirpath + "/TimeSeries_vals.txt"):
                with open(subdirpath + "/TimeSeries_vals.txt", "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        T_vals.append(float(line.strip()))
                f.close()
            else:
                #Using find_vals to extract T values from files.
                find_vals(files, [r'T_[\d]*[.][\d]+', r'T_[\d]+'], T_vals)
                # Remove T_ from T_vals and sort them in ascending order.
                T_vals = sorted([float(re.sub( "T_", "", v)) for v in T_vals])
                # Only consider unique values of T.
                T_vals =np.unique(T_vals)
            ''' Extend files with 0s at the end of "<<P(x; t)>_x>_r" for higher values of T_vals.
            # Do this using gen_missing_zero_Prelims(...)
            if (gen_missing_zero_Prelims(files, T_vals, subdirpath, "csv", 
                                          dt= dt, compare_col_labels= ["<<P(x; t)>_x>_r"]) == None):
                print("Error: Could not extend files with zeros for higher values of T in " + subdirpath)
            #'''
            # For each value of T, find Rmax of corresponding files in files.
            for t in T_vals:
                t = int(t) if float(t).is_integer() else t
                selectedfiles = glob.glob(subdirpath + "/TSERIES*_T_" + str(t) + "*.csv")

                # Using find_vals to extract R values from files.
                maxR = max([int(re.findall(r'[\d]+' , f)[-1]) for f in selectedfiles])
            
                # Save maxR to a file in subdirpath.
                try:
                    with open(subdirpath + "/maxR_T_" + str(t) + ".txt", "w") as f:
                        f.write(str(maxR+1))
                    f.close()
                except Exception as e:
                    print("Error: Could not write maxR to */" + os.path.basename(subdirpath) + "/maxR_T_" + str(t) + ".txt with error message: \n" + str(e))

                

                # Elementwise average, variance and surviving runs for each value of T in subdirpath using gen_MEAN_INDVL_Prelimsfiledata()

                df_avgtimeseries = gen_MEAN_INDVL_Prelimsfiledata(selectedfiles, subdirpath, "csv", tmax =t, dt= dt, handle_nonstandardtime_fileconflicts= "interpolate")
                try:
                    if df_avgtimeseries is not None:
                        df_avgtimeseries.to_csv(subdirpath + f"/MEAN_TSERIES_T_{t}.csv", index=False, header=True)
                except Exception as e:
                    print("Error: Could not write MEAN_TSERIES to " + subdirpath + f"/MEAN_TSERIES_T_{t}.csv with error message: \n" + str(e))
        print("\n=====================================================================================================\n") 
        print(f"Done post-processing all prelims for prefix: {pre}...")
        print("\n=====================================================================================================\n")


def main():
    
    #Find a list of immediate subdirectories in root_dir
    subdir = [os.path.basename(subdir) for subdir in glob.glob(os.path.join(root_dir, '*/'))]
    
    for pre in prefixes:

        out_dir = out_dir_noprefix;
        if pre == "":
            out_dir = out_dir_noprefix + "NoPrefix/"
        else:
            out_dir = out_dir_noprefix + pre + "/"
        
        #Subdirectories are of the form "PREFIX*_dP_{dP}_Geq_{val}"
        #Using regex and glob, create a sorted list of subdir (in ascending order of val (can be float or integer)) with prefix pre.
        
        if (Geq == "NA" and Veq == "NA"):
            subdir = sorted(glob.glob(os.path.join(root_dir, pre + "_*_dP_" + str(dP))), 
                            key=lambda x: float(re.findall(r'[\d]*[.][\d]+', x)[0]) if re.findall(r'[\d]*[.][\d]+', x) else int(re.findall(r'[\d]+', x)[-1]))
        else:

            if(Geq != "NA"):
                subdir = sorted(glob.glob(os.path.join(root_dir, pre + "_*_dP_" + str(dP) + "_Geq_" + str(Geq))), 
                        key=lambda x: float(re.findall(r'[\d]*[.][\d]+', x)[0]) if re.findall(r'[\d]*[.][\d]+', x) else int(re.findall(r'[\d]+', x)[-1]))
            elif (Veq != "NA"):
                subdir = sorted(glob.glob(os.path.join(root_dir, pre + "_*_dP_" + str(dP) + "_Veq_" + str(Veq))), 
                        key=lambda x: float(re.findall(r'[\d]*[.][\d]+', x)[0]) if re.findall(r'[\d]*[.][\d]+', x) else int(re.findall(r'[\d]+', x)[-1]))
        
        subdir_test = [os.path.basename(s) for s in subdir]
        print("Found subdirectories: " + str(subdir_test))
        #print(os.path.basename(subdir[0])[0])
        
        if pre == "":
            #If pre is empty, valid subdirectories start with a digit.
            subdir = [sub for sub in subdir if str(os.path.basename(sub))[0].isdigit()]

        print("Found subdirectories: " + str(subdir))
        if dynamic_inspect:
            input("Press Enter to continue...")
    
        for g in L:

            # In each selected subdir, all files are  of the form
            # "TimeSeries/PRELIM_TSERIES_*_G_{g}_T_{T}_dt_{dt}_a_{a_val}_*_R_{R}_*.csv"

            # First make a sorted list of a_val in each subdir.

            for sub in subdir:
                #Get all files in subdirectory Timeseries that begin with the above, but not e.
                files = glob.glob(sub + "/TimeSeries/PRELIM_*_G_" + str(g) + "_T_*.csv")
                #Find all unique values of a_val in files (in ascending order) (using regex with finding "a_{a_val}" and then removing a_).
                a_vals = []
                find_vals(files, [r'a_[\d]*[.][\d]+', r'a_[\d]+'], a_vals)
                #print("For subdir: " + sub)
                #print("Found a_vals: " + str(a_vals))
                print(len(a_vals))
                # Remove a_ from a_vals and sort them in ascending order.
                a_vals = sorted([float(re.sub( "a_", "", v)) for v in a_vals])
                # Only consider unique values of a.
                a_vals =np.unique(a_vals)

                print("For subdir: " + sub)
                print("Found a_vals: " + str(a_vals))

                # Check if out_dir exists in savedir.txt in sub.
                # If it does, provide a warning followed by a prompt to skip ("skip"/1) or continue ("continue"/2).

                
                if os.path.exists(sub + "/savedir.txt"):
                    with open(sub + "/savedir.txt", "r") as f:
                        lines = f.readlines()
                        for line in lines:
                            line.strip("\n")
                            if out_dir + "/L_" +str(g) + "_*" in line:
                                print("WARNING!!!: " + out_dir + " already exists in " + sub + "/savedir.txt")
                                if dynamic_inspect:
                                    print("Do you want to skip this directory ('skip'/1) or continue overwriting ('continue'/2)?")
                                    choice = input("Enter choice: ")
                                    if choice == "skip" or choice == "1":
                                        continue
                                    elif choice == "continue" or choice == "2":
                                        break
                                    else:
                                        print("Invalid choice. Skipping directory.")
                                        continue
                                else:
                                    print("Skipping directory.")
                                    continue
                    f.close()
                                
                for a in a_vals:
                    a = int(a) if float(a).is_integer() else a
                    savedir = "L_" + str(g) + "_a_" + str(a) + "/dP_" + str(dP) + "/Geq_"# + str(Geq) + "/TimeSeries/"
                    if Geq == "NA" and Veq != "NA":
                        savedir += str(Veq) + "/TimeSeries/"
                    else:
                        savedir += str(Geq) + "/TimeSeries/"
                    outdir = out_dir + savedir
                    # Recursively create outdir (and parent directories) if it doesn't exist.
                    Path(outdir).mkdir(parents=True, exist_ok=True)

                    t_range = find_timerange(files, a, indx_vals_t, tmin, tmax)

                    # Finding jID values in files which contain a_{a} in their name.
                    jiD_range= [ re.findall(r'jID_[\d]+' , f)[0] if re.findall(r'jID_[\d]+' , f) else "jID_-1" for f in files
                                if (re.findall(r'a_[\d]*[.][\d]+' , f) and  re.findall(r'a_[\d]*[.][\d]+' , f)[0] == "a_" + str(a))
                                or (re.findall(r'a_[\d]+' , f) and  re.findall(r'a_[\d]+' , f)[0] == "a_" + str(a))]
                    # Extract jiD values from jiD_range.
                    jiD_vals = sorted([int(re.sub("jID_", "", j)) for j in jiD_range]); 
                    # Remove duplicates from jiD_vals, and append -1 at the end to indicate possible frames with no jiD data.
                    jiD_vals = np.unique(jiD_vals); jiD_vals = np.append(jiD_vals, -1)

                    print("Processing: " + savedir)

                    for jiD in jiD_vals:    # Loop over jiD values first.
                        for t in t_range:
                            #Get all files in parent directory that begin with the above (but not subdirectories).
                            if jiD != -1:
                                selectedfiles = glob.glob(sub + "/TimeSeries/PRELIM_*_G_" + str(g) + "_T_" + str(t) +"*_a_" + str(a) + "*_jID_" + str(jiD) + "*.csv")
                            else:
                                # RECALL- A file with jiD = -1 indicates that the file does not contain unique jiD data.
                                selectedfiles = [f for f in  glob.glob(sub + "/TimeSeries/PRELIM_*_G_" + str(g) + "_T_" + str(t) +"*_a_" + str(a) + "_*.csv") 
                                                 if "jID" not in f]
                            # If selectedfiles is empty, skip the loop over t.
                            if len(selectedfiles) == 0:
                                continue
                            t_subdir = outdir
                            # Recursively create t_subdir (and parent irectories) if it doesn't exist.
                            #Path(t_subdir).mkdir(parents=True, exist_ok=True)
                            print("Processing t= " + str(t))
                            
                            # Copy all files in files to t_subdir, with new names and avoiding conflicts.
                            # If a file with the same name exists in t_subdir, rename the file to avoid conflicts.
                            # This renaming will be done by sending the conflicting file name to a function that will rename it.

                            for source_file in selectedfiles:
                                source_filename = os.path.basename(source_file)
                                # Change source_filename as follows:
                                # Source file name is of the form "PRELIM_TSERIES*_G_{g}_T_{T}_dt_{dt}_a_{a_val}_*_R_{R}_*.csv"
                                # Change it to "FRAME_T_{T}_a_{a_val}_R_{R}.csv"

                                # First extract T, a_val, R from source_filename.
                                T = t; a_val = a; R = int(re.findall(r'[\d]+' , source_filename)[-1])
                                # Find filetype of source_filename, which are the characters before the first "_"
                                filetype = source_filename.split("_")[1]
                                # Now create new filename.
                                out_filename = filetype + "_T_" + str(T) + "_a_" + str(a_val) + "_R_" + str(R) + ".csv"
                                '''
                                if re.search("GAMMA", source_filename):
                                    out_filename = "GAMMA_T_%g_a_%g" %(T, a_val) + "_R_" + str(R) + ".csv"
                                else:
                                    out_filename = "FRAME_T_%g_a_%g" %(T, a_val) + "_R_" + str(R) + ".csv"
                                '''
                                outfile = t_subdir + out_filename
                                #Check if outfile already exists.
                                #fd_source = os.open(source_file, os.O_RDONLY);
                                #print("Status of source file: " + str(os.fstat(fd_source)))
                                if os.path.exists(outfile):
                                    print("Warning: File " + out_filename + " already exists, and conflicts with source file: "
                                        + source_filename + ". Renaming file.")
                                    # First check if outfile is the same as source_file (this is done by looking at the attributes of the files (filesize etc)).
                                    # If they are the same, skip the file.
                                    for out_f in glob.glob(t_subdir + "*.csv"):

                                        if os.path.samefile(source_file, out_f):
                                            print("WARNING!!!: File " + outfile + " is the same as " + source_file + ". Skipping file.")
                                            continue

                                        fd_source = os.open(source_file, os.O_RDONLY);
                                        fd_outfile = os.open(out_f, os.O_RDONLY);
                                        if os.fstat(fd_source).st_ino == os.fstat(fd_outfile).st_ino:
                                            print("WARNING!!!: File " + outfile + " is the same as " + source_file + ". Skipping file.")
                                            continue

                                        os.close(fd_source)
                                        os.close(fd_outfile)
                                        
                                    # If outfile already exists, rename it to avoid conflicts.
                                    max_R = max([int(re.findall(r'[\d]+' , f)[-1]) for f in glob.glob(t_subdir + filetype + "_T_" + str(T) + "_a_" + str(a_val) + "*.csv")])
                                    out_filename = filetype + "_T_" + str(T) + "_a_" + str(a_val) + "_R_" + str(max_R + 1) + ".csv"
                                    outfile = t_subdir + out_filename

                                # Copy source_file to outfile.
                                try:
                                    shutil.copy2(source_file, outfile)
                                except shutil.Error as e:
                                    print("Error: Could not copy file " + source_file + "\tto:\n" + outfile + "with error message: \n" + str(e))
                                    print("Skipping file.")
                                    continue
                                except FileNotFoundError as e:
                                    print("Error: Could not find file " + source_file + "  with error message: \n" + str(e))
                                    print("Skipping file.")
                                    continue
                            # This is the end of the loop over source_file in selectedfiles.

                        # This is the end of the loop over t in t_range.

                    #Create a timeseries_vals.txt file by looking at all files matching outdir + "*.csv" and extracting the values of T.
                    test_txt_RW(outdir, "*.csv", "T", ["[\d]*[.][\d]+", "[\d]+"], savefile= "TimeSeries_vals.txt")

                # This is the end of the loop over a in a_vals.
                # Now appending list of savedir in a txt file in sub.
                
                #Open a file in sub (in append mode), and append savedir to it.
                try:
                    with open(sub + "/savedir.txt", "a") as f:
                        f.write(out_dir + "L_" + str(g) + "_*" +"\n")
                    f.close()
                except Exception as e:
                    print("Error: Could not write savedir to " + sub + "/savedir.txt with error message: \n" + str(e))
                
        print("Finished processing prefix: " + pre)
        

        # Check if out_dir exists.
        if os.path.exists(out_dir):
            print("Adding a-val.txt to " + out_dir)
            test_txt_RW(out_dir, "*/", "a", ["[\d]*[.][\d]+", "[\d]+"])

        # Read a_vals.txt in out_dir and print the values.
        if os.path.exists(out_dir + "/a_vals.txt"):
            print("Reading a_vals.txt in " + out_dir)
            with open(out_dir + "/a_vals.txt", "r") as f:
                lines = f.readlines()
                for line in lines:
                    print(line.strip())
            f.close()
    
    # This is the end of the loop over pre in prefixes.
    # Now post-process the data files in out_dir_noprefix.
    '''
    print("\n\n__________________________________________________________________________________________\n\n")
    print("Now, Post-processing data files in " + out_dir_noprefix)
    print("\n\n__________________________________________________________________________________________\n\n")
    post_process(prefixes)
    '''

set_prelims_inputs()



print("Using CPU for processing...")
# Set up a timer to measure the time taken for processing.
start_time_CPU = time.time()
if __name__ == "__main__":
    main()
    post_process(prefixes)

end_time_CPU = time.time()
elapsed_time_CPU = end_time_CPU - start_time_CPU
print(f"Time reported by CPU: {elapsed_time_CPU:.3f} seconds")