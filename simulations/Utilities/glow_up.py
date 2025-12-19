import os
import traceback
import regex as re
from pathlib import Path
import glob
import time
import shutil
import psutil
import sys
import argparse
import copy
import warnings
#from multiprocessing import Pool, cpu_count, Queue, Lock, get_start_method
import functools
from functools import lru_cache
import secrets
from collections import defaultdict, Counter
from contextlib import contextmanager



print(f"GLOW_UP: PID={os.getpid()}, __name__={__name__}, importing...")
#print(f"STACK TRACE for PID {os.getpid()}:")
#traceback.print_stack()
#print("=" * 50)

import slick as gpu
import threading
from joblib import Parallel, delayed, parallel_backend, cpu_count, Memory

if gpu.DASK_AVAILABLE:
    # Set up global Dask Client and configure it.
    #dask = gpu._dask
    #dask_distributed = gpu._dask_distributed
    dask_progress = gpu._dask_progress
    dask_delayed = gpu._dask_delayed
    dask_Client = gpu._dask_Client
    dask_LocalCluster = gpu._dask_LocalCluster
    #dask_ProgressBar = gpu._dask_ProgressBar

# Importing necessary libraries
np = gpu.np
signal = gpu.signal
fft = gpu.fft
ndimage = gpu.ndimage

interpolate = gpu.interpolate
stats = gpu.stats
asarray = gpu.asarray
asnumpy = gpu.asnumpy
is_gpu_array = gpu.is_gpu_array

# For type checking and constants - use direct CPU access
_np_base = gpu._cpu_np  # Access the underlying numpy module
_pan_base = gpu._cpu_pandas  # Access the underlying pandas module

# CPU-only modules
esda_moran = gpu.cpu_safe_import("esda.moran")
libpysal_weights = gpu.cpu_safe_import("libpysal.weights")
#from esda_moran import Moran_BV
#from libpysal_weights import Queen, KNN, lat2W

# Import specific interpolation functions
#from interpolate import InterpolatedUnivariateSpline, CubicSpline



#import pandas as pan
if gpu.GPU_AVAILABLE and not gpu.RAPIDS_AVAILABLE:
    # If GPU acceleration is available but RAPIDS is not, wrap imports to be cpu-safe.
    sklearn_cluster = gpu.cpu_safe_import("sklearn.cluster")
    sklearn_mixture = gpu.cpu_safe_import("sklearn.mixture")
    skimage_metrics = gpu.cpu_safe_import("skimage.metrics")
    sklearn_metrics = gpu.cpu_safe_import("sklearn.metrics")
    pandas = gpu.cpu_safe_import("pandas"); pan= pandas
    KMeans = sklearn_cluster.KMeans
    GaussianMixture = sklearn_mixture.GaussianMixture
    ssim = skimage_metrics.structural_similarity
    mutual_info_score = sklearn_metrics.mutual_info_score
    normalized_mutual_info_score = sklearn_metrics.normalized_mutual_info_score
    adjusted_mutual_info_score = sklearn_metrics.adjusted_mutual_info_score
else:
    # RAPIDs handles CPU fallback by default, so no need for additional handling here.
    from sklearn.cluster import KMeans
    #from sklearn.mixture import GaussianMixture
    #from skimage.metrics import structural_similarity as ssim
    #from sklearn.metrics import mutual_info_score #, normalized_mutual_info_score
    sklearn_metrics = gpu.cpu_safe_import("sklearn.metrics")
    skimage_metrics = gpu.cpu_safe_import("skimage.metrics")
    sklearn_mixture = gpu.cpu_safe_import("sklearn.mixture")
    GaussianMixture = sklearn_mixture.GaussianMixture
    ssim = skimage_metrics.structural_similarity
    mutual_info_score = sklearn_metrics.mutual_info_score
    normalized_mutual_info_score = sklearn_metrics.normalized_mutual_info_score
    adjusted_mutual_info_score = sklearn_metrics.adjusted_mutual_info_score
    #from sklearn.metrics import adjusted_mutual_info_score
    import pandas as pan

'''
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

import scipy.stats as stats
import scipy.fft as fft
from scipy.signal import find_peaks
from scipy.signal import correlate2d
from scipy.signal import fftconvolve
from scipy.interpolate import CubicSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from skimage.metrics import structural_similarity as ssim
from sklearn.metrics import mutual_info_score , normalized_mutual_info_score
from sklearn.metrics import adjusted_mutual_info_score
from esda.moran import Moran_BV
from libpysal.weights import Queen, KNN, lat2W
from collections import defaultdict

from multiprocessing import Pool, cpu_count, Queue, Lock
import functools


import secrets
'''
rng = np.random.default_rng(secrets.randbits(128))

#print_lock = Lock()

'''
This script contains a set of utility functions that can be used to read and write text files, read and write csv files,
and generally manipulate/reorganise data in given directories. It is intended to be invoked from reorganise_dir.py and other scripts
in the Utilities directory.

A detailed description of each function is provided below.
'''


@contextmanager
def managed_pool(pool_processes=None):
    from multiprocessing import Pool, cpu_count, Queue, Lock, get_start_method
    """Context manager for safe pool handling"""
    if pool_processes is None:
        pool_processes = cpu_count()
    
    pool = None
    try:
        pool = Pool(processes=pool_processes)
        yield pool
    except Exception as e:
        print(f"Pool creation failed: {e}")
        yield None
    finally:
        if pool is not None:
            pool.close()
            pool.join()

class SkipFile(Exception):
    pass

# If you need to check for both CPU and GPU arrays:
def is_array(obj):
    """Check if object is either a CPU numpy array or GPU cupy array"""
    if isinstance(obj, _np_base.ndarray):
        return True
    if gpu.GPU_AVAILABLE and gpu.is_gpu_array(obj):
        return True
    return False


# Multi-thread safe-print function using lock
'''
def safe_print(*args, **kwargs):
    print_lock.acquire()
    try:
        print(*args, **kwargs)
    finally:
        print_lock.release()
'''

# Replacing Lambda functions for defaultdict initialisation to avoid pickling issues with multiprocessing
def default_dict():
    return defaultdict(dict)

def nested_default_dict():
    return defaultdict(default_dict)

def double_nested_default_dict():
    return defaultdict(nested_default_dict)

def recursively_nested_default_dict(iter):
    if iter == 0:
        return defaultdict(dict)
    else:
        return defaultdict(recursively_nested_default_dict(iter-1))
    
def sorted_files_R(files,  ascending=True):
    sorted_files = sorted(files, key=lambda x: int(re.sub("R_", "", re.findall(r'R_[\d]+', x)[0])))
    return sorted_files if ascending else sorted_files[::-1]


#-------------------------- Defining cached functions for Dask performance --------------------------
# Libpysal weight functions are computationally expensive, so we cache them for the last 8 L vals to avoid recomputation.
@lru_cache(maxsize=16)
def cached_libpysalweights_lat2W(L):
    return libpysal_weights.lat2W(L, L)

# Set up cache directory for joblib Memory.
pysal_cache_dir = os.path.expanduser("~/.cache/libpysal_weights")
os.makedirs(pysal_cache_dir, exist_ok=True)
print(f"Using libpysal weights cache directory: {pysal_cache_dir}...."); time.sleep(1)

memory = Memory(location=pysal_cache_dir, verbose=1)

#@memory.cache
#def cached_libpysalweights_lat2W(L):
#    return libpysal_weights.lat2W(L, L)

def clear_libpysal_cache():
    print(f"[libpysal cache] Clearing cache at {pysal_cache_dir} ...")
    memory.clear(warn=False)

#-------------------------- End of cached functions for Dask performance --------------------------

'''# Summary of the function find_timerange(...)
    This function finds all T values in provided "files" (a list of files that may or may not be from the same directory)
    that have a_val = a. (t-values can be float or integer). It then extracts the "indx" smallest or largest values of T 
    ( if indx is negative/ positive respectively) from the list of T values.
    It then returns the list of T values that are in the range [min_t, max_t] (if min_t and max_t are provided).
'''
def find_timerange(files, a, indx, min_t = None, max_t = None):
    # a is the value of a in the above.
    # indx provides the number of t-values to extract from the file name.
    # Negative values of indx (say -n) indicate the n largest values of T are to be extracted.
    # Positive values of indx (say n) indicate the n smallest values of T are to be extracted.
    # Return a list of T values in files.

    # First, extract all T values from files that have a_val = a. (t-values can be float or integer).
    T_vals = [re.findall(r'T_[\d]*[.][\d]+' , f)[0] if re.findall(r'T_[\d]*[.][\d]+' , f) 
              else re.findall(r'T_[\d]+' , f)[0] for f in files 
              if (re.findall(r'a_[\d]*[.][\d]+' , f) and  re.findall(r'a_[\d]*[.][\d]+' , f)[0] == "a_" + str(a) )
              or (re.findall(r'a_[\d]+' , f) and  re.findall(r'a_[\d]+' , f)[0] == "a_" + str(a) )]
    # Remove T_ from T_vals and sort them in ascending order.
    T_vals = sorted([float(re.findall(r'[\d]*[.][\d]+', t)[0]) 
                     if re.findall(r'[\d]*[.][\d]+', t) 
                     else int(re.findall(r'[\d]+', t)[0]) for t in T_vals])
    # Only consider unique values of T.
    T_vals = sorted(list(set(T_vals)))
    # If min_t and max_t are provided, extract T values in the range [min_t, max_t].

    if np.abs(indx) > len(T_vals):
        print("Warning: indx greater than number of T values. Returning all T values.")
        return T_vals
    
    if min_t is not None:
        T_vals = [t for t in T_vals if t >= min_t]
    if max_t is not None:
        T_vals = [t for t in T_vals if t <= max_t]
    
    if np.abs(indx) > len(T_vals):
        print("Warning: indx greater than number of T values (as per constraints). Returning all T values (as per constraints).")
        return T_vals
    return T_vals[indx:] if indx < 0 else T_vals[:indx]



''' Summary of the function freedman_diaconis_bin_width(...): Compute Freedman-Diaconis bin width.'''
def freedman_diaconis_bin_width(data):
    """Compute Freedman-Diaconis bin width."""
    data = data[np.isfinite(data)]  # remove NaNs/infs
    q75, q25 = np.percentile(data, [75, 25])
    IQR = q75 - q25; n = len(data)
    if IQR == 0 or n == 0:
        return np.nan
    return 2 * IQR / np.cbrt(n)

"""Summary of the function scotts_bin_width(...): Compute Scott's bin width."""
def scotts_bin_width(data):
    #print("Using Scott's bin width...")
    #print(f"Data type: {type(data)}")
    #print(f"NP finite check: {np.isfinite(data).all()}")
    mask = np.isfinite(data)
    #print(f"IS CUPY instance data: {isinstance(data, gpu._cupy.ndarray)}")
    #print(f"Is CUPY instance mask: {isinstance(mask, gpu._cupy.ndarray)}")
    data = data[np.isfinite(data)]  # remove NaNs/infs
    std = np.std(data)
    #print(f"Calculated standard deviation: {std}")
    n = len(data)
    #print(f"Number of data points: {data.size}")
    #print(f"Calculated number of data points: {n}")
    if std == 0 or n == 0:
        return np.nan
    #print(f"Calculated Scott's bin width: {3.5 * std / np.cbrt(n)}")
    return 3.5 * std / np.cbrt(n)


'''# Summary of the function gen_interpolated_df(...)
gen_interpolated_df(df1, df2, compare_col) generates an interpolated dataframe of df2 with respect to df1, based on the values in compare_col_label.
The function first checks if compare_col_label is an existing column in both df1 and df2. If it is not, the function will return an error (None).

Next the function will do the follows:
1.) It will create a Pandas dataframe with the same columns as df2 called df_interpolate, with df_interpolate[compare_col_label] = df1[compare_col_label].
1.) It will iterate over all the columns in df2, barring compare_col_label, and will 
generate a cubic spline interpolation of each such column in df2 with respect to compare_col_label.
2.) It will use the cubic spline interpolation to evaluate the interpolated values of each column in df2 
at the values of compare_col_label in df1 (stored in df_interpolate[compare_col_label]).
3.) These interpolated values will be stored in the corresponding columns in df_interpolate.
4.) The function will return df_interpolate.
5.) If the function encounters an error while generating the cubic spline interpolation, it will return None.
'''
def gen_interpolated_df(df1, df2, compare_col_label):
    
    # Check if compare_col_label is an existing column in both df1 and df2. If it is not, return an error.
    if compare_col_label not in df1.columns or compare_col_label not in df2.columns:
        return None
    # Create a Pandas dataframe with the same columns as df2 called df_interpolate.
    df_interpolate = pan.DataFrame(columns= df2.columns)
    # Set df_interpolate[compare_col_label] = df1[compare_col_label].
    df_interpolate[compare_col_label] = df1[compare_col_label]
    interpolate_cols = [col for col in df2.columns if col != compare_col_label]
    # Iterate over all the columns in df2, barring compare_col_label.
    for col in interpolate_cols:
        # Generate a cubic spline interpolation of each such column in df2 with respect to compare_col_label.
        try:
            cs = interpolate.CubicSpline(df2[compare_col_label], df2[col])
        except Exception as e:
            print("Error: Could not generate cubic spline interpolation for column " + col + " with error message: \n" + str(e))
            return None
        # Evaluate the interpolated values of each column in df2 at the values of compare_col_label in df1.
        df_interpolate[col] = cs(df1[compare_col_label])
    return df_interpolate

'''# Summary of the function gen_zer0mean_normalised_df(...)
    This function generates a COL-SPECIFIC zero-mean normalised dataframe of df, excluding the columns in exclude_col_labels.
    The function will return the zero-mean normalised dataframe.
'''
def gen_zer0mean_normalised_df(df, exclude_col_labels= ["a_c", "x", "L"]):
    # Exclude the columns in exclude_col_labels from df.
    df = df[[col for col in df.columns if col not in exclude_col_labels]]
    # Generate zero-mean normalised columns for each column in df.
    df = df.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10) , axis=0)
    return df

def find_vals(files, reg_search_patterns, vals, dir=""):
    # Recursively iterate over reg_search_patterns and find all values in files that match the pattern.
    # Append the values to vals.
    files_copy = copy.deepcopy(files)
    if len(reg_search_patterns) == 0:
        return
    else:
        for f in files:
            vals += [re.findall(reg_search_patterns[0], f)[0] for f in files_copy if re.findall(reg_search_patterns[0], f)]
        # Remove matching values from files.
        files_copy = [f for f in files_copy if not re.findall(reg_search_patterns[0], f)]
        find_vals(files_copy, reg_search_patterns[1:], vals, dir)

''' # Summary of the function savetxt_to_dir(...)
    This function saves a dataframe or an array-like-object to a directory (savedir), with name savefilename.
    If the directory does not exist, it creates the directory. If the file already exists in the directory, 
    it renames the file to avoid conflicts.
'''
def savetxt_to_dir(data, savedir, savefilename, header = None, index = False, sep = "\t", overwrite = True, verbose = False):
    
    if not os.path.exists(savedir):
        Path(savedir).mkdir(parents=True, exist_ok=True)
    # Check if the file already exists in the directory.
    if os.path.exists(savedir + "/" + savefilename):
        # Provide a warning before overwriting the file.
        if not overwrite:
            print("Warning: File " + savedir + "/" + savefilename + " already exists. Skipping file.")
            return
    try:
        # Check if data is a dataframe or an array (or list-like object).
        if isinstance(data, _pan_base.DataFrame):
            data.to_csv(savedir + "/" + savefilename, header = header, index = index, sep = sep)
        else:
            np.savetxt(savedir + "/" + savefilename, data, delimiter = sep, fmt="%g")

        if verbose:
            with open(savedir + "/" + savefilename, "r") as f:
                lines = f.readlines()
                for line in lines:
                    print(line.strip())
                # Close the file.
                f.close()
    except Exception as e:
        print("Error: Could not write vals to " + savedir + "/" + savefilename + " with error message: \n" + str(e))
        return

''' # Summary of the function test_txt_RW(...)
    The function test_txt_RW(...) are used to read and write text files from a directory (parendir) by iterating over 
    some set of files in parendir. BY DEFAULT  test_txt_RW() iterates over all immediate subdirectories in parendir, but 
    can be used to iterate over other sets of files as well in parendir (BY CHANGING file_type).
    test_txt_RW(...) is used to read and write text files from a directory (parendir).
    If  set_of_files is provided, the function will iterate over the files in set_of_files, else it will iterate over the files in parendir as specified by file_type.
    Both functions extract all vals in files that match regex_search_patterns, with type_id as identifying prefix (if not empty)
    (such as "a" or "L" or "T" etc.) and save them to a text file in the same directory (tab delineated), with the name "savefile".
'''

def test_txt_RW(parendir, file_type="*/", type_id ="", reg_search_patterns = ["[\d]*[.][\d]+", "[\d]+"], 
                set_of_files= None,  savefile = "default", verbose= False):

    # Extract all vals in dir that match regex_search_patterns, with type_id as identifying prefix (if not empty) (such as "a" or "L" or "T" etc.)
    # and save them to a text file in the same directory (tab delineated), with the name "savefile"

    # First make a list of all immediate subdirectories in dir.
    if set_of_files is not None:
        files = set_of_files
    else:
        files = glob.glob(os.path.join(parendir, file_type))
    if type_id != "":
        reg_search_type_patterns = [type_id + "_" + reg for reg in reg_search_patterns]
    else:
        reg_search_type_patterns = reg_search_patterns
    # Recursively iterate over regex_search_patterns and find all values in files that match the pattern.
    vals=[]
    find_vals(files, reg_search_type_patterns, vals, parendir)
    #print(f"Found {type_id} vals: " + str(vals))
    # Remove type_id from vals and sort them in ascending order if type_id is not empty.
    vals = sorted([float(re.sub(type_id + "_", "", v)) for v in vals])
    vals = np.unique(vals)
    print(f"Found {type_id} vals: " + str(vals))

    if savefile == "default":
        savefile = type_id + "_vals.txt"
    savetxt_to_dir(vals, parendir, savefile, header = None, index = False, sep = "\t", overwrite = True, verbose = verbose)

''' # Summary of the function gen_MEAN_SD_COLSfiledata(...)
# This function generates the mean and standard deviation of each of the columns in the files.
# excluding exclude_col_labels. Files are standardised to have the same names of columns (but not necessarily in the same order).
# If the files have non-standardised columns (different names other than exclude_col_labels), the function will return an error.
#If nonzero is True, averages and means are calculated over only non-zero columns across all files, else over all columns.
# If add_counts is True, the function will also add a column to the dataframe next to the mean and standard deviation columns for each column-labels
# that contains the number of non-zero columns for that column-label across all files. If nonzero is False, the added columns will simply contain the number of columns
# for that column-label across all files (which will be the same as the number of files).
'''

def gen_MEAN_SD_COLSfiledata(files, pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L"], nonzero=True, add_counts= True):

    # Read the first file in files.
    if(ext == "csv"):
        df = pan.read_csv(files[0], header=0)
    else:
        try:
            df = pan.read_table(files[0], header=0)
        except Exception as e:
            print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
            return None
    
    # Store column names in col_labels. Strip elements of any leading or trailing whitespaces.
    col_labels = df.columns; col_labels = [col.strip() for col in col_labels]
    col_labels = [col for col in col_labels if col not in exclude_col_labels];
    # Remove exclude_col_labels from col_labels.

    if(add_counts):
        col_names = [ y + x + "]" for x in col_labels  for y in ["AVG[", "VAR[", "COUNT["]]
    else:
        col_names = [ y + x + "]" for x in col_labels  for y in ["AVG[", "VAR["] ]
    
    # Create a df with col_names as column names.
    df_out = pan.DataFrame(columns= col_names)
    # Pooled Data Stats will store the pooled average, standard deviation, sample size and counts for each column in col_labels across all files.
    # It will be iteratively updated as the function iterates over all files in files.
    pooled_data_stats = pan.DataFrame(columns = [y + x for x in col_labels 
            for y in ["POOLED_AVG ", "POOLED_VAR ", "POOLED_SAMPLE_SIZE ", "POOLED_COUNTS "]])
    # Initialise the pooled_data_stats to 0.
    pooled_data_stats.loc[0] = 0.0
    # Iterate over all files in files.
    for file in files:

        # Read the file.
        if(ext == "csv"):
            df = pan.read_csv(file, header=0)
        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
                return None
        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        df = df[[col for col in df.columns if col not in exclude_col_labels]]
        
        # First check if the columns in df are the same as col_labels. If not, return an error.
        if not all([col in df.columns for col in col_labels]):
            print("Error: Columns in file " + pathtodir + "/" + file + " are not standardised.")
            return None
        
        for col in col_labels:

            # Find mean, variance and sample size of each column in df.
            mean_col = df[col].mean()
            if nonzero:
                if mean_col == 0:
                    continue

            variance_col = df[col].var()
            # Store sample size in ss
            ss= len(df[col])

            # Update Pooled Data Stats iteratively.
            pooled_data_stats.loc[0, "POOLED_AVG " + col] = (pooled_data_stats.loc[0, "POOLED_AVG " + col] * pooled_data_stats.loc[0, "POOLED_SAMPLE_SIZE " + col] 
                    + mean_col * ss)/(pooled_data_stats.loc[0, "POOLED_SAMPLE_SIZE " + col] + ss)
            pooled_data_stats.loc[0, "POOLED_VAR " + col] =  (pooled_data_stats.loc[0, "POOLED_VAR " + col] * (pooled_data_stats.loc[0, "POOLED_SAMPLE_SIZE " + col] -1)
                    + variance_col * (ss -1))/(pooled_data_stats.loc[0, "POOLED_SAMPLE_SIZE " + col] + ss -2)
            pooled_data_stats.loc[0, "POOLED_SAMPLE_SIZE " + col] += ss
            pooled_data_stats.loc[0, "POOLED_COUNTS " + col] += 1

    # Done with all files in files.

    # Store the pooled_data_stats in df_out.
    for col in col_labels:
        df_out.loc[0, "AVG[" + col + "]"] = pooled_data_stats.loc[0, "POOLED_AVG " + col]
        df_out.loc[0, "VAR[" + col + "]"] = pooled_data_stats.loc[0, "POOLED_VAR " + col]
        if add_counts:
            df_out.loc[0, "COUNT[" + col + "]"] = pooled_data_stats.loc[0, "POOLED_COUNTS " + col]
    
    return df_out


'''# Summary of the function gen_MEAN_INDVL_Colsfiledata(...)
# This function generates the mean for each column for each file in files, as well as the mean across all files for each column.
# It will exclude exclude_col_labels from the mean calculations. If the files have non-standardised columns 
# (different names other than exclude_col_labels), the function will return an error.
# It will create a dataframe with the following column structure:
#R_max   AVG[{var}]_SURV   ...   AVG[{var}]_ALL   ...   AVG[{var}]_R_0   ...   AVG[{var}]_R_{R_max}
# where {var} is one of the species names NOT in exclude_col_labels.
'''
def gen_MEAN_INDVL_Colsfiledata(files, pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L"]):
    
    # Read the first file in files.
    if(ext == "csv"):
        df = pan.read_csv(files[0], header=0)
    else:
        try:
            df = pan.read_table(files[0], header=0)
        except Exception as e:
            print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
            return None
    
    # Store column names in col_labels. Strip elements of any leading or trailing whitespaces.
    col_labels = df.columns; col_labels = [col.strip() for col in col_labels]
    col_labels = [col for col in col_labels if col not in exclude_col_labels];
    # Remove exclude_col_labels from col_labels.
    tot_col_labels = len(col_labels)
    # Create a df with col_labels as column names. This will store the mean data, fot each column in col_labels for each file in files.
    # Create another df which will concatenate the mean data across all files.
    df_file = pan.DataFrame(columns= ["AVG[" + x + "]" for x in col_labels])
    df_pooled = pan.DataFrame()
    # Initialise the pooled_data_stats to 0.
    # Iterate over all files in files, sorted in ascending order of R.
    for file in sorted(files, key= lambda x: int(re.sub("R_", "", re.findall(r'R_[\d]+', x)[0]))):
        # Read the file.
        Rstr = re.findall(r'R_[\d]+', file)[0]
        R = int(re.sub("R_", "", Rstr))
        #R = int(re.search(r'R_[\d]+', file).group(0).split("_")[1])
        if(ext == "csv"):
            df = pan.read_csv(file, header=0)
        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
                return None

        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        df = df[[col for col in df.columns if col not in exclude_col_labels]]

        # First check if the columns in df are the same as col_labels. If not, return an error.
        if not all([col in df.columns for col in col_labels]):
            print("Error: Columns in file " + pathtodir + "/" + file + " are not standardised.")
            return None

        # Consolidate df_file based on R value.
        if getattr(df_file, "_data", None) is not None and df_file._data.nblocks > 100:
            df_file._consolidate_inplace()
        if getattr(df_pooled, "_data", None) is not None and df_pooled._data.nblocks > 100:
            df_pooled._consolidate_inplace()
        
        # Find mean of each column in df.
        for col in col_labels:
            mean_col = df[col].mean()
            # Store the mean in df_file.
            df_file.loc[0, "AVG[" + col + "]"] = mean_col
            # Store the mean in df_pooled.
            df_pooled.loc[0, "AVG[" + col + "]_" + Rstr] = mean_col
    
    # Done with all files in files.
    # Finally insert mean across all files for each column in col_labels in df_pooled at the start.
    for i in range(len(col_labels)):
        df_pooled.insert(i, "AVG[" + col_labels[i] + "]_ALL", df_pooled.filter(like= "AVG[" + col_labels[i] + "]_").mean(axis=1))

    # Also insert mean across all files for each column in col_labels in df_pooled that are non-zero at the start.
    for i in range(len(col_labels)):
        df_pooled.insert(i, "AVG[" + col_labels[i] + "]_SURV", df_pooled.filter(like= "AVG[" + col_labels[i] + "]_R").replace(0, np.nan).mean(axis=1))
    
    # Insert max R at the start of the dataframe.
    df_pooled.insert(0, "R", len(files))
    # Replace all NaNs with 0.
    df_pooled = df_pooled.fillna(0)
    return df_pooled

''' # Summary of the function gen_FFT_PowerSpectra(...)
    The function gen_FFT_PowerSpectra(files, pathtodir="", ext="csv", bin_mask= None, exclude_col_labels= ["a_c", "x", "L"],
    Tmin= None, Tmax = None, L="infer" verbose= False) generates the FFT power spectra of the files in files.
    It will read all files in files, and will generate the FFT power spectra for each matching column in the files for each file (replicate).
    It will first resize each column in the files to a square matrix of size L x L, where L is the maximum length of the column in the files.
    L is either provided as an argument or inferred from the maximum length of the column (as the square root of the maximum length of the column).
    If bin_mask is provided, the function will use bin_mask to bin the data in each column before generating the FFT power spectra.
    Valid bin_mask values are "read", "auto" or None. If bin_mask is "read", the function will assume the presence of a binarised mask 
    for each file in files in <path_to_dir>/BIN_MASKS  and use this mask to bin the data in each column before generating the FFT power spectra.
    If bin_mask is "auto", the function will generate a binarised mask for each file in files and use this mask to bin the data in each column 
    as well as save the binarised mask in <path_to_dir>/BIN_MASKS. If bin_mask is None, the function will not bin the data in each column before doing the FFT.
    The function will then generate the FFT power spectra for each column in the files for each file (replicate).
    It will store the FFT power spectra in a dataframe with the following column structure:
    "FREQ[{var}]", "POWER[{var}]_MEAN", "POWER[{var}]_R_0", ...., "POWER[{var}]_R_{R_max}"
    where {var} is one of the species names NOT in exclude_col_labels.
    The function will also generate the mean FFT power spectra across all files for each column in the files, which is stored in POWER[{var}]_MEAN.
'''
def gen_FFT_PowerSpectra(files, pathtodir="", ext="csv", bin_mask= None, exclude_col_labels= ["a_c", "x", "L"], 
                         Tmin= None, Tmax = None, Lval="infer", binwidth=1, verbose= False):
    
    # Read the first file in files.
    # If L is not provided, infer L from the maximum length of the column in the files.
    if(ext == "csv"):
        df = pan.read_csv(files[0], header=0)
    else:
        try:
            df = pan.read_table(files[0], header=0)
        except Exception as e:
            print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
            return None
    
    # Store column names in col_labels. Strip elements of any leading or trailing whitespaces.
    col_labels = df.columns; col_labels = [col.strip() for col in col_labels]
    col_labels = [col for col in col_labels if col not in exclude_col_labels];
    # Remove exclude_col_labels from col_labels.
    
    # This dataframe will store the FFT power spectra for each column in the files for each file (replicate)
    # This dataframe will store the mean FFT power spectra across all files for each column in the files.
    df_fft_power = pan.DataFrame()

    
    Tval = re.search(r'T_[\d]+', files[0]).group(0).split("_")[1] if re.search(r'T_[\d]+', files[0]) else re.search(r'T_[\d]*[.][\d]+', files[0]).group(0).split("_")[1]
    aval = re.search(r'a_[\d]*[.][\d]+' , files[0]).group(0).split("_")[1] if re.search(r'a_[\d]*[.][\d]+' , files[0]) else re.search(r'a_[\d]+' , files[0]).group(0).split("_")[1]
    # Check if Tmin and Tmax are provided. If they are, only calculate KDE and potential well data for Tval in the range [Tmin, Tmax].
    if Tmin is not None:
        if float(Tval) < Tmin:
            print(f"Skipping T = {Tval} and a = {aval} as it is not in the range [{Tmin}, {Tmax}].")
            return None
    if Tmax is not None:
        if float(Tval) > Tmax:
            print(f"Skipping T = {Tval} and a = {aval} as it is not in the range [{Tmin}, {Tmax}].")
            return None
    print(f"Generating FFT Power Spectra  for T = {Tval} and a = {aval}.")

    # Iterate over all files in files in ascending order of R.
    for file in sorted(files, key= lambda x: int(re.sub("R_", "", re.findall(r'R_[\d]+', x)[0]))):
        Rstr = re.findall(r'R_[\d]+', file)[0]
        #print(f"subdir: {pathtodir}, file: {file}, R: {Rstr}")
        #filename = file.split(pathtodir + "\\")[1]
        filename = os.path.basename(file)
        # Read the file.
        if(ext == "csv"):
            df = pan.read_csv(file, header=0)

        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
                return None
        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        df = df[[col for col in df.columns if col not in exclude_col_labels]]
        
        # First check if the columns in df are the same as col_labels. If not, return an error.
        if not all([col in df.columns for col in col_labels]):
            print("Error: Columns in file " + pathtodir + "/" + file + " are not standardised.")
            return None
        
        # Resize each column in the files to a square matrix of size L x L.
        if Lval == "infer":
            #L = int(np.sqrt(max([len(df[col]) for col in col_labels]))) SEVERE CUDA PERFORMANCE HIT
            L = round(asnumpy(max([len(df[col]) for col in col_labels]))**0.5)
            # Recall asnumpy() returns a numpy array from a cupy array.
            # Or if already int, string or float, it will return the same value.
        else:
            L = Lval

        if bin_mask == "read":
            try:
                # Bin mask is assumed to be in <path_to_dir>/BIN_MASKS with the same name as the file.
                # It has column names as the same as the columns in the file, with each column 
                # bearing a single threshold value (float b/w 0 and 1) for binarisation.
                BINmask = pan.read_csv(pathtodir + "/BIN_MASKS/" + filename.split("." + ext)[0] + ".txt", header="infer", sep="\t")
                BINmask.columns = [col.strip() for col in BINmask.columns]
                # Remove exclude_col_labels from bin_mask.
                BINmask = BINmask[[col for col in bin_mask.columns if col not in exclude_col_labels]]
                # First check if the columns in bin_mask are the same as col_labels. If not, return an error.
                if not all([col in bin_mask.columns for col in col_labels]):
                    print(f"ERROR: Columns in binarised mask for file: {file} are not standardised. Proceeding without binarisation.")
                    bin_mask = None
            except Exception as e:
                print("Error: Could not load binarised mask for file " + pathtodir + "/" + file + " with error message: \n" + str(e))
                print("Proceeding without binarisation....")
                bin_mask = None
        elif bin_mask == "auto" or bin_mask == "GMM":
            # We will create a binarised mask for each file in files and save it in <path_to_dir>/BIN_MASKS.
                # The binarised mask will have column names as the same as the columns in the file, with each column
                # bearing a single threshold value (float b/w 0 and 1) for binarisation.
                # This threshold value will be determined using GMM clustering with 2 clusters.
                thresholds =pan.DataFrame(columns= col_labels)
                for col in col_labels:
                    data = df[col].dropna().to_numpy(dtype=np.float64).reshape(-1, 1)
                    # Normalise the data in each column by max value in the column.
                    data /= data.max() if data.max() != 0.0 else 1.0
                    gmm = GaussianMixture(n_components=2, covariance_type='full', random_state=gpu._cpu_np.random.randint(0, 1000))
                    gmm.fit(data)
                    binned_data = gmm.predict(data)
                    # Set threshold as the mean of the means of the two clusters.
                    thresholds.loc[0, col] = float(np.mean(gmm.means_.flatten()))
                    # Next sort the centers in ascending order and reassign cluster ids accordingly.
                    sorted_indices = np.argsort(gmm.means_.flatten())
                    binned_data = sorted_indices[binned_data]
                    # Assign the binned data to the column in df.
                    df[col] = pan.Series(binned_data)

                # Print all thresholds.
                print(f"Thresholds for file {filename} are: \n" + str(thresholds))
                # Save the binarised mask in <path_to_dir>/BIN_MASKS/.
                savetxt_to_dir(thresholds, pathtodir + "/BIN_MASKS", filename.split("." + ext)[0] + ".txt", 
                               header = True, index = False, sep = "\t", overwrite = True, verbose = verbose)
                
                
                # Assign the binarised mask to bin_mask.
                BINmask = thresholds

        kbins = np.arange(1 - binwidth/2, L//2 + 1, binwidth) # Bin the data in each column before generating the FFT power spectra.
        kvals = 0.5*(kbins[1:] + kbins[:-1]) #k values at bin centers (1 , 2, 3, ... L//2)
        # Store the k values in the dataframe at the start.
        df_fft_power.insert(len(df_fft_power.columns), "k_" + Rstr, pan.Series(kvals))
        for col in col_labels:
            # Resize each column in the files to a square matrix of size L x L.
            # Bin the data in each column before generating the FFT power spectra.
            
            
            img_arr = asarray(df[col].values).astype(np.float64)  # Convert to float64 for FFT calculations.
            # Ensure img_arr is dtype float64.

            if gpu.GPU_AVAILABLE:
                print(f"Is img_arr at this point a CUPY array? {isinstance(img_arr, gpu._cupy.ndarray)}")
                print(f"Is img_arr at this point a NUMPY array? {isinstance(img_arr, gpu._cpu_np.ndarray)}")
            img_arr = np.resize(img_arr, (L, L))
            #mean_arr = img_arr.mean(); var_arr = img_arr.var()
            #print(f"Mean of {col} = {mean_arr}, Variance of {col} = {var_arr}, BEFORE NORMALISATION.")
            # Normalise the data in each column by max value in the column.
            img_arr /= img_arr.max() if img_arr.max() != 0.0 else 1.0
            mean_arr = img_arr.mean(); var_arr = img_arr.var()

            if var_arr == 0 and mean_arr == 0:
                print(f"Warning: Homogenous 0 column {col} at a = {aval}, T = {Tval}, R = {Rstr}...")
            if bin_mask == "read":
                try:
                    # Binarise the data in each column in the file using the binarised mask.
                    #img_arr = np.where(BINmask[col].to_numpy(dtype=np.float64) > img_arr, 0, img_arr)
                    # Print and check if img_arr is a numpy array or a cupy array.
                    if gpu.GPU_AVAILABLE:
                        print(f"Is img_arr CUPY array? {isinstance(img_arr, gpu._cupy.ndarray)}")
                        print(f"Is img_arr NUMPY array? {isinstance(img_arr, gpu._cpu_np.ndarray)}")
                        print(f"Is asarray(BINmask[col].values) CUPY array? {isinstance(asarray(BINmask[col].values), gpu._cupy.ndarray)}")
                        print(f"Is asarray(BINmask[col].values) NUMPY array? {isinstance(asarray(BINmask[col].values), gpu._cpu_np.ndarray)}")
                        print(f"Is BINmask[col].values a Pandas Series? {isinstance(BINmask[col].values, gpu._cpu_pandas.Series)}")
                        print(f"Is BINmask[col].values a Numpy array? {isinstance(BINmask[col].values, gpu._cpu_np.ndarray)}")
                    img_arr = np.where(asarray(BINmask[col].values) > img_arr, 0, img_arr)
                    # If bin
                except Exception as e:
                    print("Error: Could not binarise data in column " + col + " for file " + pathtodir + "/" + file + " with error message: \n" + str(e))
                    print("Proceeding without binarisation....")
            
            if verbose:
                print(f"Mean of {col} = {mean_arr}, Variance of {col} = {var_arr}, AFTER NORMALISATION.")
                
            #print(f"Mean of {col} = {mean_arr}, Variance of {col} = {var_arr}, AFTER NORMALISATION.")
            # Reshape the column to a square matrix of size L x L , and subtract the mean from each element in the column.
            

            # Get mean, var and SD of each column in df.
             #sd_col = np.sqrt(var_col)
            # First subtract the mean from each column in df.
            #df[col] = df[col].apply(lambda x: x - mean_col)
            # Find 2D FFT of each column in df
            img_fft = fft.fft2(img_arr)

            # Calculate the power spectra of each column in df by performing the following steps:
            # 1. Compute the 2D FFT of the column data.
            # 2. Flatten the 2D FFT result to a 1D array.
            # 3. Compute the Fourier amplitudes by taking the absolute value squared of the FFT result.
            # 4. Bin the Fourier amplitudes based on the radial distance from the origin in frequency space.
            # 5. Store the binned Fourier amplitudes as the power spectra for the column.
            kfreqX =fft.fftfreq(L)*L; #kfreqY = fft.fftfreq(L)*L
            kfreqX, kfreqY = np.meshgrid(kfreqX, kfreqX)
            # Calculate the radial distance from the origin in frequency space
            knrm = np.sqrt(kfreqX**2 + kfreqY**2)
            knrm = knrm.flatten(); fourier_amplitudes = np.abs(img_fft**2)
            fourier_amplitudes = fourier_amplitudes.flatten()
            
            #Ampbins, _, _ = stats.binned_statistic(knrm, fourier_amplitudes, bins=kbins, statistic='mean')
            bin_sums, _ = np.histogram(knrm, bins=kbins, weights=fourier_amplitudes)
            bin_counts, _ = np.histogram(knrm, bins=kbins)
            # Calculate the average power in each bin, avoiding division by zero.
            Ampbins = np.divide(bin_sums, bin_counts, out=np.zeros_like(bin_sums), where=bin_counts!=0)
            #WARNING: Unlike stats.binned_statistic, np.histogram does not return NaN for empty bins, but 0.
            # So we need to replace 0 with NaN in Ampbins.
            # Store the FFT power spectra for each column in the files for each file at the END of the dataframe.
            df_fft_power.insert(len(df_fft_power.columns), "POWER[" + col + "]_" + Rstr, pan.Series(Ampbins))
            # Replace all 0s in Ampbins with NaN.
            df_fft_power["POWER[" + col + "]_" + Rstr] = df_fft_power["POWER[" + col + "]_" + Rstr].replace(0, np.nan)
        # Done with all columns in col_labels for file.
    # Done with all files in files.
    
    for i in range(len(col_labels)):
        # Insert mean FFT power spectra and k-vals across all files for each column in the files in df_fft_power at the start.
        df_fft_power.insert(i, "POWER[" + col_labels[i] + "]_MEAN_ALL", df_fft_power.filter(like= "POWER[" + col_labels[i] + "]_").mean(axis=1))
    for i in range(len(col_labels)):
        # Insert mean FFT power spectra across all files for each NON-ZERO column in the files in df_fft_power at the start.
        df_fft_power.insert(i, "POWER[" + col_labels[i] + "]_MEAN_SURV", df_fft_power.filter(like= "POWER[" + col_labels[i] + "]_").replace(0, np.nan).mean(axis=1))
    for i in range(len(col_labels)):  
        df_fft_power.insert(i, "FREQ[" + col_labels[i] + "]", df_fft_power.filter(like= "k_").mean(axis=1))
            
    return df_fft_power

''' # Summary of the function quadratic_spline_roots(...)
    The function quadratic_spline_roots(spl) finds the roots of a quadratic spline function spl.
    It does so by finding the roots of the quadratic spline function in each interval [a, b] of the spline.
    Stolen from: https://stackoverflow.com/questions/50371298/find-maximum-minimum-of-a-1d-interpolated-function
'''
def quadratic_spline_roots(spl):
    roots = []
    knots = spl.get_knots()
    for a, b in zip(knots[:-1], knots[1:]):
        u, v, w = spl(a), spl((a+b)/2), spl(b)
        t = np.roots([u+w-2*v, w-u, 2*v])
        t = t[np.isreal(t) & (np.abs(t) <= 1)]
        roots.extend(t*(b-a)/2 + (b+a)/2)
    return np.array(roots)

""" # Summary of the function gen_potential_well_data(...)
The following function gen_potential_well_data(files, pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L"], evaluate_local_minima =True)
will generate the potential well data from the files in files.
It will read all files in files, and will concatenate matching columns in the files to generate the potential well data.
After concatenating the columns across all files, for each column, it will use stats.gaussian_kde to generate a kernel density estimate of the data.
It will then find the potential well for each column as potential well = -np.log10(kde(x)), where x is the data in the column.
If evaluate_local_minima is True, it will evaluate the local minima of the potential well data for each column.
It will do so by using the CubicSpline function in scipy.interpolate to find all the local minima of the potential well data.
In case of an entirely flat potential well, the function will return a warning and return local minima for that column as None.
It will then store the local minima in a different dataframe, along with the corresponding x-values.
The function will then return the potential well dataframe and (optionally) the local minima dataframes.
"""
def gen_potential_well_data(files, pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L"], Tmin= None, Tmax = None, evaluate_local_minima= True, bins=300):
    # Read the first file in files.
    if(ext == "csv"):
        df = pan.read_csv(files[0], header=0)
    else:
        try:
            df = pan.read_table(files[0], header=0)
        except Exception as e:
            print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
            return None, None
    
    # Store column names in col_labels. Strip elements of any leading or trailing whitespaces.
    col_labels = df.columns; col_labels = [col.strip() for col in col_labels];
    col_labels = [col for col in col_labels if col not in exclude_col_labels];
    # Remove exclude_col_labels from col_labels.

    kde_cols = [y + x + "]" for x in col_labels for y in ["BINS[", "KDE[", "POTENTIAL_WELL["]]
    # Create a df with col_labels as column names. This will store the KDE and potential well data.
    df_kde = pan.DataFrame(columns= kde_cols)
    if(evaluate_local_minima):
        min_cols = ["LOCAL_MINIMA[" + x + "]" for x in col_labels]
        # Create a df with col_labels as column names. This will store the local minima data.
        df_local_minima = pan.DataFrame(columns= min_cols)
    
    concat_df = pan.DataFrame()
    # Iterate over all files in files. For each file, concatenate the columns in col_labels across all files.
    # Store the concatenated data in concat_df.

    # Get T and a values from the first file.
    Tval = re.search(r'T_[\d]+', files[0]).group(0).split("_")[1] if re.search(r'T_[\d]+', files[0]) else re.search(r'T_[\d]*[.][\d]+', files[0]).group(0).split("_")[1]
    aval = re.search(r'a_[\d]*[.][\d]+' , files[0]).group(0).split("_")[1] if re.search(r'a_[\d]*[.][\d]+' , files[0]) else re.search(r'a_[\d]+' , files[0]).group(0).split("_")[1]
    # Check if Tmin and Tmax are provided. If they are, only calculate KDE and potential well data for Tval in the range [Tmin, Tmax].
    if Tmin is not None:
        if float(Tval) < Tmin:
            print(f"Skipping T = {Tval} and a = {aval} as it is not in the range [{Tmin}, {Tmax}].")
            return None, None
    if Tmax is not None:
        if float(Tval) > Tmax:
            print(f"Skipping T = {Tval} and a = {aval} as it is not in the range [{Tmin}, {Tmax}].")
            return None, None
    print(f"Generating potential well data for T = {Tval} and a = {aval}.")
    for file in files:

        #First read the file.
        if(ext == "csv"):
            df = pan.read_csv(file, header=0)
        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
                return None
        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        df = df[[col for col in df.columns if col not in exclude_col_labels]]

        # First check if the columns in df are the same as col_labels. If not, return an error.
        if not all([col in df.columns for col in col_labels]):
            print("Error: Columns in file " + pathtodir + "/" + file + " are not standardised.")
            return None, None
        
        concat_df = pan.concat([concat_df, df[col_labels]], axis=0)

    #print(concat_df.head())
    #print(concat_df.tail())
    print(concat_df.shape)
    # Done with all files in files.

    # Now generate the KDE and potential well data for each column in col_labels.
    for col in col_labels:
        is_zero = False; is_flat = False; # Flags to check if the column is entirely zero or flat. 
        #print("Evaluating column: " + col)
        #First check if the column is entirely zero. If it is, skip the column.
        if all(concat_df[col] == 0):
            print("Warning: Entirely zero column across all replicates in " + col + f" for T = {Tval} and a = {aval} ... Adjusting KDE ...")
            # Define eval_range between 0 and 1, with 300 bins. Store the KDE as 1e-10 for all values in eval_range, except 0.
            eval_range = np.linspace(0, 1, bins); kernel = np.ones(bins) * 1e-10; kernel[0] = 1; is_zero = True;
        
        # Check if the column is entirely flat, i.e. std < 1e-6. If it is, skip the column.
        elif concat_df[col].std() < 1e-6:
            print("Warning: Entirely flat column across all replicates in " + col + f"  for T = {Tval} and a = {aval} ... Adjusting KDE ...")
            # Define eval_range b/w min - 1 and max + 1, with 300 bins. Store the KDE as 1e-10 for all values in eval_range except the mean of the column.
            eval_range = np.linspace(concat_df[col].min() - 1, concat_df[col].max() + 1, bins); is_flat = True;
            kernel = np.ones(bins) * 1e-10;
            # Change the kernel value at the index bin closest to the mean of the column to 1.
            kernel[np.argmin(np.abs(eval_range - concat_df[col].mean()))] = 1;
        else:
            # Evaluate this over a certain range of x-values for each column.
            eval_range = np.linspace(concat_df[col].min(), concat_df[col].max(), bins)

            # Find the kernel density estimate of the data in col.
            try:
                kde = stats.gaussian_kde(concat_df[col]); kernel = kde.evaluate(eval_range)
            except stats.NearConstantInputWarning  or stats.ConstantInputWarning as e:
                print(f"Warning: {e} for column " + col + f" for T = {Tval} and a = {aval}. Skipping column.")
                continue

            # Report indices of kernel where kernel values are zero.
            if np.any(kernel == 0):
                #print(f"Warning: Kernel is zero for T = {Tval} and a = {aval} at indices: " + str(np.where(kernel == 0)))
                # Replace zero values in kernel with 1e-10.
                kernel[kernel == 0] = 1e-10
        
        # Find the potential well for the data in col, accounting for zero values in kernel.
        potential_well = -np.log10(kernel)
        # Store the KDE and potential well data in df_kde.
        df_kde["BINS[" + col + "]"] = pan.Series(eval_range)
        df_kde["KDE[" + col + "]"] = pan.Series(kernel)
        df_kde["POTENTIAL_WELL[" + col + "]"] = pan.Series(potential_well)

        #print(df_kde.head())
        #print(df_kde.tail())
        #print(df_kde.shape)
        #print(df_kde.describe())

        if(evaluate_local_minima):
            # Find the local minima of the potential well data in col.
            # Use CubicSpline to find the local minima.
            local_minima = []; min_pot =[];
            # local_minima will store the x-values corresponding to the local minima of the potential well.
            # min_pot will store the potential well values at the local minima (y-values).
            # First check if the potential well is entirely flat.
            if is_zero or is_flat:
                # Append the eval_range corresponding to potential well = 0 to local_minima.
                local_minima.append(eval_range[np.argmin(np.abs(potential_well))])
                min_pot.append(0)
                df_local_minima["LOCAL_MINIMA[" + col + "]"] = pan.Series(local_minima);
                df_local_minima["MIN_POT[" + col + "]"] = pan.Series(min_pot);
                continue

            cs = interpolate.InterpolatedUnivariateSpline(eval_range, potential_well, k=3)
            # Find the local minima.
            cs_1st_derivative = cs.derivative()
            cs_2nd_derivative = cs.derivative(n=2)
            extrema = quadratic_spline_roots(cs_1st_derivative)
            for x in extrema:
                if cs_2nd_derivative(x) > 0:
                    local_minima.append(x)
                    min_pot.append(cs(x))
            # If the potential well is entirely flat, return a warning and set local minima to None.
            if len(local_minima) == 0:
                print("Warning: Entirely flat potential well for column " + col + ". Local minima set to NAN.")
                local_minima.append(np.nan)
                min_pot.append(np.nan)
            # Store the list of local minima in df_local_minima.
            df_local_minima["LOCAL_MINIMA[" + col + "]"] = pan.Series(local_minima)
            df_local_minima["MIN_POT[" + col + "]"] = pan.Series(min_pot)

    # Done with all columns in col_labels.     
    return (df_kde, df_local_minima) if evaluate_local_minima else (df_kde, None)

""" Summary of the function gen_clustered_data(...)
The function gen_clustered_data() will read all files in files, and will generate clustered datsets for each file (replicate).
It will cluster the data in each file using the specified clustering algorithm (cluster_type) into n_clusters clusters.
It will save these clustered datasets as dataframes, and if evaluate_clusterfrequencies is True, it will also evaluate the cluster frequencies.
This is done by counting the size of each cluster (cluster sizes determined using scipy.ndimage.label) in each file.
The frequency of cluster sizes across all files (replicates) is then stored in a dataframe and returned.
If periodic is True, the cluster counting is done using periodic boundary conditions.
IMPORTANT: 
Cluster_type can be: "Kmeans", "GMM" (Gaussian Mixture Model)
    "Zero" (all data points above 0 are assigned to cluster 1, all data points < = 0 are assigned to cluster 0)
    "ZeroKmeans" or "ZeroGMM" (where some specified columns are zero-clustered as above, the rest are clustered using Kmeans or GMM respectively).

If "ZeroKmeans" or "ZeroGMM" is specified, the function will then look at an optional argument exclude_col_from_algorithm, which is a 
list of column names that should not be clustered using Kmeans or GMM (i.e. these columns will be zero-clustered).
If exclude_col_from_algorithm is not provided (empty list), "ZeroKmeans" or "ZeroGMM" will be identical to "Kmeans" or "GMM" respectively.
"""
def gen_clustered_data(files, pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L"], Tmin= None, Tmax = None, save_frames= True,
                        n_clusters=2, cluster_type="KMeans", evaluate_clusterfrequencies = True, periodic=True, verbose= False):
    
    # Get T and a values from the first file.
    Tval = re.search(r'T_[\d]+', files[0]).group(0).split("_")[1] if re.search(r'T_[\d]+', files[0]) else re.search(r'T_[\d]*[.][\d]+', files[0]).group(0).split("_")[1]
    aval = re.search(r'a_[\d]*[.][\d]+' , files[0]).group(0).split("_")[1] if re.search(r'a_[\d]*[.][\d]+' , files[0]) else re.search(r'a_[\d]+' , files[0]).group(0).split("_")[1]
    # Check if Tmin and Tmax are provided. If they are, only calculate KDE and potential well data for Tval in the range [Tmin, Tmax].
    if Tmin is not None:
        if float(Tval) < Tmin:
            print(f"Skipping T = {Tval} and a = {aval} as it is not in the range [{Tmin}, {Tmax}].")
            return None
    if Tmax is not None:
        if float(Tval) > Tmax:
            print(f"Skipping T = {Tval} and a = {aval} as it is not in the range [{Tmin}, {Tmax}].")
            return None
    print(f"Generating clustered data for T = {Tval} and a = {aval}.")
    
    # Read the first file in files.
    if(ext == "csv"):
        df = pan.read_csv(files[0], header=0)
    else:
        try:
            df = pan.read_table(files[0], header=0)
        except Exception as e:
            print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
            return None
    # Store column names in col_labels. Strip elements of any leading or trailing whitespaces.
    col_labels = df.columns; col_labels = [col.strip() for col in col_labels]
    col_labels = [col for col in col_labels if col not in exclude_col_labels];
    count_cols =  [y + x + "]" for x in col_labels for y in ["CLUSTER_SIZE[", "CLUSTER_COUNTS[", "CLUSTER_FREQ["]]
    
    """ Additional logic to handle cluster_type OF "ZeroKmeans" or "ZeroGMM":
    if( (cluster_type == "ZeroKmeans" or cluster_type == "ZeroGMM")):
        if exclude_col_labels_from_algo == [] or exclude_col_labels_from_algo is None:
            print(f"Warning: No columns specified to be excluded from clustering algorithm {cluster_type}. Treating as {cluster_type.replace('Zero', '')}.")
            cluster_type = cluster_type.replace("Zero", "") # If no columns are specified to be excluded from clustering algorithm, treat as Kmeans or GMM.
        else:
            # If cluster_type is "ZeroKmeans" or "ZeroGMM", then we will zero-cluster the columns in exclude_col_labels_from_algo.
            # The rest of the columns will be clustered using Kmeans or GMM respectively.
            col_labels = [col for col in col_labels if col not in exclude_col_labels_from_algo]
            print(f"Clustering algorithm {cluster_type} will be ONLY be applied to columns: {col_labels}...")
    all_col_labels = list(set(col_labels + exclude_col_labels_from_algo)) if exclude_col_labels_from_algo else col_labels
    """
    if evaluate_clusterfrequencies:
        # Create a df with col_labels as column names. This will store the cluster frequencies.
        df_cluster_freq = pan.DataFrame(columns= count_cols)

    # Default dict to store cluster sizes and counts for each column in col_labels.
    cluster_freqs = defaultdict(Counter)

    Path(pathtodir + f"/CLUST/{cluster_type}_{n_clusters}").mkdir(parents=True, exist_ok=True)

    for file in files:

        # First read the file.
        if(ext == "csv"):
            df = pan.read_csv(file, header=0)
        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
                return None
        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        df = df[[col for col in df.columns if col not in exclude_col_labels]]

        # First check if the columns in df are the same as col_labels. If not, return an error.
        if not all([col in df.columns for col in col_labels]):
            print("Error: Columns in file " + pathtodir + "/" + file + " are not standardised.")
            return None
        
        # Get L value from the first file.
        Lval = re.search(r'L_[\d]+', files[0]).group(0).split("_")[1] if re.search(r'L_[\d]+', files[0]) else re.search(r'L_[\d]*[.][\d]+', files[0]).group(0).split("_")[1]
        Lval = int(Lval) if Lval.isdigit() else round(float(Lval))

        L = round(asnumpy(df[col_labels[0]].size**0.5))
        # Recall asnumpy() returns a numpy array from a cupy array.
        assert L == Lval, f"NOTE : L value in file {file} does not match the expected L value {Lval}. Found {L}."

        # Normalise the data in each column by max value in the column.
        df[col_labels] = df[col_labels].apply(lambda x: x / x.max() if x.max() != 0.0 else 0, axis=0)

        if save_frames:
            bin_df = pan.DataFrame(columns=col_labels) # This dataframe will store the binned data for each column in the files.
        
        for col in col_labels:
            is_zero = False; is_flat = False; # Flags to check if the column is entirely zero or flat. 
            
            #print(f"Evaluating column: {col} at T = {Tval}, a = {aval} and R = {re.search(r'R_[\d]+', file).group(0).split('_')[1]}...")
            # First check if the column is entirely zero. If it is, skip the column.
            if all(df[col] == 0):
                print(f"Warning: Entirely zero column in {col} for T = {Tval} and a = {aval} ... Skipping column.")
                is_zero = True
                continue
            # Check if the column is entirely flat, i.e. std < 1e-6. If it is, skip the column.
            if df[col].std() < 1e-6:
                print(f"Warning: Entirely flat column  in {col} for T = {Tval} and a = {aval} ... Skipping column.")
                is_flat = True
                continue
            
            #Perform clustering on the data in col using the specified clustering algorithm.
            if cluster_type == "KMeans":
                test_df = pan.DataFrame({col: df[col]})
                kmeans_model = KMeans(n_clusters=n_clusters, init= 'k-means++' , max_iter= 500, 
                random_state=42).fit(test_df)
                # Get the cluster labels for each data point in col.
                cluster_labels = kmeans_model.labels_
                df["CLUSTER[" + col + "]"] =  kmeans_model.predict(test_df)
                centers = kmeans_model.cluster_centers_.flatten()
                ''' # For ZeroKmeans or ZeroGMM, we will zero-cluster the columns in exclude_col_labels_from_algo.
                elif cluster_type == "Kmeans" or cluster_type == "ZeroKmeans":
                    if col in col_labels:
                        test_df = pan.DataFrame({col: df[col]})
                        kmeans_model = KMeans(n_clusters=n_clusters, init= 'k-means++' , max_iter= 500, 
                        random_state=42).fit(test_df)
                        # Get the cluster labels for each data point in col.
                        cluster_labels = kmeans_model.labels_
                        df["CLUSTER[" + col + "]"] =  kmeans_model.predict(test_df)
                        centers = kmeans_model.cluster_centers_.flatten()
                    else:
                        # If col is in exclude_col_labels_from_algo, then we will zero-cluster the column.
                        df["CLUSTER[" + col + "]"] = df[col].apply(lambda x: 1 if x > 0 else 0)
                        centers = np.array([0, 1])
                #'''
            elif cluster_type == "GMM":
                test_df = pan.DataFrame({col: df[col]})
                gmm_model = GaussianMixture(n_components=n_clusters, covariance_type='full', random_state=42).fit(test_df)
                #Get the cluster labels for each data point in col.
                cluster_labels = gmm_model.fit_predict(test_df)
                df["CLUSTER[" + col + "]"] = cluster_labels
                centers = gmm_model.means_.flatten()
            elif cluster_type == "Zero":
                # If cluster_type is "Zero", we will assign all data points above 0 to cluster 1, and all data points <= 0 to cluster 0.
                df["CLUSTER[" + col + "]"] = df[col].apply(lambda x: 1 if x > 0 else 0)
                centers = np.array([0, 1])
            else:
                print(f"Error: Invalid cluster_type {cluster_type}. Use 'Kmeans' or 'GMM'.")
                return None
            
            ''' OLD VERSION
            palette = cluster_labels.astype(float())
            ordered_palette = np.sort(np.unique(palette))

            # Ensure that the cluster labels are in ascending order.
            sorting_order = sorted(range(len(palette)), key=lambda k: palette[k][0], reverse=True)
            if sorting_order != list(range(len(palette))):
                print(f"Warning: Cluster labels are not in ascending order for column {col} in file {file}. Reordering cluster labels.")
                for new_order, order in enumerate(sorting_order):
                    df['CLUSTER[' + col + ']'] = df['CLUSTER[' + col + ']'].replace(order, new_order + n_clusters)
                for new_order, order in enumerate(ordered_palette):
                    df['CLUSTER[' + col + ']'] = df['CLUSTER[' + col + ']'].replace(new_order + n_clusters, new_order)
            #'''

             # Reorder cluster labels so that label 0 corresponds to the lowest cluster center, etc.
            
            # Get the order that sorts the centers
            order = np.argsort(centers)
            # Map old labels to new labels
            label_map = {old: new for new, old in enumerate(asnumpy(order))}
            df['CLUSTER[' + col + ']'] = df['CLUSTER[' + col + ']'].map(label_map)

            if save_frames:
                bin_df[col] = df['CLUSTER[' + col + ']']
            

            if evaluate_clusterfrequencies:
                # First binarise df['CLUSTER[' + col + ']'] to get the cluster sizes.
                # Here all values in df['CLUSTER[' + col + ']'] are either 0 or 1, where 1 indicates the presence of a cluster.
                # All values in df['CLUSTER[' + col + ']'] > 0 are set to 1.
                #bin_df[col] = df['CLUSTER[' + col + ']'].apply(lambda x: 1 if x > 0 else 0).values.reshape((L, L))
                binary_mask = asarray(df[f"CLUSTER[{col}]"].values > 0).astype(int).reshape((L, L))
                # Count the size of each cluster in col using scipy.ndimage.label.

                non_periodic_labels, num_ind_clusters = ndimage.label(binary_mask, structure=np.ones((3, 3)))
                if periodic:
                    for x in range(L):
                        if non_periodic_labels[x, 0] > 0 and non_periodic_labels[x, -1] > 0:
                            non_periodic_labels[non_periodic_labels == non_periodic_labels[x, -1]] = non_periodic_labels[x, 0]
                            #Dealing with top and bottom borders.
                    for x in range(0, L):
                        if non_periodic_labels[0, x] > 0 and non_periodic_labels[-1, x] > 0:
                            non_periodic_labels[non_periodic_labels == non_periodic_labels[-1, x]] = non_periodic_labels[0, x]
                            #Dealing with L and R borders.

                    # Get sorted label IDs
                    label_ids = np.unique(non_periodic_labels) #Return labels in ascending order.
                    for i in  range(1, len(label_ids)):
                        if(label_ids[i] == label_ids[i-1] + 1):
                            continue
                        #Update LABEL ID to lowest extant value possible
                        non_periodic_labels[non_periodic_labels == label_ids[i]] = label_ids[i-1] + 1
                        label_ids[i] = label_ids[i-1] + 1

                Clus_ID, Clus_Size = np.unique(non_periodic_labels, return_counts=True)

                for sz, count in Counter(asnumpy(Clus_Size[1:])).items():
                    cluster_freqs[col][sz] += count


                """ # Using df approach with np.unique
                label_ids = Clus_ID
                CSize, NS = np.unique(Clus_Size[1:], return_counts=True)
                # If df_cluster_freq["CLUSTER_SIZE[" + col + "]"] is empty, initialise it with the cluster sizes.
                if df_cluster_freq[f"CLUSTER_COUNTS[{col}]"].empty:
                    df_cluster_freq[f"CLUSTER_COUNTS[{col}]"] = NS
                if df_cluster_freq[f"CLUSTER_SIZE[{col}]"].empty:
                    df_cluster_freq[f"CLUSTER_SIZE[{col}]"] = CSize
                else:
                    for i in range(len(CSize)):
                        # If CSize[i] is in df_cluster_freq["CLUSTER_SIZE[" + col + "]"], update the cluster size and cluster counts.
                        if CSize[i] in df_cluster_freq[f"CLUSTER_SIZE[{col}]"].values:
                            df_cluster_freq[f"CLUSTER_SIZE[{col}]"].iloc[df_cluster_freq[f"CLUSTER_SIZE[{col}]"].values == CSize[i]] = CSize[i]
                            df_cluster_freq[f"CLUSTER_COUNTS[{col}]"].iloc[df_cluster_freq[f"CLUSTER_SIZE[{col}]"].values == CSize[i]] += NS[i]
                        else:
                            # If CSize[i] is not in df_cluster_freq["CLUSTER_SIZE[" + col + "]"], append it to the end.
                            df_cluster_freq[f"CLUSTER_SIZE[{col}]"] = pan.concat([df_cluster_freq[f"CLUSTER_SIZE[{col}]"], pan.Series([CSize[i]])], ignore_index=True)
                            df_cluster_freq[f"CLUSTER_COUNTS[{col}]"] = pan.concat([df_cluster_freq[f"CLUSTER_COUNTS[{col}]"], pan.Series([NS[i]])], ignore_index=True)
                              
                # Sort the above two columns in df_cluster_freq by cluster size.
                df_cluster_freq[[f"CLUSTER_SIZE[{col}]", f"CLUSTER_COUNTS[{col}]"]] = df_cluster_freq[[f"CLUSTER_SIZE[{col}]", 
                        f"CLUSTER_COUNTS[{col}]"]].sort_values(by=f"CLUSTER_SIZE[{col}]").reset_index(drop=True)

                # Reset the index of the DataFrame.
                df_cluster_freq.reset_index(drop=True, inplace=True)
                 #""" 
        # End of all columns in col_labels for file.
        # Save the clustered data given by df['CLUSTER[' + col + ']'] in the file.
        # Use the file name to save the clustered data (adding CLUSTER_ prefix to the file name).
        if save_frames:
            bin_df.to_csv(pathtodir + f"/CLUST/{cluster_type}_{n_clusters}/CLUSTERED_{os.path.basename(file)}", index=False)
    # End of all files in files.
    if evaluate_clusterfrequencies:
        # Reset the index of the DataFrame.
        #df_cluster_freq.reset_index(drop=True, inplace=True)
        all_size_fits_all = sorted({sz for col in cluster_freqs for sz in cluster_freqs[col]})
        for col in col_labels:
            counts = [cluster_freqs[col].get(sz, 0) for sz in all_size_fits_all]
            df_cluster_freq[f"CLUSTER_SIZE[{col}]"] = all_size_fits_all
            df_cluster_freq[f"CLUSTER_COUNTS[{col}]"] = counts
            df_cluster_freq[f"CLUSTER_FREQ[{col}]"] = pan.Series(np.array(counts) / sum(counts)) if sum(counts) > 0 else 0
        return df_cluster_freq
    return None



def compute_mutual_information(X, Y, bins="scotts", verbose= False):
    # x, y should be 1D arrays, if not, flatten them
    X = X.flatten() if X.ndim > 1 else X; Y = Y.flatten() if Y.ndim > 1 else Y
    # If bins is numeric, use it as the number of bins. If bins is a string, set bins using Scott's method.
    #print(f"Computing MI...")
    # Print the shapes of X and Y
    #print(f"X shape: {X.shape}, Y shape: {Y.shape}")
    # Also print the means and stds of X and Y
    #print(f"X mean: {np.mean(X)}, X std: {np.std(X)}"); print(f"Y mean: {np.mean(Y)}, Y std: {np.std(Y)}")
    #print(f"X 25th percentile: {np.percentile(X, 25)}, X 75th percentile: {np.percentile(X, 75)}")
    #print(f"Y 25th percentile: {np.percentile(Y, 25)}, Y 75th percentile: {np.percentile(Y, 75)}")
    #print(f"X min: {X.min()}, X max: {X.max()}") ; print(f"Y min: {Y.min()}, Y max: {Y.max()}")

    # If bins is numeric, use it as the number of bins. If bins is a string, set bins using Scott's method.
    try:
        if isinstance(bins, str):
            
            if bins == 'freedman':
                # Use Freedman-Diaconis rule to determine the number of bins.
                bw_X = freedman_diaconis_bin_width(X); bw_Y = freedman_diaconis_bin_width(Y);
            elif bins == 'scotts':
                # Use Scott's method to determine the number of bins.
                bw_X = scotts_bin_width(X); bw_Y = scotts_bin_width(Y);
            else:
                raise TypeError(f"Invalid binning method {bins}. Use 'freedman' or 'scotts' or integer value.")

            if np.isnan(bw_X) or np.isnan(bw_Y) or bw_X <= 0 or bw_Y <= 0:
                raise ValueError("Invalid bin width calculated. Check if data is constant, empty or near constant.")
            # Both dataset partitions NEED TO have the same number of bins.
            min_bw = min(bw_X, bw_Y);
            data_range = max(X.max() - X.min(), Y.max() - Y.min())
            nbins = max(1, int(np.ceil(data_range / min_bw)))

            if( nbins > 256):
                #print(f"X bin width: {bw_X}, Y bin width: {bw_Y}, X range: {X.max() - X.min()}, Y range: {Y.max() - Y.min()}")
                #print(f"nbins: {nbins}")
                if verbose:
                    print(f"Warning: Number of bins is too high. Setting nbins to 256.") 
                #+ f"Consider using a smaller bin width or a different binning method.")
                
                nbins = 256

        elif isinstance(bins, int) or isinstance(bins, float):
            nbins = int(bins)
        else:
            raise TypeError(f"Invalid type for bins: {type(bins)}. Use 'freedman', 'scotts' or integer value.")
    except ValueError as e:
            print(f"Error: {e} for X and Y. Check if data is constant, empty or near constant. MI by definition is 0.")
            print(f"X mean: {np.mean(X)}, X std: {np.std(X)}"); print(f"Y mean: {np.mean(Y)}, Y std: {np.std(Y)}")
            print(f"X 25th percentile: {np.percentile(X, 25)}, X 75th percentile: {np.percentile(X, 75)}")
            print(f"Y 25th percentile: {np.percentile(Y, 25)}, Y 75th percentile: {np.percentile(Y, 75)}")
            print(f"X min: {X.min()}, X max: {X.max()}") ; print(f"Y min: {Y.min()}, Y max: {Y.max()}")
            return 0
    except TypeError as e:
            print(f"Error: {e} for X and Y. Setting nBins to {256}.")
            nbins = 256
        
    '''# OLD APPROACH: Compute the histogram of the joint distribution of X and Y.
    c_XY = np.histogram2d(X, Y, bins=nbins)[0]
    MI = mutual_info_score(None, None, contingency=c_XY)
    '''

    # NEW APPROACH: Use label vectors from the binned data to compute MI.
    # Digitise the data in X and Y using the histogram bin edges.
    X_digitised = np.digitize(X, bins=np.histogram_bin_edges(X, bins=nbins))
    Y_digitised = np.digitize(Y, bins=np.histogram_bin_edges(Y, bins=nbins))
    #print("Computing Mutual Information...")
    MI = mutual_info_score(X_digitised, Y_digitised)
    #'''
    #print(f"Mutual Information between X and Y: {MI} for bins = {nbins} using {bins} method.")
    return MI


def compute_normalised_mutual_information(X, Y, bins="scotts", verbose=False):
    # x, y should be 1D arrays, if not, flatten them
    X = X.flatten() if X.ndim > 1 else X; Y = Y.flatten() if Y.ndim > 1 else Y

    # If bins is numeric, use it as the number of bins. If bins is a string, set bins using Scott's method.
    try:
        if isinstance(bins, str):
            
            if bins == 'freedman':
                # Use Freedman-Diaconis rule to determine the number of bins.
                bw_X = freedman_diaconis_bin_width(X); bw_Y = freedman_diaconis_bin_width(Y);
            elif bins == 'scotts':
                # Use Scott's method to determine the number of bins.
                bw_X = scotts_bin_width(X); bw_Y = scotts_bin_width(Y);
            else:
                raise TypeError(f"Invalid binning method {bins}. Use 'freedman' or 'scotts' or integer value.")

            if np.isnan(bw_X) or np.isnan(bw_Y) or bw_X <= 0 or bw_Y <= 0:
                raise ValueError("Invalid bin width calculated. Check if data is constant, empty or near constant.")
            # Both dataset partitions NEED TO have the same number of bins.
            min_bw = min(bw_X, bw_Y);
            data_range = max(X.max() - X.min(), Y.max() - Y.min())
            nbins = max(1, int(np.ceil(data_range / min_bw)))

            if( nbins > 256):
                if verbose:
                    print(f"WARNING: Number of bins = {nbins} is too high. Setting nbins to 256.")# +
                      #f"Consider using a smaller bin width or a different binning method.")
                nbins = 256

        elif isinstance(bins, int) or isinstance(bins, float):
            nbins = int(bins)
        else:
            raise TypeError(f"Invalid type for bins: {type(bins)}. Use 'freedman', 'scotts' or integer value.")
    except ValueError as e:
            print(f"Error: {e} for X and Y. Check if data is constant, empty or near constant. MI by definition is 0.")
            print(f"X mean: {np.mean(X)}, X std: {np.std(X)}"); print(f"Y mean: {np.mean(Y)}, Y std: {np.std(Y)}")
            print(f"X 25th percentile: {np.percentile(X, 25)}, X 75th percentile: {np.percentile(X, 75)}")
            print(f"Y 25th percentile: {np.percentile(Y, 25)}, Y 75th percentile: {np.percentile(Y, 75)}")
            print(f"X min: {X.min()}, X max: {X.max()}") ; print(f"Y min: {Y.min()}, Y max: {Y.max()}")
            return 0
    except TypeError as e:
            print(f"Error: {e} for X and Y. Setting nBins to {256}.")
            nbins = 256

    X_digitised = np.digitize(X, bins=np.histogram_bin_edges(X, bins=nbins))
    Y_digitised = np.digitize(Y, bins=np.histogram_bin_edges(Y, bins=nbins))
    # Compute NMI
    NMI = normalized_mutual_info_score(X_digitised, Y_digitised)
    return NMI

def compute_adjusted_mutual_information(X, Y, bins="scotts", verbose=False):
    # x, y should be 1D arrays, if not, flatten them
    #print(f"Is X CUPY array? {isinstance(X, gpu._cupy.ndarray)}")
    #print(f"Is X NUMPY array? {isinstance(X, _np_base.ndarray)}")
    X = X.flatten() if X.ndim > 1 else X; Y = Y.flatten() if Y.ndim > 1 else Y
    #print("Working on AMI...")
    # Print the shapes of X and Y
    #print(f"X shape: {X.shape}, Y shape: {Y.shape}")
    #print(f"X 25th percentile: {np.percentile(X, 25)}, X 75th percentile: {np.percentile(X, 75)}")
    #print(f"Y 25th percentile: {np.percentile(Y, 25)}, Y 75th percentile: {np.percentile(Y, 75)}")
    #print(f"X min: {X.min()}, X max: {X.max()}") ; print(f"Y min: {Y.min()}, Y max: {Y.max()}")

    # If bins is numeric, use it as the number of bins. If bins is a string, set bins using Scott's method.
    try:
        if isinstance(bins, str):
            #print(f"Using {bins} method for binning...")
            if bins == 'freedman':
                # Use Freedman-Diaconis rule to determine the number of bins.
                bw_X = freedman_diaconis_bin_width(X); bw_Y = freedman_diaconis_bin_width(Y);
            elif bins == 'scotts':
                # Use Scott's method to determine the number of bins.
                #print(f"Is X CUPY array? {isinstance(X, gpu._cupy.ndarray)}")
                #print(f"Is X NUMPY array? {isinstance(X, _np_base.ndarray)}")
                #print(f"Is Y CUPY array? {isinstance(Y, gpu._cupy.ndarray)}")
                #print(f"Is Y NUMPY array? {isinstance(Y, _np_base.ndarray)}")
                bw_X = scotts_bin_width(X); bw_Y = scotts_bin_width(Y);
            else:
                raise TypeError(f"Invalid binning method {bins}. Use 'freedman' or 'scotts' or integer value.")

            if np.isnan(bw_X) or np.isnan(bw_Y) or bw_X <= 0 or bw_Y <= 0:
                raise ValueError("Invalid bin width calculated. Check if data is constant, empty or near constant.")
            # Both dataset partitions NEED TO have the same number of bins.
            #print("Ascertaining minimum bin width for X and Y...")
            min_bw = min(bw_X, bw_Y);
            #print(f"Getting data range for X and Y...")
            data_range = max(X.max() - X.min(), Y.max() - Y.min())
            
            nbins = max(1, int(np.ceil(data_range / min_bw)))
            #print(f"Calculating number of bins... as nbins = {nbins}.")
            if( nbins > 256):
                if verbose:
                    print(f"WARNING: Number of bins = {nbins} is too high. Setting nbins to 256") 
                      #+ f"Consider using a smaller bin width or a different binning method.")
                nbins = 256

        elif isinstance(bins, int) or isinstance(bins, float):
            nbins = int(bins)
        else:
            raise TypeError(f"Invalid type for bins: {type(bins)}. Use 'freedman', 'scotts' or integer value.")
    except ValueError as e:
            print(f"Error: {e} for X and Y. Check if data is constant, empty or near constant. MI by definition is 0.")
            print(f"X mean: {np.mean(X)}, X std: {np.std(X)}"); print(f"Y mean: {np.mean(Y)}, Y std: {np.std(Y)}")
            print(f"X 25th percentile: {np.percentile(X, 25)}, X 75th percentile: {np.percentile(X, 75)}")
            print(f"Y 25th percentile: {np.percentile(Y, 25)}, Y 75th percentile: {np.percentile(Y, 75)}")
            print(f"X min: {X.min()}, X max: {X.max()}") ; print(f"Y min: {Y.min()}, Y max: {Y.max()}")
            return 0
    except TypeError as e:
            print(f"Error: {e} for X and Y. Setting nBins to {256}.")
            nbins = 256
            print(f"STACK TRACE for PID {os.getpid()}:")
            traceback.print_stack()
            print("=" * 80)
            #exit(1)


    # Discretize the data into bins
    #print(f"Digitising data into {nbins} bins...")
    X_digitised = np.digitize(X, bins=np.histogram_bin_edges(X, bins=nbins))
    Y_digitised = np.digitize(Y, bins=np.histogram_bin_edges(Y, bins=nbins))
    # Compute AMI
    #print(f"Computing AMI for X and Y with {nbins} bins...")
    AMI = adjusted_mutual_info_score(X_digitised, Y_digitised)
    #print(f"AMI for X and Y with {nbins} bins: {AMI}")
    return AMI


def compute_BiMoronsI(X, Y):
    X = X.flatten() if X.ndim > 1 else X; Y = Y.flatten() if Y.ndim > 1 else Y
    # Compute Morons'I Global Bivariate Statistic
    #print(f"Is X CUPY array? {isinstance(X, gpu._cupy.ndarray)}")
    #print(f"Calculating L for Bi-variate Moran's I ...")
    #L = int(np.round(np.sqrt(X.shape[0]))) AGNOSTIC, BUT HAS SEVERE PERFORMANCE ISSUES
    L = round(asnumpy(X.shape[0])**0.5)  # Use round to ensure L is an integer
    #print(f"Computing Bi-variate Moran's I for L = {L} ...")
    #Spatial_Weights = libpysal_weights.lat2W(L, L)
    Spatial_Weights= cached_libpysalweights_lat2W(L)
    MoranI = esda_moran.Moran_BV(X, Y, Spatial_Weights)
    #print(f"Bi-variate Moran's I: {MoranI.I}, p-value: {MoranI.p_sim}")
    return MoranI.I, MoranI.p_sim



''' # Summary of the function compute_NCC_2DFFT(...)
This function compute_ncc_fft computes the normalized cross-correlation (NCC) between two 2D arrays x and y using FFT.
It first normalizes the arrays by subtracting the mean and dividing by the standard deviation.'''
def compute_NCC_2DFFT(x, y, zero_norm=True):
    #start = time.perf_counter()
    std_x = np.std(x)
    std_y = np.std(y)
    if std_x < 1e-8 or std_y < 1e-8:
        return None, None, None, None
    if zero_norm:
        x = (x - np.mean(x)) / std_x
        y = (y - np.mean(y)) / std_y
    #print(f"Before fftconvolve - x: {'GPU' if is_gpu_array(x) else 'CPU'}, y: {'GPU' if is_gpu_array(y) else 'CPU'}")
    # NCC computed using FFT convolution, with the y array flipped.
    corr = signal.fftconvolve(x, y[::-1, ::-1], mode='full') / (x.shape[0]*x.shape[1])
    #corr = gpu._gpu_ndimage.convolve(x, y[::-1, ::-1])/ (x.shape[0]*x.shape[1])
    #print(f"FFT operation ran on: {'GPU (cupyx.signal)' if is_gpu_array(corr) else 'CPU (scipy.signal)'}")
    max_idx = np.unravel_index(np.argmax(corr), corr.shape)
    zero_shift_idx = (y.shape[0]-1, y.shape[1]-1)

    #execution_time = time.perf_counter() - start
    #backend = 'GPU' if is_gpu_array(corr) else 'CPU'
    
    #print(f"NCC computed in {execution_time:.4f}s on {backend}")

    return corr, corr[zero_shift_idx], np.max(corr), max_idx


def timed_compute_NCC_2DFFT(x, y, zero_norm=True):

    if gpu.GPU_AVAILABLE:
        pass
        #gpu._cupy.cuda.Stream.null.synchronize()
    start = time.perf_counter()
    std_x = np.std(x)
    std_y = np.std(y)
    if std_x < 1e-8 or std_y < 1e-8:
        return None, None, None, None
    if zero_norm:
        x = (x - np.mean(x)) / std_x
        y = (y - np.mean(y)) / std_y
    #print(f"Before fftconvolve - x: {'GPU' if is_gpu_array(x) else 'CPU'}, y: {'GPU' if is_gpu_array(y) else 'CPU'}")
    # NCC computed using FFT convolution, with the y array flipped.
    corr = signal.fftconvolve(x, y[::-1, ::-1], mode='full') / (x.shape[0]*x.shape[1])
    #corr = signal.fftconvolve(x, y[::-1, ::-1], mode='full') / (x.shape[0]*x.shape[1])
    # Pad the arrays to 2*max(shape) to avoid circular convolution effects.
    #x = np.pad(x, int(x.shape[0]/2), mode='constant', constant_values=0)
    #y = np.pad(y, int(y.shape[0]/2), mode='constant', constant_values=0)
    #corr = fft.ifft2(fft.fft2(x, (2*int(x.shape[0])-1, 2*int(x.shape[1]) -1)) 
    #                 * fft.fft2(y, (2*int(y.shape[0]) -1, 2*int(y.shape[1])-1)))/ (x.shape[0]*x.shape[1])
    #corr = gpu._gpu_ndimage.convolve(x, y[::-1, ::-1])/ (x.shape[0]*x.shape[1])
    #gpu._cupy.cuda.Device(0).synchronize()  # Wait for actual GPU completion
    #print(f"FFT operation ran on: {'GPU (cupyx.signal)' if is_gpu_array(corr) else 'CPU (scipy.signal)'}")
    max_idx = np.unravel_index(np.argmax(corr), corr.shape)
    zero_shift_idx = (y.shape[0]-1, y.shape[1]-1)

    if gpu.GPU_AVAILABLE and is_gpu_array(corr):
        #gpu._cupy.cuda.Stream.null.synchronize()
        pass

    execution_time = time.perf_counter() - start
    backend = 'GPU' if is_gpu_array(corr) else 'CPU'
    
    #print(f"NCC computed in {execution_time:.4f}s on {backend}")
    #timer_arr[i] = execution_time
    # Store the execution time in the timer array at index i

    return corr, corr[zero_shift_idx], np.max(corr), max_idx, execution_time



# Define the worker function to process each time value
def process_time_value(tval, files, T, crossfiles, L, col_labels, ext, exclude_col_labels, calc_AMI, calc_MI, calc_Morans, bins, verbose):
    # Initialize nested dictionaries for this time value
    #auto_NCC_dict = defaultdict(lambda: defaultdict(dict))
    #Structure: auto_NCC_dict[key][column][NCC_key] = value
    #cross_NCC_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    # Structure: cross_NCC_dict[key][column1][column2][NCC_key] = value
    #auto_ZNCC_dict = defaultdict(lambda: defaultdict(dict))
    #cross_ZNCC_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

    # Avoiding lambda functions for better readability and to avoid pickling issues
    
    auto_NCC_dict = defaultdict(nested_default_dict)
    #auto_NCC_dict = recursively_nested_default_dict(2) # 2 levels deep
    # Structure: cross_NCC_dict[key][column1][column2][NCC_key] = value
    cross_NCC_dict = defaultdict(double_nested_default_dict)
    auto_ZNCC_dict = defaultdict(nested_default_dict)  #recursively_nested_default_dict(2) # 2 levels deep
    cross_ZNCC_dict = defaultdict(double_nested_default_dict)

    # Precompiled regex patterns for T and R values
    T_intpattern = re.compile(r'T_(\d+)\b')
    T_floatpattern = re.compile(r'T_(\d+\.\d+)\b')
    R_pattern = re.compile(r'R_(\d+)')
    
    if calc_AMI:
        #auto_AMI_dict = defaultdict(lambda: defaultdict(dict)) ; cross_AMI_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        #auto_AMI_dict = recursively_nested_default_dict(2); cross_AMI_dict = recursively_nested_default_dict(3)
        auto_AMI_dict = defaultdict(nested_default_dict); cross_AMI_dict = defaultdict(double_nested_default_dict)
    else:
        auto_AMI_dict = None; cross_AMI_dict = None
        
    if calc_MI:
        #auto_MI_dict = defaultdict(lambda: defaultdict(dict)); cross_MI_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        #auto_MI_dict = recursively_nested_default_dict(2); cross_MI_dict = recursively_nested_default_dict(3)
        auto_MI_dict = defaultdict(nested_default_dict); cross_MI_dict = defaultdict(double_nested_default_dict)
    else:
        auto_MI_dict = None; cross_MI_dict = None
    
    if calc_Morans:
        auto_Morons_dict = defaultdict(nested_default_dict); cross_Morons_dict = defaultdict(double_nested_default_dict)
    else:
        auto_Morons_dict = None; cross_Morons_dict = None

    key = (T, tval, T - float(tval), len(files))  # Common key for all files.
    #safe_print(f"Generating 2D correlation data for T0 = {T}, T1 = {tval}, delay = {T - float(tval)}, L = {L} ....")
    print(f"Generating 2D correlation data for T0 = {T}, T1 = {tval}, delay = {T - float(tval)}, L = {L} ....")
    # Iterate over all files in files, sorted in ascending order of R using regex
    
    timer_autoarr = np.zeros(10, dtype=np.float64); auto_i=0;
    #timer_crossarr = np.zeros(10, dtype=np.float64); cross_j=0;
    timer_crossarr = np.zeros(10, dtype=np.float64); cross_j=0;

    for file in sorted_files_R(files, ascending=True):
        # Read the file.
        Rstr = R_pattern.search(file).group(0)
        #Rstr = re.findall(r'R_[\d]+', file)[0]
        Lfile = int(re.search(r'L_[\d]+', file).group(0).split("_")[1]) if re.search(r'L_[\d]+', file) else None 
        filetype = re.search(r'FRAME|GAMMA', file).group(0)
        
        if Lfile != L:
            #safe_print(f"Warning: L value in file {file} is {Lfile} but expected {L}. Skipping file.")
            print(f"Warning: L value in file {file} is {Lfile} but expected {L}. Skipping file.")
            continue

        if(ext == "csv"):
            df = pan.read_csv(file, header=0)
        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                #safe_print(f"Error: Non-standard extension for file {file} with error message: \n" + str(e))
                print(f"Error: Non-standard extension for file {file} with error message: \n" + str(e))
                return auto_NCC_dict, cross_NCC_dict, auto_ZNCC_dict, cross_ZNCC_dict, auto_AMI_dict, cross_AMI_dict, auto_MI_dict, cross_MI_dict

        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        df = df[[col for col in df.columns if col not in exclude_col_labels]]

        df_Znorm = gen_zer0mean_normalised_df(df)  # Generate 0-normalised columns for df.
            
        # Get crossfiles with the same T value as tval and R value as Rstr.
        
        crossfile_tval = [x for x in crossfiles if ((m := T_intpattern.search(x)) and m.group(1) == str(tval))
                        or ((m := T_floatpattern.search(x)) and m.group(1) == str(tval))]
        """ OLD UNOPTIMISED CODE:
        crossfile_tval = [x for x in crossfiles
                        if ((m := re.findall(r'T_(\d+)', x)) and m[0] == str(tval)) 
                        or ((m := re.findall(r'T_([\d]+\.[\d]+)', x)) and m[0] == str(tval))]
        """
        crossfile_tval = [x for x in crossfile_tval if re.search(r'R_[\d]+', x).group(0) == Rstr]# and re.search(r'FRAME|GAMMA', x).group(0) == filetype]
        # Get filetype of crossfile_tval and store it in filetype_crossfile.
        
        # If no or more than 1 are found, skip the column.
        if(tval == T):
            # If tval == T, just do auto-corr with same file.
            crossfile_tval = [file]
        
        
        if len(crossfile_tval) > 2 or len(crossfile_tval) == 0:
            print(f"Warning: No or more than 2 crossfile found for T = {tval} and R = {Rstr}. Skipping file.")
            #safe_print(f"Crossfiles found: {crossfile_tval}")
            print(f"Crossfiles found: {crossfile_tval}")
            continue
        

        # Read the crossfile and zero-normalise all columns in the crossfile.
        for crossfile in crossfile_tval:
            # Get the filetype of the crossfile.
            filetype_crossfile = re.search(r'FRAME|GAMMA', crossfile).group(0)
            #timer_autoarr = np.zeros(5, dtype=np.float64); auto_i=0;
            
            
            
            if(ext == "csv"):
                df_cross = pan.read_csv(crossfile, header=0)
            else:
                try:
                    df_cross = pan.read_table(crossfile, header=0)
                except Exception as e:
                    #safe_print(f"Error: Non-standard extension for file {crossfile} with error message: \n" + str(e))
                    print(f"Error: Non-standard extension for file {crossfile} with error message: \n" + str(e))
                    return auto_NCC_dict, cross_NCC_dict, auto_ZNCC_dict, cross_ZNCC_dict, auto_AMI_dict, cross_AMI_dict, auto_MI_dict, cross_MI_dict
                    
            # Modify column names in df_cross to remove leading and trailing whitespaces.
            df_cross.columns = [col.strip() for col in df_cross.columns]
            # Remove exclude_col_labels from df_cross.
            df_cross = df_cross[[col for col in df_cross.columns if col not in exclude_col_labels]]
            df_cross_Znorm = gen_zer0mean_normalised_df(df_cross)  # Generate 0-normalised columns for df_cross.

            #safe_print(f"Calculating 2D correlation for FOCUS file {os.path.basename(file)} with CROSSfile {os.path.basename(crossfile)} at R = {Rstr}...")
            print(f"Calculating 2D correlation for FOCUS file {os.path.basename(file)} with CROSSfile {os.path.basename(crossfile)} at R = {Rstr}...")

            # First auto-NCC and auto-ZNCC: Iterate over all df.columns and calculate the NCC for each column in df.
            for col in df_Znorm.columns:
                # If filetypes are different, break the loop and proceed to cross-correlations.
                # This is to avoid cross-correlating with different filetypes.
                if filetype != filetype_crossfile:
                    break
                # Calculate the NCC for each column in df with the same column in df_cross.
                if col not in df_cross_Znorm.columns:
                    #safe_print(f"Warning: Column {col} not found in crossfile {os.path.basename(crossfile)}")
                    print(f"Warning: Column {col} not found in crossfile {os.path.basename(crossfile)}")
                    continue

                # Also check if the column is entirely zero or flat. If so, skip the column.
                if all(df_Znorm[col] == 0) or all(df_cross_Znorm[col] == 0):
                    #safe_print(f"Warning: Entirely zero column {col} in FOCUS file {os.path.basename(file)}" + 
                    print(f"Warning: Entirely zero column {col} in FOCUS file {os.path.basename(file)}" + 
                        f"or CROSSfile {os.path.basename(crossfile)} for T0 = {T}, T1 = {tval} ... Skipping column.")
                        
                    continue

                #safe_print(f"Calculating 2D AUTO Correlation for column {col}...")
                #print(f"Calculating 2D AUTO Correlation for column {col}...")

                #process = psutil.Process(os.getpid())
                #start_cpu_times = process.cpu_times()
                #start_threads = threading.active_count()
                # High-resolution timing
                #start_time = time.perf_counter()
                
                if(calc_AMI):
                    # Calculate AMI for each column in df with the same column in df_cross.
                    AMI = compute_adjusted_mutual_information(asarray(df[col].values), asarray(df_cross[col].values), bins=bins, verbose=verbose)
                    AMI_key = f"AMI{{{col}}}_{Rstr}"
                    auto_AMI_dict[key][col][AMI_key] = AMI
                    
                if(calc_MI):
                    # Calculate MI for each column in df with the same column in df_cross.
                    MI = compute_mutual_information(asarray(df[col].values), asarray(df_cross[col].values), bins=bins, verbose=verbose)
                    MI_key = f"MI{{{col}}}_{Rstr}"
                    auto_MI_dict[key][col][MI_key] = MI
                
                # Post-execution timing
                #end_time = time.perf_counter()
                #end_cpu_times = process.cpu_times()
                
                # Calculate metrics
                #wall_time = end_time - start_time
                #cpu_time = (end_cpu_times.user - start_cpu_times.user) + (end_cpu_times.system - start_cpu_times.system)
                #cpu_efficiency = cpu_time / wall_time if wall_time > 0 else 0
                
                # Print timing info
                #print(f"(PID {os.getpid()}) Wall time for AUTO AMI+MI 2D: {wall_time:.4f}s \t CPU time for AUTO AMI+MI 2D: {cpu_time:.4f}s")
                #print(f"CPU time for AUTO AMI+MI 2D: {cpu_time:.4f}s") 
                #print(f"On PID {os.getpid()} - CPU efficiency for AUTO AMI+MI 2D: {cpu_efficiency:.2f} (1.0 = perfect CPU usage)")
                #print(f"(PID {os.getpid()}) GIL impact estimate for AUTO AMI+MI 2D: {'High' if cpu_efficiency < 0.8 else 'Low'}")
                #print(f"Threads active: {threading.active_count()}")

                #process = psutil.Process(os.getpid())
                #start_cpu_times = process.cpu_times()
                #start_threads = threading.active_count()
                # High-resolution timing
                #start_time = time.perf_counter()
                
                if(calc_Morans):
                    # Calculate Morans I for each column in df with the same column in df_cross.
                    MoronsI, p_sim = compute_BiMoronsI(df_Znorm[col].values, df_cross_Znorm[col].values)
                    MoransI_key = f"BVMoransI{{{col}}}_{Rstr}"; MoransI_Pkey = f"BVMoran-P{{{col}}}_{Rstr}"
                    auto_Morons_dict[key][col][MoransI_key] = MoronsI
                    auto_Morons_dict[key][col][MoransI_Pkey] = p_sim

                # Post-execution timing
                #end_time = time.perf_counter()
                #end_cpu_times = process.cpu_times()
                
                # Calculate metrics
                #wall_time = end_time - start_time
                #cpu_time = (end_cpu_times.user - start_cpu_times.user) + (end_cpu_times.system - start_cpu_times.system)
                #cpu_efficiency = cpu_time / wall_time if wall_time > 0 else 0
                
                # Print timing info
                #print(f"(PID {os.getpid()})Wall time for AUTO MORONS 2D: {wall_time:.4f}s")
                #print(f"(PID {os.getpid()})CPU time for AUTO MORONS 2D: {cpu_time:.4f}s") 
                #print(f"(PID {os.getpid()})CPU efficiency for AUTO MORONS 2D: {cpu_efficiency:.2f} (1.0 = perfect CPU usage)")
                #print(f"(PID {os.getpid()}) GIL impact estimate for AUTO MORONS 2D: {'High' if cpu_efficiency < 0.8 else 'Low'}")
                #print(f"Threads active: {threading.active_count()}")

                x = asarray(df_Znorm[col].values).reshape(int(np.sqrt(len(df_Znorm))), -1)
                y = asarray(df_cross_Znorm[col].values).reshape(int(np.sqrt(len(df_cross_Znorm))), -1)

                NCC_key = f"NCC{{{col}}}_{Rstr}"
                NCC_index_key = f"NCC-Index{{{col}}}_{Rstr}"
                ZNCC_key = f"ZNCC{{{col}}}_{Rstr}"
                #print(f"Calculating 2D NCC for column {col}...")

                #process = psutil.Process(os.getpid())
                #start_cpu_times = process.cpu_times()
                #start_threads = threading.active_count()
                # High-resolution timing
                #start_time = time.perf_counter()

                ncc_map, zncc, peak, peak_idx = compute_NCC_2DFFT(x, y)
                #ncc_map, zncc, peak, peak_idx, exec_time = timed_compute_NCC_2DFFT(x, y)
                #timer_autoarr[auto_i % len(timer_autoarr)] = exec_time
                #auto_i += 1
                #if auto_i% len(timer_autoarr) == 0:
                #    print(f"AVG AUTO NCC 2D execution time: {np.mean(timer_autoarr):.4f}s +- {np.std(timer_autoarr):.4f}s for {len(timer_autoarr)} iterations.")

                #print("Done with AUTO NCC 2D calculation.")
                # Post-execution timing
                #end_time = time.perf_counter()
                #end_cpu_times = process.cpu_times()
                
                # Calculate metrics
                #wall_time = end_time - start_time
                #cpu_time = (end_cpu_times.user - start_cpu_times.user) + (end_cpu_times.system - start_cpu_times.system)
                #cpu_efficiency = cpu_time / wall_time if wall_time > 0 else 0
                
                # Print timing info
                #print(f"(PID {os.getpid()})Wall time for AUTO NCC 2D: {wall_time:.4f}s")
                #print(f"(PID {os.getpid()})CPU time for AUTO NCC 2D: {cpu_time:.4f}s") 
                #print(f"(PID {os.getpid()})CPU efficiency for AUTO NCC 2D: {cpu_efficiency:.2f} (1.0 = perfect CPU usage)")
                #print(f"(PID {os.getpid()})GIL impact estimate for AUTO NCC 2D: {'High' if cpu_efficiency < 0.8 else 'Low'}")
                #print(f"Threads active: {threading.active_count()}")
                if zncc is None:
                    #safe_print(f"Warning: Skipping column {col} due to flat or zero data in correlation.")
                    print(f"Warning: Skipping column {col} due to flat or zero data in correlation.")
                    continue
                    
                # Store the NCC data in auto_NCC_dict, and the ZNCC data in auto_ZNCC_dict.
                #print(f"Storing NCC and ZNCC for column {col}...")
                auto_NCC_dict[key][col][NCC_key] = peak
                auto_NCC_dict[key][col][NCC_index_key] = int(peak_idx[0] * L + peak_idx[1])  # Convert to 1D index.
                auto_ZNCC_dict[key][col][ZNCC_key] = zncc
            # End of auto-CORR loop.
            # Similar process for cross-CORR.

            #safe_print(f"Calculating CROSS Correlation b/w {filetype} and {filetype_crossfile} TYPE FILES...")
            if verbose:
                print(f"Calculating CROSS Correlation b/w {filetype} and {filetype_crossfile} TYPE FILES...")
            #gpu.init_daskworker_cuda_context() # Get CUDA contexts
            # Next cross-NCC and cross-ZNCC.
            for i, col1 in enumerate(df_Znorm.columns):
                if all(df_Znorm[col1] == 0):
                    #safe_print(f"Warning: Entirely zero column {col1} in FOCUS file {os.path.basename(file)} for T0 = {T} ... Skipping column.")
                    print(f"Warning: Entirely zero column {col1} in FOCUS file {os.path.basename(file)} for T0 = {T} ... Skipping column.")
                    continue

                for j, col2 in enumerate(df_cross_Znorm.columns):
                    if col1 == col2:
                        continue

                    # Also check if the column is entirely zero or flat. If so, skip the column.
                    if all(df_cross_Znorm[col2] == 0):
                        print(f"Warning: Entirely zero COL: {col2} in CROSSfile {os.path.basename(crossfile)} for T0 = {T}, T1 = {tval} ... Skipping column.")
                        continue

                    #safe_print(f"Calculating 2D CROSS Correlation for columns {col1} VS {col2} ...")
                    if verbose:
                        print(f"Calculating 2D CROSS Correlation for columns {col1} VS {col2} ...")

                    x = asarray(df_Znorm[col1].values).reshape(int(np.sqrt(len(df_Znorm))), -1)
                    y = asarray(df_cross_Znorm[col2].values).reshape(int(np.sqrt(len(df_cross_Znorm))), -1)

                    #process = psutil.Process(os.getpid())
                    #start_cpu_times = process.cpu_times()
                    #start_threads = threading.active_count()
                    # High-resolution timing
                    #start_time = time.perf_counter()

                    #ncc_map, zncc, peak, peak_idx = compute_NCC_2DFFT(x, y)
                    ncc_map, zncc, peak, peak_idx, exec_time = timed_compute_NCC_2DFFT(x, y)
                    timer_crossarr[cross_j % len(timer_crossarr)] = exec_time
                    cross_j += 1
                    if cross_j % len(timer_crossarr) == 0:
                        print(f"AVG CROSS NCC 2D execution time: {np.mean(timer_crossarr):.4f}s +- {np.std(timer_crossarr):.4f}s for {len(timer_crossarr)} iterations.")
                    #print("Done with CROSS NCC 2D calculation.")
                    # Post-execution timing
                    #end_time = time.perf_counter()
                    #end_cpu_times = process.cpu_times()
                    
                    # Calculate metrics
                    #wall_time = end_time - start_time
                    #cpu_time = (end_cpu_times.user - start_cpu_times.user) + (end_cpu_times.system - start_cpu_times.system)
                    #cpu_efficiency = cpu_time / wall_time if wall_time > 0 else 0
                    
                    # Print timing info
                    #print(f"(PID {os.getpid()})Wall time for CROSS NCC 2D: {wall_time:.4f}s")
                    #print(f"(PID {os.getpid()})CPU time for CROSS NCC 2D: {cpu_time:.4f}s") 
                    #print(f"(PID {os.getpid()})CPU efficiency for CROSS NCC 2D: {cpu_efficiency:.2f} (1.0 = perfect CPU usage)")
                    #print(f"(PID {os.getpid()})GIL impact estimate for CROSS NCC 2D: {'High' if cpu_efficiency < 0.8 else 'Low'}")
                    #print(f"Threads active: {threading.active_count()}")

                    NCC_index_key = f"NCC-Index{{{col1};{col2}}}_{Rstr}"
                    NCC_key = f"NCC{{{col1};{col2}}}_{Rstr}"
                    ZNCC_key = f"ZNCC{{{col1};{col2}}}_{Rstr}"

                    if zncc is None:
                        #safe_print(f"Warning: Skipping columns {col1}, {col2} due to flat or zero data.")
                        print(f"Warning: Skipping columns {col1}, {col2} due to flat or zero data.")
                        continue

                    cross_NCC_dict[key][col1][col2][NCC_key] = peak
                    cross_NCC_dict[key][col1][col2][NCC_index_key] = int(peak_idx[0] * L + peak_idx[1])
                    # Convert to 1D index.
                    cross_ZNCC_dict[key][col1][col2][ZNCC_key] = zncc

                    if(calc_AMI):
                        # Calculate AMI for each column in df with the same column in df_cross.
                        AMI = compute_adjusted_mutual_information(asarray(df[col1].values), asarray(df_cross[col2].values), bins=bins, verbose=verbose)
                        AMI_key = f"AMI{{{col1};{col2}}}_{Rstr}"
                        cross_AMI_dict[key][col1][col2][AMI_key] = AMI
                        
                    if(calc_MI):
                        # Calculate MI for each column in df with the same column in df_cross.
                        MI = compute_mutual_information(asarray(df[col1].values), asarray(df_cross[col2].values), bins=bins, verbose=verbose)
                        MI_key = f"MI{{{col1};{col2}}}_{Rstr}"
                        cross_MI_dict[key][col1][col2][MI_key] = MI

                    #process = psutil.Process(os.getpid())
                    #start_cpu_times = process.cpu_times()
                    #start_threads = threading.active_count()
                    # High-resolution timing
                    #start_time = time.perf_counter()

                    if(calc_Morans):
                        # Calculate Morans I for each column in df with the same column in df_cross.
                        MoronsI, p_sim = compute_BiMoronsI(df_Znorm[col1].values, df_cross_Znorm[col2].values)
                        MoransI_key = f"BVMoransI{{{col1};{col2}}}_{Rstr}"; MoransI_Pkey = f"BVMoran-P{{{col1};{col2}}}_{Rstr}"
                        cross_Morons_dict[key][col1][col2][MoransI_key] = MoronsI
                        cross_Morons_dict[key][col1][col2][MoransI_Pkey] = p_sim

                    # Post-execution timing
                    #end_time = time.perf_counter()
                    #end_cpu_times = process.cpu_times()
                    
                    # Calculate metrics
                    #wall_time = end_time - start_time
                    #cpu_time = (end_cpu_times.user - start_cpu_times.user) + (end_cpu_times.system - start_cpu_times.system)
                    #cpu_efficiency = cpu_time / wall_time if wall_time > 0 else 0
                    
                    # Print timing info
                    #print(f"(PID {os.getpid()})Wall time for CROSS MORAN 2D: {wall_time:.4f}s")
                    #print(f"(PID {os.getpid()})CPU time for CROSS MORAN 2D: {cpu_time:.4f}s") 
                    #print(f"(PID {os.getpid()})CPU efficiency for CROSS MORAN 2D: {cpu_efficiency:.2f} (1.0 = perfect CPU usage)")
                    #print(f"(PID {os.getpid()})GIL impact estimate for CROSS MORAN 2D: {'High' if cpu_efficiency < 0.8 else 'Low'}")
                    #print(f"Threads active: {threading.active_count()}")

    
    return auto_NCC_dict, cross_NCC_dict, auto_ZNCC_dict, cross_ZNCC_dict, auto_AMI_dict, cross_AMI_dict, auto_MI_dict, cross_MI_dict, auto_Morons_dict, cross_Morons_dict


'''# Summary of the function process_time_value_joblib_wrapper(...)
This function is a Joblib-compatible wrapper for process_time_value. 
Handles GPU arrays and returns nested defaultdicts.
It sets up a different CUDA stream for each worker if available, allowing for better isolation and performance.
It also catches exceptions and returns empty defaultdict structures on error.'''
def process_time_value_joblib_wrapper(tval, files, T, crossfiles, L, col_labels, ext,
                                    exclude_col_labels, calc_AMI, calc_MI, calc_Morans, 
                                    bins, verbose, worker_id=None):
    try:
        # Optional: Set thread-specific GPU context if available
        if gpu.GPU_AVAILABLE and worker_id is not None and not gpu.DASK_AVAILABLE:
            # Recall cupy can
            # Each worker can have its own stream for better isolation
            with gpu._cupy.cuda.Stream(non_blocking=True):
                result = process_time_value(tval, files, T, crossfiles, L, col_labels, ext,
                                          exclude_col_labels, calc_AMI, calc_MI, calc_Morans, 
                                          bins, verbose)
        else:
            result = process_time_value(tval, files, T, crossfiles, L, col_labels, ext,
                                      exclude_col_labels, calc_AMI, calc_MI, calc_Morans, 
                                      bins, verbose)
        
        return result
        
    except Exception as e:
        warnings.warn(f"Error processing time value {tval}: {e}")
        # Return empty defaultdict structure on error
        #print(f"STACK TRACE for PID {os.getpid()}:")
        #traceback.print_stack()
        #print("=" * 50)
        #exit(1)
        empty_result = (
            defaultdict(nested_default_dict), #auto_NCC_dict
            defaultdict(double_nested_default_dict), #cross_NCC_dict
            defaultdict(nested_default_dict), #auto_ZNCC_dict
            defaultdict(double_nested_default_dict), #cross_ZNCC_dict
            defaultdict(nested_default_dict) if calc_AMI else None, #auto_AMI_dict
            defaultdict(double_nested_default_dict) if calc_AMI else None, #cross_AMI_dict
            defaultdict(lambda: defaultdict(dict)) if calc_MI else None, #auto_MI_dict
            defaultdict(lambda: defaultdict(lambda: defaultdict(dict))) if calc_MI else None, #cross_MI_dict
            defaultdict(lambda: defaultdict(dict)) if calc_Morans else None, #auto_Morons_dict
            defaultdict(lambda: defaultdict(lambda: defaultdict(dict))) if calc_Morans else None #cross_Morons_dict
        )
        return empty_result


''' # Summary of the function gen_2DCorr_data(...)
This function gen_NCC_data will generate auto- and cross-correlation data for the files in files.
These files are named in the format: FRAME_T_{T}_a_{a_val}_R_{R}.csv or GAMMA_T_{T}_a_{a_val}_R_{R}.csv.
By iterating over the files in increasing values of R, the column-specific time auto-NCC is generated by auto-correlation of each column in the file,
with a matching column in crossfiles (if provided), with the same R value as the file.
The cross-NCC on the other hand consists of two steps:
1. Pure spatial cross-correlation (at the same time): Cross-correlation of the data in each included column with the data  in another column in the same file.
2. Spatio-temporal cross-correlation of the data in each included column with the data in another column in a different file (crossfiles) with the same R value as the file.
The function will generate the NCC data for each column in the files for each file (replicate).
The auto-NCC (for best shifts with highest peaks) data is stored in a df with the following column structure:
"t0", "t1", "t-delay", "Rmax", "AVG[NCC[{var}]], ...., "NCC[{var}]_R_0", "NCC-Index[{var}]_R_0" ...., "NCC[{var}]_R_{R_max}", "NCC-Index[{var}]_R_{R_max}"
where {var} is one of the species names NOT in exclude_col_labels, Rmax is the maximum R value in the files, t0 is the time of the first frame, t1 is the time of the compared frame, 
and t-delay is the time delay between the two frames.
The cross-NCC (for best shifts with highest peaks) data is stored in a df with the following column structure:
"t0", "t1", "t-delay", "Rmax", "AVG[NCC[{var1}][{var2}]]", ...., "AVG[NCC[{var1}][{varN}]]" , "AVG[NCC[{var2}][{var3}]]", ...., "AVG[NCC[{vari}][{vark}]]" , ...., AVG[NCC[{varN-1}][{varN}]],
"NCC[{var1}][{var2}]_R_0" , "NCC-Index[{var1}][{var2}]_R_0" ...., "NCC[{var1}][{var2}]_R_{R_max}", "NCC-Index[{var1}][{var2}]_R_{R_max}" ...., 
"NCC[{varN-1}][{varN}]_R_{R_max}", "NCC-Index[{varN-1}][{varN}]_R_{R_max}" ...., "NCC[{varN-1}][{varN}]_R_{R_max}", "NCC-Index[{varN-1}][{varN}]_R_{R_max}"
where {var1}, {var2}, ... {varN} are the species names NOT in exclude_col_labels, Rmax is the maximum R value in the files, t0 is the time of the first frame, t1 is the time of the compared frame,
and t-delay is the time delay between the two frames, and i < k (the order of columns in the files).
The function will also generate the mean NCC data across all files for each column in the files, which is stored in AVG[NCC[{var}]] and AVG[NCC[{var1}][{var2}]] for auto- and cross-NCC respectively.

Similarly, computes the Adjusted Mutual Information (AMI) and Mutual Information (MI) for each column in the files for each file (replicate)
in the same way as NCC, and stores them in the same column structure as above.

Also the 2D Pearson correlation coefficient in each case corresponds to the zero-shift NCC, store these values in auto-ZNCC and cross-ZNCC dfs with the same column structure as above.

NOTE: ZNCC is the zero-normalised cross-correlation, so columns in the files are normalised to zero mean and unit variance before calculating the ZNCC.

INPUTS:
files: List of files which are named in the format: FRAME_T_{T}_a_{a_val}_R_{R}.csv or GAMMA_T_{T}_a_{a_val}_R_{R}.csv.
       These are the focus files for which the auto- and cross-correlation data is to be generated.
crossfiles: List of all crossfiles which are named in the format: FRAME_T_{T}_a_{a_val}_R_{R}.csv or GAMMA_T_{T}_a_{a_val}_R_{R}.csv.
    These are the potential files against which auto (over time) and cross-correlation is to be generated.
pathtodir: Path to the directory where the files are located.
ext: Extension of the files. Default is "csv".
exclude_col_labels: List of column labels that are excluded from the analysis, if they exist. Default is ["a_c", "x", "L", "GAM[P(x; t)]"].
calc_AMI, calc_MI, calc_Morans: Boolean values to indicate whether to calculate AMI, MI and Morans I respectively. 
    Default is True for AMI, False for MI and Morans I.
bins: Binning method to use for calculating AMI and MI. Default is "scotts" (Scott's method), 
    BUT can also be set to "freedman" (Freedman-Diaconis rule) or an integer value for the number of bins.
verbose: Boolean value to indicate whether to print verbose output. Default is False.
ncores: Number of cores to use for parallel processing. Default is 1.
    If ncores > 1, the function will use multiprocessing to parallelise the computation of NCC data across files.
'''


def gen_2DCorr_data(files, T, matchT=[], crossfiles=[], pathtodir="", ext="csv", exclude_col_labels= ["a_c", "x", "L", "GAM[P(x; t)]"], 
                    calc_AMI = True, calc_MI = False, calc_Morans=False, bins="scotts", verbose=False, ncores=1):

    # Store unique column names present in files in col_labels. Strip elements of any leading or trailing whitespaces.
    col_labels = [] #L = 0;
    # Iterate over other files in files, and store unique column names in col_labels.
    for file in files:
        if(ext == "csv"):
            df = pan.read_csv(file, header=0)
        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                print("Error: Non-standard extension for file " + pathtodir +"/" + file + " with error message: \n" + str(e))
                return None, None, None, None, None, None, None, None
        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        col_labels.extend([col for col in df.columns if col not in exclude_col_labels])

        if file == files[0]:
            # Get L from the number of rows in the first file.
            #L = int(np.sqrt(len(df))) #int(np.sqrt(len(df[df.columns[0]])))
            L = round(asnumpy(len(df)**0.5))
            print(f"Found L = {L}.")
    
    # Remove duplicates and exclude_col_labels from col_labels.
    col_labels = list(set(col_labels))
    col_labels = [col for col in col_labels if col not in exclude_col_labels]

    auto_NCC_cols = ["t0", "t1", "t-delay", "Rmax"] + ["AVG[NCC{" + col + "}]" for col in col_labels] + [y + col + "}_R_" + str(i) for col 
                                                                                                      in col_labels for y in ["NCC{", "NCC-Index{"] for i in range(0, len(files))]
    cross_NCC_cols = ["t0", "t1", "t-delay", "Rmax"] + ["AVG[NCC{" + col1 + ";" + col2 + "}]" for col1 
                    in col_labels for col2 in col_labels if col1 != col2] + [y + col1 + ";" + col2 + "}_R_" + str(i) 
                    for col1 in col_labels for col2 in col_labels if col1 != col2 for y in ["NCC{", "NCC-Index{"]  for i in range(0, len(files))]

    df_auto_NCC = pan.DataFrame(columns=auto_NCC_cols)
    df_cross_NCC = pan.DataFrame(columns=cross_NCC_cols)

    df_auto_ZNCC = pan.DataFrame(columns=auto_NCC_cols) # Stores ZNCC (NCC(0,0)) OR 2D Pearson correlation coefficient for each column in the files.
    df_cross_ZNCC = pan.DataFrame(columns=cross_NCC_cols)

    if calc_AMI:
        # Drop all columns containing "NCC-Index" from df_auto_AMI, and replace all occurrences of "NCC" with "AMI".
        auto_AMI_cols = [col.replace("NCC", "AMI") for col in auto_NCC_cols if "NCC-Index" not in col]
        cross_AMI_cols = [col.replace("NCC", "AMI") for col in cross_NCC_cols if "NCC-Index" not in col]
        df_auto_AMI = pan.DataFrame(columns=auto_AMI_cols)
        df_cross_AMI = pan.DataFrame(columns=cross_AMI_cols)
    else:
        df_auto_AMI = None
        df_cross_AMI = None

    if calc_MI:
        # Drop all columns containing "NCC-Index" from df_auto_MI, and replace all occurrences of "NCC" with "MI".
        auto_MI_cols = [col.replace("NCC", "MI") for col in auto_NCC_cols if "NCC-Index" not in col]
        cross_MI_cols = [col.replace("NCC", "MI") for col in cross_NCC_cols if "NCC-Index" not in col]
        df_auto_MI = pan.DataFrame(columns=auto_MI_cols)
        df_cross_MI = pan.DataFrame(columns=cross_MI_cols)
    else:
        df_auto_MI = None
        df_cross_MI = None

    if calc_Morans:
        # Drop all columns containing "NCC-Index" from df_auto_MoransI, and replace all occurrences of "NCC" with "MoransI".
        auto_MoransI_cols = [col.replace("NCC", "BVMoransI") for col in auto_NCC_cols if "NCC-Index" not in col]
        cross_MoransI_cols = [col.replace("NCC", "BVMoransI") for col in cross_NCC_cols if "NCC-Index" not in col]
        df_auto_MoransI = pan.DataFrame(columns=auto_MoransI_cols)
        df_cross_MoransI = pan.DataFrame(columns=cross_MoransI_cols)
    else:
        df_auto_MoransI = None
        df_cross_MoransI = None

    # Initialize empty dictionaries to collect results from parallel processing
    all_auto_NCC_dict = defaultdict(lambda: defaultdict(dict))
    all_cross_NCC_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    all_auto_ZNCC_dict = defaultdict(lambda: defaultdict(dict))
    all_cross_ZNCC_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    
    if calc_AMI:
        all_auto_AMI_dict = defaultdict(lambda: defaultdict(dict))
        all_cross_AMI_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    else:
        all_auto_AMI_dict = None
        all_cross_AMI_dict = None
        
    if calc_MI:
        all_auto_MI_dict = defaultdict(lambda: defaultdict(dict))
        all_cross_MI_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    else:
        all_auto_MI_dict = None
        all_cross_MI_dict = None

    if calc_Morans:
        all_auto_MoransI_dict = defaultdict(lambda: defaultdict(dict))
        all_cross_MoransI_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    else:
        all_auto_MoransI_dict = None
        all_cross_MoransI_dict = None

    # Ensure ncores doesn't exceed the number of time values or available CPU cores
    available_cores = cpu_count()
    ncores = min(ncores, len(matchT), available_cores)
    print(f"Using {ncores} cores out of {available_cores} available cores to process {len(matchT)} time values")
    

    
    if gpu.DASK_AVAILABLE:
        #dask_loccluster = dask_LocalCluster(n_workers=ncores, threads_per_worker=1, memory_limit='2GB')
        mem_limit = "3GB"
        daskclient = dask_Client(processes=True, n_workers=ncores, threads_per_worker=1, memory_limit='3GB',
            worker_kwargs={'death_timeout': 60,  'preload': []}, timeout=60)
        # Each dask worker gets single thread with its own CUDA context!
        print(f"Using Dask LocalCluster with {ncores} workers and 1 thread per worker, with {mem_limit} memory limit per worker.")
        if gpu.GPU_AVAILABLE:
            daskclient.register_worker_plugin(gpu.setup_daskworker_gpu_context)
        # Assigns local CUDA context to each worker.
        joblib_backend = 'dask'
        print("Dashboard URL: ", daskclient.dashboard_link)
        # Pause for 5s
        time.sleep(8)
    else:
        # Use threading backend for GPU compatibility if GPU acceleration is available
        daskclient = None
        joblib_backend = 'threading' if gpu.GPU_AVAILABLE else 'loky'

    print(f"Using Joblib backend: {joblib_backend}")
    # Sort matchT in descending order (so t-delay increases over rows)
    matchT = sorted(matchT, reverse=True)
    '''#OLD MULTIPROCESSING CODE
    # Prepare worker function with partial to pass all the constant parameters
    worker_func = functools.partial(
        process_time_value, files=files,
        T=T,
        crossfiles=crossfiles, L=L,
        col_labels=col_labels, ext=ext,
        exclude_col_labels=exclude_col_labels,
        calc_AMI=calc_AMI, calc_MI=calc_MI, calc_Morans=calc_Morans,
        bins=bins, verbose=verbose
    )
    '''
    # Execute the parallel processing
    if ncores > 1:
        ''' OLD MULTIPROCESSING CODE
        with Pool(processes=ncores) as pool:
        #with managed_pool(pool_processes=ncores) as pool:
            results = pool.map(worker_func, matchT)
        '''
        #with joblib.Parallel(n_jobs=ncores, backend=joblib_backend) as parallel:
        #    results = parallel(joblib.delayed(worker_func)(tval) for tval in matchT)
        # End of parallel processing

        # Execute parallel processing with specified backend
        with parallel_backend(joblib_backend, n_jobs=ncores):
            results = Parallel(
                batch_size="auto",
                verbose=10 if True else 0 #verbose else 0  # Joblib verbosity
            )([
                    delayed(process_time_value_joblib_wrapper)(
                    tval, files, T, crossfiles, L, col_labels, ext,
                    exclude_col_labels, calc_AMI, calc_MI, calc_Morans, bins, verbose,
                    worker_id=i % ncores if ncores > 0 else 0
                    ) for i, tval in enumerate(matchT)
            ])
        # End of parallel processing

        #If Dask client is used, unregister the worker plugin
        if gpu.DASK_AVAILABLE and daskclient is not None:
            # Wait for all tasks to complete
            daskclient.run_on_workers(lambda: None)  # Synchronization point
            #daskclient.unregister_worker_plugin(gpu.setup_daskworker_gpu_context)
            daskclient.close()
            print("Dask client closed and worker plugin unregistered.")

        # Try removing pysal_cache if it exists
        if os.path.exists(pysal_cache_dir):
            try:
                shutil.rmtree(pysal_cache_dir)
                print(f"Removed pysal cache directory: {pysal_cache_dir}")
            except Exception as e:
                print(f"Error removing pysal cache directory {pysal_cache_dir}: {e}")

        # Collect and merge results from all processes
        for result in results:
            auto_NCC_dict, cross_NCC_dict, auto_ZNCC_dict, cross_ZNCC_dict, auto_AMI_dict, cross_AMI_dict, auto_MI_dict, cross_MI_dict, auto_MoransI_dict, cross_MoransI_dict = result
            
            # Merge dictionaries
            for key, col_data in auto_NCC_dict.items():
                for col, rep_data in col_data.items():
                    all_auto_NCC_dict[key][col].update(rep_data)
                    
            for key, col_data in cross_NCC_dict.items():
                for col1, col_data_2 in col_data.items():
                    for col2, rep_data in col_data_2.items():
                        all_cross_NCC_dict[key][col1][col2].update(rep_data)
                        
            for key, col_data in auto_ZNCC_dict.items():
                for col, rep_data in col_data.items():
                    all_auto_ZNCC_dict[key][col].update(rep_data)
                    
            for key, col_data in cross_ZNCC_dict.items():
                for col1, col_data_2 in col_data.items():
                    for col2, rep_data in col_data_2.items():
                        all_cross_ZNCC_dict[key][col1][col2].update(rep_data)
            
            if calc_AMI:
                for key, col_data in auto_AMI_dict.items():
                    for col, rep_data in col_data.items():
                        all_auto_AMI_dict[key][col].update(rep_data)
                        
                for key, col_data in cross_AMI_dict.items():
                    for col1, col_data_2 in col_data.items():
                        for col2, rep_data in col_data_2.items():
                            all_cross_AMI_dict[key][col1][col2].update(rep_data)
            
            if calc_MI:
                for key, col_data in auto_MI_dict.items():
                    for col, rep_data in col_data.items():
                        all_auto_MI_dict[key][col].update(rep_data)
                        
                for key, col_data in cross_MI_dict.items():
                    for col1, col_data_2 in col_data.items():
                        for col2, rep_data in col_data_2.items():
                            all_cross_MI_dict[key][col1][col2].update(rep_data)

            if calc_Morans:
                for key, col_data in auto_MoransI_dict.items():
                    for col, rep_data in col_data.items():
                        all_auto_MoransI_dict[key][col].update(rep_data)
                        
                for key, col_data in cross_MoransI_dict.items():
                    for col1, col_data_2 in col_data.items():
                        for col2, rep_data in col_data_2.items():
                            all_cross_MoransI_dict[key][col1][col2].update(rep_data)
    else:
        # Sequential processing if ncores is 1
        for tval in matchT:
            #result = worker_func(tval)
            result = process_time_value_joblib_wrapper(tval, files, T, crossfiles, L, col_labels, ext,
                    exclude_col_labels, calc_AMI, calc_MI, calc_Morans, bins, verbose)
            auto_NCC_dict, cross_NCC_dict, auto_ZNCC_dict, cross_ZNCC_dict, auto_AMI_dict, cross_AMI_dict, auto_MI_dict, cross_MI_dict, auto_MoransI_dict, cross_MoransI_dict = result
            
            # Merge dictionaries
            for key, col_data in auto_NCC_dict.items():
                for col, rep_data in col_data.items():
                    all_auto_NCC_dict[key][col].update(rep_data)
                    
            for key, col_data in cross_NCC_dict.items():
                for col1, col_data_2 in col_data.items():
                    for col2, rep_data in col_data_2.items():
                        all_cross_NCC_dict[key][col1][col2].update(rep_data)
                        
            for key, col_data in auto_ZNCC_dict.items():
                for col, rep_data in col_data.items():
                    all_auto_ZNCC_dict[key][col].update(rep_data)
                    
            for key, col_data in cross_ZNCC_dict.items():
                for col1, col_data_2 in col_data.items():
                    for col2, rep_data in col_data_2.items():
                        all_cross_ZNCC_dict[key][col1][col2].update(rep_data)
            
            if calc_AMI:
                for key, col_data in auto_AMI_dict.items():
                    for col, rep_data in col_data.items():
                        all_auto_AMI_dict[key][col].update(rep_data)
                        
                for key, col_data in cross_AMI_dict.items():
                    for col1, col_data_2 in col_data.items():
                        for col2, rep_data in col_data_2.items():
                            all_cross_AMI_dict[key][col1][col2].update(rep_data)
            
            if calc_MI:
                for key, col_data in auto_MI_dict.items():
                    for col, rep_data in col_data.items():
                        all_auto_MI_dict[key][col].update(rep_data)
                        
                for key, col_data in cross_MI_dict.items():
                    for col1, col_data_2 in col_data.items():
                        for col2, rep_data in col_data_2.items():
                            all_cross_MI_dict[key][col1][col2].update(rep_data)

            if calc_Morans:
                for key, col_data in auto_MoransI_dict.items():
                    for col, rep_data in col_data.items():
                        all_auto_MoransI_dict[key][col].update(rep_data)
                        
                for key, col_data in cross_MoransI_dict.items():
                    for col1, col_data_2 in col_data.items():
                        for col2, rep_data in col_data_2.items():
                            all_cross_MoransI_dict[key][col1][col2].update(rep_data)

            

    # Construct the dataframes from the merged dictionaries
    # Sort keys by t-delay to ensure rows are ordered by increasing t-delay
    sorted_keys_auto_NCC = sorted(all_auto_NCC_dict.keys(), key=lambda x: x[2])  # Sort by t-delay (index 2)
    
    # Construct dataframes for auto-NCC
    for key in sorted_keys_auto_NCC:
        col_data = all_auto_NCC_dict[key]
        row = {"t0": key[0], "t1": key[1], "t-delay": key[2], "Rmax": key[3]}
        for var, rep_var_data in col_data.items():
            NCC_rep_vals = [val for repkey, val in rep_var_data.items() if repkey.startswith("NCC{")]
            if len(NCC_rep_vals) > 0:
                row[f"AVG[NCC{{{var}}}]"] = np.mean(NCC_rep_vals)
            row.update(rep_var_data)
        df_auto_NCC = pan.concat([df_auto_NCC, pan.DataFrame([row])], ignore_index=True)

    
    # Construct dataframes for auto-ZNCC
    sorted_keys = sorted(all_auto_ZNCC_dict.keys(), key=lambda x: x[2])  # Sort by t-delay (index 2)
    for key in sorted_keys:
        if key not in all_auto_ZNCC_dict:
            continue
        col_data = all_auto_ZNCC_dict[key]
        row = {"t0": key[0], "t1": key[1], "t-delay": key[2], "Rmax": key[3]}
        for var, rep_var_data in col_data.items():
            ZNCC_rep_vals = [val for repkey, val in rep_var_data.items() if repkey.startswith("ZNCC{")]
            if len(ZNCC_rep_vals) > 0:
                row[f"AVG[ZNCC{{{var}}}]"] = np.mean(ZNCC_rep_vals)
            row.update(rep_var_data)
        df_auto_ZNCC = pan.concat([df_auto_ZNCC, pan.DataFrame([row])], ignore_index=True)

    # Construct dataframes for cross-NCC
    sorted_keys = sorted(all_cross_NCC_dict.keys(), key=lambda x: x[2])  # Sort by t-delay (index 2)
    for key in sorted_keys:
        if key not in all_cross_NCC_dict:
            continue
        col_data = all_cross_NCC_dict[key]
        row = {"t0": key[0], "t1": key[1], "t-delay": key[2], "Rmax": key[3]}
        for var1, col_data_2 in col_data.items():
            for var2, rep_var1var2_data in col_data_2.items():
                NCC_rep_vals = [val for repkey, val in rep_var1var2_data.items() if repkey.startswith("NCC{")]
                if len(NCC_rep_vals) > 0:
                    row[f"AVG[NCC{{{var1};{var2}}}]"] = np.mean(NCC_rep_vals)
                row.update(rep_var1var2_data)
        df_cross_NCC = pan.concat([df_cross_NCC, pan.DataFrame([row])], ignore_index=True)

    # Construct dataframes for cross-ZNCC
    sorted_keys = sorted(all_cross_ZNCC_dict.keys(), key=lambda x: x[2])  # Sort by t-delay (index 2)
    for key in sorted_keys:
        if key not in all_cross_ZNCC_dict:
            continue
        col_data = all_cross_ZNCC_dict[key]
        row = {"t0": key[0], "t1": key[1], "t-delay": key[2], "Rmax": key[3]}
        for var1, col_data_2 in col_data.items():
            for var2, rep_var1var2_data in col_data_2.items():
                ZNCC_rep_vals = [val for repkey, val in rep_var1var2_data.items() if repkey.startswith("ZNCC{")]
                if len(ZNCC_rep_vals) > 0:
                    row[f"AVG[ZNCC{{{var1};{var2}}}]"] = np.mean(ZNCC_rep_vals)
                row.update(rep_var1var2_data)
        df_cross_ZNCC = pan.concat([df_cross_ZNCC, pan.DataFrame([row])], ignore_index=True)

    # Construct dataframes for AMI and MI if calculated
    if calc_AMI and all_auto_AMI_dict is not None:
        df_auto_AMI = pan.DataFrame(columns=auto_AMI_cols)
        df_cross_AMI = pan.DataFrame(columns=cross_AMI_cols)
        
        
        # Construct dataframes for auto-AMI
        sorted_keys = sorted(all_auto_AMI_dict.keys(), key=lambda x: x[2])  # Sort by t-delay (index 2)
        for key in sorted_keys:
            col_data = all_auto_AMI_dict[key]
            row = {"t0": key[0], "t1": key[1], "t-delay": key[2], "Rmax": key[3]}
            for var, rep_var_data in col_data.items():
                AMI_rep_vals = [val for repkey, val in rep_var_data.items() if repkey.startswith("AMI{")]
                if len(AMI_rep_vals) > 0:
                    row[f"AVG[AMI{{{var}}}]"] = np.mean(AMI_rep_vals)
                row.update(rep_var_data)
            df_auto_AMI = pan.concat([df_auto_AMI, pan.DataFrame([row])], ignore_index=True)
        
        # Construct dataframes for cross-AMI
        sorted_keys = sorted(all_cross_AMI_dict.keys(), key=lambda x: x[2])  # Sort by t-delay (index 2)
        for key in sorted_keys:
            col_data = all_cross_AMI_dict[key]
            row = {"t0": key[0], "t1": key[1], "t-delay": key[2], "Rmax": key[3]}
            for var1, col_data_2 in col_data.items():
                for var2, rep_var1var2_data in col_data_2.items():
                    AMI_rep_vals = [val for repkey, val in rep_var1var2_data.items() if repkey.startswith("AMI{")]
                    if len(AMI_rep_vals) > 0:
                        row[f"AVG[AMI{{{var1};{var2}}}]"] = np.mean(AMI_rep_vals)
                    row.update(rep_var1var2_data)
            df_cross_AMI = pan.concat([df_cross_AMI, pan.DataFrame([row])], ignore_index=True)

    
    if calc_MI and all_auto_MI_dict is not None:
        df_auto_MI = pan.DataFrame(columns=auto_MI_cols)
        df_cross_MI = pan.DataFrame(columns=cross_MI_cols)
        
        # Construct dataframes for auto-MI
        sorted_keys = sorted(all_auto_MI_dict.keys(), key=lambda x: x[2])
        for key in sorted_keys:
            col_data = all_auto_MI_dict[key]
            row = {"t0": key[0], "t1": key[1], "t-delay": key[2], "Rmax": key[3]}
            for var, rep_var_data in col_data.items():
                MI_rep_vals = [val for repkey, val in rep_var_data.items() if repkey.startswith("MI{")]
                if len(MI_rep_vals) > 0:
                    row[f"AVG[MI{{{var}}}]"] = np.mean(MI_rep_vals)
                row.update(rep_var_data)
            df_auto_MI = pan.concat([df_auto_MI, pan.DataFrame([row])], ignore_index=True)

        # Construct dataframes for cross-MI
        sorted_keys = sorted(all_cross_MI_dict.keys(), key=lambda x: x[2])
        for key in sorted_keys:
            col_data = all_cross_MI_dict[key]
            row = {"t0": key[0], "t1": key[1], "t-delay": key[2], "Rmax": key[3]}
            for var1, col_data_2 in col_data.items():
                for var2, rep_var1var2_data in col_data_2.items():
                    MI_rep_vals = [val for repkey, val in rep_var1var2_data.items() if repkey.startswith("MI{")]
                    if len(MI_rep_vals) > 0:
                        row[f"AVG[MI{{{var1};{var2}}}]"] = np.mean(MI_rep_vals)
                    row.update(rep_var1var2_data)
            df_cross_MI = pan.concat([df_cross_MI, pan.DataFrame([row])], ignore_index=True)


    if calc_Morans and all_auto_MoransI_dict is not None:
        df_auto_MoransI = pan.DataFrame(columns=auto_MoransI_cols)
        df_cross_MoransI = pan.DataFrame(columns=cross_MoransI_cols)
        
        # Construct dataframes for auto-MoransI
        sorted_keys = sorted(all_auto_MoransI_dict.keys(), key=lambda x: x[2])
        for key in sorted_keys:
            col_data = all_auto_MoransI_dict[key]
            row = {"t0": key[0], "t1": key[1], "t-delay": key[2], "Rmax": key[3]}
            for var, rep_var_data in col_data.items():
                MoransI_rep_vals = [val for repkey, val in rep_var_data.items() if repkey.startswith("BVMoransI{")]
                if len(MoransI_rep_vals) > 0:
                    row[f"AVG[BVMoransI{{{var}}}]"] = np.mean(MoransI_rep_vals)
                row.update(rep_var_data)
            df_auto_MoransI = pan.concat([df_auto_MoransI, pan.DataFrame([row])], ignore_index=True)

        # Construct dataframes for cross-MoransI
        sorted_keys = sorted(all_cross_MoransI_dict.keys(), key=lambda x: x[2])
        for key in sorted_keys:
            col_data = all_cross_MoransI_dict[key]
            row = {"t0": key[0], "t1": key[1], "t-delay": key[2], "Rmax": key[3]}
            for var1, col_data_2 in col_data.items():
                for var2, rep_var1var2_data in col_data_2.items():
                    MoransI_rep_vals = [val for repkey, val in rep_var1var2_data.items() if repkey.startswith("BVMoransI{")]
                    if len(MoransI_rep_vals) > 0:
                        row[f"AVG[BVMoransI{{{var1};{var2}}}]"] = np.mean(MoransI_rep_vals)
                    row.update(rep_var1var2_data)
            df_cross_MoransI = pan.concat([df_cross_MoransI, pan.DataFrame([row])], ignore_index=True)
        


    # Drop all empty columns from the dataframes
    df_auto_NCC.dropna(axis=1, how='all', inplace=True)
    df_cross_NCC.dropna(axis=1, how='all', inplace=True)
    df_auto_ZNCC.dropna(axis=1, how='all', inplace=True)
    df_cross_ZNCC.dropna(axis=1, how='all', inplace=True)
    if df_auto_AMI is not None:
        df_auto_AMI.dropna(axis=1, how='all', inplace=True)
    if df_cross_AMI is not None:
        df_cross_AMI.dropna(axis=1, how='all', inplace=True)
    if df_auto_MI is not None:
        df_auto_MI.dropna(axis=1, how='all', inplace=True)
    if df_cross_MI is not None:
        df_cross_MI.dropna(axis=1, how='all', inplace=True)
    if df_auto_MoransI is not None:
        df_auto_MoransI.dropna(axis=1, how='all', inplace=True)
    if df_cross_MoransI is not None:
        df_cross_MoransI.dropna(axis=1, how='all', inplace=True)

    # Warn if duplicate columns exist in the dataframes
    if df_auto_NCC.columns.duplicated().any():
        print("WARNING!!!!! : Duplicate columns found in auto-NCC dataframe. These are: " + str(df_auto_NCC.columns[df_auto_NCC.columns.duplicated()]))
    if df_cross_NCC.columns.duplicated().any():
        print("WARNING!!!!!: Duplicate columns found in cross-NCC dataframe. These are: " + str(df_cross_NCC.columns[df_cross_NCC.columns.duplicated()]))
    if df_auto_ZNCC.columns.duplicated().any():
        print("WARNING!!!!!: Duplicate columns found in auto-ZNCC dataframe. These are: " + str(df_auto_ZNCC.columns[df_auto_ZNCC.columns.duplicated()]))
    if df_cross_ZNCC.columns.duplicated().any():
        print("WARNING!!!!!: Duplicate columns found in cross-ZNCC dataframe. These are: " + str(df_cross_ZNCC.columns[df_cross_ZNCC.columns.duplicated()]))
    if df_auto_AMI is not None and df_auto_AMI.columns.duplicated().any():
        print("WARNING!!!!!: Duplicate columns found in auto-AMI dataframe. These are: " + str(df_auto_AMI.columns[df_auto_AMI.columns.duplicated()]))
    if df_cross_AMI is not None and df_cross_AMI.columns.duplicated().any():
        print("WARNING!!!!!: Duplicate columns found in cross-AMI dataframe. These are: " + str(df_cross_AMI.columns[df_cross_AMI.columns.duplicated()]))
    if df_auto_MI is not None and df_auto_MI.columns.duplicated().any():
        print("WARNING!!!!!: Duplicate columns found in auto-MI dataframe. These are: " + str(df_auto_MI.columns[df_auto_MI.columns.duplicated()]))
    if df_cross_MI is not None and df_cross_MI.columns.duplicated().any():
        print("WARNING!!!!!: Duplicate columns found in cross-MI dataframe. These are: " + str(df_cross_MI.columns[df_cross_MI.columns.duplicated()]))
    if df_auto_MoransI is not None and df_auto_MoransI.columns.duplicated().any():
        print("WARNING!!!!!: Duplicate columns found in auto-MoransI dataframe. These are: " + str(df_auto_MoransI.columns[df_auto_MoransI.columns.duplicated()]))
    if df_cross_MoransI is not None and df_cross_MoransI.columns.duplicated().any():
        print("WARNING!!!!!: Duplicate columns found in cross-MoransI dataframe. These are: " + str(df_cross_MoransI.columns[df_cross_MoransI.columns.duplicated()]))
    

    return df_auto_NCC, df_cross_NCC, df_auto_ZNCC, df_cross_ZNCC, df_auto_AMI, df_cross_AMI, df_auto_MI, df_cross_MI, df_auto_MoransI, df_cross_MoransI
    
    


# =============================== FUNCTIONS SPECIFIC TO FILES PROCESSING ============================

def rename(dir, keyword, new_keyword):
    # Rename all files in dir (including the entire filepath) that contain keyword to new_keyword.
    # This function is used to rename files that conflict with other files in the directory.
    # This function is called by main() when a file with the same name already exists in the new directory.
    # The function will rename the file to avoid conflicts.
    # The function will rename the file by replacing keyword with new_keyword in the filename.
    # The function will rename the file by adding a number to the end of the filename if new_keyword is already present in the filename.

    files = glob.glob(dir + "\\*" + keyword + "*")
    for file in files:
        new_file = file.replace(keyword, new_keyword)
        if os.path.exists(new_file):
            new_file = new_file.replace(new_keyword, new_keyword + "_1")
        os.rename(file, new_file)


# =============================== FUNCTIONS SPECIFIC TO PRELIM FILES PROCESSING ============================

# Computes and returns ln-distributed points from t=10^0 to log10(t_max)
    # (the latter rounded down to 1 decimal place)
def logarithm10_time_bins(t_max, dt, size=0):
    
    
    if t_max <= 1:
        return [0.0]  # No points recorded if t_max < e^0

    tend_one = np.log10(t_max*1.01)
    power = np.arange(0.0, tend_one, 0.04)  # 25 time-stamps per decade

    t_fin = np.floor( np.power(10, power) / dt) * dt

    t_fin[t_fin <= 1.0] += dt  # Ensure t_fin is strictly above 1.0

    logspaced = np.unique(t_fin).tolist()

    # Removing time-stamps that are too close to each other (within dt/2.0)
    
    n = len(logspaced)
    for i in range(n - 1, 0, -1):
        if logspaced[i] < logspaced[i - 1] + dt / 2.0:
            logspaced.pop(i)

    
    # Add 0.0 to the list of time-stamps if it is not already present at the start
    if logspaced[0] > 0.0:
        logspaced.insert(0, 0.0)
    
    print(f"Generated {len(logspaced)} logarithmically spaced time-stamps from 0 to {tend_one} (log10({t_max})) with dt = {dt}.")
    #print(f"Time-stamps: {logspaced}")
    if size > 0 and len(logspaced) != size:
        print(f"Warning! The number of time-stamps generated {len(logspaced)} is not equal to the specified size = {size}. Returning the first  {size} time-stamps.")
        return logspaced[:size]
    else:
        return logspaced

# Computes and returns ln-distributed points from t=10^0 to log10(t_max)
# for small dt values (dt < 0.05). 
def logarithm10_time_bins_smalldt(t_max, dt, dt_factor =1, linear_windows =[(600, 1500, 50)], size=0):
    
    dt_effective = dt * dt_factor
    
    if t_max <= 1:
        return [0.0]  # No points recorded if t_max < e^0

    tend_one = np.log10(t_max*1.01)
    power = np.arange(0.0, tend_one, 0.04)  # 25 time-stamps per decade

    t_fin = np.floor( np.power(10, power) / dt_effective) * dt_effective

    t_fin[t_fin <= 1.0] += dt_effective  # Ensure t_fin is strictly above 1.0

    logspaced = np.unique(t_fin).tolist()

    # Removing time-stamps that are too close to each other (within dt/2.0)
    
    n = len(logspaced)
    for i in range(n - 1, 0, -1):
        if logspaced[i] < logspaced[i - 1] + dt / 2.0:
            logspaced.pop(i)

    
    # Add 0.0 to the list of time-stamps if it is not already present at the start
    if logspaced[0] > 0.0:
        logspaced.insert(0, 0.0)

    # For linear windows, first generate linearly spaced time-stamps, then remove the time-stamps that are already present in logspaced.
    # And finally, add the linearly spaced time-stamps to logspaced at the appropriate positions.
    for window in linear_windows:
        t_start, t_end, dt_lin = window
        lin_spaced = np.arange(t_start, t_end, dt_lin).tolist()
        # Remove time-stamps from logspaced between t_start and t_end.
        logspaced = [t for t in logspaced if t < t_start or t > t_end]
        # Add the linearly spaced time-stamps to logspaced at the appropriate positions.
        for t in lin_spaced:
            for i in range(len(logspaced) - 1):
                if logspaced[i] < t < logspaced[i + 1]:
                    logspaced.insert(i + 1, t)
                    break

    # Sort the time-stamps in logspaced to ensure they are in ascending order.
    logspaced.sort()

    
    print(f"Generated {len(logspaced)} logarithmically spaced time-stamps from 0 to {tend_one} (log10({t_max})) with dt = {dt}.")
    #print(f"Time-stamps: {logspaced}")
    if size > 0 and len(logspaced) != size:
        print(f"Warning! The number of time-stamps generated {len(logspaced)} is not equal to the specified size = {size}. Returning the first  {size} time-stamps.")
        return logspaced[:size]
    else:
        return logspaced    


''' # Summary of function get_1D_HarmonicFreq_Prelimsdata(...)
This function performs FFT-based harmonic frequency analysis on time series data, identifying spectral peaks
and computing averaged harmonic frequencies across experimental replicates.

## Processing Steps:
1. Validates X values are linearly spaced (required for FFT)
2. Performs FFT on each Y column not in exclude_Y_cols
3. Generates output dataframe with structure: k(X), FFT{Ycol1} ... FFT{YcolEnd}
4. If report_maxima=True:
   - Fits cubic splines to FFT data and finds local maxima using derivative analysis
   - Handles flat/noisy signals robustly with noise-aware detection
   - Constrains peak search to original frequency domain
   - For each base variable <vari>, averages peaks across all replicate columns <vari>_R_<digit>
   - Creates peak_df with columns: FFT-PEAK[{col}], FFT-PEAK-VALS[{col}], FFT-PEAK-REPAVG{<vari>}, FFT-PEAK-VALS-REPAVG{<vari>}

## Returns:
- df_fft: DataFrame with frequency domain data
- peak_df: DataFrame with peak analysis (if report_maxima=True), otherwise None

## Input Requirements:
- df: DataFrame with time series data
- X: Either a string (column name or index name) OR an array-like object (np.array, list) representing the independent variable
  * If string: Must exist as a column in df or as an index/MultiIndex level when X_as_index=True
  * If array: Must be linearly spaced, 1D, and have the same length as df rows
- exclude_Y_cols: List of column names to exclude from FFT analysis
- X_as_index: Boolean, set True if X refers to an index level rather than a column
- report_maxima: Boolean, whether to perform peak detection and harmonic averaging
- maxima_finder: Method for peak detection, either "cubic_spline" or "find_peaks" (scipy.signal)

## DataFrame Structure Requirements:
Input dataframe should have structure:
X, AVG[<var1>] ... AVG[<varN>], <var1>_R_{i}, ... <var1>_R_{M}, ..., <varN>_R_{k}, ... <varN>_R_{l}
where <var1>, <var2>, ... are the dependent variables, and <var1>_R_{i}, ... are the replicate columns.
NOTE: Replicate column names MUST be in the format <var>_R_<int> to be recognised as replicates by the function.

IMPORTANT: Error Conditions:
- X array length mismatch with DataFrame
- Non-linearly spaced X values
- Missing column/index references
- Insufficient data for spline fitting
'''
def get_1D_HarmonicFreq_Prelimsdata(df, pathtodir="", X="t", exclude_Y_cols=[], X_as_index= True, skip_row=0, 
            rel_threshold =0.05, abs_threshold=None, report_maxima=True, maxima_finder="cubic_spline"):

    df_timeseries = df.copy(); # Copy the dataframe to avoid modifying the original one
    df_fft = pan.DataFrame() # Create an empty dataframe to store the FFT results
    peak_df = pan.DataFrame() if report_maxima else None # Create an empty dataframe to store the peak results if report_maxima is True
    df_timeseries.columns = [col.strip() for col in df_timeseries.columns] # Remove leading and trailing whitespaces from column names
    print(f"gen_1D_HarmonicFreq_Prelimsdata: Processing {len(df_timeseries)} rows and {len(df_timeseries.columns)} columns.")
    # Check first if X is an array-like structure (np.array, list, df[col], etc.) or a string.
    if isinstance(X, (np.ndarray, list)):
        # If X is an array-like structure, check if it has the same length as the dataframe.
        Xvalues = np.array(X); Xvalues = Xvalues.flatten() # Flatten the array to ensure it is 1D
        if len(Xvalues) != len(df_timeseries):
            print(f"Error: Length of provided X-array ({len(X)}) does not match the length of the dataframe ({len(df_timeseries)}).")
            return None, None
        # If it matches, set Xvalues to X.
        df_timeseries["FFT_X"] = pan.Series(Xvalues)  # Add X as a column to the dataframe for FFT processing.
        X = "FFT_X"  # Set X to the new column name for further processing.
    elif isinstance(X, str):
        # If X_as_index, check if df is MultiIndex, and if so, set X as provided X label, else if single index, set X as the single index.
        if X in df_timeseries.columns:
            Xvalues = asarray(df_timeseries[X].values)
        elif X in df_timeseries.index.names:
            df_timeseries = df_timeseries.reset_index(level=X) #if isinstance(df_timeseries.index, pan.MultiIndex) else df_timeseries.reset_index()
            Xvalues = asarray(df_timeseries[X].values)
        else:
            print(f"Error: X={X} is not in the columns, index or MultiIndex levels. Please provide a valid X label or pass as np.linspace.")
            return None, None

    print(f"gen_1D_HarmonicFreq_Prelimsdata: X-values are {Xvalues} with length {len(Xvalues)}.")
    # Drop skipped rows if skip_row > 0
    if skip_row > 0:
        print(f"Skipping first {skip_row} rows from the dataframe...")
        df_timeseries = df_timeseries.iloc[skip_row:].reset_index(drop=True)
    # Next remove exclude_Y_cols from df_timeseries (and add X to the column if it is a string).
    Y_cols = [col for col in df_timeseries.columns if col not in exclude_Y_cols and col != X and col not in df_timeseries.index.names]
    if isinstance(X, str):
        df_timeseries = df_timeseries[[X] + Y_cols]
    else:
        df_timeseries = df_timeseries[Y_cols]  

    # Replace missing values with 0.0
    df_timeseries.fillna(0.0, inplace=True)
    print(f"Xvals-diff-std: {np.std(np.diff(Xvalues))}")
    # Also check if Xvalues is linearly spaced (i.e. std of the differences is < 1)
    if np.std(np.diff(Xvalues)) > 1:
        print(f"Error: Detected X-values are not linearly spaced. Please provide a valid X array reference.")
        #print(f"X-values: {Xvalues}")
        #print(f"X-diff: {np.diff(Xvalues)}")
        #print(f"X-diff-std: {np.std(np.diff(Xvalues))}")
        #return None, None
        # Try interpolating df_timeseries to a linear space, using gen_interpolated_df()
        print("WARNING: Interpolating time values in file to generate a standardised time range.")
        # First generate a linear space of X values using max and min values of df_timeseries[X]
        Xvalues = np.linspace(np.min(Xvalues), np.max(Xvalues), len(df_timeseries))
        #print(f"Xvals: {Xvalues} with shape {Xvalues.shape} and length {len(Xvalues)}.")
        print(isinstance(Xvalues, np.ndarray))
        print(isinstance(gpu.to_cpu(Xvalues), _np_base.ndarray))
        # Save Xvalues as column "X" in a new temp dataframe called df_temp
        df_temp = pan.DataFrame({X: Xvalues})
        df_timeseries = gen_interpolated_df(df_temp, df_timeseries, X)
        if df is None:
            print(f"Error: Interpolation failed at {pathtodir}, with compare_label={X}")
            return None, None
        Xvalues = asarray(df_timeseries[X].values)  # Update Xvalues to the new interpolated values
        print(f"Interpolation successful. New Xvalues length: {len(Xvalues)}, with STD of differences: {np.std(np.diff(Xvalues))}")


    #df_fft[f"k({X})"] = np.fft.fftfreq(len(Xvalues), d=(Xvalues[1] - Xvalues[0]))[:len(Xvalues)//2]
    try:
        df_fft[f"k({X})"] = pan.Series(fft.fftfreq(len(Xvalues), d=(Xvalues[1] - Xvalues[0]))[:len(Xvalues)//2])
    except Exception as e:
        print(f"Error: Failed to compute FFT frequencies for X={X}. Exception: {e}")
        return None, None
    # Get the positive frequencies corresponding to the FFT result

    # Perform FFT on each Y column, storing the results in df_fft
    for col in Y_cols:
        # Perform FFT on the Y column
        
        #fft_result = np.fft.fft(df_timeseries[col].values)[:len(Xvalues)//2]
        fft_result = fft.fft(asarray(df_timeseries[col].values))[:len(Xvalues)//2]
        print(f"Performed FFT on column '{col}' with {len(fft_result)} frequency components.")
        # Store the positive frequencies in df_fft
        df_fft[f"FFT{{{col}}}"] = pan.Series(np.abs(fft_result))  # Store the absolute values of the FFT result
        if report_maxima:
            # Find the local maxima corresponding to first 5 peaks (if any) in the FFT result
            harm_peaks = np.nan*np.zeros(10); harm_peaks_fftvals = np.nan*np.zeros(10)
            # Stores the harmonic freq peaks and their corresponding FFT values 
            # (upto 10 peaks in descending order of magnitude)

            # First check if the FFT result is flat (or nearly flat) by seeing if fft_result[1:] is all zeros
            # Recall FFT of a flat signal is a Dirac delta-function at 0 frequency.
            df_fft_DC_max = df_fft.loc[1:, f"FFT{{{col}}}"].max()
            df_fft_DC_STD = df_fft.loc[1:, f"FFT{{{col}}}"].std(); df_fft_DC_mean = df_fft.loc[1:, f"FFT{{{col}}}"].mean()
            is_flat= ( df_fft_DC_max < 1e-6 ) or\
                ( df_fft_DC_STD < 1e-6 and df_fft_DC_mean < 1e-6 ) or np.all(np.isclose(fft_result[1:], 0.0))
                
            if is_flat:
                print("Warning: Flat or nearly flat column seen in:" + col 
                      + f" with DC Mode {df_fft_DC_max} and Mean {df_fft_DC_mean} ... Adjusting Peak to 0 ...")
                harm_peaks[0] = 0.0; harm_peaks_fftvals[0] = np.abs(fft_result[0])
                # Recall FFT of a flat signal is a Dirac delta-function at 0 frequency.
                peak_df[f"FFT-PEAK[{col}]"] = pan.Series(harm_peaks);
                peak_df[f"FFT-PEAK-VALS[{col}]"] = pan.Series(harm_peaks_fftvals)
                continue
            
            # Set a 5% threshold for peak detection to avoid noise
            #max_signal = np.max(df_fft[f"FFT{{{col}}}"]); thresh = 0.05*max_signal
            max_signal = np.max(df_fft.loc[1:, f"FFT{{{col}}}"]); 
            thresh = rel_threshold*max_signal if abs_threshold is None else max(rel_threshold*max_signal, abs_threshold)

            # Next fit a cubic spline to the FFT result.
            if maxima_finder == "cubic_spline":
                cs = interpolate.InterpolatedUnivariateSpline(df_fft[f"k({X})"], df_fft[f"FFT{{{col}}}"], k=3)
                # Find the local minima.
                cs_1st_derivative = cs.derivative()
                cs_2nd_derivative = cs.derivative(n=2)
                #quadratic_spline_roots() finds the roots of the quadratic spline function 
                # (1st derivative of cubic spline cs) in each interval [a, b] of the spline.
                extrema = quadratic_spline_roots(cs_1st_derivative)
                local_maxima = sorted([x for x in extrema if cs_2nd_derivative(x) < 0], key=lambda x: cs(x), reverse=True)
                # All local maxima sorting in descending order of their FFT values (given by negative 2nd derivative).

                # Next remove all local maxima that are below the threshold or not in the range of the FFT result
                local_maxima = [x for x in local_maxima if cs(x) > thresh and x >= df_fft[f"k({X})"].min() and x <= df_fft[f"k({X})"].max()]
                maxima_data = [(freq, cs(freq)) for freq in local_maxima]

            elif maxima_finder == "find_peaks":
                # Use scipy's find_peaks to find the local maxima
                peaks, _ = signal.find_peaks(df_fft.loc[1:, f"FFT{{{col}}}"], height=thresh); peaks +=1;
                peak_frequencies = asarray(df_fft[f"k({X})"].iloc[asnumpy(peaks)].values)
                peak_fft_values = asarray(df_fft[f"FFT{{{col}}}"].iloc[asnumpy(peaks)].values)
                # Create (frequency, value) pairs and sort
                maxima_data = sorted(zip(peak_frequencies, peak_fft_values), 
                                    key=lambda pair: pair[1], reverse=True)
                # Filter maxima_data to remove peaks with frequencies less than the second frequency (df_fft[f"k({X})"].iloc[1])
                #maxima_data = [pair for pair in maxima_data if pair[0] >= df_fft[f"k({X})"].iloc[1]]
            else:
                print(f"Error: Invalid maxima_finder '{maxima_finder}'. Use 'cubic_spline' or 'find_peaks'.")
                return None, None
            # Assign these values to harm_peaks and harm_peaks_fftvals
            for i, (peakfreq, peakval) in enumerate(maxima_data):
                if i >= len(harm_peaks):
                    break
                if peakfreq <= df_fft[f"k({X})"].iloc[1]:
                    continue # Skip frequencies less than the second frequency
                harm_peaks[i] = peakfreq
                harm_peaks_fftvals[i] = peakval
                

            peak_df[f"FFT-PEAK[{col}]"] = pan.Series(harm_peaks);
            peak_df[f"FFT-PEAK-VALS[{col}]"] = pan.Series(harm_peaks_fftvals)
    
    # Now calculate the mean harmonic peak frequencies for each {vari} by averaging over the local maxima for each <vari>_R_{*} columns
    if report_maxima:
        # Find all columns that match the pattern <basevar>_R_<digit>
        replicate_pattern = re.compile(r'^(.+)_R_(\d+)$')
        base_varcols =[ replicate_pattern.match(col).group(1) for col in Y_cols if replicate_pattern.match(col) ]
        # Get the unique base variable names
        base_varcols = list(set(base_varcols))
        for base_varcol in base_varcols:
                # Get all columns that start with the base variable name
                rep_cols = [c for c in Y_cols if c.startswith(base_varcol)]
                # Calculate the mean harmonic peak frequencies by averaging over all rep_cols
                peak_df[f"FFT-PEAK-REPAVG{{{base_varcol}}}"] = peak_df[[f"FFT-PEAK[{c}]" for c in rep_cols]].mean(axis=1)
                peak_df[f"FFT-PEAK-VALS-REPAVG{{{base_varcol}}}"] = peak_df[[f"FFT-PEAK-VALS[{c}]" for c in rep_cols]].mean(axis=1)
                # Also calculate the variance of the harmonic peak frequencies by averaging over all rep_cols
                peak_df[f"FFT-PEAK-REPVAR{{{base_varcol}}}"] = peak_df[[f"FFT-PEAK[{c}]" for c in rep_cols]].var(axis=1)
                peak_df[f"FFT-PEAK-VALS-REPVAR{{{base_varcol}}}"] = peak_df[[f"FFT-PEAK-VALS[{c}]" for c in rep_cols]].var(axis=1)
        
        return df_fft, peak_df
    
    return df_fft, None

''' # Summary of the function gen_1D_HarmonicFreq_PrelimsWrapper(...)
This function is a wrapper for gen_1D_HarmonicFreq_Prelimsdata that processes multiple files in a directory.
It reads each file, extracts the time series data, and applies the harmonic frequency analysis.
index_cols is a list of columns to be used as index for the dataframe (either numeric values indicating iloc index or strings indicating column names).
If index_cols has more than one element, the dataframes will be read as MultiIndex dataframes.
'''
def gen_1D_HarmonicFreq_PrelimsWrapper(files, pathtodir, ext="csv", X="t", exclude_Y_cols=[], index_cols = [], X_as_index=True, skip_row=0,
            rel_threshold =0.05, abs_threshold=None, report_maxima=True, maxima_finder="cubic_spline"):
    for file in files:
        try:
            if ext == "csv":
                df = pan.read_csv(file, header=0, index_col=index_cols if len(index_cols) > 0 else None)
            else:
                df = pan.read_table(file, header=0, index_col=index_cols if len(index_cols) > 0 else None)
        except Exception as e:
            print(f"Error: Failed to read file {file} with error message: {e}")
            continue
        
        # Modify column names in df to remove leading and trailing whitespaces.
        df.columns = [col.strip() for col in df.columns]
        
        # Call gen_1D_HarmonicFreq_Prelimsdata on the dataframe
        df_fft, peak_df = get_1D_HarmonicFreq_Prelimsdata(df, pathtodir, X=X, exclude_Y_cols=exclude_Y_cols, X_as_index=X_as_index,
                                                            skip_row=skip_row, rel_threshold=rel_threshold, abs_threshold=abs_threshold,
                                                            report_maxima=report_maxima, maxima_finder=maxima_finder)
        if df_fft is None:
            print(f"Error: Failed to process file {file}. Skipping this file.")
            continue
        
        savedir = pathtodir if pathtodir else os.path.dirname(file)
        savedir += "/FFT/" if not savedir.endswith("/") else "FFT/"
        # Create the directory if it does not exist
        Path(savedir).mkdir(parents=True, exist_ok=True)
        output_file = os.path.join(savedir, f"FFTSig_{os.path.basename(file)}")
        try:
            if ext == "csv":
                df_fft.to_csv(output_file, index=False, sep=",", header=True)
            else:
                df_fft.to_csv(output_file, index=False, header=True)
        except Exception as e:
            print(f"Error: Failed to save FFT data for file {file} with error message: {e}")
            continue
        
        if report_maxima and peak_df is not None:
            try:
                peak_output_file = os.path.join(savedir, f"HarmonicPeaks_{os.path.basename(file)}")
                if ext == "csv":
                    peak_df.to_csv(peak_output_file, index=False, sep=",", header=True)
                else:
                    peak_df.to_table(peak_output_file, index=False, header=True)
                print(f"Saved harmonic peaks data for file {file} to {peak_output_file}.")
            except Exception as e:
                print(f"Error: Failed to save harmonic peaks data for file {file} with error message: {e}")
                continue


'''# Summary of the function gen_missing_zero_Prelims(...)
# This function combs through the files in files, which are in the directory pathtodir.
These files follow the naming convention "/TSERIES*_T_{Tval}_a_{aval}_R_{Rval}.csv" where Tval, aval and Rval are the values of T, a and R respectively.
The list of Tvals is given by T_vals. This function will do the following:
1. It will go through all the files in files in descending order of T_vals (starting from the second largest Tval).
2. If at least the last two values in the columns corresponding to compare_col_labels are zero, it will "extend" the matched files over time (for all larger Tvals (say Tval_larger1, Tval_larger2, ...))
by creating new files with the same {aval}, but with {Tval_largeri} replaced by the larger Tval {Tval}, and {Rval} replaced by Rmax +1, where Rmax is the largest Rval for the files matching the given {Tval_largeri}
This "extension" will occur by copying the matched files, and appending 0 values to all columns in the copied files for the new time-stamps.
3. If the last two values in the columns corresponding to compare_col_labels are not zero, it will skip the file.
4. For the "extended" files, the "t" column (time-stamps) will be generated using logarithmically spaced time-stamps from 0 to {Tval_largeri}, with time-stamps spaced by dt.

Finally, if the files have non-standardised columns (a label in include_col_labels is not present in the column names of the files), the function will return an error.
'''
def gen_missing_zero_Prelims(files, pathtodir, T_vals, ext="csv", dt =None, compare_col_labels= ["<<P(x; t)>_x>_r" ] ):
    
    T_vals_desc = np.sort(T_vals)[::-1]
    for Tval in T_vals_desc[1:]:
        # Going from second largest Tval to the smallest Tval.
        # Find all files in files that have Tval.
        t = int(Tval) if float(Tval).is_integer() else Tval
        selectedfiles = glob.glob(pathtodir + "/TSERIES*_T_" + str(t) + "*.csv")
        maxR = len(selectedfiles)
        for file in selectedfiles:
            try:
                Rval = int(re.search(r'R_[\d]+', file).group(0).split("_")[1])
                # Read the file.
                if(ext == "csv"):
                    df = pan.read_csv(file, header=0)
                else:
                    try:
                        df = pan.read_table(file, header=0)
                    except Exception as e:
                        print("Error: Non-standard extension for file " + pathtodir +"/" + file + " with error message: \n" + str(e))
                        return None
                # Modify column names in df to remove leading and trailing whitespaces.
                df.columns = [col.strip() for col in df.columns]
                # First check if the column names in compare_col_labels are present in df.columns. If not, return an error.
                if not all([col in df.columns for col in compare_col_labels]):
                    print("Error: Columns in file " + pathtodir + "/" + file + " are not standardised.")
                    return None
                # Check if the last two values in the columns corresponding to compare_col_labels are zero.
                for col in compare_col_labels:
                    if not all(df[col].iloc[-2:] == 0):
                        raise SkipFile
                print(f"Extending file {file} at T = {Tval} over time.")
                # If the last two values in the columns corresponding to compare_col_labels are zero, "extend" the matched files over time.
                # This is done by creating new files with the same {aval}, but with {Tval_larger} replaced by the larger Tval {Tval}, and {Rval} replaced by Rmax +1.
                # The "extension" will occur by copying the matched files, and appending 0 values to all columns in the copied files for the new time-stamps.
                # The "extended" files will have the "t" column (time-stamps) generated using logarithmically spaced time-stamps from 0 to {Tval_larger}, with time-stamps spaced by dt.
                for Tval_larger in T_vals_desc:
                    if Tval_larger <= Tval:
                        continue
                    # Find all files in files that have Tval_larger.
                    t_larger = int(Tval_larger) if float(Tval_larger).is_integer() else Tval_larger
                    selectedfiles_larger = glob.glob(pathtodir + "/TSERIES*_T_" + str(t_larger) + "*.csv")
                    if len(selectedfiles_larger) != 0:
                        df_larger = pan.read_csv(selectedfiles_larger[0], header=0)
                    Rmax = len(selectedfiles_larger)
                    # Find the largest Rval for the files matching the given {Tval_larger}.
                    Rmax = max([int(re.search(r'R_[\d]+', file).group(0).split("_")[1]) for file in selectedfiles_larger])
                    # Create a new file with the same {aval}, but with {Tval_larger} replaced by the larger Tval {Tval}, and {Rval} replaced by Rmax +1.
                    new_file = file.replace("T_" + str(t), "T_" + str(t_larger)).replace("R_" + str(Rval), "R_" + str(Rmax + 1))
                    new_df = df.copy()
                    # Generate the "t" column (time-stamps) using logarithmically spaced time-stamps from 0 to {Tval_larger}, with time-stamps spaced by dt.
                    if not all(df["t"] == 0) or dt == None:
                        new_df["t"] = logarithm10_time_bins(Tval_larger, dt, len(df))
                    else:
                        new_df["t"] = logarithm10_time_bins(Tval_larger, dt)
                    # Append 0 values to all columns in the copied files for the new time-stamps.
                    for col in new_df.columns:
                        if col != "t":
                            new_df[col] = new_df[col].append(pan.Series([0]*(len(new_df) - len(df))), ignore_index=True)
                    # Write the new file to the directory.
                    new_df.to_csv(new_file, index=False)
                    print(f"Extended file {file} at T = {Tval} over time to {new_file} at T = {Tval_larger} and R = {Rmax + 1}.")
                    # Generate the "t" column (time-stamps) using logarithmically spaced time-stamps from 0 to {Tval_larger}, with time-stamps spaced by dt.
            except SkipFile:
                print(f"Skipping file {file} at T = {Tval} as the last two values in the columns corresponding to compare_col_labels are not zero.")
                continue
            except Exception as e:
                print(f"Error: {e} for file {file} at T = {Tval}. Skipping file.")
                
    return 1

'''# Summary of the function gen_MEAN_INDVL_Colsfiledata(...)
# This function generates the mean for each column for each file in files, as well as the mean across all files for each column.
# It will only include include_col_labels from the mean calculations. If the files have non-standardised columns 
# (a label in include_col_labels is not present in the column names of the files), the function will return an error.
# It will create a dataframe with the following column structure:
#R_max  t   AVG[{var}]_SURV   ...   SD[{var}]_SURV  ... R_SURV[{var}]   ...   AVG[{var}]_ALL   ...   SD[{var}]_ALL
# where {var} is one of the species names in include_col_labels (excluding "t").
# The optional argument dt is used to generate a logarithmically spaced time range from 0 to t_max, with time-stamps spaced by dt,
# if the time values in the files are all zero.
# The optional argument handle_nonstandard_time is used to handle non-standard time value conflicts in the files.
# It has THREE VALID OPTIONS:
# 1. If handle_nonstandard_time = "error", the function will return an error if the time values in the files are not standardised.
# 2. If handle_nonstandard_time = "skip", the function will skip files with non-standardised time values.
# 3. If handle_nonstandard_time = "interpolate", the function will interpolate (using cubic splines) the time values in the files to generate a standardised time range.
This will be done by calling the function gen_INTERPOLATED_filedata(...) with the files and include_col_labels.
'''
def gen_MEAN_INDVL_Prelimsfiledata(files, pathtodir="", ext="csv", tmax =None, dt =None, include_col_labels= ["t",  "<<P(x; t)>_x>_r" , 
                                        "<<G(x; t)>_x>_r", "<<Pr(x; t)>_x>_r", "<<W(x; t)>_x>_r", "<<O(x; t)>_x>_r"], handle_nonstandardtime_fileconflicts = "error"):
    
    # Read the first file in files.
    if(ext == "csv"):
        df = pan.read_csv(files[0], header=0)
    else:
        try:
            df = pan.read_table(files[0], header=0)
        except Exception as e:
            print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
            return None
    
    # Store column names in col_labels. Strip elements of any leading or trailing whitespaces.
    df = df.reset_index(drop=True)
    df.columns = [col.strip() for col in df.columns]
    col_labels = [col for col in df.columns if col in include_col_labels and col != "t"];

    # Create a df with col_labels as column names. This will store the mean data, fot each column in col_labels for each file in files.
    df_file = pan.DataFrame(columns= ["Rmax", "t"] + ["AVG[" + x + "]_ALL" for x in col_labels] + ["VAR[" + x + "]_ALL" for x in col_labels] 
                            + ["AVG[" + x + "]_SURV" for x in col_labels] + ["VAR[" + x + "]_SURV" for x in col_labels] + ["R_SURV[" + x + "]" for x in col_labels])
    
    pooled_df = pan.DataFrame()

    if not all(df["t"] == 0) or dt == None or tmax == None:
        pooled_df["t"] = df["t"] # Store the time values in pooled_df.
    else:
        # Use a custom time range if the time values are all zero.
        if dt > 0.05:
            pooled_df["t"] = logarithm10_time_bins(tmax, dt, len(df))
        else:
            pooled_df["t"] = logarithm10_time_bins_smalldt(tmax, dt, dt_factor =1, linear_windows =[(600, 1500, 50)], size= len(df))


    maxR = len(files)
    # Iterate over all files in files. For each file, concatenate the columns in col_labels across all files.
    # Store the concatenated data in concat_df.
    for file in files:

        #First read the file.
        if(ext == "csv"):
            try:
                df = pan.read_csv(file, header=0)
            except Exception as e:
                print("Error: Failed to read file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
                return None
        else:
            try:
                df = pan.read_table(file, header=0)
            except Exception as e:
                print("Error: Non-standard extension for file " + pathtodir +"/" + files[0] + " with error message: \n" + str(e))
                return None
        # Modify column names in df to remove leading and trailing whitespaces.
        df = df.reset_index(drop=True)
        df.columns = [col.strip() for col in df.columns]
        # Remove exclude_col_labels from df.
        df = df[[col for col in df.columns if col in include_col_labels]]

        # First check if the columns in df are the same as col_labels. If not, return an error.
        if not all([col in df.columns for col in col_labels]):
            print("Error: Columns in file " +  file + " are not standardised.")
            return None
        
        # Also check if the time values in the files are the same. If not, return an error.
        #if not all(df["t"] == pooled_df["t"]) and not all(df["t"] == 0):
        if not (df["t"].reset_index(drop=True).equals(pooled_df["t"].reset_index(drop=True))) and not all(df["t"] == 0):
            print("Error: Time values in file " +  file + " are not standardised.")
            if handle_nonstandardtime_fileconflicts == "error":
                return None
            elif handle_nonstandardtime_fileconflicts == "skip":
                print("Skipping file " + pathtodir + "/" + file + " due to non-standardised time values.")
                continue
            elif handle_nonstandardtime_fileconflicts == "interpolate":
                print("WARNING: Interpolating time values in file to generate a standardised time range.")
                df = gen_interpolated_df(pooled_df, df, "t")
                if df is None:
                    print("Error: Interpolation failed for file " +  file + ". Skipping file.")
                    continue
                else:
                    print("Interpolation successful for file.")

                    #First set all values in df < abs(1e-6) to 0 (to avoid negative values in the interpolated data).
                    df = df.applymap(lambda x: 0 if x < 1e-6 and x > -1e-6 else x)
                    # Save the interpolated file to disk as pathtodir + "/INTERPOLATED_" + file.
                    # First get the filename by removing the pathtodir from file.
                    filename = file.split(pathtodir + "/")[1]
                    df.to_csv(pathtodir + "/INTERPOLATED_" + filename, index=False)
            else:
                print("Error: Invalid option for handle_nonstandard_time. Please choose from 'error', 'skip' or 'interpolate'.")
                return None
            
        
        Rstr = re.findall(r'R_[\d]+', file)[0]
        
        for col in col_labels:

            # Store column values in pooled_df
            pooled_df[col + "_" + Rstr] = df[col]

    # Done with all files in files.
    # Now generate element-wise mean and standard deviation for each column in col_labels, across all files.
    # Store the mean and standard deviation in df_file.

    for col in col_labels:

        # Find mean and standard deviation of each column in pooled_df.
        mean_col = pooled_df.filter(like= col).mean(axis=1)
        variance_col = pooled_df.filter(like= col).var(axis=1)
        # Store the mean and standard deviation in df_file.
        df_file["AVG[" + col + "]_ALL"] = mean_col
        df_file["VAR[" + col + "]_ALL"] = variance_col
        
        # Surviving data is calculated as the element-wise mean and standard deviation of the column across all files for each non-zero element.
        # Store the mean and standard deviation in df_file.
        mean_col_surv = pooled_df.filter(like= col).replace(0, np.nan).mean(axis=1)
        variance_col_surv = pooled_df.filter(like= col).replace(0, np.nan).var(axis=1)
        df_file["AVG[" + col + "]_SURV"] = mean_col_surv
        df_file["VAR[" + col + "]_SURV"] = variance_col_surv

        # Store the number of non-zero elements in the column across all files in df_file.
        R_surv_col = pooled_df.filter(like= col).replace(0, np.nan).count(axis=1)
        df_file["R_SURV[" + col + "]"] = R_surv_col
        # Store the Rmax in df_file, as a column with the same value for all rows.
        df_file["Rmax"] = [maxR for i in range(len(df_file))]
        df_file["t"] = pooled_df["t"]
    
    return df_file
        
