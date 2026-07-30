'''
# This Python script accepts a CSV file as input, normalises the data b/w 0 and 1 (grayscale), 
# and then clusters the data using KMeans clustering. The original image is optionally displayed 
#in a subplot next to the clustered image and is saved as a PNG file. 
#Finally the clustered matrix is saved as a csv file to the output directory.

'''

import pandas as pan
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import cv2
import glob
import seaborn as sea
import regex as re
import secrets
from pathlib import Path
from scipy.ndimage import gaussian_filter
from scipy.ndimage import label

from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler

rng = np.random.default_rng(secrets.randbits(128))

def rgb_to_hex(r, g, b):
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)

def hex_to_rgb(hex):
    hex = hex.lstrip('#')
    return tuple(int(hex[i:i+2], 16) for i in (0, 2, 4))

def grayscale_to_hex(value):
    value = int(value * 255)
    return rgb_to_hex(value, value, value)

def normalise_data(data):
    # Normalise the data to be between 0 and 1
    scaler = MinMaxScaler()
    data = scaler.fit_transform(data)
    return data

def cluster_data(data, n_clusters, s, rep_no):
    # Cluster the data using KMeans
    kmeans_model = KMeans(n_clusters=n_clusters, init= 'k-means++' , 
                          n_init= 1, max_iter= 500, random_state= rng.integers(0, 1000)).fit(data)
    return kmeans_model#.labels_

def plot_clustered_data(data_df, L, n_clusters, output_dir, filename, **kwargs):
    # Plot the clustered data
    fig, ax = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    # Extract the species number and other parameters from kwargs
    s = kwargs.get("Sp_No", "NA")
    R_val = kwargs.get("R", "NA"); a_val = kwargs.get("a", "NA"); dP_val = kwargs.get("dP", "NA"); t_val = kwargs.get("t", "NA")

    # Plot the original data
    pos = ax[0].imshow(data_df[f'Original Value_S_{s}_R_{R_val}'].to_numpy().reshape(L, L), cmap='gray', aspect='auto')
    ax[0].set_title('Original Data')

    # Set colourbar for the original data
    cbar1 = fig.colorbar(pos, ax=ax[0]); cbar1.minorticks_on()

    
    try:
        # First try to get the hex_list of the cluster centers from kwargs (kwargs["hex_list"] if present)
        hex_list = kwargs["hex_list"]
    except KeyError:
        # If not present, use default grayscale values
        f_list =np.linspace(0, 1, num=n_clusters)
        hex_list = [grayscale_to_hex(value) for value in f_list] 

    #print(hex_list)
    # Get cmap from the hex_list
    # Convert hexadecimal colors to RGB tuples
    colors_rgb = [tuple(int(hex.lstrip('#')[i:i+2], 16) / 255.0 for i in (0, 2, 4)) for hex in hex_list]

    # Create a custom colormap using the specified colors
    cmap_custom = LinearSegmentedColormap.from_list('custom_colormap', colors_rgb, N=len(hex_list))

        
    # Plot the clustered data
    pos = ax[1].imshow(data_df[f'Ordn_cluster_S_{s}_R_{R_val}'].to_numpy().reshape(L, L), cmap= cmap_custom, aspect='auto')
    ax[1].set_title('Clustered Data')
    cbar2 = fig.colorbar(pos, ax=ax[1]); cbar2.minorticks_on()

    #Unwrap the rest of the kwargs and use it as the title of the overall plot
    title = r"$\rho_%s(t= %s)$, $R = $" %(s, t_val) + str(a_val) + r", $dP = $" + str(dP_val) + r", $t = $" + str(t_val) + r", $n = $" + str(R_val) 
    plt.suptitle(title)
    output_path = os.path.join(output_dir, f'Rho{s}_{filename}.png')
    plt.savefig(output_path)
    if( a_val % 0.003 == 0.0):
        plt.show()
    plt.close()

def save_clustered_data(clustered_datadf, output_dir, filename, S_list=[1,2], Rep_List= [0]):
    # Save the clustered data to a CSV file
    output_path = os.path.join(output_dir, f'{filename}.csv')
    
    # Extract the cluster labels for each species and save them to the CSV file
    for s in S_list:
        for rep_no in Rep_List:
            clustered_datadf[f'n_cluster_S_{s}_R_{rep_no}'] = clustered_datadf[f'n_cluster_S_{s}_R_{rep_no}'].astype(int)
            clustered_datadf[f'Ordn_cluster_S_{s}_R_{rep_no}'] = clustered_datadf[f'Ordn_cluster_S_{s}_R_{rep_no}'].astype(int)
            clustered_datadf[f'cluster_center_S_{s}_R_{rep_no}'] = clustered_datadf[f'cluster_center_S_{s}_R_{rep_no}'].apply(lambda x: x[0])
    # Save only the columns that are required
    clipped_df = clustered_datadf[[col for col in clustered_datadf.columns if col.split('_')[0] not in ['hex']]]
    clipped_df.to_csv(output_path, index=False)

def generate_clusters(input_dir, output_png_dir, output_file_dir, n_clusters, L, dP, Geq, a_vals =[], S_list=[1,2], verbose=False):
    # Read the input CSV file
    # Check if kwargs["a"] is present
    # The input directory has a "a_vals.txt" file which stores a decimal number on each line. Extract these numbers as a list.

    # Check if a_vals is empty list.
    if len(a_vals) == 0:
        try:
            with open(os.path.join(input_dir, "a_vals.txt"), "r") as a_file:
                a_vals = [float(line.strip()) for line in a_file]
        except FileNotFoundError:
            print("a_vals.txt not found in the input directory. Terminating the program.")
            return
    print(f"List of a values: {a_vals}")
    
    for a in a_vals:
        parendir = input_dir + "L_" + str(L) + "_a_" + str(a) + "/dP_" + str(dP) + "/Geq_" + str(Geq) + "/"

        # Check if the parent directory exists
        if not os.path.exists(parendir):
            print(f"Parent directory {parendir} not found. Skipping...")
            continue
        
        # In parendir, we have files in sub-directories of the form T_{t}/FRAME_T_{t}_a_{a}_R_{R}.csv 
        #where t is a time value (can be a decimal) and R is a replicate number (integers >=0).
        # We first extract the t values from the subdirectories and sort them in ascending order.
        # Then we iterate over each replicate file for each t value and cluster the data (extracting the replicate number from the filename).

        # Extract the t values from the subdirectories
        t_vals = []
        t_vals = [re.findall(r'T_[\d]*[.][\d]+' , subdir)[0] if re.findall(r'T_[\d]*[.][\d]+' , subdir) 
                    else re.findall(r'T_[\d]+' , subdir)[0] for subdir in os.listdir(parendir) 
                    if os.path.isdir(os.path.join(parendir, subdir))]
        t_vals = sorted(list(set([float(re.findall(r'[\d]*[.][\d]+', t)[0]) 
                     if re.findall(r'[\d]*[.][\d]+', t) 
                     else int(re.findall(r'[\d]+', t)[0]) for t in t_vals])))
        # Only consider unique values of T.
        #t_vals = sorted(list(set(t_vals)))
        # Remove t_values that are less than 30000
        t_vals = [t for t in t_vals if t >= 30000]  
        print(t_vals)
        # Iterate over each t value
        for t in t_vals:
            # Iterate over each replicate file for each t value
            # Create empty dataframe to store the clustered data
            data_df = pan.DataFrame()
            max_rep = 0 #Stores maximum replicates encountered in each cell
            rep_list=[] #Stores list of replicates in each T directory.

            for file in glob.glob(parendir + "T_%g/*.csv" %(t)): # This is the replicate loop
                print(file)
                max_rep += 1
                for s in S_list: # This is the species loop
                    data = np.genfromtxt(file, delimiter=',', skip_header=1, usecols=s+2)
                    filename = Path(file).stem

                    #Extract the replicate number from the filename (R_{R}) (R is an integer >=0, last number in the filename)
                    rep_no = int(re.findall(r'\d+', filename)[-1])

                    if rep_no not in rep_list:
                        rep_list.append(rep_no)

                    # Normalise the data
                    data = normalise_data(data.reshape(-1,1))

                    #data = gaussian_filter(data.reshape(-1,L), sigma=0.75)

                    print(f"Processing at a = {a}, t = {t}, Species = {s}, Replicate = {rep_no}, dP = {dP}, Geq = {Geq},")

                    #Flattening the data to 1D array
                    data = data.flatten()
                    # Convert np array to pandas dataframe
                    #data_df = pan.DataFrame({f"Original Value_S_{s}_R_{rep_no}": data})
                    data_df[f"Original Value_S_{s}_R_{rep_no}"] = data

                    # Cluster the data
                    kmeans_model = cluster_data(data.reshape(-1,1), n_clusters, s, rep_no)

                    #Assigning the cluster labels to the data_df as a new column "n_cluster"
                    data_df[f'n_cluster_S_{s}_R_{rep_no}'] = kmeans_model.labels_
                    data_df[f'hex_code_S_{s}_R_{rep_no}'] = data_df[f"Original Value_S_{s}_R_{rep_no}"].apply(grayscale_to_hex)
                    data_df[f'cluster_center_S_{s}_R_{rep_no}'] = data_df[f'n_cluster_S_{s}_R_{rep_no}'].apply(lambda x: kmeans_model.cluster_centers_[x])

                    palette = kmeans_model.cluster_centers_.astype(float).flatten()
                    ord_palette = palette.copy();
                    ord_palette.sort()

                    #Order the cluster centers in ascending order in ord_palette and save indices of th

                    #Order the cluster labels in data_df by the order of ord_palette
                    data_df[f'Ordn_cluster_S_{s}_R_{rep_no}'] = data_df[f'n_cluster_S_{s}_R_{rep_no}'].apply(lambda x: np.where(ord_palette == palette[x])[0][0])
                    #print(data_df[f'n_cluster_S_{s}_R_{rep_no}'][:10])
                    #print(data_df[f'Ordn_cluster_S_{s}_R_{rep_no}'][:10])
                    if(verbose):
                        
                        print(f"Cluster centers for Sp {s}, Rep {rep_no}:\t"); print(palette)
                        print(f"Ordered Cluster centers for Sp {s}, Rep {rep_no}:\t"); print(ord_palette)
                        for each_cluster in range(0, n_clusters):
                            print('Cluster ', each_cluster)
                            #sea.palplot(sea.color_palette(list(data_df[data_df[f'n_cluster_S_{s}_R_{rep_no}'] == each_cluster][f'hex_code_S_{s}_R_{rep_no}'][:6])))
                            #plt.show(); plt.close()
                    
                    hex_palette = [grayscale_to_hex(value) for value in ord_palette]

                    print(hex_palette)
                
                    if(verbose):
                        # Plot the clustered data
                        estraz_args = {"Sp_No": s, "R": rep_no, "a": a, "dP": dP, "t": t, "Geq": Geq, "hex_list": hex_palette}
                        plot_clustered_data(data_df, L, n_clusters, output_png_dir, filename, **estraz_args)

                # End of s-loop
                if(verbose):
                    if(rep_no == 2):
                        print(f"At a = {a}, t = {t}, dP = {dP}, Geq = {Geq}")
                        print(data_df.head())
                        print(data_df.describe())
                        #print(data_df.info())
                        print(data_df.shape)
                        #print(data_df.columns)    
            # End of file-loop (replicate loop), save the clustered data to a CSV file
            rep_list.sort()
            filename = f"CLUSTERED_{n_clusters}_FRAME_T_{t}_a_{a}_dP_{dP}_Geq_{Geq}_R_{max_rep}"
            save_clustered_data(data_df, output_file_dir, filename, S_list, rep_list)
        # End of t-loop
    # End of a-loop

'''
Analyse the clusters generated by the generate_clusters function. This function reads the clustered data from a given input directory,
(which is the output_file_dir from the generate_clusters function), and generates a series of plots and statistics to analyse the clusters.

The function receives the following inputs:
1. input_dir: The input directory where the clustered data is stored.
2. output_png_dir: The output directory where the plots generated by this function will be saved.
3. output_file_dir: The output directory where the statistics generated by this function will be saved.
4. focus_clusters: A list of cluster numbers (IDs) to focus on.
5. L: The size of the lattice.
6. dP: The value of dP.
7. Geq: The value of Geq.
8. a_vals: A list of a values to consider. This is optional, and if not provided, the function will read the a_vals.txt file from the input directory.
9. S_list: A list of species numbers to consider. Default is [1, 2].
10. verbose: A boolean flag to print verbose output. Default is False.

The function will do the following:
1. Iterate over each a value and each t value, and read the clustered data from the input directory.
2. On reading the clustered data, it will iterate over each species and focus cluster, and do the following:
3. Iterate over each replicate, and do the following:
    a. Mask the clustered data to only consider the focus cluster (non-focus clusters will be set to 0).
    b. Use scipy.ndimage.label to label the focus clusters, and extract the number and size of each cluster.
    c. Create an array with columns containing the counts of each cluster size, normalised by total number of clusters (for a given species and replicate).
    d. Create (or append to) another array with columns containing the counts of each cluster size, normalised by total number of clusters (for a given species and replicate) for all replicates.
e. Create a plot of the cluster size distribution for each species and replicate.
f. Create a plot of the cluster size distribution for all replicates for each species.
g. Save the cluster size distribution data to a CSV file in the output_file_dir directory.
h. Save the plots to the output_png_dir directory.
'''


def analyse_clusters(input_dir, output_png_dir, output_file_dir, n_clusters, L, dP, Geq, maxR, focus_clusters =[0], a_vals =[], S_list=[1,2], verbose=False):

    if len(a_vals) == 0:
        try:
            with open(os.path.join(input_dir, "a_vals.txt"), "r") as a_file:
                a_vals = [float(line.strip()) for line in a_file]
        except FileNotFoundError:
            print("a_vals.txt not found in the input directory. Terminating the program.")
            return
    print(f"List of a values: {a_vals}")
    t_vals = []
    
    for a in a_vals:
        parendir = input_dir
        # Check if the parent directory exists
        if not os.path.exists(parendir):
            print(f"Parent directory {parendir} not found. Skipping...")
            continue
        
        # In parendir, we have files in sub-directories of the form T_{t}/FRAME_T_{t}_a_{a}_R_{R}.csv 
        #where t is a time value (can be a decimal) and R is a replicate number (integers >=0).
        # We first extract the t values from the subdirectories and sort them in ascending order.
        # Then we iterate over each replicate file for each t value and cluster the data (extracting the replicate number from the filename).

        # Get all files in the input directory (parendir) that match the pattern "CLUSTERED_{n_clusters}_FRAME_T_{t}_a_{a}_dP_{dP}_Geq_{Geq}_R_{R}.csv"

        files = glob.glob(parendir + f"CLUSTERED_{n_clusters}_FRAME_T_*_a_{a}_*.csv")
        t_vals = [re.findall(r'T_[\d]*[.][\d]+' , f)[0] if re.findall(r'T_[\d]*[.][\d]+' , f) 
                    else re.findall(r'T_[\d]+' , f)[0] for f in files]
        t_vals = sorted(list(set([float(re.findall(r'[\d]*[.][\d]+', t)[0]) 
                     if re.findall(r'[\d]*[.][\d]+', t) 
                     else int(re.findall(r'[\d]+', t)[0]) for t in t_vals])))
        t_vals = [t for t in t_vals if t >= 45000]  
        print(t_vals)
        
        # Iterate over each t value
        for t in t_vals:
            # All replicates for given a and t are stored in a single unique file
            # CLUSTERED_{n_clusters}_FRAME_{t}_*_a_{a}_dP_{dP}_Geq_{Geq}_R_{maxR}.csv" 

            if maxR == 0:
                # Get the maximum replicate number for the given a and t values
                file = glob.glob(parendir + f"CLUSTERED_{n_clusters}_FRAME_T_{t}_a_{a}_dP_{dP}_Geq_{Geq}_R_*.csv")[0]
                maxR = re.findall(r'[\d]+' , file)[-1]
                print(f"Max Replicates: {maxR}")
            
            filename = f"CLUSTERED_{n_clusters}_FRAME_T_{t}_a_{a}_dP_{dP}_Geq_{Geq}_R_{maxR}.csv"
            file = os.path.join(parendir, filename)
            

            # Read the clustered data from the input directory
            data_df = pan.read_csv(file, delimiter=',', header= 0)
            replicates_df = pan.DataFrame()

            for s in S_list:
                for cluster_type in focus_clusters:
                    for rep_no in range(0, maxR):
                        
                        # If there are no features at all for the given species and replicate, skip to the next replicate
                        if np.count_nonzero(data_df.columns.str.contains(f'Ordn_cluster_S_{s}_R_{rep_no}')) == 0:
                            if(verbose):
                                print(f"For Species {s}, the Replicate {rep_no}, is dead. Skipping...")
                            continue

                        if(verbose):
                            print("\n=====================================================================================================\n")
                            print(f"Processing at a = {a}, t = {t}, Species = {s}, Replicate = {rep_no}, dP = {dP}, Geq = {Geq}, Focus Cluster = {cluster_type}")
                            print("=====================================================================================================\n")

                        # Mask the clustered data to only consider the focus cluster (non-focus clusters will be set to 0)
                        mask = data_df[f'Ordn_cluster_S_{s}_R_{rep_no}'] == cluster_type
                        masked_data = data_df[f'Ordn_cluster_S_{s}_R_{rep_no}'].copy()
                        print(masked_data.head())
                        if(cluster_type == 0):
                            masked_data[mask] = 1
                            if(verbose):
                                print(f"Focus Cluster Type {cluster_type} is zero-labelled, so relabelled to 1.")
                        masked_data[~mask] = 0

                        print(masked_data.head())
                        
                        # Use scipy.ndimage.label to label the focus clusters
                        non_periodic_labels, num_clusters = label(masked_data.to_numpy().reshape(L, L))
                        
                        p,q = non_periodic_labels.shape
                        if verbose:
                            
                            print(f"Non- Periodic Labelled Frame For Species {s}, Replicate {rep_no}, \& Focus Cluster Type  {cluster_type}:\t")
                            print("NON-PERIODIC NON ZERO COUNT:\t" +str(np.count_nonzero(non_periodic_labels)))
                            print("Number of non-zero clusters:\t" +str(num_clusters)) #


                        #CREDITS FOR PERIODIC CONNECTIONS TO: 
                        #https://stackoverflow.com/questions/55953353/how-to-specify-a-periodic-connection-for-features-of-scipy-ndimage-label

                        #Now to create cluster labels that obey periodic boundary conditions.

                        ordered_periodic_labels = non_periodic_labels.copy()

                        for x in range(0, p):
                            if ordered_periodic_labels[x, 0] > 0 and ordered_periodic_labels[x, -1] > 0:
                                ordered_periodic_labels[ordered_periodic_labels == ordered_periodic_labels[x, -1]] = ordered_periodic_labels[x, 0]
                                #Dealing with top and bottom borders.
                        for x in range(0, q):
                            if ordered_periodic_labels[0, x] > 0 and ordered_periodic_labels[-1, x] > 0:
                                ordered_periodic_labels[ordered_periodic_labels == ordered_periodic_labels[-1, x]] = ordered_periodic_labels[0, x]
                                #Dealing with L and R borders.

                        label_ids = np.unique(ordered_periodic_labels.copy()) #Return labels in ascending order.

                        if verbose:
                            print(f"BEFORE REDUCTION, Periodic Labelled Frame For Species {s}, Replicate {rep_no}, \& Focus Cluster Type {cluster_type}:\t")
                            print("PERIODIC NON ZERO COUNT:\t" +str(np.count_nonzero(ordered_periodic_labels)))
                            print("Number of non-zero PERIODIC clusters:\t" +str(len(label_ids)) + " with last label being " + str(label_ids[-1])) 
                            print("PERIODIC LABEL IDS:\t" +str(label_ids[1:int(len(label_ids)/12)+1]) + " ....")

                        for i in  range(1, len(label_ids)):
                            if(label_ids[i] == label_ids[i-1] + 1):
                                continue
                            #Update LABEL ID to lowest extant value possible
                            ordered_periodic_labels[ordered_periodic_labels == label_ids[i]] = label_ids[i-1] + 1
                            label_ids[i] = label_ids[i-1] + 1
                        
                        if verbose:
                            print(f"AFTER REDUCTION, Periodic Labelled Frame For Species {s}, Replicate {rep_no}, \& Focus Cluster Type  {cluster_type}:\t")
                            print("PERIODIC NON ZERO COUNT:\t" +str(np.count_nonzero(ordered_periodic_labels)))
                            print("Number of non-zero PERIODIC clusters:\t" +str(len(label_ids)) + " with last label being " + str(label_ids[-1]))  
                            print("PERIODIC LABEL IDS:\t" +str(label_ids[1:int(len(label_ids)/12)+1]) + " ....")

                        # Here, we have the periodic labels. Now we can proceed to extract the number and size of each cluster.
                        periodic_Clus_ID, periodic_Clus_Size = np.unique(ordered_periodic_labels, return_counts=True)
                        #Here periodic_Clus_ID contains the cluster IDs (of each unique cluster) and periodic_Clus_Size contains the size of each unique cluster.

                        #Next, we need to remove the zeroth cluster, which is the background cluster.
                        # Also to construct frequency distribution of cluster sizes, we need to count all clusters that have the same size.
                        CSize, NS = np.unique(periodic_Clus_Size[1:], return_counts=True)

                        NS = NS/(len(label_ids)-1) #Normalising by total number of clusters.

                        # Create (or append to) another array with columns containing the counts of each cluster size, normalised by total number of clusters for all replicates

                        replicates_df[f'sizeC_S_{s}_R_{rep_no}_fC_{cluster_type}'] = pan.Series(CSize)
                        replicates_df[f'NS_S_{s}_R_{rep_no}_fC_{cluster_type}'] = pan.Series(NS)
                        

                        if(False):
                            outdir = output_png_dir + f"Sp_{s}/Indiv_Replicates"
                            #Make the output directory for images if it does not exist
                            Path(outdir).mkdir(parents=True, exist_ok=True)
                            # Create a plot of the cluster size distribution for each species and replicate and focus cluster
                            fig, ax = plt.subplots(1, 1, figsize=(8, 6))
                            sea.scatterplot(data=replicates_df, x=f'sizeC_S_{s}_R_{rep_no}_fC_{cluster_type}', y=f'NS_S_{s}_R_{rep_no}_fC_{cluster_type}', ax=ax, s=2, alpha=0.5)
                            sea.lineplot(data=replicates_df, x=f'sizeC_S_{s}_R_{rep_no}_fC_{cluster_type}', y=f'NS_S_{s}_R_{rep_no}_fC_{cluster_type}', ax=ax, alpha=0.2)
                            ax.set_xlabel('Cluster Size (s)')
                            ax.set_ylabel(r'Normalised Frequency $N_{%s}(T_{eq} = %s)$' %(s, t))
                            ax.set_title(r'$N_{%s}(T_{eq} = %g)$ vs $s$ for Focus $ = %g$, $ R = %g$, Rep $ = %g$, $dP = %g$' %(s, t, cluster_type, a, rep_no, dP))
                            plt.yscale('log'); plt.xscale('log');
                            output_path = os.path.join(outdir, f'Ns_fC_{cluster_type}_S_{s}_T_{t}_a_{a}_R_{rep_no}.png')
                            plt.savefig(output_path)
                            plt.close()
                            #plt.show(); plt.close()
                    # End of rep-loop

                    # Create a plot of the cluster size distribution for all replicates for each species and focus cluster
                    # Use seaborn to create a scatter plot followed by a line plot, with a different hue for each replicate data.
                    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
                    colours = sea.color_palette("husl", maxR)
                    for rep_no in range(0, maxR):
                        #Plot only non
                        ax.scatter(replicates_df[f'sizeC_S_{s}_R_{rep_no}_fC_{cluster_type}'], replicates_df[f'NS_S_{s}_R_{rep_no}_fC_{cluster_type}'], s=2, alpha=0.3, label=f'Replicate {rep_no}', color=colours[rep_no])
                        sea.lineplot(data=replicates_df, x=f'sizeC_S_{s}_R_{rep_no}_fC_{cluster_type}', y=f'NS_S_{s}_R_{rep_no}_fC_{cluster_type}', ax=ax, alpha=0.2, color=colours[rep_no])
                    ax.set_xlabel('Cluster Size (s)')
                    ax.set_ylabel(r'Normalised Frequency $N_{%s}(T_{eq} = %s)$' %(s, t))
                    #ax.set_ylim(ymax= 1.0)
                    #ax.set_xlim(xmax= 1000)
                    header_focus = "Lowest Act. Cluster" if cluster_type == 0 else "Highest Act. Cluster" if cluster_type == n_clusters-1 else f"Med Act Cluster {cluster_type}"

                    ax.set_title(r'$N_{%s}(T_{eq} = %g)$ vs $s$, %s, $ R = %g$, $dP = %g$, $n = %g$' %(s, t, header_focus, a, dP, maxR))
                    try:
                        plt.yscale('log'); plt.xscale('log');

                        outdir = output_png_dir + f"Sp_{s}"
                        #Make the output directory for images if it does not exist
                        Path(outdir).mkdir(parents=True, exist_ok=True)
                        output_path = os.path.join(outdir, f'Ns_fC_{cluster_type}_S_{s}_T_{t}_a_{a}_R_{maxR}.png')
                        plt.savefig(output_path, dpi= 500)
                        #plt.show(); plt.close()
                        plt.close()
                    except ValueError:
                        plt.close()
                        print("ValueError: math domain error. Skipping...")
                        continue
                
                # End of cluster_type-loop
            # End of s-loop
            # Save the cluster size distribution data to a CSV file in the output_file_dir directory
            outdir = output_file_dir
            Path(outdir).mkdir(parents=True, exist_ok=True)
            output_path = os.path.join(outdir, f'ClusterFreqDistr_T_{t}_a_{a}_R_{maxR}.csv')
            replicates_df.to_csv(output_path, index=False)
        
        # End of t-loop
    # End of a-loop

    # Finally create a video for each species and focus clusters using the images generated in the output_png_dir directory
    # The video goes over all the a-values for a given value of t and focus cluster.
    # The video is saved in the output_png_dir directory using cv2.videowriter.

    for s in S_list:
        for cluster_type in focus_clusters:
            images = []
            pngdir = output_png_dir + f"Sp_{s}"
            savedir = output_png_dir + f"Sp_{s}/Videos"
            Path(savedir).mkdir(parents=True, exist_ok=True)
            t_max = t_vals[-1]
            for a in a_vals:
                filename = f'Ns_fC_{cluster_type}_S_{s}_T_{t_max}_a_{a}_R_{maxR}.png'
                if os.path.exists(os.path.join(pngdir, filename)):
                    images.append(cv2.imread(os.path.join(pngdir, filename)))
            
            video_name = f'ClusterFreqDistr_S_{s}_fC_{cluster_type}_T_{t_max}_R_{maxR}.mp4'
            output_path = os.path.join(savedir, video_name)
            height, width, layers = images[0].shape

            video = cv2.VideoWriter(output_path, cv2.VideoWriter_fourcc(*'mp4v'), fps=1, frameSize=(width, height))

            for image in images:
                video.write(image)
            
            cv2.destroyAllWindows()
            video.release()






                    




'''
Conditions for generating clusters through the generate_clusters function (Unsupervised Clustering):
n_clusters = 4
input_dir = r'StdParam_20_100_EQ/DiC/'
output_png_dir = r'../../Images/Clustered_Images/' + input_dir + f"Clus_{n_clusters}/"
output_file_dir = r'../Data/Clusters/' + input_dir + f"Clus_{n_clusters}/"
'''

# Conditions for analysing clusters through the analyse_clusters function.
n_clusters = 4
input_dir = f'../Data/Clusters/StdParam_20_100_EQ/DiC/Gauss0.75_Clus_{n_clusters}/'
output_png_dir = f'../../Images/ClusterAggrStats/StdParam_20_100_EQ/DiC/Gauss0.75_nClus_{n_clusters}/'
output_file_dir = f'../Data/StdParam_20_100_EQ/DiC/AggrStats/Gauss0.75_nClus_{n_clusters}/'


#Make the output directory for images if it does not exist
Path(output_png_dir).mkdir(parents=True, exist_ok=True)
#Make the output directory for files if it does not exist
Path(output_file_dir).mkdir(parents=True, exist_ok=True)

sea.set_theme(style='whitegrid')
L =128; dP= 30000; Geq =4.802; Sp_list = [1, 2]; ver= True
maxR = 0; focal_clusters = [0, n_clusters-1]

analyse_clusters(input_dir, output_png_dir, output_file_dir, n_clusters, L, dP, Geq, maxR=5, focus_clusters=focal_clusters, S_list=Sp_list, verbose=ver)
#input_dir, output_png_dir, output_file_dir, n_clusters, L, dP, Geq, maxR, focus_clusters =[0], a_vals =[], S_list=[1,2], verbose=False
#generate_clusters(input_dir, output_png_dir, output_file_dir, n_clusters, L, dP, Geq, S_list=Sp_list, verbose=ver)


