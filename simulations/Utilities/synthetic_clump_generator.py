import numpy as np
import pandas as pan
import scipy.ndimage as ndimage
import regex as re
import matplotlib.pyplot as plt
import os
import glob
import random

SPB =3
prefix ="HEXBLADE";
L = 256; r =40;
output_dir = f"../Input/Rietkerk/{SPB}Sp/{prefix}/L_{L}_a_0/"

# Create the output directory if it does not exist
os.makedirs(output_dir, exist_ok=True)


hex_length = {"  P(x; t)": 20}
gauss_2sd_radius = {"  P(x; t)": 4}
gauss_amp = {"  P(x; t)": (10000, 10000)} #Tuple of (max, min) amplitudes
min_val = {"  P(x; t)": 0.0}
percent_missing = {"  P(x; t)": 0.1}
retain_cluster = {"  P(x; t)": 0}
gradient_type = {"  P(x; t)": None} #Valid arguments: "linear_X", "linear_-X", "linear_Y", "linear_-Y", "radial", "Random" or None
radial_gradient_origin = {"  P(x; t)": (L//2, L//2)}



def create_gaussian_nonperiodic(L, center, amplitude, sigma, min_val):
    """Create a 2D Gaussian centered at `center` in a grid of size `L`."""
    y, x = np.indices((L, L))
    y0, x0 = center
    # Calculate the distance with periodic boundary conditions
    gaussian = (amplitude) * np.exp( -(np.power(x - x0, 2)+ np.power(y - y0,2)) / (2*np.power(sigma,2)))
    #gaussian = (amplitude)/(2*np.power(sigma,2)*np.pi) * np.exp( -(np.power(x - x0, 2)+ np.power(y - y0,2)) / (2*np.power(sigma,2)))
    gaussian[gaussian < min_val] = min_val
    return gaussian

def create_gaussian_periodic(L, center, amplitude, sigma, min_val):
    """Create a 2D Gaussian centered at `center` in a grid of size `L`."""
    y, x = np.indices((L, L))
    y0, x0 = center
    # Calculate the distance with periodic boundary conditions
    dx = np.minimum(np.abs(x - x0), L - np.abs(x - x0))
    dy = np.minimum(np.abs(y - y0), L - np.abs(y - y0))
    gaussian = (amplitude/2*np.power(sigma,2)*np.pi) * np.exp( -(np.power(dx, 2)+ np.power(dy,2)) / (2*np.power(sigma,2)))
    gaussian[gaussian < min_val] = min_val
    return gaussian

''' SUMMARY OF tile_hex_grid_basic()
    Create a hexagonal grid pattern.
    If gradient is not None, the hexagon centers are assigned values based on the gradient.
    gradient has the following valid arguments: "linear_+x", "linear_-x", "linear_+y", "linear_-y", "radial" or None.
    radial_gradient_origin is the origin of the gradient (only used if gradient is set to "radial").
    The function returns a grid of size L x L with hexagons assigned values based on the gradient.
    "linear" gradients are linearly increasing/decreasing in the (+/-) x or y direction.
    gradient_range is a tuple (min, max) that specifies the range of the gradient.
    "radial" gradient is a radial gradient with the origin at radial_gradient_origin.
'''
def tile_hex_grid_basic(L, hex_length, gradient=None, gradient_range=(0.1, 1), radial_gradient_origin=(0, 0)):
    """Create a hexagonal grid pattern."""
    grid = np.zeros((L, L))
    hex_height = np.sqrt(3) * hex_length
    hex_width = 2 * hex_length
    
    for i in range(int(hex_height), L, int(hex_height)):
        for j in range(int(hex_width / 2), L, int(hex_width)):
            if (i // int(hex_height)) % 2 == 1:
                center = (i, j)
                #print("On-tile", i, j)
            else:
                #print("Off-tile", i, j, int(hex_width / 2))
                center = (i, j + int(hex_width / 2))
            #grid[center] = 1 #FOR TRUE HEXAGON CENTERS
            grid[i%L, j%L] = 1 # FOR SIMPLIFIED (SQUARE) HEXAGON ("HEXBLADE") CENTERS

    if gradient != "Random" and gradient is not None:
        print("Gradient: ", gradient)
        # Assign values to hexagon centers defined previously based on the gradient
        y, x = np.indices((L, L))
        if gradient == "linear_X":
            # Linear gradient in the +x direction, with max values assigned to the rightmost hexagon centres (with max x-coordinate)
            dx = np.abs(x - np.max(np.argwhere(grid > 0)[:, 1])); max_dx = np.max(dx)
            linear_gradient = np.interp(dx, (0, max_dx), gradient_range)
            # Assign the linear gradient to the hexagon centers
            grid[grid > 0] = linear_gradient[grid > 0]

        elif gradient == "linear_-X":
            # Linear gradient in the -x direction, with max values assigned to the leftmost hexagon centres (with min x-coordinate)
            dx = np.abs(x - np.min(np.argwhere(grid > 0)[:, 1])); max_dx = np.max(dx);
            linear_gradient = np.interp(dx, (0, max_dx), gradient_range)
            grid[grid > 0] = linear_gradient[grid > 0]
            '''
            descending_centers = np.sort(np.argwhere(grid > 0), axis=0)[::-1]
            min_x = np.min(descending_centers[:, 1]); max_x = np.max(descending_centers[:, 1])
            for center in descending_centers:
                grid[center[0], center[1]] = np.interp(center[1], (min_x, max_x), gradient_range)
            '''
        elif gradient == "linear_Y":
            # Linear gradient in the +y direction, with max values assigned to the topmost hexagon centres (with max y-coordinate)
            dy = np.abs(y - np.max(np.argwhere(grid > 0)[:, 0])); max_dy = np.max(dy)
            linear_gradient = np.interp(dy, (0, max_dy), gradient_range)
            grid[grid > 0] = linear_gradient[grid > 0]
            '''
            ascending_centers = np.sort(np.argwhere(grid > 0), axis=0)
            max_y = np.max(ascending_centers[:, 0]); min_y = np.min(ascending_centers[:, 0])
            for center in ascending_centers:
                grid[center[0], center[1]] = np.interp(center[0], (min_y, max_y), gradient_range)
            '''
        elif gradient == "linear_-Y":
            # Linear gradient in the -y direction, with max values assigned to the bottommost hexagon centres (with min y-coordinate)
            dy = np.abs(y - np.min(np.argwhere(grid > 0)[:, 0])); max_dy = np.max(dy)
            linear_gradient = np.interp(dy, (0, max_dy), gradient_range)
            grid[grid > 0] = linear_gradient[grid > 0]
            '''
            descending_centers = np.sort(np.argwhere(grid > 0), axis=0)[::-1]
            min_y = np.min(descending_centers[:, 0]); max_y = np.max(descending_centers[:, 0])
            for center in descending_centers:
                grid[center[0], center[1]] = np.interp(center[0], (min_y, max_y), gradient_range)
            '''
        elif gradient == "radial":
            # Radial gradient with gradient origin at radial_gradient_origin
            # Gradient values decrease radially from the origin
            y, x = np.indices((L, L))
            dx = np.minimum(np.abs(x - radial_gradient_origin[1]), L - np.abs(x - radial_gradient_origin[1]))
            dy = np.minimum(np.abs(y - radial_gradient_origin[0]), L - np.abs(y - radial_gradient_origin[0]))
            radial_gradient = np.interp(np.sqrt(np.power(dx, 2) + np.power(dy, 2)), (0, L/np.sqrt(2)), (gradient_range[1], gradient_range[0]))
            # Assign the radial gradient to the hexagon centers
            grid[grid > 0] = radial_gradient[grid > 0]
        else:
            raise ValueError("Invalid gradient argument. Valid arguments: 'linear_X', 'linear_-X', 'linear_+Y', 'linear_-Y', 'radial', 'Random' or None")
    
    return grid

def plot_hex_grid(L, hex_length, gradient=None, gradient_range=(0.1, 1), radial_gradient_origin=(0, 0)):
    """Plot the hexagonal grid pattern."""
    grid = tile_hex_grid_basic(L, hex_length, gradient, gradient_range, radial_gradient_origin)
    # Plotting hexagon centers
    plt.imshow(grid, cmap='hot', interpolation='nearest')
    # Next add the hexagon patterns
    plt.show()
    plt.close()



def generate_grid(L, hex_length, gauss_2sd_radius, gauss_amp_range, retain_cluster=0, percent_missing=0, min_val=0, gradient=None, radial_gradient_origin=(0, 0)):
    if gauss_2sd_radius >= hex_length:
        raise ValueError("gauss_2sd_radius must be less than hex_length")
    
    if len(gauss_amp_range) != 2:
        raise ValueError("gauss_amp_range must be a tuple of length 2")
    
    # Get the Gaussian SD.
    gauss_sd = gauss_2sd_radius / 2.0
    
    grid = np.zeros((L, L))
    if gauss_amp_range[0] < gauss_amp_range[1]:
        gradient_range = (gauss_amp_range[0], gauss_amp_range[1])
    else:
        gradient_range = (0.1, 1)
    hex_grid = tile_hex_grid_basic(L, hex_length, gradient, gradient_range, radial_gradient_origin=radial_gradient_origin)
    #plot_hex_grid(L, hex_length, gradient, gradient_range, radial_gradient_origin=radial_gradient_origin)
    # Identify hexagon centers (non-zero values in hex_grid)
    centers = np.argwhere(hex_grid > 0)
    num_hex = len(centers)
    # If retain_cluster is set, randomly select a cluster of hexagons to retain, equal to retain_cluster
    if retain_cluster > 0 and retain_cluster <= num_hex:
        selected_centers = random.sample(list(centers), retain_cluster)
    else:
        # Randomly select (1 - percent_missing) fraction of hexagon centers
        selected_centers = random.sample(list(centers), int((1 - percent_missing) * num_hex))
    
    for center in selected_centers:
        print(center)
        print(hex_grid[center[0],center[1]], gradient)
        if(gradient is not None and gradient != "Random"):
            gauss_amp = hex_grid[center[0],center[1]]
        elif(gauss_amp_range[0] < gauss_amp_range[1] and gradient == "Random"):
            gauss_amp = np.random.uniform(gauss_amp_range[0], gauss_amp_range[1])
        else:
            gauss_amp = gauss_amp_range[0]

        print(gauss_amp)
        
        gaussian = create_gaussian_nonperiodic(L, center, gauss_amp, gauss_sd, min_val)
        grid += gaussian
    
    grid[grid < min_val] = min_val
    return grid

def plot_grid(grid, title=""):
    # Add a colorbar to the plot
    plt.imshow(grid, cmap='hot', interpolation='nearest')
    cbar = plt.colorbar()
    plt.title(title)
    plt.show()
    plt.close()

def save_grid(df, replicate, savedir, frame_subdir ="SEP_{sep}_WID_{wid}/{grad}_GR/MSF_{missing}/MAMP_{mamp}_MIMP_{mimp}/", overwrite=False):
    """Save the grid to CSV, ensuring no file overwrite."""
    # Creating the directory structure for the frame.
    first_key = list(df.keys())[2]
    save_frame_dir = os.path.join(savedir, frame_subdir.format(sep = hex_length[first_key], missing = percent_missing[first_key], grad = gradient_type[first_key],
            wid = gauss_2sd_radius[first_key], mamp = gauss_amp[first_key][0], mimp= gauss_amp[first_key][1], minval = min_val[first_key]))
    os.makedirs(save_frame_dir, exist_ok=True)
    
    existing_files = glob.glob(os.path.join(save_frame_dir, 'FRAME_T_0_a_0_R_*.csv'))
    existing_replicates = [int(re.findall(r'[\d]+' , f)[-1]) for f in existing_files]
    
    if replicate in existing_replicates and not overwrite:
        replicate = max(existing_replicates) + 1

    
    
    filename = os.path.join(save_frame_dir, f"FRAME_T_0_a_0_R_{replicate}.csv")
    df.to_csv(filename, index=False)

def main(hex_length, gauss_2sd_radius, L, r, gauss_amp, min_val, output_dir, percent_missing, gradient_type, radial_gradient_origin):
    for replicate in range(0, r):
        # Create a pandas dataframe with columns: "a_c", "x" followed by keys in hex_length
        # Each column should be a 1D array of length L x L
        df = pan.DataFrame(columns=["a_c", "x"] + list(hex_length.keys()))
        #Store a 0 array of length L x L in the "a_c" column
        df["a_c"] = np.zeros(L*L); df["x"] = np.arange(L*L); 
        first_key = list(hex_length.keys())[0]
        for key in hex_length:
            grid = generate_grid(L, hex_length[key], gauss_2sd_radius[key], gauss_amp[key], percent_missing = percent_missing[key], min_val = min_val[key], 
                                 retain_cluster = retain_cluster[key], gradient = gradient_type[key], radial_gradient_origin = radial_gradient_origin[key])
            #plot_grid(grid, key)
            #Flatten the grid and store it in the respective column
            df[key] = grid.flatten()
            print("Max value: ", np.max(grid))
            #grid = df[key].values.reshape(L, L)]
            #plot_grid(grid, key)
        print(gradient_type)
        if(gradient_type[first_key] == None or gauss_amp[first_key][0] == gauss_amp[first_key][1]):
            framesubdir = "SEP_{sep}_WID_{wid}/MSF_{missing}/MAMP_{mamp}_MINVAL_{minval}/"
        else:
            framesubdir = "SEP_{sep}_WID_{wid}/{grad}_GR/MSF_{missing}/MAMP_{mamp}_MIMP_{mimp}/"
        print(framesubdir)
        save_grid(df, replicate, output_dir,frame_subdir = framesubdir, overwrite=False)

if __name__ == "__main__":
    # Example usage
    main(hex_length, gauss_2sd_radius, L, r, gauss_amp, min_val, output_dir, percent_missing, gradient_type, radial_gradient_origin)
