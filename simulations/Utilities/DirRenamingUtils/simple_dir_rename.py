import os
from pathlib import Path
import re
import glob


''' This script is designed to rename files within specified directories (given by a list of dir paths)
    by replacing file names which match a specific pattern (e.g., 'a_1.728_' or 'a_2e-06_') with a new name that 
    where the matched pattern is replaced by a new user provided pattern (e.g., 'a_1.728_scaled_' or 'a_0.000002_scaled_'),
    while preserving the rest of the file name. The script uses regular expressions to identify the pattern in the file names and performs the renaming operation accordingly.
'''

def rename_files_in_directory(dir_path, old_pattern, new_pattern, user_confirmation=True):

    # Use glob to find all files in the directory that match the old pattern
    files_to_rename = glob.glob(os.path.join(dir_path, f"*{old_pattern}*"))

    print("\n====================================================================\n")
    print(f"⚠️ Found {len(files_to_rename)} files to rename in '{dir_path}' matching the pattern '{old_pattern}'.")
    # Print immediate parent directory of the files to be renamed for user confirmation
    parent_dir = os.path.basename(os.path.normpath(dir_path))
    print(f"In Parent Directory: '{parent_dir}'...")
    print("\n====================================================================\n")

    for file_path in files_to_rename:
        # Extract the file name from the path
        file_name = os.path.basename(file_path)

        # Replace the old pattern with the new pattern in the file name
        new_file_name = file_name.replace(old_pattern, new_pattern)

        # Create the new file path
        new_file_path = os.path.join(dir_path, new_file_name)

        # Rename the file
        if user_confirmation:
            print(f"⚠️ Renaming: {file_name} -> {new_file_name}.\t Proceed? (Y/y/N/n): ", end="")
            user_input = input()
            if user_input.lower() != 'y':
                print("⚠️ Skipping renaming.")
                continue
        try:
            os.rename(file_path, new_file_path)
            print(f"✅ Renamed: {file_name} -> {new_file_name}")
        except Exception as e:
            print(f"❌ Error renaming {file_path} to {new_file_path}: {e}")


if __name__ == "__main__":
    # Example usage
    dir_paths = [
        "../../Data/DP/Frames/Stochastic/3Sp_mG_20_mP_100/DsB6L9-UA0A0-1UNI_0.00172-0.00173_dP_1_Geq_0/",
    ]
    old_pattern = "_a_1.728_"  # Example old pattern to be replaced
    new_pattern = "_a_1.7275_"  # Example new pattern to replace with
    user_confirmation = False  # Set to False to skip user confirmation for each file
    for dir_path in dir_paths:
        rename_files_in_directory(dir_path, old_pattern, new_pattern, user_confirmation)
