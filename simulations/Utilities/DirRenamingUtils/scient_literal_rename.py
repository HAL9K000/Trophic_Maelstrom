import os
import re
from pathlib import Path

def convert_scientific_to_decimal(match, scale_factor):
    """Converts a scientific literal (e.g., 2e-06) to decimal, scaled by a factor."""
    scientific_literal = match.group(1)  # Extract the scientific notation value (e.g., 2e-06)
    decimal_value = float(scientific_literal) * scale_factor  # Convert to decimal and scale
    #Decimal value is rounded to 5 decimal places if necessary
    decimal_value = round(decimal_value, 5); decimal_value = "{:g}".format(decimal_value)
    return f"a_{decimal_value}"

def rename_with_scientific_to_decimal(path, scale_factor):
    """
    Renames directories or files if their name contains a scientific literal in the form a_{scientific_literal}.
    """
    # Regex pattern to match 'a_{scientific_literal}'
    pattern = r'a_([+-]?\d+(?:\.\d+)?[eE][+-]?\d+)'

    # Check if the path's name contains a scientific literal and needs renaming
    if re.search(pattern, path.name):
        # Replace the scientific literal in the name with its decimal equivalent
        new_name = re.sub(pattern, lambda match: convert_scientific_to_decimal(match, scale_factor), path.name)
        new_path = path.with_name(new_name)

        # Rename the file or directory
        print(f"Renaming: {path} -> {new_path}")
        path.rename(new_path)
        return new_path  # Return the new path for further renaming if necessary
    return path  # If no change, return the original path

def process_directory(directory_path, scale_factor):
    """Recursively processes all directories and files in the directory to rename if necessary."""
    path = Path(directory_path)

    # Traverse the directory tree in reverse order to ensure sub-directories are processed before their parents
    for item in sorted(path.rglob('*'), key=lambda p: len(p.parts), reverse=True):
        if item.is_file() or item.is_dir():
            rename_with_scientific_to_decimal(item, scale_factor)

    # Finally, process the root directory itself
    rename_with_scientific_to_decimal(path, scale_factor)

if __name__ == "__main__":
    # Input: Directory path and scale factor
    directory_path = "../../Data/DP/Prelims/Stochastic/2Sp/"
    scale_factor = float(input("Enter the scaling factor (e.g., 1000): "))

    if os.path.isdir(directory_path):
        process_directory(directory_path, scale_factor)
    else:
        print(f"Error: The provided path '{directory_path}' is not a valid directory.")