    in_par_prelims_file="$1"

# Convert time gap from hours to seconds
time_gap_seconds=$(( time_gap_hours * 3600 ))

# Function to log errors
log_error() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') ERROR: $1" >> error_log.txt
}

# Function to log updates
log_update() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') UPDATE: $1" >> error_log.txt
}
    
    while IFS= read -r line || [[ -n "$line" ]]; do
        # Extract input parameters from the line (assuming well-formed lines)
        indir=$(echo "$line" | grep -oP '(?<=--indir )[^ ]+')
        outdir=$(echo "$line" | grep -oP '(?<=--outdir )[^ ]+')
        
        if [[ -z "$indir" || -z "$outdir" ]]; then
            log_error "Invalid line: $line"
            continue
        fi

        # Remove trailing whitespaces and backslash (/) from outdir and indir
        outdir=$(echo "$outdir" | sed 's/[[:space:]]*\/$//' | sed 's/[[:space:]]*$//')
        indir=$(echo "$outdir" | sed 's/[[:space:]]*\/$//' | sed 's/[[:space:]]*$//')

                # Run the Python script
        python3.11 reorganise_prelims_dir.py $line
        if [[ $? -ne 0 ]]; then
            log_error "Python script reorganise_prelims_dir.py failed for line: $line"
            continue
        fi

    done < "$in_par_prelims_file"