#!/bin/bash

# This script is used to prepare a job array script for Amarel Cluster.
# It accepts two arguements: the input file and the job_name.
# The input file is a tab delineated txt file with the following columns:
#dt	t_max	L	R	a_st	a_end	div	dP	PREFIXES,
# where the eight arguement on each line (corresponding to div) 
# gives the --cpus-per-task parameter for each job in the array.
# The job_name is the name of the job array script to be created.
# Finally the number of lines in the input file (minus the header) is used to determine the number of jobs in the array.
# Additionally, both the number of jobs and the number of cpus-per-task can be provided as optional arguements,
# in which case they will override the values from the input file.
# After reading the necessary parameters, the actual script will be read using a here document to the submission script $(job_name)_array.sh

# Check if the number of arguments is at least 2
if [ $# -lt 3 ]; then
    echo "Usage: $0 <path/to/init_file.txt> <job_name> <SPB> [Optional <init =0> <n> <k>]"
    echo "Where: n is --cpus-per-task and k is the number of jobs in the array."
    exit 1
fi

#Check if the input file exists
if [ ! -f $1 ]; then
    echo "Error: Input file not found. Are you sure the path is correct?. Exiting."
    exit 1
fi


# Extract n and k from the arguments if they are provided
if [ $# -ge 6 ]; then
    k=$6
else
    # Read k from the init file (number of lines in the file - 1)
    k=$(wc -l < $1)
fi

if [ $# -ge 5 ]; then
    n=$5
else
    # Read n as the maximum value of div (seventh column) in the init file
    n=$(awk 'NR>1 {print $7}' $1 | sort -nr | head -n 1)
fi
if [ $# -ge 4 ]; then
    init=$4
else
    init=0
fi

# All the parameters are provided
echo "k = $k, n = $n", init = $init


# Generate the SLURM script with the correct array range
# Generate the SLURM script with the correct array range
cat << EOF > ${2}_array.sh
#!/bin/bash
#SBATCH --job-name=${2}     
#SBATCH --partition=p_deenr_1      
#SBATCH --requeue  
#SBATCH --priority=100
#SBATCH --array=0-$((k))%5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$((n))      # IDEALLY THIS SHOULD BE DYNAMICALLY SET, VARYING FOR EACH JOB.
#SBATCH --mem-per-cpu=8G          # Memory per CPU core
#SBATCH --time=14-00:00:00      # Total run time limit (DD-HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH --mail-type=END
#SBATCH --mail-user=xXx@uni.gov

# Path to the input file containing the arguments
init_file=$1
# Read the corresponding line of arguments based on SLURM_ARRAY_TASK_ID
read -r p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 <<< \$(sed -n "\$((SLURM_ARRAY_TASK_ID+2))p" \$init_file)
CPUS_PER_TASK=\$p7
#Strip whitespaces and newlines from p12
p12=\$(echo \$p12 | tr -d '[:space:]')
prefix=\$p12
#Remove trailing whitespaces and new lines from p9 if SpB is 1 else remove trailing whitespaces and new lines from p10
if [ $3 -eq 1 ]; then
    p9=\$(echo \$p9 | tr -d '[:space:]')
    prefix=\$p11
elif [ $3 -eq 2 ]; then
    p11=\$(echo \$p11 | tr -d '[:space:]')
    prefix=\$p11
else
    p12=\$(echo \$p12 | tr -d '[:space:]')
    prefix=\$p12
fi
if [ $init -eq 2 ]; then
    # Check if p13 is not empty
    if [ -z "\$p13" ]; then
        p12=$(echo \$p12 | tr -d '[:space:]')
    else
        p13=$(echo \$p13 | tr -d '[:space:]')
    fi
fi

echo \$p1 \$p2 \$p3 \$p4 \$p5 \$p6 \$p7 \$p8 \$p9 \$p10 \$p11 \$p12 \$p13
# Dynamically set the number of CPUs per task

# Move to the directory where the source files are located
cd ..

# Use the variables p1, p2, p3, p4, p5, p6, p7, p8, p9 to compile the source files:
if [ $init -ne 2 ]; then
    g++ -O3 -march=native -DSPB=${3} -DINIT=${init} multiSPDP.cpp order_${3}stocDP_unity.cpp -L\${HOME}/Installs/fftw-3.3.10/lib -lfftw3_threads -lfftw3 -lm -fopenmp -o ama_${2}_\${SLURM_ARRAY_TASK_ID}.out -std=c++23
else
    g++ -O3 -march=native -DSPB=${3} -DINIT=${init} multiSPDP.cpp order_${3}stocDP_burnin.cpp -L\${HOME}/Installs/fftw-3.3.10/lib -lfftw3_threads -lfftw3 -lm -fopenmp -o ama_${2}_\${SLURM_ARRAY_TASK_ID}.out -std=c++23
fi

# Run the compiled program with the input parameters
./ama_${2}_\${SLURM_ARRAY_TASK_ID}.out \$p1 \$p2 \$p3 \$p4 \$p5 \$p6 \$p7 \$p8 \$p9 \$p10 \$p11 \$p12 \$p13 &> std_${2}_\${SLURM_ARRAY_TASK_ID}.txt

// Finally, move back to the submission directory
cd \$SLURM_SUBMIT_DIR
EOF

sbatch ${2}_array.sh

# Make the script executable
# Print the name of the generated script
echo "Generated script: ${2}_array.sh"

# End of script