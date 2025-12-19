#!/bin/bash

# This script is used in conjunction with the Linux utility screen and nvcc (CUDA compiler).
# It is used to run a series of jobs in parallel, each in a separate screen session.
# It accepts upto 5 arguments: the input file, SpB (the number of biotic species), Init_Types (type of frame initialisation), 
# Mv_Invariance (type of movement invariance), and ARCH-SM (CUDA architecture and streaming multiprocessor).
# The input file is a tab delineated txt file with the following columns:
#dt	t_max	L	R	a_st	a_end	div	dP	PREFIXES
# where the seventh arguement on each line (corresponding to div) gives the number of threads.
# Each line in the input file corresponds to a separate job, which will be run in a separate screen session.
# An array of screen session names is provided, and the script will check if each screen session is running before starting a new one.
# If the screen session is not running, the script will start a new screen session and run the job.
# If the screen session is running, the script will skip to the next screen session name in the array.

screen_names=(
    "hongkong" "tokyo" "kyoto" "manilla" 
    "osaka" "busan" "edo" "macao"
    "yokohama" "shizuoka" "nanjing" "chiba"
    "taipei" "shanghai" "nagoya" "hanoi" 
    "shenzhen" "seoul" "incheon" "okinawa"
    "kobe" "sapporo" "fukuoka" "sendai"
    "casablanca" "tangiers" "tunis" "fes" 
    "maghreb" "marrakesh" "rabat" "algiers" 
    "medea" "gibraltar" "dakar" "atlas" 
    "kairouan" "tripoli" "cyrene" "sahara"
    "tuareg" "tombouctou" "djemila" "timbuktu")

# Check if the number of arguments is at least 1
if [ $# -lt 2 ]; then
    echo "Usage: $0 <Path/to/init_file.txt> [Optional: <SpB> <Init_Types> <Mv_Invariance> <ARCH-SM>]"
    echo " where SpB is the number of biotic species (default is 3), 
            Init_Types is type of frame initialisation (default is 1),
            Mv_Invariance is the type of movement invariance (default is 0 = Distance Inv, 1 = Velocity Inv) and \n
            ARCH-SM is the virtual architecture and streaming multiprocessor of CUDA-compatible device (default is left to nvcc to decide)."

    exit 1
fi

# Check if the screen utility is installed.
if ! command -v screen &> /dev/null; then
    echo "Error: screen utility is not installed. Please install screen and try again."
    exit 1
fi

#Check if the input file exists
if [ ! -f $1 ]; then
    echo "Error: Input file not found. Are you sure the path is correct?. Exiting."
    exit 1
fi

# Check if there is a CUDA-compatible device
if ! command -v nvcc &> /dev/null; then
    echo "Error: nvcc is not installed. Please install CUDA and try again, or try running the non-CUDA version of the script."
    exit 1
fi

# If the number of arguments is less than 3, set the default value of SpB to 3
if [ $# -lt 3 ]; then
    SPB=3
else
    SPB=$2
fi

# If the number of arguments is less than 4, set the default value of Init_Types to 1
if [ $# -lt 4 ]; then
    INIT=1
else
    INIT=$3
fi

# If the number of arguments is less than 5, set the default value of Mv_Invariance to 0
if [ $# -lt 5 ]; then
    MV_INVARIANCE=0
else
    MV_INVARIANCE=$4
fi

# If the number of arguments is greater than 5, set the ARCH_SM flag to the specified value, else set it to an empty string
if [ $# -ge 5 ]; then
    ARCH_SM=$5
fi



# Store the current directory
curr_dir=$(pwd)

#Get timestamp at the initialisation of the script in format DD{Month}YYYY
timestamp=$(date +"%d%b%Y")
echo "Starting screen sessions for jobs in $1 with SpB=$2 and Init_Types=$3 at $timestamp ...."

# First, create the object files for the cpp scripts and for the CUDA scripts
#if [ $INIT -ne 2 ]; then
    # Compile order_${SPB}stoc_unity_rietkerk.cpp
#    bash -c "cd .. ; g++-14 -DSPB=${SPB} -DINIT=${INIT} -DMV_INVARIANCE=${MV_INVARIANCE} -DBARRACUDA -c order_${SPB}stoc_unity_rietkerk.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include"
#else
#    bash -c "cd .. ; g++-14 -DSPB=${SPB} -DINIT=${INIT} -DMV_INVARIANCE=${MV_INVARIANCE} -DBARRACUDA -c order_${SPB}stoc_burnin_rietkerk.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include"
#fi
#echo "Created object file for order_${SPB}stoc ... script"
# Next, compile the main cpp script
#bash -c "cd .. ; g++-14 -DSPB=${SPB} -DINIT=${INIT} -DMV_INVARIANCE=${MV_INVARIANCE} -DBARRACUDA -c rietkerk_bjork_FKE.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include"
#echo "Created object file for rietkerk_bjork_FKE ... script"

#Next create object files for the CUDA scripts, with ARCH-SM flag if specified, else no flag
if [ -z "$ARCH_SM" ]; then
    bash -c "cd .. ; nvcc -O3 --use_fast_math -c -o kernel_FKE.o  kernel_advdiff_FKE2D.cu -std=c++20;"
else
    bash -c "cd .. ; nvcc -arch=sm_${ARCH_SM} -O3 --use_fast_math -c -o kernel_FKE.o  kernel_advdiff_FKE2D.cu -std=c++20;"
fi
echo "Created object file for kernel_FKE ... script"

# Function that checks if starts a new screen session and runs the current job after first checking if a screen with the same name is already running.
# If a screen session with the same name is running, it will increment the screen_idx_mod variable to use the next screen name in the array.
# It will also issue a warning message to the user.
# The function takes the screen session name, the current screen_index, i as inputs, read the parameters from the input file, and starts the screen session.
# It will then execute the compiled program with the input parameters in that screen and finally return the updated screen_index.
start_screen(){
    local file_path=$1
    local spb=$2
    local screen_index=$3
    local i=$4
    local init=$5
    local mvinv=$6
    local arch_sm=$7
    local prefix=""
    # If screen_index exceeds the length of the screen_names array, issue an error message and exit
    if [ $((screen_index)) -ge ${#screen_names[@]} ]; then
        echo "Error: Insufficient screen names provided. Exiting."
        exit 1
    fi
    read -r p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 <<< $(sed -n "$((i+2))p" $1)
    #Remove trailing whitespaces and new lines from p9
    if [ $spb -eq 1 ]; then
        p9=$(echo $p9 | tr -d '[:space:]')
        prefix=$p9

    elif [ $spb -eq 2 ]; then
        p11=$(echo $p11 | tr -d '[:space:]')
        prefix=$p11
    else
        p12=$(echo $p12 | tr -d '[:space:]')
        prefix=$p12
    fi
    if [ $init -eq 2 ]; then
        # Check if p13 is not empty
        if [ -z "$p13" ]; then
            p12=$(echo $p12 | tr -d '[:space:]')
        else
            p13=$(echo $p13 | tr -d '[:space:]')
        fi
    fi
    # Check if mvinv is not empty, if it is, set it to 0 (distance invariance)
    if [ -z "$mvinv" ]; then
        mvinv=0
    fi
    
    if ! screen -list | grep -q "${screen_names[$screen_index]}"; then
        # Start a new screen session and run the job
        # First creating common variables
        # Common variables
        local screen_name="${screen_names[$screen_index]}"
        local main_obj="${screen_name}_rietkerk_bjork_FKE.o"
        local executable="${screen_name}_cuFKE${spb}.out"
        local compile_flags="-DSPB=${spb} -DINIT=${init} -DMV_INVARIANCE=${mvinv} -DBARRACUDA"
        local opt_flags="-O3 -march=native -fPIE -fopenmp -std=c++20"
        local cuda_flags="-I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart -lcudadevrt -Wl,-rpath=/usr/local/cuda/lib64"
        # Define source file and object file based on init type
        local src_file
        local order_obj
        local params
        #Done by linking and compiling the object files, then running the executable
        if [ $init -ne 2 ]; then
            src_file="order_${spb}stoc_unity_rietkerk.cpp"
            order_obj="${screen_name}_order_${spb}stoc_unity.o"
            params="$p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12"
            #screen -dmS ${screen_names[$screen_index]} bash -c "cd .. ; nvcc -arch=sm_${arch_sm} -O3 -o ${screen_names[$screen_index]}_cuFKE${spb}.out rietkerk_bjork_FKE.o order_${spb}stoc_unity_rietkerk.o kernel_FKE.o -Xcompiler -fopenmp -O3 -march=native -std=c++20 -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart; ./${screen_names[$screen_index]}_cuFKE${spb}.out $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12 &> stderr_${screen_names[$screen_index]}.txt; cd $curr_dir"
            #screen -dmS ${screen_names[$screen_index]} bash -c "cd .. ; g++-14 -DSPB=${spb} -DINIT=${init} -DMV_INVARIANCE=${mvinv} -DBARRACUDA -o ${screen_names[$screen_index]}_order_${spb}stoc_unity.o -c order_${spb}stoc_unity_rietkerk.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include; g++-14 -DSPB=${spb} -DINIT=${init} -DMV_INVARIANCE=${mvinv} -DBARRACUDA -o ${screen_names[$screen_index]}_rietkerk_bjork_FKE.o -c rietkerk_bjork_FKE.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include; g++-14 -o ${screen_names[$screen_index]}_cuFKE${spb}.out ${screen_names[$screen_index]}_rietkerk_bjork_FKE.o ${screen_names[$screen_index]}_order_${spb}stoc_unity_rietkerk.o kernel_FKE.o -fopenmp -O3 -march=native -std=c++20 -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart -lcudadevrt -Wl,-rpath=/usr/local/cuda/lib64; ./${screen_names[$screen_index]}_cuFKE${spb}.out $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12 &> stderr_${screen_names[$screen_index]}.txt; cd $curr_dir"
        else
            src_file="order_${spb}stoc_burnin_rietkerk.cpp"
            order_obj="${screen_name}_order_${spb}stoc_burnin.o"
            params="$p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12 $p13"
            #screen -dmS ${screen_names[$screen_index]} bash -c "cd .. ; nvcc -arch=sm_${arch_sm} -O3 -o ${screen_names[$screen_index]}_cuFKE${spb}.out rietkerk_bjork_FKE.o order_${spb}stoc_burnin_rietkerk.o kernel_FKE.o -Xcompiler -fopenmp -O3 -march=native -std=c++20 -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart; ./${screen_names[$screen_index]}_cuFKE${spb}.out $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12 $p13 &> stderr_${screen_names[$screen_index]}.txt; cd $curr_dir"
            #screen -dmS ${screen_names[$screen_index]} bash -c "cd .. ; g++-14 -DSPB=${spb} -DINIT=${init} -DMV_INVARIANCE=${mvinv} -DBARRACUDA -o ${screen_names[$screen_index]}_order_${spb}stoc_burnin.o -c order_${spb}stoc_burnin_rietkerk.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include; g++-14 -DSPB=${spb} -DINIT=${init} -DMV_INVARIANCE=${mvinv} -DBARRACUDA -o ${screen_names[$screen_index]}_rietkerk_bjork_FKE.o -c rietkerk_bjork_FKE.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include; g++-14 -o ${screen_names[$screen_index]}_cuFKE${spb}.out ${screen_names[$screen_index]}_rietkerk_bjork_FKE.o ${screen_names[$screen_index]}_order_${spb}stoc_burnin_rietkerk.o kernel_FKE.o -fopenmp -O3 -march=native -std=c++20 -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart -lcudadevrt -Wl,-rpath=/usr/local/cuda/lib64; ./${screen_names[$screen_index]}_cuFKE${spb}.out $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12 $p13 &> stderr_${screen_names[$screen_index]}.txt; cd $curr_dir"
        fi

        # Start the screen session, compile and run the scripts
        screen -dmS $screen_name bash -c "
        cd .. || exit 1
        # Compile source files
        g++-14 $compile_flags -o $order_obj -c $src_file $opt_flags -I/usr/local/cuda/include && \
        g++-14 $compile_flags -o $main_obj -c rietkerk_bjork_FKE.cpp $opt_flags -I/usr/local/cuda/include && \
        # Link executable
        g++-14 -o $executable $main_obj $order_obj kernel_FKE.o \
        -fopenmp -O3 -march=native -std=c++20 \
        $cuda_flags && \
        
        # Run executable
        ./$executable $params &> stderr_${screen_name}_${timestamp}.txt
        
        cd $curr_dir"
        screen_index=$((screen_index+1))
    else
        #echo "Screen session ${screen_names[$screen_index]} is already running. Using the next screen name."
        screen_index=$((screen_index+1))
        # Recursively call the function with the updated screen_idx
        screen_index=$(start_screen $file_path $spb $screen_index $i $init $mvinv $arch_sm)
    fi
    echo "$screen_index"
}

# Get the number of screen sessions from the number of lines in the input file
num_screens=$(wc -l < $1)
screen_idx=0

# Read the input file and start a screen session for each job
for ((i=0; i<num_screens; i++)); do
    screen_idx=$(start_screen $1 $2 $screen_idx $i $3 $4 $5)
    #echo "Screen index: $screen_idx"
    LINE=$(sed -n "$((i+2))p" $1)
    echo "Starting screen session ${screen_names[$((screen_idx-1))]} for job $i with parameters: $LINE"
done


#g++-14 -DSPB=3 -DINIT=0 -DMV_INVARIANCE=0 -DBARRACUDA -c order_3stoc_unity_rietkerk.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include
#g++-14 -DSPB=3 -DINIT=0 -DMV_INVARIANCE=0 -DBARRACUDA -c rietkerk_bjork_FKE.cpp -O3 -march=native -fPIE -fopenmp -std=c++20 -I/usr/local/cuda/include
# nvcc -arch=sm_89 -O3 -c -o kernel_FKE.o  kernel_advdiff_FKE2D.cu -std=c++20
#g++-14 -o cornwall_cuFKE3.out rietkerk_bjork_FKE.o order_3stoc_unity_rietkerk.o kernel_FKE.o -fopenmp -O3 -march=native -std=c++20 -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart -lcudadevrt -Wl,-rpath=/usr/local/cuda/lib64