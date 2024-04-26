#!/bin/bash

# file       test.sh.cpp
# author    Pavel Kratochvil
#           xkrato61@vutbr.cz
# version   2023
# date      26 April  2020

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <file_path> <iteration_count>"
    exit 1
fi

# get the script arguments
file_path=$1
iteration_count=$2

# Check if the iteration count is a non-negative integer
if ! [[ "$iteration_count" =~ ^[0-9]+$ ]]; then
    echo "Error: Iteration count must be a non-negative integer."
    exit 1
fi

# Check if the file exists
if [ ! -f "$file_path" ]; then
    echo "Error: File '$file_path' not found."
    exit 1
fi

# get the physical core count on mac or linux machine
os=$(uname)
if [ "$os" = "Darwin" ]; then
    # macOS
    physical_cores=$(sysctl -n hw.physicalcpu)
elif [ "$os" = "Linux" ]; then
    # Linux
    cores_per_socket=$(lscpu | grep 'Core(s) per socket:' | awk '{print $4}')
    n_sockets=$(lscpu | grep 'Socket(s):' | awk '{print $2}')
    physical_cores=$((cores_per_socket * n_sockets))
else
    echo "Unsupported operating system: $os"
    exit 1
fi

# get the input file lines
input_file_lines=$(wc -l < "$file_path" | awk '{print $1}')

# choose the minimum for optimum utilization
np=$input_file_lines
if [ "$physical_cores" -lt "$input_file_lines" ]; then
    np=$physical_cores
fi

# Compile and run the program
mpicxx -std=c++20 "life.cpp" -o life && mpiexec -np $np ./life "$file_path" $iteration_count