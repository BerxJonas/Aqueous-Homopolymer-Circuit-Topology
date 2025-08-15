#!/bin/bash

# Define the values for the fourth argument
last_args=(240 250 260 270 280 290 298 310 320 330 340 350 360)

# Loop over each value of the last argument
for last_arg in "${last_args[@]}"
do
    # Loop over third argument from 0.1 to 1.0 in steps of 0.1
    for i in $(seq 0.1 0.1 1.0)
    do
        #echo "Running with arguments: 4 30 $i $last_arg"
        echo "Running with arguments: 1 30 $i $last_arg"
        #./circuit_topology_analysis.exe 4 30 $i $last_arg
        ./circuit_topology_analysis.exe 1 30 $i $last_arg
    done
done

