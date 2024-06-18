#! /usr/bin/bash

num="8"

for i in {1..5}
do
    dir="/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/${num}k_particles/exp${i}"
    cd $dir
    mkdir analysis_results
    mkdir hole_analysis_results
    cd "/home/felipe/HYDROGELS/Hydrogel-Simulation/analysis"
    python3 analysis.py "${dir}/dump_shearing_${num}k_${i}.lammpstrj"
done
