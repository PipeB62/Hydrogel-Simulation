#! /usr/bin/bash

num="16"

for i in {1..5}
do

    dir="/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/${num}k_particles/exp${i}"
    cd $dir
    mkdir analysis_results
    cd "/home/felipe/HYDROGELS/Hydrogel-Simulation/analysis"
    julia fractal_dimension.jl "${dir}/dump_shearing_${num}k_${i}.lammpstrj" $N_l

done