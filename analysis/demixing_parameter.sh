#! /usr/bin/bash

num="16"

N_ls=("2" "4" "5" "7" "10" )
for i in {1..5}
do
    for N_l in ${N_ls[@]}
    do
        dir="/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/${num}k_particles/exp${i}"
        cd $dir
        mkdir analysis_results
        cd "/home/felipe/HYDROGELS/Hydrogel-Simulation/analysis"
        julia demixing_parameter.jl "${dir}/dump_shearing_${num}k_${i}.lammpstrj" $N_l
    done
done