#! /usr/bin/bash

num="8"

for i in {1..5}
do
    dir="/media/felipe/Hydrogel_sim_experiments/SizeExperiments/${num}k_particles/exp${i}"
    cd dir
    mkdir analysis_results
    mkdir hole_analysis_results