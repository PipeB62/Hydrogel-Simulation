#! /usr/bin/bash

#Dar como input L (longitud de caja), numero de xl, numero de monomeros, lambda (tres cuerpos), sigma (patches), l (longitud centro-patch), shear_rate, directorio, y apellido del archivo
#Ejemplo:
#bash create_simulation.sh 51.08 1120 14480 15.0 0.4 0.3612 0.01 directorio apellido

L=$1
xl_num=$2
mon_num=$3
lam=$4
sigma=$5
l=$6
shear_rate=$7
directorio=$8
apellido=$9

mkdir "$directorio/input_data"

cd input_data
python3 patchy_particles.py $l $directorio
python3 patches.py $sigma $directorio
python3 3b.py $lam $sigma $directorio

cd ../self_assembly
python3 create_input_selfassembly_v2.py $L $xl_num $mon_num $lam $sigma $l $directorio $apellido

cd ../shearing
python3 create_input_shearing.py $L $shear_rate $directorio $apellido

touch "$directorio/PARAMETERS.txt"
printf "L = $L \nxl_num = $xl_num \nmon_num = $mon_num \nlambda (3b) = $lam \nsigma (patches) = $sigma \nl (c-p) = $l \nshear rate = $shear_rate" > "$directorio/PARAMETERS.txt"
