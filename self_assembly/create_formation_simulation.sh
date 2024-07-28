#! /usr/bin/bash

#Dar como input L (longitud de caja), numero de xl, numero de monomeros, lambda (tres cuerpos), sigma (patches), l (longitud centro-patch), directorio, y apellido del archivo

L=$1
xl_num=$2
mon_num=$3
lam=$4
sigma=$5
l=$6
directorio=$7
apellido=$8

mkdir "$directorio/input_data"
python3 create_input_selfassembly2.py $L $xl_num $mon_num $lam $sigma $l $directorio $apellido
cd ../input_data
python3 patchy_particles.py $l $directorio
python3 patches.py $sigma $directorio
python3 3b.py $lam $sigma $directorio

touch "$directorio/PARAMETERS.txt"
printf "L = $L \n xl_num = $xl_num \n mon_num = $mon_num \n lambda (3b) = $lam \n sigma (patches) = $sigma \n l (c-p) = $l" > "$directorio/PARAMETERS.txt"
