#! /usr/bin/bash

directory=$1 #Directorio donde se van a generar los archivos necesarios para la simulacion

#Parametros del experimento
L=51.08
N_xl=1120
N_mon=14880
lam=15.0
sigma=0.4
l_cp=0.3612
shear_rate=(0.01 0.001 0.0001)

#Crear todas las carpetas
mkdir $directory/FullExperiment1 
mkdir $directory/FullExperiment1/Formation 
mkdir $directory/FullExperiment1/dGamma1
mkdir $directory/FullExperiment1/dGamma2
mkdir $directory/FullExperiment1/dGamma3
mkdir $directory/FullExperiment1/input_data
for i in {1..10}
do 
    mkdir $directory/FullExperiment1/Formation/exp$i
    mkdir $directory/FullExperiment1/dGamma1/exp$i
    mkdir $directory/FullExperiment1/dGamma2/exp$i
    mkdir $directory/FullExperiment1/dGamma3/exp$i
done

#Crear input_data
cd input_data
python3 patchy_particles.py $l_cp $directory/FullExperiment1/input_data
python3 patches.py $sigma $directory/FullExperiment1/input_data
python3 3b.py $lam $sigma $directory/FullExperiment1/input_data

#Crear inputs para la formacion
#L,xl_num,mon_num,l_cp,savedir,inputfilesdir,last_name
cd ../self_assembly
for i in {1..10}
do
    python3 create_input_selfassembly_v2.py $L $N_xl $N_mon $l_cp $directory/FullExperiment1/Formation/exp$i $directory/FullExperiment1/input_data $i
done

#Crear inputs para el shear 
#L,shear_rate,savedir,system_formation,inputfilesdir,last_name
cd ../shearing
for k in {1..3} 
do
    for i in {1..10}
    do
        python3 create_input_shearing.py $L ${shear_rate[$k-1]} $directory/FullExperiment1/dGamma$k/exp$i $directory/FullExperiment1/Formation/exp$i/system_formation_$i.data $directory/FullExperiment1/input_data dGamma${k}_${i}
    done
done

echo "L = $L" >> $directory/FullExperiment1/PARAMETERS.txt 
echo "xl_num = $N_xl" >> $directory/FullExperiment1/PARAMETERS.txt
echo "mon_num = $N_mon" >> $directory/FullExperiment1/PARAMETERS.txt
echo "lambda (3b) = $lam" >> $directory/FullExperiment1/PARAMETERS.txt
echo "sigma (patches) = $sigma" >> $directory/FullExperiment1/PARAMETERS.txt
echo "l (c-p) = $l_cp" >> $directory/FullExperiment1/PARAMETERS.txt 
echo "shear rate = ${shear_rate[@]}" >> $directory/FullExperiment1/PARAMETERS.txt