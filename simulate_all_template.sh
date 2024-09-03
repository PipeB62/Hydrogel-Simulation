dir=$(pwd)

for i in {1..10}
do
	cd $dir/exp$i
	echo "#!/bin/bash" >> simulate.sge
	echo "#$ -cwd" >> simulate.sge
	echo "#$ -j y" >> simulate.sge
	echo "#$ -S /bin/bash" >> simulate.sge
	echo "#$ -N Name_$i" >> simulate.sge
	echo ". /etc/profile.d/modules.sh" >> simulate.sge
	echo "module load gcc/8.3.0" >> simulate.sge
	echo "module load openmpi/gcc/64/1.10.1" >> simulate.sge
	echo "#$ -pe mpich 4" >> simulate.sge
	echo "mpirun -np 4 /mnt/MD1200B/cferreiro/fbenavides/lammps-2Aug2023/src/lmp_mpi -in input_name_$i.lammps" >> simulate.sge
	
	qsub simulate.sge
	sleep 2 

done