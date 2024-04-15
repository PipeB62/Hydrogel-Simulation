#Esto es para que ovito y matplotlib funcionen correctamente
import os
os.environ['OVITO_GUI_MODE'] = '1'

#Importar librerias
from ovito.io import import_file
from ovito.modifiers import ClusterAnalysisModifier, CreateBondsModifier
import numpy as np
import matplotlib.pyplot as plt 

#Version 
version = input('Version: ')

#Importar 1 dump de formacion 
input_filepath= "/home/pipe/lammps/code/HYDROGELS/selfassembly2/dump_v"+version+".lammpstrj"

node = import_file(input_filepath)

#Obtener numero de iteraciones
iter_num = node.source.num_frames

x = np.arange(1,iter_num)

#Modificador para crear bonds patch-patch
Rlim = 0.5 #Cutoff para crear bonds patch-patch
create_bonds1 = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise) 
create_bonds1.set_pairwise_cutoff(2,2,Rlim)
create_bonds1.set_pairwise_cutoff(2,4,Rlim)
node.modifiers.append(create_bonds1) #Agregar modificador para crear bonds

#Calcular numero de bonds en cada iteracion
bonds_v_t = [] 
for i in range(1,iter_num):
    data=node.compute(i)
    bonds_v_t.append(data.particles.bonds.count) 


#Modificador para crear bonds centro-patch
Rlim2 = 0.5
create_bonds2 = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise) 
create_bonds2.set_pairwise_cutoff(1,2,Rlim2)
create_bonds2.set_pairwise_cutoff(3,4,Rlim2)
node.modifiers.append(create_bonds2) #Agregar modificador para crear bonds

#Modificador para numero de clusters
clusters = ClusterAnalysisModifier(neighbor_mode=ClusterAnalysisModifier.NeighborMode.Bonding,
                                sort_by_size = True)
node.modifiers.append(clusters)

#Calcular numero de clusters en cada iteracion
clusternum_v_t = []
bigclustersz_v_t = []
for i in range(1,iter_num):
    data=node.compute(i)
    clusternum_v_t.append(data.attributes["ClusterAnalysis.cluster_count"])
    bigclustersz_v_t.append(data.attributes["ClusterAnalysis.largest_size"])

pot_eng_mat = np.loadtxt("/home/pipe/lammps/code/HYDROGELS/selfassembly2/pot_ene_v"+version+".dat")
pot_eng = pot_eng_mat[:,1]

plt.figure()
plt.plot(x,bonds_v_t)
plt.title('Bonds')
plt.xlabel('Time')
plt.ylabel('Number of bonds')

plt.figure()
plt.plot(x,clusternum_v_t)
plt.title('Number of clusters')
plt.title('Clusters')
plt.xlabel('Time')
plt.ylabel('Number of clusters')

plt.figure()
plt.plot(x,bigclustersz_v_t)
plt.title('Biggest cluster size')
plt.xlabel('Time')
plt.ylabel('Number of atoms in biggest cluster')

'''
plt.figure()
plt.plot(x[1:],pot_eng)
plt.title('Potential Energy')
plt.xlabel('Time')
plt.ylabel('Potential Energy')
'''
plt.show()