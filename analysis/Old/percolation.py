
#Importar librerias
from ovito.io import import_file
from ovito.modifiers import ClusterAnalysisModifier, CreateBondsModifier, ReplicateModifier
import numpy as np

version = input('Version: ')

input_filepath="/home/pipe/lammps/code/HYDROGELS/selfassembly2/system_v"+version+".data"

node = import_file(input_filepath)

node2 = import_file(input_filepath)
node2.modifiers.append(ReplicateModifier(adjust_box = True, num_x = 3, num_y = 3, num_z = 3))

#Obtener numero de iteraciones
iter_num = node.source.num_frames

#Modificador para crear bonds patch-patch
Rlim = 0.5 #Cutoff para crear bonds patch-patch
create_bonds1 = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise) 
create_bonds1.set_pairwise_cutoff(2,2,Rlim)
create_bonds1.set_pairwise_cutoff(2,4,Rlim)
node.modifiers.append(create_bonds1) #Agregar modificador para crear bonds
node2.modifiers.append(create_bonds1) #Agregar modificador para crear bonds

#Modificador para crear bonds centro-patch
Rlim2 = 0.45
create_bonds2 = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise) 
create_bonds2.set_pairwise_cutoff(1,2,Rlim2)
create_bonds2.set_pairwise_cutoff(3,4,Rlim2)
node.modifiers.append(create_bonds2) #Agregar modificador para crear bonds
node2.modifiers.append(create_bonds2) #Agregar modificador para crear bonds

#Modificador para numero de clusters
clusters = ClusterAnalysisModifier(neighbor_mode=ClusterAnalysisModifier.NeighborMode.Bonding,
                                   sort_by_size = True)
node.modifiers.append(clusters)
node2.modifiers.append(clusters)

print('Cluster size:',node.compute().attributes['ClusterAnalysis.largest_size'])
print('Cluster number:',node.compute().attributes['ClusterAnalysis.cluster_count'])

print('Big cluster size:',node2.compute().attributes['ClusterAnalysis.largest_size'])
print('Big cluster number:',node2.compute().attributes['ClusterAnalysis.cluster_count'])




