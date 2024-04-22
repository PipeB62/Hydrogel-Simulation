#Esto es para que ovito y matplotlib funcionen correctamente
import os
os.environ['OVITO_GUI_MODE'] = '1'

#Importar librerias
from ovito.io import import_file
from ovito.modifiers import ClusterAnalysisModifier, CreateBondsModifier
import numpy as np
import matplotlib.pyplot as plt 

input_filefolder = input("Folder directory: ")
input_filenames = []
legend = []
for i in range(3):
    file = input("Filename (type end to stop): ")
    if file == "end":
        break
    else:
        input_filenames.append(file)
        input('Legend: ')

colors = ['red','green','blue']
markers = ['|','_','.']

plt.figure('bonds_v_t')
plt.title('Bonds')

plt.figure('clusternum_v_t')
plt.title('Number of clusters')

plt.figure('clustersize_v_t')
plt.title('Biggest cluster size')

for k,file in enumerate(input_filenames):
    node = import_file(input_filefolder+"/"+file)

    #Obtener numero de iteraciones
    iter_num = node.source.num_frames

    x = np.arange(iter_num-1)

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
    #Rlim2 = 0.37 #para cuerpos rigidos
    Rlim2 = 0.5 #para bonds
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

    plt.figure('bonds_v_t')
    plt.plot(x,bonds_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')

    plt.figure('clusternum_v_t')
    plt.plot(x,clusternum_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')

    plt.figure('clustersize_v_t')
    plt.plot(x,bigclustersz_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')

plt.figure('bonds_v_t')
plt.legend(legend)
plt.xlabel('Time')
plt.ylabel('Number of bonds')

plt.figure('clusternum_v_t')
plt.legend(legend)
plt.xlabel('Time')
plt.ylabel('Number of clusters')

plt.figure('clustersize_v_t')
plt.legend(legend)
plt.xlabel('Time')
plt.ylabel('Number of atoms in biggest cluster')

plt.show()












