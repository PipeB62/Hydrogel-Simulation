'''
En este archivo se encuentran las funciones count_bonds_and_clusters() y non_affine_sq_disp()
count_bonds_and_clusters() calcula el numero de bonds, numero de clusters y tamaño del cluster más grande con respecto al tiempo. 
non_affine_sq_disp() calcula el non affine squared displacement del sistema con respecto al tiempo. Toma como referencia el frame 0. 

El main() da la opcion de especificar los directorios de dumps arbitrarios para realizar este analisis y hacer graficas. 
'''
def count_bonds_and_clusters(dumpdir):
    from ovito.io import import_file
    from ovito.modifiers import ClusterAnalysisModifier, CreateBondsModifier
    import numpy as np

    node = import_file(dumpdir)

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

    return bonds_v_t,clusternum_v_t,bigclustersz_v_t,x

def non_affine_sq_disp(dumpdir):
    from ovito.io import import_file
    from ovito.modifiers import SelectTypeModifier,AtomicStrainModifier, DeleteSelectedModifier, UnwrapTrajectoriesModifier
    from ovito.pipeline import ReferenceConfigurationModifier
    import numpy as np

    node = import_file(dumpdir)

    #Obtener numero de iteraciones
    iter_num = node.source.num_frames

    x = np.arange(iter_num-1)

    #Seleccionar patches
    selectpatches = SelectTypeModifier(operate_on = "particles",
                                        property = "Particle Type",
                                        types = {2,4})
    node.modifiers.append(selectpatches)

    #Eliminar patches
    deletepatches = DeleteSelectedModifier()
    node.modifiers.append(deletepatches)

    #Unwrap trajectories
    unwrap = UnwrapTrajectoriesModifier()
    node.modifiers.append(unwrap)

    #Non affine square displacement
    atomic_strain1 = AtomicStrainModifier(output_nonaffine_squared_displacements=True,
                                        affine_mapping = ReferenceConfigurationModifier.AffineMapping.ToReference,
                                        use_frame_offset = False,
                                        minimum_image_convention = False)
    
    node.modifiers.append(atomic_strain1)

    Dsq_v_t = []
    for i in range(1,iter_num):
        data=node.compute(i)
        per_particle_dsq = data.particles["Nonaffine Squared Displacement"]
        Dsq_v_t.append(sum(per_particle_dsq)/len(per_particle_dsq))

    return Dsq_v_t,x

def main():
    #Esto es para que ovito y matplotlib funcionen correctamente
    import os
    os.environ['OVITO_GUI_MODE'] = '1'

    import numpy as np
    import matplotlib.pyplot as plt 
    import sys

    input_filedir = sys.argv[1:]

    colors = ['red','green','blue']
    markers = ['|','_','.']

    fig1,ax1 = plt.subplots()
    ax1.set_title('Bonds')
    ax1.legend(input_filedir)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Number of bonds')

    fig2,ax2 = plt.subplots()
    ax2.set_title('Number of clusters')
    ax2.legend(input_filedir)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Number of clusters')

    fig3,ax3 = plt.subplots()
    ax3.set_title('Biggest cluster size')
    ax3.legend(input_filedir)
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Number of atoms in biggest cluster')

    fig4,ax4 = plt.subplots()
    ax4.set_title('Non affine sq displacement')
    ax4.legend(input_filedir)
    ax4.set_xlabel('Time')
    ax4.set_ylabel('$D^2$')  

    for k,file in enumerate(input_filedir):

        bonds_v_t,clusternum_v_t,bigclustersz_v_t,x1 = count_bonds_and_clusters(file)
        D_sq_v_t,x2 = non_affine_sq_disp(file)
        
        ax1.plot(x1,bonds_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')
        ax2.plot(x1,clusternum_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')
        ax3.plot(x1,bigclustersz_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')
        ax4.plot(x2,D_sq_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')

    plt.show()

if __name__=="__main__":
    main()














