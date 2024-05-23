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
    import json

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
    
    dir = dumpdir.split('/')
    dir[-1] = "analysis_results"
    savedir = '/'.join(dir)

    print('Esribiendo resultados en archivos json')
    with open (f'{savedir}/analysis_bonds_v_t.json','w') as f:
        json.dump(bonds_v_t,f)
    with open (f'{savedir}/analysis_clusternum_v_t.json','w') as f:
        json.dump(clusternum_v_t,f)
    with open (f'{savedir}/analysis_bigclustersz_v_t.json','w') as f:
        json.dump(bigclustersz_v_t,f)
    x = x.tolist()
    with open (f'{savedir}/analysis_x1.json','w') as f:
        json.dump(x,f)
    print('FIN')

    #return bonds_v_t,clusternum_v_t,bigclustersz_v_t,x

def non_affine_sq_disp(dumpdir):
    from ovito.io import import_file
    from ovito.modifiers import SelectTypeModifier,AtomicStrainModifier, DeleteSelectedModifier, UnwrapTrajectoriesModifier
    from ovito.pipeline import ReferenceConfigurationModifier
    import numpy as np
    import json

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

    dir = dumpdir.split('/')
    dir[-1] = "analysis_results"
    savedir = '/'.join(dir)

    print('Esribiendo resultados en archivos json')
    with open (f'{savedir}/analysis_Dsq_v_t.json','w') as f:
        json.dump(Dsq_v_t,f)
    x = x.tolist()
    with open (f'{savedir}/analysis_x2.json','w') as f:
        json.dump(x,f)
    print('FIN')

    #return Dsq_v_t,x

def main():
    import sys

    dumpdir = sys.argv[1]

    print('Iniciando analisis')
    count_bonds_and_clusters(dumpdir)
    non_affine_sq_disp(dumpdir)
    print('FIN')

if __name__=="__main__":
    main()














