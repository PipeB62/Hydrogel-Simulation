'''
En este archivo se encuentran las funciones count_bonds_and_clusters() y non_affine_sq_disp()
count_bonds_and_clusters() calcula el numero de bonds, numero de clusters y tamaño del cluster más grande con respecto al tiempo. 
non_affine_sq_disp() calcula el non affine squared displacement del sistema con respecto al tiempo. Toma como referencia el frame 0. 

El main() da la opcion de especificar los directorios de dumps arbitrarios para realizar este analisis y hacer graficas. 
'''

def write_json(savedir,name,data):
    import json 
    with open (f'{savedir}/{name}.json','w') as f:
        json.dump(data,f)

def load_data(data):
    import numpy as np
    a = np.loadtxt(data)
    return a.tolist()


def count_bonds_and_clusters(dumpdir,savedir,savelm):
    from ovito.io import import_file
    from ovito.modifiers import ClusterAnalysisModifier, CreateBondsModifier
    import numpy as np
    import json

    print("Inicio bonds y clusters")

    node = import_file(dumpdir)

    #Obtener numero de iteraciones
    iter_num = node.source.num_frames

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
    Rlim2 = 0.52 #para bonds
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
    
    #dir = dumpdir.split('/')
    #dir[-1] = "analysis_results"
    #savedir = '/'.join(dir)

    print('Esribiendo resultados en archivos json')
    write_json(savedir,f"bondnum_{savelm}",bonds_v_t)
    write_json(savedir,f"clusternum_{savelm}",clusternum_v_t)
    write_json(savedir,f"bigclustersz_{savelm}",bigclustersz_v_t)
    print('FIN')

    #return bonds_v_t,clusternum_v_t,bigclustersz_v_t,x

def non_affine_sq_disp(dumpdir,savedir):
    from ovito.io import import_file
    from ovito.modifiers import SelectTypeModifier,AtomicStrainModifier, DeleteSelectedModifier, UnwrapTrajectoriesModifier
    from ovito.pipeline import ReferenceConfigurationModifier
    import numpy as np
    import json

    print("Inicio non-affine sq displacement")

    node = import_file(dumpdir)

    #Obtener numero de iteraciones
    iter_num = node.source.num_frames

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

    #dir = dumpdir.split('/')
    #dir[-1] = "analysis_results"
    #savedir = '/'.join(dir)

    print('Esribiendo resultados en archivos json')
    write_json(savedir,"Dsq",Dsq_v_t)
    print('FIN')

    #return Dsq_v_t,x

def strain(dumpdir,savedir,savelm):
    import json

    print("Inicio strain")

    #obtener numero de frames
    frame_num = 0
    with open(dumpdir,"r") as dump:
        for line in dump.readlines():
            if line == "ITEM: TIMESTEP\n":
                frame_num+=1
    
    dump = open(dumpdir,"r")
    strain = []
    jumps=0
    for frame in range(frame_num):
        for i in range(3):
            dump.readline()
        atom_num = int(dump.readline())
        dump.readline()
        xmin,xmax,xy = [float(x) for x in dump.readline().split()]
        dump.readline()
        zmin,zmax,xz = [float(x) for x in dump.readline().split()]
        dump.readline()
        for i in range(atom_num):
            dump.readline()
        L = zmax-zmin
        s = xy/L+jumps   
        if len(strain)>0:
            sa = strain[-1]
        else:
            sa = 0    
        if s-sa<0:
            jumps+=1
            s=s+1
        strain.append(s)
    dump.close()

    #dir = dumpdir.split('/')
    #dir[-1] = "analysis_results"
    #savedir = '/'.join(dir)

    print('Esribiendo resultados en archivos json')
    write_json(savedir,f"strain_{savelm}",strain)
    print("FIN")

def stress(presdir,savedir):

    print("Inicio stress")

    stress = [-row[4] for row in load_data(presdir)]

    print("Escribiendo resultados en archivos json")
    write_json(savedir,f"stress",stress)
    print("FIN")

def pot_ene_formation(pot_ene_dir,savedir,savelm):
    print("Inicio pot_ene_formation")

    pe = [row[1] for row in load_data(pot_ene_dir)]
    frames = [row[0] for row in load_data(pot_ene_dir)]

    print("Escribiendo resultados en archivos json")
    write_json(savedir,f"pot_ene_formation_{savelm}",pe)
    write_json(savedir,f"pot_ene_formation_frames_{savelm}",frames)
    print("FIN")

def main():
    import sys

    #dir = sys.argv[1]
    #dumpdir = sys.argv[1]
    presdir = sys.argv[1]
    savedir = sys.argv[2]
    #pot_ene = sys.argv[4]

    #dumpdir = "/".join([dir,dump])

    #print('Iniciando analisis')
    #count_bonds_and_clusters(dumpdir,savedir,savelm)
    #non_affine_sq_disp(dumpdir,savedir)
    #strain(dumpdir,savedir,savelm)
    stress(presdir,savedir)
    #pot_ene_formation(dir,pot_ene)
    print('FIN')

if __name__=="__main__":
    main()














