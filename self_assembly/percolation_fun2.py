def read_system_data(datafile):
    import numpy as np

    #Inicializar variables
    N_atoms = 0
    L = 0
    centers_coords = [] #x,y,z
    centers_ids = [] #atomid,moleculeid,type
    patches_coords = [] #x,y,z
    patches_ids = [] #atomid,moleculeid,type

    with open(datafile, "r") as data:
        for line in data.readlines():
            aa = line.split()
            if len(aa)==2 and aa[1]=="atoms": #Leer numero de atomos
                N_atoms = int(aa[0])
            elif len(aa)==4 and aa[2] =="xlo": #Leer longitud de caja
                L = float(aa[1])-float(aa[0])
            elif len(aa)==10: #Leer info de atomos
                if aa[2]=="1" or aa[2]=="3": #Centros
                    centers_ids.append([int(k) for k in aa[:3]])
                    centers_coords.append([float(k) for k in aa[4:7]])
                else: #Patches
                    patches_ids.append([int(k) for k in aa[:3]])
                    patches_coords.append([float(k) for k in aa[4:7]])
            if len(centers_coords)+len(patches_coords)==N_atoms and N_atoms>0: #Terminar cuando se lean todos los atomos
                break
    
    #Convertir a numpy array
    centers_coords = np.array(centers_coords,dtype=np.float64)
    centers_ids = np.array(centers_ids,dtype=np.int64)
    patches_coords = np.array(patches_coords,dtype=np.float64)
    patches_ids = np.array(patches_ids,dtype=np.int64)

    return N_atoms,L,centers_ids,centers_coords,patches_ids,patches_coords

def periodic_data(ids,coords,L):
    import numpy as np

    #Crear vectores de traslacion
    tr = []
    for q in [-1,0,1]:
        for p in [-1,0,1]:
            for o in [-1,0,1]:
                tr.append(np.array([L*o,L*p,L*q],dtype=np.float64))
    
    N_atoms = np.shape(coords)[0] #Obtener numero de atomos
    ids_periodic = np.zeros((len(tr)*N_atoms,3),dtype=np.int64) #Inicializar matriz para id periodicos
    coords_periodic = np.zeros((len(tr)*N_atoms,3),dtype=np.float64) #Inicializar matriz para coordenadas periodicas

    #Llenar matrices con datos periodicos
    for i in range(len(tr)):
        ids_periodic[N_atoms*i:N_atoms*(i+1),:] = ids 
        coords_periodic[N_atoms*i:N_atoms*(i+1),:] = coords + tr[i]

    return ids_periodic,coords_periodic 

def get_neighbors(c_id,centers_ids,patches_ids, patches_ids_periodic, patches_neigh_indexes):
    import numpy as np

    mol_id = centers_ids[centers_ids[:,0]==c_id][:,1] #Encontrar mol id del centro
    
    mol_patches_indexes = np.where(patches_ids[:,1]==mol_id)[0].tolist() #Indices de los patches correspondientes a molecula actual

    #Indices de los patches periodicos vecinos
    p_neighs_ix=[] 
    for p_ix in mol_patches_indexes:
        p_neighs_ix.extend(patches_neigh_indexes[p_ix])

    #Molecule ids de vecinos
    neigh_mol=[]
    for ii in p_neighs_ix:
        cpnmol = patches_ids_periodic[ii,1]
        if cpnmol != mol_id:
            neigh_mol.append(cpnmol)

    #Id del centro de los vecinos
    neighs = []
    for ii in neigh_mol:
        neighs.append(centers_ids[centers_ids[:,1]==ii][:,0][0])
    
    return neighs
        
def percolation(directory,sigma):
    import numpy as np
    from scipy.spatial import KDTree

    datafile=f"{directory}/system_check.data"

    _,L,centers_ids,_,patches_ids,patches_coords = read_system_data(datafile) #Leer info del system.data
    #N_patches = np.shape(patches_ids)[0] #Obtener numero de patches
    N_centers = np.shape(centers_ids)[0] #Obtener numero de centros (Numero de particulas)

    patches_ids_periodic,patches_coords_periodic = periodic_data(patches_ids,patches_coords,L) #Obtener datos periodicos

    patches_tree = KDTree(patches_coords) #Crear KDTree con datos de patches
    patches_tree_periodic = KDTree(patches_coords_periodic) #Crear KDTree con datos periodicos

    a = 1.5
    r_c = a*sigma
    patches_neigh_indexes = patches_tree.query_ball_tree(patches_tree_periodic,r_c) #Obtener indices de vecinos para cada patch

    centerlist = centers_ids[:,0].tolist() #Obtener lista de todos los centros
    clustersz = []  #Inicializar lista para almacenar tamaños de cluster
    
    #BFS para encontrar clusters
    while len(centerlist)>0:
        queue = [centerlist[0]]
        visited = [centerlist[0]]
        centerlist.pop(0)
        while len(queue)>0:
            c_center = queue[0]
            neighs = get_neighbors(c_center,centers_ids,patches_ids, patches_ids_periodic, patches_neigh_indexes)
            for nn in neighs:
                if nn not in visited:
                    queue.append(nn)
                    visited.append(nn)
                    centerlist.remove(nn)
            queue.pop(0)
            #print(len(queue))
        clustersz.append(len(visited))

    #clusternum = len(clustersz)
    max_cluster_sz = max(clustersz)
    percentage = max_cluster_sz/N_centers

    return percentage

#percolation("/media/felipe/Files/Hydrogel_sim_experiments/SelfAssemby/exp2",0.2)





