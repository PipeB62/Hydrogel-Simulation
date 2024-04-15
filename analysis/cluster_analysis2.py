import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
from scipy.spatial import KDTree

#Leer system.data de lammps
with open("system.data", "r") as search:
    for num, line in enumerate(search, 1):
        if 'atoms' in line:
            atomnum = int(line.split()[0])
        if 'Atoms' in line:
            atoms_line = num
        elif 'Velocities' in line:
            velocities_line = num
        elif 'Bonds' in line: 
            bonds_line = num
        elif 'Angles' in line:
            angles_line = num

#Leer dump.lammpstrj de lammps para extraer linea de inicio de cada iteracion
with open("dump.lammpstrj", "r") as search:
    tablestart = []
    boxboundstart=[]
    for num, line in enumerate(search, 1):
        if 'ITEM: ATOMS id type x y z' in line:
            tablestart.append(num)
        if 'ITEM: BOX BOUNDS xy xz yz pp pp pp' in line:
            boxboundstart.append(num)

#Importar tabla con info de molecule_ids. Esa info es fija, por lo que no va dentro de los ciclos
atomsfull = pd.read_csv("system.data",sep=' ',skiprows=atoms_line,header=None,nrows=atomnum).to_numpy() 

iter_num = len(tablestart) #Numero de frames en el dump
clusternum_v_t = [] #Inicalizar lista para almacenar numero de clusters en cada frame
clustersize_v_t = [] #Inicalizar lista para almacenar tamaño del cluster principal en cada frame
bonds_v_t = [] #Inicializar lista para almacenar numero de bonds en cada frame
for iter in range(iter_num):
    
    atoms = pd.read_csv("dump.lammpstrj",sep=' ',skiprows=int(tablestart[iter]),header=None,nrows=atomnum).to_numpy() #Importar lista de atomos en iteracion actual
    boxinfo = pd.read_csv("dump.lammpstrj",sep=' ',skiprows=int(boxboundstart[iter]),header=None,nrows=3).to_numpy() #Importar informacion de caja
    xls = atoms[np.logical_or(atoms[:,1]==2,atoms[:,1]==4)] #Extraer info de cross-linkers
    xlsnum = np.size(xls[:,0]) #Obtener numero de cross-linkers

    #Extraer info de la boundig box actual. Cambia cada iteracion por la deformacion
    xlo_bound = boxinfo[0,0] #Extraer xlo_bound
    xhi_bound = boxinfo[0,1] #Extraer xhi_bound
    xy = boxinfo[0,2] #Extraer xy
    ylo_bound = boxinfo[1,0] #Extraer ylo_bound
    yhi_bound = boxinfo[1,1] #Extraer yhi_bound
    xz = boxinfo[1,2] #Extraer xz
    zlo_bound = boxinfo[2,0] #Extraer zlo_bound
    zhi_bound = boxinfo[2,1] #Extraer zhi_bound
    yz = boxinfo[2,2] #Extraer yz
    #Obtener parametros de trclinic box a partir de los parametros de la bounding box
    xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
    xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
    ylo = ylo_bound - min(0.0,yz)
    yhi = yhi_bound - max(0.0,yz)
    zlo = zlo_bound
    zhi = zhi_bound

    #Crear vectores eje y almacenar desplazamientos en una misma matriz
    a = np.array([xhi-xlo,0,0]) 
    b = np.array([xy,yhi-ylo,0])
    c = np.array([xz,yz,zhi-zlo])
    periodic_disp = np.zeros([26,3])
    periodic_disp[0,:] = a
    periodic_disp[1,:] = -a
    periodic_disp[2,:] = b
    periodic_disp[3,:] = -b
    periodic_disp[4,:] = c
    periodic_disp[5,:] = -c
    periodic_disp[6,:] = a + b
    periodic_disp[7,:] = a - b
    periodic_disp[8,:] = a + c
    periodic_disp[9,:] = a - c
    periodic_disp[10,:] = -a + b
    periodic_disp[11,:] = -a - b
    periodic_disp[12,:] = -a + c
    periodic_disp[13,:] = -a - c
    periodic_disp[14,:] = c + b
    periodic_disp[15,:] = c - b
    periodic_disp[16,:] = -c + b
    periodic_disp[17,:] = -c - b
    periodic_disp[18,:] = a + b + c
    periodic_disp[19,:] = -a + b + c
    periodic_disp[20,:] = a - b + c
    periodic_disp[21,:] = a + b - c
    periodic_disp[22,:] = -a - b + c
    periodic_disp[23,:] = -a + b - c
    periodic_disp[24,:] = a - b - c
    periodic_disp[25,:] = -a - b - c

    Rlim = 0.5 #Limite de distancia para definir un bond

    #Crear lista de xls incluyendo copias periodicas
    xls_periodic = np.zeros([27*xlsnum,5]) #Inicializar matriz
    xls_periodic[0:xlsnum,:] = xls #Copia central
    for j in range(26): # 26 copias para cada cara, arista, esquina con su respectivo desplazamiento
        xls_periodic[(j+1)*xlsnum:(j+2)*xlsnum,0] = xls[:,0]
        xls_periodic[(j+1)*xlsnum:(j+2)*xlsnum,1] = xls[:,1]
        xls_periodic[(j+1)*xlsnum:(j+2)*xlsnum,2:] = xls[:,2:] + periodic_disp[j,:]

    #Obtener vecinos utilizando KDTree
    tree1 = KDTree(xls[:,2:]) #Crear KDTree de cross linkers (solo posiciones reales)
    tree2 = KDTree(xls_periodic[:,2:]) #Crear KDTree de cross linkers con copias periodicas
    neigh_indexes = tree1.query_ball_tree(tree2,Rlim) #Encontrar vecinos de cada cross linker
    

    neighslist = [] #inicializar lista de bonds
    for j in range(xlsnum): #Iterar para todos los cross linkers
        current_id = xls[j,0].tolist() #Obtener id de cross linker actual
        vecinos = xls_periodic[neigh_indexes[j],0].tolist() #Extraer ids de vecinos
        if current_id in vecinos: #Si el id actual está en la lista de vecinos
            vecinos.remove(current_id) #Quitar id actual de lista de vecinos
        neighslist.append(vecinos)
    
    xlsorder = xls[:,0].tolist() #Guardar el 'orden' de los ids de los cross-linkers como aparecen en el dump

    #Calcular numero de bonds metodo 1
    neighslist_noempty = [x for x in neighslist if x != []] #Quitar listas vacias (cross linkers sin vecinos)
    flat_list = [x for xs in neighslist_noempty for x in xs ]
    bondsnum = len(flat_list)/2 #Contar el numero de bonds
    bonds_v_t.append(bondsnum) #Guardar numero de bonds en iteracion actual

    # --- Aqui comienza codigo para detectar clusters ---
    molnum = max(atomsfull[:,1]).astype(int) #Obtener el numero de moleculas en el sistema 
    finalclusterlist=[] #lista para almacenar las listas de molecule_ids de cada cluster
    finalclustersize=[] #lista para almacenar el tamaño de cada cluster

    intlist = np.arange(1,molnum+1) #Crear lista de enteros 
    intlistnz = intlist #inicializar lista intlistnz
    clusternum = 1 #inicializar while
    while clusternum<=molnum:
        
        initmol = intlistnz[0] #molecule_id para inicializar algoritmo
        clusterlist = np.array([initmol]) #Inicalizar lista con molecule_ids pertenecientes al cluster
        for i in range(molnum):

            current_mol = clusterlist[i] #Molecula actual
            molecule = atomsfull[atomsfull[:,1]==current_mol] #Encontrar atomos que pertenecen a molecula actual
            #Identificar cross-linkers de molecula actual
            molxls = molecule[np.logical_or(molecule[:,2]==2,molecule[:,2]==4)] 
            molxls = molxls[:,0] 
            #Identificar centro de molecula actual
            center = molecule[np.logical_or(molecule[:,2]==1,molecule[:,2]==3)]
            center = center[:,0]

            index1 = 0
            bonded_molecule = [] #Lista para almacenar moleculas unidas a molecula actual
            for j in molxls: #Encontrar atomos unidos a los cross-linkers
                xlindex = xlsorder.index(j)
                b1b2 = neighslist[xlindex]
                if len(b1b2)>3:
                    print('WARNING: Rlim demasiado grande')
                if np.size(b1b2)==0: #Si no hay atomos distintos al centro unidos a ese cross linker
                    continue
                else: #Si hay atomos unidos al xl
                    for k in range(len(b1b2)):
                        aa = atomsfull[atomsfull[:,0]==b1b2[k]] #Buscar atomo en lista de atomos
                        bonded_molecule.append(aa[0,1]) #Identificar a que molecula corresponde y agregar a lista de moleculas unidas
            bonded_molecule = np.array(bonded_molecule)

            #Agregar moleculas a la lista del cluster asegurandose de no repetir
            sizebef = np.size(clusterlist) #obtener tamaño del clusterlist antes agregar moleculas
            for j in range(np.size(bonded_molecule)): #Iterar por todos las moleculas en bonded_molecule
                if bonded_molecule[j] in clusterlist: #Si la molecula ya se encuentra en el clusterlist, no agregar
                    continue
                else:
                    clusterlist = np.append(clusterlist,bonded_molecule[j]) #Si la molecula no se encuentra en el clusterlist, agregarla
            sizeaft = np.size(clusterlist) #Obtener tamaño del cluster despues de agregar moleculas

            #Si estamos en el ultimo elemento de la lista del cluster y no incrementa el tamaño, terminar
            if sizebef == sizeaft and i == sizebef-1: 
                #print(clusterlist,'Tamaño del cluster:',sizeaft) #Imprimir lista de cluster y su tamaño
                finalclusterlist.append(clusterlist) #Guardar lista de cluster
                finalclustersize.append(sizeaft) #Guardar tamaño de cluster
                break
        
        #Actualizar intlist quitando todas las moleculas pertenecientes al cluster actual
        intlistrm = np.sort(clusterlist)-1 #Crear lista con indices de moleculas del cluster
        intlist[intlistrm.astype(int)] = 0 #Remover moleculas del cluster de la lista de enteros
        intlistnz = intlist[intlist!=0] #Quitar valores nulos de la lista para reiniciar algoritmo
        
        if np.size(intlistnz) == 0: #Cuando se agoten todos los atomos terminar
            break

        clusternum+=1

    #print('Numero de clusters:',clusternum) #Imprimir el numero de clusters
    clusternum_v_t.append(clusternum) #Guardar numero de clusters de iteracion actual
    clustersize_v_t.append(max(finalclustersize)) #Guardar tamaño del cluster mas grande de iteracion actual

#Graficas
time = np.arange(1,len(clusternum_v_t)+1)
plt.scatter(time,clusternum_v_t)
plt.title('Número de clusters')
plt.xlabel('Tiempo')
plt.ylabel('Número de clusters')
plt.figure()
plt.scatter(time,clustersize_v_t)
plt.title('Tamaño del cluster más grande')
plt.xlabel('Tiempo')
plt.ylabel('Número de partículas del cluster')
plt.figure()
plt.scatter(time,bonds_v_t)
plt.title('Número de bonds')
plt.xlabel('Tiempo')
plt.ylabel('Número de bonds')
plt.show()


