
def hole_analysys_box(dumpdir):
    import numpy as np
    from random import random
    from scipy.spatial import KDTree

    #obtener numero de frames
    frame_num = 0
    with open(dumpdir,"r") as dump:
        for line in dump.readlines():
            if line == "ITEM: TIMESTEP\n":
                frame_num+=1
    
    calc_frames_num = 100
    calc_frames_step = int(frame_num/calc_frames_num)
    calc_frames = np.arange(0,frame_num,calc_frames_step)
    print('Calculated frames:',calc_frames)
    
    dump = open(dumpdir,"r") #abrir dump
    distr = [] #inicializar lista para guardar radios de huecos
    ave = [] #inicializar lista para guardar promedio de radio de huecos
    hole_num = [] #inicializar lista para guardar numero de huecos
    for frame in range(frame_num):
        dump.readline() #saltar texto ITEM: TIMESTEP
        timestep = int(dump.readline()) #leer iteracion actual
        dump.readline() #saltar texto ITEM: NUMBER OF ATOMS
        n = int(dump.readline()) #leer numero de atomos
        dump.readline() #saltar texto ITEM: BOX BOUNDS xy xz yz pp pp pp
        xmin,xmax,xy = dump.readline().split() #Leer coords de caja en x. Leer xy
        ymin,ymax,yz = dump.readline().split() #Leer coords de caja en y. Leer yz
        zmin,zmax,xz = dump.readline().split() #Leer coords de caja en z. Leer xz
        #convertir a float
        xmin = float(xmin)
        xmax = float(xmax)
        xy = float(xy) 
        ymin = float(ymin)
        ymax = float(ymax)
        yz = float(yz)  
        zmin = float(zmin)
        zmax = float(zmax)
        xz = float(xz)
        L = zmax-zmin #Longitud de la caja sin deformar
        dump.readline() #saltar texto ITEM: ATOMS id type x y z
        #Leer info de atomos y quitar patches
        r = []
        for i in range(n):
            line = dump.readline().split()
            if line[1] == "2" or line[1] == "4":
                continue
            else:
                r.append([float(ii) for k,ii in enumerate(line) if k>1])
        #Analisis de huecos
        if frame in calc_frames:
            print('Frame:',frame)
            r = np.array(r) #Convertir a numpy array
            n_centers = np.size(r,0)

            T = (L/2)*np.ones((n_centers,3)) #Vector de traslacion para poner origen en esquina de la caja
            r = r+T #generar traslacion

            #Asegurar que todas las particulas esten dentro de la caja (condiciones periodicas). 
            r[:,0] = r[:,0] - (xy/L)*r[:,1] #shear inverso
            for i in range(n_centers):
                for j in range(3):
                    if r[i,j] >= L:
                        r[i,j] -= L
                    elif r[i,j] < 0:
                        r[i,j] +=L
            r[:,0] = r[:,0] + (xy/L)*r[:,1] #shear 

            #Mapeo a caja central
            for i in range(n_centers):
                if r[i,0] >= L:
                    r[i,0] -= L
                elif r[i,0] < 0:
                    r[i,0] += L
            '''
            #Revisar
            for i in range(n_centers):
                if r[i,0]==L or r[i,1]==L or r[i,2]==L:
                    print(L)
                    print(f'particula {i} en orilla. coordenadas {r[i,0]} {r[i,1]} {r[i,2]}')
            '''
            #Generar rran 
            rrand_num = 100000
            rran = np.random.rand(rrand_num,3)
            rran = L*rran
         
            #Busqueda de distancia a primer atomo para cada rran usando KDTree
            #print('start search')
            atomtree = KDTree(r,boxsize=L) #generar KDTree
            neighs = atomtree.query(rran,k=1) #Realizar busqueda
            #print('search done. initialize cleaning')
            current_distr = neighs[0].tolist() #Convertir a lista de python
            particle_rad = 2**(1/6)
            for i in range(len(current_distr)): 
                if current_distr[i] < particle_rad: #Si la distancia es menor a 0.5, esta "dentro" de una particula
                    current_distr[i] = None #Quitar valor
                else:
                    current_distr[i] -= particle_rad #Si la distancia es mayor a 0.5, esta fuera de un esfera y representa un "hueco"
            current_distr = [i for i in current_distr if i is not None] #Quitar Nones 
            #print('cleaning done')
            distr.append(current_distr) #Guardar distancias a primer vecino en distr
            ave.append(np.average(current_distr)) #Guardar promedios
            hole_num.append(len(current_distr)/rrand_num) #Guardar numero de huecos

    dump.close()   
    
    return distr,ave,hole_num,calc_frames

def main():
    import sys
    import matplotlib.pyplot as plt
    import numpy as np

    dumpdir = sys.argv[1]
    distr,ave,hole_num,calc_frames = hole_analysys_box(dumpdir)

    bin = np.linspace(0,7,50)

    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()

    ax1.hist(distr[0],bins=bin,histtype='step',label='Primer frame')
    ax1.hist(distr[-1],bins=bin,histtype='step',label='Ultimo frame')
    ax1.legend()
    ax1.set_title('Distribuciones de huecos')
    ax1.set_xlabel('Tamaño de hueco')
    ax1.set_ylabel('Numero de huecos')

    ax2.scatter(calc_frames,ave)
    ax2.set_title("Radio promedio de hueco")
    ax2.set_xlabel('Frame')
    ax2.set_ylabel('Radio')

    ax3.scatter(calc_frames,hole_num)
    ax3.set_title("Porcentaje de huecos")
    ax3.set_xlabel('Frame')
    ax3.set_ylabel('Porcentaje de huecos')
    plt.show()

    

if __name__=="__main__":
    main()