
def hole_analysis_pores(dumpdir):
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

def hole_analysis_lines(dumpdir):
    import numpy as np
    from random import random
    from scipy.spatial import KDTree

    #obtener numero de frames
    frame_num = 0
    with open(dumpdir,"r") as dump:
        for line in dump.readlines():
            if line == "ITEM: TIMESTEP\n":
                frame_num+=1
    
    calc_frames_num = 50
    calc_frames_step = int(frame_num/calc_frames_num)
    calc_frames = np.arange(0,frame_num,calc_frames_step)

    print('Calculated frames:',calc_frames)
    
    dump = open(dumpdir,"r") #abrir dump
    distr = [] #inicializar lista para guardar radios de huecos
    ave = [] #inicializar lista para guardar promedio de radio de huecos
    hole_num = [] #inicializar lista para guardar numero de huecos
    sample_distr = []
    prob_distr = []
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

            particle_rad = 2**(1/6)
            sample_size = 50000

            current_distr = []
            current_sample_distr = []

            for ss in range(sample_size):
                #Generar linea 
                #metodo1
                line_points = L*np.random.rand(2,3)
                l21 = line_points[1,:]-line_points[0,:]
                l_len = np.linalg.norm(l21)
                lu21 = l21/l_len
                current_sample_distr.append(l_len)

                #metodo2
                '''
                i = np.random.randint(n_centers)
                j = i
                while i == j:
                    j = np.random.randint(n_centers)
                line_points = np.zeros((2,3))
                check = 0
                while check == 0:
                    eu1 = np.random.random_sample((1,3))
                    eu1 = eu1/np.linalg.norm(eu1)
                    eu2 = np.random.random_sample((1,3))
                    eu2 = eu2/np.linalg.norm(eu2)
                    line_points[0,:] = r[i,:] + 0.5*eu1
                    line_points[1,:] = r[j,:] + 0.5*eu2
                    if np.any(np.logical_or(line_points>L,line_points<0)):
                        check = 0
                    else:
                        check = 1
                l21 = line_points[1,:]-line_points[0,:]
                l_len = np.linalg.norm(l21)
                lu21 = l21/l_len
                current_sample_distr.append(l_len)
                '''

                #Obtener vectores que inician en las orillas de la linea y terminan en cada particula
                ri1 = r - line_points[0,:] 
                ri2 = r - line_points[1,:] 

                #Obtener el producto punto para determinar si estan "dentro de la linea"
                ri1_d_lu21 = np.matmul(ri1,np.transpose(lu21))
                ri2_d_lu12 = np.matmul(ri2,np.transpose(-lu21))
                aa = np.logical_and(ri1_d_lu21>0,ri2_d_lu12>0)

                #Medir distancia perpendicular a linea y determinar si hay colision
                g = ri1[aa,:] - np.outer(ri1_d_lu21[aa],lu21)
                g_norm = np.linalg.norm(g,axis=1)

                if np.any(g_norm<particle_rad):
                    continue
                else:
                    current_distr.append(l_len)

            distr.append(current_distr)
            ave.append(np.average(current_distr)) #Guardar promedios
            hole_num.append(len(current_distr)/sample_size) #Guardar numero de huecos
            sample_distr.append(current_sample_distr)
            
            sample_hist,bin = np.histogram(current_sample_distr,bins=50)
            hole_hist,trsh = np.histogram(current_distr,bins=bin)
            prob_distr.append(np.divide(hole_hist,sample_hist))
    
    return distr,ave,hole_num,calc_frames,sample_distr,prob_distr,bin



def main():
    import sys
    import matplotlib.pyplot as plt
    import numpy as np

    dumpdir = sys.argv[1]

    #distr,ave,hole_num,calc_frames = hole_analysis_pores(dumpdir)
    distr,ave,hole_num,calc_frames,sample_distr,prob_distr,bin1 = hole_analysis_lines(dumpdir)
    big_r = [max(k) for k in distr]

    bin = np.linspace(0,90,90)

    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    fig4,ax4 = plt.subplots()
    fig5,ax5 = plt.subplots()
    fig6,ax6 = plt.subplots()

    ax1.hist(distr[0],bins=bin,histtype='step',label='Primer frame')
    ax1.hist(distr[-1],bins=bin,histtype='step',label='Ultimo frame')
    ax1.legend()
    ax1.set_title('Distribuciones de huecos')
    #ax1.set_xlabel('Tamaño de hueco')
    ax1.set_xlabel('Longitud de hueco')
    ax1.set_ylabel('Numero de huecos')

    ax2.scatter(calc_frames,ave)
    #ax2.set_title("Radio promedio de hueco")
    ax2.set_title('Longitud promedio de hueco')
    ax2.set_xlabel('Frame')
    ax2.set_ylabel('Radio')

    ax3.scatter(calc_frames,hole_num)
    ax3.set_title("Porcentaje de huecos")
    ax3.set_xlabel('Frame')
    ax3.set_ylabel('Porcentaje de huecos')

    ax4.scatter(calc_frames,big_r)
    #ax4.set_title("Radio del hueco mas grande")
    ax4.set_title("Longitud del hueco mas grande")
    ax4.set_xlabel('Frame')
    ax4.set_ylabel('RadioS')

    ax5.hist(sample_distr[0],bins=bin,histtype='step',label='Primer frame')
    ax5.hist(sample_distr[-1],bins=bin,histtype='step',label='Ultimo frame')
    ax5.legend()
    ax5.set_title('Distribucion de muestreo')
    ax5.set_xlabel('Longitud de linea')
    ax5.set_ylabel('Numero de lineas')

    ax6.stairs(prob_distr[0],bin1,label='Primer frame')
    ax6.stairs(prob_distr[-1],bin1,label='Ultimo frame')
    ax6.legend()
    ax6.set_title('Distribucion de probabilidad de huecos')
    ax6.set_xlabel('Longitud de linea')
    ax6.set_ylabel('Probabilidad')
    
    plt.show()
    

if __name__=="__main__":
    main()