
def hole_analysis_pores(dumpdir):
    import numpy as np
    from random import random
    from scipy.spatial import KDTree
    import json 

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

    dir = dumpdir.split('/')
    dir[-1] = "hole_analysis_results"
    savedir = '/'.join(dir)

    print('Esribiendo resultados en archivos json')
    with open (f'{savedir}/hole_analysis_pores_distr.json','w') as f:
        json.dump(distr,f)
    with open (f'{savedir}/hole_analysis_pores_ave.json','w') as f:
        json.dump(ave,f)
    with open (f'{savedir}/hole_analysis_pores_hole_num.json','w') as f:
        json.dump(hole_num,f)
    calc_frames = calc_frames.tolist()
    with open (f'{savedir}/hole_analysis_pores_calc_frames.json','w') as f:
        json.dump(calc_frames,f)
    print('FIN')

    #return distr,ave,hole_num,calc_frames



def hole_analysis_lines(dumpdir):
    import numpy as np
    from random import random
    from scipy.spatial import KDTree
    import json

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
            sample_size = 100000

            current_distr = []
            current_sample_distr = []
            
            print('start for')
            for ss in range(sample_size):
                if ss%1000 ==0:
                    print(ss/1000)
                #Generar linea 
                line_points = L*np.random.rand(2,3)
                l21 = line_points[1,:]-line_points[0,:]
                l_len = np.linalg.norm(l21)
                lu21 = l21/l_len
                current_sample_distr.append(l_len)
                
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
                

            print('end for')

            distr.append(current_distr)
            ave.append(np.average(current_distr)) #Guardar promedios
            hole_num.append(len(current_distr)/sample_size) #Guardar numero de huecos
            sample_distr.append(current_sample_distr)
            
            sample_hist,bin = np.histogram(current_sample_distr,bins=50)
            hole_hist,trsh = np.histogram(current_distr,bins=bin)
            prob_distr.append(np.divide(hole_hist,sample_hist).tolist())
    
    dir = dumpdir.split('/')
    dir[-1] = 'hole_analysis_results'
    savedir = '/'.join(dir)

    print('Esribiendo resultados en archivos json')
    with open (f'{savedir}/hole_analysis_lines_distr.json','w') as f:
        json.dump(distr,f)
    with open (f'{savedir}/hole_analysis_lines_ave.json','w') as f:
        json.dump(ave,f)
    with open (f'{savedir}/hole_analysis_lines_hole_num.json','w') as f:
        json.dump(hole_num,f)
    calc_frames = calc_frames.tolist()
    with open (f'{savedir}/hole_analysis_lines_calc_frames.json','w') as f:
        json.dump(calc_frames,f)
    with open (f'{savedir}/hole_analysis_lines_sample_distr.json','w') as f:
        json.dump(sample_distr,f)
    with open (f'{savedir}/hole_analysis_lines_prob_distr.json','w') as f:
        json.dump(prob_distr,f)
    bin = bin.tolist()
    with open (f'{savedir}/hole_analysis_lines_bin.json','w') as f:
        json.dump(bin,f)
    print('FIN')

    #return distr,ave,hole_num,calc_frames,sample_distr,prob_distr,bin



def main():
    import sys
    import matplotlib.pyplot as plt
    import numpy as np

    dumpdir = sys.argv[1]

    #print('Iniciando analisis con metodo de poros')
    #hole_analysis_pores(dumpdir)
    #print('FIN')

    print('Iniciando analisis con metodo de lineas')
    hole_analysis_lines(dumpdir)
    print('FIN')

if __name__=="__main__":
    main()