'''
Script para hacer el analisis y exportar graficas del experimento de tamaño
Se realizaron 5 experimentos para diferentes tamaños de sistema -> 2k, 4k, 8k y 16k particulas
Se obtienen los promedios de los 5 experimentos de cada tamaño y se comparan los resultados
El objetivo es determinar el tamaño adecuado para los experimentos subsecuentes

Al correr el script es necesario dar como argumento el directorio de la carpeta SizeExperiments. 

La carpeta SizeExperiments debe contener una carpeta para cada tamaño de sistema llamada /Ak_particles, donde A = 2,4,8,16.
Las carpetas #k_particles deben contener una carpeta para cada experimento llamada /expB, donde B = 1,2,3,4,5
Las carpetas expB deben contener el dump del experimento dump_shearing_Ak_B.lammpstrj
También deben contener el archivo con la informacion de las mediciones de presion presion_ave_shearing_Ak_B.data

'''

def main():
    #Esto es para que ovito y matplotlib funcionen correctamente
    import os
    os.environ['OVITO_GUI_MODE'] = '1'

    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    from analysis import count_bonds_and_clusters,non_affine_sq_disp

    #El argumento debe ser el directorio SizeExperiments

    directory = sys.argv[1]
    pk = [2,4,8,16] #Lista con numeros de particulas (x10^3)

    pres_average = [] #Lista para guardar promedios de presion por iteracion
    pres_std = [] #Lista para guardar desviaciones estandar promedio 
    bonds_average = [] #Lista para guardar promedios de presion por iteracion
    clusters_average = [] #Lista para guardar promedios de clusters por iteracion
    for n in pk: #iterar numero de particulas
        print(f'{n}k start')
        pres_data = [] #Lista para almacenar datos de presion de cada experimento. Se reinicia al cambiar numero de particulas
        bonds_data = [] #Lista para almacenar datos de bonds de cada experimento. Se reinicia al cambiar numero de particulas
        clusters_data=[] #Lista para almacenar datos de clusters cada experimento. Se reinicia al cambiar numero de particulas
        for i in range(5): #iterar cada experimento
            print(f'{n}k exp{i} start')
            #---Stress Strain---#
            data = np.loadtxt(f"{directory}/{n}k_particles/exp{i+1}/presion_ave_shearing_{n}k_{i+1}.data") #Extraer datos de presion de archivo para experimento i
            pres_data.append(data[:,4]) #Exportar componente xy de presion y agregar a la lista pres_data

            #---Bonds,Clusters---#
            bonds_v_t,clusternum_v_t,_,_ = count_bonds_and_clusters(f"{directory}/{n}k_particles/exp{i+1}/dump_shearing_{n}k_{i+1}.lammpstrj")
            bonds_data.append(bonds_v_t) #Agregar experimento actual a bonds_data
            clusters_data.append(clusternum_v_t) #Agregar experimento actual a clusters_data
            print(f'{n}k exp{i} done')

        np_pres_data = np.array(pres_data) #convertir lista pres_data en numpy array
        pres_average.append(-1*np.average(np_pres_data,axis=0)) #obtener promedio de presion en cada iteracion y guardar en pres_average
        pres_std.append(np.average(np.std(np_pres_data,axis=0))) #obtener desviacion estandar promedio y guardar en pres_std

        np_bonds_data = np.array(bonds_data) #convertir lista bonds_data en numpy array
        bonds_average.append(np.average(np_bonds_data,axis=0)) #obtener promedio de bonds en cada iteracion y guardar en bonds_average

        np_clusters_data = np.array(clusters_data) #convertir lista clusters_data en numpy array
        clusters_average.append(np.average(np_clusters_data,axis=0)) #obtener promedio de clusters en cada iteracion y guardar en bonds_average

        print(f'{n}k done')

    print('Doing graphs')

    # --- Graficas ---
    x = np.linspace(0,2,len(pres_average[0])) #generar eje x para graficas

    #Grafica stress-strain promedio
    fig1,ax1 = plt.subplots()
    for i in range(4):
        ax1.plot(x,pres_average[i])
    ax1.set_title("Stress-Strain. $\dot{\gamma}=0.01$")
    ax1.set_xlabel("$\gamma$")
    ax1.set_ylabel("$\sigma_{xy}$")
    ax1.legend(("2k","4k","8k","16k"))

    #Grafica stress-strain promedio escala log-log
    fig2,ax2 = plt.subplots()
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    for i in range(4):
        ax2.scatter(x,pres_average[i])
    ax2.set_title("Stress-Strain. $\dot{\gamma}=0.01$. log-log")
    ax2.set_xlabel("$\gamma$")
    ax2.set_ylabel("$\sigma_{xy}$")
    ax2.legend(("2k","4k","8k","16k"))

    '''
    #grafica desviacion estandar stress promedio
    fig3,ax3 = plt.subplots()
    ax3.scatter(pk,pres_std)
    ax3.set_title("Average standard deviation of stress between experiments")
    ax3.set_xlabel("Number of particles ($10^3$)")
    ax3.set_ylabel("Average standard deviation")

    #Grafica bonds promedio 
    for i,n in enumerate(pk):
        plt.figure()
        plt.plot(x,bonds_average[i])
        plt.title("Bonds. %sk particles. $\dot{\gamma}=0.01$" % n )
        plt.xlabel("$\gamma$")
        plt.ylabel("$Bonds$")
    '''
    #Grafica superposicion stres-strain con bonds 16k 
    fig4,ax4 = plt.subplots()
    ax4.plot(x,pres_average[3],color='blue')          
    ax4.set_ylabel('$\sigma_{xy}$', color='blue')
    ax4.tick_params(axis='y', colors='blue')
    ax5 = ax4.twinx()                  
    ax5.plot(x,bonds_average[3],color='red')
    ax5.set_ylabel('Bonds', color='red')       
    ax5.tick_params(axis='y', colors='red')

    plt.show()

if __name__=="__main__":
    main()





