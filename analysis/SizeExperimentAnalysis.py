#Esto es para que ovito y matplotlib funcionen correctamente
import os
os.environ['OVITO_GUI_MODE'] = '1'

import numpy as np
import matplotlib.pyplot as plt
from ovito.io import import_file, export_file
from ovito.modifiers import ClusterAnalysisModifier, CreateBondsModifier

directory = input("SizeExperiment Directory: ")
pk = [2,4,8,16] #Lista con numeros de particulas (x10^3)

pres_average = [] #Lista para guardar promedios de presion por iteracion
pres_std = [] #Lista para guardar desviaciones estandar promedio 
bonds_average = [] #Lista para guardar promedios de presion por iteracion
for n in pk: #iterar numero de particulas
    pres_data = [] #Lista para almacenar datos de cada experimento. Se reinicia al cambiar numero de particulas
    bonds_data = [] #Lista para almacenar datos de cada experimento. Se reinicia al cambiar numero de particulas
    for i in range(5): #iterar cada experimento
        #---Stress Strain---#
        data = np.loadtxt(f"{directory}/{n}k_particles/exp{i+1}/presion_ave_shearing_{n}k_{i+1}.data") #Extraer datos de presion de archivo para experimento i
        pres_data.append(data[:,4]) #Exportar componente xy de presion y agregar a la lista pres_data

        #---Bonds---#
        node = import_file(f"{directory}/{n}k_particles/exp{i+1}/dump_shearing_{n}k_{i+1}.lammpstrj") #Extraer datos de presion de archivo para experimento i
        iter_num = node.source.num_frames
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
        bonds_data.append(bonds_v_t) #Agregar experimento actual a bonds_data

    np_pres_data = np.array(pres_data) #convertir lista pres_data en numpy array
    pres_average.append(-1*np.average(np_pres_data,axis=0)) #obtener promedio de presion en cada iteracion y guardar en pres_average
    pres_std.append(np.average(np.std(np_pres_data,axis=0))) #obtener desviacion estandar promedio y guardar en pres_std

    np_bonds_data = np.array(bonds_data) #convertir lista pres_data en numpy array
    bonds_average.append(np.average(np_bonds_data,axis=0)) #obtener promedio de bonds en cada iteracion y guardar en bonds_average

# --- Graficas ---
x = np.linspace(0,2,len(pres_average[0])) #generar eje x para graficas

#Grafica stress-strain promedio 
plt.figure()
for i in range(4):
    plt.plot(x,pres_average[i])
plt.title("Stress-Strain. $\dot{\gamma}=0.01$")
plt.xlabel("$\gamma$")
plt.ylabel("$\sigma_{xy}$")
plt.legend(("2k","4k","8k","16k"))
#grafica desviacion estandar stress promedio
plt.figure()
plt.scatter(pk,pres_std)
plt.title("Average standard deviation of stress between experiments")
plt.xlabel("Number of particles ($10^3$)")
plt.ylabel("Average standard deviation")
'''
#Grafica bonds promedio 
for i,n in enumerate(pk):
    plt.figure()
    plt.plot(x,bonds_average[i])
    plt.title("Bonds. %sk particles. $\dot{\gamma}=0.01$" % n )
    plt.xlabel("$\gamma$")
    plt.ylabel("$Bonds$")
'''
#Grafica superposicion stres-strain con bonds 16k 
fig1 = plt.figure()
#fig1.subplots_adjust(right=0.15)
ax1 = fig1.add_subplot()           
ax1.plot(x,pres_average[3],color='blue')          
ax1.set_ylabel('$\sigma_{xy}$', color='blue')
ax1.tick_params(axis='y', colors='blue')
ax2 = ax1.twinx()                  
ax2.plot(x,bonds_average[3],color='red')
ax2.set_ylabel('Bonds', color='red')       
ax2.tick_params(axis='y', colors='red')

plt.show()






