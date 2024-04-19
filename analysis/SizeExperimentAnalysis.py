import numpy as np
import matplotlib.pyplot as plt

directory = input("SizeExperiment Directory: ")

pk = [2,4,8] #Lista con numeros de particulas (x10^3)
pres_average = [] #Lista para guardar promedios de presion por iteracion
pres_std = [] #Lista para guardar desviaciones estandar promedio 
for n in pk: #iterar numero de particulas
    pres_data = [] #Lista para almacenar datos de cada experimento. Se reinicia al cambiar numero de particulas
    for i in range(5): #iterar cada experimento
        data = np.loadtxt(f"{directory}/{n}k_particles/exp{i+1}/presion_ave_shearing_{n}k_{i+1}.data") #Extraer datos de presion de archivo para experimento i
        pres_data.append(data[:,4]) #Exportar componente xy de presion y agregar a la lista pres_data
    np_pres_data = np.array(pres_data) #convertir lista pres_data en numpy array
    pres_average.append(-1*np.average(np_pres_data,axis=0)) #obtener promedio de presion en cada iteracion y guardar en pres_average
    pres_std.append(np.average(np.std(np_pres_data,axis=0))) #obtener desviacion estandar promedio y guardar en pres_std


# --- Graficas ---
x = np.linspace(0,2,len(pres_average[0])) #generar eje x para graficas
#Grafica stress-strain promedio 
plt.figure()
plt.plot(x,pres_average[0])
plt.plot(x,pres_average[1])
plt.plot(x,pres_average[2])
plt.title("Stress-Strain. $\dot{\gamma}=0.01$")
plt.xlabel("$\gamma$")
plt.ylabel("$\sigma_{xy}$")
plt.legend(("2k","4k","8k"))
#grafica desviacion estandar promedio
plt.figure()
plt.scatter(pk,pres_std)
plt.title("Average standard deviation between experiments")
plt.xlabel("Number of particles ($10^3$)")
plt.ylabel("Average standard deviation")

plt.show()






