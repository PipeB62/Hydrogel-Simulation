import numpy as np
import matplotlib.pyplot as plt 

#Version
version = input('Version: ')

#Importar dumps de shearing
input_filefolder = "/home/pipe/lammps/code/HYDROGELS/selfassembly2/shearing/estres_v"+version+"/"
input_filenames = ["presion_ave_1em4.dat","presion_ave_1em3.dat","presion_ave_1em2.dat"]

input_filenames.remove("presion_ave_1em4.dat")
#input_filenames.remove("presion_ave_1em3.dat")

colors = ['red','green','blue']
markers = ['x','o','s']
x = np.linspace(0,2,num=200)
#x = np.linspace(0,1,num=100)

plt.figure()
for k,filename in enumerate(input_filenames):
    input_filepath = input_filefolder+filename
    data = np.loadtxt(input_filepath)

    #sigma_xy = data[1:-1,3] #(1)
    sigma_xy = -data[:,4] #(2)

    plt.plot(x,sigma_xy,c=colors[k],marker=markers[k],linestyle='none',fillstyle='none')

plt.legend(['1e-4','1e-3','1e-2'],title='$\dot{\gamma}$')
plt.xlabel('$\gamma$')
plt.ylabel('$\sigma_{xy}$')
plt.title('Stress-Strain')
plt.show()
