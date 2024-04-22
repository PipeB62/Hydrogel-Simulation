import numpy as np
import matplotlib.pyplot as plt 


#Importar dumps de shearing
input_filefolder = input("Folder directory: ")
input_filenames = []
legend = []
for i in range(3):
    file = input("Filename (type end to stop): ")
    if file == "end":
        break
    else:
        input_filenames.append(file)
        input('Legend: ')


colors = ['red','green','blue']
markers = ['x','o','s']
x = np.linspace(0,2,num=200)

plt.figure()
for k,filename in enumerate(input_filenames):
    input_filepath = input_filefolder+"/"+filename
    data = np.loadtxt(input_filepath)
    sigma_xy = -data[-200:,4]
    plt.plot(x,sigma_xy,c=colors[k],marker=markers[k],linestyle='none',fillstyle='none')

plt.legend(legend)

plt.xlabel('$\gamma$')
plt.ylabel('$\sigma_{xy}$')
plt.title('Stress-Strain')
plt.show()