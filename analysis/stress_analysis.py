'''
Este script recibe como argumento una carpeta y archivos de presion contenidos en dicha carpeta
Realiza una grafica comparando las curvas stress-strain para los archivos de presion proporcionados

'''

def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import sys

    #Importar dumps de shearing
    input_filefolder = sys.argv[1]
    input_filenames = sys.argv[2:]

    colors = ['red','green','blue']
    markers = ['x','o','s']
    x = np.linspace(0,2,num=200)

    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    for k,filename in enumerate(input_filenames):
        input_filepath = input_filefolder+"/"+filename
        data = np.loadtxt(input_filepath)
        sigma_xy = -data[-200:,4]
        ax1.plot(x,sigma_xy,c=colors[k],marker=markers[k],linestyle='none',fillstyle='none')
        ax2.plot(x,sigma_xy,c=colors[k],marker=markers[k],linestyle='none',fillstyle='none')

    #plt.legend(legend)

    ax1.set_xlabel('$\gamma$')
    ax2.set_xlabel('$\gamma$')
    ax1.set_ylabel('$\sigma_{xy}$')
    ax2.set_ylabel('$\sigma_{xy}$')
    ax1.set_title('Stress-Strain')
    ax2.set_title('Stress-Strain. log-log')

    plt.show()

if __name__=="__main__":
    main()
