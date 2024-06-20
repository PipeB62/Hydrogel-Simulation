def load_json(dir,name):
    import json
    with open (f'{dir}/{name}.json') as f:
        data = json.load(f)
    if isinstance(data,str):
        data = json.loads(data)
    return(data)

def derivative(data,step):
    a = []
    for i in range(1,len(data)):
        a.append((data[i]-data[i-1])/step)
    return a

#Grafica stress-strain con los datos de analisis comparados verticalmente
def figura1():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio 
    dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop

    #importar datos de estres y strain promediados
    stress_v_t_ave = [-i for i in load_json(dir,"stress_v_t_ave")]
    strain = load_json(dir,"strain")
    for i,x in enumerate(strain):
        if x < 0:
            strain[i] = x+1
    strain.pop(0) #Quitar el frame 0 del strain
    #Obtener el yield strain. Lugar de maximo estres
    yield_ix = stress_v_t_ave.index(max(stress_v_t_ave))
    yield_strain = strain[yield_ix]

    #derivada del estres con respecto al strain
    dsigma_dgamma = derivative(stress_v_t_ave,0.01)
    maxderix = dsigma_dgamma.index(max(dsigma_dgamma))
    max_der_strain = strain[maxderix+1]

    #Crear figura
    fig = plt.figure(1,figsize=(5,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(5,4)
    gs.update(wspace=0.2,hspace=0)

    #Crear subplot1 y graficar
    subplot1 = fig.add_subplot(gs[0,0:4])
    plt.plot(strain[:100],stress_v_t_ave[:100],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    plt.ylabel('$\sigma_{xy}$')

    #Limites
    plt.xlim((0,1))
    plt.ylim((-0.01,0.1))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    plt.yticks(np.arange(0,0.11,0.03))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #crear subplot 5 y graficar 
    subplot5 = fig.add_subplot(gs[1,0:4])
    plt.plot(strain[1:100],dsigma_dgamma[0:99],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\frac{d \sigma_{xy}}{d \gamma}$")

    #Limites
    plt.xlim((0,1))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de bonds
    bonds_v_t_ave = load_json(dir,"bonds_v_t_ave")
    bonds_v_t_ave.pop(0)

    #crear subplot 2 y graficar
    subplot2 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain[:100],bonds_v_t_ave[:100],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel('$N_b$')

    #Limites
    plt.xlim((0,1))
    plt.ylim((21300,21700))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de clusters 
    clusternum_v_t_ave = load_json(dir,"clusternum_v_t_ave")
    clusternum_v_t_ave.pop(0)

    #crear subplot 3 y graficar
    subplot3 = fig.add_subplot(gs[3,0:4])
    plt.plot(strain[:100],clusternum_v_t_ave[:100],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel('$N_c$')

    #Limites
    plt.xlim((0,1))
    plt.ylim((-2,28))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de non affine squared displacement 
    dsq_v_t_ave = load_json(dir,"dsq_v_t_ave")
    dsq_v_t_ave.pop(0)

    #crear subplot 4 y graficar
    subplot4 = fig.add_subplot(gs[4,0:4])
    plt.plot(strain[:100],dsq_v_t_ave[:100],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel('$D^2$')
    plt.xlabel(r'$\gamma$')

    #Limites
    plt.xlim((0,1))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)

    plt.savefig('D:/ExperimentData/SizeExperiment/graficas/stress_strain_analysis.pdf',dpi=300,bbox_inches='tight')
    
    plt.show()

#Grafica de curvas stress strain de diferentes tamaños de sistema
def figura2():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio 

    strain = load_json("D:/ExperimentData/SizeExperiment/16k_particles/averaged_data","strain")
    for i,x in enumerate(strain):
        if x < 0:
            strain[i] = x+1
    strain.pop(0) #Quitar el frame 0 del strain

    allsigma = []
    for i in [2,4,8,16]:
        dir = f"D:/ExperimentData/SizeExperiment/{i}k_particles/averaged_data"
        allsigma.append([-k for k in load_json(dir,"stress_v_t_ave")])
    
    fig,ax = plt.subplots(figsize=(5,5))

    ax.plot(strain[:100],allsigma[0][:100],marker='o',markersize=4,mfc='w',alpha = 0.9,linestyle='-',label="2k particles")
    ax.plot(strain[:100],allsigma[1][:100],marker='^',markersize=4,mfc='w',alpha = 0.9,linestyle='-',label="4k particles")
    ax.plot(strain[:100],allsigma[2][:100],marker='D',markersize=4,mfc='w',alpha = 0.9,linestyle='-',label="8k particles")
    ax.plot(strain[:100],allsigma[3][:100],marker='s',markersize=4,mfc='w',alpha = 0.9,linestyle='-',label="16k particles")

    ax.legend()

    ax.set_xlabel("$\gamma$")
    ax.set_ylabel("$\sigma_{xy}$")

    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    x1,x2,y1,y2 = 0, 0.2, -0.0025, 0.005
    axins = ax.inset_axes(
        [0.55, 0.1, 0.4, 0.3],
        xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
    axins.plot(strain[:20],allsigma[0][:20],marker='o',markersize=4,mfc='w',alpha = 0.9,linestyle='-')
    axins.plot(strain[:20],allsigma[1][:20],marker='^',markersize=4,mfc='w',alpha = 0.9,linestyle='-')
    axins.plot(strain[:20],allsigma[2][:20],marker='D',markersize=4,mfc='w',alpha = 0.9,linestyle='-')
    axins.plot(strain[:20],allsigma[3][:20],marker='s',markersize=4,mfc='w',alpha = 0.9,linestyle='-')
    axins.minorticks_on()
    axins.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    axins.tick_params(direction='in',which='major',length=5,right=True,top=True)
    axins.set_xticks(np.arange(0,0.21,0.2))
    axins.set_yticks(np.arange(0,0.0051,0.005))
    ax.indicate_inset_zoom(axins, edgecolor="black")

    plt.savefig('D:/ExperimentData/SizeExperiment/graficas/stress_strain_sizes.pdf',dpi=300,bbox_inches='tight')
    plt.show()

def figura3():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio 

    strain = load_json("D:/ExperimentData/SizeExperiment/16k_particles/averaged_data","strain")
    for i,x in enumerate(strain):
        if x < 0:
            strain[i] = x+1
    strain.pop(0) #Quitar el frame 0 del strain

    allsigma = []
    for i in [2,4,8,16]:
        dir = f"D:/ExperimentData/SizeExperiment/{i}k_particles/averaged_data"
        allsigma.append([-k for k in load_json(dir,"stress_v_t_ave")])
    
    fig,ax = plt.subplots(figsize=(5,5))

    ax.plot(strain[:50],allsigma[0][:50],marker='o',markersize=4,mfc='w',alpha = 0.9,linestyle='--',label="2k particles")
    ax.plot(strain[:50],allsigma[1][:50],marker='^',markersize=4,mfc='w',alpha = 0.9,linestyle='--',label="4k particles")
    ax.plot(strain[:50],allsigma[2][:50],marker='D',markersize=4,mfc='w',alpha = 0.9,linestyle='--',label="8k particles")
    ax.plot(strain[:50],allsigma[3][:50],marker='s',markersize=4,mfc='w',alpha = 0.9,linestyle='--',label="16k particles")

    ax.legend()

    ax.set_xlabel("$\gamma$")
    ax.set_ylabel("$\sigma_{xy}$")

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    plt.savefig('D:/ExperimentData/SizeExperiment/graficas/stress_strain_loglog.pdf',dpi=300,bbox_inches='tight')
    plt.show()



def main():
    #figura1()
    #figura2()
    figura3()

if __name__ == '__main__':
    main()