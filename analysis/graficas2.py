def load_json(dir,name):
    import json
    with open (f'{dir}/{name}.json') as f:
        data = json.load(f)
    if isinstance(data,str):
        data = json.loads(data)
    return(data)

def load_data(dir,name):
    import numpy as np
    a = np.loadtxt(f"{dir}/{name}")
    return a.tolist()

def derivative(data,step):
    a = []
    for i in range(1,len(data)):
        a.append((data[i]-data[i-1])/step)
    return a

def derivative2(data,step):
    a = []
    for i in range(1,len(data)-1):
        a.append((data[i+1]-data[i-1])/(2*step))
    return a

#Grafica stress-strain con los datos de analisis comparados verticalmente
def stress_strain_basic_analysis():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/averaged_data" #pc escritorio 
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop

    #importar datos de estres y strain promediados
    stress_v_t_ave = load_json(dir,"analysis_stress_ave")
    stress_v_t_std = load_json(dir,"analysis_stress_std")
    strain = load_json(dir,"analysis_strain_ave")
    strain.pop(0) #Quitar el frame 0 del strain

    #Obtener el yield strain. Lugar de maximo estres
    yield_ix = stress_v_t_ave.index(max(stress_v_t_ave))
    yield_strain = strain[yield_ix]

    #derivada del estres con respecto al strain
    dsigma_dgamma = derivative2(stress_v_t_ave,0.01)
    maxderix = dsigma_dgamma.index(max(dsigma_dgamma))
    max_der_strain = strain[maxderix+1]
    max_der_strain = 0.5

    #Crear figura
    fig = plt.figure(1,figsize=(8,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(5,4)
    gs.update(wspace=0.2,hspace=0)

    #Crear subplot1 y graficar
    subplot1 = fig.add_subplot(gs[0,0:4])
    plt.plot(strain[:],stress_v_t_ave[:])
    plt.errorbar(strain[0::30],stress_v_t_ave[0::30],yerr=stress_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    plt.ylabel('$\sigma_{xy}$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((-0.01,0.1))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.yticks(np.arange(0,0.11,0.03))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #crear subplot 2 y graficar 
    subplot2 = fig.add_subplot(gs[1,0:4])
    plt.plot(strain[1:-1],dsigma_dgamma[:],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\frac{d \sigma_{xy}}{d \gamma}$")

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de bonds
    bonds_v_t_ave = load_json(dir,"analysis_bonds_v_t_ave")
    bonds_v_t_std = load_json(dir,"analysis_bonds_v_t_std")

    #crear subplot 3 y graficar
    subplot3 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain[:],bonds_v_t_ave[:])
    plt.errorbar(strain[0::30],bonds_v_t_ave[0::30],yerr=bonds_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel('$N_b$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((21300,21700))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de clusters 
    clusternum_v_t_ave = load_json(dir,"analysis_clusternum_v_t_ave")
    clusternum_v_t_std = load_json(dir,"analysis_clusternum_v_t_std")

    #crear subplot 4 y graficar
    subplot4 = fig.add_subplot(gs[3,0:4])
    plt.plot(strain[:],clusternum_v_t_ave[:])
    plt.errorbar(strain[0::30],clusternum_v_t_ave[0::30],yerr=clusternum_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel('$N_c$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((-2,28))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de tamaño del cluster mas grande
    bigclustersz_v_t_ave = load_json(dir,"analysis_bigclustersz_v_t_ave")
    bigclustersz_v_t_std = load_json(dir,"analysis_bigclustersz_v_t_std")

    #crear subplot 4 y graficar
    subplot4 = fig.add_subplot(gs[4,0:4])
    plt.plot(strain[:],bigclustersz_v_t_ave[:])
    plt.errorbar(strain[0::30],bigclustersz_v_t_ave[0::30],yerr=bigclustersz_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel('Cluster Size')
    plt.xlabel(r'$\gamma$')

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #plt.savefig('D:/ExperimentData/SizeExperiment/graficas/stress_strain_analysis.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/Graficas/stress_strain_basic_analysis.pdf',dpi=300,bbox_inches='tight')

    plt.show()

def stress_strain_basic_analysis_30k():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/averaged_data" #pc escritorio 
    dir30k = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/SizeTest/exp1/analysis_results"
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop

    #importar datos de estres y strain promediados
    stress_v_t_ave = load_json(dir,"analysis_stress_ave")
    stress_v_t_std = load_json(dir,"analysis_stress_std")
    strain = load_json(dir,"analysis_strain_ave")
    strain.pop(0) #Quitar el frame 0 del strain

    stress30k = load_json(dir30k,"analysis_stress")
    strain30k = load_json(dir30k,"analysis_strain")
    strain30k.pop(0) #Quitar el frame 0 del strain

    #Obtener el yield strain. Lugar de maximo estres
    yield_ix = stress_v_t_ave.index(max(stress_v_t_ave))
    yield_strain = strain[yield_ix]

    #derivada del estres con respecto al strain
    dsigma_dgamma = derivative2(stress_v_t_ave,0.01)
    maxderix = dsigma_dgamma.index(max(dsigma_dgamma))
    max_der_strain = strain[maxderix+1]
    max_der_strain = 0.5

    #Crear figura
    fig = plt.figure(1,figsize=(8,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(5,4)
    gs.update(wspace=0.2,hspace=0)

    #Crear subplot1 y graficar
    subplot1 = fig.add_subplot(gs[0,0:4])
    plt.plot(strain[:],stress_v_t_ave[:],label="16k")
    plt.errorbar(strain[0::30],stress_v_t_ave[0::30],yerr=stress_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    plt.plot(strain30k[:],stress30k[:],label="32k")
    plt.legend()

    #Labels
    plt.ylabel('$\sigma_{xy}$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((-0.01,0.1))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.yticks(np.arange(0,0.11,0.03))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #crear subplot 2 y graficar 
    subplot2 = fig.add_subplot(gs[1,0:4])
    plt.plot(strain[1:-1],dsigma_dgamma[:],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\frac{d \sigma_{xy}}{d \gamma}$")

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de bonds
    bonds_v_t_ave = load_json(dir,"analysis_bonds_v_t_ave")
    bonds_v_t_std = load_json(dir,"analysis_bonds_v_t_std")
    maxbonds1 = max(bonds_v_t_ave)
    bonds_v_t_ave = [i/maxbonds1 for i in bonds_v_t_ave]
    bonds_v_t_std = [i/maxbonds1 for i in bonds_v_t_std]

    bonds30k = load_json(dir30k,"analysis_bonds_v_t")
    maxbonds2 = max(bonds30k)
    bonds30k = [i/maxbonds2 for i in bonds30k]

    #crear subplot 3 y graficar
    subplot3 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain[:],bonds_v_t_ave[:],label="16k")
    plt.errorbar(strain[0::30],bonds_v_t_ave[0::30],yerr=bonds_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    plt.plot(strain30k[:],bonds30k[:],label="32k")
    plt.legend()

    #Labels 
    plt.ylabel('$N_b$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((21300,21700))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de clusters 
    clusternum_v_t_ave = load_json(dir,"analysis_clusternum_v_t_ave")
    clusternum_v_t_std = load_json(dir,"analysis_clusternum_v_t_std")
    maxcluster1 = max(clusternum_v_t_ave)
    clusternum_v_t_ave = [i/maxcluster1 for i in clusternum_v_t_ave]
    clusternum_v_t_std = [i/maxcluster1 for i in clusternum_v_t_std]

    clusternum30k = load_json(dir30k,"analysis_clusternum_v_t")
    maxcluster2 = max(clusternum30k)
    clusternum30k = [i/maxcluster2 for i in clusternum30k]

    #crear subplot 4 y graficar
    subplot4 = fig.add_subplot(gs[3,0:4])
    plt.plot(strain[:],clusternum_v_t_ave[:],label="16k")
    plt.errorbar(strain[0::30],clusternum_v_t_ave[0::30],yerr=clusternum_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    plt.plot(strain30k[:],clusternum30k[:],label="32k")
    plt.legend()

    #Labels 
    plt.ylabel('$N_c$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((-2,28))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de tamaño del cluster mas grande
    bigclustersz_v_t_ave = load_json(dir,"analysis_bigclustersz_v_t_ave")
    bigclustersz_v_t_std = load_json(dir,"analysis_bigclustersz_v_t_std")
    maxclsz1 = max(bigclustersz_v_t_ave)
    bigclustersz_v_t_ave = [i/maxclsz1 for i in bigclustersz_v_t_ave]
    bigclustersz_v_t_std = [i/maxclsz1 for i in bigclustersz_v_t_std]

    bigclustersz30k = load_json(dir30k,"analysis_bigclustersz_v_t")
    maxclsz2 = max(bigclustersz30k)
    bigclustersz30k = [i/maxclsz2 for i in bigclustersz30k]

    #crear subplot 4 y graficar
    subplot4 = fig.add_subplot(gs[4,0:4])
    plt.plot(strain[:],bigclustersz_v_t_ave[:],label="16k")
    plt.errorbar(strain[0::30],bigclustersz_v_t_ave[0::30],yerr=bigclustersz_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    plt.plot(strain30k[:],bigclustersz30k[:],label="32k")
    plt.legend()

    #Labels 
    plt.ylabel('Cluster Size')
    plt.xlabel(r'$\gamma$')

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #plt.savefig('D:/ExperimentData/SizeExperiment/graficas/stress_strain_analysis.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('/media/felipe/KINGSTON/Hydrogel_sim_experiments/SizeTest/exp1/Graficas/stress_strain_basic_analysis.pdf',dpi=300,bbox_inches='tight')

    plt.show()

def stress_strain_xl_distance():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/averaged_data" #pc escritorio 
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop

    #importar datos de estres y strain promediados
    stress_v_t_ave = load_json(dir,"analysis_stress_ave")
    stress_v_t_std = load_json(dir,"analysis_stress_std")
    strain = load_json(dir,"analysis_strain_ave")
    strain_nz = strain[1:]

    calcframes = load_json(dir,"xldistance_calcframes_ave")
    c = 0
    strain_cf = []
    for i,s in enumerate(strain):
        if i == calcframes[c]-1:
            c+=1
            strain_cf.append(s)

    #Obtener el yield strain. Lugar de maximo estres
    yield_ix = stress_v_t_ave.index(max(stress_v_t_ave))
    yield_strain = strain[yield_ix]
    yield_strain = 0.9 #sobreescribir al tanteo

    #derivada del estres con respecto al strain
    dsigma_dgamma = derivative2(stress_v_t_ave,0.01)
    maxderix = dsigma_dgamma.index(max(dsigma_dgamma))
    max_der_strain = strain[maxderix+1]
    max_der_strain = 0.5 #sobreescribir al tanteo

    #Crear figura
    fig = plt.figure(1,figsize=(8,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(5,4)
    gs.update(wspace=0.2,hspace=0)

    #Crear subplot1 y graficar
    subplot1 = fig.add_subplot(gs[0,0:4])
    plt.plot(strain_nz[:],stress_v_t_ave[:])
    plt.errorbar(strain_nz[0::30],stress_v_t_ave[0::30],yerr=stress_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    plt.ylabel('$\sigma_{xy}$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((-0.01,0.1))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.yticks(np.arange(0,0.11,0.03))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #crear subplot 2 y graficar 
    subplot2 = fig.add_subplot(gs[1,0:4])
    plt.plot(strain_nz[1:-1],dsigma_dgamma[:],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\frac{d \sigma_{xy}}{d \gamma}$")

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de xl coord
    mean_xl_coord_v_t_ave = load_json(dir,"mean_xl_coord_v_t_ave")
    mean_xl_coord_v_t_std = load_json(dir,"mean_xl_coord_v_t_std")

    #crear subplot 3 y graficar
    subplot3 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain_cf[:],mean_xl_coord_v_t_ave[:])
    plt.errorbar(strain_cf[::3],mean_xl_coord_v_t_ave[::3],yerr=mean_xl_coord_v_t_std[::3],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r'$\langle C_{xl} \rangle$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((21300,21700))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de xl distance
    mean_xl_distance_v_t_ave = load_json(dir,"mean_xl_distance_v_t_ave")
    mean_xl_distance_v_t_std = load_json(dir,"mean_xl_distance_v_t_std")

    #crear subplot 4 y graficar
    subplot4 = fig.add_subplot(gs[3,0:4])
    plt.plot(strain_cf[:],mean_xl_distance_v_t_ave[:])
    plt.errorbar(strain_cf[::3],mean_xl_distance_v_t_ave[::3],yerr=mean_xl_distance_v_t_std[::3],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r'$\langle D_{xl} \rangle$')
    plt.xlabel(r'$\gamma$')
    #Limites
    plt.xlim((0,10))
    #plt.ylim((-2,28))

    #Ticks
    #plt.xticks(np.arange(0,10.1,0.5))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #imporar std_xldistance
    std_xl_distance_v_t_ave = load_json(dir,"std_xl_distance_v_t_ave")
    std_xl_distance_v_t_std = load_json(dir,"std_xl_distance_v_t_std")

    #crear subplot 5 y graficar
    subplot5 = fig.add_subplot(gs[4,0:4])
    plt.plot(strain_cf[:],std_xl_distance_v_t_ave[:])
    plt.errorbar(strain_cf[::3],std_xl_distance_v_t_ave[::3],yerr=std_xl_distance_v_t_std[::3],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r'$\sigma \left( D_{xl} \right)$')
    plt.xlabel(r'$\gamma$')
    #Limites
    plt.xlim((0,10))
    #plt.ylim((-2,28))

    #Ticks
    #plt.xticks(np.arange(0,10.1,0.5))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    #plt.tick_params(labelbottom=False)

    #plt.savefig('D:/ExperimentData/SizeExperiment/graficas/stress_strain_xldistance.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/Graficas/stress_strain_xldistance.pdf',dpi=300,bbox_inches='tight')
    
    plt.show()

def stress_strain_xl_distance_30k():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/averaged_data" #pc escritorio 
    dir30k = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/SizeTest/exp1/analysis_results"
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop

    #importar datos de estres y strain promediados
    stress_v_t_ave = load_json(dir,"analysis_stress_ave")
    stress_v_t_std = load_json(dir,"analysis_stress_std")
    strain = load_json(dir,"analysis_strain_ave")
    strain_nz = strain[1:]

    stress30k = load_json(dir30k,"analysis_stress")
    strain30k = load_json(dir30k,"analysis_strain")
    strain30k_nz = strain30k[1:]

    calcframes = load_json(dir,"xldistance_calcframes_ave")
    c = 0
    strain_cf = []
    for i,s in enumerate(strain):
        if i == calcframes[c]-1:
            c+=1
            strain_cf.append(s)

    #Obtener el yield strain. Lugar de maximo estres
    yield_ix = stress_v_t_ave.index(max(stress_v_t_ave))
    yield_strain = strain[yield_ix]

    #derivada del estres con respecto al strain
    dsigma_dgamma = derivative2(stress_v_t_ave,0.01)
    maxderix = dsigma_dgamma.index(max(dsigma_dgamma))
    max_der_strain = strain[maxderix+1]
    max_der_strain = 0.5

    #Crear figura
    fig = plt.figure(1,figsize=(8,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(4,4)
    gs.update(wspace=0.2,hspace=0)

    #Crear subplot1 y graficar
    subplot1 = fig.add_subplot(gs[0,0:4])
    plt.plot(strain_nz[:],stress_v_t_ave[:],label="16k")
    plt.errorbar(strain_nz[0::30],stress_v_t_ave[0::30],yerr=stress_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    plt.plot(strain30k_nz,stress30k,label="32k")
    plt.legend()

    #Labels
    plt.ylabel('$\sigma_{xy}$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((-0.01,0.1))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.yticks(np.arange(0,0.11,0.03))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #crear subplot 2 y graficar 
    subplot2 = fig.add_subplot(gs[1,0:4])
    plt.plot(strain_nz[1:-1],dsigma_dgamma[:],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\frac{d \sigma_{xy}}{d \gamma}$")

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de xl coord
    mean_xl_coord_v_t_ave = load_json(dir,"mean_xl_coord_v_t_ave")
    mean_xl_coord_v_t_std = load_json(dir,"mean_xl_coord_v_t_std")

    mean_xl_coord30k = load_json(dir30k,"mean_xl_coord_v_t")

    #crear subplot 3 y graficar
    subplot3 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain_cf[:],mean_xl_coord_v_t_ave[:],linestyle="--",label="16k")
    plt.errorbar(strain_cf[:],mean_xl_coord_v_t_ave[:],yerr=mean_xl_coord_v_t_std[:],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    plt.plot(strain_cf[:],mean_xl_coord30k[:],linestyle="--",label="32k")
    plt.legend()

    #Labels 
    plt.ylabel(r'$\langle C_{xl} \rangle$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((21300,21700))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de xl distance
    mean_xl_distance_v_t_ave = load_json(dir,"mean_xl_distance_v_t_ave")
    mean_xl_distance_v_t_std = load_json(dir,"mean_xl_distance_v_t_std")

    mean_xl_distance30k = load_json(dir30k,"mean_xl_distance_v_t")

    #crear subplot 4 y graficar
    subplot4 = fig.add_subplot(gs[3,0:4])
    plt.plot(strain_cf[:],mean_xl_distance_v_t_ave[:],linestyle="--",label="16k")
    plt.errorbar(strain_cf[:],mean_xl_distance_v_t_ave[:],yerr=mean_xl_distance_v_t_std[:],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    plt.plot(strain_cf[:],mean_xl_distance30k[:],linestyle="--",label="32k")
    plt.legend()

    #Labels 
    plt.ylabel(r'$\langle D_{xl} \rangle$')
    plt.xlabel(r'$\gamma$')
    #Limites
    plt.xlim((0,10))
    #plt.ylim((-2,28))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #plt.savefig('D:/ExperimentData/SizeExperiment/graficas/stress_strain_xldistance.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('/media/felipe/KINGSTON/Hydrogel_sim_experiments/SizeTest/exp1/Graficas/stress_strain_xldistance_30k.pdf',dpi=300,bbox_inches='tight')
    
    plt.show()

def stress_strain_loglog():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/averaged_data" #usb

    stress = load_json(dir,"analysis_stress_ave")
    stress_std = load_json(dir,"analysis_stress_std")
    strain = load_json(dir,"analysis_strain_ave")
    strain.pop(0)
    
    fig,ax = plt.subplots(figsize=(5,5))

    ax.errorbar(strain[:],stress[:],yerr=stress_std[:],marker='o',markersize=4,mfc='w',alpha = 0.9,linestyle='--',capsize=1)

    ax.set_xlabel("$\gamma$")
    ax.set_ylabel("$\sigma_{xy}$")

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #fig.savefig("E:/Hydrogel_sim_experiments/Test1/exp1/Graficas/stress_strain_loglog_1.pdf",dpi=300,bbox_inches='tight')
    fig.savefig('/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/Graficas/stress_strain_loglog.pdf',dpi=300,bbox_inches='tight')

    plt.show()

#Grafica de curvas stress strain de diferentes tamaños de sistema
def stress_strain_sz():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio 
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop

    strain = load_json(dir,"strain")
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

def pot_ene_formation():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SelfAssemby/exp1_mp" #pc escritorio 
    #dir = "E:/Hydrogel_sim_experiments/Test1/exp1"
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/averaged_data"

    fig1,ax1 = plt.subplots(figsize=(5,5))

    pot_ene = load_json(dir,"pot_ene_formation_ave")
    pot_ene_std = load_json(dir,"pot_ene_formation_std")
    frames = load_json(dir,"pot_ene_formation_frames_ave")
    time = [float(i)*0.002 for i in frames]

    ax1.plot(time,pot_ene)
    ax1.errorbar(time[::30],pot_ene[::30],yerr=pot_ene_std[::30],capsize=2,linestyle="None",marker="s",mfc="w",markersize=3)

    #Labels
    ax1.set_ylabel('Potential Energy')
    ax1.set_xlabel('Time')

    #Ticks
    #ax1.set_xticks(np.arange(0,9))
    #ax.set_yticks(np.arange(0,1.31,0.4))
    ax1.minorticks_on()
    ax1.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax1.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig1.savefig("/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/Graficas/pot_ene_formation.pdf",dpi=300,bbox_inches='tight')

    plt.show()

#Grafica de curvas stress strain en escala log-log
def stress_strain_loglog_size():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/averaged_data"

    strain = load_json(dir,"analysis_strain_ave")
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

#Grafica stress-strain con los datos de huecos (poros) comparados verticalmente
def figura4():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_hole_data" #pc escritorio 
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop
    dir2 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio

    #importar datos de estres y strain promediados
    stress_v_t_ave = [-i for i in load_json(dir2,"stress_v_t_ave")]
    strain = load_json(dir2,"strain")
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

    #crear subplot 2 y graficar 
    subplot2 = fig.add_subplot(gs[1,0:4])
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

    #importar datos de tamaño promedio de hueco
    hole_sz = load_json(dir,"av_hole_sz_pore_ave")
    hole_sz.pop(0)

    #crear subplot 3 y graficar 
    subplot3 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain[:100],hole_sz[:100],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\langle r_{p} \rangle$")

    #Limites
    plt.xlim((0,1))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de hueco mas grande
    max_hole_sz = load_json(dir,"max_hole_sz_pore_ave")
    max_hole_sz.pop(0)

    #crear subplot 4 y graficar 
    subplot4 = fig.add_subplot(gs[3,0:4])
    plt.plot(strain[:100],max_hole_sz[:100],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\max(r_{p}) $")

    #Limites
    plt.xlim((0,1))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de huecos
    hole_num = load_json(dir,"hole_num_pore_ave")
    hole_num.pop(0)

    #crear subplot 5 y graficar 
    subplot5 = fig.add_subplot(gs[4,0:4])
    plt.plot(strain[:100],hole_num[:100],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$N_{p}/N_{tot}$")
    plt.xlabel(r'$\gamma$')

    #Limites
    plt.xlim((0,1))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    #plt.tick_params(labelbottom=False)

    #plt.savefig('D:/ExperimentData/SizeExperiment/graficas/stress_strain_loglog.pdf',dpi=300,bbox_inches='tight')
    plt.savefig("/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/hole_pore_analysis.pdf",dpi=300,bbox_inches='tight')
    plt.show()

#Grafica stress-strain con los datos de huecos (lines) comparados verticalmente
def figura5():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_hole_data" #pc escritorio 
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop
    dir2 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio

    #importar datos de estres y strain promediados
    stress_v_t_ave = [-i for i in load_json(dir2,"stress_v_t_ave")]
    strain = load_json(dir2,"strain")
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

    #crear subplot 2 y graficar 
    subplot2 = fig.add_subplot(gs[1,0:4])
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

    #importar datos de tamaño promedio de hueco
    hole_sz = load_json(dir,"av_hole_sz_line_ave")
    hole_sz.pop(0)

    #crear subplot 3 y graficar 
    strain2 = np.linspace(0,1,25)
    subplot3 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain2,hole_sz[:25],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\langle L_{h} \rangle$")

    #Limites
    plt.xlim((0,1))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de hueco mas grande
    max_hole_sz = load_json(dir,"max_hole_sz_line_ave")
    max_hole_sz.pop(0)

    #crear subplot 4 y graficar 
    subplot4 = fig.add_subplot(gs[3,0:4])
    plt.plot(strain2,max_hole_sz[:25],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\max(L_{h}) $")

    #Limites
    plt.xlim((0,1))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de huecos
    hole_num = load_json(dir,"hole_num_line_ave")
    hole_num.pop(0)

    #crear subplot 5 y graficar 
    subplot5 = fig.add_subplot(gs[4,0:4])
    plt.plot(strain2,hole_num[:25],marker='s',markersize=2,mfc='w',linestyle='-')
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$N_{h}/N_{tot}$")
    plt.xlabel(r'$\gamma$')

    #Limites
    plt.xlim((0,1))

    #Ticks
    plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    #plt.tick_params(labelbottom=False)

    #plt.savefig('D:/ExperimentData/SizeExperiment/graficas/stress_strain_loglog.pdf',dpi=300,bbox_inches='tight')
    plt.savefig("/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/hole_line_analysis.pdf",dpi=300,bbox_inches='tight')
    plt.show()

#Distribuciones de huecos (poros) 
def figura6():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_hole_data" #pc escritorio
    dir2 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data"
    dir3 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/exp5/hole_analysis_results"

    distrs = load_json(dir,"hole_distr_pore_ave")
    bbin = load_json(dir,"binpore")
    strain = load_json(dir2,"strain")
    for i,x in enumerate(strain):
        if x < 0:
            strain[i] = x+1
    frames = load_json(dir3,"hole_analysis_pores_calc_frames")

    xbin = []
    for i in range(1,len(bbin)):
        xbin.append((bbin[i]+bbin[i-1])/2)
    xbin = np.array(xbin)

    fig,ax = plt.subplots(figsize=(5,5))

    ax.stairs(distrs[0], bbin, label=f"$\gamma = {np.round(strain[frames[0]],decimals=2)}$",color="b")
    ax.plot(xbin,distrs[0],linestyle="None",marker="s",mfc="w",color="b",markersize=3)
    ax.stairs(distrs[35], bbin, label=f"$\gamma = {np.round(strain[frames[30]],decimals=2)}$",color="g")
    ax.plot(xbin,distrs[35],linestyle="None",marker="s",mfc="w",color="g",markersize=3)
    ax.stairs(distrs[50], bbin, label=f"$\gamma = {np.round(strain[frames[50]],decimals=2)}$",color="y")
    ax.plot(xbin,distrs[50],linestyle="None",marker="s",mfc="w",color="y",markersize=3)

    ax.set_xlim(0,4)
    ax.set_xlabel("$r_p$")
    ax.set_ylabel("$N_p$")

    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    ax.legend()

    fig.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/hole_pore_distr.pdf',dpi=300,bbox_inches='tight')

    plt.show()

#Distribuciones de huecos (lineas) 
def figura7():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_hole_data" #pc escritorio
    dir2 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data"
    dir3 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/exp5/hole_analysis_results"

    distrs = load_json(dir,"hole_distr_line_ave")
    bbin = load_json(dir,"binline")

    xbin = []
    for i in range(1,len(bbin)):
        xbin.append((bbin[i]+bbin[i-1])/2)
    xbin = np.array(xbin)

    strain = load_json(dir2,"strain")
    for i,x in enumerate(strain):
        if x < 0:
            strain[i] = x+1
    frames = [int(i) for i in load_json(dir3,"hole_analysis_line_calc_frames")]

    fig,ax = plt.subplots(figsize=(5,5))

    ax.stairs(distrs[0],bbin,color='b')
    ax.plot(xbin,distrs[0],marker='s',markersize=3,color='b',linestyle="None",mfc="w",alpha=0.8,label=f"$\gamma = {np.round(strain[frames[0]-1],decimals=2)}$")

    ax.stairs(distrs[16],bbin,color='y')
    ax.plot(xbin,distrs[16],marker='s',markersize=3,color='y',linestyle="None",mfc="w",alpha=0.8,label=f"$\gamma = {np.round(strain[frames[12]-1],decimals=2)}$")

    ax.stairs(distrs[23],bbin,color='g')
    ax.plot(xbin,distrs[23],marker='s',markersize=3,color='g',linestyle="None",mfc="w",alpha=0.8,label=f"$\gamma = {np.round(strain[frames[23]-1],decimals=2)}$")

    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)


    ax.set_xlabel("$L_{h}$")
    ax.set_ylabel("$N_h$")

    ax.legend()

    fig.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/hole_line_distr.pdf',dpi=300,bbox_inches='tight')

    plt.show()

def divide_hists(h1,h2):
    h1_h2 = []
    for i in range(len(h1)):
        if h2[i]>0:
            h1_h2.append(h1[i]/h2[i])
        else:
            h1_h2.append(0)
    return h1_h2

#Probabilidades de huecos (lineas)
def figura8():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_hole_data" #pc escritorio
    dir2 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data"
    dir3 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/exp5/hole_analysis_results"

    distrs = load_json(dir,"hole_distr_line_ave")
    bbin = load_json(dir,"binline")

    xbin = []
    for i in range(1,len(bbin)):
        xbin.append((bbin[i]+bbin[i-1])/2)
    xbin = np.array(xbin)

    strain = load_json(dir2,"strain")
    for i,x in enumerate(strain):
        if x < 0:
            strain[i] = x+1
    frames = [int(i) for i in load_json(dir3,"hole_analysis_line_calc_frames")]
    sample_distr = load_json(dir,"sample_distr_line_ave")

    fig,ax = plt.subplots(2,1,figsize=(5,5))

    ax[0].stairs(divide_hists(distrs[0],sample_distr[0]),bbin,color="b")
    ax[0].plot(xbin,divide_hists(distrs[0],sample_distr[0]),linestyle="None",marker="s",markersize=3,alpha=0.8,color="b",mfc="w",label=f"$\gamma = {np.round(strain[frames[0]-1],decimals=2)}$")

    ax[0].stairs(divide_hists(distrs[16],sample_distr[16]),bbin,color="y")
    ax[0].plot(xbin,divide_hists(distrs[16],sample_distr[16]),linestyle="None",marker="s",markersize=3,alpha=0.8,color="y",mfc="w",label=f"$\gamma = {np.round(strain[frames[12]-1],decimals=2)}$")

    ax[0].stairs(divide_hists(distrs[23],sample_distr[23]),bbin,color="g")
    ax[0].plot(xbin,divide_hists(distrs[23],sample_distr[23]),linestyle="None",marker="s",markersize=3,alpha=0.8,color="g",mfc="w",label=f"$\gamma = {np.round(strain[frames[23]-1],decimals=2)}$")

    ax[0].set_ylim([0,1])
    
    #ax[0].set_ylabel("$P(L)$")
    ax[0].minorticks_on()
    ax[0].tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax[0].tick_params(direction='in',which='major',length=5,right=True,top=True)
    ax[0].tick_params(labelbottom=False)

    ax[0].legend()

    ax[1].stairs(divide_hists(distrs[0],sample_distr[0]),bbin,color="b")
    ax[1].plot(xbin,divide_hists(distrs[0],sample_distr[0]),linestyle="None",marker="s",markersize=3,alpha=0.8,color="b",mfc="w",label=f"$\gamma = {np.round(strain[frames[0]-1],decimals=2)}$")

    ax[1].stairs(divide_hists(distrs[16],sample_distr[16]),bbin,color="y")
    ax[1].plot(xbin,divide_hists(distrs[16],sample_distr[16]),linestyle="None",marker="s",markersize=3,alpha=0.8,color="y",mfc="w",label=f"$\gamma = {np.round(strain[frames[12]-1],decimals=2)}$")

    ax[1].stairs(divide_hists(distrs[23],sample_distr[23]),bbin,color="g")
    ax[1].plot(xbin,divide_hists(distrs[23],sample_distr[23]),linestyle="None",marker="s",markersize=3,alpha=0.8,color="g",mfc="w",label=f"$\gamma = {np.round(strain[frames[23]-1],decimals=2)}$")

    ax[1].set_yscale("log")
    ax[1].set_ylim([0,1])
    ax[1].set_xlabel("$L_h$")
    ax[1].minorticks_on()
    ax[1].tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax[1].tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.ylabel("$P(L)$")

    fig.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/hole_prob_line.pdf',dpi=300,bbox_inches='tight')

    plt.show()

#Ajuste de probabilidades de huecos a modelo exponencial
def ajuste8():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 
    from matplotlib.ticker import NullFormatter

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_hole_data" #pc escritorio
    dir2 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data"
    dir3 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/exp5/hole_analysis_results"

    distrs = load_json(dir,"hole_distr_line_ave")
    bbin = load_json(dir,"binline")
    strain = load_json(dir2,"strain")
    for i in range(len(strain)-1):
        dif = strain[i+1]-strain[i]
        if dif < 0:
            strain[i+1:] = [x+1 for x in strain[i+1:]]
       
    frames = [int(i) for i in load_json(dir3,"hole_analysis_line_calc_frames")]
    sample_distr = load_json(dir,"sample_distr_line_ave")

    N_p = 16000
    V_p = (4/3)*np.pi*(2**(-5/6))**3
    V_c = (2*25.54)**3
    P_0 = 1-N_p*(V_p)/(V_c)
    print(P_0)

    xbin = []
    for i in range(1,len(bbin)):
        xbin.append((bbin[i]+bbin[i-1])/2)
    xbin = np.array(xbin)

    a = []
    rsq = []
    erco = []
    for i in range(int(len(distrs))):
        pr = np.array(divide_hists(distrs[i],sample_distr[i]))
        for j in range(len(pr)):
            if pr[j] == 0:
                iz = j
                break
        
        #Calcular a
        #f = -np.log(pr[:iz])
        f = np.log(np.divide(P_0,pr[:iz]))
        aa = np.dot(f[:iz],xbin[:iz])/np.dot(xbin[:iz],xbin[:iz])
        a.append(aa)

        #Calcular rsq
        #fit = np.exp(-aa*xbin)
        fit = P_0*np.exp(-aa*xbin)
        ssres = np.sum(np.square(pr-fit))
        sstot = np.sum(np.square(pr-np.mean(pr)))
        rsq.append(1-(ssres/sstot))
        #rsq.append(sstot)

        #error en cola
        Lco = 60
        for j,x in enumerate(xbin):
            if x > Lco:
                coix = j
                break
        erco.append(np.sum(pr[coix:iz]-fit[coix:iz])/(iz-coix))
        

    fig1,ax1 = plt.subplots(figsize=(5,5))
    calcstrain = [strain[i-1] for k,i in enumerate(frames) if k < len(a)]
    ax1.plot(calcstrain,a,marker = 's',linestyle="None")
    ax1.set_xlabel("$\gamma$")
    ax1.set_ylabel("$a(\gamma)$")
    ax1.minorticks_on()
    ax1.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax1.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #fig1.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/a(gamma).pdf',dpi=300,bbox_inches='tight')

    fig2,ax2 = plt.subplots(figsize=(5,5))

    #pf0 = 0
    pf0=0
    ax2.plot(xbin,divide_hists(distrs[pf0],sample_distr[pf0]), label=f"$\gamma = {np.round(strain[frames[pf0]-1],decimals=2)}$",color='b',linestyle="None",marker="s",alpha=0.5,markersize=3)
    ax2.plot(xbin,[P_0*np.e**(-a[pf0]*x) for x in xbin],color="b")

    #pf1 = 16
    pf1 = 30
    ax2.plot(xbin,divide_hists(distrs[pf1],sample_distr[pf1]), label=f"$\gamma = {np.round(strain[frames[pf1]-1],decimals=2)}$",color='y',linestyle="None",marker="s",alpha=0.5,markersize=3)
    ax2.plot(xbin,[P_0*np.e**(-a[pf1]*x) for x in xbin],color="y")

    #pf2 = 24
    pf2 = 46
    ax2.plot(xbin,divide_hists(distrs[pf2],sample_distr[pf2]), label=f"$\gamma = {np.round(strain[frames[pf2]-1],decimals=2)}$",color='g',linestyle="None",marker="s",alpha=0.5,markersize=3)
    ax2.plot(xbin,[P_0*np.e**(-a[pf2]*x) for x in xbin],color="g")

    ax2.minorticks_on()
    ax2.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax2.tick_params(direction='in',which='major',length=5,right=True,top=True)
    ax2.set_xlabel("$L_{h}$")
    ax2.set_ylabel("$P(L_{h})$")
    ax2.set_yscale("log")
    ax2.legend()

    x1,x2,y1,y2 = 30, 40, 9*10**(-3), 5*10**(-2)
    axins = ax2.inset_axes(
        [0.1, 0.1, 0.4, 0.3],
        xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
    axins.plot(xbin,divide_hists(distrs[pf0],sample_distr[pf0]),color='b',linestyle="None",marker="s")
    axins.plot(xbin,[P_0*np.e**(-a[pf0]*x) for x in xbin],color="b")
    axins.plot(xbin,divide_hists(distrs[pf1],sample_distr[pf1]),color='y',linestyle="None",marker="s")
    axins.plot(xbin,[P_0*np.e**(-a[pf1]*x) for x in xbin],color="y")
    axins.plot(xbin,divide_hists(distrs[pf2],sample_distr[pf2]),color='g',linestyle="None",marker="s")
    axins.plot(xbin,[P_0*np.e**(-a[pf2]*x) for x in xbin],color="g")
    #axins.minorticks_on()
    axins.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    axins.tick_params(direction='in',which='major',length=5,right=True,top=True)
    axins.set_yscale("log")
    axins.yaxis.set_major_formatter(NullFormatter())
    axins.yaxis.set_minor_formatter(NullFormatter())
    
    ax2.indicate_inset_zoom(axins, edgecolor="black")

    #fig2.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/prob_line_comp_ajuste.pdf',dpi=300,bbox_inches='tight')
    fig2.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/prob_line_comp_ajuste_hasta2.pdf',dpi=300,bbox_inches='tight')

    fig3,ax3 = plt.subplots(figsize=(5,5))
    ax3.plot(calcstrain,rsq,marker = 's',linestyle="None")
    ax3.set_xlabel("$\gamma$")
    ax3.set_ylabel("$R^2$")
    ax3.minorticks_on()
    ax3.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax3.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #fig3.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/rsq_a.pdf',dpi=300,bbox_inches='tight')

    fig4,ax4 = plt.subplots(figsize=(5,5))
    ax4.plot(calcstrain,erco,marker = 's',linestyle="--")
    ax4.set_xlabel("$\gamma$")
    ax4.set_ylabel("$MSD(L>60)$")
    ax4.minorticks_on()
    ax4.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax4.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig4.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/msd_l60.pdf',dpi=300,bbox_inches='tight')

    plt.show()

def demixing_param():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio 
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop

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


    phi2 = load_json(dir,"demixing_param_2_ave")
    phi2_std = load_json(dir,"demixing_param_2_std")
    phi4 = load_json(dir,"demixing_param_4_ave")
    phi4_std = load_json(dir,"demixing_param_4_std")
    phi5 = load_json(dir,"demixing_param_5_ave")
    phi5_std = load_json(dir,"demixing_param_5_std")
    phi7 = load_json(dir,"demixing_param_7_ave")
    phi7_std = load_json(dir,"demixing_param_7_std")
    phi10 = load_json(dir,"demixing_param_10_ave")
    phi10_std = load_json(dir,"demixing_param_10_std")

    #Crear subplot1 y graficar
    subplot1 = fig.add_subplot(gs[0,0:4])
    #subplot1.plot(np.linspace(0,2,50),phi2[:50],marker='s',markersize=2,mfc='w',linestyle='-')
    subplot1.errorbar(np.linspace(0,2,50),phi2[:50],yerr=phi2_std[:50],marker='s',markersize=3,mfc='w',linestyle='-', capsize=2)
    subplot1.axvline(x=yield_strain,linestyle='--',c='black')
    subplot1.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    subplot1.set_ylabel('$\psi_{2}$')

    #Limits
    subplot1.set_ylim((-0.0002,0.0042))
 
    #Ticks
    subplot1.set_xticks(np.arange(0,2.1,0.4))
    subplot1.set_yticks(np.arange(0,0.0041,0.001))
    
    subplot1.tick_params(direction='in',length=5,right=True,top=True)
    subplot1.minorticks_on()
    subplot1.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    subplot1.tick_params(direction='in',which='major',length=5,right=True,top=True)
    
    #Crear subplot2 y graficar
    subplot2 = fig.add_subplot(gs[1,0:4])
    #subplot2.plot(np.linspace(0,2,50),phi4[:50],marker='s',markersize=2,mfc='w',linestyle='-')
    subplot2.errorbar(np.linspace(0,2,50),phi4[:50],yerr=phi4_std[:50],marker='s',markersize=3,mfc='w',linestyle='-', capsize=2)
    subplot2.axvline(x=yield_strain,linestyle='--',c='black')
    subplot2.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    subplot2.set_ylabel('$\psi_{4}$')

    #Limits
    subplot2.set_ylim((-0.005,0.12))

    #Ticks
    subplot2.set_xticks(np.arange(0,2.1,0.4))
    subplot2.set_yticks(np.arange(0,0.12,0.03))
    #subplot2.tick_params(direction='in',length=5,right=True,top=True)
    subplot2.minorticks_on()
    subplot2.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    subplot2.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot3 y graficar
    subplot3 = fig.add_subplot(gs[2,0:4])
    #subplot3.plot(np.linspace(0,2,50),phi5[:50],marker='s',markersize=2,mfc='w',linestyle='-')
    subplot3.errorbar(np.linspace(0,2,50),phi5[:50],yerr=phi5_std[:50],marker='s',markersize=3,mfc='w',linestyle='-', capsize=2)
    subplot3.axvline(x=yield_strain,linestyle='--',c='black')
    subplot3.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    subplot3.set_ylabel('$\psi_{5}$')
    
    #Limites
    subplot3.set_ylim((-0.01,0.35))

    #Ticks
    subplot3.set_xticks(np.arange(0,2.1,0.4))
    subplot3.set_yticks(np.arange(0,0.31,0.1))
    #subplot2.tick_params(direction='in',length=5,right=True,top=True)
    subplot3.minorticks_on()
    subplot3.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    subplot3.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot4 y graficar
    subplot4 = fig.add_subplot(gs[3,0:4])
    #subplot4.plot(np.linspace(0,2,50),phi7[:50],marker='s',markersize=2,mfc='w',linestyle='-')
    subplot4.errorbar(np.linspace(0,2,50),phi7[:50],yerr=phi7_std[:50],marker='s',markersize=3,mfc='w',linestyle='-', capsize=2)
    subplot4.axvline(x=yield_strain,linestyle='--',c='black')
    subplot4.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    subplot4.set_ylabel('$\psi_{7}$')
    
    #Limites
    subplot4.set_ylim((-0.05,1.5))

    #Ticks
    subplot4.set_xticks(np.arange(0,2.1,0.4))
    subplot4.set_yticks(np.arange(0,1.31,0.4))
    #subplot2.tick_params(direction='in',length=5,right=True,top=True)
    subplot4.minorticks_on()
    subplot4.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    subplot4.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot5 y graficar
    subplot4 = fig.add_subplot(gs[4,0:4])
    #subplot4.plot(np.linspace(0,2,50),phi10[:50],marker='s',markersize=2,mfc='w',linestyle='-')
    subplot4.errorbar(np.linspace(0,2,50),phi10[:50],yerr=phi10_std[:50],marker='s',markersize=3,mfc='w',linestyle='-', capsize=2)
    subplot4.axvline(x=yield_strain,linestyle='--',c='black')
    subplot4.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    subplot4.set_ylabel('$\psi_{10}$')
    subplot4.set_xlabel('$\gamma$')
    
    #Limites
    subplot4.set_ylim((1,7))

    #Ticks
    subplot4.set_xticks(np.arange(0,2.1,0.4))
    subplot4.set_yticks(np.arange(1,6.1,2))
    #subplot2.tick_params(direction='in',length=5,right=True,top=True)
    subplot4.minorticks_on()
    subplot4.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    subplot4.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/demixing_parameter.pdf',dpi=300,bbox_inches='tight')

    plt.show()

def fractal_dimension():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio 
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop

    dir2 = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/exp1/analysis_results"

    #importar datos de estres y strain promediados
    stress_v_t_ave = [-i for i in load_json(dir,"stress_v_t_ave")]
    strain = load_json(dir,"strain")
    for i in range(len(strain)-1):
        dif = strain[i+1]-strain[i]
        if dif < 0:
            strain[i+1:] = [x+1 for x in strain[i+1:]]
    xstrain = strain.copy()
    strain.pop(0) #Quitar el frame 0 del strain
    #Obtener el yield strain. Lugar de maximo estres
    yield_ix = stress_v_t_ave.index(max(stress_v_t_ave))
    yield_strain = strain[yield_ix]

    #derivada del estres con respecto al strain
    dsigma_dgamma = derivative(stress_v_t_ave,0.01)
    maxderix = dsigma_dgamma.index(max(dsigma_dgamma))
    max_der_strain = strain[maxderix+1]

    #Importar frames calculados
    calc_frames = [int(i-1) for i in load_json(dir2,"calc_frames_fractal_dimension")]
    print(calc_frames)
    calc_strain = []
    for i in calc_frames:
        calc_strain.append(xstrain[i])

    #Crear figura
    fig,ax = plt.subplots(figsize=(5,5))

    D_ave = load_json(dir,"D_ave")
    D_std = load_json(dir,"D_std")

    ax.errorbar(calc_strain,D_ave,yerr=D_std,marker='s',markersize=3,mfc='w',linestyle='-', capsize=2)
    ax.axvline(x=yield_strain,linestyle='--',c='black')
    ax.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    ax.set_ylabel('$D$')
    ax.set_xlabel('$\gamma$')
    
    #Limites

    #Ticks
    #ax.set_xticks(np.arange(0,2.1,0.4))
    #ax.set_yticks(np.arange(0,1.31,0.4))
    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/Graficas/fractal_dimension.pdf',dpi=300,bbox_inches='tight')

    plt.show()

def distr_vecinos():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    mon = [1, 462, 14417, 0, 0, 0, 0, 0, 0, 0]
    xl = [0, 6, 102, 434, 578, 0, 0, 0, 0, 0]
    x = [0,1,2,3,4,5,6,7,8,9]

    fig1,ax1 = plt.subplots(figsize=(5,5))
    ax1.plot(x,mon,marker='s',markersize=8,mfc='b',alpha = 0.9,linestyle='None',label="Monomers")
    fig2,ax2 = plt.subplots(figsize=(5,5))
    ax2.plot(x,xl,marker='o',markersize=8,mfc='b',alpha = 0.9,linestyle='None',label="Cross-linkers")

    #Labels
    ax1.set_ylabel('Number of monomers')
    ax1.set_xlabel('Number of neighbors')
    ax2.set_ylabel('Number of cross-linkers')
    ax2.set_xlabel('Number of neighbors')
    #Limites

    #Ticks
    ax1.set_xticks(np.arange(0,9))
    ax2.set_xticks(np.arange(0,9))
    #ax.set_yticks(np.arange(0,1.31,0.4))
    #ax.minorticks_on()
    #ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax1.tick_params(direction='in',which='major',length=5,right=True,top=True)
    ax2.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig1.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SelfAssemby/Results/HistMon.pdf',dpi=300,bbox_inches='tight')
    fig2.savefig('/media/felipe/Files/Hydrogel_sim_experiments/SelfAssemby/Results/HistXl.pdf',dpi=300,bbox_inches='tight')

    plt.show()

def pot_ene_formation_se():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SelfAssemby/exp1_mp" #pc escritorio 
    dir = "E:/Hydrogel_sim_experiments/Test1/exp1"

    fig1,ax1 = plt.subplots(figsize=(5,5))

    pot_ene = [row[1] for row in load_data(dir,"pot_ene_formation_1.data")]
    frames = [row[0] for row in load_data(dir,"pot_ene_formation_1.data")]
    time = [float(i)*0.002 for i in frames]

    ax1.plot(time,pot_ene)

    #Labels
    ax1.set_ylabel('Potential Energy')
    ax1.set_xlabel('Time')

    #Ticks
    #ax1.set_xticks(np.arange(0,9))
    #ax.set_yticks(np.arange(0,1.31,0.4))
    ax1.minorticks_on()
    ax1.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax1.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig1.savefig("E:/Hydrogel_sim_experiments/Test1/exp1/Graficas/pot_ene_formation_1.pdf",dpi=300,bbox_inches='tight')

    plt.show()

#Grafica stress-strain con los datos de analisis comparados verticalmente
def stress_strain_basic_analysis_se():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio 
    #dir = "E:/Hydrogel_sim_experiments/Test1/exp1/analysis_results" #usb
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/exp4/analysis_results"

    #importar datos de estres y strain promediados
    stress_v_t = load_json(dir,"analysis_stress")
    strain = load_json(dir,"analysis_strain")
    strain.pop(0) #Quitar el frame 0 del strain
    #Obtener el yield strain. Lugar de maximo estres
    yield_ix = stress_v_t.index(max(stress_v_t))
    yield_strain = strain[yield_ix]

    #derivada del estres con respecto al strain
    dsigma_dgamma = derivative(stress_v_t,0.01)
    maxderix = dsigma_dgamma.index(max(dsigma_dgamma))
    max_der_strain = strain[maxderix+1]

    #Crear figura
    fig = plt.figure(1,figsize=(5,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(5,4)
    gs.update(wspace=0.2,hspace=0)

    #Crear subplot1 y graficar
    subplot1 = fig.add_subplot(gs[0,0:4])
    plt.plot(strain,stress_v_t)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    plt.ylabel('$\sigma_{xy}$')

    #Limites
    plt.xlim((0,10))
    plt.ylim((-0.01,0.1))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.yticks(np.arange(0,0.11,0.03))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #crear subplot 2 y graficar 
    subplot2 = fig.add_subplot(gs[1,0:4])
    plt.plot(strain[1:],dsigma_dgamma)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\frac{d \sigma_{xy}}{d \gamma}$")

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de bonds
    bonds_v_t = load_json(dir,"analysis_bonds_v_t")

    #crear subplot 3 y graficar
    subplot3 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain,bonds_v_t)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel('$N_b$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((21300,21700))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de clusters 
    clusternum_v_t = load_json(dir,"analysis_clusternum_v_t")

    #crear subplot 4 y graficar
    subplot4 = fig.add_subplot(gs[3,0:4])
    plt.plot(strain,clusternum_v_t)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel('$N_c$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((-2,28))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de tamaño del cluster mas grande
    bigclustersz_v_t = load_json(dir,"analysis_bigclustersz_v_t")

    #crear subplot 5 y graficar
    subplot5 = fig.add_subplot(gs[4,0:4])
    plt.plot(strain,bigclustersz_v_t)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel('Cluster size')
    plt.xlabel(r'$\gamma$')

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #fig.savefig("E:/Hydrogel_sim_experiments/Test1/exp1/Graficas/BasicAnalysis_1.pdf",dpi=300,bbox_inches='tight')
    
    plt.show()

#stress-strain log-log
def stress_strain_loglog_se():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio
    #dir = "D:/ExperimentData/SizeExperiment/16k_particles/averaged_data" #laptop
    dir = "E:/Hydrogel_sim_experiments/Test1/exp1/analysis_results" #usb

    sigma = load_json(dir,"analysis_stress")
    strain = load_json(dir,"analysis_strain")
    strain.pop(0)
    
    fig,ax = plt.subplots(figsize=(5,5))

    ax.plot(strain[:100],sigma[:100],marker='o',markersize=4,mfc='w',alpha = 0.9,linestyle='--')

    ax.set_xlabel("$\gamma$")
    ax.set_ylabel("$\sigma_{xy}$")

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig("E:/Hydrogel_sim_experiments/Test1/exp1/Graficas/stress_strain_loglog_1.pdf",dpi=300,bbox_inches='tight')
    plt.show()

#Grafica stress-strain con los datos de xldistance comparados verticalmente
def xldistance_se():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio 
    #dir = "E:/Hydrogel_sim_experiments/Test1/exp2/analysis_results" #usb
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/exp1/analysis_results"

    #importar datos de estres y strain promediados
    stress_v_t = load_json(dir,"analysis_stress")
    strain = load_json(dir,"analysis_strain")
    strain_nz = strain[1:]

    calcframes = load_json(dir,"xldistance_calcframes")
    c = 0
    strain_cf = []
    for i,s in enumerate(strain):
        if i == calcframes[c]-1:
            c+=1
            strain_cf.append(s)

    #Obtener el yield strain. Lugar de maximo estres
    yield_ix = stress_v_t.index(max(stress_v_t))
    yield_strain = strain[yield_ix]

    #derivada del estres con respecto al strain
    dsigma_dgamma = derivative(stress_v_t,0.01)
    maxderix = dsigma_dgamma.index(max(dsigma_dgamma))
    max_der_strain = strain[maxderix+1]

    #Crear figura
    fig = plt.figure(1,figsize=(5,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(5,4)
    gs.update(wspace=0.2,hspace=0)

    #Crear subplot1 y graficar
    subplot1 = fig.add_subplot(gs[0,0:4])
    plt.plot(strain_nz,stress_v_t)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    plt.ylabel('$\sigma_{xy}$')

    #Limites
    plt.xlim((0,10))
    plt.ylim((-0.01,0.1))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.yticks(np.arange(0,0.11,0.03))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #crear subplot 2 y graficar 
    subplot2 = fig.add_subplot(gs[1,0:4])
    plt.plot(strain_nz[1:],dsigma_dgamma)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\frac{d \sigma_{xy}}{d \gamma}$")

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de numero de xl_distance
    xldistance_v_t = load_json(dir,"mean_xl_distance_v_t")
    std_xldistance_v_t = load_json(dir,"std_xl_distance_v_t")

    #crear subplot 3 y graficar
    subplot3 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain_cf,xldistance_v_t,marker="s",markersize=3,mfc="w",linestyle="--")
    #plt.errorbar(strain_cf,xldistance_v_t,yerr=std_xldistance_v_t)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r'$\langle D_{xl} \rangle$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((21300,21700))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #crear subplot 3 y graficar
    subplot3 = fig.add_subplot(gs[3,0:4])
    plt.plot(strain_cf,std_xldistance_v_t,marker="s",markersize=3,mfc="w",linestyle="--")
    #plt.errorbar(strain_cf,xldistance_v_t,yerr=std_xldistance_v_t)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r'$\sigma\left(D_{xl}\right)$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((21300,21700))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de coordinacion promedio de cross linkers
    xl_coord_v_t = load_json(dir,"mean_xl_coord_v_t")

    #crear subplot 5 y graficar
    subplot5 = fig.add_subplot(gs[4,0:4])
    plt.plot(strain_cf,xl_coord_v_t,marker="s",markersize=3,mfc="w",linestyle="--")
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r'$\langle C_{xl} \rangle$')
    plt.xlabel("$\gamma$")

    #Limites
    plt.xlim((0,10))
    plt.ylim((1.9,3.1))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    plt.yticks(np.arange(2,3.1,0.5))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    #plt.tick_params(labelbottom=False)

    #fig.savefig("E:/Hydrogel_sim_experiments/Test1/exp1/Graficas/xldistance_1.pdf",dpi=300,bbox_inches='tight')
    
    plt.show()

def selected_xl_distances():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio
    #dir = "E:/Hydrogel_sim_experiments/Test1/exp1/analysis_results" #usb
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/exp1/analysis_results"

    selected_avdistances = load_json(dir,"selected_avdistances_v_t")
    strain = load_json(dir,"analysis_strain")
    calcframes = load_json(dir,"xldistance_calcframes")

    strain_cf = []
    c = 0
    for i,s in enumerate(strain):
        if i == calcframes[c]-1:
            c+=1
            strain_cf.append(s)

    
    
    fig,ax = plt.subplots(figsize=(5,5))

    for i in range(11,12):
        data = [row[i] for row in selected_avdistances]
        ax.plot(strain_cf,data,marker='o',markersize=4,mfc='w',alpha = 0.9,linestyle='--')

    ax.set_xlabel("$\gamma$")
    #ax.set_ylabel("$\sigma_{xy}$")

    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #fig.savefig("E:/Hydrogel_sim_experiments/Test1/exp1/Graficas/stress_strain_loglog_1.pdf",dpi=300,bbox_inches='tight')
    plt.show()

def chain_analysis_se():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio 
    #dir = "E:/Hydrogel_sim_experiments/Test1/exp2/analysis_results" #usb
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/exp1/analysis_results"

    #importar datos de estres y strain promediados
    stress_v_t = load_json(dir,"analysis_stress")
    strain = load_json(dir,"analysis_strain")
    strain_nz = strain[1:]

    calcframes = load_json(dir,"chain_analysis_calcframes")
    c = 0
    strain_cf = []
    for i,s in enumerate(strain):
        if i == calcframes[c]-1:
            c+=1
            strain_cf.append(s)

    #Obtener el yield strain. Lugar de maximo estres
    yield_ix = stress_v_t.index(max(stress_v_t))
    yield_strain = strain[yield_ix]

    #derivada del estres con respecto al strain
    dsigma_dgamma = derivative(stress_v_t,0.01)
    maxderix = dsigma_dgamma.index(max(dsigma_dgamma))
    max_der_strain = strain[maxderix+1]

    #Crear figura
    fig = plt.figure(1,figsize=(5,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(5,4)
    gs.update(wspace=0.2,hspace=0)

    #Crear subplot1 y graficar
    subplot1 = fig.add_subplot(gs[0,0:4])
    plt.plot(strain_nz,stress_v_t)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    plt.ylabel('$\sigma_{xy}$')

    #Limites
    plt.xlim((0,10))
    plt.ylim((-0.01,0.07))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.yticks(np.arange(0,0.11,0.03))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #crear subplot 2 y graficar 
    subplot2 = fig.add_subplot(gs[1,0:4])
    plt.plot(strain_nz[1:],dsigma_dgamma)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\frac{d \sigma_{xy}}{d \gamma}$")

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de curvatura
    curvature_v_t = load_json(dir,"mean_curvature_v_t")

    #crear subplot 3 y graficar
    subplot3 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain_cf,curvature_v_t,marker="s",markersize=3,mfc="w",linestyle="--")
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r'$\langle \kappa \rangle$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((21300,21700))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #Importar datos de numero de cadenas
    dangling_chains_num_v_t = load_json(dir,"dangling_chains_num_v_t")
    linked_chains_num_v_t = load_json(dir,"linked_chains_num_v_t")

    #crear subplot 4 y graficar
    subplot4 = fig.add_subplot(gs[3,0:4])
    plt.plot(strain_cf,dangling_chains_num_v_t,marker="s",markersize=3,mfc="w",linestyle="--",label="Dangling")
    plt.plot(strain_cf,linked_chains_num_v_t,marker="s",markersize=3,mfc="w",linestyle="--",label="Linked")
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')
    plt.legend()

    #Labels 
    plt.ylabel(r'$ N_{chains}$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((21300,21700))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de longitud de cadenas
    mean_dangling_chain_length_v_t = load_json(dir,"mean_dangling_chain_length_v_t")
    mean_linked_chain_length_v_t = load_json(dir,"mean_linked_chain_length_v_t")

    #crear subplot 5 y graficar
    subplot5 = fig.add_subplot(gs[4,0:4])
    plt.plot(strain_cf,mean_dangling_chain_length_v_t,marker="s",markersize=3,mfc="w",linestyle="--",label="Dangling")
    plt.plot(strain_cf,mean_linked_chain_length_v_t,marker="s",markersize=3,mfc="w",linestyle="--",label="Linked")
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')
    plt.legend()

    #Labels 
    plt.ylabel(r'$\langle L_{chain} \rangle$')
    plt.xlabel("$\gamma$")

    #Limites
    plt.xlim((0,10))
    #plt.ylim((1.9,3.1))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.yticks(np.arange(2,3.1,0.5))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    #plt.tick_params(labelbottom=False)

    #fig.savefig("E:/Hydrogel_sim_experiments/Test1/exp1/Graficas/xldistance_1.pdf",dpi=300,bbox_inches='tight')
    
    plt.show()

def chain_analysis():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    #Directorios con datos
    #dir = "/media/felipe/Files/Hydrogel_sim_experiments/SizeExperiments/16k_particles/averaged_data" #pc escritorio 
    #dir = "E:/Hydrogel_sim_experiments/Test1/exp2/analysis_results" #usb
    dir = "/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/averaged_data"

    #importar datos de estres y strain promediados
    stress_v_t_ave = load_json(dir,"analysis_stress_ave")
    stress_v_t_std = load_json(dir,"analysis_stress_std")
    strain = load_json(dir,"analysis_strain_ave")
    strain_nz = strain[1:]

    calcframes = load_json(dir,"chain_analysis_calcframes_ave")
    c = 0
    strain_cf = []
    for i,s in enumerate(strain):
        if i == calcframes[c]-1:
            c+=1
            strain_cf.append(s)

    #Obtener el yield strain. Lugar de maximo estres
    yield_ix = stress_v_t_ave.index(max(stress_v_t_ave))
    yield_strain = strain[yield_ix]

    #derivada del estres con respecto al strain
    dsigma_dgamma = derivative(stress_v_t_ave,0.01)
    maxderix = dsigma_dgamma.index(max(dsigma_dgamma))
    max_der_strain = strain[maxderix+1]
    max_der_strain = 0.5

    #Crear figura
    fig = plt.figure(1,figsize=(5,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(5,4)
    gs.update(wspace=0.2,hspace=0)

    #Crear subplot1 y graficar
    subplot1 = fig.add_subplot(gs[0,0:4])
    plt.plot(strain_nz,stress_v_t_ave)
    plt.errorbar(strain_nz[0::30],stress_v_t_ave[0::30],yerr=stress_v_t_std[0::30],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels
    plt.ylabel('$\sigma_{xy}$')

    #Limites
    plt.xlim((0,10))
    plt.ylim((-0.01,0.07))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.yticks(np.arange(0,0.11,0.03))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #crear subplot 2 y graficar 
    subplot2 = fig.add_subplot(gs[1,0:4])
    plt.plot(strain_nz[1:],dsigma_dgamma)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r"$\frac{d \sigma_{xy}}{d \gamma}$")

    #Limites
    plt.xlim((0,10))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de curvatura
    curvature_v_t_ave = load_json(dir,"mean_curvature_v_t_ave")
    curvature_v_t_std = load_json(dir,"mean_curvature_v_t_std")

    #crear subplot 3 y graficar
    subplot3 = fig.add_subplot(gs[2,0:4])
    plt.plot(strain_cf,curvature_v_t_ave,marker="s",markersize=3,mfc="w",linestyle="--")
    plt.errorbar(strain_cf[::2],curvature_v_t_ave[::2],yerr=curvature_v_t_std[::2],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)
    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')

    #Labels 
    plt.ylabel(r'$\langle \kappa \rangle$')

    #Limites
    plt.xlim((0,10))
    #plt.ylim((21300,21700))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #Importar datos de numero de cadenas
    dangling_chains_num_v_t_ave = load_json(dir,"dangling_chains_num_v_t_ave")
    dangling_chains_num_v_t_std = load_json(dir,"dangling_chains_num_v_t_std")
    linked_chains_num_v_t_ave = load_json(dir,"linked_chains_num_v_t_ave")
    linked_chains_num_v_t_std = load_json(dir,"linked_chains_num_v_t_std")

    #crear subplot 4 y graficar
    subplot4 = fig.add_subplot(gs[3,0:4])

    plt.plot(strain_cf,dangling_chains_num_v_t_ave,marker="s",markersize=3,mfc="w",linestyle="--",label="Dangling")
    plt.errorbar(strain_cf[::2],dangling_chains_num_v_t_ave[::2],yerr=dangling_chains_num_v_t_std[::2],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)

    plt.plot(strain_cf,linked_chains_num_v_t_ave,marker="s",markersize=3,mfc="w",linestyle="--",label="Linked")
    plt.errorbar(strain_cf[::2],linked_chains_num_v_t_ave[::2],yerr=linked_chains_num_v_t_std[::2],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)

    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')
    plt.legend(loc=4)

    #Labels 
    plt.ylabel(r'$ N_{chains}$')

    #Limites
    plt.xlim((0,10))
    plt.ylim((0,1800))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    plt.yticks(np.arange(200,1800,800))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    plt.tick_params(labelbottom=False)

    #importar datos de longitud de cadenas
    mean_dangling_chain_length_v_t_ave = load_json(dir,"mean_dangling_chain_length_v_t_ave")
    mean_dangling_chain_length_v_t_std = load_json(dir,"mean_dangling_chain_length_v_t_std")
    mean_linked_chain_length_v_t_ave = load_json(dir,"mean_linked_chain_length_v_t_ave")
    mean_linked_chain_length_v_t_std = load_json(dir,"mean_linked_chain_length_v_t_std")

    #crear subplot 5 y graficar
    subplot5 = fig.add_subplot(gs[4,0:4])
    plt.plot(strain_cf,mean_dangling_chain_length_v_t_ave,marker="s",markersize=3,mfc="w",linestyle="--",label="Dangling")
    plt.errorbar(strain_cf[::2],mean_dangling_chain_length_v_t_ave[::2],yerr=mean_dangling_chain_length_v_t_std[::2],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)

    plt.plot(strain_cf,mean_linked_chain_length_v_t_ave,marker="s",markersize=3,mfc="w",linestyle="--",label="Linked")
    plt.errorbar(strain_cf[::2],mean_linked_chain_length_v_t_ave[::2],yerr=mean_linked_chain_length_v_t_std[::2],marker='s',markersize=2,mfc='w',linestyle='None',capsize=2)

    plt.axvline(x=yield_strain,linestyle='--',c='black')
    plt.axvline(x=max_der_strain,linestyle='--',c='red')
    plt.legend()

    #Labels 
    plt.ylabel(r'$\langle L_{chain} \rangle$')
    plt.xlabel("$\gamma$")

    #Limites
    plt.xlim((0,10))
    plt.ylim((5.7,8.2))

    #Ticks
    #plt.xticks(np.arange(0,1.1,0.2))
    #plt.yticks(np.arange(2,3.1,0.5))
    #plt.tick_params(direction='in',length=5,right=True,top=True)
    plt.minorticks_on()
    plt.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    plt.tick_params(direction='in',which='major',length=5,right=True,top=True)
    #plt.tick_params(labelbottom=False)

    fig.savefig("/media/felipe/KINGSTON/Hydrogel_sim_experiments/Test1/Graficas/chain_analysis.pdf",dpi=300,bbox_inches='tight')
    
    plt.show()

def main():

    #stress_strain_basic_analysis()
    #stress_strain_basic_analysis_30k()
    #stress_strain_xl_distance()
    #stress_strain_xl_distance_30k()
    #stress_strain_loglog()
    #demixing_param()
    #fractal_dimension()
    #pot_ene_formation()
    #distr_vecinos()
    #pot_ene_formation_se()
    #xldistance_se()
    #selected_xl_distances()
    #chain_analysis_se()
    chain_analysis()

if __name__ == '__main__':
    main()