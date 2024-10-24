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

def stress_strain():  
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    shear_rates = [r"$\dot{\gamma}=1\times 10^{-2}$",r"$\dot{\gamma}=1\times 10^{-3}$",r"$\dot{\gamma}=1\times 10^{-4}$"]
    colors = ["#1f77b4","#ff7f0e","#2ca02c"]

    stress_ave = []
    stress_err = []
    strain = []
    for i in range(3):
        dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/averaged_data"
        stress_ave.append(load_json(dir,"stress_ave"))
        stress_err.append(load_json(dir,"stress_err"))
        strain.append(load_json(dir,"strain_ave")[1:])
        

    fig,ax = plt.subplots(figsize=(5,5))

    for i in range(3):
        ax.plot(strain[i],stress_ave[i],label=shear_rates[i],color=colors[i])
        ax.errorbar(strain[i][::30],stress_ave[i][::30],yerr=stress_err[i][::30],marker='o',markersize=4,mfc='w',alpha = 0.9,linestyle='None',capsize=1,color=colors[i])

    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$\sigma_{xy}$")

    ax.legend()

    #ax.set_xscale('log')
    #ax.set_yscale('log')

    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/stress_strain.pdf",dpi=300,bbox_inches='tight')

    plt.show()

def stress_strain_loglog():  
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    shear_rates = [r"$\dot{\gamma}=1\times 10^{-2}$",r"$\dot{\gamma}=1\times 10^{-3}$",r"$\dot{\gamma}=1\times 10^{-4}$"]
    colors = ["#1f77b4","#ff7f0e","#2ca02c"]

    stress_ave = []
    stress_err = []
    strain = []
    for i in range(3):
        dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/averaged_data"
        stress_ave.append(load_json(dir,"stress_ave"))
        stress_err.append(load_json(dir,"stress_err"))
        strain.append(load_json(dir,"strain_ave")[1:])
        

    fig,ax = plt.subplots(figsize=(5,5))

    for i in range(3):
        ax.plot(strain[i],stress_ave[i],label=shear_rates[i],color=colors[i])
        ax.errorbar(strain[i][::30],stress_ave[i][::30],yerr=stress_err[i][::30],marker='o',markersize=4,mfc='w',alpha = 0.9,linestyle='None',capsize=1,color=colors[i])

    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$\sigma_{xy}$")

    ax.legend()

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/stress_strain_loglog.pdf",dpi=300,bbox_inches='tight')

    plt.show()

def basic_analysis():

    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    #Crear figura
    fig = plt.figure(1,figsize=(8,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(5,4)
    gs.update(wspace=0.2,hspace=0)

    shear_rates = [r"$\dot{\gamma}=1\times 10^{-2}$",r"$\dot{\gamma}=1\times 10^{-3}$",r"$\dot{\gamma}=1\times 10^{-4}$"]
    colors = ["#1f77b4","#ff7f0e","#2ca02c"]

    #Importar datos
    stress_ave = []
    stress_err = []
    strain = []
    ddtstress = []

    bondnum_ave=[]
    bondnum_err=[]

    clusternum_ave=[]
    clusternum_err=[]

    bigclustersz_ave=[]
    bigclustersz_err=[]

    for i in range(3):
        dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/averaged_data"
        stress_ave.append(load_json(dir,"stress_ave"))
        stress_err.append(load_json(dir,"stress_err"))
        strain.append(load_json(dir,"strain_ave"))
        ddtstress.append(derivative(load_json(dir,"stress_ave"),0.01))
        bondnum_ave.append(load_json(dir,"bondnum_ave"))
        bondnum_err.append(load_json(dir,"bondnum_err"))
        clusternum_ave.append(load_json(dir,"clusternum_ave"))
        clusternum_err.append(load_json(dir,"clusternum_err"))
        bigclustersz_ave.append(load_json(dir,"bigclustersz_ave"))
        bigclustersz_err.append(load_json(dir,"bigclustersz_err"))

    #Crear subplot
    ax1 = fig.add_subplot(gs[0,0:4])

    for i in range(3):
        ax1.plot(strain[i][1:],stress_ave[i],label=shear_rates[i],color=colors[i])

    ax1.legend()
    ax1.set_ylabel(r"$\sigma_{xy}$")

    ax1.minorticks_on()
    ax1.set_xticks(np.arange(0,20.1,5))
    ax1.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax1.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax2 = fig.add_subplot(gs[1,0:4])

    for i in range(3):
        ax2.plot(strain[i][2:],ddtstress[i],label=shear_rates[i],color=colors[i])

    ax2.set_ylabel(r"$\frac{d}{d\gamma}\sigma_{xy}$")

    ax2.minorticks_on()
    ax2.set_xticks(np.arange(0,20.1,5))
    ax2.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax2.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax3 = fig.add_subplot(gs[2,0:4])

    for i in range(3):
        ax3.plot(strain[i][1:],bondnum_ave[i],label=shear_rates[i],color=colors[i])
        ax3.errorbar(strain[i][1::60],bondnum_ave[i][::60],yerr=bondnum_err[i][::60],linestyle="None",capsize=2,color=colors[i])

    ax3.set_ylabel("$N_{bonds}$")

    ax3.minorticks_on()
    ax3.set_xticks(np.arange(0,20.1,5))
    ax3.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax3.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax4 = fig.add_subplot(gs[3,0:4])

    for i in range(3):
        ax4.plot(strain[i][1:],clusternum_ave[i],label=shear_rates[i],color=colors[i])
        ax4.errorbar(strain[i][1::60],clusternum_ave[i][::60],yerr=clusternum_err[i][::60],linestyle="None",capsize=2,color=colors[i])

    ax4.set_ylabel("$N_{clusters}$")

    ax4.minorticks_on()
    ax4.set_xticks(np.arange(0,20.1,5))
    ax4.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax4.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax4 = fig.add_subplot(gs[4,0:4])

    for i in range(3):
        ax4.plot(strain[i][1:],bigclustersz_ave[i],label=shear_rates[i],color=colors[i])
        ax4.errorbar(strain[i][1::60],bigclustersz_ave[i][::60],yerr=bigclustersz_err[i][::60],linestyle="None",capsize=2,color=colors[i])

    ax4.set_ylabel("$Max Cluster Size$")

    ax4.minorticks_on()
    ax4.set_xticks(np.arange(0,20.1,5))
    ax4.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax4.tick_params(direction='in',which='major',length=5,right=True,top=True)

    ax4.set_xlabel(r"$\gamma$")

    fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/basic_analysis.pdf",dpi=300,bbox_inches='tight')

    plt.show()

def xldistance():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    #Crear figura
    fig = plt.figure(1,figsize=(8,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(5,4)
    gs.update(wspace=0.2,hspace=0)

    shear_rates = [r"$\dot{\gamma}=1\times 10^{-2}$",r"$\dot{\gamma}=1\times 10^{-3}$",r"$\dot{\gamma}=1\times 10^{-4}$"]
    colors = ["#1f77b4","#ff7f0e","#2ca02c"]

    #Importar datos
    stress_ave = []
    stress_err = []
    strain = []
    ddtstress = []

    mean_xl_coordination_ave=[]
    mean_xl_coordination_err=[]

    std_xl_coordination_ave=[]
    std_xl_coordination_err=[]

    mean_xl_distance_ave=[]
    mean_xl_distance_err=[]

    std_xl_distance_ave=[]
    std_xl_distance_err=[]

    xldistance_calcframes=[]

    for i in range(3):
        dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/averaged_data"

        stress_ave.append(load_json(dir,"stress_ave"))
        stress_err.append(load_json(dir,"stress_err"))
        strain.append(load_json(dir,"strain_ave"))
        ddtstress.append(derivative(load_json(dir,"stress_ave"),0.01))

        mean_xl_coordination_ave.append(load_json(dir,"mean_xl_coordination_ave"))
        mean_xl_coordination_err.append(load_json(dir,"mean_xl_coordination_err"))

        std_xl_coordination_ave.append(load_json(dir,"std_xl_coordination_ave"))
        std_xl_coordination_err.append(load_json(dir,"std_xl_coordination_err"))

        mean_xl_distance_ave.append(load_json(dir,"mean_xl_distance_ave"))
        mean_xl_distance_err.append(load_json(dir,"mean_xl_distance_err"))

        std_xl_distance_ave.append(load_json(dir,"std_xl_distance_ave"))
        std_xl_distance_err.append(load_json(dir,"std_xl_distance_err"))

        xldistance_calcframes.append(load_json(dir,"xldistance_calcframes_ave"))
    
    strain_cf=[]
    for i in range(3):
        ss_cf=[]
        for k in xldistance_calcframes[i]:
            ss_cf.append(strain[i][int(k)-1])    
        strain_cf.append(ss_cf)


    #Crear subplot
    ax1 = fig.add_subplot(gs[0,0:4])

    for i in range(3):
        ax1.plot(strain[i][1:],stress_ave[i],label=shear_rates[i],color=colors[i])

    ax1.legend()
    ax1.set_ylabel(r"$\sigma_{xy}$")

    ax1.minorticks_on()
    ax1.set_xticks(np.arange(0,20.1,5))
    ax1.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax1.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax2 = fig.add_subplot(gs[1,0:4])

    for i in range(3):
        ax2.plot(strain_cf[i],mean_xl_coordination_ave[i],label=shear_rates[i],color=colors[i])
        ax2.errorbar(strain_cf[i][::10],mean_xl_coordination_ave[i][::10],yerr=mean_xl_coordination_err[i][::10],linestyle="None",capsize=2,color=colors[i])

    ax2.set_ylabel(r"$\langle C_{xl} \rangle$")

    ax2.minorticks_on()
    ax2.set_xticks(np.arange(0,20.1,5))
    ax2.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax2.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax3 = fig.add_subplot(gs[2,0:4])

    for i in range(3):
        ax3.plot(strain_cf[i],std_xl_coordination_ave[i],label=shear_rates[i],color=colors[i])
        ax3.errorbar(strain_cf[i][::10],std_xl_coordination_ave[i][::10],yerr=std_xl_coordination_err[i][::10],linestyle="None",capsize=2,color=colors[i])

    ax3.set_ylabel(r"$\sigma(C_{xl})$")

    ax3.minorticks_on()
    ax3.set_xticks(np.arange(0,20.1,5))
    ax3.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax3.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax4 = fig.add_subplot(gs[3,0:4])

    for i in range(3):
        ax4.plot(strain_cf[i],mean_xl_distance_ave[i],label=shear_rates[i],color=colors[i])
        ax4.errorbar(strain_cf[i][::10],mean_xl_distance_ave[i][::10],yerr=mean_xl_distance_err[i][::10],linestyle="None",capsize=2,color=colors[i])

    ax4.set_ylabel(r"$\langle L_{xl} \rangle$")

    ax4.minorticks_on()
    ax4.set_xticks(np.arange(0,20.1,5))
    ax4.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax4.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax5 = fig.add_subplot(gs[4,0:4])

    for i in range(3):
        ax5.plot(strain_cf[i],std_xl_distance_ave[i],label=shear_rates[i],color=colors[i])
        ax5.errorbar(strain_cf[i][::10],std_xl_distance_ave[i][::10],yerr=std_xl_distance_err[i][::10],linestyle="None",capsize=2,color=colors[i])

    ax5.set_ylabel(r"$\sigma(L_{xl})$")
    ax5.set_xlabel(r"$\gamma$")

    ax5.minorticks_on()
    ax5.set_xticks(np.arange(0,20.1,5))
    ax5.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax5.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/xldistance.pdf",dpi=300,bbox_inches='tight')

    plt.show()

def chain_analysis():

    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    #Crear figura
    fig = plt.figure(1,figsize=(8,10))

    #Gridspec para subplots
    gs = gridspec.GridSpec(6,4)
    gs.update(wspace=0.2,hspace=0)

    shear_rates = [r"$\dot{\gamma}=1\times 10^{-2}$",r"$\dot{\gamma}=1\times 10^{-3}$",r"$\dot{\gamma}=1\times 10^{-4}$"]
    colors = ["#1f77b4","#ff7f0e","#2ca02c"]

    #Importar datos
    stress_ave = []
    stress_err = []
    strain = []
    ddtstress = []

    dangling_chains_num_ave=[]
    dangling_chains_num_err=[]

    linked_chains_num_ave=[]
    linked_chains_num_err=[]

    mean_curvature_ave=[]
    mean_curvature_err=[]

    mean_total_chain_length_ave=[]
    mean_total_chain_length_err=[]

    std_total_chain_length_ave=[]
    std_total_chain_length_err=[]

    chain_analysis_calcframes=[]

    for i in range(3):
        dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/averaged_data"

        stress_ave.append(load_json(dir,"stress_ave"))
        stress_err.append(load_json(dir,"stress_err"))
        strain.append(load_json(dir,"strain_ave"))
        ddtstress.append(derivative(load_json(dir,"stress_ave"),0.01))

        dangling_chains_num_ave.append(load_json(dir,"dangling_chains_num_ave"))
        dangling_chains_num_err.append(load_json(dir,"dangling_chains_num_err"))

        linked_chains_num_ave.append(load_json(dir,"linked_chains_num_ave"))
        linked_chains_num_err.append(load_json(dir,"linked_chains_num_err"))

        mean_curvature_ave.append(load_json(dir,"mean_curvature_ave"))
        mean_curvature_err.append(load_json(dir,"mean_curvature_err"))

        mean_total_chain_length_ave.append(load_json(dir,"mean_total_chain_length_ave"))
        mean_total_chain_length_err.append(load_json(dir,"mean_total_chain_length_err"))

        std_total_chain_length_ave.append(load_json(dir,"std_total_chain_length_ave"))
        std_total_chain_length_err.append(load_json(dir,"std_total_chain_length_err"))

        chain_analysis_calcframes.append(load_json(dir,"chain_analysis_calcframes_ave"))
    
    strain_cf=[]
    for i in range(3):
        ss_cf=[]
        for k in chain_analysis_calcframes[i]:
            ss_cf.append(strain[i][int(k)-1])    
        strain_cf.append(ss_cf)


    #Crear subplot
    ax1 = fig.add_subplot(gs[0,0:4])

    for i in range(3):
        ax1.plot(strain[i][1:],stress_ave[i],label=shear_rates[i],color=colors[i])

    ax1.legend(loc="right")
    ax1.set_ylabel(r"$\sigma_{xy}$")

    ax1.minorticks_on()
    ax1.set_xticks(np.arange(0,20.1,5))
    ax1.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax1.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax2 = fig.add_subplot(gs[1,0:4])

    for i in range(3):
        ax2.plot(strain_cf[i],dangling_chains_num_ave[i],label=shear_rates[i],color=colors[i])
        ax2.errorbar(strain_cf[i][::10],dangling_chains_num_ave[i][::10],yerr=dangling_chains_num_err[i][::10],linestyle="None",capsize=2,color=colors[i])

    ax2.set_ylabel(r"$N_{dangling}$")

    ax2.minorticks_on()
    ax2.set_xticks(np.arange(0,20.1,5))
    ax2.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax2.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax3 = fig.add_subplot(gs[2,0:4])

    for i in range(3):
        ax3.plot(strain_cf[i],linked_chains_num_ave[i],label=shear_rates[i],color=colors[i])
        ax3.errorbar(strain_cf[i][::10],linked_chains_num_ave[i][::10],yerr=linked_chains_num_err[i][::10],linestyle="None",capsize=2,color=colors[i])

    ax3.set_ylabel(r"$N_{linked}$")

    ax3.minorticks_on()
    ax3.set_xticks(np.arange(0,20.1,5))
    ax3.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax3.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax4 = fig.add_subplot(gs[3,0:4])

    for i in range(3):
        ax4.plot(strain_cf[i],mean_curvature_ave[i],label=shear_rates[i],color=colors[i])
        ax4.errorbar(strain_cf[i][::10],mean_curvature_ave[i][::10],yerr=mean_curvature_err[i][::10],linestyle="None",capsize=2,color=colors[i])

    ax4.set_ylabel(r"$\kappa$")

    ax4.set_ylim((2.2,4.3))

    ax4.minorticks_on()
    ax4.set_xticks(np.arange(0,20.1,5))
    ax4.set_yticks(np.arange(2.5,4.1,0.5))
    ax4.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax4.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax5 = fig.add_subplot(gs[4,0:4])

    for i in range(3):
        ax5.plot(strain_cf[i],mean_total_chain_length_ave[i],label=shear_rates[i],color=colors[i])
        ax5.errorbar(strain_cf[i][::10],mean_total_chain_length_ave[i][::10],yerr=mean_total_chain_length_err[i][::10],linestyle="None",capsize=2,color=colors[i])

    ax5.set_ylabel(r"$\langle L_{chain} \rangle$")

    ax5.minorticks_on()
    ax5.set_xticks(np.arange(0,20.1,5))
    ax5.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax5.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax6 = fig.add_subplot(gs[5,0:4])

    for i in range(3):
        ax6.plot(strain_cf[i],std_total_chain_length_ave[i],label=shear_rates[i],color=colors[i])
        ax6.errorbar(strain_cf[i][::10],std_total_chain_length_ave[i][::10],yerr=std_total_chain_length_err[i][::10],linestyle="None",capsize=2,color=colors[i])

    ax6.set_ylabel(r"$\sigma ( L_{chain} )$")
    ax6.set_xlabel(r"$\gamma$")

    ax6.minorticks_on()
    ax6.set_xticks(np.arange(0,20.1,5))
    ax6.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax6.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/chain_analysis.pdf",dpi=300,bbox_inches='tight')

    plt.show()
    
def stress_shearrate():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    shear_rates = [r"$\dot{\gamma}=1\times 10^{-2}$",r"$\dot{\gamma}=1\times 10^{-3}$",r"$\dot{\gamma}=1\times 10^{-4}$"]
    colors = ["#1f77b4","#ff7f0e","#2ca02c"]

    stress_ave = []
    stress_err = []

    for i in range(10):
        dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/averaged_data"
        stress_ave.append(load_json(dir,"stress_ave"))
        stress_err.append(load_json(dir,"stress_err"))
    
    strain = load_json("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma1/averaged_data","strain_ave")[1:]

    index=0
    for i,s in enumerate(strain):
        if s>=15:
            index=i
            break
    
    eq_stress = []
    eq_stress_err = []
    for i in range(10):
        eq_stress.append(np.mean(stress_ave[i][index:]))
        err = 0
        for k in stress_err[i][index:]:
            err+=k**2
        eq_stress_err.append((1/len(stress_err[i][index:]))*np.sqrt(err))

    shearrate = [0.01,0.001,0.0001,0.00125,0.0025,0.00375,0.005,0.00625,0.0075,0.00875]

    order = [2,1,3,4,5,6,7,8,9,0]
    shearrate_ordered = []
    eq_stress_ordered = []
    eq_stress_err_ordered = []
    for i in order:
        shearrate_ordered.append(shearrate[i])
        eq_stress_ordered.append(eq_stress[i])
        eq_stress_err_ordered.append(eq_stress_err[i])

    fig,ax = plt.subplots(figsize=(5,5))

    ax.plot(shearrate_ordered,eq_stress_ordered,linestyle="--")
    ax.errorbar(shearrate_ordered,eq_stress_ordered,yerr=eq_stress_err_ordered,marker='s',markersize=6,mfc='w',alpha = 0.9,linestyle='None',capsize=1)

    ax.set_xlabel(r"$\dot{\gamma}$")
    ax.set_ylabel(r"$\sigma_{xy}$")

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim(0,0.011)
    #ax.set_ylim(0,0.0125)

    ax.minorticks_on()
    #ax.set_yticks(np.arange(0.001,0.0121,0.002))
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/stress_shearrate.pdf",dpi=300,bbox_inches='tight')

    plt.show()

def nematic_order_parameter():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    #Crear figura
    fig = plt.figure(1,figsize=(8,8))

    #Gridspec para subplots
    gs = gridspec.GridSpec(2,4)
    gs.update(wspace=0.2,hspace=0)

    shear_rates = [r"$\dot{\gamma}=1\times 10^{-2}$",r"$\dot{\gamma}=1\times 10^{-3}$",r"$\dot{\gamma}=1\times 10^{-4}$"]
    colors = ["#1f77b4","#ff7f0e","#2ca02c"]

    #Importar datos
    stress_ave = []
    stress_err = []
    strain = []

    S_ave = []
    S_err = []

    calcframes=[]

    for i in range(3):
        dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/averaged_data"

        stress_ave.append(load_json(dir,"stress_ave"))
        stress_err.append(load_json(dir,"stress_err"))
        strain.append(load_json(dir,"strain_ave"))

        S_ave.append(load_json(dir,"nematic_order_parameter_ave"))
        S_err.append(load_json(dir,"nematic_order_parameter_err"))

        calcframes.append(load_json(dir,"nematic_order_parameter_calcframes_ave"))


    strain_cf=[]
    for i in range(3):
        ss_cf=[]
        for k in calcframes[i]:
            ss_cf.append(strain[i][int(k)-1])    
        strain_cf.append(ss_cf)


    #Crear subplot
    ax1 = fig.add_subplot(gs[0,0:4])

    for i in range(3):
        ax1.plot(strain[i][1:],stress_ave[i],label=shear_rates[i],color=colors[i])

    ax1.legend()
    ax1.set_ylabel(r"$\sigma_{xy}$")

    ax1.minorticks_on()
    ax1.set_xticks(np.arange(0,20.1,5))
    ax1.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax1.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #Crear subplot
    ax2 = fig.add_subplot(gs[1,0:4])

    for i in range(3):
        ax2.plot(strain_cf[i],S_ave[i],label=shear_rates[i],color=colors[i])
        ax2.errorbar(strain_cf[i][::10],S_ave[i][::10],yerr=S_err[i][::10],linestyle="None",capsize=2,color=colors[i])


    ax2.set_ylabel("$S$")
    ax2.set_xlabel(r"$\gamma$")

    ax2.set_ylim((0,0.62))

    ax2.minorticks_on()
    ax2.set_xticks(np.arange(0,20.1,5))
    ax2.set_yticks(np.arange(0,0.61,0.1))
    ax2.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax2.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/nematic_order_parameter.pdf",dpi=300,bbox_inches='tight')

    plt.show()

def extension():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    shear_rates = [r"$\dot{\gamma}=1\times 10^{-2}$",r"$\dot{\gamma}=1\times 10^{-3}$",r"$\dot{\gamma}=1\times 10^{-4}$"]
    colors = ["#1f77b4","#ff7f0e","#2ca02c"]

    stress_ave = []
    stress_err = []
    strain = []

    for i in range(3):
        dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/averaged_data"
        stress_ave.append(load_json(dir,"stress_ave"))
        stress_err.append(load_json(dir,"stress_err"))
        strain.append(load_json(dir,"strain_ave")[1:])

    stress_ext_ave = []
    stress_ext_err = []
    strain_ext = []    
    dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma1/averaged_data"
    stress_ext_ave.append(load_json(dir,"stress_extension1_ave"))
    stress_ext_err.append(load_json(dir,"stress_extension1_err"))
    strain_ext.append(load_json(dir,"strain_extension1_ave")[1:]) 

    for i,ss in enumerate(strain_ext):
        strain_ext[i]=[x+20 for x in ss]

    index=0
    for i,s in enumerate(strain[0]):
        if s>=15:
            index=i
            break
    
    eq_stress = []
    eq_stress_err = []
    for i in range(3):
        eq_stress.append(np.mean(stress_ave[i][index:]))
        err = 0
        for k in stress_err[i][index:]:
            err+=k**2
        eq_stress_err.append((1/len(stress_err[i][index:]))*np.sqrt(err))

    fig,ax = plt.subplots(figsize=(5,5))

    for i in range(2):
        ax.plot(strain[i],stress_ave[i],label=shear_rates[i],color=colors[i])
        ax.errorbar(strain[i][::30],stress_ave[i][::30],yerr=stress_err[i][::30],marker='o',markersize=4,mfc='w',alpha = 0.9,linestyle='None',capsize=1,color=colors[i])

    ax.plot(strain_ext[0],stress_ext_ave[0],color=colors[0])
    ax.errorbar(strain_ext[0][::30],stress_ext_ave[0][::30],yerr=stress_ext_err[0][::30],marker='o',markersize=4,mfc='w',alpha = 0.9,linestyle='None',capsize=1,color=colors[0])

    ax.hlines(eq_stress[1],15,26,linestyle="--",color="red")

    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$\sigma_{xy}$")

    ax.legend()

    #ax.set_xscale('log')
    #ax.set_yscale('log')

    ax.minorticks_on()
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/stress_strain_extension1.pdf",dpi=300,bbox_inches='tight')

    plt.show()

def structure_factor():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    shear_rates = [r"$\dot{\gamma}=1\times 10^{-2}$",r"$\dot{\gamma}=1\times 10^{-3}$",r"$\dot{\gamma}=1\times 10^{-4}$"]
    colors = ["#1f77b4","#ff7f0e","#2ca02c"]

    #Crear figura
    fig = plt.figure(1,figsize=(8,4))

    #Gridspec para subplots
    gs = gridspec.GridSpec(1,1)
    gs.update(wspace=0.2,hspace=0.5)

    #ax = fig.add_subplot(gs[0,0:4])

    #Importar datos 

    dir = f"E:/Hydrogel_sim_experiments/FullExperiment1/dGamma1/exp1/analysis_results"
    S_vt = load_json(dir,"structure_factor_vt_x")
    k = load_json(dir,"structure_factor_vt_kvalues_x")

    ax = fig.add_subplot(gs[0,0])
    for i in [0,4]:
        ax.plot(k,S_vt[i],marker="s",linestyle="--",label=f"{i}")

    ax.legend()

    #fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/velocity_profile_dGamma3_exp1.pdf",dpi=300,bbox_inches='tight')

    plt.show()

def ISF():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    deltat = [0.001*1000,0.001*10000,0.001*100000]
    titles = [r"$\dot{\gamma}=10^{-2}$",r"$\dot{\gamma}=10^{-3}$",r"$\dot{\gamma}=10^{-4}$"]

    ISF_xy=[]
    ISF_xz=[]
    ISF_yz=[]
    time_xy=[]
    time_xz=[]
    time_yz=[]

    for i in range(3):
        #dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/exp1/analysis_results"
        dir = f"E:/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/exp1/analysis_results"
        ISF_xy.append(load_json(dir,"ISF_xy"))
        frames_xy = load_json(dir,"ISF_calcframes_xy")
        time_xy.append([deltat[i]*j for j in frames_xy])

        ISF_xz.append(load_json(dir,"ISF_xz"))
        frames_xz = load_json(dir,"ISF_calcframes_xz")
        time_xz.append([deltat[i]*j for j in frames_xz])

        ISF_yz.append(load_json(dir,"ISF_yz"))
        frames_yz = load_json(dir,"ISF_calcframes_yz")
        time_yz.append([deltat[i]*j for j in frames_yz])


    #Crear figura
    fig = plt.figure(1,figsize=(8,10))

    #Gridspec para subplots
    gs = gridspec.GridSpec(3,4)
    gs.update(wspace=0.2,hspace=0.7)

    for i in range(3):

        ax = fig.add_subplot(gs[i,0:4])

        ax.set_title(titles[i])

        ax.plot(time_xy[i],ISF_xy[i],marker="s",markersize=2,label="xy")
        ax.plot(time_xz[i],ISF_xz[i],marker="s",markersize=2,label="xz")
        ax.plot(time_yz[i],ISF_yz[i],marker="s",markersize=2,label="yz")

        ax.set_xlabel("$t$")
        ax.set_ylabel("$F$")

        ax.legend()

        ax.set_xscale('log')
        #ax.set_yscale('log')

        #ax.set_xlim(0,0.011)
        #ax.set_ylim(0,0.0125)

        ax.minorticks_on()
        #ax.set_yticks(np.arange(0.001,0.0121,0.002))
        ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
        ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    #fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/stress_shearrate.pdf",dpi=300,bbox_inches='tight')

    plt.show()

def vel_distr():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    dir = f"E:/Hydrogel_sim_experiments/FullExperiment1/dGamma1/exp1/analysis_results"
    all_velocities = load_json(dir,"all_velocities")

    fig,ax = plt.subplots(figsize=(5,5))

    ax.hist(all_velocities,bins=np.arange(0,30,0.5),density=True)
    ax.set_yscale("log")

    plt.show()

def velocity_profile():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.size'] = 14

    dir = f"E:/Hydrogel_sim_experiments/FullExperiment1/dGamma1/exp1/analysis_results"
    velocity_profile = load_json(dir,"velocity_profile")
    velocity_profile_bins = load_json(dir,"velocity_profile_bins")

    fig,ax = plt.subplots(figsize=(5,5))
    
    for i in range(10):
        ax.plot(velocity_profile_bins,velocity_profile[i],marker="s")

    plt.show()


def main():

    #stress_strain()
    #stress_strain_loglog()
    #basic_analysis()
    #xldistance()
    #chain_analysis()
    #stress_shearrate()
    #nematic_order_parameter()
    #extension()
    #structure_factor()
    #ISF()
    #vel_distr()
    velocity_profile()

if __name__=="__main__":
    main()