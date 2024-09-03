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

    ax.set_xlabel("$\gamma$")
    ax.set_ylabel("$\sigma_{xy}$")

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

    ax.set_xlabel("$\gamma$")
    ax.set_ylabel("$\sigma_{xy}$")

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
    ax1.set_ylabel("$\sigma_{xy}$")

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
    ax1.set_ylabel("$\sigma_{xy}$")

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
    ax1.set_ylabel("$\sigma_{xy}$")

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
    strain = []
    for i in range(3):
        dir = f"/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/dGamma{i+1}/averaged_data"
        stress_ave.append(load_json(dir,"stress_ave"))
        stress_err.append(load_json(dir,"stress_err"))
        strain.append(load_json(dir,"strain_ave")[1:])

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

    shearrate = [0.01,0.001,0.0001]

    fig,ax = plt.subplots(figsize=(5,5))

    ax.plot(shearrate,eq_stress,linestyle="--")
    ax.errorbar(shearrate,eq_stress,yerr=eq_stress_err,marker='s',markersize=6,mfc='w',alpha = 0.9,linestyle='None',capsize=1)

    ax.set_xlabel("$\dot{\gamma}$")
    ax.set_ylabel("$\sigma_{xy}$")

    #ax.set_xscale('log')
    #ax.set_yscale('log')

    ax.set_xlim(0,0.011)
    ax.set_ylim(0,0.0125)

    ax.minorticks_on()
    ax.set_yticks(np.arange(0.001,0.0121,0.002))
    ax.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/stress_shearrate.pdf",dpi=300,bbox_inches='tight')

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
    ax1.set_ylabel("$\sigma_{xy}$")

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
    ax2.set_xlabel("$\gamma$")

    ax2.set_ylim((0,0.62))

    ax2.minorticks_on()
    ax2.set_xticks(np.arange(0,20.1,5))
    ax2.set_yticks(np.arange(0,0.61,0.1))
    ax2.tick_params(direction='in',which='minor',length=2,right=True,top=True)
    ax2.tick_params(direction='in',which='major',length=5,right=True,top=True)

    fig.savefig("/media/felipe/Files/Hydrogel_sim_experiments/FullExperiment1/Graficas/nematic_order_parameter.pdf",dpi=300,bbox_inches='tight')

    plt.show()



def main():

    #stress_strain()
    #stress_strain_loglog()
    #basic_analysis()
    #xldistance()
    #chain_analysis()
    #stress_shearrate()
    nematic_order_parameter()

if __name__=="__main__":
    main()