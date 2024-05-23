def load_json(dir,name):
    import json
    with open (f'{dir}/hole_analysis_line_{name}.json') as f:
        content = json.load(f)
    data = json.loads(content)
    return(data)

def graficas_analysis_base(dirs):
    import numpy as np
    import matplotlib.pyplot as plt 
    import json

    colors = ['red','green','blue']
    markers = ['|','_','.']

    fig1,ax1 = plt.subplots()
    ax1.set_title('Bonds')
    ax1.legend(dirs)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Number of bonds')

    fig2,ax2 = plt.subplots()
    ax2.set_title('Number of clusters')
    ax2.legend(dirs)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Number of clusters')

    fig3,ax3 = plt.subplots()
    ax3.set_title('Biggest cluster size')
    ax3.legend(dirs)
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Number of atoms in biggest cluster')

    fig4,ax4 = plt.subplots()
    ax4.set_title('Non affine sq displacement')
    ax4.legend(dirs)
    ax4.set_xlabel('Time')
    ax4.set_ylabel('$D^2$')  

    fig5,ax5 = plt.subplots(3,1,sharex=True)

    
    for k,file in enumerate(dirs):

        #leer archivos json
        with open (f'{file}/analysis_bonds_v_t.json') as f:
            bonds_v_t = json.load(f)
        with open (f'{file}/analysis_clusternum_v_t.json') as f:
            clusternum_v_t = json.load(f)
        with open (f'{file}/analysis_bigclustersz_v_t.json') as f:
            bigclustersz_v_t = json.load(f)
        with open (f'{file}/analysis_x1.json') as f:
            x1 = json.load(f)
        with open (f'{file}/analysis_Dsq_v_t.json') as f:
            D_sq_v_t = json.load(f)
        with open (f'{file}/analysis_x2.json') as f:
            x2 = json.load(f)
        
        ax1.plot(x1,bonds_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')
        ax2.plot(x1,clusternum_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')
        ax3.plot(x1,bigclustersz_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')
        ax4.plot(x2,D_sq_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')
        ax5[0].plot(x1,bonds_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')
        ax5[1].plot(x1,clusternum_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')
        ax5[2].plot(x1,D_sq_v_t,c=colors[k],marker=markers[k],fillstyle='none',linestyle='None')

    plt.show()

def graficas_formation_analysis(dir):
    pass

def graficas_hole_analysis_pore(dir):
    import json
    import numpy as np
    import matplotlib.pyplot as plt

    #leer archivos json
    with open (f'{dir}/hole_analysis_pores_distr.json') as f:
        distr = json.load(f)
    with open (f'{dir}/hole_analysis_pores_ave.json') as f:
        ave = json.load(f)
    with open (f'{dir}/hole_analysis_pores_hole_num.json') as f:
        hole_num = json.load(f)
    with open (f'{dir}/hole_analysis_pores_calc_frames.json') as f:
        calc_frames = json.load(f)

    big_r = [max(k) for k in distr]

    bin = np.linspace(0,7,50)

    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    fig4,ax4 = plt.subplots()

    ax1.hist(distr[0],bins=bin,histtype='step',label='Primer frame')
    ax1.hist(distr[-1],bins=bin,histtype='step',label='Ultimo frame')
    ax1.legend()
    ax1.set_title('Distribuciones de huecos')
    ax1.set_xlabel('Tamaño de hueco')
    ax1.set_ylabel('Numero de huecos')

    ax2.scatter(calc_frames,ave)
    ax2.set_title("Radio promedio de hueco")
    ax2.set_xlabel('Frame')
    ax2.set_ylabel('Radio')

    ax3.scatter(calc_frames,hole_num)
    ax3.set_title("Porcentaje de huecos")
    ax3.set_xlabel('Frame')
    ax3.set_ylabel('Porcentaje de huecos')

    ax4.scatter(calc_frames,big_r)
    ax4.set_title("Radio del hueco mas grande")
    ax4.set_xlabel('Frame')
    ax4.set_ylabel('RadioS')
    
    plt.show()

def graficas_hole_analysis_line(dir):
    import json
    import numpy as np
    import matplotlib.pyplot as plt

    print('Cargando datos ...')
    distr = load_json(dir,"distr")
    ave = load_json(dir,"ave")
    calc_frames = load_json(dir,"calc_frames")
    hole_num = load_json(dir,"hole_num")

    """
    with open (f'{dir}/hole_analysis_line_sample_distr.json') as f:
        sample_distr = json.load(f)
    with open (f'{dir}/hole_analysis_line_prob_distr.json') as f:
        prob_distr = json.load(f)
    with open (f'{dir}/hole_analysis_line_bin.json') as f:
        bin1 = json.load(f)
    """
    print('Generando graficas')
    big_r = [max(k) for k in distr]

    bin = np.linspace(0,90,90)

    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    fig4,ax4 = plt.subplots()
    fig5,ax5 = plt.subplots()
    fig6,ax6 = plt.subplots()

    ax1.hist(distr[0],bins=bin,histtype='step',label='Primer frame')
    ax1.hist(distr[-1],bins=bin,histtype='step',label='Ultimo frame')
    ax1.legend()
    ax1.set_title('Distribuciones de huecos')
    ax1.set_xlabel('Longitud de hueco')
    ax1.set_ylabel('Numero de huecos')
    
    ax2.scatter(calc_frames,ave)
    ax2.set_title('Longitud promedio de hueco')
    ax2.set_xlabel('Frame')
    ax2.set_ylabel('Radio')

    ax3.scatter(calc_frames,hole_num)
    ax3.set_title("Porcentaje de huecos")
    ax3.set_xlabel('Frame')
    ax3.set_ylabel('Porcentaje de huecos')

    """
    ax4.scatter(calc_frames,big_r)
    ax4.set_title("Longitud del hueco mas grande")
    ax4.set_xlabel('Frame')
    ax4.set_ylabel('RadioS')

    ax5.hist(sample_distr[0],bins=bin,histtype='step',label='Primer frame')
    ax5.hist(sample_distr[-1],bins=bin,histtype='step',label='Ultimo frame')
    ax5.legend()
    ax5.set_title('Distribucion de muestreo')
    ax5.set_xlabel('Longitud de linea')
    ax5.set_ylabel('Numero de lineas')

    ax6.stairs(prob_distr[0],bin1,label='Primer frame')
    ax6.stairs(prob_distr[-1],bin1,label='Ultimo frame')
    ax6.legend()
    ax6.set_title('Distribucion de probabilidad de huecos')
    ax6.set_xlabel('Longitud de linea')
    ax6.set_ylabel('Probabilidad')
    """
    plt.show()

def main():
    import sys 
    dirs = sys.argv[1:]

    #graficas_analysis_base(dirs)

    #graficas_hole_analysis_pore(dirs[0])

    graficas_hole_analysis_line(dirs[0]) 



if __name__ == '__main__':
    main()