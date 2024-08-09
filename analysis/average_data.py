def load_json(dir,name):
    import json
    with open (f'{dir}/{name}.json') as f:
        data = json.load(f)
    if isinstance(data,str):
        data = json.loads(data)
    return(data)

def write_json(savedir,name,data):
    import json 
    with open (f'{savedir}/{name}.json','w') as f:
        json.dump(data,f)

def load_data(dir,name):
    import numpy as np
    a = np.loadtxt(f"{dir}/{name}")
    return a.tolist()

def average(data):
    import numpy as np

    exp_num = len(data)
    iter = len(data[0])
    ave_data = []
    std = []
    for t in range(iter):
        ave = 0
        for exp in range(exp_num):
            ave += data[exp][t]
        ave = ave/exp_num
        ave_data.append(ave)

        var = 0
        for exp in range(exp_num):
            var += (data[exp][t]-ave)**2
        var = var/(exp_num-1)
        std.append(np.sqrt(var))

    return ave_data, std

def average_hist(data):
    exp_num = len(data)
    iter_num = len(data[0])
    col_num = len(data[0][0])
    ave_data=[]
    for iter in range(iter_num):
        h = []
        for col in range(col_num):
            ave = 0
            for exp in range(exp_num):
                ave += data[exp][iter][col]
            ave = ave/exp_num
            h.append(ave)
        ave_data.append(h)
    return ave_data

def average_general():
    import sys 

    expdir = sys.argv[1]
    savedir = sys.argv[2]
    dataname = sys.argv[3]
    expnum = int(sys.argv[4])

    data = []

    for exp in range(1,expnum+1):
        dir = f"{expdir}/exp{exp}/analysis_results"
        data.append(load_json(dir,dataname))

    data_ave, data_std = average(data)

    write_json(savedir,f"{dataname}_ave",data_ave)
    write_json(savedir,f"{dataname}_std",data_std)
    
def average_basic_analysis():
    import sys 

    expdir = sys.argv[1]

    bonds_v_t = []
    dsq_v_t = []
    clusternum_v_t = []
    bigclustersz_v_t = []
    stress_v_t = []

    for exp in range(1,5+1):
        dir = f"{expdir}/exp{exp}/analysis_results"
        dir2 = f"{expdir}/exp{exp}"
        bonds_v_t.append(load_json(dir,"analysis_bonds_v_t"))
        dsq_v_t.append(load_json(dir,"analysis_Dsq_v_t"))
        clusternum_v_t.append(load_json(dir,"analysis_clusternum_v_t"))
        bigclustersz_v_t.append(load_json(dir,"analysis_bigclustersz_v_t"))
        stress_v_t.append([row[4] for row in load_data(dir2,f"presion_ave_shearing_16k_{exp}.data")])

    bonds_v_t_ave, bonds_v_t_std = average(bonds_v_t)
    dsq_v_t_ave, dsq_v_t_std = average(dsq_v_t)
    clusternum_v_t_ave, clusternum_v_t_std = average(clusternum_v_t)
    bigclustersz_v_t_ave, bigclustersz_v_t_std = average(bigclustersz_v_t)
    stress_v_t_ave, stress_v_t_std = average(stress_v_t)

    savedir = f"{expdir}/averaged_data"

    write_json(savedir,"bonds_v_t_ave",bonds_v_t_ave)
    write_json(savedir,"dsq_v_t_ave",dsq_v_t_ave)
    write_json(savedir,"clusternum_v_t_ave",clusternum_v_t_ave)
    write_json(savedir,"bigclustersz_v_t_ave",bigclustersz_v_t_ave)
    write_json(savedir,"stress_v_t_ave",stress_v_t_ave)

    write_json(savedir,"bonds_v_t_std",bonds_v_t_std)
    write_json(savedir,"dsq_v_t_std",dsq_v_t_std)
    write_json(savedir,"clusternum_v_t_std",clusternum_v_t_std)
    write_json(savedir,"bigclustersz_v_t_std",bigclustersz_v_t_std)
    write_json(savedir,"stress_v_t_std",stress_v_t_std)

def average_hole_analysis():
    import sys 
    import numpy as np

    expdir = sys.argv[1]
    expnum = 5

    av_hole_sz_pore = []
    max_hole_sz_pore = []
    hole_distr_pore = []
    hole_num_pore = []
    binpore = np.linspace(0,10,51)

    av_hole_sz_line = []
    max_hole_sz_line = []
    hole_distr_line = []
    sample_distr_line = []
    hole_num_line = []
    sample_distr = []
    binline = np.linspace(0,90,101)

    for exp in range(1,expnum+1):
        dir = f"{expdir}/exp{exp}/hole_analysis_results"

        av_hole_sz_pore.append(load_json(dir,"hole_analysis_pores_ave"))
        av_hole_sz_line.append(load_json(dir,"hole_analysis_line_ave"))

        hole_num_pore.append(load_json(dir,"hole_analysis_pores_hole_num"))
        hole_num_line.append(load_json(dir,"hole_analysis_line_hole_num"))

        current_distr_pore = load_json(dir,"hole_analysis_pores_distr")
        current_distr_line = load_json(dir,"hole_analysis_line_distr")
        current_sample_distr = load_json(dir,"hole_analysis_line_sample_distr")

        histspore = []
        for cdistr in current_distr_pore:
            h,b = np.histogram(cdistr,bins=binpore)
            histspore.append(h.tolist())
        hole_distr_pore.append(histspore)

        histsline = []
        for cdistr in current_distr_line:
            h,b = np.histogram(cdistr,bins=binline)
            histsline.append(h.tolist())
        hole_distr_line.append(histsline)

        sample_hists = []
        for cdistr in current_sample_distr:
            h,b = np.histogram(cdistr,bins=binline)
            sample_hists.append(h.tolist())
        sample_distr.append(sample_hists)

        maxpore = [max(i) for i in current_distr_pore]
        maxline = [max(i) for i in current_distr_line]
        max_hole_sz_pore.append(maxpore)
        max_hole_sz_line.append(maxline)
    

    av_hole_sz_pore_ave, av_hole_sz_pore_std = average(av_hole_sz_pore)
    av_hole_sz_line_ave, av_hole_sz_line_std = average(av_hole_sz_line)
    hole_num_pore_ave, hole_num_pore_std = average(hole_num_pore)
    hole_num_line_ave, hole_num_line_std = average(hole_num_line)
    max_hole_sz_pore_ave, max_hole_sz_pore_std = average(max_hole_sz_pore)
    max_hole_sz_line_ave, max_hole_sz_line_std = average(max_hole_sz_line)
    hole_distr_pore_ave = average_hist(hole_distr_pore)
    hole_distr_line_ave = average_hist(hole_distr_line)
    sample_distr_line_ave = average_hist(sample_distr)
   

    savedir = f"{expdir}/averaged_hole_data"
    write_json(savedir,"av_hole_sz_pore_ave",av_hole_sz_pore_ave)
    write_json(savedir,"av_hole_sz_line_ave",av_hole_sz_line_ave)
    write_json(savedir,"hole_num_pore_ave",hole_num_pore_ave)
    write_json(savedir,"hole_num_line_ave",hole_num_line_ave)
    write_json(savedir,"max_hole_sz_pore_ave",max_hole_sz_pore_ave)
    write_json(savedir,"max_hole_sz_line_ave",max_hole_sz_line_ave)
    write_json(savedir,"hole_distr_pore_ave",hole_distr_pore_ave)
    write_json(savedir,"hole_distr_line_ave",hole_distr_line_ave)
    write_json(savedir,"sample_distr_line_ave",sample_distr_line_ave)
    write_json(savedir,"binpore",binpore.tolist())
    write_json(savedir,"binline",binline.tolist())

    write_json(savedir,"av_hole_sz_pore_std",av_hole_sz_pore_std)
    write_json(savedir,"av_hole_sz_line_std",av_hole_sz_line_std)
    write_json(savedir,"hole_num_pore_std",hole_num_pore_std)
    write_json(savedir,"hole_num_line_std",hole_num_line_std)
    write_json(savedir,"max_hole_sz_pore_std",max_hole_sz_pore_std)
    write_json(savedir,"max_hole_sz_line_std",max_hole_sz_line_std)

def average_demixing_param():

    import sys 
    import numpy as np

    expdir = sys.argv[1]
    expnum = 5

    phi2 = []
    phi4 = []
    phi5 = []
    phi7 = []
    phi10 = []

    for exp in range(1,expnum+1):
        dir = f"{expdir}/exp{exp}/analysis_results"
        phi2.append(load_json(dir,"demixing_param_2"))
        phi4.append(load_json(dir,"demixing_param_4"))
        phi5.append(load_json(dir,"demixing_param_5"))
        phi7.append(load_json(dir,"demixing_param_7"))
        phi10.append(load_json(dir,"demixing_param_10"))

    phi2_ave, std_phi2 = average(phi2)
    phi4_ave, std_phi4 = average(phi4)
    phi5_ave, std_phi5 = average(phi5)
    phi7_ave, std_phi7 = average(phi7)
    phi10_ave, std_phi10 = average(phi10)

    savedir = f"{expdir}/averaged_data"

    write_json(savedir,"demixing_param_2_ave",phi2_ave)
    write_json(savedir,"demixing_param_4_ave",phi4_ave)
    write_json(savedir,"demixing_param_5_ave",phi5_ave)
    write_json(savedir,"demixing_param_7_ave",phi7_ave)
    write_json(savedir,"demixing_param_10_ave",phi10_ave)

    write_json(savedir,"demixing_param_2_std",std_phi2)
    write_json(savedir,"demixing_param_4_std",std_phi4)
    write_json(savedir,"demixing_param_5_std",std_phi5)
    write_json(savedir,"demixing_param_7_std",std_phi7)
    write_json(savedir,"demixing_param_10_std",std_phi10)

def average_fractal_dimension():

    import sys 
    import numpy as np

    expdir = sys.argv[1]
    expnum = 5

    D = []

    for exp in range(1,expnum+1):
        dir = f"{expdir}/exp{exp}/analysis_results"
        D.append(load_json(dir,"fractal_dimension"))

    D_ave, D_std = average(D)

    savedir = f"{expdir}/averaged_data"

    write_json(savedir,"D_ave",D_ave)

    write_json(savedir,"D_std",D_std)


def main():
    average_general()

if __name__ == "__main__":
    main()
