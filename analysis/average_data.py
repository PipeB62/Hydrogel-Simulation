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
    exp_num = len(data)
    iter = len(data[0])
    ave_data = []
    for t in range(iter):
        ave = 0
        for exp in range(exp_num):
            ave += data[exp][t]
        ave = ave/exp_num
        ave_data.append(ave)
    return ave_data

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

def average_basic_analysis():
    import sys 

    expdir = sys.argv[1]

    bonds_v_t = []
    dsq_v_t = []
    clusternum_v_t = []
    bigclustersz_v_t = []
    stress_v_t = []

    for exp in range(1,5):
        dir = f"{expdir}/exp{exp}/analysis_results"
        dir2 = f"{expdir}/exp{exp}"
        bonds_v_t.append(load_json(dir,"analysis_bonds_v_t"))
        dsq_v_t.append(load_json(dir,"analysis_Dsq_v_t"))
        clusternum_v_t.append(load_json(dir,"analysis_clusternum_v_t"))
        bigclustersz_v_t.append(load_json(dir,"analysis_bigclustersz_v_t"))
        stress_v_t.append([row[4] for row in load_data(dir2,f"presion_ave_shearing_2k_{exp}.data")])

    bonds_v_t_ave = average(bonds_v_t)
    dsq_v_t_ave = average(dsq_v_t)
    clusternum_v_t_ave = average(clusternum_v_t)
    bigclustersz_v_t_ave = average(bigclustersz_v_t)
    stress_v_t_ave = average(stress_v_t)

    savedir = f"{expdir}/averaged_data"
    write_json(savedir,"bonds_v_t_ave",bonds_v_t_ave)
    write_json(savedir,"dsq_v_t_ave",dsq_v_t_ave)
    write_json(savedir,"clusternum_v_t_ave",clusternum_v_t_ave)
    write_json(savedir,"bigclustersz_v_t_ave",bigclustersz_v_t_ave)
    write_json(savedir,"stress_v_t_ave",stress_v_t_ave)

def average_hole_analysis():
    import sys 
    import numpy as np

    expdir = sys.argv[1]

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
    binline = np.linspace(0,90,101)

    for exp in range(1,5):
        dir = f"{expdir}/exp{exp}/hole_analysis_results"

        av_hole_sz_pore.append(load_json(dir,"hole_analysis_pores_ave"))
        av_hole_sz_line.append(load_json(dir,"hole_analysis_line_ave"))

        hole_num_pore.append(load_json(dir,"hole_analysis_pores_hole_num"))
        hole_num_line.append(load_json(dir,"hole_analysis_line_hole_num"))

        current_distr_pore = load_json(dir,"hole_analysis_pores_distr")
        current_distr_line = load_json(dir,"hole_analysis_line_distr")

        histspore = []
        for cdistr in current_distr_pore:
            h,b = np.histogram(cdistr,bins=binpore)
            histspore.append(h)
        hole_distr_pore.append(histspore)

        histsline = []
        for cdistr in current_distr_line:
            h,b = np.histogram(cdistr,bins=binline)
            histsline.append(h)
        hole_distr_line.append(histsline)

        maxpore = [max(i) for i in current_distr_pore]
        maxline = [max(i) for i in current_distr_line]
        max_hole_sz_pore.append(maxpore)
        max_hole_sz_line.append(maxline)
    
    sample_distr_line,b = np.histogram(load_json(dir,"hole_analysis_line_sample_distr"),bins=binline)

    av_hole_sz_pore_ave = average(av_hole_sz_pore)
    av_hole_sz_line_ave = average(av_hole_sz_line)
    hole_num_pore_ave = average(hole_num_pore)
    hole_num_line_ave = average(hole_num_line)
    max_hole_sz_pore_ave = average(max_hole_sz_pore)
    max_hole_sz_line_ave = average(max_hole_sz_line)
    hole_distr_pore_ave = average_hist(hole_distr_pore)
    hole_distr_line_ave = average_hist(hole_distr_line)
   

    savedir = f"{expdir}/averaged_hole_data"
    write_json(savedir,"av_hole_sz_pore_ave",av_hole_sz_pore_ave)
    write_json(savedir,"av_hole_sz_line_ave",av_hole_sz_line_ave)
    write_json(savedir,"hole_num_pore_ave",hole_num_pore_ave)
    write_json(savedir,"hole_num_line_ave",hole_num_line_ave)
    write_json(savedir,"max_hole_sz_pore_ave",max_hole_sz_pore_ave)
    write_json(savedir,"max_hole_sz_line_ave",max_hole_sz_line_ave)
    write_json(savedir,"hole_distr_pore_ave",hole_distr_pore_ave)
    write_json(savedir,"hole_distr_line_ave",hole_distr_line_ave)
    write_json(savedir,"sample_distr_line",sample_distr_line)
    write_json(savedir,"binpore",binpore)
    write_json(savedir,"binline",binline)

def main():
    average_hole_analysis()

if __name__ == "__main__":
    main()
