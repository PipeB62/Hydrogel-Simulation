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

def main():
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

if __name__ == "__main__":
    main()
