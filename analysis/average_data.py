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

    expdir = sys.argv[1] #directorio que contiene todos los experimentos
    savedir = sys.argv[2] #directorio donde guardar datos promediados
    dataname = sys.argv[3] #nombre del archivo json que contiene los datos en cada experimento. (sin el .json)
    expnum = int(sys.argv[4]) #numero de experimentos

    data = []

    for exp in range(1,expnum+1):
        dir = f"{expdir}/exp{exp}/analysis_results"
        data.append(load_json(dir,dataname))

    data_ave, data_std = average(data)

    write_json(savedir,f"{dataname}_ave",data_ave)
    write_json(savedir,f"{dataname}_err",data_std)
    
def main():
    average_general()

if __name__ == "__main__":
    main()
