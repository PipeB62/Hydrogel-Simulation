import numpy as np
import os 
import sys

#Dar como input lambda, sigma y el directorio del experimento en ese orden

def V_3_f(r,sigma,r_min,r_c):
    if r < r_min:
        V_bond = 1
    else:
        V_bond = -2*((sigma**4)/(2*r**4) -1)*np.exp(sigma/(r-r_c) +2)
    return V_bond

def drV_3_f(r,sigma,r_min,r_c):
    if r < r_min:
        drV_3 = 0
    else:
        drV_3 = 2*((2*sigma**4)/r**5 + ((sigma**4)/(2*r**4) -1)*((sigma)/(r-r_c)**2))*np.exp(sigma/(r-r_c) +2)
    return drV_3

'''
#Obtener direccion de carpeta input_data
dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path_split = dir_path.split("/")
dir_path_split[-1]="input_data"
inputfilesdir = "/".join(dir_path_split)
'''

N = 12
a = 1.5
sigma = float(sys.argv[2])
r_c = sigma*a
r_min = sigma
lam = float(sys.argv[1])
eps = 1

inputfilesdir = f"{sys.argv[3]}/input_data"

minR = 0.01
maxR = r_c-0.001

r = np.linspace(minR,maxR,N)
angles = np.linspace(180/(4*N),180*(1-1/(4*N)),2*N)

str_3btable = f"3b_{lam}.table"
str_3b2table = f"3b_2_{lam}.table"

# ---Tabla 1---
data = np.zeros((N*N*(N+1),11))
data[:,0] = np.arange(1,N*N*(N+1)+1)

index = 0
for i1 in range(N):
    r_ij = r[i1]
    #print('rij '+str(r_ij))
    for i2 in range(i1,N):
        r_ik = r[i2]
        #print('rik '+str(r_ik))
        for i3 in range(2*N):
            theta = angles[i3]

            data[index,1:4] = [r_ij,r_ik,theta]
            data[index,4] = -lam*eps*V_3_f(r_ik,sigma,r_min,r_c)*drV_3_f(r_ij,sigma,r_min,r_c)
            data[index,5] = -lam*eps*V_3_f(r_ij,sigma,r_min,r_c)*drV_3_f(r_ik,sigma,r_min,r_c)
            data[index,6] = lam*eps*V_3_f(r_ik,sigma,r_min,r_c)*drV_3_f(r_ij,sigma,r_min,r_c)
            data[index,7] = 0
            data[index,8] = lam*eps*V_3_f(r_ij,sigma,r_min,r_c)*drV_3_f(r_ik,sigma,r_min,r_c)
            data[index,9] = 0
            data[index,10] = lam*eps*V_3_f(r_ij,sigma,r_min,r_c)*V_3_f(r_ik,sigma,r_min,r_c)
            index += 1
            #print(index)

file = open(f"{inputfilesdir}/{str_3btable}",'a') #abrir archivo en modo append
file.truncate(0) #borrar contenido actual
file.write('# Tabulated threebody potential\n\n'+
           'ENTRY1\n'+
           'N '+str(N)+' rmin '+str(minR)+' rmax '+str(maxR)+'\n\n')
np.savetxt(file,data,delimiter=' ',newline='\n',comments='',fmt=' '.join(['%i'] + ['%1.4f']*10))
file.close()

# ---Tabla 2---
data2 = np.zeros((2*N*N*N,11))
data2[:,0] = np.arange(1,2*N*N*N+1)

index = 0
for i1 in range(N):
    r_ij = r[i1]
    #print('rij '+str(r_ij))
    for i2 in range(N):
        r_ik = r[i2]
        #print('rik '+str(r_ik))
        for i3 in range(2*N):
            theta = angles[i3]

            data2[index,1:4] = [r_ij,r_ik,theta]
            data2[index,4] = -lam*eps*V_3_f(r_ik,sigma,r_min,r_c)*drV_3_f(r_ij,sigma,r_min,r_c)
            data2[index,5] = -lam*eps*V_3_f(r_ij,sigma,r_min,r_c)*drV_3_f(r_ik,sigma,r_min,r_c)
            data2[index,6] = lam*eps*V_3_f(r_ik,sigma,r_min,r_c)*drV_3_f(r_ij,sigma,r_min,r_c)
            data2[index,7] = 0
            data2[index,8] = lam*eps*V_3_f(r_ij,sigma,r_min,r_c)*drV_3_f(r_ik,sigma,r_min,r_c)
            data2[index,9] = 0
            data2[index,10] = lam*eps*V_3_f(r_ij,sigma,r_min,r_c)*V_3_f(r_ik,sigma,r_min,r_c)
            index += 1
            #print(index)

file = open(f"{inputfilesdir}/{str_3b2table}",'a') #abrir archivo en modo append
file.truncate(0) #borrar contenido actual
file.write('# Tabulated threebody potential\n\n'+
           'ENTRY1\n'+
           'N '+str(N)+' rmin '+str(minR)+' rmax '+str(maxR)+'\n\n')
np.savetxt(file,data2,delimiter=' ',newline='\n',comments='',fmt=' '.join(['%i'] + ['%1.4f']*10))
file.close()

#tablepath = "/home/felipe/HYDROGELS/Hydrogel-Simulation/input_data/"
tablepath = inputfilesdir

#3b file
file = open(f"{inputfilesdir}/threebody_{lam}.3b",'a')
file.truncate(0) #borrar contenido actual
for i in range(8):
    if i == 0:
        s1 = 'Mon\n'
        s2 = 'Mon\n'
        s3 = 'Mon\n'
        table = tablepath+'/'+str_3btable+'\n'
    elif i == 1:
        s1 = 'Xl\n'
        s2 = 'Mon\n'
        s3 = 'Mon\n'
        table = tablepath+'/'+str_3btable+'\n'
    elif i == 2: 
        s1 = 'Mon\n'
        s2 = 'Xl\n'
        s3 = 'Mon\n'
        table = tablepath+'/'+str_3b2table+'\n'
    elif i == 3:
        s1 = 'Mon\n'
        s2 = 'Mon\n'
        s3 = 'Xl\n'
        table = tablepath+'/'+str_3b2table+'\n'
    elif i == 4:
        s1 = 'Mon\n'
        s2 = 'Xl\n'
        s3 = 'Xl\n'
        table = tablepath+'/'+str_3btable+'\n'
    elif i == 5: 
        s1 = 'Xl\n'
        s2 = 'Mon\n'
        s3 = 'Xl\n'
        table = tablepath+'/'+str_3b2table+'\n'
    elif i == 6:
        s1 = 'Xl\n'
        s2 = 'Xl\n'
        s3 = 'Mon\n'
        table = tablepath+'/'+str_3b2table+'\n'
    elif i == 7:
        s1 = 'Xl\n'
        s2 = 'Xl\n'
        s3 = 'Xl\n'
        table = tablepath+'/'+str_3btable+'\n'
    file.write(s1+
            s2+
            s3+
            str(r_c-0.001)+'\n'+
            table+
            'ENTRY1\n'+
            'linear\n'+
            str(N)+'\n') 
file.close()




