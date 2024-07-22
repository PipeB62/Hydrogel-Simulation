
#Dar como input lambda, sigma y el directorio del experimento en ese orden

def V_3_f(r,sigma,r_min,r_c):
    import numpy as np

    if r < r_min:
        V_3 = 1
    elif r >= r_min and r <= r_c:
        V_3 = -2*((sigma**4)/(2*r**4) -1)*np.exp(sigma/(r-r_c) +2)
    return V_3

def V_tb_f(r_ij,r_ik,sigma,r_min,r_c,lam,eps):
    import numpy as np

    V_threebody = lam*eps*V_3_f(r_ij,sigma,r_min,r_c)*V_3_f(r_ik,sigma,r_min,r_c)

    return V_threebody

def main():

    import numpy as np
    import os 
    import sys
    
    N = 12
    a = 1.5
    sigma = float(sys.argv[2])
    r_c = sigma*a
    r_min = sigma
    lam = float(sys.argv[1])
    eps = 1

    h = 0.0001

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

                dVtb_drij = (V_tb_f(r_ij+(h/2),r_ik,sigma,r_min,r_c,lam,eps) - V_tb_f(r_ij-(h/2),r_ik,sigma,r_min,r_c,lam,eps))/h
                dVtb_drik = (V_tb_f(r_ij,r_ik+(h/2),sigma,r_min,r_c,lam,eps) - V_tb_f(r_ij,r_ik-(h/2),sigma,r_min,r_c,lam,eps))/h

                data[index,1:4] = [r_ij,r_ik,theta]
                data[index,4] = -(1/r_ij)*(dVtb_drij) #fi1
                data[index,5] = -(1/r_ik)*(dVtb_drik) #fi2
                data[index,6] = (1/r_ij)*(dVtb_drij) #fj1
                data[index,7] = 0  #fj2
                data[index,8] = (1/r_ik)*(dVtb_drik) #fk1
                data[index,9] = 0 #fk2
                data[index,10] = V_tb_f(r_ij,r_ik,sigma,r_min,r_c,lam,eps) #e
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

                dVtb_drij = (V_tb_f(r_ij+(h/2),r_ik,sigma,r_min,r_c,lam,eps) - V_tb_f(r_ij-(h/2),r_ik,sigma,r_min,r_c,lam,eps))/h
                dVtb_drik = (V_tb_f(r_ij,r_ik+(h/2),sigma,r_min,r_c,lam,eps) - V_tb_f(r_ij,r_ik-(h/2),sigma,r_min,r_c,lam,eps))/h

                data2[index,1:4] = [r_ij,r_ik,theta]
                data2[index,4] = -(1/r_ij)*(dVtb_drij) #fi1
                data2[index,5] = -(1/r_ik)*(dVtb_drik) #fi2
                data2[index,6] = (1/r_ij)*(dVtb_drij) #fj1
                data2[index,7] = 0  #fj2
                data2[index,8] = (1/r_ik)*(dVtb_drik) #fk1
                data2[index,9] = 0 #fk2
                data2[index,10] = V_tb_f(r_ij,r_ik,sigma,r_min,r_c,lam,eps) #e
                index += 1
                

    file = open(f"{inputfilesdir}/{str_3b2table}",'a') #abrir archivo en modo append
    file.truncate(0) #borrar contenido actual
    file.write('# Tabulated threebody potential\n\n'+
            'ENTRY1\n'+
            'N '+str(N)+' rmin '+str(minR)+' rmax '+str(maxR)+'\n\n')
    np.savetxt(file,data2,delimiter=' ',newline='\n',comments='',fmt=' '.join(['%i'] + ['%1.4f']*10))
    file.close()

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


if __name__ == "__main__":
    main()

