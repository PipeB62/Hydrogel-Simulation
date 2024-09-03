from random import randint
import os 
import sys


#Pedir al usuario parametros de simulacion
L,shear_rate,max_strain,savedir,init_system,inputfilesdir,last_name = sys.argv[1:]
L = float(L)
shear_rate = float(shear_rate)

#Generar semillas
seeds = []
for i in range(10):
    seeds.append(randint(100000, 999999))

#Parametros fijos
save_every = 100000
timestep = 0.001
damp = 0.1 
temp = 0.05 #Temperatura objetivo

#max_strain = 20
max_strain=int(max_strain)
delta_gamma = 0.01

flip_strain = 0.5*L
delta_gamma_distance = delta_gamma*L
n = int(max_strain/delta_gamma)
relaxation_time = delta_gamma/shear_rate
shear_every = int(relaxation_time/timestep)
shear_iter_num = n*shear_every 
n_pres_av = int(shear_every/2)

#Definir unidades, condiciones de frontera, tipos de atomos. Newton on necesario para potencial de tres cuerpos. No guardar log en archivo
inicializacion = (f"log log_shearing_{last_name}.lammps \n\n" 
                  "units lj \n" 
                  "boundary p p p \n" 
                  "atom_style full \n" 
                  "newton on \n\n")

#Definir interacciones. centro-centro = lj. xl-xl = patches.sw. threebody para xl-xl-xl
pair_definitions = ("pair_style hybrid/overlay lj/cut 1.122462 sw threebody off threebody/table\n"
                    f"pair_coeff * * threebody/table {inputfilesdir}/threebody.3b NULL Mon NULL Xl\n"
                    f"pair_coeff 2 4 sw {inputfilesdir}/patches.sw NULL Mon NULL Xl\n"
                    f"pair_coeff 2 2 sw {inputfilesdir}/patches.sw NULL Mon NULL Xl\n"
                    "pair_coeff 4 4 none\n"
                    "pair_coeff 1 1 lj/cut 1.0 1.0\n"
                    "pair_coeff 1 2 none\n"
                    "pair_coeff 1 3 lj/cut 1.0 1.0\n"
                    "pair_coeff 2 3 none\n"
                    "pair_coeff 3 3 lj/cut 1.0 1.0\n"
                    "pair_coeff 3 4 none\n"
                    "pair_coeff 1 4 none\n\n")

#Definir bonds y angles armonicos para uniones centro-patch
bonds_angles = ("bond_style harmonic \n"
                "angle_style harmonic \n\n")

#importar y cambiar caja
#read_system = (f"read_data {directory}/system_formation_{last_name}.data extra/special/per/atom 10 \n" 
#               "change_box all triclinic \n\n") 

read_system = (f"read_data {init_system} extra/special/per/atom 10 \n" 
               "change_box all triclinic \n\n") 

#Computes dentro de lammps
pressure = ("compute presion all pressure NULL pair bond \n"
           f"fix presion_ave all ave/time 1 {n_pres_av} {shear_every} c_presion[*] file {savedir}/presion_ave_shearing_{last_name}.data \n\n")

#Visualizacion
visualization = (f"dump mydmp all atom {shear_every} {savedir}/dump_shearing_{last_name}.lammpstrj \n"
                  "dump_modify mydmp scale no \n\n"
                 f"thermo {save_every} \n"
                  "thermo_style custom step temp xy c_presion[4] \n\n")

#Definir dinamica nve
fix_nve = (f"timestep {timestep} \n"
           "fix mynve all nve \n\n")

#Definir termostato de langevin. 
langevin = (f"fix thermostat all langevin {temp} {temp} {damp} {seeds[2]} \n\n")

shear = (f"variable n_loop loop {n} \n"
         "variable xy equal xy \n" 
         f"variable flip_strain equal {flip_strain} \n"
         "label runloop \n"
         f"change_box all xy delta {delta_gamma_distance} remap \n"
         'if "${xy} == ${flip_strain}" then "change_box all xy final -${flip_strain}" \n'
         f"run {shear_every} \n"
         "next n_loop \n"
         "jump SELF runloop \n\n")

write_data = f"write_data {savedir}/system_shearing_{last_name}.data"

input_file = [inicializacion,bonds_angles,read_system,pair_definitions,pressure,visualization,fix_nve,langevin,
              shear,write_data]
with open(f"{savedir}/input_shearing_{last_name}.lammps","a") as f:  
    f.truncate(0)
    f.writelines(input_file)


