from random import randint

#Pedir al usuario parametros de simulacion
L = input("L = ")
xl_num = input("Number of cross-linkers: ")
mon_num = input("Number of monomers: ")
directory = input("File directory: ")
last_name = input("File name locator: ")

#Generar semillas
seeds = []
for i in range(10):
    seeds.append(randint(100000, 999999))

#Parametros fijos
save_every = 10000
timestep = 0.002
temp = 0.05
damp = 0.5
iter_num_1_3 = 500000
itern_num_2 = 100000

#Definir unidades, condiciones de frontera, tipos de atomos. Newton on necesario para potencial de tres cuerpos. No guardar log en archivo
inicializacion = ("log none \n\n" 
                  "units lj \n" 
                  "boundary p p p \n" 
                  "atom_style full \n" 
                  "newton on \n\n")

#Definir caja 
sim_box = (f"region sim_box block -{L} {L} -{L} {L} -{L} {L} \n"
           "create_box 4 sim_box &\n"
           "   bond/types 1 &\n"
           "   angle/types 2 &\n"
           "   extra/bond/per/atom 10 &\n"
           "   extra/angle/per/atom 6 &\n"
           "   extra/special/per/atom 10\n\n")

#Crear atomos y definir masas
atom_types = ("mass 1 1.0 \n" #Monomer
              "mass 2 0.01 \n" #Patches
              "mass 3 1.0 \n" #Xlinker
              "mass 4 0.01 \n\n") #Xlinker patches

#Definir interacciones. centro-centro = lj. xl-xl = patches.sw. threebody para xl-xl-xl
pair_definitions = ("pair_style hybrid/overlay lj/cut 1.122462 sw threebody off threebody/table\n"
                    "pair_coeff * * threebody/table ../input_data/threebody.3b NULL Mon NULL Xl\n"
                    "pair_coeff 2 4 sw ../input_data/patches.sw NULL Mon NULL Xl\n"
                    "pair_coeff 2 2 sw ../input_data/patches.sw NULL Mon NULL Xl\n"
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
                "bond_coeff 1 100.0 0.361231024154685 \n\n"
                "angle_style harmonic \n"
                "angle_coeff 1 100.0 180 \n"
                "angle_coeff 2 100.0 109.4712 \n\n")

#Generar moleculas
spawn = (f"region spawn_box block -{L} {L} -{L} {L} -{L} {L} \n\n" 
         "molecule 1 ../input_data/monomer.mol \n" 
         "molecule 2 ../input_data/xlinker.mol \n\n"
         f"create_atoms 0 random {xl_num} {seeds[0]} spawn_box mol 2 {seeds[1]} overlap 1.13 maxtry 5000 \n"
         f"create_atoms 0 random {mon_num} {seeds[2]} spawn_box mol 1 {seeds[3]} overlap 1.13 maxtry 5000 \n\n")

#Computes dentro de lammps
computes = ("compute pot_ene all pe \n\n")

#Visualizacion
visualization = (f"dump mydmp all atom {save_every} {directory}/dump_formation_{last_name}.lammpstrj \n"
                 "dump_modify mydmp scale no \n"
                 f"thermo {save_every} \n\n")

#Minimizar energia para evitar particulas superpuestas
minimize_energy = "minimize 1.0e-4 1.0e-6 1000 100000 \n\n"

#Guardar energia potencial en archivo
fix_save_poten = f"fix save_pot_ene all ave/time 1 1 {save_every} c_pot_ene file {directory}/pot_ene_formation_{last_name}.data \n\n" 

#Definir dinamica NVE 
fix_nve = (f"timestep {timestep} \n"
           "fix mynve all nve \n\n" )

#Percolation_fun de python
percolation_check = ("variable check_percolation python percolation \n"
                     "python percolation return v_check_percolation format i file percolation_fun.py \n\n")

#fase 1: subir temperatura
fase1 = (f"fix thermostat all langevin 0.0 {temp} {damp} {seeds[4]} \n" 
         f"run {iter_num_1_3} start 0 \n\n")

#fase 2: formacion
fase2 = (f"fix thermostat all langevin {temp} {temp} {damp} {seeds[5]} \n"
         "label formation_loop \n"
         f"run {itern_num_2} start 0 \n"
         "write_data system_check.data \n"
         "print ${check_percolation} \n"
         'if "${check_percolation} == 1" then "jump SELF bajar_temp" else "jump SELF formation_loop" \n\n')

#fase 3: bajar temperatura
fase3 = ("label bajar_temp \n"
         f"fix thermostat all langevin {temp} 0.0 {damp} {seeds[6]} \n"
         f"run {iter_num_1_3} start 0 \n\n")

#Guardar sistema final
write_data = f"write_data {directory}/system_formation_{last_name}.data"


#Crear y escribir. Reescribir si el archivo ya existe
input_file = [inicializacion,sim_box,atom_types,pair_definitions,bonds_angles,spawn,computes,visualization,
              minimize_energy,fix_save_poten,fix_nve,percolation_check,fase1,fase2,fase3,write_data]
with open(f"{directory}/input_formation_{last_name}.lammps","a") as f:  
    f.truncate(0)
    f.writelines(input_file)



