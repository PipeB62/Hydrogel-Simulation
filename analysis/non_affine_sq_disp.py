import os
os.environ['OVITO_GUI_MODE'] = '1'

#Importar librerias
from ovito.io import import_file, export_file
from ovito.modifiers import SelectTypeModifier,AtomicStrainModifier, DeleteSelectedModifier
from ovito.pipeline import ReferenceConfigurationModifier
import matplotlib.pyplot as plt
import numpy as np

#Version
version = input("Version: ")

#Importar dumps de shearing
input_filefolder = "/home/pipe/lammps/code/HYDROGELS/selfassembly2/shearing/dumps_v"+version+"/"
input_filenames = ["dump_v"+version+"_1em4.lammpstrj","dump_v"+version+"_1em3.lammpstrj","dump_v"+version+"_1em2.lammpstrj"]

input_filenames.remove("dump_v"+version+"_1em4.lammpstrj")
input_filenames.remove("dump_v"+version+"_1em3.lammpstrj")

colors = ['red','green','blue']
markers = ['x','o','s']
x = np.linspace(0,2,num=200)

plt.figure(1)
plt.title('Average non-affine square displacement. Sliding reference')
plt.figure(2)
plt.title('Average non-affine square displacement. Frame 0 reference')

for k,filename in enumerate(input_filenames):
    input_filepath = input_filefolder+filename
    node = import_file(input_filepath)

    #Obtener numero de iteraciones
    iter_num = node.source.num_frames

    #Seleccionar patches
    selectpatches = SelectTypeModifier(operate_on = "particles",
                                        property = "Particle Type",
                                        types = {2,4})
    node.modifiers.append(selectpatches)

    #Eliminar patches
    deletepatches = DeleteSelectedModifier()
    node.modifiers.append(deletepatches)

    #Non affine square displacement
    atomic_strain1 = AtomicStrainModifier(output_nonaffine_squared_displacements=True,
                                        affine_mapping = ReferenceConfigurationModifier.AffineMapping.ToCurrent,
                                        use_frame_offset = True,
                                        minimum_image_convention = True)
    atomic_strain2 = AtomicStrainModifier(output_nonaffine_squared_displacements=True,
                                        affine_mapping = ReferenceConfigurationModifier.AffineMapping.ToCurrent,
                                        minimum_image_convention = True)
    
    node.modifiers.append(atomic_strain1)

    Dsq_v_t = []
    for i in range(1,iter_num):
        data=node.compute(i)
        per_particle_dsq = data.particles["Nonaffine Squared Displacement"]
        Dsq_v_t.append(sum(per_particle_dsq)/len(per_particle_dsq))

    plt.figure(1)
    plt.plot(x,Dsq_v_t,c=colors[k])

    del(node.modifiers[2])
    node.modifiers.append(atomic_strain2)

    Dsq_v_t = []
    for i in range(1,iter_num):
        data=node.compute(i)
        per_particle_dsq = data.particles["Nonaffine Squared Displacement"]
        Dsq_v_t.append(sum(per_particle_dsq)/len(per_particle_dsq))

    plt.figure(2)
    plt.plot(x,Dsq_v_t,c=colors[k],marker=markers[k],linestyle='none',fillstyle='none')

plt.figure(1)
plt.legend(['1e-4','1e-3','1e-2'],title='$\dot{\gamma}$')
plt.xlabel('$\gamma$')
plt.ylabel('$D^2_{min}$')

plt.figure(2)
plt.legend(['1e-4','1e-3','1e-2'],title='$\dot{\gamma}$')
plt.xlabel('$\gamma$')
plt.ylabel('$D^2_{min}$')

plt.show()