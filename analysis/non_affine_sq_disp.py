import os
os.environ['OVITO_GUI_MODE'] = '1'

#Importar librerias
from ovito.io import import_file, export_file
from ovito.modifiers import SelectTypeModifier,AtomicStrainModifier, DeleteSelectedModifier, UnwrapTrajectoriesModifier
from ovito.pipeline import ReferenceConfigurationModifier
import matplotlib.pyplot as plt
import numpy as np

input_filefolder = input("Folder directory: ")
input_filenames = []
legend = []
for i in range(3):
    file = input("Filename (type end to stop): ")
    if file == "end":
        break
    else:
        input_filenames.append(file)
        input('Legend: ')

colors = ['red','green','blue']
markers = ['x','o','s']
x = np.linspace(0,2,num=200)


for k,filename in enumerate(input_filenames):
    input_filepath = input_filefolder+'/'+filename
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

    #Unwrap trajectories
    unwrap = UnwrapTrajectoriesModifier()
    node.modifiers.append(unwrap)

    #Non affine square displacement
    atomic_strain1 = AtomicStrainModifier(output_nonaffine_squared_displacements=True,
                                        affine_mapping = ReferenceConfigurationModifier.AffineMapping.ToReference,
                                        use_frame_offset = False,
                                        minimum_image_convention = False)
    
    node.modifiers.append(atomic_strain1)

    Dsq_v_t = []
    for i in range(1,iter_num):
        data=node.compute(i)
        per_particle_dsq = data.particles["Nonaffine Squared Displacement"]
        Dsq_v_t.append(sum(per_particle_dsq)/len(per_particle_dsq))

    plt.figure(1)
    plt.plot(x,Dsq_v_t,c=colors[k])


plt.figure(1)
plt.title('Average non-affine square displacement. Frame 0 reference')
plt.legend(legend)
plt.xlabel('$\gamma$')
plt.ylabel('$D^2_{min}$')

plt.show()