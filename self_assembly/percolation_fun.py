def percolation(directory,sigma):
    from ovito.io import import_file
    from ovito.modifiers import ClusterAnalysisModifier, CreateBondsModifier, ReplicateModifier
    import numpy as np


    input_filepath=f"{directory}/system_check.data"

    node = import_file(input_filepath)

    node.modifiers.append(ReplicateModifier(adjust_box = True, num_x = 3, num_y = 3, num_z = 3))

    #Modificador para crear bonds patch-patch
    a = 1.5
    Rlim = sigma*a #Cutoff para crear bonds patch-patch
    create_bonds1 = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise) 
    create_bonds1.set_pairwise_cutoff(2,2,Rlim)
    create_bonds1.set_pairwise_cutoff(2,4,Rlim)
    node.modifiers.append(create_bonds1) #Agregar modificador para crear bonds

    #Modificador para numero de clusters
    clusters = ClusterAnalysisModifier(
        neighbor_mode=ClusterAnalysisModifier.NeighborMode.Bonding,
        sort_by_size = True)
    node.modifiers.append(clusters)

    data = node.compute()
    cluster_table = data.tables["clusters"]
    biggest_cluster_sz = cluster_table["Cluster Size"][0]
    n_particles = data.particles.count

    percentage = biggest_cluster_sz/n_particles
    
    return percentage
    
#print(percolation("/media/felipe/Files/Hydrogel_sim_experiments/SelfAssemby/exp0",0.25))





