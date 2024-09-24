using JSON
using StaticArrays
using Plots
using LinearAlgebra
using NearestNeighbors

include("analysisModule.jl")
using .analysisModule

function main()

    #dir = ARGS[1]
    #dumpp = ARGS[2]
    #system = ARGS[3]

    #dumpdir = dir * "/" * dumpp
    #systemdir = dir * "/" * system

    dumpdir = ARGS[1]
    systemdir = ARGS[2]
    savedir = ARGS[3]

    println("Inicio chain analysis")

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = 100
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 1:calc_frames_step:frame_num
    #calc_frames = [1]

    bondnum,bonds = read_system_bonds(systemdir) #Leer bonds del archivo system.data

    dump = open(dumpdir,"r") #Abrir dump

    mean_linked_chain_length_v_t = Vector{Float64}()
    mean_dangling_chain_length_v_t = Vector{Float64}()
    mean_total_chain_length_v_t = Vector{Float64}()
    std_total_chain_length_v_t = Vector{Float64}()
    dangling_chains_num_v_t = Vector{Float64}()
    linked_chains_num_v_t = Vector{Float64}()
    mean_curvature_v_t = Vector{Float64}()
    loopnum_v_t = Vector{Float64}()

    for frame in 1:frame_num

        timestep,n,L,xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

        if frame in calc_frames

            #Triclinic box vectors
            a = SVector{3,Float64}([L,0,0])
            b = SVector{3,Float64}([xy,L,0])
            c = SVector{3,Float64}([0,0,L])

            #println("Frame: ",frame)
            centers_coords,centers_id,patches_coords,patches_id = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
            centers_coords = fix_boundaries(centers_coords,xy,L) #Asegurar que todas las particulas esten dentro de la caja
            patches_coords = fix_boundaries(patches_coords,xy,L) #Asegurar que todas las particulas esten dentro de la caja
            patches_id_periodic, patches_coords_periodic = periodic_data(patches_id,patches_coords,a,b,c) #Obtener datos con copias periodicas de los patches
            patches_tree = KDTree(patches_coords_periodic) #Crear KDTree de patches con condiciones de frontera periodicas
            #centers_id_periodic, centers_coords_periodic = periodic_data(patches_id,patches_coords,L) #Obtener datos con copias periodicas de los centros

            xl_list = get_xl_list(centers_id)
            r_c = 0.5 #Radio cutoff para determinar vecinos patch-patch
            
            loops = 0
            dangling_chains = 0
            linked_chain_length_histogram = zeros(Int64,100)
            chain_lengths_total = Vector{Float64}()
            chain_lengths_linked = Vector{Float64}()
            chain_lengths_dangling = Vector{Float64}()
            curvature = Vector{Float64}()
            for c_xl in xl_list
                #print(c_xl, " ")

                chain_starts = get_neigbors(c_xl,patches_id,patches_coords,patches_id_periodic,patches_tree,bonds,r_c)

                for mon_0 in chain_starts
                    
                    xlEnd = -1
                    queue = Vector{Int64}() #Cola
                    visited = Vector{Int64}() #Nodos visitados
                    push!(queue,mon_0)
                    push!(visited,mon_0)

                    safe = 1
                    while length(queue)>0

                        c_node = first(queue) #Obtener nodo actual con el primer elemento del queue
                        popfirst!(queue) #Quitar primer elemento del queue
                        c_neighs = get_neigbors(c_node,patches_id,patches_coords,patches_id_periodic,patches_tree,bonds,r_c) #obtener vecinos del nodo actual
                        if c_node == mon_0
                            deleteat!(c_neighs, findall(x->x==c_xl, c_neighs)) #qutar c_xl de la lista de vecinos para el primer nodo
                        end
                        #BFS
                        for cn in c_neighs
                            if get_type(cn,centers_id) == "mon" && !(cn in visited)
                                push!(queue,cn)
                                push!(visited,cn)
                            elseif get_type(cn,centers_id) == "xl" && !(cn in visited) 
                                xlEnd=cn
                                if cn == c_xl
                                    loops+=1 
                                    #println("loop ",c_xl," ",mon_0)
                                end
                            end
                        end
                        safe+=1
                        if safe > 100_000
                            println("SAFE")
                            break
                        end
                    end

                    s_chain = 2^(1/6)*length(visited)#Longitud de la cadena. +2 para contar los cross linkers

                    push!(chain_lengths_total,s_chain) 

                    if xlEnd==-1 #Si no termino en un cross linker (dangling)
                        dangling_chains+=1
                        push!(chain_lengths_dangling,s_chain)
                    else #Si la cadena une dos cross linkers
                        push!(chain_lengths_linked,s_chain)
                        linked_chain_length_histogram[length(visited)]+=1

                        if length(visited)>3
                            kappa = get_curvature(visited,centers_coords,centers_id)
                            push!(curvature,kappa)
                        end
                    end
                end
            end

            linked_chain_length_histogram = linked_chain_length_histogram./2 #Las cadenas linked se estan contando doble

            #println("Loops: ",loops)

            #println("Dangling_chains: ",dangling_chains)
            push!(dangling_chains_num_v_t,dangling_chains)

            #println("Linked chains: ",sum(linked_chain_length_histogram))
            push!(linked_chains_num_v_t,sum(linked_chain_length_histogram))

            #println("Linked chain length histogram: ",linked_chain_length_histogram)

            #println("Mean curvature: ",mean(curvature))
            push!(mean_curvature_v_t,mean(curvature))

            #println("Mean total chain length: ", mean(chain_lengths_total))
            c_mean_total_chain_length, c_std_total_chain_length = mean_std(chain_lengths_total)
            push!(mean_total_chain_length_v_t,c_mean_total_chain_length)
            push!(std_total_chain_length_v_t,c_std_total_chain_length)

            #println("Mean dangling chain length: ", mean(chain_lengths_dangling))
            push!(mean_dangling_chain_length_v_t,mean(chain_lengths_dangling))

            #println("Mean linked chain length: ", mean(chain_lengths_linked))
            push!(mean_linked_chain_length_v_t,mean(chain_lengths_linked))

            push!(loopnum_v_t,loops/2)

        else
            for _ in 1:n
                readline(dump)
            end
        end
    end
    close(dump)

    #savedir = dir * "/" * "analysis_results"

    write_json(savedir,loopnum_v_t,"loop_num")
    write_json(savedir,dangling_chains_num_v_t,"dangling_chains_num")
    write_json(savedir,linked_chains_num_v_t,"linked_chains_num")
    write_json(savedir,mean_curvature_v_t,"mean_curvature")
    write_json(savedir,mean_total_chain_length_v_t,"mean_total_chain_length")
    write_json(savedir,std_total_chain_length_v_t,"std_total_chain_length")
    write_json(savedir,mean_dangling_chain_length_v_t,"mean_dangling_chain_length")
    write_json(savedir,mean_linked_chain_length_v_t,"mean_linked_chain_length")
    write_json(savedir,calc_frames,"chain_analysis_calcframes")
    
end

main()