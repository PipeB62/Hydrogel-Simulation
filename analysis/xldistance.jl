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

    println("Inicio xl distance")

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = 100
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 1:calc_frames_step:frame_num
    #calc_frames = [1]

    bondnum,bonds = read_system_bonds(systemdir) #Leer bonds del archivo system.data

    dump = open(dumpdir,"r") #Abrir dump

    mean_coordination = Vector{Float64}()
    std_coordination = Vector{Float64}()
    mean_distance = Vector{Float64}()
    std_distance = Vector{Float64}()

    #selected_xls = [2391,1631,2446,1481,1466,56,156,3201,2986,2911,2706,1901]
    #selected_avdistances_v_t = Vector{Vector{Float64}}()
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

            xl_coordinations = Vector{Int64}()
            av_distances = Vector{Float64}()
            #selected_avdistances = Vector{Float64}(undef,length(selected_xls))
            for c_xl in xl_list
                #print(c_xl, " ")
                queue = Vector{Int64}() #Cola
                visited = Vector{Int64}() #Nodos visitados
                push!(queue,c_xl)
                push!(visited,c_xl)
                
                linked_xls = Vector{Int64}()

                safe = 1
                while length(queue)>0

                    c_node = first(queue) #Obtener nodo actual con el primer elemento del queue
                    popfirst!(queue) #Quitar primer elemento del queue
                    c_neighs = get_neigbors(c_node,patches_id,patches_coords,patches_id_periodic,patches_tree,bonds,r_c)
                    for cn in c_neighs
                        if get_type(cn,centers_id) == "mon" && !(cn in visited)
                            push!(queue,cn)
                            push!(visited,cn)
                        elseif get_type(cn,centers_id) == "xl" && !(cn in visited) && !(cn in linked_xls)
                            push!(linked_xls,cn)
                            push!(visited,cn)
                        end
                    end
                    safe+=1
                    if safe > 100_000
                        println("SAFE")
                        break
                    end
                end

                #=
                for i in eachindex(selected_xls)
                    if c_xl == selected_xls[i]
                        av_distance = 0
                        for l_xl in linked_xls
                            av_distance += real_distance(c_xl,l_xl,centers_coords,centers_id,a,b,c)
                        end
                        av_distance = av_distance/length(linked_xls)
                        selected_avdistances[i] = av_distance
                    end
                end
                =#

                if length(linked_xls)>0
                    av_distance = 0
                    for l_xl in linked_xls
                        av_distance += real_distance(c_xl,l_xl,centers_coords,centers_id,a,b,c)
                        #push!(av_distances,real_distance(c_xl,l_xl,centers_coords,centers_id,L))
                    end
                    av_distance = av_distance/length(linked_xls)
                    push!(xl_coordinations,length(linked_xls))
                    push!(av_distances,av_distance)
                else
                    push!(xl_coordinations,0)
                    push!(av_distances,0.0)
                end
            end
            
            #println("Mean coordination: ",mean(xl_coordinations))
            c_mean_coordinations,c_std_coordinations = mean_std(xl_coordinations)
            push!(mean_coordination,c_mean_coordinations)
            push!(std_coordination,c_std_coordinations)

            c_mean_distances,c_std_distances = mean_std(av_distances)
            #println("Mean distance: ",c_mean_distances)
            #println("StD distance: ",c_std_distances)
            push!(mean_distance,c_mean_distances)
            push!(std_distance,c_std_distances)

            #push!(selected_avdistances_v_t,selected_avdistances)

        else
            for _ in 1:n
                readline(dump)
            end
        end
    end
    close(dump)

    #savedir = dir * "/" * "analysis_results"
    #println("Coordination: ",mean_coordination)
    write_json(savedir,mean_coordination,"mean_xl_coordination")
    #println("Distances: ",mean_distance)
    write_json(savedir,mean_distance,"mean_xl_distance")

    write_json(savedir,std_distance,"std_xl_distance")

    write_json(savedir,std_coordination,"std_xl_coordination")

    write_json(savedir,calc_frames,"xldistance_calcframes")

    #write_json(savedir,selected_avdistances_v_t,"selected_avdistances")

end

main()

