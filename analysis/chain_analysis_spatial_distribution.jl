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
    xyzdir = ARGS[3]

    println("Inicio chain analysis")

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = frame_num/2
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 1:calc_frames_step:frame_num
    #calc_frames = [1]

    bondnum,bonds = read_system_bonds(systemdir) #Leer bonds del archivo system.data

    dump = open(dumpdir,"r") #Abrir dump
    xyzfile = open(xyzdir,"w")

    for frame in 1:frame_num

        timestep,n,L,xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

        if frame in calc_frames

            #Triclinic box vectors
            a = SVector{3,Float64}([L,0,0])
            b = SVector{3,Float64}([xy,L,0])
            c = SVector{3,Float64}([0,0,L])

            println("Frame: ",frame)
            centers_coords,centers_id,patches_coords,patches_id = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
            centers_coords = fix_boundaries(centers_coords,xy,L) #Asegurar que todas las particulas esten dentro de la caja
            patches_coords = fix_boundaries(patches_coords,xy,L) #Asegurar que todas las particulas esten dentro de la caja
            patches_id_periodic, patches_coords_periodic = periodic_data(patches_id,patches_coords,a,b,c) #Obtener datos con copias periodicas de los patches
            patches_tree = KDTree(patches_coords_periodic) #Crear KDTree de patches con condiciones de frontera periodicas
            
            free_centers = [centers_id[i][1] for i in eachindex(centers_id)]

            xl_list = get_xl_list(centers_id)
            r_c = 0.5 #Radio cutoff para determinar vecinos patch-patch
            
            write(xyzfile, "$(length(centers_id))\n")
            write(xyzfile,"Lattice=\"$L 0 0 $xy $L 0 0 0 $L\" Origin=\"0.0 0.0 0.0\" Properties=type:I:1:pos:R:3\n")
            
            for c_xl in xl_list
                #print(c_xl, " ")
                
                if c_xl in free_centers
                    xlcoords = get_coords(c_xl,centers_id,centers_coords)
                    write(xyzfile,"5 $(xlcoords[1]) $(xlcoords[2]) $(xlcoords[3])\n")
                    remove!(free_centers,c_xl)
                end

                chain_starts = get_neigbors(c_xl,patches_id,patches_coords,patches_id_periodic,patches_tree,bonds,r_c)

                for mon_0 in chain_starts
                    if mon_0 in free_centers
                        xlEnd = -1
                        loop=0
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
                                        loop+=1 
                                    end
                                end
                            end
                            safe+=1
                            if safe > 100_000
                                println("SAFE")
                                break
                            end
                        end
    
                        type=0
                        if xlEnd==-1 #Si no termino en un cross linker (dangling)
                            type=1
                        elseif loop ==1 #Si la cadena es un loop
                            type=3
                        else #Si la cadena es linkeada
                            type=2
                        end

                        for cc in visited
                            if cc in free_centers
                                ccoords = get_coords(cc,centers_id,centers_coords)
                                write(xyzfile,"$type $(ccoords[1]) $(ccoords[2]) $(ccoords[3])\n")
                                remove!(free_centers,cc)
                            end
                        end
                    end
                end
            end

            for cc in free_centers
                ccoords = get_coords(cc,centers_id,centers_coords)
                write(xyzfile,"4 $(ccoords[1]) $(ccoords[2]) $(ccoords[3])\n")
            end

        else
            for _ in 1:n
                readline(dump)
            end
        end
    end
    close(dump)

end

main()