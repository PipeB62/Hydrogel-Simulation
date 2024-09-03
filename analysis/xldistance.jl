using JSON
using StaticArrays
using Plots
using LinearAlgebra
using NearestNeighbors

using Profile
using BenchmarkTools

function cdot(u1::SVector{3,Float64}, u2::SVector{3,Float64}) 
    sum(u1.*u2)
end

function norm(u::SVector{3,Float64})
    return sqrt(cdot(u,u))
end

function mean(v::Vector{Float64})
    return sum(v)/length(v)
end

function mean(v::Vector{Int64})
    return sum(v)/length(v)
end

function mean_std(v::Vector{Float64})
    mean = sum(v)/length(v)
    var = 0 
    for a in v
        var += (a-mean)^2
    end
    var = var/(length(v)-1)
    std = sqrt(var)

    return mean,std
end

function mean_std(v::Vector{Int64})
    mean = sum(v)/length(v)
    var = 0 
    for a in v
        var += (a-mean)^2
    end
    var = var/(length(v)-1)
    std = sqrt(var)

    return mean,std
end

function read_frame_num(dumpdir)
    dump = open(dumpdir,"r") 
    frame_num = 0
    for line in readlines(dump)
        if line == "ITEM: TIMESTEP"
            frame_num += 1
        end
    end
    close(dump)
    return frame_num 
end

#Lee la infromacion del frame del dump
function read_dump_info(dump)
    readline(dump) #saltar texto ITEM: TIMESTEP
    timestep = readline(dump) #leer iteracion actual
    readline(dump) #saltar texto ITEM: NUMBER OF ATOMS
    n = parse(Int64,readline(dump)) #leer numero de atomos
    readline(dump) #saltar texto ITEM: BOX BOUNDS xy xz yz pp pp pp
    #Leer coords de caja en x. Leer xy
    xline = [parse(Float64,ii) for ii in split(readline(dump)," ")] 
    xy = xline[3]
    #Leer coords de caja en y. Leer yz
    yline = [parse(Float64,ii) for ii in split(readline(dump)," ")]
    L = yline[2]-yline[1]
    readline(dump) #Leer coords de caja en z. Leer xz
    readline(dump) #saltar texto ITEM: ATOMS id type x y z

    return timestep,n,L,xy
end

#Lee las coordenads de las particulas del dump para un frame
function read_dump_particles(dump,n::Int64)
    centers_coords = Vector{SVector{3,Float64}}()
    patches_coords = Vector{SVector{3,Float64}}()
    centers_id = Vector{SVector{2,Int64}}()
    patches_id = Vector{Int64}()
    for _ in 1:n
        line = split(readline(dump)," ")
        aa = SVector{3,Float64}(parse.(Float64,line[3:5])) 
        if line[2] == "2" || line[2] == "4"
            bb = parse(Int64,line[1])
            push!(patches_coords,aa)
            push!(patches_id,bb)
        else
            cc = SVector{2,Float64}(parse.(Float64,line[1:2]))
            push!(centers_coords,aa)
            push!(centers_id,cc)
        end
    end
    return centers_coords,centers_id,patches_coords,patches_id
end

#Lee la informacion de los bonds del system.data
function read_system_bonds(datadir)

    jumplines = 1
    bondnum = 0
    open(datadir,"r") do data

        for line in readlines(data)
            jumplines+=1
            text = split(line," ")
            if line == "Bonds"
                break
            end
            if length(text)==2 
                if text[2] == "bonds"
                    bondnum = parse(Int64,text[1])
                end
            end
        end

    end

    bonds = Array{Float64}(undef,bondnum,2)

    open(datadir,"r") do data

        for _ in 1:jumplines
            readline(data)
        end

        for i in 1:bondnum
            aa = [parse(Float64,a) for a in split(readline(data)," ")]
            bonds[i,:] = aa[3:end] 
        end
    end     

    return bondnum,bonds

end

#Obtiene la lista de ids de cross linkers 
function get_xl_list(centers_id)
    xl_list = Vector{Int64}()
    for i in eachindex(centers_id)
        if centers_id[i][2] == 3
            push!(xl_list,centers_id[i][1])
        end
    end
    return xl_list
end

#Aplica una traslacion diagonal para mover el centro a la esquina de la caja
function traslacion(r::Vector{SVector{3,Float64}},T::Float64)
    for i in eachindex(r)    
        r[i] = r[i] + SVector{3,Float64}(T,T,T)
    end
    return r
end

#Aplica el shear a todas las particulas
function shear(r::Vector{SVector{3,Float64}},xy::Float64,L::Float64)
    M = SMatrix{3,3,Float64}(1,0,0,xy/L,1,0,0,0,1)
    for i in eachindex(r)
        r[i] = M*r[i]
    end
    return r
end

#Asegura que todas las particulas esten dentro de la caja
function wrap_boundaries(r::Vector{SVector{3,Float64}},L::Float64)
    for i in eachindex(r)
        ri = r[i]
        t = MVector{3,Float64}(0.0,0.0,0.0)
        for j in 1:3
            if ri[j] >= L 
                t[j] = -L 
            elseif ri[j] < 0
                t[j] = L 
            end
        end
        r[i] = ri + t
    end
    return r
end

#Hace el shear inverso, asegura que todas las particulas esten dentro de la caja y regresa el shear.
function fix_boundaries(r::Vector{SVector{3,Float64}},xy::Float64,L::Float64)
    r = traslacion(r,L/2) #Mover origen a orilla de la caja
    r = shear(r,-xy,L) #Aplicar el shear inverso
    r = wrap_boundaries(r,L) #Asegurar que todas las particulas esten dentro de la caja
    r = shear(r,xy,L) #Revertir el shear
    return r
end

#Crear datos con copias periodicas (centros)
function periodic_data(id::Vector{SVector{2,Int64}},coords::Vector{SVector{3,Float64}},a::SVector{3,Float64},b::SVector{3,Float64},c::SVector{3,Float64})

    coords_periodic = Vector{SVector{3,Float64}}()
    id_periodic = Vector{SVector{2,Int64}}()
    
    for q in -1:1
        for p in -1:1
            for o in -1:1
                tr = o*a + p*b + q*c
                for i in eachindex(coords)
                    push!(coords_periodic,coords[i]+tr)
                    push!(id_periodic,id[i])
                end
            end
        end
    end

    #=
    a=1
    for n in 0:26
        println(id_periodic[a+n*length(coords)],coords_periodic[a+n*length(coords)])
    end
    =#

    return id_periodic, coords_periodic

end

#Crear datos con copias periodicas (patches)
function periodic_data(id::Vector{Int64},coords::Vector{SVector{3,Float64}},a::SVector{3,Float64},b::SVector{3,Float64},c::SVector{3,Float64})

    coords_periodic = Vector{SVector{3,Float64}}()
    id_periodic = Vector{Int64}()
    
    for q in -1:1
        for p in -1:1
            for o in -1:1
                tr = tr = o*a + p*b + q*c
                for i in eachindex(coords)
                    push!(coords_periodic,coords[i]+tr)
                    push!(id_periodic,id[i])
                end
            end
        end
    end

    #=
    a=1
    for n in 0:26
        println(id_periodic[a+n*length(coords)],coords_periodic[a+n*length(coords)])
    end
    =#

    return id_periodic, coords_periodic

end

#Obtiene los vecinos unidos por interaccion patch-patch de una particula dada
function get_neigbors(c_id::Int64,
    patches_id::Vector{Int64},
    patches_coords::Vector{SVector{3,Float64}},
    patches_id_periodic::Vector{Int64},
    #patches_coords_periodic::Vector{SVector{3,Float64}},
    patches_tree,
    bonds::Array{Float64},
    r_c::Float64)

    #Encontrar patches de particula
    #println("c_id: ",c_id)
    indexes1 = findall(item -> item == c_id, bonds[:,1])
    indexes2 = findall(item -> item == c_id, bonds[:,2])
    c_patches = bonds[indexes1,2]
    append!(c_patches,bonds[indexes2,1])
    #println("c_patches: ",c_patches)

    #Usar kd-tree para encontrar vecinos del patch (patches unidos)
    neigh_patches = Vector{Int64}()

    for aa in eachindex(c_patches)

        cpatch = c_patches[aa]
        point = SVector{3,Float64}(0,0,0)
        for i in eachindex(patches_id)
            if patches_id[i] == cpatch
                point = patches_coords[i]
                break
            end
        end
        neigh_patches_ix = inrange(patches_tree,point,r_c)
        #println(neigh_patches_ix)
        for i in neigh_patches_ix
            if patches_id_periodic[i] != cpatch
                push!(neigh_patches,patches_id_periodic[i])
            end
        end    

    end
    #println("neigh_patches: ",neigh_patches)

    #Encontrar centros correspondiente a los patches vecinos
    neighs = Vector{Int64}()
    for cnp in neigh_patches
        indexes1 = findall(item -> item == cnp, bonds[:,1])
        indexes2 = findall(item -> item == cnp, bonds[:,2])
        c_centers = bonds[indexes1,2]
        append!(c_centers,bonds[indexes2,1])
        append!(neighs,c_centers)
    end

    #println("neighs: ",neighs)

    return neighs

end

function get_type(c_id::Int64,centers_id::Vector{SVector{2,Int64}})
    numtype = 0
    for i in eachindex(centers_id)
        if centers_id[i][1] == c_id
            numtype =  centers_id[i][2]
            break
        end
    end
    if numtype == 1
        return "mon"
    elseif numtype == 3
        return "xl"
    end
end

function real_distance(c_id1::Int64,
    c_id2::Int64,
    centers_coords::Vector{SVector{3,Float64}},
    centers_id::Vector{SVector{2,Int64}},
    a::SVector{3,Float64},
    b::SVector{3,Float64},
    c::SVector{3,Float64})

    r1 = SVector{3,Float64}(0,0,0)
    for i in eachindex(centers_id)
        if centers_id[i][1] == c_id1
            r1 = centers_coords[i]
            break
        end
    end

    r2 = SVector{3,Float64}(0,0,0)
    for i in eachindex(centers_id)
        if centers_id[i][1] == c_id2
            r2 = centers_coords[i]
            break
        end
    end

    ds = Vector{Float64}(undef,27)
    aa = 1
    for q in -1:1
        for p in -1:1
            for o in -1:1
                tr = o*a + p*b + q*c
                r2t = r2 + tr
                ds[aa] = norm(r1-r2t)
                aa+=1
            end
        end
    end

    return minimum(ds)

end

#Guarda los datos en un archivo json
function write_json(savedir,data,name)
    json_string = JSON.json(data)
    open("$(savedir)/$(name).json","w") do f
        JSON.print(f, json_string)
    end
end

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

