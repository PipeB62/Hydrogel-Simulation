using JSON
using StaticArrays
using Plots
using LinearAlgebra
using NearestNeighbors

using Profile
using BenchmarkTools

using DelimitedFiles

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
    #r = wrap_boundaries(r,L) #Mapeo a caja central
    return r
end

#Crear datos con copias periodicas (centros)
function periodic_data(id::Vector{SVector{2,Int64}},coords::Vector{SVector{3,Float64}},L::Float64)

    coords_periodic = Vector{SVector{3,Float64}}()
    id_periodic = Vector{SVector{2,Int64}}()
    
    for q in -1:1
        for p in -1:1
            for o in -1:1
                tr = SVector{3,Float64}([o*L,p*L,q*L])
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
function periodic_data(id::Vector{Int64},coords::Vector{SVector{3,Float64}},L::Float64)

    coords_periodic = Vector{SVector{3,Float64}}()
    id_periodic = Vector{Int64}()
    
    for q in -1:1
        for p in -1:1
            for o in -1:1
                tr = SVector{3,Float64}([o*L,p*L,q*L])
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
        #patches_tree = KDTree(patches_coords_periodic)
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
    L::Float64)

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
                tr = SVector{3,Float64}([o*L,p*L,q*L])
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

    dir = ARGS[1]
    dumpp = ARGS[2]

    dumpdir = dir * "/" * dumpp

    frame_num = read_frame_num(dumpdir)

    calc_frame = 1

    dump = open(dumpdir,"r") #Abrir dump

    for frame in 1:frame_num

        timestep,n,L,xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

        if frame == calc_frame
            println("Frame: ",frame)
            println("xy = ",xy)
            centers_coords,centers_id,patches_coords,patches_id = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
            
            ix = 1 
            for ii in centers_id
                if ii[1]==48225
                    break
                else
                    ix+=1
                end
            end
            println(centers_id[ix])
            println(centers_coords[ix])

            centers_coords_x = [i[1] for i in centers_coords]
            centers_coords_y = [i[2] for i in centers_coords]
            centers_coords_z = [i[3] for i in centers_coords]
            patches_coords_x = [i[1] for i in patches_coords]
            patches_coords_y = [i[2] for i in patches_coords]
            patches_coords_z = [i[3] for i in patches_coords]

            centers_id_w = Int64[i[1] for i in centers_id]

            open(dir*"/particles_pre.xyz", "w") do file
                writedlm(file, n)
                write(file,"Lattice=\"$L 0 0 $xy $L 0 0 0 $L\" Origin=\"$(-L/2) $(-L/2) $(-L/2)\" Properties=type:I:1:id:I:1:pos:R:3\n")
                
                writedlm(file, Union{Int64,Float64}[ones(Int,length(centers_coords)) centers_id_w centers_coords_x centers_coords_y centers_coords_z])
                writedlm(file, Union{Int64,Float64}[2*ones(Int,length(patches_coords)) patches_id patches_coords_x patches_coords_y patches_coords_z])
            end

            centers_coords = fix_boundaries(centers_coords,xy,L) #Mapeo a caja central
            patches_coords = fix_boundaries(patches_coords,xy,L) #Mapeo a caja central

            centers_coords_x = [i[1] for i in centers_coords]
            centers_coords_y = [i[2] for i in centers_coords]
            centers_coords_z = [i[3] for i in centers_coords]
            patches_coords_x = [i[1] for i in patches_coords]
            patches_coords_y = [i[2] for i in patches_coords]
            patches_coords_z = [i[3] for i in patches_coords]

            open(dir*"/particles_post.xyz", "w") do file
                writedlm(file, n)
                write(file,"Lattice=\"$L 0 0 $xy $L 0 0 0 $L\" Properties=type:I:1:id:I:1:pos:R:3\n")
                
                writedlm(file, Union{Int64,Float64}[ones(Int,length(centers_coords)) centers_id_w centers_coords_x centers_coords_y centers_coords_z])
                writedlm(file, Union{Int64,Float64}[2*ones(Int,length(patches_coords)) patches_id patches_coords_x patches_coords_y patches_coords_z])
            end
            
            break

        else
            for _ in 1:n
                readline(dump)
            end
        end
    end
    close(dump)

end

main()




