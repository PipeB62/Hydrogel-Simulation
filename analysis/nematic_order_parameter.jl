using JSON
using StaticArrays
using Plots
using LinearAlgebra
using NearestNeighbors


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

function get_coords(c_id::Int64,centers_id::Vector{SVector{2,Int64}},centers_coords::Vector{SVector{3,Float64}})
    index=0
    for i in eachindex(centers_id)
        if centers_id[i][1] == c_id
            index = i
            break
        end
    end
    coords = centers_coords[index]
    return coords
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

function get_curvature(chainids::Vector{Int64},
    centers_coords::Vector{SVector{3,Float64}},
    centers_id::Vector{SVector{2,Int64}})

    R = Vector{SVector{3,Float64}}(undef,length(chainids))

    for k in eachindex(chainids)
        c_id = chainids[k]
        for i in eachindex(centers_id)
            if centers_id[i][1] == c_id
                R[k] = centers_coords[i]
                break
            end
        end
    end

    t = Vector{SVector{3,Float64}}(undef,length(chainids)-1)
    for i in 1:length(R)-1
        t[i]=(R[i+1]-R[i])/norm(R[i+1]-R[i])
    end

    kappa_ds = Vector{Float64}(undef,length(chainids)-2)
    for i in 1:length(t)-1
        kappa_ds[i]=sqrt(2*(1-cdot(t[i],t[i+1])))
    end

    kappa = sum(kappa_ds)

end

#Guarda los datos en un archivo json
function write_json(savedir,data,name)
    json_string = JSON.json(data)
    open("$(savedir)/$(name).json","w") do f
        JSON.print(f, json_string)
    end
end

function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

function get_orientation_vector(c_id::Int64,
    bonds::Array{Float64},
    patches_id::Vector{Int64},
    patches_coords::Vector{SVector{3,Float64}},
    a::SVector{3,Float64},
    b::SVector{3,Float64},
    c::SVector{3,Float64})

    #Encontrar patches de particula
    indexes1 = findall(item -> item == c_id, bonds[:,1])
    indexes2 = findall(item -> item == c_id, bonds[:,2])
    c_patches = bonds[indexes1,2]
    append!(c_patches,bonds[indexes2,1])

    #Extraer posicion de los patches
    v1 = SVector{3,Float64}([0,0,0])
    v2 = SVector{3,Float64}([0,0,0])
    check=0
    for i in eachindex(patches_id)
        if patches_id[i]==c_patches[1]
            v1=patches_coords[i]
            check+=1
        elseif patches_id[i]==c_patches[2]
            v2=patches_coords[i]
            check+=1
        end

        if check==2
            break
        end
    end

    #Crear vector de orientacion unitario m considerando fronteras periodicas
    m_all = Vector{SVector{3,Float64}}(undef,27)
    aa = 1
    for q in -1:1
        for p in -1:1
            for o in -1:1
                tr = o*a + p*b + q*c
                v2t = v2 + tr
                m_all[aa] = v2t-v1
                aa+=1
            end
        end
    end

    m = argmin(x->norm(x),m_all)
    m = m./norm(m)
    
    return m
end

function main()

    dumpdir = ARGS[1]
    systemdir = ARGS[2]
    savedir = ARGS[3]

    println("Inicio parametro de orden nematico")

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = frame_num/10
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 1:calc_frames_step:frame_num
    write_json(savedir,calc_frames,"nematic_order_parameter_calcframes")

    _,bonds = read_system_bonds(systemdir) #Leer bonds del archivo system.data

    dump = open(dumpdir,"r") #Abrir dump

    S_vt = Vector{Float64}()

    for frame in 1:frame_num

        _,n,L,xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

        if frame in calc_frames

            #Triclinic box vectors
            a = SVector{3,Float64}([L,0,0])
            b = SVector{3,Float64}([xy,L,0])
            c = SVector{3,Float64}([0,0,L])

            #println("Frame: ",frame)
            _,centers_id,patches_coords,patches_id = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
            patches_coords = fix_boundaries(patches_coords,xy,L) #Asegurar que todas las particulas esten dentro de la caja
            
            M = Array{Float64}(undef,3,3)
            nmons=0
            for i in eachindex(centers_id)
                c_id = centers_id[i][1]
                c_type = centers_id[i][2]
                if c_type==1
                    nmons+=1
                    m = get_orientation_vector(c_id,bonds,patches_id,patches_coords,a,b,c)
                    for alpha in 1:3
                        for beta in 1:3
                            M[alpha,beta]+=m[alpha]*m[beta]
                        end
                    end
                end
            end
            M = M./nmons

            Q = (3/2)*M - (1/2)*I

            eig = eigen(Q)
            index = findfirst(x->x==maximum(eig.values),eig.values)

            S = eig.values[index]
            n_vec = eig.vectors[:,index]

            println("S = ",S)
            push!(S_vt,S)

            #println("n = ",n_vec)

        else
            for _ in 1:n
                readline(dump)
            end
        end


    end
    close(dump)

    write_json(savedir,S_vt,"nematic_order_parameter")

end


main()