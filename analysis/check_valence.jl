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

function read_system_atoms(datadir)

    jumplines = 1
    atomnum = 0
    open(datadir,"r") do data

        for line in readlines(data)
            jumplines+=1
            text = split(line," ")
            if line == "Atoms # full"
                break
            end
            if length(text)==2 
                if text[2] == "atoms"
                    atomnum = parse(Int64,text[1])
                end
            end
        end

    end

    centers_id = Vector{SVector{2,Int64}}()
    centers_coords = Vector{SVector{3,Float64}}()
    patches_id = Vector{Int64}()
    patches_coords = Vector{SVector{3,Float64}}()

    open(datadir,"r") do data

        for _ in 1:jumplines
            readline(data)
        end

        for _ in 1:atomnum
            aa = split(readline(data)," ")
            aa_id = [parse(Int64,aa[k]) for k in 1:length(aa) if k == 1 || k == 3]
            aa_coord = [parse(Float64,aa[k]) for k in 1:length(aa) if k>4 && k<8]
            if aa_id[2] == 1 || aa_id[2] == 3 #Centros
                push!(centers_id,SVector{2,Int64}(aa_id))
                push!(centers_coords,SVector{3,Float64}(aa_coord))
            elseif aa_id[2] == 2 || aa_id[2] == 4 #Patches
                push!(patches_id,aa_id[1])
                push!(patches_coords,SVector{3,Float64}(aa_coord))
            end

        end
    end     

    return centers_id, centers_coords, patches_id, patches_coords

end

function read_system_L(datadir)
    L = 0
    open(datadir,"r") do data

        for line in readlines(data)
            text = split(line," ")
            if length(text)==4 
                if text[3] == "xlo"
                    L = abs(parse(Float64,text[2]) - parse(Float64,text[1]))
                end
            end
        end

    end

    return L
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
    r = wrap_boundaries(r,L) #Mapeo a caja central
    return r
end

#Crear datos con copias periodicas 
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

function main()

    systemdir = ARGS[1]
    sigma = parse(Float64,ARGS[2])

    bondnum,bonds = read_system_bonds(systemdir) #Leer bonds del archivo system.data
    centers_id, centers_coords, patches_id, patches_coords = read_system_atoms(systemdir)
    L = read_system_L(systemdir)

    patches_id_periodic, patches_coords_periodic = periodic_data(patches_id,patches_coords,L) #Obtener datos con copias periodicas de los patches
    patches_tree = KDTree(patches_coords_periodic) #Crear KDTree de patches con condiciones de frontera periodicas

    r_c = 1.5*sigma

    neighbor_distr_mon = Vector{Int64}()
    neighbor_distr_xl = Vector{Int64}()
    aglomerations = 0
    aglom_id = Vector{Int64}()
    for i in eachindex(centers_id)
        neighnum = length(get_neigbors(centers_id[i][1],patches_id,patches_coords,patches_id_periodic,patches_tree,bonds,r_c))
        if centers_id[i][2] == 1
            push!(neighbor_distr_mon,neighnum)
            if neighnum>2
                aglomerations += 1
                push!(aglom_id,centers_id[i][1])
            end
        elseif centers_id[i][2] == 3 
            
            push!(neighbor_distr_xl,neighnum)
            if neighnum>4
                aglomerations += 1
                push!(aglom_id,centers_id[i][1])
            end
        end
        #print(aglomerations," ")
    end

    println()
    println("Aglomerations: ",aglomerations)
    
    histogram_mon = Vector{Int64}(undef,10)
    histogram_xl = Vector{Int64}(undef,10)
    for i in 1:10
        histogram_mon[i] = count(a->(a==i-1),neighbor_distr_mon)
        histogram_xl[i] = count(a->(a==i-1),neighbor_distr_xl)
    end
    println("Distribucion de vecinos de monomeros: ",histogram_mon)
    println("Distribucion de vecinos de cross-linkers: ",histogram_xl)
    if aglomerations>10
        println("Aglomerated particles: ",aglom_id[1:10])
    elseif aglomerations>0
        println("Aglomerated particles: ",aglom_id[1:end])
    end

end

main()