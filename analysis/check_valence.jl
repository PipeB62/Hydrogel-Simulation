using JSON
using StaticArrays
using Plots
using LinearAlgebra
using NearestNeighbors

include("analysisModule.jl")
using .analysisModule

function main()

    systemdir = ARGS[1]
    sigma = parse(Float64,ARGS[2])

    bondnum,bonds = read_system_bonds(systemdir) #Leer bonds del archivo system.data
    centers_id, centers_coords, patches_id, patches_coords = read_system_atoms(systemdir)
    L = read_system_L(systemdir)

    patches_id_periodic, patches_coords_periodic = periodic_data(patches_id,patches_coords,L) #Obtener datos con copias periodicas de los patches
    patches_tree = KDTree(patches_coords_periodic) #Crear KDTree de patches con condiciones de frontera periodicas

    r_c = 1.3*sigma

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