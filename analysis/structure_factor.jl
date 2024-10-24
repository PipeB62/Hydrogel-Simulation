using JSON
using StaticArrays
using Plots
using LinearAlgebra
using NearestNeighbors

include("analysisModule.jl")
using .analysisModule

function main()

    dumpdir = ARGS[1]
    savedir = ARGS[2]
    direccion = ARGS[3]

    println("Inicio Factor de estructura")

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = 10
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 2:calc_frames_step:frame_num

    calc_frames = [2,12,22,32,42,52]

    dump = open(dumpdir,"r") #Abrir dump

    #Leer primer frame (primer frame de referencia)
    timestep,n,L,r_xy = read_dump_info(dump) #Leer informacion del dump en el frame actual
    for _ in 1:n
        readline(dump)
    end

    #Caja central
    a_0 = SVector{3,Float64}([L,0,0])
    b_0 = SVector{3,Float64}([0,L,0])
    c_0 = SVector{3,Float64}([0,0,L])

    #sampleo de vectores k
    n_max = 100
    k_vectors = Vector{SVector{3,Float64}}()
    kvalues = (2*pi/L)*(1:n_max)

    for ni in 1:n_max 
        if direccion == "x"
            kvec = (2*pi/L)*SVector{3,Float64}([ni,0,0])
        elseif direccion == "y"
            kvec = (2*pi/L)*SVector{3,Float64}([0,ni,0])
        elseif direccion == "z"
            kvec = (2*pi/L)*SVector{3,Float64}([0,0,ni])
        end
        push!(k_vectors,kvec)
    end
    N_k = length(k_vectors)
    #println("N_k: ",N_k)

    S_vt = Vector{Vector{Float64}}()
    p_xy = r_xy
    flipcount = 0
    for frame in 2:frame_num

        timestep,n,L,c_xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

        #Obtener el tilt real (sin flip)
        if abs(c_xy-p_xy)>1
            flipcount+=1
        end
        real_xy =L*flipcount+c_xy

        S = Vector{Float64}()
        if frame in calc_frames

            print(real_xy," ")

            #Triclinic box vectors
            a = SVector{3,Float64}([L,0,0])
            b = SVector{3,Float64}([c_xy,L,0])
            c = SVector{3,Float64}([0,0,L])

            c_centers_coords,_,_,_ = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
            c_centers_coords = fix_boundaries(c_centers_coords,L,a,b,c) #Asegurar que todas las particulas esten dentro de la caja
            c_centers_coords = shear(c_centers_coords,-real_xy,L) #Revertir shear
            c_centers_coords = wrap_boundaries(c_centers_coords,a_0,b_0,c_0) #Mapeo a caja central

            N_p = length(c_centers_coords) #Numero de particulas (centros)

            for kvec in k_vectors
                rePhi = 0
                imPhi = 0
                for i in eachindex(c_centers_coords)
                    ri = c_centers_coords[i] #posicion en tiempo t

                    rePhi += cos(cdot(kvec,ri))
                    imPhi += sin(cdot(kvec,ri))
                end
                push!(S,(1/N_p)*(rePhi^2 + imPhi^2)) #Calcular
            end

            push!(S_vt,S)

        else
            for _ in 1:n
                readline(dump)
            end
        end

        p_xy = c_xy

    end

    close(dump)

    write_json(savedir,S_vt,"structure_factor_vt_$(direccion)")
    write_json(savedir,kvalues,"structure_factor_vt_kvalues_$(direccion)")

    #=
    plot(kvalues,S_vt[1],show=true)
    plot!(kvalues,S_vt[2],show=true)
    readline()
    =#

end

main()