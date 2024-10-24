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
    remove_shear = ARGS[3]

    println("Inicio velocity profile analysis")

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = 10
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 2:calc_frames_step:frame_num

    dump = open(dumpdir,"r") #Abrir dump

    #Leer primer frame (primer frame de referencia)

    timestep,n,L,r_xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

    #Caja central
    a_0 = SVector{3,Float64}([L,0,0])
    b_0 = SVector{3,Float64}([0,L,0])
    c_0 = SVector{3,Float64}([0,0,L])

    #Triclinic box vectors
    a = SVector{3,Float64}([L,0,0])
    b = SVector{3,Float64}([r_xy,L,0])
    c = SVector{3,Float64}([0,0,L])

    r_centers_coords,r_centers_id,_,_ = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
    r_centers_coords = fix_boundaries(r_centers_coords,L,a,b,c) #Asegurar que todas las particulas esten dentro de la caja
    if remove_shear=="yes"
        println("removing shear")
        r_centers_coords = shear(r_centers_coords,-r_xy,L) #Reverir shear
        r_centers_coords = wrap_boundaries(r_centers_coords,a_0,b_0,c_0)
    end

    n_bins = 50
    binstep = L/n_bins 
    bins = binstep/2:binstep:L

    dt = 0.001 * 1000

    velocity_profile_vt = Vector{Vector{Float64}}()
    p_xy = r_xy
    flipcount = 0
    for frame in 2:frame_num

        timestep,n,L,c_xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

        #Obtener el tilt real (sin flip)
        if abs(c_xy-p_xy)>1
            flipcount+=1
        end
        real_xy =L*flipcount+c_xy

        #Triclinic box vectors
        a = SVector{3,Float64}([L,0,0])
        b = SVector{3,Float64}([c_xy,L,0])
        c = SVector{3,Float64}([0,0,L])

        c_centers_coords,c_centers_id,_,_ = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
        c_centers_coords = fix_boundaries(c_centers_coords,L,a,b,c) #Asegurar que todas las particulas esten dentro de la caja
        if remove_shear=="yes"
            c_centers_coords = shear(c_centers_coords,-real_xy,L)
            c_centers_coords = wrap_boundaries(c_centers_coords,a_0,b_0,c_0)
        end

        centers_list = [c_centers_id[i][1] for i in eachindex(c_centers_id)]
        #println(length(centers_list))

        velocity_profile = zeros(Float64,n_bins)
        particles_per_bin = zeros(Int64,n_bins)
        for i in eachindex(centers_list)
            cid = centers_list[i]
            current_coords = c_centers_coords[i]

            prev_coords = SVector{3,Float64}([0,0,0])
            for j in eachindex(r_centers_id)
                rid = r_centers_id[j][1]
                if rid==cid
                    prev_coords = r_centers_coords[j]
                    break 
                end
            end
            
            ds_all = Vector{SVector{3,Float64}}(undef,27)
            aa = 1
            for q in -1:1
                for p in -1:1
                    for o in -1:1
                        tr = o*a + p*b + q*c
                        current_coords_t = current_coords + tr
                        ds_all[aa] = current_coords_t-prev_coords
                        aa+=1
                    end
                end
            end

            ds = argmin(x->norm(x),ds_all)
            v = ds./dt

            if frame in calc_frames
                deltas = [abs(current_coords[2]-k) for k in bins]
                bin_index = argmin(deltas)
                velocity_profile[bin_index] += v[1]
                particles_per_bin[bin_index] += 1
            end
        end

        if frame in calc_frames 
            println("Frame: ",frame)
            velocity_profile = velocity_profile./particles_per_bin
            push!(velocity_profile_vt,velocity_profile)
        end

        #Frame actual ahora es el de referencia
        r_centers_coords = c_centers_coords
        r_centers_id = c_centers_id
        p_xy = c_xy

    end
    close(dump)

    write_json(savedir,velocity_profile_vt,"velocity_profile")
    write_json(savedir,bins,"velocity_profile_bins")
    write_json(savedir,calc_frames,"velocity_profile_calcframes")

end

main()