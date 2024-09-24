using JSON
using StaticArrays
using Plots
using LinearAlgebra
using NearestNeighbors

include("analysisModule.jl")
using .analysisModule

function main()


    dumpdir = ARGS[1]
    xyzdir = ARGS[2]
    remove_shear = ARGS[3]

    println("Inicio velocity analysis")

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = frame_num/100
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 2:calc_frames_step:frame_num

    dump = open(dumpdir,"r") #Abrir dump
    xyzfile = open(xyzdir,"w")

    #Leer primer frame (primer frame de referencia)

    timestep,n,L,r_xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

    r_centers_coords,r_centers_id,_,_ = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
    r_centers_coords = fix_boundaries(r_centers_coords,r_xy,L) #Asegurar que todas las particulas esten dentro de la caja
    if remove_shear=="yes"
        println("removing shear")
        r_centers_coords = shear(r_centers_coords,-r_xy,L)
    end

    dt = 0.001

    for frame in 2:frame_num

        timestep,n,L,c_xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

        #Triclinic box vectors
        if remove_shear=="yes"
            a = SVector{3,Float64}([L,0,0])
            b = SVector{3,Float64}([0,L,0])
            c = SVector{3,Float64}([0,0,L])
        else
            a = SVector{3,Float64}([L,0,0])
            b = SVector{3,Float64}([c_xy,L,0])
            c = SVector{3,Float64}([0,0,L])
        end

        c_centers_coords,c_centers_id,_,_ = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
        c_centers_coords = fix_boundaries(c_centers_coords,c_xy,L) #Asegurar que todas las particulas esten dentro de la caja
        if remove_shear=="yes"
            c_centers_coords = shear(c_centers_coords,-c_xy,L)
        end

        if frame in calc_frames
            print(frame," ")
            write(xyzfile, "$(length(c_centers_id))\n")
            if remove_shear=="yes"
                write(xyzfile,"Lattice=\"$L 0 0 0 $L 0 0 0 $L\" Origin=\"0.0 0.0 0.0\" Properties=id:I:1:pos:R:3:velo:R:3\n") #revirtiendo shear
            else
                write(xyzfile,"Lattice=\"$L 0 0 $c_xy $L 0 0 0 $L\" Origin=\"0.0 0.0 0.0\" Properties=id:I:1:pos:R:3:velo:R:3\n") #sin reverir shear
            end
        end

        centers_list = [c_centers_id[i][1] for i in eachindex(c_centers_id)]

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
                write(xyzfile,"$cid $(current_coords[1]) $(current_coords[2]) $(current_coords[3]) $(v[1]) $(v[2]) $(v[3])\n")
            end
            
        end

        #Frame actual ahora es el de referencia
        r_centers_coords = c_centers_coords
        r_centers_id = c_centers_id

    end
    close(dump)

end

main()