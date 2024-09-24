using JSON
using StaticArrays
using Plots
using LinearAlgebra
using NearestNeighbors

include("analysisModule.jl")
using .analysisModule


function main()

    dumpdir = ARGS[1]
    systemdir = ARGS[2]
    xyzdir = ARGS[3]

    println("Inicio local order parameter")

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = 10
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 1:calc_frames_step:frame_num

    _,bonds = read_system_bonds(systemdir) #Leer bonds del archivo system.data

    dump = open(dumpdir,"r") #Abrir dump
    xyzfile = open(xyzdir,"w")

    rho_ij = 5.0

    for frame in 1:frame_num

        timestep,n,L,xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

        #Triclinic box vectors
        a = SVector{3,Float64}([L,0,0])
        b = SVector{3,Float64}([xy,L,0])
        c = SVector{3,Float64}([0,0,L])

        centers_coords,centers_id,patches_coords,patches_id = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
        centers_coords = fix_boundaries(centers_coords,xy,L) #Asegurar que todas las particulas esten dentro de la caja

        if frame in calc_frames
            print(frame," ")
            write(xyzfile, "$(length(centers_id))\n")
            write(xyzfile,"Lattice=\"$L 0 0 $xy $L 0 0 0 $L\" Origin=\"0.0 0.0 0.0\" Properties=id:I:1:pos:R:3:s:R:1\n") #sin reverir shear
            
            for i in eachindex(centers_id)

                s = 0
                n_ij = 0

                i_id = centers_id[i][1]
                i_coords = centers_coords[i]
                itype = centers_id[i][2]

                if itype==1
                    for j in eachindex(centers_id)
                        j_id = centers_id[j][1]
                        j_coords = centers_coords[j]
                        jtype = centers_id[j][2]

                        if !(i_id==j_id) && jtype==1
                            delta = real_distance(i_coords,j_coords,a,b,c)
                            if delta<rho_ij
                                i_director = get_orientation_vector(i_id,bonds,patches_id,patches_coords,a,b,c)
                                j_director = get_orientation_vector(j_id,bonds,patches_id,patches_coords,a,b,c)
                                s += (3/2)*abs(cdot(i_director,j_director))-1/2
                                n_ij += 1
                            end
                        end
                    end
                    s = s/n_ij
                end

                write(xyzfile,"$i_id $(i_coords[1]) $(i_coords[2]) $(i_coords[3]) $s\n")
            end
        end
    end
    close(dump)

end

main()