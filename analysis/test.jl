using JSON
using StaticArrays
using Plots
using LinearAlgebra
using NearestNeighbors

using Profile
using BenchmarkTools

using DelimitedFiles

include("analysisModule.jl")
using .analysisModule

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




