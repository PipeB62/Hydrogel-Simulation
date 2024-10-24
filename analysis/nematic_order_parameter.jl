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
            patches_coords = fix_boundaries(patches_coords,L,a,b,c) #Asegurar que todas las particulas esten dentro de la caja
            
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