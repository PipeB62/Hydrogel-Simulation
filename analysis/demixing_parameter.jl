using JSON
using StaticArrays
using Plots

include("analysisModule.jl")
using .analysisModule

function main()

    dumpdir = ARGS[1] 
    N_l = parse(Int64,ARGS[2])

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = 50
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 1:calc_frames_step:frame_num
    
    N_particles = 16000
    L = 2*25.54
    rho_0 = N_particles/(L^3)

    dump = open(dumpdir,"r")

    phi = Vector{Float64}()

    for frame in 1:frame_num

        timestep,n,L,xy = read_dump_info(dump)

        if frame in calc_frames

            println("Frame:","$(frame) ")

            r = read_dump_coords(dump,n) #Leer posiciones de particulas del dump
            r = fix_boundaries(r,xy,L) #Asegurar que todas las particulas esten dentro de la caja (condiciones periodicas) y mapeo a caja central. 

            push!(phi,demixing_param(r,L,N_l,rho_0))

        else
            for _ in 1:n
                readline(dump)
            end
        end

    end

    #display(scatter(1:length(phi),phi,show=true))
    #readline()

    dir = split(dumpdir,"/")
    dir[end] = "analysis_results"
    savedir = join(dir,"/")
    write_json(savedir,phi,"demixing_param_$(N_l)")

end

main()




