using JSON
using StaticArrays
using Plots
using LinearAlgebra

include("analysisModule.jl")
using .analysisModule

function main()

    dumpdir = ARGS[1] 

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = 50
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 1:calc_frames_step:frame_num
    
    N_particles = 16000
    L = 2*25.54

    dump = open(dumpdir,"r")

    D = Vector{Float64}()

    for frame in 1:frame_num

        timestep,n,L,xy = read_dump_info(dump)

        if frame in calc_frames

            println("Frame:","$(frame) ")

            r = read_dump_coords(dump,n) #Leer posiciones de particulas del dump
            r = fix_boundaries(r,xy,L) #Asegurar que todas las particulas esten dentro de la caja (condiciones periodicas) y mapeo a caja central. 

            push!(D,fractal_dimension(r,L))

        else
            for _ in 1:n
                readline(dump)
            end
        end

    end

    #display(scatter(1:length(D),D,show=true))
    #readline()

    dir = split(dumpdir,"/")
    dir[end] = "analysis_results"
    savedir = join(dir,"/")
    write_json(savedir,D,"fractal_dimension")
    write_json(savedir,calc_frames,"calc_frames_fractal_dimension")

end

main()

function test()
    dumpdir = ARGS[1] 

    frame_num = read_frame_num(dumpdir)
  
    N_particles = 16000
    L = 2*25.54

    dump = open(dumpdir,"r")

    for frame in 1:frame_num

        timestep,n,L,xy = read_dump_info(dump)

        if frame == 17

            println("Frame:","$(frame) ")

            r = read_dump_coords(dump,n) #Leer posiciones de particulas del dump
            r = fix_boundaries(r,xy,L) #Asegurar que todas las particulas esten dentro de la caja (condiciones periodicas) y mapeo a caja central. 

            fractal_dimension(r,L)

        else
            for _ in 1:n
                readline(dump)
            end
        end

    end

end

#test()