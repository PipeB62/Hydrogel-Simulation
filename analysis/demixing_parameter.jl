using JSON
using StaticArrays
using Plots

include("analysisModule.jl")
using .analysisModule

function main()

    dumpdir = ARGS[1] 
    savedir = ARGS[2]
    N_l = parse(Int64,ARGS[3])

    frame_num = read_frame_num(dumpdir)

    dump = open(dumpdir,"r")

    calc_frames_num = 50
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 2:calc_frames_step:frame_num

    #Leer primer frame (referencia)
    timestep,n,L,r_xy = read_dump_info(dump)
    r_centers_coords,_,_,_ = read_dump_particles(dump,n) #Leer posiciones de particulas del dump

    #Caja central
    a_0 = SVector{3,Float64}([L,0,0])
    b_0 = SVector{3,Float64}([0,L,0])
    c_0 = SVector{3,Float64}([0,0,L])
    
    N_particles = length(r_centers_coords)
    rho_0 = N_particles/(L^3)

    phi = Vector{Float64}()
    p_xy = r_xy
    flipcount=0
    for frame in 2:frame_num   

        timestep,n,L,c_xy = read_dump_info(dump)

        #Obtener el tilt real (sin flip)
        if abs(c_xy-p_xy)>1
            flipcount+=1
        end
        real_xy = L*flipcount+c_xy

        #Triclinic box vectors
        a = SVector{3,Float64}([L,0,0])
        b = SVector{3,Float64}([c_xy,L,0])
        c = SVector{3,Float64}([0,0,L])

        if frame in calc_frames

            println("Frame:","$(frame) ")

            c_centers_coords,_,_,_ = read_dump_particles(dump,n) #Leer posiciones de particulas del dump
            c_centers_coords = fix_boundaries(c_centers_coords,L,a,b,c) #Asegurar que todas las particulas esten dentro de la caja
            c_centers_coords = shear(c_centers_coords,-real_xy,L) #Revertir shear
            c_centers_coords = wrap_boundaries(c_centers_coords,a_0,b_0,c_0) #Mapeo a caja central

            push!(phi,demixing_param(c_centers_coords,L,N_l,rho_0))

        else
            for _ in 1:n
                readline(dump)
            end
        end

        p_xy = c_xy

    end

    scatter(1:length(phi),phi,show=true)
    readline()

    #write_json(savedir,phi,"demixing_param_$(N_l)")

end

main()




