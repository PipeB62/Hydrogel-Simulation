using JSON
using StaticArrays
using Plots
using LinearAlgebra

function read_frame_num(dumpdir)
    dump = open(dumpdir,"r") 
    frame_num = 0
    for line in readlines(dump)
        if line == "ITEM: TIMESTEP"
            frame_num += 1
        end
    end
    close(dump)
    return frame_num 
end

#Lee la infromacion del frame del dump
function read_dump_info(dump)
    readline(dump) #saltar texto ITEM: TIMESTEP
    timestep = readline(dump) #leer iteracion actual
    readline(dump) #saltar texto ITEM: NUMBER OF ATOMS
    n = parse(Int64,readline(dump)) #leer numero de atomos
    readline(dump) #saltar texto ITEM: BOX BOUNDS xy xz yz pp pp pp
    #Leer coords de caja en x. Leer xy
    xline = [parse(Float64,ii) for ii in split(readline(dump)," ")] 
    xy = xline[3]
    #Leer coords de caja en y. Leer yz
    yline = [parse(Float64,ii) for ii in split(readline(dump)," ")]
    L = yline[2]-yline[1]
    readline(dump) #Leer coords de caja en z. Leer xz
    readline(dump) #saltar texto ITEM: ATOMS id type x y z

    return timestep,n,L,xy
end

#Lee las coordenads de las particulas del dump para un frame
function read_dump_coords(dump,n::Int64)
    r = Vector{SVector{3,Float64}}()
    for _ in 1:n
        line = split(readline(dump)," ")
        if line[2] == "2" || line[2] == "4"
            continue
        else
            ri = @SVector [parse(Float64,line[ii]) for ii in 3:5]
            push!(r,ri)
        end
    end
    return r
end

#Aplica una traslacion diagonal para mover el centro a la esquina de la caja
function traslacion(r::Vector{SVector{3,Float64}},T::Float64)
    for i in eachindex(r)    
        r[i] = r[i] + SVector{3,Float64}(T,T,T)
    end
    return r
end

#Aplica el shear a todas las particulas
function shear(r::Vector{SVector{3,Float64}},xy::Float64,L::Float64)
    M = SMatrix{3,3,Float64}(1,0,0,xy/L,1,0,0,0,1)
    for i in eachindex(r)
        r[i] = M*r[i]
    end
    return r
end

#Asegura que todas las particulas esten dentro de la caja
function wrap_boundaries(r::Vector{SVector{3,Float64}},L::Float64)
    for i in eachindex(r)
        ri = r[i]
        t = MVector{3,Float64}(0.0,0.0,0.0)
        for j in 1:3
            if ri[j] >= L 
                t[j] = -L 
            elseif ri[j] < 0
                t[j] = L 
            end
        end
        r[i] = ri + t
    end
    return r
end

#Hace el shear inverso, asegura que todas las particulas esten dentro de la caja y regresa el shear.
function fix_boundaries(r::Vector{SVector{3,Float64}},xy::Float64,L::Float64)
    r = traslacion(r,L/2) #Mover origen a orilla de la caja
    r = shear(r,-xy,L) #Aplicar el shear inverso
    r = wrap_boundaries(r,L) #Asegurar que todas las particulas esten dentro de la caja
    r = shear(r,xy,L) #Revertir el shear
    r = wrap_boundaries(r,L) #Mapeo a caja central
    return r
end

#Guarda los datos en un archivo json
function write_json(savedir,data,name)
    json_string = JSON.json(data)
    open("$(savedir)/$(name).json","w") do f
        JSON.print(f, json_string)
    end
end

function fractal_dimension(r::Vector{SVector{3,Float64}},L::Float64)

    N_ls = Vector{Int64}(2:15)

    N_b = Vector{Int64}(undef,length(N_ls))
    eps = Vector{Float64}(undef,length(N_ls))

    i=1
    for N_l in N_ls

        step = L/N_l 
        spots = step/2:step:L

        n_particles_box = zeros(Int64,N_l^3)

        for ri in r
            boxlocator = Vector{Int64}(undef,3)
            for j in 1:3
                boxlocator[j] = findmin([abs(k-ri[j]) for k in spots])[2]
            end
            loc = boxlocator[1] + N_l*(boxlocator[2]-1) + N_l^2 * (boxlocator[3]-1)
            n_particles_box[loc] += 1
        end

        empty_boxes = 0
        for aa in n_particles_box
            if aa == 0
                empty_boxes += 1
            end
        end

        N_b[i] = length(n_particles_box)-empty_boxes
        eps[i] = step

        i+=1

    end

    Y = [log(N) for N in N_b]
    A = ones(length(N_ls),2)
    A[:,1] = [log(1/e) for e in eps]

    X_ls = inv(A'*A)*A'*Y

    #test
    #display(scatter([log(1/e) for e in eps],Y,show=true))
    #xx = [log(1/e) for e in eps]
    #yy = X_ls[1].*xx .+ X_ls[2]
    #display(plot!(xx,yy,show=true))
    #readline()
    #test

    return abs(X_ls[1])

end

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