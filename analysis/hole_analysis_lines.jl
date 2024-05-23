using JSON

function dot(u1, u2) 
    return sum(u1.*u2)
end

function norm(u)
    return sqrt(sum(u.^2))
end

function mean(v::Vector{Float64})
    return sum(v)/length(v)
end

function test_overlap(r::Array{Float64}, rho::Float64, sample_sz::Int64, L::Float64)
    current_distr = Array{Float64}([])
    for k in 1:sample_sz
        l::Array{Float64} = L.*rand(3,2)
        ul::Array{Float64} = (l[:,1]-l[:,2])./norm(l[:,1]-l[:,2])
        overlap = false
        d_perp::Float64 = 0.0
        for i in 1:size(r,2)
            rvec = r[:,i]
            d_pll2::Float64 = dot(ul, rvec-l[:,2])
            d_pll1::Float64 = dot(-ul, rvec-l[:,1])
            if ((d_pll2>0) && (d_pll1>0))        
                d_perp = norm((rvec-l[:,2]) - ul*dot(ul, rvec-l[:,2])) 
                if d_perp<rho
                    overlap = true
                    break
                end
            end
        end
        if overlap == false
            push!(current_distr,d_perp)
        end
    end
    return current_distr
end

#Obtiene el numero de frames en el dump
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
    r = Array{Float64}(undef,3,0)
    for i in 1:n
        line = split(readline(dump)," ")
        if line[2] == "2" || line[2] == "4"
            continue
        else
            ri = Vector{Float64}([parse(Float64,line[ii]) for ii in 3:5])
            r = hcat(r,ri)
        end
    end
    return r
end

#Aplica una traslacion diagonal para mover el centro a la esquina de la caja
function traslacion(r::Array{Float64},T::Float64)
    for i in 1:3
        for j in 1:size(r)[2]
            r[i,j] = r[i,j] + T
        end
    end
    return r
end

#Aplica el shear a todas las particulas
function shear(r::Array{Float64},xy::Float64,L::Float64)
    for i in 1:size(r)[2]
        r[1,i] = r[1,i] + (xy/L)*r[2,i]
    end
    return r
end

#Asegura que todas las particulas esten dentro de la caja
function wrap_boundaries(r::Array{Float64},L::Float64)
    for i in 1:3
        for j in 1:size(r)[2]
            if r[i,j] >= L 
                r[i,j] -= L 
            elseif r[i,j] < 0
                r[i,j] += L 
            end
        end
    end
    return r
end

#Hace el shear inverso, asegura que todas las particulas esten dentro de la caja y regresa el shear.
function fix_boundaries(r::Array{Float64},xy::Float64,L::Float64)
    r1 = traslacion(r,L/2)
    r2 = shear(r1,-xy,L)
    r3 = wrap_boundaries(r2,L)
    r4 = shear(r3,xy,L)
    r5 = wrap_boundaries(r4,L)
    return r5
end

#Guarda los datos en un archivo json
function write_json(savedir,data,name)
    json_string = JSON.json(data)
    open("$(savedir)/hole_analysis_line_$(name).json","w") do f
        JSON.print(f, json_string)
    end
end

#obtener numero de frames
dumpdir = ARGS[1] 

frame_num = read_frame_num(dumpdir)

calc_frames_num = 50
calc_frames_step = floor(frame_num/calc_frames_num)
calc_frames = 1:calc_frames_step:frame_num

sample_sz = 100000
rho = 2^(1/6)

dump = open(dumpdir,"r")

distr = Vector{Vector{Float64}}([])
ave = Vector{Float64}([])
hole_num = Vector{Float64}([])
for frame in 1:frame_num

    timestep,n,L,xy = read_dump_info(dump)

    #Analisis de huecos
    if frame in calc_frames
        println("$(frame) ")

        #Leer info de atomos y quitar patches
        r = read_dump_coords(dump,n)

        #Asegurar que todas las particulas esten dentro de la caja (condiciones periodicas) y mapeo a caja central. 
        r = fix_boundaries(r,xy,L)

        #Obtener distribucion de frame actual y guardar en distr
        current_distr = test_overlap(r,rho,sample_sz,L)
        push!(distr, current_distr)
        push!(ave, mean(current_distr))
        push!(hole_num, length(current_distr)/sample_sz)

    else
        for i in 1:n
            readline(dump)
        end
    end

end

dir = split(dumpdir,"/")
dir[end] = "hole_analysis_results"
savedir = join(dir,'/')

write_json(savedir,distr,"distr")
write_json(savedir,ave,"ave")
write_json(savedir,calc_frames,"calc_frames")
write_json(savedir,hole_num,"hole_num")


close(dump)