using JSON
using StaticArrays

function cdot(u1::SVector{3,Float64}, u2::SVector{3,Float64}) 
    sum(u1.*u2)
end

function norm(u::SVector{3,Float64})
    return sqrt(cdot(u,u))
end

function mean(v::Vector{Float64})
    return sum(v)/length(v)
end

@views function check_overlap(r::SVector{3,Float64},l1::SVector{3,Float64},l2::SVector{3,Float64},rho::Float64)
    ul = SVector{3,Float64}((l1-l2)/norm(l1-l2))
    ri1 = SVector{3,Float64}(r-l1)
    ri2 = SVector{3,Float64}(r-l2)
    d_pll2::Float64 = cdot(ul, ri2)
    d_pll1::Float64 = cdot(-ul, ri1)
    overlap = false
    if ((d_pll2>0) && (d_pll1>0))        
        d_perp = norm(ri2 - ul*d_pll2) 
        if d_perp<rho
            overlap = true
        end
    end
    return overlap
end

@views function line_distr(r::Vector{SVector{3,Float64}},sample_sz::Int64,L::Float64,rho::Float64)
    distr = Vector{Float64}()
    sample_distr = Vector{Float64}(undef,sample_sz)
    for k in 1:sample_sz
        l1 = L* @SVector rand(Float64,3)
        l2 = L* @SVector rand(Float64,3)
        sample_distr[k] = norm(l1-l2)
        overlap = true
        for i in 1:size(r,1)
            overlap = check_overlap(r[i],l1,l2,rho)
            if overlap
                break
            end
        end
        if !overlap 
            push!(distr,norm(l1-l2))
        end
    end
    return distr,sample_distr
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

sample_sz = 1_000_000
rho = 2^(-5/6)

dump = open(dumpdir,"r")

distr = Vector{Vector{Float64}}()
ave = Vector{Float64}()
hole_num = Vector{Float64}()
sample_distr = Vector{Vector{Float64}}()
for frame in 1:frame_num

    timestep,n,L,xy = read_dump_info(dump)

    #Analisis de huecos
    if frame in calc_frames
        println("Frame:","$(frame) ")

        #Leer info de atomos y quitar patches
        println("Leyendo datos...")
        r = read_dump_coords(dump,n)

        #Asegurar que todas las particulas esten dentro de la caja (condiciones periodicas) y mapeo a caja central. 
        r = fix_boundaries(r,xy,L)

        #Obtener distribucion de frame actual y guardar en distr
        println("Calculando...")
        current_distr,current_sample_distr = line_distr(r,sample_sz,L,rho)
        push!(distr,current_distr)
        push!(ave,mean(current_distr))
        push!(hole_num,length(current_distr)/sample_sz)
        push!(sample_distr,current_sample_distr)
        println("Fin")
    else
        for _ in 1:n
            readline(dump)
        end
    end

end

dir = split(dumpdir,"/")
dir[end] = "hole_analysis_results"
savedir = join(dir,"/")

write_json(savedir,distr,"distr")
write_json(savedir,ave,"ave")
write_json(savedir,calc_frames,"calc_frames")
write_json(savedir,hole_num,"hole_num")
write_json(savedir,sample_distr,"sample_distr")


close(dump)