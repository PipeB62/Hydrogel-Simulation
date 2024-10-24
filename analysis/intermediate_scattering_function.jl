using JSON
using StaticArrays
using LinearAlgebra
using NearestNeighbors

include("analysisModule.jl")
using .analysisModule

function main()

    dumpdir = ARGS[1]
    savedir = ARGS[2]
    plano = ARGS[3]

    println("Inicio ISF")

    frame_num = read_frame_num(dumpdir)

    calc_frames_num = frame_num/2
    calc_frames_step = floor(frame_num/calc_frames_num)
    calc_frames = 2:calc_frames_step:frame_num

    dump = open(dumpdir,"r") #Abrir dump

    #Leer primer frame (primer frame de referencia)

    timestep,n,L,r_xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

    #Caja central
    a_0 = SVector{3,Float64}([L,0,0])
    b_0 = SVector{3,Float64}([0,L,0])
    c_0 = SVector{3,Float64}([0,0,L])

    #Triclinic box vectors
    a = SVector{3,Float64}([L,0,0])
    b = SVector{3,Float64}([r_xy,L,0])
    c = SVector{3,Float64}([0,0,L])

    #Referencia. t=0
    r_centers_coords,r_centers_id,_,_ = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
    r_centers_coords = fix_boundaries(r_centers_coords,L,a,b,c) #Asegurar que todas las particulas esten dentro de la caja
    r_centers_coords = shear(r_centers_coords,-r_xy,L) #Revertir shear
    r_centers_coords = wrap_boundaries(r_centers_coords,a_0,b_0,c_0) #Mapear a caja central

    N_p = length(r_centers_coords) #Numero de particulas (centros)

    #Vectores normales al plano de interes

    if plano == "xy"
        n_plano = SVector{3,Float64}([0,0,1]) #Plano xy
    elseif plano == "xz"
        n_plano = SVector{3,Float64}([0,1,0]) #Plano xz
    elseif plano == "yz"
        n_plano = SVector{3,Float64}([1,0,0]) #Plano yz
    else
        n_plano = SVector{3,Float64}([0,0,0]) #Sin proyectar
    end

    P_op = I-outer_product(n_plano,n_plano) #Operador de proyeccion sobre plano

    #sampleo de vectores k
    k = (2*pi/L) * 10
    delta_k = 0.01*k 
    k_max = k + delta_k/2
    k_min = k - delta_k/2
    n_max = ceil(L*k_max/(2*pi))
    #println(k_min," ",k," ",k_max)
    #println("n_max: ",n_max)

    k_vectors = Vector{SVector{3,Float64}}()

    for n_1 in 0:n_max 
        for n_2 in 0:n_max 
            for n_3 in 0:n_max
                kvec = (2*pi/L)*SVector{3,Float64}([n_1,n_2,n_3])
                kvec_norm = norm(kvec)
                if kvec_norm<=k_max && kvec_norm>=k_min && cdot(kvec,n_plano)==0
                    push!(k_vectors,kvec)
                end
            end
        end
    end
    N_k = length(k_vectors)
    println("N_k: ",N_k)

    F = Vector{Float64}() #|F(k,t)| como funcion del tiempo
    p_xy = r_xy
    flipcount = 0
    for frame in 2:frame_num

        timestep,n,L,c_xy = read_dump_info(dump) #Leer informacion del dump en el frame actual

        #Obtener el tilt real (sin flip)
        if abs(c_xy-p_xy)>1
            flipcount+=1
        end
        real_xy =L*flipcount+c_xy

        if frame in calc_frames

            println("Frame: ",frame)
            #println("c_xy: ",c_xy)
            #println("real_xy: ",real_xy)

            #Triclinic box vectors
            a = SVector{3,Float64}([L,0,0])
            b = SVector{3,Float64}([c_xy,L,0])
            c = SVector{3,Float64}([0,0,L])

            c_centers_coords,c_centers_id,_,_ = read_dump_particles(dump,n) #Leer coordenadas y id del dump en el frame actual
            c_centers_coords = fix_boundaries(c_centers_coords,L,a,b,c) #Asegurar que todas las particulas esten dentro de la caja
            c_centers_coords = shear(c_centers_coords,-real_xy,L) #Revertir shear
            c_centers_coords = wrap_boundaries(c_centers_coords,a_0,b_0,c_0) #Mapeo a caja central

            Fkt = 0
            test = 0
            for kvec in k_vectors

                reF = 0
                imF = 0
                for i in eachindex(c_centers_id)
                    iid = c_centers_id[i][1] #id de particula

                    rt = c_centers_coords[i] #posicion en tiempo t
                    r0 = get_coords(iid,r_centers_id,r_centers_coords) #posicion en tiempo 0

                    delta_r_all = Vector{SVector{3,Float64}}(undef,27)
                    aa = 1
                    for q in -1:1
                        for p in -1:1
                            for o in -1:1
                                tr = o*a_0 + p*b_0 + q*c_0
                                rt_tr = rt + tr
                                delta_r_all[aa] = rt_tr-r0
                                aa+=1
                            end
                        end
                    end

                    delta_r = argmin(x->norm(x),delta_r_all) 
                    delta_r_ab = SVector{3,Float64}(P_op*delta_r)

                    #=
                    if iid==r_centers_id[1832][1] && test==0
                        println("rt ",rt)
                        println("ro ",r0)
                        println("deltar ",delta_r_ab)
                    end
                    =#
                    reF += cos(cdot(kvec,delta_r_ab))
                    imF += sin(cdot(kvec,delta_r_ab))
                end
                Fkt += (1/N_p)*sqrt(reF^2 + imF^2) #Calcular |F(kvec,t)|

                test+=1
            end
            Fkt = Fkt/N_k #Promediar sobre valores de F para todos los vectores k

            push!(F,Fkt)

        else
            for _ in 1:n
                readline(dump)
            end
        end

        p_xy = c_xy

    end

    close(dump)

    write_json(savedir,F,"ISF_$(plano)")
    write_json(savedir,calc_frames,"ISF_calcframes_$(plano)")

end

main()