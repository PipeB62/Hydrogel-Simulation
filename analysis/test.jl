using Profile
using StaticArrays
using BenchmarkTools

function cdot(u1::SVector{3,Float64}, u2::SVector{3,Float64}) 
    sum(u1.*u2)
end

function norm(u::SVector{3,Float64})
    return sqrt(cdot(u,u))
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
    for _ in 1:sample_sz
        l1 = L* @SVector rand(Float64,3)
        l2 = L* @SVector rand(Float64,3)
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
    return distr
end

L = 54.0
rho = 2^(1/6)
sample_sz = 1000000
r = [L*rand(SVector{3, Float64}) for i in 1:16000]

@btime line_distr($r,$sample_sz,$L,$rho)
#@code_warntype line_distr(r,sample_sz,L,rho)
#@profview line_distr(r,sample_sz,L,rho)

#l1 = L* @SVector rand(Float64,3)
#l2 = L* @SVector rand(Float64,3)

#@code_warntype check_overlap(r[1],l1,l2,rho)




