using Profile

function dot(u1, u2) 
    return sum(u1.*u2)
end

function norm(u)
    return sqrt(sum(u.^2))
end

function test_overlap(r::Array{Float64}, rho::Float64, sample_sz::Int64, L::Float64)
    current_distr = Vector{Float64}([])
    l = Array{Float64,2}(undef,3,2)
    ul = Vector{Float64}(undef,3)
    @views for k in 1:sample_sz
        l .= L.*rand(3,2)
        ul .= (l[:,1]-l[:,2])./norm(l[:,1]-l[:,2])
        overlap = false
        d_perp::Float64 = 0.0
        for i in 1:size(r,2)
            d_pll2::Float64 = dot(ul, r[:,i]-l[:,2])
            d_pll1::Float64 = dot(-ul, r[:,i]-l[:,1])
            if ((d_pll2>0) && (d_pll1>0))        
                d_perp = norm((r[:,i]-l[:,2]) - ul*dot(ul, r[:,i]-l[:,2])) 
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

L = 54.0
r = L.*rand(3,16000)
rho = 2^(1/6)
sample_sz = 10000

test_overlap(r,rho,sample_sz,L)

@time a = test_overlap(r,rho,sample_sz,L)
println(length(a))

#@profile test_overlap(r,rho,sample_sz,L)
#Profile.print()