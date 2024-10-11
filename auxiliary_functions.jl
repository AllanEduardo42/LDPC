################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Auxiliary functions

using Random

"""For each column j of matrix H, find the indices i where H[i,j] = 1, and 
return a vector "vn2cn" where vn2cn[j] is a vector containing the
indices i where H[i,j] = 1"""
function make_vn2cn_list(H::BitMatrix)

    N = size(H,2)
    vn2cn = Vector{Vector{Int}}()
    for n in 1:N
        @inbounds push!(vn2cn,findall(x -> x == true, H[:,n]))
    end

    return vn2cn

end

"""For each row i of matrix H, find the indices j where H[i,j] = 1, and 
return a vector "cn2vn" where cn2vn[i] is a vector containing the
indices j where H[i,j] = 1"""
function make_cn2vn_list(H::BitMatrix)

    M = size(H,1)
    cn2vn = Vector{Vector{Int}}()
    for m in 1:M
        @inbounds push!(cn2vn,findall(x -> x == true, H[m,:]))
    end

    return cn2vn

end

function 
    init_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )
    
    for n in eachindex(vn2cn)
        for m in vn2cn[n]
            @inbounds Lq[n,m] = Lf[n]
        end
    end
end

function
    init_Lq!(
        Lq::Array{<:AbstractFloat,3},
        Lf::Matrix{<:AbstractFloat},
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )
   
    for n in eachindex(vn2cn)
        for m in vn2cn[n]
            @inbounds Lq[n,m,1] = Lf[n,1]
            @inbounds Lq[n,m,2] = Lf[n,2]
        end
    end

end

function
    received_signal!(
        t::Vector{<:AbstractFloat},
        noise::Vector{<:AbstractFloat},
        σ::AbstractFloat,
        u::Vector{<:AbstractFloat},
        rng_noise::AbstractRNG)

    randn!(rng_noise,noise)
    @fastmath noise .*= σ
    @fastmath t .= u .+ noise

end

function
    resetfactors!(
        Factors::Matrix{<:AbstractFloat},
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )
    
    for n in eachindex(vn2cn)
        for m in vn2cn[n]
            Factors[m,n] = 1.0
        end
    end
end