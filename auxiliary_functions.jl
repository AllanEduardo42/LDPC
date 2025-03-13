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
    return make_vn2cn_list(H,N)

end

function make_vn2cn_list(H::BitMatrix, N::Integer)

    vn2cn = Vector{Vector{Int}}()
    @inbounds for n in 1:N
        push!(vn2cn,findall(isone, H[:,n]))
    end

    return vn2cn

end

"""For each row i of matrix H, find the indices j where H[i,j] = 1, and 
return a vector "cn2vn" where cn2vn[i] is a vector containing the
indices j where H[i,j] = 1"""
function make_cn2vn_list(H::BitMatrix)

    M = size(H,1)
   
    return make_cn2vn_list(H,M)

end

function make_cn2vn_list(H::BitMatrix, M::Integer)

    cn2vn = Vector{Vector{Int}}()
    @inbounds for m in 1:M
        push!(cn2vn,findall(isone, H[m,:]))
    end

    return cn2vn

end

function 
    init_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        M::Integer,
        N::Integer
    )

    m = 0
    @inbounds for ml in 0:N:M*(N-1)
        m += 1
        for n in cn2vn[m]
            Lq[ml+n] = Lf[n]
        end
    end
end

function
    init_Lq!(
        Lq::Array{<:AbstractFloat,3},
        Lf::Matrix{<:AbstractFloat},
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )
   
    @inbounds for n in eachindex(vn2cn)
        for m in vn2cn[n]
            Lq[n,m,1] = Lf[n,1]
            Lq[n,m,2] = Lf[n,2]
        end
    end

end

function
    received_signal!(
        signal::AbstractArray{<:AbstractFloat},
        noise::Vector{<:AbstractFloat},
        σ::AbstractFloat,
        u::Vector{<:AbstractFloat},
        rng_noise::AbstractRNG)

    randn!(rng_noise,noise)
    @fastmath noise .*= σ
    @fastmath signal .= u .+ noise

end

function
    resetmatrix!(
        X::Matrix{<:Real},
        M::Integer,
        N::Integer,
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        value::Real
    )    
    
    n = 0
    @inbounds for nl in 0:M:N*(M-1)
        n += 1
        for m in vn2cn[n]
            X[nl+m] = value
        end
    end
end