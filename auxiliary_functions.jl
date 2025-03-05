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
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )
    
    @inbounds for n in eachindex(vn2cn)
        for m in vn2cn[n]
            Lq[n,m] = Lf[n]
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
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        value::Real
    )    
    
    @fastmath @inbounds begin
        for n in 1:N
            nl = LinearIndices(X)[1,n]-1
            for m in vn2cn[n]
                X[nl+m] = value
            end
        end
    end
end