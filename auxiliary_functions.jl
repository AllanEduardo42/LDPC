################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Auxiliary functions

"""For each column j of matrix H, find the indices i where H[i,j] = 1, and 
return a vector "vn2cn" where vn2cn[j] is a vector containing the
indices i where H[i,j] = 1"""
function make_vn2cn_list(H::BitMatrix)

    N = size(H,2)
    vn2cn = Vector{Vector{Int}}()
    for vn in 1:N
        push!(vn2cn,findall(x -> x == true, H[:,vn]))
    end

    return vn2cn

end

"""For each row i of matrix H, find the indices j where H[i,j] = 1, and 
return a vector "cn2vn" where cn2vn[i] is a vector containing the
indices j where H[i,j] = 1"""
function make_cn2vn_list(H::BitMatrix)

    M = size(H,1)
    cn2vn = Vector{Vector{Int}}()
    for cn in 1:M
        push!(cn2vn,findall(x -> x == true, H[cn,:]))
    end

    return cn2vn

end

function 
    init_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )
    
    for vn in eachindex(vn2cn)
        for cn in vn2cn[vn]
            @inbounds Lq[vn,cn] = Lf[vn]
        end
    end
end

function
    init_Lq!(
        Lq::Array{<:AbstractFloat,3},
        Lf::Matrix{<:AbstractFloat},
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )
   
    for vn in eachindex(vn2cn)
        for cn in vn2cn[vn]
            @inbounds Lq[vn,cn,1] = Lf[vn,1]
            @inbounds Lq[vn,cn,2] = Lf[vn,2]
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
    reset_factors!(
        Factors::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )

    for cn in eachindex(cn2vn)
        for vn in cn2vn[cn]
            Factors[cn,vn] = 1.0
        end
    end
end