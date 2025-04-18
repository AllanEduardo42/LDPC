################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Auxiliary functions

using Random

"""For each column j of matrix H, find the indices i where H[i,j] = 1, and 
return a vector "vn2cn" where vn2cn[j] is a vector containing the
indices i where H[i,j] = 1"""
function make_vn2cn_list(H::Matrix{Bool})

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
function make_cn2vn_list(H::Matrix{Bool})

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
            Lq[m,n] = Lf[n]
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
            Lq[m,n,1] = Lf[n,1]
            Lq[m,n,2] = Lf[n,2]
        end
    end

end

function
    received_signal!(
        signal::AbstractArray{<:AbstractFloat},
        noise::Vector{<:AbstractFloat},
        σ::AbstractFloat,
        u::Vector{<:AbstractFloat},
        rng_noise::AbstractRNG,
        noisetest::Nothing
    )

    @fastmath begin
        randn!(rng_noise,noise)
        noise .*= σ
        signal .= u .+ noise
    end

end

function
    received_signal!(
        signal::AbstractArray{<:AbstractFloat},
        noise::Vector{<:AbstractFloat},
        σ::AbstractFloat,
        u::Vector{<:AbstractFloat},
        rng_noise::AbstractRNG,
        noisetest::Vector{<:AbstractFloat}
    )

    @fastmath begin
        noisetest .*= σ
        signal .= u .+ noisetest
    end

end

function
    resetmatrix!(
        X::Matrix{<:Real},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        value::Real
    )    
    
    @inbounds for n in eachindex(vn2cn)
        for m in vn2cn[n]
            X[m,n] = value
        end
    end
end

function
    resetmatrix!(
        X::Array{<:Real,3},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        value::Real
    )    
    
    @inbounds for n in eachindex(vn2cn)
        for m in vn2cn[n]
            X[m,n,1] = value
            X[m,n,2] = value
        end
    end
end

function generate_message!(
    msg::Vector{Bool},
    rgn_msg::AbstractRNG,
    ::Nothing
)

    rand!(rgn_msg,msg,Bool)

end

function generate_message!(
    msg::Vector{Bool},
    ::AbstractRNG,
    msgtest::Vector{Bool}
)

    msg .= msgtest

end

function print_test(
    text::String,
    array::Vector{Bool}

)    
    println()
    println("$text (L = $(length(array))):")
    for i in eachindex(array)
        print(Int(array[i]))
        if i%80 == 0
            println()
        end
    end
    println()
end 