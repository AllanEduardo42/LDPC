################################################################################
# Allan Eduardo Feitosa
# 19 set 2024
# Function to calculate the prior probabilities

function
    calc_Lf!(
        Lf::Vector{<:AbstractFloat},
        t::Vector{<:AbstractFloat},
        σ²::AbstractFloat
    )

    for i in eachindex(t)
        @inbounds @fastmath Lf[i] = -2*t[i]/σ²
    end

end

function
    calc_f!(
        f::Matrix{<:AbstractFloat},
        t::Vector{<:AbstractFloat},
        σ²::AbstractFloat
    )

    k = 1/sqrt(2*π*σ²)

    for i in eachindex(t)
        @inbounds @fastmath f[i,1] = k*exp(-(t[i]+1)^2/(2*σ²))
        @inbounds @fastmath f[i,2] = k*exp(-(t[i]-1)^2/(2*σ²))
    end

    normalize!(f)


end

"""Normalize a binary probability distribution function"""
function normalize!(f::Matrix{<:AbstractFloat})

    N = size(f,1)
    
    for n in 1:N
        @inbounds @fastmath α = f[n,1] + f[n,2]
        @inbounds @fastmath f[n,1] = f[n,1]/α
        @inbounds @fastmath f[n,2] = f[n,2]/α
    end

end