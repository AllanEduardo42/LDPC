################################################################################
# Allan Eduardo Feitosa
# 19 set 2024
# Function to calculate the prior probabilities

function
    calc_Lf!(
        Lf::Vector{<:AbstractFloat},
        signal::Vector{<:AbstractFloat},
        σ²::AbstractFloat
    )

    for i in eachindex(signal)
        @inbounds @fastmath Lf[i] = -2*signal[i]/σ²
    end

end

function
    calc_Lf!(
        f::Matrix{<:AbstractFloat},
        signal::Vector{<:AbstractFloat},
        σ²::AbstractFloat
    )

    k = 1/sqrt(2*π*σ²)

    for i in eachindex(signal)
        @inbounds @fastmath f[i,1] = k*exp(-(signal[i]+1)^2/(2*σ²))
        @inbounds @fastmath f[i,2] = k*exp(-(signal[i]-1)^2/(2*σ²))
    end

    normalize!(f)

end

function normalize!(f::Matrix{<:AbstractFloat})

    N = size(f,1)
    
    for n in 1:N
        @inbounds @fastmath α = f[n,1] + f[n,2]
        @inbounds @fastmath f[n,1] = f[n,1]/α
        @inbounds @fastmath f[n,2] = f[n,2]/α
    end

end