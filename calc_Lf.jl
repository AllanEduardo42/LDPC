################################################################################
# Allan Eduardo Feitosa
# 19 set 2024
# Function to calculate the prior probabilities

function
    calc_Lf!(
        Lf::AbstractVector{<:AbstractFloat},
        signal::Vector{<:AbstractFloat},
        σ²::AbstractFloat,
        bptype::String
    )

    @fastmath @inbounds for i in eachindex(signal)
        Lf[i] = -2*signal[i]/σ²
    end
    if bptype == "TABL"
        # scale for table
        Lf .*= SIZE_per_RANGE
    end

end

# MKAY
function
    calc_Lf!(
        f::Matrix{<:AbstractFloat},
        signal::Vector{<:AbstractFloat},
        σ²::AbstractFloat
    )

    @fastmath k = 1/sqrt(2*π*σ²)

    @inbounds @fastmath  for i in eachindex(signal)
        f[i,1] = k*exp(-(signal[i]+1)^2/(2*σ²))
        f[i,2] = k*exp(-(signal[i]-1)^2/(2*σ²))
    end

    normalize!(f)

end

function normalize!(f::Matrix{<:AbstractFloat})

    N = size(f,1)
    
    @inbounds @fastmath for n in 1:N
        α = f[n,1] + f[n,2]
        f[n,1] = f[n,1]/α
        f[n,2] = f[n,2]/α
    end

end