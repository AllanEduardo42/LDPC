################################################################################
# Allan Eduardo Feitosa
# 19 set 2024
# Function to calculate the prior probabilities

function
    calc_Lf!(
        Lf::Vector{Float64},
        twoZc::Int,
        signal::Vector{Float64},
        σ²::Float64
    )

    @fastmath @inbounds begin
        for i in eachindex(signal)
            Lf[twoZc+i] = -2*signal[i]/σ²
        end
    end

end

# MKAY
function
    calc_Lf!(
        f::Matrix{Float64},
        twoZc::Int,
        signal::Vector{Float64},
        σ²::Float64
    )

    @fastmath @inbounds begin

        k = 1/sqrt(2*π*σ²)

        for i in eachindex(signal)
            f[twoZc+i,1] = k*exp(-(signal[i]+1)^2/(2*σ²))
            f[twoZc+i,2] = k*exp(-(signal[i]-1)^2/(2*σ²))
        end

        normalize!(f)
    end

end

function normalize!(f::AbstractMatrix{Float64})
    
    @inbounds @fastmath for n in axes(f,1)
        α = f[n,1] + f[n,2]
        f[n,1] = f[n,1]/α
        f[n,2] = f[n,2]/α
    end

end